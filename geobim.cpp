/********************************************************************************
 *                                                                              *
 * This file is part of TUDelft Esri GEOBIM.                                    *
 *                                                                              *
 * License: APACHE                                                              *
 *                                                                              *
 ********************************************************************************/

 // Example invocations:
 // -o --entities=IfcWall;IfcSlab;IfcWindow;IfcDoor x.ifc y 0.01 0.1
 // --entities=IfcWall;IfcSlab x y 0.01 0.1

#define ENSURE_2ND_OP_NARROWER

#include <ifcgeom/kernels/cgal/CgalKernel.h>
#include <ifcgeom/schema_agnostic/IfcGeomFilter.h>
#include <ifcgeom/schema_agnostic/IfcGeomIterator.h>
#include <ifcconvert/validation_utils.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <CGAL/minkowski_sum_3.h>
#include <CGAL/Nef_nary_union_3.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/boost/graph/copy_face_graph.h>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>

#include <nlohmann/json.hpp>

#include <chrono>

template <class Poly_B, class Poly_A>
typename std::enable_if<std::is_same<Poly_A, Poly_B>::value>::type poly_copy(Poly_B& poly_b, const Poly_A& poly_a) {
	poly_b = poly_a;
}

template <class Poly_B, class Poly_A>
typename std::enable_if<!std::is_same<Poly_A, Poly_B>::value>::type poly_copy(Poly_B& poly_b, const Poly_A& poly_a) {
	poly_b.clear();
	CGAL::copy_face_graph(poly_a, poly_b);
}

template <class Kb, class Ka>
typename std::enable_if<std::is_same<Ka, Kb>::value>::type transformation_copy(CGAL::Aff_transformation_3<Kb>& trsf_b, const CGAL::Aff_transformation_3<Ka>& trsf_a) {
	trsf_b = trsf_a;
}

template <class Kb, class Ka>
typename std::enable_if<!std::is_same<Ka, Kb>::value>::type transformation_copy(CGAL::Aff_transformation_3<Kb>& trsf_b, const CGAL::Aff_transformation_3<Ka>& trsf_a) {
	CGAL::NT_converter<
		typename Ka::RT, 
		typename Kb::RT> converter;
	trsf_b = CGAL::Aff_transformation_3<Kb>(
		converter(trsf_a.hm(0, 0)),
		converter(trsf_a.hm(0, 1)),
		converter(trsf_a.hm(0, 2)),
		converter(trsf_a.hm(0, 3)),
		converter(trsf_a.hm(1, 0)),
		converter(trsf_a.hm(1, 1)),
		converter(trsf_a.hm(1, 2)),
		converter(trsf_a.hm(1, 3)),
		converter(trsf_a.hm(2, 0)),
		converter(trsf_a.hm(2, 1)),
		converter(trsf_a.hm(2, 2)),
		converter(trsf_a.hm(2, 3)),
		converter(trsf_a.hm(3, 3))
	);
}


namespace po = boost::program_options;

class timer_class {
	typedef std::map<std::string, double> timings_map;
	timings_map timings_;

public:
	// @todo something like a scoped_stopwatch using RAII
	class stopwatch {
		typedef std::chrono::time_point<std::chrono::steady_clock> timepoint;
		std::string key;
		timepoint start_, stop_;
		timings_map& timings;

	public:
		stopwatch(timings_map& timings) : timings(timings) {}
		stopwatch& start(const std::string& s) {
			key = s;
			start_ = std::chrono::steady_clock::now();
			return *this;
		}
		void stop() {
			stop_ = std::chrono::steady_clock::now();
			timings[key] += std::chrono::duration<double>(stop_ - start_).count();
		}
	};

	stopwatch measure(const std::string& s) {
		return stopwatch(timings_).start(s);
	}

	void print(std::ostream& s) const {
		size_t max_key_length = 0;
		for (auto& p : timings_) {
			if (p.first.size() > max_key_length) {
				max_key_length = p.first.size();
			}
		}
		for (auto& p : timings_) {
			s << p.first << std::string(max_key_length - p.first.size(), ' ') << ":" << p.second << std::endl;
		}
	}
};

static timer_class timer;

struct abstract_writer {
	std::array<Kernel_::Point_3, 3> points_from_facet(cgal_shape_t::Facet_handle f) {
		return {
				f->facet_begin()->vertex()->point(),
				f->facet_begin()->next()->vertex()->point(),
				f->facet_begin()->next()->next()->vertex()->point()
		};
	}

	std::array<CGAL::Simple_cartesian<double>::Point_3, 3> points_from_facet(CGAL::Polyhedron_3<CGAL::Simple_cartesian<double>>::Facet_handle f) {
		return {
				f->facet_begin()->vertex()->point(),
				f->facet_begin()->next()->vertex()->point(),
				f->facet_begin()->next()->next()->vertex()->point()
		};
	}

	std::array<Kernel_::Point_3, 3> points_from_facet(std::list<cgal_shape_t::Facet_handle>::iterator f) {
		return points_from_facet(*f);
	}
};

// OBJ writer for CGAL facets paired with a style
struct simple_obj_writer : public abstract_writer {
	int group_id = 1;
	int vertex_count = 1;
	std::ofstream obj, mtl;
	ifcopenshell::geometry::taxonomy::colour GRAY;

	simple_obj_writer(const std::string& fn_prefix)
		: obj((fn_prefix + ".obj").c_str())
		, mtl((fn_prefix + ".mtl").c_str())
		, GRAY(0.6, 0.6, 0.6)
	{
		obj << "mtllib " << fn_prefix << ".mtl\n";
	}
	
	template <typename It>
	void operator()(const ifcopenshell::geometry::taxonomy::style* style, It begin, It end) {
		const auto& diffuse = *(style ? style->diffuse.get_value_or(GRAY).components : GRAY.components);

		obj << "g group-" << group_id << "\n";
		obj << "usemtl m" << group_id << "\n";
		mtl << "newmtl m" << group_id << "\n";
		mtl << "kd " << diffuse(0) << " " << diffuse(1) << " " << diffuse(2) << "\n";

		group_id++;

		for (auto it = begin; it != end; ++it) {
			auto points = points_from_facet(it);
			for (int i = 0; i < 3; ++i) {
				obj << "v "
					<< points[i].cartesian(0) << " "
					<< points[i].cartesian(1) << " "
					<< points[i].cartesian(2) << "\n";
			}
			obj << "f "
				<< (vertex_count + 0) << " "
				<< (vertex_count + 1) << " "
				<< (vertex_count + 2) << "\n";
			vertex_count += 3;
		}
	}
};

struct city_json_writer : public abstract_writer {
	ifcopenshell::geometry::taxonomy::colour GRAY;

	using json = nlohmann::json;

	std::string filename;

	std::vector<std::array<double, 3>> vertices;
	std::vector<std::vector<std::vector<std::vector<int>>>> boundaries;
	std::vector<std::vector<int>> boundary_materials;

	json materials;

	city_json_writer(const std::string& fn_prefix)
		: filename(fn_prefix + ".json")
		, materials(json::array())
		, GRAY(0.6, 0.6, 0.6)
	{
		// assumes one solid.
		boundaries.emplace_back();
		boundary_materials.emplace_back();
	}

	template <typename It>
	void operator()(const ifcopenshell::geometry::taxonomy::style* style, It begin, It end) {
		const auto& diffuse = *(style ? style->diffuse.get_value_or(GRAY).components : GRAY.components);

		json material = json::object();
		material["name"] = "material-" + boost::lexical_cast<std::string>(materials.size());
		material["diffuseColor"] = std::array<double, 3>{diffuse(0), diffuse(1), diffuse(2)};
		material["specularColor"] = std::array<double, 3>{0., 0., 0.};
		material["shininess"] = 0.;
		material["isSmooth"] = false;
		materials.push_back(material);

		for (auto it = begin; it != end; ++it) {
			auto points = points_from_facet(it);
			std::vector<int> faces;
			for (int i = 0; i < 3; ++i) {
				faces.push_back(vertices.size());
				vertices.push_back({{
					CGAL::to_double(points[i].cartesian(0)),
					CGAL::to_double(points[i].cartesian(1)),
					CGAL::to_double(points[i].cartesian(2))
				}});
			}
			boundaries.front().push_back({ faces });
			boundary_materials.front().push_back(materials.size() - 1);
		}
	}

	void finalize() {
		json city;

		city["type"] = "CityJSON";
		city["version"] = "1.0";
		city["extensions"] = json::object();
		city["metadata"]["referenceSystem"] = "urn:ogc:def:crs:EPSG::2355";
		city["vertices"] = vertices;
		city["appearance"]["materials"] = materials;

		auto& building1 = city["CityObjects"]["id-1"];
		building1["type"] = "Building";
		building1["geographicalExtent"] = std::array<double, 6>{0, 0, 0, 1, 1, 1};
		/*
		building1["attributes"]["measuredHeight"] = 22.3;
		building1["attributes"]["roofType"] = "gable";
		building1["attributes"]["owner"] = "Elvis Presley";
		*/

		json geom = json::object();
		geom["type"] = "Solid";
		geom["lod"] = 2;
		geom["boundaries"] = boundaries;
		geom["material"]["diffuse"]["values"] = boundary_materials;
		building1["geometry"].push_back(geom);

		std::ofstream(filename.c_str()) << city;
	}

	~city_json_writer() {
		finalize();
	}
};

// Global settings that are derived from the command line invocation
// parsed by Boost.ProgramOptions
struct geobim_settings {
	std::string input_filename, output_filename;
	std::vector<double> radii;
	bool apply_openings, apply_openings_posthoc, debug, exact_segmentation, minkowski_triangles;
	ifcopenshell::geometry::settings settings;
	boost::optional<std::set<std::string>> entity_names;
	bool entity_names_included;
	IfcParse::IfcFile* file;
};

// Generated for every representation *item* in the IFC file
struct shape_callback_item {
	IfcUtil::IfcBaseEntity* src;
	std::string id, type, geom_reference;
	cgal_placement_t transformation;
	cgal_shape_t polyhedron;
	boost::optional<ifcopenshell::geometry::taxonomy::style> style;
	boost::optional<Eigen::Vector3d> wall_direction;
	std::list<shape_callback_item*> openings;

	bool to_nef_polyhedron(CGAL::Nef_polyhedron_3<Kernel_>& nef) {
		auto T1 = timer.measure("ifc_element_to_nef");
		nef = ifcopenshell::geometry::utils::create_nef_polyhedron(polyhedron);
		T1.stop();

		if (nef.is_empty() || !nef.is_simple()) {
			nef.clear();
			return false;
		}

		return true;
	}
};

// Prototype of a context to which processed shapes will be fed
struct execution_context {
	virtual void operator()(shape_callback_item&) = 0;
};

#include <CGAL/Polygon_mesh_processing/self_intersections.h>

struct opening_collector : public execution_context {
	std::list<shape_callback_item> list;
	std::map<IfcUtil::IfcBaseEntity*, IfcUtil::IfcBaseEntity*> opening_to_elem;
	std::multimap<IfcUtil::IfcBaseEntity*, shape_callback_item*> map;

	opening_collector(IfcParse::IfcFile* f) {
		auto rels = f->instances_by_type("IfcRelVoidsElement");
		if (!rels) {
			return;
		}
		for (auto& rel : *rels) {
			auto be = (IfcUtil::IfcBaseEntity*) ((IfcUtil::IfcBaseEntity*)rel)->get_value<IfcUtil::IfcBaseClass*>("RelatingBuildingElement");
			auto op = (IfcUtil::IfcBaseEntity*) ((IfcUtil::IfcBaseEntity*)rel)->get_value<IfcUtil::IfcBaseClass*>("RelatedOpeningElement");
			opening_to_elem.insert({ op, be });
		}
	}

	void operator()(shape_callback_item& item) {
		auto opit = opening_to_elem.find(item.src);
		if (opit != opening_to_elem.end()) {
			list.push_back(item);
			map.insert({ opit->second, &list.back() });
		}
	}
};

struct debug_writer : public execution_context {
	void operator()(shape_callback_item& item) {
		simple_obj_writer obj("debug-" + boost::lexical_cast<std::string>(item.src->data().id()));
		obj(nullptr, item.polyhedron.facets_begin(), item.polyhedron.facets_end());
	}
};

// State that is relevant accross the different radii for which the
// program is executed. Mostly related for the mapping back to
// semantics from the original IFC input.
template <typename TreeKernel>
struct global_execution_context : public execution_context {
	typedef CGAL::Polyhedron_3<TreeKernel> TreeShapeType;
	typedef CGAL::AABB_face_graph_triangle_primitive<TreeShapeType, CGAL::Default, CGAL::Tag_false> Primitive;
	typedef CGAL::AABB_traits<TreeKernel, Primitive> AAbbTraits;
	typedef CGAL::AABB_tree<AAbbTraits> AAbbTree;

	AAbbTree tree;

	std::list<ifcopenshell::geometry::taxonomy::style> styles;
	std::map<std::pair<double, std::pair<double, double>>, typename decltype(styles)::iterator> diffuse_to_style;

	// A reference is kept to the original shapes in a std::list.
	// Later an aabb tree is used map eroded triangle centroids
	// back to the original elements to preserve semantics.
	std::list<TreeShapeType> triangulated_shape_memory;
	std::map<typename TreeShapeType::Facet_handle, typename decltype(styles)::const_iterator> facet_to_style;

#ifdef GEOBIM_DEBUG
	simple_obj_writer obj_;
#endif

	global_execution_context()
#ifdef GEOBIM_DEBUG
		: obj_("debug-tree")
#endif
	{
		// style 0 is for elements without style annotations
		styles.emplace_back();
	}

	void operator()(shape_callback_item& item) {
		size_t style_idx = item.style ? styles.size() : 0;

		// Group taxonomy::styles based on diffuse colour since we
		// do not have an equality operator on it.
		typename decltype(styles)::iterator sit = styles.begin();
		if (item.style && item.style->diffuse) {
			const auto& cc = *item.style->diffuse->components;
			auto c = std::make_pair(cc(0), std::make_pair(cc(1), cc(2)));
			auto it = diffuse_to_style.find(c);
			if (it == diffuse_to_style.end()) {
				styles.push_back(*item.style);
				sit = --styles.end();
				diffuse_to_style.insert({ c, sit });
			} else {
				sit = it->second;
			}
		}
		/*
		auto p = item.polyhedron;
		// Apply transformation
		for (auto &vertex : vertices(p)) {
			vertex->point() = vertex->point().transform(item.transformation);
		}
		*/

		TreeShapeType tree_polyhedron;
		poly_copy(tree_polyhedron, item.polyhedron);

		typename TreeKernel::Aff_transformation_3 transformation;
		transformation_copy(transformation, item.transformation);

		std::transform(
			tree_polyhedron.points_begin(), tree_polyhedron.points_end(), 
			tree_polyhedron.points_begin(), transformation);

		triangulated_shape_memory.push_back(tree_polyhedron);
		CGAL::Polygon_mesh_processing::triangulate_faces(triangulated_shape_memory.back());
		tree.insert(faces(triangulated_shape_memory.back()).first, faces(triangulated_shape_memory.back()).second, triangulated_shape_memory.back());

#ifdef GEOBIM_DEBUG
		obj_(nullptr, triangulated_shape_memory.back().facets_begin(), triangulated_shape_memory.back().facets_end());
#endif

		for (auto& f : faces(triangulated_shape_memory.back())) {
			typename TreeShapeType::Facet_handle F = f;
			facet_to_style.insert({ F, sit });
		}
	}

	void finalize() {
		tree.build();
		tree.accelerate_distance_queries();
	}

	typedef std::vector<std::pair<
		ifcopenshell::geometry::taxonomy::style*,
		std::list<cgal_shape_t::Facet_handle>>> segmentation_return_type;

	segmentation_return_type segment(const cgal_shape_t& input) {
		segmentation_return_type result(styles.size());
		auto it = styles.begin();
		for (size_t i = 0; i < styles.size(); ++i, ++it) {
			result[i].first = i ? &*it : nullptr;
		}

		CGAL::Cartesian_converter<Kernel_,TreeKernel> converter;

		for (auto &f : faces(input)) {
			auto O = CGAL::centroid(
				converter(f->facet_begin()->vertex()->point()),
				converter(f->facet_begin()->next()->vertex()->point()),
				converter(f->facet_begin()->next()->next()->vertex()->point())
			);
			auto pair = tree.closest_point_and_primitive(O);
			typename TreeShapeType::Face_handle F = pair.second.first;
			auto it = facet_to_style.find(F);
			if (it != facet_to_style.end()) {
				int sid = std::distance(styles.cbegin(), it->second);
				result[sid].second.push_back(f);
			}
		}

		return result;
	}
};

// State (polyhedra mostly) that are relevant only for one radius
struct capturing_execution_context : public execution_context {
	std::list<shape_callback_item> items;

	void operator()(shape_callback_item& item) {
		items.push_back(item);
	}
};

template <class HDS> 
class PolyFromMesh : public CGAL::Modifier_base<HDS> {
private:

	std::list<cgal_point_t> points_;
	std::vector<std::vector<int>> indices_;

public:
	PolyFromMesh(const std::list<cgal_point_t>& points, const std::vector<std::vector<int>>& indices)
		: points_(points)
		, indices_(indices) {}

	void operator()(HDS& hds) {
		CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, false);

		B.begin_surface(points_.size(), indices_.size());

		for (auto& p : points_) {
			B.add_vertex(p);
		}

		for (auto& fs : indices_) {
			B.begin_facet();
			for (auto& i : fs) {
				B.add_vertex_to_facet(i);
			}
			B.end_facet();
		}

		B.end_surface();
	}
};

#include <CGAL/exceptions.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>

// State (polyhedra mostly) that are relevant only for one radius
struct radius_execution_context : public execution_context {
	double radius;
	CGAL::Nef_nary_union_3< CGAL::Nef_polyhedron_3<Kernel_> > union_collector;
	CGAL::Nef_polyhedron_3<Kernel_> padding_cube, padding_cube_2, boolean_result, exterior, bounding_box, complement, complement_padded;
	cgal_shape_t polyhedron, polyhedron_exterior;
	enum extract_component { INTERIOR, EXTERIOR };
	bool minkowski_triangles_;

	radius_execution_context(double r, bool narrower = false, bool minkowski_triangles = false) : radius(r), minkowski_triangles_(minkowski_triangles) {
		{
			auto polycube = ifcopenshell::geometry::utils::create_cube(r);
			padding_cube = ifcopenshell::geometry::utils::create_nef_polyhedron(polycube);
		}

#ifdef ENSURE_2ND_OP_NARROWER
		if (narrower) {
			// double r2 = boost::math::float_advance(r, +5);
			double r2 = r + 1e-7;
			std::cout << r << " -> " << r2 << std::endl;
			auto polycube = ifcopenshell::geometry::utils::create_cube(r2);
			padding_cube_2 = ifcopenshell::geometry::utils::create_nef_polyhedron(polycube);
		} else // -> ...
#endif
		padding_cube_2 = padding_cube;
	}

	IfcUtil::IfcBaseEntity* previous_src = nullptr;
	std::string previous_geom_ref;
	CGAL::Nef_nary_union_3< CGAL::Nef_polyhedron_3<Kernel_> > per_product_collector;
	cgal_placement_t last_place;

	void operator()(shape_callback_item& item) {
		if (item.src != previous_src) {
			if (item.geom_reference == previous_geom_ref) {
				std::cout << "Reusing padded geometry for " << item.src->data().toString() << std::endl;
				auto product = per_product_collector.get_union();
				product.transform(last_place.inverse());
				product.transform(item.transformation);
				union_collector.add_polyhedron(product);
				return;
			}
			per_product_collector = CGAL::Nef_nary_union_3< CGAL::Nef_polyhedron_3<Kernel_> >();
		}
		
		CGAL::Polyhedron_3<CGAL::Epick> poly_triangulated;
		poly_copy(poly_triangulated, item.polyhedron);
		if (!CGAL::Polygon_mesh_processing::triangulate_faces(poly_triangulated)) {
			std::cerr << "unable to triangulate all faces" << std::endl;
			return;
		}

		std::vector<
			std::pair<
			boost::graph_traits<CGAL::Polyhedron_3<CGAL::Epick>>::face_descriptor,
			boost::graph_traits<CGAL::Polyhedron_3<CGAL::Epick>>::face_descriptor>> self_intersections;
		CGAL::Polygon_mesh_processing::self_intersections(poly_triangulated, std::back_inserter(self_intersections));

		previous_src = item.src;
		previous_geom_ref = item.geom_reference;
		last_place = item.transformation;
		
		CGAL::Nef_polyhedron_3<Kernel_> item_nef, result;
		bool item_nef_succeeded;
		if (!(item_nef_succeeded = item.to_nef_polyhedron(item_nef))) {
			std::cerr << "no nef for product" << std::endl;
		}

		bool result_set = false;
		
		if (!(minkowski_triangles_ || !item_nef_succeeded || !self_intersections.empty())) {
			item_nef.transform(item.transformation);

			auto T0 = timer.measure("minkowski_sum");
			try {
				result = CGAL::minkowski_sum_3(item_nef, padding_cube);
				result_set = true;
			} catch (CGAL::Failure_exception&) {
				std::cerr << "Minkowski on volume failed, retrying with individual triangles" << std::endl;
			}
			T0.stop();
		} 
		
		if (!result_set) {
			auto T2 = timer.measure("self_intersection_handling");

			if (self_intersections.size()) {
				std::cerr << self_intersections.size() << " self-intersections for product" << std::endl;
			}

			CGAL::Nef_nary_union_3< CGAL::Nef_polyhedron_3<Kernel_> > accum;

			for (auto &face : faces(poly_triangulated)) {

				if (!face->is_triangle()) {
					std::cout << "Warning: non-triangular face!" << std::endl;
					continue;
				}

				CGAL::Polyhedron_3<CGAL::Epick>::Halfedge_around_facet_const_circulator current_halfedge = face->facet_begin();
				CGAL::Point_3<CGAL::Epick> points[3];

				int i = 0;
				do {
					points[i] = current_halfedge->vertex()->point();
					++i;
					++current_halfedge;
				} while (current_halfedge != face->facet_begin());

				double A = std::sqrt(CGAL::to_double(CGAL::Triangle_3<CGAL::Epick>(points[0], points[1], points[2]).squared_area()));
				if (A < 1.e-5) {
					std::cout << "Skipping triangle with area " << A << std::endl;
					continue;
				}

				cgal_shape_t T;
				CGAL::Cartesian_converter<CGAL::Epick, CGAL::Epeck> C;
				T.make_triangle(C(points[0]), C(points[1]), C(points[2]));

				CGAL::Nef_polyhedron_3<Kernel_> Tnef(T);

				CGAL::Nef_polyhedron_3<Kernel_> padded = CGAL::minkowski_sum_3(Tnef, padding_cube);
				accum.add_polyhedron(padded);
			}

			result = accum.get_union();
			result.transform(item.transformation);

			T2.stop();
		}

		auto T1 = timer.measure("opening_handling");
		if (item.wall_direction && item.openings.size()) {
			static const Eigen::Vector3d Zax(0, 0, 1);
			// @todo derive from model.
			// @todo since many walls will be parallel we can cache these polyhedrons
			static const double EPS = 1.e-5;

			auto Yax = Zax.cross(*item.wall_direction).normalized();

			auto x0 = *item.wall_direction * -radius;
			auto x1 = *item.wall_direction * +radius;

			auto y0 = Yax * -(radius + EPS);
			auto y1 = Yax * +(radius + EPS);

			Kernel_::Point_3 X0(x0(0), x0(1), x0(2));
			Kernel_::Point_3 X1(x1(0), x1(1), x1(2));

			Kernel_::Point_3 Y0(y0(0), y0(1), y0(2));
			Kernel_::Point_3 Y1(y1(0), y1(1), y1(2));

			Kernel_::Point_3 Z0(0, 0, -radius);
			Kernel_::Point_3 Z1(0, 0, +radius);

			// CGAL::Nef_polyhedron_3<Kernel_> X(CGAL::Segment_3<Kernel_>(X0, X1));
			CGAL::Nef_polyhedron_3<Kernel_> Y(CGAL::Segment_3<Kernel_>(Y0, Y1));
			// CGAL::Nef_polyhedron_3<Kernel_> Z(CGAL::Segment_3<Kernel_>(Z0, Z1));
			// auto ZX = CGAL::minkowski_sum_3(X, Z);

			// manual minkowski sum...
			CGAL::Polyhedron_3<Kernel_> zx;
			std::list<cgal_point_t> zx_points{ {
				CGAL::ORIGIN + ((X0 - CGAL::ORIGIN) + (Z0 - CGAL::ORIGIN)),
				CGAL::ORIGIN + ((X1 - CGAL::ORIGIN) + (Z0 - CGAL::ORIGIN)),
				CGAL::ORIGIN + ((X1 - CGAL::ORIGIN) + (Z1 - CGAL::ORIGIN)),
				CGAL::ORIGIN + ((X0 - CGAL::ORIGIN) + (Z1 - CGAL::ORIGIN))
			} };
			std::vector<std::vector<int>> zx_idxs{ {{{0,1,2,3}}} };
			PolyFromMesh<cgal_shape_t::HDS> m(zx_points, zx_idxs);
			zx.delegate(m);
			CGAL::Nef_polyhedron_3<Kernel_> ZX(zx);

			CGAL::Nef_nary_union_3< CGAL::Nef_polyhedron_3<Kernel_> > opening_union;
			for (auto& op : item.openings) {
				auto bounds = create_bounding_box(op->polyhedron);				
				CGAL::Nef_polyhedron_3<Kernel_> opening_nef;
				if (!op->to_nef_polyhedron(opening_nef)) {
					std::cerr << "no nef for opening" << std::endl;
					continue;
				}
				opening_nef.transform(op->transformation);
				bounds.transform(op->transformation);

				auto temp = bounds - opening_nef;
				temp = CGAL::minkowski_sum_3(ZX, temp);
				temp = bounds - temp;
				temp = CGAL::minkowski_sum_3(temp, Y);

#ifdef GEOBIM_DEBUG
				simple_obj_writer x("opening-" + op->id);
				x(nullptr, item.polyhedron.facets_begin(), item.polyhedron.facets_end());
				x(nullptr, op->polyhedron.facets_begin(), op->polyhedron.facets_end());
				CGAL::Polyhedron_3<Kernel_> temp2;
				temp.convert_to_polyhedron(temp2);
				x(nullptr, temp2.facets_begin(), temp2.facets_end());
#endif

				result -= temp;
			}
		}
		T1.stop();
		
		union_collector.add_polyhedron(result);
		per_product_collector.add_polyhedron(result);
	}

	// Extract the exterior component of a CGAL Polyhedron
	cgal_shape_t extract(const cgal_shape_t& input, extract_component component) const {
		if (input.facets_begin() == input.facets_end()) {
			throw std::runtime_error("Empty input operand to extract()");
		}

		// Input is going to be mutated, so make a copy first
		cgal_shape_t input_copy = input;

		Kernel_::FT max_x = -1e9;
		cgal_shape_t::Vertex_handle left_most_vertex;
		for (auto it = input_copy.edges_begin(); it != input_copy.edges_end(); ++it) {
			auto vx = it->vertex()->point().cartesian(0);
			if (vx > max_x) {
				max_x = vx;
				left_most_vertex = it->vertex();
			}
		}

		
		if (component == EXTERIOR) {
			std::set<cgal_shape_t::Facet_handle> empty;
			auto connected = connected_faces(left_most_vertex->halfedge()->facet(), empty);

			std::set<cgal_shape_t::Facet_handle> outer_facets(connected.begin(), connected.end());

			while (std::distance(input_copy.facets_begin(), input_copy.facets_end()) > outer_facets.size()) {
				for (auto it = input_copy.facets_begin(); it != input_copy.facets_end(); ++it) {
					if (outer_facets.find(&*it) == outer_facets.end()) {
						input_copy.erase_connected_component(it->halfedge());
						break;
					}
				}
			}
		} else {
			input_copy.erase_connected_component(left_most_vertex->halfedge());
		}

		return input_copy;
	}

	// Create a bounding box (six-faced Nef poly) around a CGAL Polyhedron
	CGAL::Nef_polyhedron_3<Kernel_> create_bounding_box(const cgal_shape_t& input) const {
		// @todo we can probably use
		// CGAL::Nef_polyhedron_3<Kernel_> nef(CGAL::Nef_polyhedron_3<Kernel_>::COMPLETE)
		// Implementation detail: there is always an implicit box around the Nef, even when open or closed.
		// Never mind: The Minkowski sum cannot operate on unbounded inputs...

		// Create the complement of the Nef by subtracting from its bounding box,
		// see: https://github.com/tudelft3d/ifc2citygml/blob/master/off2citygml/Minkowski.cpp#L23
		auto bounding_box = CGAL::Polygon_mesh_processing::bbox(input);
		Kernel_::Point_3 bbmin(bounding_box.xmin(), bounding_box.ymin(), bounding_box.zmin());
		Kernel_::Point_3 bbmax(bounding_box.xmax(), bounding_box.ymax(), bounding_box.zmax());
		Kernel_::Vector_3 d(radius, radius, radius);
		bbmin = CGAL::ORIGIN + ((bbmin - CGAL::ORIGIN) - d);
		bbmax = CGAL::ORIGIN + ((bbmax - CGAL::ORIGIN) + d);
		cgal_shape_t poly_box = ifcopenshell::geometry::utils::create_cube(bbmin, bbmax);
		return ifcopenshell::geometry::utils::create_nef_polyhedron(poly_box);
	}

	// Completes the boolean union, extracts exterior and erodes padding radius
	void finalize() {
		{
			auto T = timer.measure("nef_boolean_union");
			// @todo spatial sorting?
			boolean_result = union_collector.get_union();
			T.stop();
		}

		auto T2 = timer.measure("result_nef_processing");
		polyhedron = ifcopenshell::geometry::utils::create_polyhedron(boolean_result);

		{
			simple_obj_writer tmp_debug("debug-after-boolean");
			tmp_debug(nullptr, polyhedron.facets_begin(), polyhedron.facets_end());
		}

		// @todo Wasteful: remove interior on Nef?
		polyhedron_exterior = extract(polyhedron, EXTERIOR);

		{
			simple_obj_writer tmp_debug("debug-exterior");
			tmp_debug(nullptr, polyhedron_exterior.facets_begin(), polyhedron_exterior.facets_end());
		}

		exterior = ifcopenshell::geometry::utils::create_nef_polyhedron(polyhedron_exterior);
		bounding_box = create_bounding_box(polyhedron);

		complement = bounding_box - exterior;
		complement.extract_regularization();
		// @nb padding cube is potentially slightly larger to result in a thinner result
		// then another radius for comparison.
		complement_padded = CGAL::minkowski_sum_3(complement, padding_cube_2);
		complement_padded.extract_regularization();
		T2.stop();

		{
			auto T = timer.measure("result_nef_to_poly");
#if 0
			// @todo I imagine this operation is costly, we can also convert the padded complement to
			// polyhedron, and remove the connected component that belongs to the bbox, then reverse
			// the remaining poly to point to the interior?
			exterior -= complement_padded;
#else
			// Rougly twice as fast as the complexity is half (box complexity is negligable).

			// Re above: extracting the interior shell did not prove to be reliable even with
			// the undocumented function convert_inner_shell_to_polyhedron(). Therefore we
			// subtract from the padded box as that will have lower complexity than above.
			// Mark_bounded_volumes on the completement also did not work.
			exterior = bounding_box - complement_padded;
#endif
			T.stop();
		}

		exterior.extract_regularization();

		if (exterior.is_simple()) {
			auto T1 = timer.measure("result_nef_to_poly");
			polyhedron_exterior = ifcopenshell::geometry::utils::create_polyhedron(exterior);
			T1.stop();

			auto vol = CGAL::Polygon_mesh_processing::volume(polyhedron_exterior);
			std::cout << "Volume with radius " << radius << " is " << vol << std::endl;
		} else {
			CGAL::convert_nef_polyhedron_to_polygon_mesh(exterior, polyhedron_exterior);
			std::cout << "Result with radius " << radius << " is not manifold" << std::endl;
		}
	}
};

#ifdef ENSURE_2ND_OP_NARROWER
#define MAKE_OP2_NARROWER -2e-7
#else
#define MAKE_OP2_NARROWER
#endif

struct radius_comparison {
	struct hollow_solid {
		typedef CGAL::Nef_polyhedron_3<Kernel_> nef;

		nef bbox, complement, complement_padded, inner, hollow, cube;

		double D;

		hollow_solid(radius_execution_context& a, double d) {
			D = d;
			bbox = a.create_bounding_box(a.polyhedron_exterior);
			complement = bbox - a.exterior;
			complement.extract_regularization();
			auto polycube = ifcopenshell::geometry::utils::create_cube(d);
			cube = ifcopenshell::geometry::utils::create_nef_polyhedron(polycube);
			complement_padded = CGAL::minkowski_sum_3(complement, cube);
			complement_padded.extract_regularization();
			inner = bbox - complement_padded;
			inner.extract_regularization();
			{
				auto inner_poly = ifcopenshell::geometry::utils::create_polyhedron(inner);
				simple_obj_writer obj("debug-inner-" + boost::lexical_cast<std::string>(D));
				obj(nullptr, inner_poly.facets_begin(), inner_poly.facets_end());
			}
			hollow = a.exterior - inner;
			hollow.extract_regularization();
		}
	};

	typedef CGAL::Nef_polyhedron_3<Kernel_> nef;
	typedef CGAL::Polyhedron_3<Kernel_> poly;

	nef difference_nef;
	poly difference_poly;
	hollow_solid A, B;

	radius_comparison(radius_execution_context& a, radius_execution_context& b, double d)
//		: A(a, d), B(b, boost::math::float_advance(d, -10))
		: A(a, d), B(b, d MAKE_OP2_NARROWER)
	{
		difference_nef = B.hollow - A.hollow;
		difference_nef.extract_regularization();
		difference_poly = ifcopenshell::geometry::utils::create_polyhedron(difference_nef);
	}
};

// Parse the command line settings and do basic initialization
int parse_command_line(geobim_settings& settings, int argc, char** argv) {
	std::string entities;

	typedef po::command_line_parser command_line_parser;
	po::options_description options("Command line options");
	options.add_options()
		("help,h", "display usage information")
		("version,v", "display version information")
		("debug,d", "more verbose log messages")
		("openings,o", "whether to process opening subtractions")
		("timings,t", "print timings after execution")
		//("exact-segmentation,e", "use exact kernel in proximity tree for semantic association (not recommended)")
		("openings-posthoc,O", "whether to process opening subtractions posthoc")
		("minkowski-triangles,T", "force minkowski sum on individual triangles, slow..")
		("entities,e", new po::typed_value<std::string, char>(&entities), "semicolon separated list of IFC entities to include or exclude")
		("exclude,x", "entities are to be excluded")
		("input-file", new po::typed_value<std::string, char>(&settings.input_filename), "input IFC file")
		("output-file", new po::typed_value<std::string, char>(&settings.output_filename), "output OBJ file")
		("radii", boost::program_options::value<std::vector<std::string>>()->multitoken());

	po::positional_options_description positional_options;
	positional_options.add("input-file", 1);
	positional_options.add("output-file", 1);
	positional_options.add("radii", -1);

	po::variables_map vmap;
	try {
		po::store(command_line_parser(argc, argv).
			options(options).positional(positional_options).run(), vmap);
	} catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
	}

	po::notify(vmap);

	if (vmap.count("version")) {
		return 0;
	} else if (vmap.count("help")) {
		return 0;
	} else if (!vmap.count("input-file")) {
		return 1;
	}

	if (!vmap.count("output-file")) {
		settings.output_filename = settings.input_filename;
	}

	// using Kernel_::FT creates weird segfaults, probably due to how ifopsh constructs the box, coordinates components cannot be shared?
	if (vmap.count("radii")) {
		const auto& radii_str = vmap["radii"].as<std::vector<std::string>>();
		std::transform(radii_str.begin(), radii_str.end(), std::back_inserter(settings.radii), [](const std::string& s) {
			// radii.push_back(CGAL::Gmpq(argv[i]));
			return boost::lexical_cast<double>(s);
		});
		std::sort(settings.radii.begin(), settings.radii.end());
	}

	if (settings.radii.empty()) {
		// Binary search will be used.
	}

	settings.apply_openings = vmap.count("openings");
	settings.apply_openings_posthoc = vmap.count("openings-posthoc");
	settings.exact_segmentation = vmap.count("exact-segmentation");
	settings.debug = vmap.count("debug");
	settings.minkowski_triangles = vmap.count("minkowski-triangles");

	settings.file = new IfcParse::IfcFile(settings.input_filename);

	if (!settings.file->good()) {
		std::cerr << "[Error] Unable to parse input file '" << settings.input_filename << "'";
		return 1;
	}

	settings.settings.set(ifcopenshell::geometry::settings::USE_WORLD_COORDS, false);
	settings.settings.set(ifcopenshell::geometry::settings::WELD_VERTICES, false);
	settings.settings.set(ifcopenshell::geometry::settings::SEW_SHELLS, true);
	settings.settings.set(ifcopenshell::geometry::settings::CONVERT_BACK_UNITS, true);
	settings.settings.set(ifcopenshell::geometry::settings::DISABLE_TRIANGULATION, true);
	settings.settings.set(ifcopenshell::geometry::settings::DISABLE_OPENING_SUBTRACTIONS, !settings.apply_openings);

	if (vmap.count("entities")) {
		std::vector<std::string> tokens;
		boost::split(tokens, vmap["entities"].as<std::string>(), boost::is_any_of(";"));
		if (!tokens.empty()) {
			settings.entity_names.emplace();
			settings.entity_names->clear();
			settings.entity_names->insert(tokens.begin(), tokens.end());
		}
		settings.entity_names_included = !vmap.count("exclude");
	}

	return 0;
}

// Interprets IFC geometries by means of IfcOpenShell CGAL and
// pass result to callback
template <typename Fn>
int process_geometries(geobim_settings& settings, Fn& fn) {

	// @todo static now to prevent being gc'ed
	static opening_collector all_openings(settings.file);

	// Capture all openings beforehand, they are later assigned to the
	// building elements.
	auto opening_settings = settings;
	opening_settings.entity_names.emplace();
	opening_settings.entity_names->insert("IfcOpeningElement");
	opening_settings.entity_names_included = true;
	if (settings.apply_openings_posthoc && settings.entity_names != opening_settings.entity_names) {
		process_geometries(opening_settings, all_openings);
	}

	std::vector<ifcopenshell::geometry::filter_t> filters;
	if (settings.entity_names) {
		if (!settings.entity_names_included) {
			settings.entity_names->insert("IfcSpace");
			settings.entity_names->insert("IfcOpeningElement");
		}
		filters.push_back(IfcGeom::entity_filter(settings.entity_names_included, false, *settings.entity_names));
	} else {
		filters.push_back(IfcGeom::entity_filter(false, false, {"IfcSpace", "IfcOpeningElement"}));
	}
	
	ifcopenshell::geometry::Iterator context_iterator("cgal", settings.settings, settings.file, filters);

	auto T = timer.measure("ifc_geometry_processing");
	if (!context_iterator.initialize()) {
		return 1;
	}
	T.stop();

	size_t num_created = 0;

	auto axis_settings = settings.settings;
	axis_settings.set(ifcopenshell::geometry::settings::EXCLUDE_SOLIDS_AND_SURFACES, true);
	axis_settings.set(ifcopenshell::geometry::settings::INCLUDE_CURVES, true);
	auto geometry_mapper = ifcopenshell::geometry::impl::mapping_implementations().construct(settings.file, axis_settings);

	for (;; ++num_created) {
		bool has_more = true;
		if (num_created) {
			auto T0 = timer.measure("ifc_geometry_processing");
			has_more = context_iterator.next();
			T0.stop();
		}
		ifcopenshell::geometry::NativeElement* geom_object = nullptr;
		if (has_more) {
			geom_object = context_iterator.get_native();
		}
		if (!geom_object) {
			break;
		}
		
		if (geom_object->guid() == "3es57B9Kr3nxL4uBITV$0e") {
			std::cout << "NOTICE Skipping: " << geom_object->product()->data().toString() << std::endl;
			continue;
		}

		std::cout << "Processing: " << geom_object->product()->data().toString() << std::endl;
		
		const auto& n = *geom_object->transformation().data().components;
		const cgal_placement_t element_transformation (
			n(0, 0), n(0, 1), n(0, 2), n(0, 3),
			n(1, 0), n(1, 1), n(1, 2), n(1, 3),
			n(2, 0), n(2, 1), n(2, 2), n(2, 3));

		boost::optional<Eigen::Vector3d> wall_direction;


		std::list<shape_callback_item*> openings;

		auto p = all_openings.map.equal_range(geom_object->product());
		for (auto it = p.first; it != p.second; ++it) {
			openings.push_back(it->second);
		}

		if (settings.apply_openings_posthoc && geom_object->product()->declaration().is("IfcWall")) {
			auto T2 = timer.measure("wall_axis_handling");
			auto item = geometry_mapper->map(geom_object->product());
			if (item) {
				typedef ifcopenshell::geometry::taxonomy::collection cl;
				typedef ifcopenshell::geometry::taxonomy::loop l;
				typedef ifcopenshell::geometry::taxonomy::edge e;
				typedef ifcopenshell::geometry::taxonomy::point3 p;
				auto edge = (e*)((l*)((cl*)((cl*)item)->children[0])->children[0])->children[0];
				if (edge->basis == nullptr && edge->start.which() == 0 && edge->end.which() == 0) {
					const auto& p0 = boost::get<p>(edge->start);
					const auto& p1 = boost::get<p>(edge->end);
					Eigen::Vector4d P0 = p0.components->homogeneous();
					Eigen::Vector4d P1 = p1.components->homogeneous();
					auto V0 = n * P0;
					auto V1 = n * P1;
					std::cout << "Axis " << V0(0) << " " << V0(1) << " " << V0(2) << " -> "
						<< V1(0) << " " << V1(1) << " " << V1(2) << std::endl;
					wall_direction = (V1 - V0).head<3>().normalized();
				}
				T2.stop();
			}
		}

		for (auto& g : geom_object->geometry()) {
			auto s = ((ifcopenshell::geometry::CgalShape*) g.Shape())->shape();
			const auto& m = *g.Placement().components;
			
			const cgal_placement_t part_transformation (
				m(0, 0), m(0, 1), m(0, 2), m(0, 3),
				m(1, 0), m(1, 1), m(1, 2), m(1, 3),
				m(2, 0), m(2, 1), m(2, 2), m(2, 3));

			// Apply transformation
			for (auto &vertex : vertices(s)) {
				vertex->point() = vertex->point().transform(part_transformation);
			}

			boost::optional<ifcopenshell::geometry::taxonomy::style> opt_style;
			if (g.hasStyle()) {
				opt_style = g.Style();
			}

			shape_callback_item item {
				geom_object->product(),
				geom_object->guid(),
				geom_object->type(),
				geom_object->geometry().id(),
				element_transformation,
				s,
				opt_style,
				wall_direction,
				openings
			};
				
			fn(item);

			std::cout << "Processed: " << geom_object->product()->data().toString() << " part: #" << geom_object->geometry().id() << std::endl;
		}
		
		std::cout << "Progress: " << context_iterator.progress() << std::endl;

	}

	return 0;
}

// A structure for recieving processed shapes simply defers to a vector of contexts
struct shape_callback {
	std::vector<execution_context*> contexts;

	void operator()(shape_callback_item& item) {
		for (auto& c : contexts) {
			(*c)(item);
		}
	}
};

template <typename It>
double initialize_radius_context_and_get_volume_with_cache(It first, It second, double radius) {
	static std::map<double, double> cache_;

	auto it = cache_.find(radius);
	if (it != cache_.end()) {
		std::cout << "Used cache for R=" << radius << " V=" << it->second << std::endl;
		return it->second;
	}

	radius_execution_context rec(radius, false);
	// Unfortunately for_each() is by value so needs to be wrapped in a lambda with side effects
	std::for_each(first, second, [&rec](auto& v) {
		rec(v);
	});
	rec.finalize();
	double V = CGAL::to_double(CGAL::Polygon_mesh_processing::volume(rec.polyhedron_exterior));

	std::cout << "Calculated for R=" << radius << " V=" << V << std::endl;

	cache_.insert(it, { radius , V });
	return V;
}

template <typename It>
double binary_search(It first, It second, std::pair<double, double> range) {
	double abc[3];
	std::tie(abc[0], abc[2]) = range;

	std::cout << "Testing " << abc[0] << " and " << abc[2] << std::endl;

	auto a_vol = initialize_radius_context_and_get_volume_with_cache(first, second, abc[0]);
	auto b_vol = initialize_radius_context_and_get_volume_with_cache(first, second, abc[2]);

	if (a_vol * 1.1 < b_vol) {
		if ((abc[2] - abc[0]) < 1.e-4) {
			std::cout << "Terminating search at " << abc[0] << " and " << abc[2] << std::endl;
			return (abc[0] + abc[2]) / 2.;
		}

		abc[1] = (abc[0] + abc[2]) / 2.;
		for (int i = 1; i >= 0; --i) {
			auto r = binary_search(first, second, { abc[i], abc[i+1] });
			if (r != abc[i + 1]) {
				return r;
			}
		}
	}

	return abc[2];
}

int main(int argc, char** argv) {
	Logger::SetOutput(&std::cerr, &std::cerr);
	Logger::Verbosity(Logger::LOG_NOTICE);

	geobim_settings settings;
	parse_command_line(settings, argc, argv);

	global_execution_context<CGAL::Epick> global_context;
	global_execution_context<Kernel_> global_context_exact;

	shape_callback callback;
	if (settings.exact_segmentation) {
		callback.contexts.push_back(&global_context_exact);
	} else {
		callback.contexts.push_back(&global_context);
	}

#ifdef GEOBIM_DEBUG
	callback.contexts.push_back(new debug_writer);
#endif

	if (settings.radii.empty()) {
		auto cec = new capturing_execution_context;
		callback.contexts.push_back(cec);
		process_geometries(settings, callback);
		auto R = binary_search(cec->items.begin(), cec->items.end(), { 1.e-3, 0.2 });
		std::cout << "Largest gap found with R ~ " << (R * 2.) << std::endl;
	} else {
		std::vector<radius_execution_context> radius_contexts;
		bool first = true;
		for (double r : settings.radii) {
			// 2nd is narrower (depending on ifdef above, appears to be necessary).
			radius_contexts.emplace_back(r, !first, settings.minkowski_triangles);
			first = false;
		}

		for (auto& c : radius_contexts) {
			callback.contexts.push_back(&c);
		}

		process_geometries(settings, callback);

		std::cout << "done processing geometries" << std::endl;
		delete settings.file;

		auto T1 = timer.measure("semantic_segmentation");
		if (settings.exact_segmentation) {
			global_context_exact.finalize();
		} else {
			global_context.finalize();
		}
		T1.stop();

		for (auto& c : radius_contexts) {
			c.finalize();

			auto T0 = timer.measure("semantic_segmentation");
			global_execution_context<Kernel_>::segmentation_return_type style_facet_pairs;
			if (settings.exact_segmentation) {
				style_facet_pairs = global_context_exact.segment(c.polyhedron_exterior);
			} else {
				style_facet_pairs = global_context.segment(c.polyhedron_exterior);
			}
			T0.stop();

			city_json_writer write_city(settings.output_filename + boost::lexical_cast<std::string>(c.radius));
			for (auto& p : style_facet_pairs) {
				write_city(p.first, p.second.begin(), p.second.end());
			}

			simple_obj_writer write_obj(settings.output_filename + boost::lexical_cast<std::string>(c.radius));
			for (auto& p : style_facet_pairs) {
				write_obj(p.first, p.second.begin(), p.second.end());
			}
		}

		auto T2 = timer.measure("difference_overlay");
		auto it = radius_contexts.begin();
		for (auto jt = it + 1; jt != radius_contexts.end(); ++it, ++jt) {
			radius_comparison difference(*it, *jt, 0.001);
			simple_obj_writer obj("difference-"
				+ boost::lexical_cast<std::string>(it->radius) + "-"
				+ boost::lexical_cast<std::string>(jt->radius));
			obj(nullptr, difference.difference_poly.facets_begin(), difference.difference_poly.facets_end());
		}
		T2.stop();

	}

	timer.print(std::cout);
}
