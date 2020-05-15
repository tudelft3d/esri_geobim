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

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>

#include <chrono>

// Can be used to convert polyhedron from exact to inexact and vice-versa
template <class Polyhedron_input,
	class Polyhedron_output>
	struct Copy_polyhedron_to
	: public CGAL::Modifier_base<typename Polyhedron_output::HalfedgeDS> {
	Copy_polyhedron_to(const Polyhedron_input& in_poly)
		: in_poly(in_poly) {}

	void operator()(typename Polyhedron_output::HalfedgeDS& out_hds) {
		typedef typename Polyhedron_output::HalfedgeDS Output_HDS;
		typedef typename Polyhedron_input::HalfedgeDS Input_HDS;

		CGAL::Polyhedron_incremental_builder_3<Output_HDS> builder(out_hds);

		typedef typename Polyhedron_input::Vertex_const_iterator Vertex_const_iterator;
		typedef typename Polyhedron_input::Facet_const_iterator  Facet_const_iterator;
		typedef typename Polyhedron_input::Halfedge_around_facet_const_circulator HFCC;

		CGAL::Cartesian_converter<
			Polyhedron_input::Traits::Kernel,
			Polyhedron_output::Traits::Kernel> converter;

		builder.begin_surface(in_poly.size_of_vertices(),
			in_poly.size_of_facets(),
			in_poly.size_of_halfedges());

		for (Vertex_const_iterator
			vi = in_poly.vertices_begin(), end = in_poly.vertices_end();
			vi != end; ++vi) {
			builder.add_vertex(converter(vi->point()));
		}

		typedef CGAL::Inverse_index<Vertex_const_iterator> Index;
		Index index(in_poly.vertices_begin(), in_poly.vertices_end());

		for (Facet_const_iterator
			fi = in_poly.facets_begin(), end = in_poly.facets_end();
			fi != end; ++fi) {
			HFCC hc = fi->facet_begin();
			HFCC hc_end = hc;
			builder.begin_facet();
			do {
				builder.add_vertex_to_facet(index[hc->vertex()]);
				++hc;
			} while (hc != hc_end);
			builder.end_facet();
		}
		builder.end_surface();
	} // end operator()(..)
	private:
		const Polyhedron_input& in_poly;
}; // end Copy_polyhedron_to<>

template <class Poly_B, class Poly_A>
typename std::enable_if<std::is_same<Poly_A, Poly_B>::value>::type poly_copy(Poly_B& poly_b, const Poly_A& poly_a) {
	poly_b = poly_a;
}

template <class Poly_B, class Poly_A>
typename std::enable_if<!std::is_same<Poly_A, Poly_B>::value>::type poly_copy(Poly_B& poly_b, const Poly_A& poly_a) {
	poly_b.clear();
	Copy_polyhedron_to<Poly_A, Poly_B> modifier(poly_a);
	poly_b.delegate(modifier);
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

// OBJ writer for CGAL facets paired with a style
struct simple_obj_writer {
	int group_id = 1;
	int vertex_count = 1;
	std::ofstream obj, mtl;
	ifcopenshell::geometry::taxonomy::colour BLACK;

	simple_obj_writer(const std::string& fn_prefix)
		: obj((fn_prefix + ".obj").c_str())
		, mtl((fn_prefix + ".mtl").c_str()) {
		obj << "mtllib " << fn_prefix << ".mtl\n";
		BLACK.components = Eigen::Vector3d(0., 0., 0.);
	}

	std::array<Kernel_::Point_3, 3> points_from_facet(cgal_shape_t::Facet_handle f) {
		return {
				f->facet_begin()->vertex()->point(),
				f->facet_begin()->next()->vertex()->point(),
				f->facet_begin()->next()->next()->vertex()->point()
		};
	}

	std::array<Kernel_::Point_3, 3> points_from_facet(std::list<cgal_shape_t::Facet_handle>::iterator f) {
		return points_from_facet(*f);
	}
	
	template <typename It>
	void operator()(const ifcopenshell::geometry::taxonomy::style* style, It begin, It end) {
		auto diffuse = style ? style->diffuse.get_value_or(BLACK).components : BLACK.components;

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

// Global settings that are derived from the command line invocation
// parsed by Boost.ProgramOptions
struct geobim_settings {
	std::string input_filename, output_filename;
	std::vector<double> radii;
	bool apply_openings, apply_openings_posthoc, debug, exact_segmentation;
	ifcopenshell::geometry::settings settings;
	std::set<std::string> entity_names;
	IfcParse::IfcFile* file;
};

// Generated for every representation *item* in the IFC file
struct shape_callback_item {
	IfcUtil::IfcBaseEntity* src;
	std::string id, type;
	cgal_shape_t polyhedron;
	CGAL::Nef_polyhedron_3<Kernel_> nef_polyhedron;
	boost::optional<ifcopenshell::geometry::taxonomy::style> style;
	boost::optional<Eigen::Vector3d> wall_direction;
	std::list<shape_callback_item*> openings;
};

// Prototype of a context to which processed shapes will be fed
struct execution_context {
	virtual void operator()(shape_callback_item&) = 0;
};

struct opening_collector : public execution_context {
	std::list<shape_callback_item> list;
	std::map<IfcUtil::IfcBaseEntity*, IfcUtil::IfcBaseEntity*> opening_to_elem;
	std::multimap<IfcUtil::IfcBaseEntity*, shape_callback_item*> map;

	opening_collector(IfcParse::IfcFile* f) {
		auto rels = f->instances_by_type("IfcRelVoidsElement");
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
	std::map<std::pair<double, std::pair<double, double>>, decltype(styles)::iterator> diffuse_to_style;

	// A reference is kept to the original shapes in a std::list.
	// Later an aabb tree is used map eroded triangle centroids
	// back to the original elements to preserve semantics.
	std::list<TreeShapeType> triangulated_shape_memory;
	std::map<typename TreeShapeType::Facet_handle, decltype(styles)::const_iterator> facet_to_style;

	global_execution_context() {
		// style 0 is for elements without style annotations
		styles.emplace_back();
	}

	void operator()(shape_callback_item& item) {
		size_t style_idx = item.style ? styles.size() : 0;

		// Group taxonomy::styles based on diffuse colour since we
		// do not have an equality operator on it.
		decltype(styles)::iterator sit = styles.begin();
		if (item.style && item.style->diffuse) {
			auto cc = item.style->diffuse->components;
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

		TreeShapeType tree_polyhedron;
		poly_copy(tree_polyhedron, item.polyhedron);

		triangulated_shape_memory.push_back(tree_polyhedron);
		CGAL::Polygon_mesh_processing::triangulate_faces(triangulated_shape_memory.back());
		tree.insert(faces(triangulated_shape_memory.back()).first, faces(triangulated_shape_memory.back()).second, triangulated_shape_memory.back());
		for (auto& f : faces(triangulated_shape_memory.back())) {
			TreeShapeType::Facet_handle F = f;
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
			auto it = facet_to_style.find(pair.second);
			if (it != facet_to_style.end()) {
				int sid = std::distance(styles.cbegin(), it->second);
				result[sid].second.push_back(f);
			}
		}

		return result;
	}
};

// State (polyhedra mostly) that are relevant only for one radius
struct radius_execution_context : public execution_context {
	double radius;
	CGAL::Nef_nary_union_3< CGAL::Nef_polyhedron_3<Kernel_> > union_collector;
	CGAL::Nef_polyhedron_3<Kernel_> padding_cube, boolean_result, exterior, bounding_box, complement, complement_padded;
	cgal_shape_t polyhedron, polyhedron_exterior;
	enum extract_component { INTERIOR, EXTERIOR };

	radius_execution_context(double r) : radius(r) {
		auto polycube = ifcopenshell::geometry::utils::create_cube(r);
		padding_cube = ifcopenshell::geometry::utils::create_nef_polyhedron(polycube);
	}

	void operator()(shape_callback_item& item) {
		CGAL::Nef_polyhedron_3<Kernel_> result;

		auto T0 = timer.measure("minkowski_sum");
		result = CGAL::minkowski_sum_3(item.nef_polyhedron, padding_cube);
		T0.stop();

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

			CGAL::Nef_polyhedron_3<Kernel_> X(CGAL::Segment_3<Kernel_>(X0, X1));
			CGAL::Nef_polyhedron_3<Kernel_> Y(CGAL::Segment_3<Kernel_>(Y0, Y1));
			CGAL::Nef_polyhedron_3<Kernel_> Z(CGAL::Segment_3<Kernel_>(Z0, Z1));
			auto ZX = CGAL::minkowski_sum_3(X, Z);

			CGAL::Nef_nary_union_3< CGAL::Nef_polyhedron_3<Kernel_> > opening_union;
			for (auto& op : item.openings) {
				auto bounds = create_bounding_box(op->polyhedron);
				auto temp = bounds - op->nef_polyhedron;
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
	}

	// Extract the exterior component of a CGAL Polyhedron
	cgal_shape_t extract(const cgal_shape_t& input, extract_component component) const {
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
			boolean_result = union_collector.get_union();
			T.stop();
		}

		auto T2 = timer.measure("result_nef_processing");
		polyhedron = ifcopenshell::geometry::utils::create_polyhedron(boolean_result);
		polyhedron_exterior = extract(polyhedron, EXTERIOR);
		exterior = ifcopenshell::geometry::utils::create_nef_polyhedron(polyhedron_exterior);
		bounding_box = create_bounding_box(polyhedron);
		complement = bounding_box - exterior;
		complement.extract_regularization();
		complement_padded = CGAL::minkowski_sum_3(complement, padding_cube);
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

		auto T1 = timer.measure("result_nef_to_poly");
		polyhedron_exterior = ifcopenshell::geometry::utils::create_polyhedron(exterior);
		T1.stop();

		auto vol = CGAL::Polygon_mesh_processing::volume(polyhedron_exterior);
		std::cout << "Volume with radius " << radius << " is " << vol << std::endl;
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
		("exact-segmentation,e", "use exact kernel in proximity tree for semantic association (not recommended)")
		("openings-posthoc,O", "whether to process opening subtractions posthoc")
		("entities,e", new po::typed_value<std::string, char>(&entities), "semicolon separated list of IFC entities to include")
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
	const auto& radii_str = vmap["radii"].as<std::vector<std::string>>();
	std::transform(radii_str.begin(), radii_str.end(), std::back_inserter(settings.radii), [](const std::string& s) {
		// radii.push_back(CGAL::Gmpq(argv[i]));
		return boost::lexical_cast<double>(s);
	});

	if (settings.radii.empty()) {
		settings.radii.push_back(0.01);
	}

	settings.apply_openings = vmap.count("openings");
	settings.apply_openings_posthoc = vmap.count("openings-posthoc");
	settings.exact_segmentation = vmap.count("exact-segmentation");
	settings.debug = vmap.count("debug");

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

	settings.entity_names = { "IfcWall", "IfcSlab" };
	if (vmap.count("entities")) {
		std::vector<std::string> tokens;
		boost::split(tokens, vmap["entities"].as<std::string>(), boost::is_any_of(";"));
		if (!tokens.empty()) {
			settings.entity_names.clear();
			settings.entity_names.insert(tokens.begin(), tokens.end());
		}
	}

	return 0;
}

// Interprets IFC geometries by means of IfcOpenShell CGAL and
// pass result to callback
template <typename Fn>
int process_geometries(geobim_settings& settings, Fn& fn) {

	opening_collector all_openings(settings.file);

	// Capture all openings beforehand, they are later assigned to the
	// building elements.
	auto opening_settings = settings;
	opening_settings.entity_names = { "IfcOpeningElement" };
	if (settings.apply_openings_posthoc && settings.entity_names != opening_settings.entity_names) {
		process_geometries(opening_settings, all_openings);
	}

	std::vector<ifcopenshell::geometry::filter_t> filters = {
		IfcGeom::entity_filter(true, false, settings.entity_names)
	};

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

		const auto& n = geom_object->transformation().data().components;
		const cgal_placement_t trsf2(
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
			typedef ifcopenshell::geometry::taxonomy::collection cl;
			typedef ifcopenshell::geometry::taxonomy::loop l;
			typedef ifcopenshell::geometry::taxonomy::edge e;
			typedef ifcopenshell::geometry::taxonomy::point3 p;
			auto edge = (e*)((l*)((cl*)((cl*)item)->children[0])->children[0])->children[0];
			if (edge->basis == nullptr && edge->start.which() == 0 && edge->end.which() == 0) {
				auto p0 = boost::get<p>(edge->start);
				auto p1 = boost::get<p>(edge->end);
				Eigen::Vector4d P0;
				P0 << p0.components, 1.;
				Eigen::Vector4d P1;
				P1 << p1.components, 1.;
				auto V0 = n * P0;
				auto V1 = n * P1;
				std::cout << "Axis " << V0(0) << " " << V0(1) << " " << V0(2) << " -> "
					<< V1(0) << " " << V1(1) << " " << V1(2) << std::endl;
				wall_direction = (V1 - V0).head<3>().normalized();
			}
			T2.stop();
		}

		for (auto& g : geom_object->geometry()) {
			auto s = ((ifcopenshell::geometry::CgalShape*) g.Shape())->shape();
			const auto& m = g.Placement().components;
			
			const cgal_placement_t trsf(
				m(0, 0), m(0, 1), m(0, 2), m(0, 3),
				m(1, 0), m(1, 1), m(1, 2), m(1, 3),
				m(2, 0), m(2, 1), m(2, 2), m(2, 3));

			// Apply transformation
			for (auto &vertex : vertices(s)) {
				vertex->point() = vertex->point().transform(trsf).transform(trsf2);
			}

			boost::optional<ifcopenshell::geometry::taxonomy::style> opt_style;
			if (g.hasStyle()) {
				opt_style = g.Style();
			}

			auto T1 = timer.measure("ifc_element_to_nef");
			CGAL::Nef_polyhedron_3<Kernel_> part_nef = ifcopenshell::geometry::utils::create_nef_polyhedron(s);
			T1.stop();

			if (part_nef.is_empty() || !part_nef.is_simple()) {
				continue;
			}

			fn(shape_callback_item{
				geom_object->product(),
				geom_object->guid(),
				geom_object->type(),
				s,
				part_nef,
				opt_style,
				wall_direction,
				openings });

			std::cout << "Processed: " << geom_object->product()->data().toString() << std::endl;
		}
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

int main(int argc, char** argv) {
	geobim_settings settings;
	parse_command_line(settings, argc, argv);

	global_execution_context<CGAL::Simple_cartesian<double>> global_context;
	global_execution_context<Kernel_> global_context_exact;

	std::vector<radius_execution_context> radius_contexts;
	for (double r : settings.radii) {
		radius_contexts.emplace_back(r);
	}

	shape_callback callback;
	if (settings.exact_segmentation) {
		callback.contexts.push_back(&global_context_exact);
	} else {
		callback.contexts.push_back(&global_context);
	}
	for (auto& c : radius_contexts) {
		callback.contexts.push_back(&c);
	}

	process_geometries(settings, callback);

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

		simple_obj_writer write_obj(settings.output_filename + boost::lexical_cast<std::string>(c.radius));
		for (auto& p : style_facet_pairs) {
			write_obj(p.first, p.second.begin(), p.second.end());
		}
	}

	timer.print(std::cout);
}
