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

typedef CGAL::AABB_face_graph_triangle_primitive<cgal_shape_t, CGAL::Default, CGAL::Tag_false> Primitive;
typedef CGAL::AABB_traits<Kernel_, Primitive> AAbbTraits;
typedef CGAL::AABB_tree<AAbbTraits> AAbbTree;

namespace po = boost::program_options;

// Global settings that are derived from the command line invocation
// parsed by Boost.ProgramOptions
struct geobim_settings {
	std::string input_filename, output_filename;
	std::vector<double> radii;
	bool apply_openings, debug;
	ifcopenshell::geometry::settings settings;
	std::set<std::string> entity_names;
	IfcParse::IfcFile* file;
};

// Generated for every representation *item* in the IFC file
struct shape_callback_item {
	std::string id, type;
	cgal_shape_t polyhedron;
	CGAL::Nef_polyhedron_3<Kernel_> nef_polyhedron;
	boost::optional<ifcopenshell::geometry::taxonomy::style> style;
};

// Prototype of a context to which processed shapes will be fed
struct execution_context {
	virtual void operator()(shape_callback_item&) = 0;
};

// State that is relevant accross the different radii for which the
// program is executed. Mostly related for the mapping back to
// semantics from the original IFC input.
struct global_execution_context : public execution_context {
	AAbbTree tree;

	std::list<ifcopenshell::geometry::taxonomy::style> styles;
	std::map<std::pair<double, std::pair<double, double>>, decltype(styles)::iterator> diffuse_to_style;

	// A reference is kept to the original shapes in a std::list.
	// Later an aabb tree is used map eroded triangle centroids
	// back to the original elements to preserve semantics.
	std::list<cgal_shape_t> triangulated_shape_memory;
	std::map<cgal_shape_t::Facet_handle, decltype(styles)::const_iterator> facet_to_style;

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

		triangulated_shape_memory.push_back(item.polyhedron);
		CGAL::Polygon_mesh_processing::triangulate_faces(triangulated_shape_memory.back());
		tree.insert(faces(triangulated_shape_memory.back()).first, faces(triangulated_shape_memory.back()).second, triangulated_shape_memory.back());
		for (auto& f : faces(triangulated_shape_memory.back())) {
			cgal_shape_t::Facet_handle F = f;
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

		for (auto &f : faces(input)) {
			auto O = CGAL::centroid(
				f->facet_begin()->vertex()->point(),
				f->facet_begin()->next()->vertex()->point(),
				f->facet_begin()->next()->next()->vertex()->point()
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
		union_collector.add_polyhedron(CGAL::minkowski_sum_3(item.nef_polyhedron, padding_cube));
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
		boolean_result = union_collector.get_union();
		polyhedron = ifcopenshell::geometry::utils::create_polyhedron(boolean_result);
		polyhedron_exterior = extract(polyhedron, EXTERIOR);
		exterior = ifcopenshell::geometry::utils::create_nef_polyhedron(polyhedron_exterior);
		bounding_box = create_bounding_box(polyhedron);
		complement = bounding_box - exterior;
		complement.extract_regularization();
		complement_padded = CGAL::minkowski_sum_3(complement, padding_cube);
		complement_padded.extract_regularization();

		auto start = std::chrono::steady_clock::now();
#if 0
		// @todo I imagine this operation is costly, we can also convert the padded complement to
		// polyhedron, and remove the connected component that belongs to the bbox, then reverse
		// the remaining poly to point to the interior?
		exterior -= complement_padded;
#else
		// Marginally faster.

		// Re above: extracting the interior shell did not prove to be reliable even with
		// the undocumented function convert_inner_shell_to_polyhedron(). Therefore we
		// subtract from the padded box as that will have lower complexity than above.
		// Mark_bounded_volumes on the completement also did not work.
		auto box_padded = CGAL::minkowski_sum_3(bounding_box, padding_cube);
		auto exterior = box_padded - complement_padded;
#endif
		auto end = std::chrono::steady_clock::now();
		std::cout << "conversion back took " << std::chrono::duration<double, std::milli>(end-start).count() << " ms" << endl;

		exterior.extract_regularization();
		polyhedron_exterior = ifcopenshell::geometry::utils::create_polyhedron(exterior);

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
int process_geometries(geobim_settings& settings, Fn fn) {
	std::vector<ifcopenshell::geometry::filter_t> filters = {
		IfcGeom::entity_filter(true, false, settings.entity_names)
	};

	ifcopenshell::geometry::Iterator context_iterator("cgal", settings.settings, settings.file, filters);

	if (!context_iterator.initialize()) {
		return 1;
	}

	size_t num_created = 0;

	for (;; ++num_created) {
		bool has_more = true;
		if (num_created) {
			has_more = context_iterator.next();
		}
		ifcopenshell::geometry::NativeElement* geom_object = nullptr;
		if (has_more) {
			geom_object = context_iterator.get_native();
		}
		if (!geom_object) {
			break;
		}

		for (auto& g : geom_object->geometry()) {
			auto s = ((ifcopenshell::geometry::CgalShape*) g.Shape())->shape();
			const auto& m = g.Placement().components;
			const auto& n = geom_object->transformation().data().components;

			const cgal_placement_t trsf(
				m(0, 0), m(0, 1), m(0, 2), m(0, 3),
				m(1, 0), m(1, 1), m(1, 2), m(1, 3),
				m(2, 0), m(2, 1), m(2, 2), m(2, 3));

			const cgal_placement_t trsf2(
				n(0, 0), n(0, 1), n(0, 2), n(0, 3),
				n(1, 0), n(1, 1), n(1, 2), n(1, 3),
				n(2, 0), n(2, 1), n(2, 2), n(2, 3));

			// Apply transformation
			for (auto &vertex : vertices(s)) {
				vertex->point() = vertex->point().transform(trsf).transform(trsf2);
			}

			boost::optional<ifcopenshell::geometry::taxonomy::style> opt_style;
			if (g.hasStyle()) {
				opt_style = g.Style();
			}

			CGAL::Nef_polyhedron_3<Kernel_> part_nef = ifcopenshell::geometry::utils::create_nef_polyhedron(s);

			if (part_nef.is_empty() || !part_nef.is_simple()) {
				continue;
			}

			fn(shape_callback_item{
				geom_object->guid(),
				geom_object->type(),
				s,
				part_nef,
				opt_style });

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

// OBJ writer for CGAL facets paired with a style
struct simple_obj_writer {
	int group_id = 1;
	int vertex_count = 1;
	std::ofstream obj, mtl;
	ifcopenshell::geometry::taxonomy::colour BLACK;

	simple_obj_writer(const std::string& fn_prefix)
		: obj((fn_prefix + ".obj").c_str())
		, mtl((fn_prefix + ".mtl").c_str())
	{
		obj << "mtllib " << fn_prefix << ".mtl\n";
		BLACK.components = Eigen::Vector3d(0., 0., 0.);
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
			auto& f = *it;
			Kernel_::Point_3 points[] = {
				f->facet_begin()->vertex()->point(),
				f->facet_begin()->next()->vertex()->point(),
				f->facet_begin()->next()->next()->vertex()->point()
			};
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

int main(int argc, char** argv) {
	geobim_settings settings;
	parse_command_line(settings, argc, argv);

	global_execution_context global_context;
	std::vector<radius_execution_context> radius_contexts;
	for (double r : settings.radii) {
		radius_contexts.emplace_back(r);
	}

	shape_callback callback;
	callback.contexts.push_back(&global_context);
	for (auto& c : radius_contexts) {
		callback.contexts.push_back(&c);
	}

	process_geometries(settings, callback);

	global_context.finalize();

	for (auto& c : radius_contexts) {
		c.finalize();
		auto style_facet_pairs = global_context.segment(c.polyhedron_exterior);

		simple_obj_writer write_obj(settings.output_filename + boost::lexical_cast<std::string>(c.radius));
		for (auto& p : style_facet_pairs) {
			write_obj(p.first, p.second.begin(), p.second.end());
		}
	}
}
