/********************************************************************************
 *                                                                              *
 * This file is part of TUDelft Esri GEOBIM.                                    *
 *                                                                              *
 * License: APACHE                                                              *
 *                                                                              *
 ********************************************************************************/

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

typedef CGAL::AABB_face_graph_triangle_primitive<cgal_shape_t, CGAL::Default, CGAL::Tag_false> Primitive;
typedef CGAL::AABB_traits<Kernel_, Primitive> AAbbTraits;
typedef CGAL::AABB_tree<AAbbTraits> AAbbTree;


int main (int argc, char** argv) {

	typedef CGAL::Box_intersection_d::Box_with_handle_d<double, 3, size_t> Box;

	if (argc < 4) {
		std::cerr << "[Error] Usage: '" << argv[0] << " [-d] [-e=IfcWall;IfcSlab;...] [-o] <ifc_file> <out_obj_file> <dilate_radius_0 .. n>'" << std::endl;
		return 1;
	}

	// using Kernel_::FT creates weird segfaults, probably due to how ifopsh constructs the box, coordinates components cannot be shared?
	// std::vector<Kernel_::FT> radii;
	std::vector<double> radii;
	for (int i = 3; i < argc; ++i) {
		// radii.push_back(CGAL::Gmpq(argv[i]));
		radii.push_back(boost::lexical_cast<double>(argv[i]));
	}

	std::string fn = argv[1];
	std::string ofn = argv[2];
	IfcParse::IfcFile f(fn);

	if (!f.good()) {
		std::cerr << "[Error] Unable to parse input file '" << fn << "'";
		return 1;
	}

	ifcopenshell::geometry::settings settings;
	settings.set(ifcopenshell::geometry::settings::USE_WORLD_COORDS, false);
	settings.set(ifcopenshell::geometry::settings::WELD_VERTICES, false);
	settings.set(ifcopenshell::geometry::settings::SEW_SHELLS, true);
	settings.set(ifcopenshell::geometry::settings::CONVERT_BACK_UNITS, true);
	settings.set(ifcopenshell::geometry::settings::DISABLE_TRIANGULATION, true);
	settings.set(ifcopenshell::geometry::settings::DISABLE_OPENING_SUBTRACTIONS, true);

	std::vector<ifcopenshell::geometry::filter_t> no_openings_and_spaces = {
		IfcGeom::entity_filter(false, false, {"IfcOpeningElement", "IfcSpace"})
	};

	std::vector<ifcopenshell::geometry::filter_t> walls_and_slabs = {
		IfcGeom::entity_filter(true, false, {"IfcWall", "IfcSlab"})
	};

	ifcopenshell::geometry::Iterator context_iterator("cgal", settings, &f, walls_and_slabs);

	if (!context_iterator.initialize()) {
		return 1;
	}

	size_t num_created = 0;

	std::vector< CGAL::Nef_nary_union_3< CGAL::Nef_polyhedron_3<Kernel_> > > all_union(radii.size());
	std::vector< CGAL::Nef_polyhedron_3<Kernel_> > cubes;

	AAbbTree tree;

	std::list<ifcopenshell::geometry::taxonomy::style> styles;
	std::map<std::pair<double, std::pair<double, double>>, decltype(styles)::iterator> diffuse_to_style;
	// style 0 is for elements without style annotations
	styles.emplace_back();
	std::vector<Box> facet_boxes_with_styles;
	std::list<cgal_shape_t> triangulated_shape_memory;
	std::map<cgal_shape_t::Facet_handle, decltype(styles)::const_iterator> facet_to_style;


	for (auto& r : radii) {
		auto polycube = ifcopenshell::geometry::utils::create_cube(r);
		cubes.push_back(ifcopenshell::geometry::utils::create_nef_polyhedron(polycube));
	}

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

			size_t style_idx = g.hasStyle() ? styles.size() : 0;

			// Group taxonomy::styles based on diffuse colour since we
			// do not have an equality operator on it.
			decltype(styles)::iterator sit = styles.begin();
			if (g.hasStyle() && g.Style().diffuse) {
				auto cc = (*g.Style().diffuse).components;
				auto c = std::make_pair(cc(0), std::make_pair(cc(1), cc(2)));
				auto it = diffuse_to_style.find(c);
				if (it == diffuse_to_style.end()) {
					styles.push_back(g.Style());
					sit = --styles.end();
					diffuse_to_style.insert({ c, sit });
				} else {
					sit = it->second;
				}
			}

			// A reference is kept to the original shapes in a std::list.
			// Later an aabb tree is used map eroded triangle centroids
			// back to the original elements to preserve semantics.
			triangulated_shape_memory.push_back(s);
			CGAL::Polygon_mesh_processing::triangulate_faces(triangulated_shape_memory.back());
			tree.insert(faces(triangulated_shape_memory.back()).first, faces(triangulated_shape_memory.back()).second, triangulated_shape_memory.back());
			for (auto& f : faces(triangulated_shape_memory.back())) {
				cgal_shape_t::Facet_handle F = f;
				facet_to_style.insert({ F, sit });
			}

			CGAL::Nef_polyhedron_3<Kernel_> part_nef = ifcopenshell::geometry::utils::create_nef_polyhedron(s);
			
			if (part_nef.is_empty()) {
				// std::wcout << "not simple" << std::endl;
				continue;
			}

			if (!part_nef.is_simple()) {
				// std::wcout << "not simple" << std::endl;
				continue;
			}

			auto uit = all_union.begin();
			auto cit = cubes.begin();
			for (; uit != all_union.end(); ++uit, ++cit) {
				uit->add_polyhedron(CGAL::minkowski_sum_3(part_nef, (*cit)));
			}

			std::cout << geom_object->product()->data().toString() << std::endl;
		}
	}

	auto uit = all_union.begin();
	auto rit = radii.begin();
	auto cit = cubes.begin();

	for (; uit != all_union.end(); ++uit, ++rit, ++cit) {

		std::stringstream ss;
		ss << (*rit);
		std::string radius_string = ss.str();
		std::string fn1 = ofn + "-full-" + radius_string + ".off";
		std::string fn2 = ofn + "-dilated-" + radius_string + ".off";
		std::string fn3 = ofn + "-" + radius_string + ".off";

		auto result = uit->get_union();

		auto poly = ifcopenshell::geometry::utils::create_polyhedron(result);

		{
			std::ofstream fs(fn1.c_str());
			fs.precision(12);
			fs << poly;
		}

		Kernel_::FT maxx = -1e9;
		cgal_shape_t::Vertex_handle left_most_vertex;
		for (auto it = poly.edges_begin(); it != poly.edges_end(); ++it) {
			auto vx = it->vertex()->point().cartesian(0);
			if (vx > maxx) {
				maxx = vx;
				left_most_vertex = it->vertex();
			}
		}

		std::set<cgal_shape_t::Facet_handle> empty;
		auto connected = connected_faces(left_most_vertex->halfedge()->facet(), empty);

		std::set<cgal_shape_t::Facet_handle> outer_facets(connected.begin(), connected.end());

		while (std::distance(poly.facets_begin(), poly.facets_end()) > outer_facets.size()) {
			for (auto it = poly.facets_begin(); it != poly.facets_end(); ++it) {
				if (outer_facets.find(&*it) == outer_facets.end()) {
					poly.erase_connected_component(it->halfedge());
					break;
				}
			}
		}

		{
			std::ofstream fs(fn2.c_str());
			fs.precision(12);
			fs << poly;
		}

		auto exterior_nef = ifcopenshell::geometry::utils::create_nef_polyhedron(poly);

		// Create the complement of the Nef by subtracting from its bounding box,
		// see: https://github.com/tudelft3d/ifc2citygml/blob/master/off2citygml/Minkowski.cpp#L23
		auto bounding_box = CGAL::Polygon_mesh_processing::bbox(poly);
		Kernel_::Point_3 bbmin (bounding_box.xmin(), bounding_box.ymin(), bounding_box.zmin());
		Kernel_::Point_3 bbmax (bounding_box.xmax(), bounding_box.ymax(), bounding_box.zmax());
		Kernel_::Vector_3 d(*rit, *rit, *rit);
		bbmin = CGAL::ORIGIN + ((bbmin - CGAL::ORIGIN) - d);
		bbmax = CGAL::ORIGIN + ((bbmax - CGAL::ORIGIN) + d);
		auto bbpoly = ifcopenshell::geometry::utils::create_cube(bbmin, bbmax);
		auto bbnef = ifcopenshell::geometry::utils::create_nef_polyhedron(bbpoly);
		bbnef.extract_regularization();
		auto complement = bbnef - exterior_nef;
		complement.extract_regularization();

		{
			auto complement_poly = ifcopenshell::geometry::utils::create_polyhedron(complement);
			std::ofstream fs("padded");
			fs.precision(12);
			fs << complement_poly;
		}

		auto padded = CGAL::minkowski_sum_3(complement, (*cit));
		padded.extract_regularization();

		{
			auto padded_poly = ifcopenshell::geometry::utils::create_polyhedron(padded);
			std::ofstream fs("padded");
			fs.precision(12);
			fs << padded_poly;
		}

		// @todo I imagine this operation is costly, we can also convert the padded complement to
		// polyhedron, and remove the connected component that belongs to the bbox, then reverse
		// the remaining poly to point to the interior?
		// or: extract_interior?
		exterior_nef -= padded;
		exterior_nef.extract_regularization();
		auto eroded_poly = ifcopenshell::geometry::utils::create_polyhedron(exterior_nef);

		{
			std::ofstream fs(fn3.c_str());
			fs.precision(12);
			fs << eroded_poly;
		}
		
		tree.build();
		tree.accelerate_distance_queries();

		std::vector< std::list<cgal_shape_t::Facet_handle> > facets_by_style(styles.size());

		for (auto &f : faces(eroded_poly)) {
			auto O = CGAL::centroid(
				f->facet_begin()->vertex()->point(),
				f->facet_begin()->next()->vertex()->point(),
				f->facet_begin()->next()->next()->vertex()->point()
			);
			auto pair = tree.closest_point_and_primitive(O);
			auto it = facet_to_style.find(pair.second);
			if (it != facet_to_style.end()) {
				int sid = std::distance(styles.cbegin(), it->second);
				facets_by_style[sid].push_back(f);
			}
		}

		ifcopenshell::geometry::taxonomy::colour BLACK;
		BLACK.components = Eigen::Vector3d(0., 0., 0.);

		{
			std::ofstream fs("colors.obj");
			std::ofstream fs2("colors.mtl");
			fs << "mtllib colors.mtl\n";
			int g = 1;
			int vs = 1;
			auto sit = styles.begin();
			for (auto& li : facets_by_style) {
				fs << "g group-" << g << "\n";
				fs << "usemtl m" << g << "\n";
				fs2 << "newmtl m" << g << "\n";
				auto diff = sit->diffuse.get_value_or(BLACK).components;
				fs2 << "kd " << diff(0) << " " << diff(1) << " " << diff(2) << "\n";
				g++;
				for (auto& f : li) {
					Kernel_::Point_3 points[] = {
						f->facet_begin()->vertex()->point(),
						f->facet_begin()->next()->vertex()->point(),
						f->facet_begin()->next()->next()->vertex()->point()
					};
					for (int i = 0; i < 3; ++i) {
						fs << "v "
							<< points[i].cartesian(0) << " "
							<< points[i].cartesian(1) << " "
							<< points[i].cartesian(2) << "\n";
					}
					fs << "f "
						<< (vs + 0) << " "
						<< (vs + 1) << " "
						<< (vs + 2) << "\n";
					vs += 3;
				}
				++sit;
			}
		}

	}
}
