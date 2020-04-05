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

int main (int argc, char** argv) {
	if (argc != 3) {
		std::cerr << "[Error] Usage: '" << argv[0] << " <ifc_file> <out_obj_file>'";
		return 1;
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

	CGAL::Nef_polyhedron_3<Kernel_> all_union;

	auto polycube = ifcopenshell::geometry::utils::create_cube(0.01);
	auto cube = ifcopenshell::geometry::utils::create_nef_polyhedron(polycube);

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

		/*
		std::stringstream ss;
		ss << geom_object->product()->data().toString();
		auto sss = ss.str();
		std::wcout << sss.c_str() << std::endl;
		*/

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

			CGAL::Nef_polyhedron_3<Kernel_> part_nef = ifcopenshell::geometry::utils::create_nef_polyhedron(s);
			
			if (part_nef.is_empty()) {
				// std::wcout << "not simple" << std::endl;
				continue;
			}

			if (!part_nef.is_simple()) {
				// std::wcout << "not simple" << std::endl;
				continue;
			}

			all_union += CGAL::minkowski_sum_3(part_nef, cube);

			std::cout << geom_object->product()->data().toString() << std::endl;
		}
	}

	auto poly = ifcopenshell::geometry::utils::create_polyhedron(all_union);

	{
		std::ofstream fs("poly.off");
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

	/*
	std::set<cgal_shape_t::Halfedge_handle> outer_edges;
	for (auto& f : outer_facets) {
		auto begin = f->facet_begin();
		auto it = begin;
		do {
			outer_edges.insert(it);
		} while (++it != begin);
	}
	*/

	while (std::distance(poly.facets_begin(), poly.facets_end()) > outer_facets.size()) {
		for (auto it = poly.facets_begin(); it != poly.facets_end(); ++it) {
			if (outer_facets.find(&*it) == outer_facets.end()) {
				poly.erase_connected_component(it->halfedge());
				break;
			}
		}
	}

	{
		std::ofstream fs(ofn.c_str());
		fs.precision(12);
		fs << poly;
	}
}
