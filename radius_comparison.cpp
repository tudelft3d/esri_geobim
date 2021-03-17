#include "radius_comparison.h"

#include "writer.h"
#include <CGAL/minkowski_sum_3.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>

CGAL::Nef_polyhedron_3<Kernel_> create_bounding_box(const cgal_shape_t & input, double radius);

radius_comparison::hollow_solid::hollow_solid(radius_execution_context & a, double d) {
	D = d;

	{
		cgal_shape_t pl;
		simple_obj_writer obj("hollow-input-" + a.radius_str);
		obj(nullptr, a.polyhedron_exterior.facets_begin(), a.polyhedron_exterior.facets_end());
	}

	bbox = create_bounding_box(a.polyhedron_exterior, a.radius);
	auto exterior = ifcopenshell::geometry::utils::create_nef_polyhedron(a.polyhedron_exterior);

	complement = bbox - exterior;
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
	hollow = exterior - inner;
	hollow.extract_regularization();
}

#define MAKE_OP2_NARROWER -2e-7

radius_comparison::radius_comparison(radius_execution_context & a, radius_execution_context & b, double d)
	: A(a, d), B(b, d + MAKE_OP2_NARROWER) {

	{
		cgal_shape_t pl;
		CGAL::convert_nef_polyhedron_to_polygon_mesh(A.hollow, pl);
		CGAL::Polygon_mesh_processing::triangulate_faces(pl);
		simple_obj_writer obj("debug-hollow-a");
		obj(nullptr, pl.facets_begin(), pl.facets_end());
	}

	{
		cgal_shape_t pl;
		CGAL::convert_nef_polyhedron_to_polygon_mesh(B.hollow, pl);
		CGAL::Polygon_mesh_processing::triangulate_faces(pl);
		simple_obj_writer obj("debug-hollow-b");
		obj(nullptr, pl.facets_begin(), pl.facets_end());
	}

	difference_nef = B.hollow - A.hollow;
	difference_nef.extract_regularization();
	// difference_poly = ifcopenshell::geometry::utils::create_polyhedron(difference_nef);
	CGAL::convert_nef_polyhedron_to_polygon_mesh(difference_nef, difference_poly);
	CGAL::Polygon_mesh_processing::triangulate_faces(difference_poly);
}
