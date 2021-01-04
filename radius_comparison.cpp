#include "radius_comparison.h"

#include "writer.h"
#include <CGAL/minkowski_sum_3.h>

radius_comparison::hollow_solid::hollow_solid(radius_execution_context & a, double d) {
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

#define MAKE_OP2_NARROWER -2e-7

radius_comparison::radius_comparison(radius_execution_context & a, radius_execution_context & b, double d)
	: A(a, d), B(b, d + MAKE_OP2_NARROWER) {
	difference_nef = B.hollow - A.hollow;
	difference_nef.extract_regularization();
	difference_poly = ifcopenshell::geometry::utils::create_polyhedron(difference_nef);
}
