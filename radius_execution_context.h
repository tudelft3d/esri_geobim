#ifndef RADIUS_EXECUTION_CONTEXT_H
#define RADIUS_EXECUTION_CONTEXT_H

#include "processing.h"
#include "utils.h"

#include <CGAL/Nef_nary_union_3.h>

// State (polyhedra mostly) that are relevant only for one radius
struct radius_execution_context : public execution_context {
	double radius;
	CGAL::Nef_nary_union_3< CGAL::Nef_polyhedron_3<Kernel_> > union_collector;
	CGAL::Nef_polyhedron_3<Kernel_> padding_cube, padding_cube_2, boolean_result, exterior, bounding_box, complement, complement_padded;
	cgal_shape_t polyhedron, polyhedron_exterior;
	enum extract_component { INTERIOR, EXTERIOR };
	bool minkowski_triangles_;

	radius_execution_context(double r, bool narrower = false, bool minkowski_triangles = false);

	IfcUtil::IfcBaseEntity* previous_src = nullptr;
	std::string previous_geom_ref;
	CGAL::Nef_nary_union_3< CGAL::Nef_polyhedron_3<Kernel_> > per_product_collector;
	cgal_placement_t last_place;

	void operator()(shape_callback_item& item);

	// Extract the exterior component of a CGAL Polyhedron
	cgal_shape_t extract(const cgal_shape_t& input, extract_component component) const;

	// Create a bounding box (six-faced Nef poly) around a CGAL Polyhedron
	CGAL::Nef_polyhedron_3<Kernel_> create_bounding_box(const cgal_shape_t& input) const;

	// Completes the boolean union, extracts exterior and erodes padding radius
	void finalize();
};

#endif
