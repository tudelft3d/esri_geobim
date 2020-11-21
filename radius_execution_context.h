#ifndef RADIUS_EXECUTION_CONTEXT_H
#define RADIUS_EXECUTION_CONTEXT_H

#include "processing.h"
#include "utils.h"

#include <CGAL/Nef_nary_union_3.h>

template <typename T>
class lazy_nary_union {
	std::list<T> a;
	CGAL::Nef_nary_union_3<T> b, empty;
	T c;
	int state = 0;

public:
	void add_polyhedron(const T& t) {
		if (state == 0) {
			a.push_back(t);
		} else {
			throw std::runtime_error("");
		}
	}

	const T& get_union() {
		if (state == 0) {
			++state;
			for (auto& p : a) {
				b.add_polyhedron(p);
			}
			++state;
			if (!a.empty()) {
				c = b.get_union();
			}
		}
		return c;
	}

	void clear() {
		state = 0;
		a.clear();
		b = empty;
	}
};

// State (polyhedra mostly) that are relevant only for one radius
struct radius_execution_context : public execution_context {
	double radius;
	CGAL::Nef_nary_union_3< CGAL::Nef_polyhedron_3<Kernel_> > union_collector;
	CGAL::Nef_polyhedron_3<Kernel_> padding_cube, padding_cube_2, boolean_result, exterior, bounding_box, complement, complement_padded;
	cgal_shape_t polyhedron_exterior;
	enum extract_component { INTERIOR, EXTERIOR, LARGEST_AREA, SECOND_LARGEST_AREA };
	bool minkowski_triangles_, no_erosion_, empty_;

	radius_execution_context(double r, bool narrower = false, bool minkowski_triangles = false, bool no_erosion = false);

	IfcUtil::IfcBaseEntity* previous_src = nullptr;
	std::string previous_geom_ref;
	lazy_nary_union<CGAL::Nef_polyhedron_3<Kernel_> > per_product_collector;
	cgal_placement_t last_place;

	void operator()(shape_callback_item& item);

	// Extract the exterior component of a CGAL Polyhedron
	cgal_shape_t extract(const cgal_shape_t& input, extract_component component) const;
	void extract_in_place(cgal_shape_t& input, extract_component component) const;

	// Create a bounding box (six-faced Nef poly) around a CGAL Polyhedron
	CGAL::Nef_polyhedron_3<Kernel_> create_bounding_box(const cgal_shape_t& input) const;

	// Completes the boolean union, extracts exterior and erodes padding radius
	void finalize();

	bool empty() const { return empty_; }
};

#endif
