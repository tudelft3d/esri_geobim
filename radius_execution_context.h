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

#include <bitset>

struct radius_settings : std::bitset<4> {
	enum V {
		NARROWER, MINKOWSKI_TRIANGLES, NO_EROSION, SPHERE
	};
	radius_settings& set(V v, bool b) {
		(*this)[(int) v] = b;
		return *this;
	}
	bool get(V v) const {
		return (*this)[(int) v];
	}
};

// State (polyhedra mostly) that are relevant only for one radius
struct radius_execution_context : public execution_context {
	radius_settings settings_;
	std::string radius_str;
	double radius, ico_edge_length;
	CGAL::Nef_polyhedron_3<Kernel_> exterior;
	cgal_shape_t polyhedron_exterior;
	enum extract_component { INTERIOR, EXTERIOR, LARGEST_AREA, SECOND_LARGEST_AREA };
	bool minkowski_triangles_, no_erosion_, empty_;

	// @todo this sets ico_edge_length
	CGAL::Nef_polyhedron_3<Kernel_> construct_padding_volume_(const boost::optional<double>& = boost::none);

	boost::optional<size_t> threads_;

	radius_execution_context(const std::string& radius, radius_settings=radius_settings());
	radius_execution_context(const radius_execution_context&) = delete;
	radius_execution_context& operator=(const radius_execution_context&) = delete;


	IfcUtil::IfcBaseEntity* previous_src = nullptr;
	std::string previous_geom_ref;
	// lazy_nary_union<CGAL::Nef_polyhedron_3<Kernel_> > per_product_collector;
	cgal_placement_t last_place;

	std::vector< std::future<void> > threadpool_;
	std::list< std::pair<CGAL::Nef_polyhedron_3<Kernel_>, CGAL::Nef_polyhedron_3<Kernel_> >  > padding_volumes_;
	// Extra indirection for easy swaps
	// Not used currently
	std::vector< std::pair<CGAL::Nef_polyhedron_3<Kernel_>, CGAL::Nef_polyhedron_3<Kernel_> >* > padding_vol_ptr_;

	void set_threads(size_t n);
	
	struct geometry_reference {
		IfcUtil::IfcBaseEntity* target;
		cgal_placement_t inverse, own;
	};
	typedef std::list< CGAL::Nef_polyhedron_3<Kernel_> > result_list_t;

	std::map<std::string, IfcUtil::IfcBaseEntity*> first_product_for_geom_id;
	std::map<IfcUtil::IfcBaseEntity*, geometry_reference> reused_products;
	std::map<IfcUtil::IfcBaseEntity*, cgal_placement_t> placements;
	std::map<IfcUtil::IfcBaseEntity*, result_list_t> product_geometries;

	// + map for reused items
	// + finalize() does insertion in n-ary_bool
	void operator()(shape_callback_item* item);

	// Extract the exterior component of a CGAL Polyhedron
	cgal_shape_t extract(const cgal_shape_t& input, extract_component component) const;
	void extract_in_place(cgal_shape_t& input, extract_component component) const;

	// Create a bounding box (six-faced Nef poly) around a CGAL Polyhedron
	// CGAL::Nef_polyhedron_3<Kernel_> create_bounding_box(const cgal_shape_t& input) const;

	// Completes the boolean union, extracts exterior and erodes padding radius
	void finalize();

	bool empty() const { return empty_; }
};

#endif
