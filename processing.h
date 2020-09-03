#ifndef PROCESSING_H
#define PROCESSING_H

#include "timer.h"
#include "settings.h"

#define ENSURE_2ND_OP_NARROWER

#ifdef ENSURE_2ND_OP_NARROWER
#define MAKE_OP2_NARROWER -2e-7
#else
#define MAKE_OP2_NARROWER
#endif

#include <ifcgeom/kernels/cgal/CgalKernel.h>
#include <ifcgeom/schema_agnostic/IfcGeomFilter.h>

#include <string>

// Generated for every representation *item* in the IFC file
struct shape_callback_item {
	IfcUtil::IfcBaseEntity* src;
	std::string id, type, geom_reference;
	cgal_placement_t transformation;
	cgal_shape_t polyhedron;
	boost::optional<ifcopenshell::geometry::taxonomy::style> style;
	boost::optional<Eigen::Vector3d> wall_direction;
	std::list<shape_callback_item*> openings;

	bool to_nef_polyhedron(CGAL::Nef_polyhedron_3<Kernel_>& nef);
};

// Prototype of a context to which processed shapes will be fed
struct execution_context {
	virtual void operator()(shape_callback_item&) = 0;
};

struct debug_writer : public execution_context {
	void operator()(shape_callback_item& item);
};

// An execution context that stores processed items from the file.
struct capturing_execution_context : public execution_context {
	std::list<shape_callback_item> items;

	void operator()(shape_callback_item& item) {
		items.push_back(item);
	}
};

// A structure for recieving processed shapes simply defers to a vector of contexts
struct shape_callback {
	std::vector<execution_context*> contexts;

	void operator()(shape_callback_item& item) {
		for (auto& c : contexts) {
			(*c)(item);
		}
	}
};

// Interprets IFC geometries by means of IfcOpenShell CGAL and
// pass result to callback
int process_geometries(geobim_settings& settings, const std::function <void(shape_callback_item&)>& fn);

#endif
