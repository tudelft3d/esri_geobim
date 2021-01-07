#ifndef CONTEXT_H
#define CONTEXT_H

#include <ifcgeom/kernels/cgal/CgalKernel.h>

#include <boost/optional.hpp>

#include <list>
#include <string>

// Generated for every representation *item* in the IFC file
struct shape_callback_item {
	IfcUtil::IfcBaseEntity* src;
	std::string id, type, part_reference, geom_reference;
	cgal_placement_t transformation;
	cgal_shape_t polyhedron;
	const ifcopenshell::geometry::taxonomy::style* style;
	boost::optional<Eigen::Vector3d> wall_direction;
	std::list<shape_callback_item*> openings;

	bool to_nef_polyhedron(CGAL::Nef_polyhedron_3<Kernel_>& nef);
};

// Prototype of a context to which processed shapes will be fed
struct execution_context {
	virtual void operator()(shape_callback_item&) = 0;
};

#endif