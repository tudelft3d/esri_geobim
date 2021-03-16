#ifndef PROCESSING_H
#define PROCESSING_H

#include "timer.h"
#include "context.h"
#include "settings.h"
#include "opening_collector.h"

#include <ifcgeom/kernels/cgal/CgalKernel.h>
#include <ifcgeom/schema_agnostic/IfcGeomFilter.h>
#include <ifcgeom/schema_agnostic/IfcGeomIterator.h>

#include <string>

template <typename K>
struct non_manifold_polyhedron {
	typedef CGAL::Point_3<K> P;
	std::vector<P> points;
	std::vector<std::vector<size_t>> indices;
};

// Light weight representation to be stored in global exec context
struct rgb : public std::pair<double, std::pair<double, double>> {
	rgb(double r, double g, double b)
		: std::pair<double, std::pair<double, double>>(r, { g, b }) {}
	rgb(const Eigen::Vector3d& v)
		: std::pair<double, std::pair<double, double>>(v(0), { v(1), v(2) }) {}

	double& r() { return this->first; }
	const double& r() const { return this->first; }

	double& g() { return this->second.first; }
	const double& g() const { return this->second.first; }

	double& b() { return this->second.second; }
	const double& b() const { return this->second.second; }
};

// Light weight representation to be stored in global exec context
struct item_info {
	// @nb we don't store a reference to the ifcopenshell entity instance so the files can be freed from memory
	// we can store a const reference to the ifcopenshell latebound schema type names.
	const std::string& entity_type;
	std::string guid;
	rgb* diffuse;
};

struct debug_writer : public execution_context {
	void operator()(shape_callback_item* item);
};

// An execution context that stores processed items from the file.
struct capturing_execution_context : public execution_context {
	std::list<shape_callback_item*> items;

	void operator()(shape_callback_item* item) {
		items.push_back(item);
	}

	template <typename Fn>
	void run(Fn fn) {
		for (auto& i : items) {
			fn(i);
		}
	}
};

// A structure for recieving processed shapes simply defers to a vector of contexts
struct shape_callback {
	std::vector<execution_context*> contexts;

	void operator()(shape_callback_item* item) {
		for (auto& c : contexts) {
			(*c)(item);
		}
	}
};

// Interprets IFC geometries by means of IfcOpenShell CGAL and
// pass result to callback
struct process_geometries {
	geobim_settings settings;
	opening_collector all_openings;
	std::list<ifcopenshell::geometry::Iterator*> iterators;

	process_geometries(geobim_settings&);
	~process_geometries();
	int operator()(const std::function<void(shape_callback_item*)>&);
};

#endif
