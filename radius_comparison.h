#ifndef RADIUS_COMPARISON_H
#define RADIUS_COMPARISON_H

#include "radius_execution_context.h"

#include <ifcgeom/kernels/cgal/CgalKernel.h>

#include <CGAL/Polygon_mesh_processing/measure.h>

// A comparison between two exterior-shell polyhedra to identify gaps
// in a facade that are filled by a larger dilation radius but left 
// open when using the smaller radius.
struct radius_comparison {
	struct hollow_solid {
		typedef CGAL::Nef_polyhedron_3<Kernel_> nef;

		nef bbox, complement, complement_padded, inner, hollow, cube;

		double D;

		hollow_solid(radius_execution_context& a, double d);
	};

	typedef CGAL::Nef_polyhedron_3<Kernel_> nef;
	typedef CGAL::Polyhedron_3<Kernel_> poly;

	nef difference_nef;
	poly difference_poly;
	hollow_solid A, B;

	radius_comparison(radius_execution_context& a, radius_execution_context& b, double d);
};

template <typename It>
double initialize_radius_context_and_get_volume_with_cache(It first, It second, const std::string& radius) {
	static std::map<std::string, double> cache_;

	auto it = cache_.find(radius);
	if (it != cache_.end()) {
		std::cout << "Used cache for R=" << radius << " V=" << it->second << std::endl;
		return it->second;
	}

	radius_execution_context rec(radius);
	// Unfortunately for_each() is by value so needs to be wrapped in a lambda with side effects
	std::for_each(first, second, [&rec](auto& v) {
		rec(v);
	});
	rec.finalize();
	double V = CGAL::to_double(CGAL::Polygon_mesh_processing::volume(rec.polyhedron_exterior));

	std::cout << "Calculated for R=" << radius << " V=" << V << std::endl;

	cache_.insert(it, { radius , V });
	return V;
}

template <typename It>
std::string binary_search(It first, It second, std::pair<std::string, std::string> range) {
	std::string abc[3];
	std::tie(abc[0], abc[2]) = range;

	std::cout << "Testing " << abc[0] << " and " << abc[2] << std::endl;

	auto a_vol = initialize_radius_context_and_get_volume_with_cache(first, second, abc[0]);
	auto b_vol = initialize_radius_context_and_get_volume_with_cache(first, second, abc[2]);

	if (a_vol * 1.1 < b_vol) {
		double a = boost::lexical_cast<double>(abc[0]);
		double c = boost::lexical_cast<double>(abc[2]);
		auto b = (a + c) / 2.;

		abc[1] = boost::lexical_cast<std::string>(b);

		if ((c - a) < 1.e-4) {
			std::cout << "Terminating search at " << abc[0] << " and " << abc[2] << std::endl;
			return boost::lexical_cast<std::string>((a + c) / 2.);
		}

		for (int i = 1; i >= 0; --i) {
			auto r = binary_search(first, second, { abc[i], abc[i + 1] });
			if (r != abc[i + 1]) {
				return r;
			}
		}
	}

	return abc[2];
}

#endif
