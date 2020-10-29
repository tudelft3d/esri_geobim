/********************************************************************************
 *                                                                              *
 * This file is part of TUDelft Esri GEOBIM.                                    *
 *                                                                              *
 * License: APACHE                                                              *
 *                                                                              *
 ********************************************************************************/

#include "timer.h"
#include "writer.h"
#include "settings.h"
#include "processing.h"
#include "radius_execution_context.h"
#include "global_execution_context.h"
#include "radius_comparison.h"

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>

int main(int argc, char** argv) {
	Logger::SetOutput(&std::cerr, &std::cerr);
	Logger::Verbosity(Logger::LOG_NOTICE);

	geobim_settings settings;
	parse_command_line(settings, argc, argv);

	// global_execution_context<CGAL::Simple_cartesian<double>> global_context;
	// Epick seems to have better performance?
	global_execution_context<CGAL::Epick> global_context;
	global_execution_context<Kernel_> global_context_exact;

	shape_callback callback;
	if (settings.exact_segmentation) {
		callback.contexts.push_back(&global_context_exact);
	} else {
		callback.contexts.push_back(&global_context);
	}

#ifdef GEOBIM_DEBUG
	callback.contexts.push_back(new debug_writer);
#endif

	if (settings.radii.empty()) {
		auto cec = new capturing_execution_context;
		callback.contexts.push_back(cec);
		process_geometries(settings, callback);
		auto R = binary_search(cec->items.begin(), cec->items.end(), { 1.e-3, 0.2 });
		std::cout << "Largest gap found with R ~ " << (R * 2.) << std::endl;
	} else {
		std::vector<radius_execution_context> radius_contexts;
		bool first = true;
		for (double r : settings.radii) {
			// 2nd is narrower (depending on ifdef above, appears to be necessary).
			radius_contexts.emplace_back(r, !first, settings.minkowski_triangles, settings.no_erosion);
			first = false;
		}

		for (auto& c : radius_contexts) {
			callback.contexts.push_back(&c);
		}

		process_geometries(settings, callback);

		std::cout << "done processing geometries" << std::endl;

		for (auto& f : settings.file) {
			delete f;
		}		

		auto T1 = timer::measure("semantic_segmentation");
		if (settings.exact_segmentation) {
			global_context_exact.finalize();
		} else {
			global_context.finalize();
		}
		T1.stop();

		for (auto& c : radius_contexts) {
			if (c.empty()) {
				continue;
			}

			c.finalize();

			auto T0 = timer::measure("semantic_segmentation");
			global_execution_context<Kernel_>::segmentation_return_type style_facet_pairs;
			if (settings.exact_segmentation) {
				style_facet_pairs = global_context_exact.segment(c.polyhedron_exterior);
			} else {
				style_facet_pairs = global_context.segment(c.polyhedron_exterior);
			}
			T0.stop();

			city_json_writer write_city(settings.output_filename + boost::lexical_cast<std::string>(c.radius));
			for (auto& p : style_facet_pairs) {
				write_city(p.first, p.second.begin(), p.second.end());
			}

			simple_obj_writer write_obj(settings.output_filename + boost::lexical_cast<std::string>(c.radius));
			for (auto& p : style_facet_pairs) {
				write_obj(p.first, p.second.begin(), p.second.end());
			}
		}

		auto T2 = timer::measure("difference_overlay");
		auto it = radius_contexts.begin();
		for (auto jt = it + 1; jt != radius_contexts.end(); ++it, ++jt) {
			radius_comparison difference(*it, *jt, 0.001);
			simple_obj_writer obj("difference-"
				+ boost::lexical_cast<std::string>(it->radius) + "-"
				+ boost::lexical_cast<std::string>(jt->radius));
			obj(nullptr, difference.difference_poly.facets_begin(), difference.difference_poly.facets_end());
		}
		T2.stop();

	}

	timer::print(std::cout);
}
