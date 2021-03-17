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

namespace {
	template <typename W>
	struct output_writer {
		W& w;

		template <typename R>
		void operator()(R& r) {
			for (auto& p : r) {
				w(p.first, p.second.begin(), p.second.end());
			}
		}
	};
}

int main(int argc, char** argv) {
	Logger::SetOutput(&std::cerr, &std::cerr);
	Logger::Verbosity(Logger::LOG_NOTICE);

	geobim_settings settings;
	parse_command_line(settings, argc, argv);

	// global_execution_context<CGAL::Simple_cartesian<double>> global_context;
	// Epick seems to have better performance?
	global_execution_context<CGAL::Epick> global_context;
	global_execution_context<Kernel_> global_context_exact;

	shape_callback callback_global;
	if (settings.exact_segmentation) {
		callback_global.contexts.push_back(&global_context_exact);
	} else {
		callback_global.contexts.push_back(&global_context);
	}

#ifdef GEOBIM_DEBUG
	callback.contexts.push_back(new debug_writer);
#endif

	std::unique_ptr<process_geometries> p;

	if (settings.radii.empty()) {
		// @todo
		/*
		auto cec = new capturing_execution_context;
		callback.contexts.push_back(cec);
		p = std::make_unique<process_geometries>(settings);
		(*p)(std::ref(callback));
		auto R = binary_search(cec->items.begin(), cec->items.end(), { "0.001", "0.2" });
		std::cout << "Largest gap found with R / 2 ~ " << R << std::endl;
		*/
	} else {
		std::vector<std::unique_ptr<radius_execution_context>> radius_contexts;
		bool first = true;
		for (auto& r : settings.radii) {
			// 2nd is narrower (depending on ifdef above, appears to be necessary).
			radius_contexts.push_back(std::make_unique<radius_execution_context>(r, radius_settings()
				.set(radius_settings::NARROWER, !first)
				.set(radius_settings::MINKOWSKI_TRIANGLES, settings.minkowski_triangles)
				.set(radius_settings::NO_EROSION, settings.no_erosion)
				.set(radius_settings::SPHERE, settings.spherical_padding)));
			first = false;
			if (settings.threads) {
				radius_contexts.back()->set_threads(*settings.threads);
			}
		}

		capturing_execution_context cec;

		p = std::make_unique<process_geometries>(settings);
		(*p)(std::ref(cec));

		{
			shape_callback callback;
			for (auto& c : radius_contexts) {
				callback.contexts.push_back(&*c);
			}

			cec.run(std::ref(callback));
		}

		cec.run(std::ref(callback_global));

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
			if (c->empty()) {
				continue;
			}

			c->finalize();
		}

		p.reset();

		for (auto& c : radius_contexts) {
			auto T0 = timer::measure("semantic_segmentation");

			std::list<item_info*> all_infos;

			if (settings.exact_segmentation) {
				all_infos = global_context_exact.all_item_infos();
			}
			else {
				all_infos = global_context.all_item_infos();
			}

			// pop the first 'empty' info
			all_infos.pop_front();

			city_json_writer write_city(settings.output_filename + c->radius_str);
			simple_obj_writer write_obj(settings.output_filename + c->radius_str);
			external_element_collector write_elem(settings.output_filename + ".external", all_infos);

			boost::variant<
				global_execution_context<Kernel_>::segmentation_return_type,
				global_execution_context<Kernel_>::segmentation_return_type_2
			> style_facet_pairs;

			if (settings.spherical_padding) {
				if (settings.exact_segmentation) {
					style_facet_pairs = global_context_exact.segment(c->polyhedron_exterior_nm);
				} else {
					style_facet_pairs = global_context.segment(c->polyhedron_exterior_nm);
				}
				write_city.point_lookup = &c->polyhedron_exterior_nm.points;
				write_obj.point_lookup = &c->polyhedron_exterior_nm.points;
				write_elem.point_lookup = &c->polyhedron_exterior_nm.points;
			} else if (settings.exact_segmentation) {
				style_facet_pairs = global_context_exact.segment(c->polyhedron_exterior);
			} else {
				style_facet_pairs = global_context.segment(c->polyhedron_exterior);
			}

			T0.stop();

			{	
				output_writer<city_json_writer> vis{ write_city };
				boost::apply_visitor(vis, style_facet_pairs);
				write_city.finalize();
			}

			{
				output_writer<simple_obj_writer> vis{ write_obj };
				boost::apply_visitor(vis, style_facet_pairs);
				write_obj.finalize();
			}

			{
				output_writer<external_element_collector> vis{ write_elem };
				boost::apply_visitor(vis, style_facet_pairs);
				write_elem.finalize();
			}
		}

		auto T2 = timer::measure("difference_overlay");
		auto it = radius_contexts.begin();
		for (auto jt = it + 1; jt != radius_contexts.end(); ++it, ++jt) {
			radius_comparison difference(**it, **jt, 0.001);
			simple_obj_writer obj("difference-"
				+ (*it)->radius_str + "-"
				+ (*jt)->radius_str);
			obj(nullptr, difference.difference_poly.facets_begin(), difference.difference_poly.facets_end());
		}
		T2.stop();

	}

	timer::print(std::cout);
}
