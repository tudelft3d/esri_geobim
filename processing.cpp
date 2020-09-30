#include <ifcconvert/validation_utils.h>
#include <ifcgeom/schema_agnostic/IfcGeomIterator.h>

#include "processing.h"
#include "writer.h"
#include "opening_collector.h"

bool shape_callback_item::to_nef_polyhedron(CGAL::Nef_polyhedron_3<Kernel_>& nef) {
	auto T1 = timer::measure("ifc_element_to_nef");
	nef = ifcopenshell::geometry::utils::create_nef_polyhedron(polyhedron);
	T1.stop();

	if (nef.is_empty() || !nef.is_simple()) {
		nef.clear();
		return false;
	}

	return true;
}

void debug_writer::operator()(shape_callback_item& item) {
	simple_obj_writer obj("debug-" + boost::lexical_cast<std::string>(item.src->data().id()));
	obj(nullptr, item.polyhedron.facets_begin(), item.polyhedron.facets_end());
}

// Interprets IFC geometries by means of IfcOpenShell CGAL and
// pass result to callback
int process_geometries(geobim_settings & settings, const std::function<void(shape_callback_item&)>& fn) {

	// @todo static now to prevent being gc'ed
	static opening_collector all_openings(settings.file);

	// Capture all openings beforehand, they are later assigned to the
	// building elements.
	auto opening_settings = settings;
	opening_settings.entity_names.emplace();
	opening_settings.entity_names->insert("IfcOpeningElement");
	opening_settings.entity_names_included = true;
	if (settings.apply_openings_posthoc && settings.entity_names != opening_settings.entity_names) {
		process_geometries(opening_settings, std::ref(all_openings));
	}

	std::vector<ifcopenshell::geometry::filter_t> filters;
	if (settings.entity_names) {
		if (!settings.entity_names_included) {
			settings.entity_names->insert("IfcSpace");
			settings.entity_names->insert("IfcOpeningElement");
		}
		filters.push_back(IfcGeom::entity_filter(settings.entity_names_included, false, *settings.entity_names));
	} else {
		filters.push_back(IfcGeom::entity_filter(false, false, { "IfcSpace", "IfcOpeningElement" }));
	}

	for (auto f : settings.file) {

		ifcopenshell::geometry::Iterator context_iterator("cgal", settings.settings, f, filters);

		auto T = timer::measure("ifc_geometry_processing");
		if (!context_iterator.initialize()) {
			return 1;
		}
		T.stop();

		size_t num_created = 0;

		auto axis_settings = settings.settings;
		axis_settings.set(ifcopenshell::geometry::settings::EXCLUDE_SOLIDS_AND_SURFACES, true);
		axis_settings.set(ifcopenshell::geometry::settings::INCLUDE_CURVES, true);
		auto geometry_mapper = ifcopenshell::geometry::impl::mapping_implementations().construct(f, axis_settings);

		for (;; ++num_created) {
			bool has_more = true;
			if (num_created) {
				auto T0 = timer::measure("ifc_geometry_processing");
				has_more = context_iterator.next();
				T0.stop();
			}
			ifcopenshell::geometry::NativeElement* geom_object = nullptr;
			if (has_more) {
				geom_object = context_iterator.get_native();
			}
			if (!geom_object) {
				break;
			}

			if (geom_object->guid() == "3es57B9Kr3nxL4uBITV$0e") {
				std::cout << "NOTICE Skipping: " << geom_object->product()->data().toString() << std::endl;
				continue;
			}

			std::cout << "Processing: " << geom_object->product()->data().toString() << std::endl;

			const auto& n = geom_object->transformation().data().ccomponents();
			const cgal_placement_t element_transformation(
				n(0, 0), n(0, 1), n(0, 2), n(0, 3),
				n(1, 0), n(1, 1), n(1, 2), n(1, 3),
				n(2, 0), n(2, 1), n(2, 2), n(2, 3));

			boost::optional<Eigen::Vector3d> wall_direction;


			std::list<shape_callback_item*> openings;

			auto p = all_openings.map.equal_range(geom_object->product());
			for (auto it = p.first; it != p.second; ++it) {
				openings.push_back(it->second);
			}

			if (settings.apply_openings_posthoc && geom_object->product()->declaration().is("IfcWall")) {
				auto T2 = timer::measure("wall_axis_handling");
				auto item = geometry_mapper->map(geom_object->product());
				if (item) {
					typedef ifcopenshell::geometry::taxonomy::collection cl;
					typedef ifcopenshell::geometry::taxonomy::loop l;
					typedef ifcopenshell::geometry::taxonomy::edge e;
					typedef ifcopenshell::geometry::taxonomy::point3 p;
					auto edge = (e*)((l*)((cl*)((cl*)item)->children[0])->children[0])->children[0];
					if (edge->basis == nullptr && edge->start.which() == 0 && edge->end.which() == 0) {
						const auto& p0 = boost::get<p>(edge->start);
						const auto& p1 = boost::get<p>(edge->end);
						Eigen::Vector4d P0 = p0.ccomponents().homogeneous();
						Eigen::Vector4d P1 = p1.ccomponents().homogeneous();
						auto V0 = n * P0;
						auto V1 = n * P1;
						std::cout << "Axis " << V0(0) << " " << V0(1) << " " << V0(2) << " -> "
							<< V1(0) << " " << V1(1) << " " << V1(2) << std::endl;
						wall_direction = (V1 - V0).head<3>().normalized();
					}
					T2.stop();
				}
			}

			for (auto& g : geom_object->geometry()) {
				auto s = ((ifcopenshell::geometry::CgalShape*) g.Shape())->shape();
				const auto& m = g.Placement().ccomponents();

				const cgal_placement_t part_transformation(
					m(0, 0), m(0, 1), m(0, 2), m(0, 3),
					m(1, 0), m(1, 1), m(1, 2), m(1, 3),
					m(2, 0), m(2, 1), m(2, 2), m(2, 3));

				// Apply transformation
				for (auto &vertex : vertices(s)) {
					vertex->point() = vertex->point().transform(part_transformation);
				}

				const ifcopenshell::geometry::taxonomy::style* opt_style;
				if (g.hasStyle()) {
					opt_style = &g.Style();
				}

				shape_callback_item item{
					geom_object->product(),
					geom_object->guid(),
					geom_object->type(),
					geom_object->geometry().id(),
					element_transformation,
					s,
					opt_style,
					wall_direction,
					openings
				};

				fn(item);

				std::cout << "Processed: " << geom_object->product()->data().toString() << " part: #" << geom_object->geometry().id() << std::endl;
			}

			std::cout << "Progress: " << context_iterator.progress() << std::endl;

		}

	}

	return 0;
}
