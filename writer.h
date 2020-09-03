#ifndef WRITER_H
#define WRITER_H

#include <ifcgeom/kernels/cgal/CgalKernel.h>

#include <nlohmann/json.hpp>

#include <array>
#include <fstream>

// Abstract writer class that takes triangular facets.
struct abstract_writer {
	std::array<Kernel_::Point_3, 3> points_from_facet(cgal_shape_t::Facet_handle f) {
		return {
				f->facet_begin()->vertex()->point(),
				f->facet_begin()->next()->vertex()->point(),
				f->facet_begin()->next()->next()->vertex()->point()
		};
	}

	std::array<CGAL::Simple_cartesian<double>::Point_3, 3> points_from_facet(CGAL::Polyhedron_3<CGAL::Simple_cartesian<double>>::Facet_handle f) {
		return {
				f->facet_begin()->vertex()->point(),
				f->facet_begin()->next()->vertex()->point(),
				f->facet_begin()->next()->next()->vertex()->point()
		};
	}

	std::array<Kernel_::Point_3, 3> points_from_facet(std::list<cgal_shape_t::Facet_handle>::iterator f) {
		return points_from_facet(*f);
	}
};

// OBJ writer for CGAL facets paired with a style
struct simple_obj_writer : public abstract_writer {
	int group_id = 1;
	int vertex_count = 1;
	std::ofstream obj, mtl;
	ifcopenshell::geometry::taxonomy::colour GRAY;

	simple_obj_writer(const std::string& fn_prefix)
		: obj((fn_prefix + ".obj").c_str())
		, mtl((fn_prefix + ".mtl").c_str())
		, GRAY(0.6, 0.6, 0.6) {
		obj << "mtllib " << fn_prefix << ".mtl\n";
	}

	template <typename It>
	void operator()(const ifcopenshell::geometry::taxonomy::style* style, It begin, It end) {
		const auto& diffuse = *(style ? style->diffuse.get_value_or(GRAY).components : GRAY.components);

		obj << "g group-" << group_id << "\n";
		obj << "usemtl m" << group_id << "\n";
		mtl << "newmtl m" << group_id << "\n";
		mtl << "kd " << diffuse(0) << " " << diffuse(1) << " " << diffuse(2) << "\n";

		group_id++;

		for (auto it = begin; it != end; ++it) {
			auto points = points_from_facet(it);
			for (int i = 0; i < 3; ++i) {
				obj << "v "
					<< points[i].cartesian(0) << " "
					<< points[i].cartesian(1) << " "
					<< points[i].cartesian(2) << "\n";
			}
			obj << "f "
				<< (vertex_count + 0) << " "
				<< (vertex_count + 1) << " "
				<< (vertex_count + 2) << "\n";
			vertex_count += 3;
		}
	}
};

struct city_json_writer : public abstract_writer {
	ifcopenshell::geometry::taxonomy::colour GRAY;

	using json = nlohmann::json;

	std::string filename;

	std::vector<std::array<double, 3>> vertices;
	std::vector<std::vector<std::vector<std::vector<int>>>> boundaries;
	std::vector<std::vector<int>> boundary_materials;

	json materials;

	city_json_writer(const std::string& fn_prefix)
		: filename(fn_prefix + ".json")
		, materials(json::array())
		, GRAY(0.6, 0.6, 0.6) {
		// assumes one solid.
		boundaries.emplace_back();
		boundary_materials.emplace_back();
	}

	template <typename It>
	void operator()(const ifcopenshell::geometry::taxonomy::style* style, It begin, It end) {
		const auto& diffuse = *(style ? style->diffuse.get_value_or(GRAY).components : GRAY.components);

		json material = json::object();
		material["name"] = "material-" + boost::lexical_cast<std::string>(materials.size());
		material["diffuseColor"] = std::array<double, 3>{diffuse(0), diffuse(1), diffuse(2)};
		material["specularColor"] = std::array<double, 3>{0., 0., 0.};
		material["shininess"] = 0.;
		material["isSmooth"] = false;
		materials.push_back(material);

		for (auto it = begin; it != end; ++it) {
			auto points = points_from_facet(it);
			std::vector<int> faces;
			for (int i = 0; i < 3; ++i) {
				faces.push_back(vertices.size());
				vertices.push_back({ {
					CGAL::to_double(points[i].cartesian(0)),
					CGAL::to_double(points[i].cartesian(1)),
					CGAL::to_double(points[i].cartesian(2))
				} });
			}
			boundaries.front().push_back({ faces });
			boundary_materials.front().push_back(materials.size() - 1);
		}
	}

	void finalize() {
		json city;

		city["type"] = "CityJSON";
		city["version"] = "1.0";
		city["extensions"] = json::object();
		city["metadata"]["referenceSystem"] = "urn:ogc:def:crs:EPSG::2355";
		city["vertices"] = vertices;
		city["appearance"]["materials"] = materials;

		auto& building1 = city["CityObjects"]["id-1"];
		building1["type"] = "Building";
		building1["geographicalExtent"] = std::array<double, 6>{0, 0, 0, 1, 1, 1};
		/*
		building1["attributes"]["measuredHeight"] = 22.3;
		building1["attributes"]["roofType"] = "gable";
		building1["attributes"]["owner"] = "Elvis Presley";
		*/

		json geom = json::object();
		geom["type"] = "Solid";
		geom["lod"] = 2;
		geom["boundaries"] = boundaries;
		geom["material"]["diffuse"]["values"] = boundary_materials;
		building1["geometry"].push_back(geom);

		std::ofstream(filename.c_str()) << city;
	}

	~city_json_writer() {
		finalize();
	}
};

#endif
