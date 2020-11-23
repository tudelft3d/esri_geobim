#include "radius_execution_context.h"
#include "writer.h"

#include <ifcconvert/validation_utils.h>

#include <CGAL/exceptions.h>
#include <CGAL/minkowski_sum_3.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <boost/foreach.hpp>

radius_execution_context::radius_execution_context(double r, bool narrower, bool minkowski_triangles, bool no_erosion)
	: radius(r)
	, minkowski_triangles_(minkowski_triangles) 
	, no_erosion_(no_erosion)
	, empty_(true)
{
	{
		auto polycube = ifcopenshell::geometry::utils::create_cube(r);
		padding_cube = ifcopenshell::geometry::utils::create_nef_polyhedron(polycube);
	}

#ifdef ENSURE_2ND_OP_NARROWER
	if (narrower) {
		// double r2 = boost::math::float_advance(r, +5);
		double r2 = r + 1e-7;
		std::cout << r << " -> " << r2 << std::endl;
		auto polycube = ifcopenshell::geometry::utils::create_cube(r2);
		padding_cube_2 = ifcopenshell::geometry::utils::create_nef_polyhedron(polycube);
	} else // -> ...
#endif
		padding_cube_2 = padding_cube;
}

namespace {
	void minkowski_sum_triangles(const CGAL::Polyhedron_3<CGAL::Epick>& poly_triangulated, CGAL::Nef_polyhedron_3<Kernel_>& padding_cube, CGAL::Nef_polyhedron_3<Kernel_>& result) {
		CGAL::Nef_nary_union_3< CGAL::Nef_polyhedron_3<Kernel_> > accum;

		for (auto &face : faces(poly_triangulated)) {

			if (!face->is_triangle()) {
				std::cout << "Warning: non-triangular face!" << std::endl;
				continue;
			}

			CGAL::Polyhedron_3<CGAL::Epick>::Halfedge_around_facet_const_circulator current_halfedge = face->facet_begin();
			CGAL::Point_3<CGAL::Epick> points[3];

			int i = 0;
			do {
				points[i] = current_halfedge->vertex()->point();
				++i;
				++current_halfedge;
			} while (current_halfedge != face->facet_begin());

			double A = std::sqrt(CGAL::to_double(CGAL::Triangle_3<CGAL::Epick>(points[0], points[1], points[2]).squared_area()));
			if (A < (1.e-5 * 1.e-5 * 0.5)) {
				std::cout << "Skipping triangle with area " << A << std::endl;
				continue;
			}

			cgal_shape_t T;
			CGAL::Cartesian_converter<CGAL::Epick, CGAL::Epeck> C;
			T.make_triangle(C(points[0]), C(points[1]), C(points[2]));

			CGAL::Nef_polyhedron_3<Kernel_> Tnef(T);

			CGAL::Nef_polyhedron_3<Kernel_> padded = CGAL::minkowski_sum_3(Tnef, padding_cube);
			accum.add_polyhedron(padded);
		}

		result = accum.get_union();
	}
}


void radius_execution_context::operator()(shape_callback_item& item) {
	if (item.src != previous_src) {
		if (item.geom_reference == previous_geom_ref) {
			std::cout << "Reusing padded geometry for " << item.src->data().toString() << std::endl;
			auto product = per_product_collector.get_union();
			product.transform(last_place.inverse());
			product.transform(item.transformation);
			union_collector.add_polyhedron(product);
			return;
		}
		// per_product_collector = CGAL::Nef_nary_union_3< CGAL::Nef_polyhedron_3<Kernel_> >();
		per_product_collector.clear();
	}

	CGAL::Polyhedron_3<CGAL::Epick> poly_triangulated;
	util::copy::polyhedron(poly_triangulated, item.polyhedron);
	if (!CGAL::Polygon_mesh_processing::triangulate_faces(poly_triangulated)) {
		std::cerr << "unable to triangulate all faces" << std::endl;
		return;
	}

	std::vector<
		std::pair<
		boost::graph_traits<CGAL::Polyhedron_3<CGAL::Epick>>::face_descriptor,
		boost::graph_traits<CGAL::Polyhedron_3<CGAL::Epick>>::face_descriptor>> self_intersections;
	CGAL::Polygon_mesh_processing::self_intersections(poly_triangulated, std::back_inserter(self_intersections));

	previous_src = item.src;
	previous_geom_ref = item.geom_reference;
	last_place = item.transformation;

	CGAL::Nef_polyhedron_3<Kernel_> item_nef, result;
	bool item_nef_succeeded;
	if (!(item_nef_succeeded = item.to_nef_polyhedron(item_nef))) {
		std::cerr << "no nef for product" << std::endl;
	}

	bool result_set = false;
	bool failed = false;

	if (!(minkowski_triangles_ || !item_nef_succeeded || !self_intersections.empty())) {
		item_nef.transform(item.transformation);

		auto T0 = timer::measure("minkowski_sum");
		try {
			CGAL::Nef_polyhedron_3<Kernel_>* item_nef_copy = new CGAL::Nef_polyhedron_3<Kernel_>(item_nef);
			result = CGAL::minkowski_sum_3(*item_nef_copy, padding_cube);
			// So this is funky, we got segfaults in the destructor when exceptions were
			// raised, so we only delete when minkowski (actually the convex_decomposition)
			// succeed. Otherwise, we just have to incur some memory leak.
			// @todo report this to cgal.
			delete item_nef_copy;
			result_set = true;
		} catch (CGAL::Failure_exception&) {
			failed = true;
			std::cerr << "Minkowski on volume failed, retrying with individual triangles" << std::endl;
		}
		T0.stop();
	}

	double max_triangle_area = 0.;

	for (auto &face : faces(poly_triangulated)) {

		if (!face->is_triangle()) {
			std::cout << "Warning: non-triangular face!" << std::endl;
			continue;
		}

		CGAL::Polyhedron_3<CGAL::Epick>::Halfedge_around_facet_const_circulator current_halfedge = face->facet_begin();
		CGAL::Point_3<CGAL::Epick> points[3];

		int i = 0;
		do {
			points[i] = current_halfedge->vertex()->point();
			++i;
			++current_halfedge;
		} while (current_halfedge != face->facet_begin());

		double A = std::sqrt(CGAL::to_double(CGAL::Triangle_3<CGAL::Epick>(points[0], points[1], points[2]).squared_area()));

		if (A > max_triangle_area) {
			max_triangle_area = A;
		}
	}

	if (!result_set && (poly_triangulated.size_of_facets() > 1000 || max_triangle_area < 1.e-5)) {

		if (poly_triangulated.size_of_facets() > 1000) {
			std::cerr << "Too many individual triangles, using bounding box" << std::endl;
		} else {
			std::cerr << "Max triangle area is " << max_triangle_area << ", using bounding box" << std::endl;
		}
		auto bb = CGAL::Polygon_mesh_processing::bbox(item.polyhedron);
		cgal_point_t lower(bb.min(0) - radius, bb.min(1) - radius, bb.min(2) - radius);
		cgal_point_t upper(bb.max(0) + radius, bb.max(1) + radius, bb.max(2) + radius);
		auto bbpl = ifcopenshell::geometry::utils::create_cube(lower, upper);
		result = ifcopenshell::geometry::utils::create_nef_polyhedron(bbpl);

	} else if (!result_set) {
		
		auto T2 = timer::measure("self_intersection_handling");

		if (self_intersections.size()) {
			std::cerr << self_intersections.size() << " self-intersections for product" << std::endl;
		}

		minkowski_sum_triangles(poly_triangulated, padding_cube, result);
		result.transform(item.transformation);

		T2.stop();
	}

	auto T1 = timer::measure("opening_handling");
	if (item.wall_direction && item.openings.size()) {
		static const Eigen::Vector3d Zax(0, 0, 1);
		// @todo derive from model.
		// @todo since many walls will be parallel we can cache these polyhedrons
		static const double EPS = 1.e-5;

		auto Yax = Zax.cross(*item.wall_direction).normalized();

		auto x0 = *item.wall_direction * -radius;
		auto x1 = *item.wall_direction * +radius;

		auto y0 = Yax * -(radius + EPS);
		auto y1 = Yax * +(radius + EPS);

		Kernel_::Point_3 X0(x0(0), x0(1), x0(2));
		Kernel_::Point_3 X1(x1(0), x1(1), x1(2));

		Kernel_::Point_3 Y0(y0(0), y0(1), y0(2));
		Kernel_::Point_3 Y1(y1(0), y1(1), y1(2));

		Kernel_::Point_3 Z0(0, 0, -radius);
		Kernel_::Point_3 Z1(0, 0, +radius);

		// CGAL::Nef_polyhedron_3<Kernel_> X(CGAL::Segment_3<Kernel_>(X0, X1));
		CGAL::Nef_polyhedron_3<Kernel_> Y(CGAL::Segment_3<Kernel_>(Y0, Y1));
		// CGAL::Nef_polyhedron_3<Kernel_> Z(CGAL::Segment_3<Kernel_>(Z0, Z1));
		// auto ZX = CGAL::minkowski_sum_3(X, Z);

		// manual minkowski sum...
		CGAL::Polyhedron_3<Kernel_> zx;
		std::list<cgal_point_t> zx_points{ {
			CGAL::ORIGIN + ((X0 - CGAL::ORIGIN) + (Z0 - CGAL::ORIGIN)),
			CGAL::ORIGIN + ((X1 - CGAL::ORIGIN) + (Z0 - CGAL::ORIGIN)),
			CGAL::ORIGIN + ((X1 - CGAL::ORIGIN) + (Z1 - CGAL::ORIGIN)),
			CGAL::ORIGIN + ((X0 - CGAL::ORIGIN) + (Z1 - CGAL::ORIGIN))
		} };
		std::vector<std::vector<int>> zx_idxs{ {{{0,1,2,3}}} };
		util::PolyFromMesh<cgal_shape_t::HDS> m(zx_points, zx_idxs);
		zx.delegate(m);
		CGAL::Nef_polyhedron_3<Kernel_> ZX(zx);

		CGAL::Nef_nary_union_3< CGAL::Nef_polyhedron_3<Kernel_> > opening_union;
		for (auto& op : item.openings) {

			const auto& xdir = *item.wall_direction;
			double min_dot = +std::numeric_limits<double>::infinity();
			double max_dot = -std::numeric_limits<double>::infinity();
			double min_z = +std::numeric_limits<double>::infinity();
			double max_z = -std::numeric_limits<double>::infinity();

			for (const auto& v : vertices(op->polyhedron)) {
				auto p = v->point();
				p = p.transform(op->transformation);
				Eigen::Vector3d vv(
					CGAL::to_double(p.cartesian(0)),
					CGAL::to_double(p.cartesian(1)),
					CGAL::to_double(p.cartesian(2))
				);
				double d = xdir.dot(vv);
				if (d < min_dot) {
					min_dot = d;
				}
				if (d > max_dot) {
					max_dot = d;
				}
				if (vv.z() < min_z) {
					min_z = vv.z();
				}
				if (vv.z() > max_z) {
					max_z = vv.z();
				}
			}

			// These are basically workarounds for a bug in
			// nary_union<T> which is also used on minkowski of concave operands.
			// It segaults on getting front() of an empty queue.
			// with the patch https://patch-diff.githubusercontent.com/raw/CGAL/cgal/pull/4768.patch
			// (applied now by the IfcOpenShell build script)
			// these fixes are not necessary anymore.
			if ((max_dot - min_dot) < radius * 2) {
				std::cerr << "Opening too narrow to have effect after incorporating radius, skipping" << std::endl;
				continue;
			}

			if ((max_z - min_z) < radius * 2) {
				std::cerr << "Opening too narrow to have effect after incorporating radius, skipping" << std::endl;
				continue;
			}

			auto bounds = create_bounding_box(op->polyhedron);
			CGAL::Nef_polyhedron_3<Kernel_> opening_nef;
			if (!op->to_nef_polyhedron(opening_nef)) {
				std::cerr << "no nef for opening" << std::endl;
				continue;
			}
			opening_nef.transform(op->transformation);
			bounds.transform(op->transformation);

			auto temp = bounds - opening_nef;
			temp = CGAL::minkowski_sum_3(ZX, temp);
			temp = bounds - temp;
			temp = CGAL::minkowski_sum_3(temp, Y);

#ifdef GEOBIM_DEBUG
			simple_obj_writer x("opening-" + op->id);
			x(nullptr, item.polyhedron.facets_begin(), item.polyhedron.facets_end());
			x(nullptr, op->polyhedron.facets_begin(), op->polyhedron.facets_end());
			CGAL::Polyhedron_3<Kernel_> temp2;
			temp.convert_to_polyhedron(temp2);
			x(nullptr, temp2.facets_begin(), temp2.facets_end());
#endif

			result -= temp;
		}
	}
	T1.stop();

	union_collector.add_polyhedron(result);
	per_product_collector.add_polyhedron(result);

	empty_ = false;
}

namespace {
	double bbox_diagonal(const CGAL::Bbox_3& b) {
		return std::sqrt(
			(b.xmax() - b.xmin()) * (b.xmax() - b.xmin()) +
			(b.ymax() - b.ymin()) * (b.ymax() - b.ymin()) +
			(b.zmax() - b.zmin()) * (b.zmax() - b.zmin())
		);
	}
}


namespace {
	// Return pair in map with largest value
	template<typename T, typename U>
	std::pair<T, U> map_max_value(const std::map<T, U>& x) {
		return *std::max_element(x.begin(), x.end(), [](const std::pair<T, U>& p1, const std::pair<T, U>& p2) {
			return p1.second < p2.second;
		});
	}

	template<size_t N, typename T, typename U>
	std::pair<T, U> nth_largest_value_from_map(const std::map<T, U>& x) {
		std::vector<std::pair<U, T>> vec;
		for (auto& y : x) {
			vec.push_back({ y.second, y.first });
		}
		std::sort(vec.begin(), vec.end());
		std::reverse(vec.begin(), vec.end());
		if (N > vec.size()) {
			throw std::runtime_error("not enough components");
		}
		auto it = vec.begin() + N;
		return { it->second, it->first };
	}

}

// Extract the exterior component of a CGAL Polyhedron
void radius_execution_context::extract_in_place(cgal_shape_t& input, extract_component component) const {
	if (input.facets_begin() == input.facets_end()) {
		throw std::runtime_error("Empty input operand to extract()");
	}

	// split_connected_components() is introduced in CGAL 5 :(

	typedef boost::graph_traits<cgal_shape_t>::face_descriptor face_descriptor;
	typedef boost::graph_traits<cgal_shape_t>::vertex_descriptor vertex_descriptor;
	typedef boost::graph_traits<cgal_shape_t>::vertex_iterator vertex_iterator;

	boost::property_map<cgal_shape_t, boost::face_external_index_t>::type fim
		= get(boost::face_external_index, input);

	boost::vector_property_map<size_t,
		boost::property_map<cgal_shape_t, boost::face_external_index_t>::type>
		fsm(fim);

	auto ffim = CGAL::Polygon_mesh_processing::parameters::face_index_map(fim);

	CGAL::Polygon_mesh_processing::connected_components(
		input, 
		fsm,
		ffim);

	face_descriptor largest_component_facet;
	std::map<size_t, size_t> component_sizes;
	std::map<size_t, double> component_areas;

	BOOST_FOREACH(face_descriptor f, faces(input)) {
		auto idx = fsm[f];
		component_sizes[idx] ++;
		component_areas[idx] += CGAL::to_double(CGAL::Polygon_mesh_processing::area(std::vector<face_descriptor>{f}, input));
		if (component_sizes.rbegin()->first == idx) {
			largest_component_facet = f;
		}
	}

	for (auto& p : component_sizes) {
		std::cout << "component " << p.first << " has " << p.second << " and area " << component_areas[p.first] << std::endl;
	}

	typedef std::map<vertex_descriptor, std::size_t>   Internal_vertex_map;
	typedef boost::associative_property_map<Internal_vertex_map>   Vertex_index_map;
	Internal_vertex_map internal_vertex_index_map;
	Vertex_index_map vertex_index_map(internal_vertex_index_map);
	vertex_iterator vb, ve;
	std::size_t counter = 0;
	for (boost::tie(vb, ve) = vertices(input); vb != ve; ++vb, ++counter) {
		put(vertex_index_map, *vb, counter);
	}

	std::vector<size_t> components_to_keep;

	if (component == LARGEST_AREA) {
		components_to_keep.push_back(map_max_value(component_areas).first);
	} else if (component == SECOND_LARGEST_AREA) {
		components_to_keep.push_back(nth_largest_value_from_map<1>(component_areas).first);
	}

	CGAL::Polygon_mesh_processing::keep_connected_components(
		input,
		components_to_keep,
		fsm,
		CGAL::Polygon_mesh_processing::parameters::vertex_index_map(vertex_index_map));
}

cgal_shape_t radius_execution_context::extract(const cgal_shape_t& input, extract_component component) const {
	auto input_copy = input;
	extract_in_place(input_copy, component);
	return input_copy;
}

// Create a bounding box (six-faced Nef poly) around a CGAL Polyhedron

CGAL::Nef_polyhedron_3<Kernel_> radius_execution_context::create_bounding_box(const cgal_shape_t & input) const {
	// @todo we can probably use
	// CGAL::Nef_polyhedron_3<Kernel_> nef(CGAL::Nef_polyhedron_3<Kernel_>::COMPLETE)
	// Implementation detail: there is always an implicit box around the Nef, even when open or closed.
	// Never mind: The Minkowski sum cannot operate on unbounded inputs...

	// Create the complement of the Nef by subtracting from its bounding box,
	// see: https://github.com/tudelft3d/ifc2citygml/blob/master/off2citygml/Minkowski.cpp#L23
	auto bounding_box = CGAL::Polygon_mesh_processing::bbox(input);
	Kernel_::Point_3 bbmin(bounding_box.xmin(), bounding_box.ymin(), bounding_box.zmin());
	Kernel_::Point_3 bbmax(bounding_box.xmax(), bounding_box.ymax(), bounding_box.zmax());
	Kernel_::Vector_3 d(radius, radius, radius);
	bbmin = CGAL::ORIGIN + ((bbmin - CGAL::ORIGIN) - d);
	bbmax = CGAL::ORIGIN + ((bbmax - CGAL::ORIGIN) + d);
	cgal_shape_t poly_box = ifcopenshell::geometry::utils::create_cube(bbmin, bbmax);
	return ifcopenshell::geometry::utils::create_nef_polyhedron(poly_box);
}

// Completes the boolean union, extracts exterior and erodes padding radius

void radius_execution_context::finalize() {
	{
		auto T = timer::measure("nef_boolean_union");
		// @todo spatial sorting?
		boolean_result = union_collector.get_union();
		T.stop();
	}

	if (no_erosion_) {
		exterior = boolean_result;
	} else {

		auto T2 = timer::measure("result_nef_processing");
		polyhedron_exterior = ifcopenshell::geometry::utils::create_polyhedron(boolean_result);

		{
			simple_obj_writer tmp_debug("debug-after-boolean");
			tmp_debug(nullptr, polyhedron_exterior.facets_begin(), polyhedron_exterior.facets_end());
		}

		// @todo Wasteful: remove interior on Nef?
		extract_in_place(polyhedron_exterior, LARGEST_AREA);

		{
			simple_obj_writer tmp_debug("debug-exterior");
			tmp_debug(nullptr, polyhedron_exterior.facets_begin(), polyhedron_exterior.facets_end());
		}

#if 0

		exterior = ifcopenshell::geometry::utils::create_nef_polyhedron(polyhedron_exterior);
		bounding_box = create_bounding_box(polyhedron_exterior);

		complement = bounding_box - exterior;
		complement.extract_regularization();
		// @nb padding cube is potentially slightly larger to result in a thinner result
		// then another radius for comparison.
		complement_padded = CGAL::minkowski_sum_3(complement, padding_cube_2);
		complement_padded.extract_regularization();
		T2.stop();

		{
			auto T = timer::measure("result_nef_to_poly");
#if 0
			// @todo I imagine this operation is costly, we can also convert the padded complement to
			// polyhedron, and remove the connected component that belongs to the bbox, then reverse
			// the remaining poly to point to the interior?
			exterior -= complement_padded;
#else
			// Rougly twice as fast as the complexity is half (box complexity is negligable).

			// Re above: extracting the interior shell did not prove to be reliable even with
			// the undocumented function convert_inner_shell_to_polyhedron(). Therefore we
			// subtract from the padded box as that will have lower complexity than above.
			// Mark_bounded_volumes on the completement also did not work.
			exterior = bounding_box - complement_padded;
#endif
			T.stop();
		}

		exterior.extract_regularization();

		if (exterior.is_simple()) {
			auto T1 = timer::measure("result_nef_to_poly");
			polyhedron_exterior = ifcopenshell::geometry::utils::create_polyhedron(exterior);
			T1.stop();

			auto vol = CGAL::Polygon_mesh_processing::volume(polyhedron_exterior);
			std::cout << "Volume with radius " << radius << " is " << vol << std::endl;
		} else {
			CGAL::convert_nef_polyhedron_to_polygon_mesh(exterior, polyhedron_exterior);
			std::cout << "Result with radius " << radius << " is not manifold" << std::endl;
		}

#else

		CGAL::Polyhedron_3<CGAL::Epick> poly_triangulated;
		util::copy::polyhedron(poly_triangulated, polyhedron_exterior);
		if (!CGAL::Polygon_mesh_processing::triangulate_faces(poly_triangulated)) {
			std::cerr << "unable to triangulate all faces" << std::endl;
			return;
		}
		minkowski_sum_triangles(poly_triangulated, padding_cube, exterior);
		if (exterior.is_simple()) {
			polyhedron_exterior = ifcopenshell::geometry::utils::create_polyhedron(exterior);
		} else {
			CGAL::convert_nef_polyhedron_to_polygon_mesh(exterior, polyhedron_exterior);
		}

		extract_in_place(polyhedron_exterior, SECOND_LARGEST_AREA);
#endif

	}

	std::cout << "exterior poly num facets: " << polyhedron_exterior.size_of_facets() << std::endl;
}

