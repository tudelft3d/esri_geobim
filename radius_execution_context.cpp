#include "radius_execution_context.h"
#include "writer.h"

#include <ifcconvert/validation_utils.h>

#include <CGAL/exceptions.h>
#include <CGAL/minkowski_sum_3.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>

#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <boost/foreach.hpp>

#include <CGAL/boost/graph/helpers.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/cluster_point_set.h>
#include <CGAL/compute_average_spacing.h>

#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

static bool ENSURE_2ND_OP_NARROWER = true;
static double MAKE_OP2_NARROWER = ENSURE_2ND_OP_NARROWER ? -2e-7 : 0.0;

namespace {
	template <typename K>
	struct average_point {
		CGAL::Vector_3<K> accum;
		size_t n = 0;

		void add(const CGAL::Point_3<K>& p) {
			accum += p - CGAL::ORIGIN;
			n += 1;
		}

		operator CGAL::Point_3<K>() const {
			return CGAL::ORIGIN + accum / n;
		}
	};

	// @todo this should not be based on euclidean distance, but rather distance over edge/face as it will now close holes which is not the intention
	template <typename T>
	bool cluster_vertices(T& s, double r) {
		typedef typename T::Traits::Kernel K;
		typedef CGAL::Point_3<K> P;

		std::vector<P> points;
		std::vector<std::vector<size_t>> indices;
		CGAL::Polygon_mesh_processing::polygon_mesh_to_polygon_soup(s, points, indices);

		std::map<P, size_t> clusters;
		boost::associative_property_map<std::map<P, size_t>> clusters_map(clusters);

		std::cout << "Clustering with radius " << r << std::endl;

		int n = CGAL::cluster_point_set(points, clusters_map, CGAL::parameters::neighbor_radius(r));
		std::cout << n << " clusters" << std::endl;

		std::vector<average_point<K>> new_points_accum(n);
		for (auto& p : clusters) {
			new_points_accum[p.second].add(p.first);
		}

		std::vector<P> new_points;
		for (auto& p : new_points_accum) {
			new_points.push_back(p);
		}

		std::vector<std::vector<size_t>> new_indices;

		for (auto& x : indices) {
			std::vector<size_t> transformed;
			std::transform(x.begin(), x.end(), std::back_inserter(transformed), [&clusters, &points](size_t i) {
				return clusters.find(points[i])->second;
			});
			std::set<size_t> transformed_unique(transformed.begin(), transformed.end());
			if (transformed_unique.size() == transformed.size()) {
				new_indices.push_back(transformed);
			}
		}

		std::map<std::set<size_t>, size_t> triangle_use;
		for (auto& x : new_indices) {
			std::set<size_t> s(x.begin(), x.end());
			triangle_use[s] ++;
		}

		auto it = new_indices.end();
		while (it > new_indices.begin())
		{
			it--;
			std::set<size_t> s(it->begin(), it->end());
			if (triangle_use[s] == 2) {
				it = new_indices.erase(it);
				std::cerr << "erasing" << std::endl;
			}
		}

		std::ofstream fs2("points.txt");
		fs2 << "[\n";
		for (auto& p : new_points) {
			fs2 << "[" << p.cartesian(0) << "," << p.cartesian(1) << "," << p.cartesian(2) << "],\n";
		}
		fs2.seekp(-3, std::ios_base::cur);
		fs2 << "]";
		fs2.close();
		
		std::ofstream fs("debug.txt");

		std::map<std::pair<size_t, size_t>, size_t> edge_use;
		auto add_use = [&edge_use](size_t a, size_t b) {
			if (a > b) std::swap(a, b);
			edge_use[{a, b}]++;
		};
		fs << "[\n";
		for (auto& tri : new_indices) {
			fs << "[" << tri[0] << "," << tri[1] << "," << tri[2] << "],\n";
			add_use(tri[0], tri[1]);
			add_use(tri[1], tri[2]);
			add_use(tri[2], tri[0]);
		}
		fs.seekp(-3, std::ios_base::cur);
		fs << "]";
		fs.close();

		for (auto& p : edge_use) {
			if (p.second != 2) {
				std::cerr << "non-manifold: " << p.first.first << " " << p.first.second << std::endl;
				return false;
				std::cerr << "attach debugger and press key" << std::endl;
				std::cin.get();
			}
		}

		std::cerr << "removed points: " << CGAL::Polygon_mesh_processing::remove_isolated_points_in_polygon_soup(new_points, new_indices) << std::endl;

		auto v = CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(new_indices);
		std::cerr << "valid: " << v << std::endl;
		if (!v) {
			CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(new_indices);
		}

		T new_poly;
		CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(new_points, new_indices, new_poly);

		/*
		for (auto it = new_poly.vertices_begin(); it != new_poly.vertices_begin(); ++it) {
			if (it != it->halfedge()->vertex()) {
				auto jt = std::find(new_points.begin(), new_points.end(), it->point());
				auto NN = std::distance(new_points.begin(), jt);
				std::cin.get();
				std::cerr << NN << std::endl;
			}
		}
		*/

		s = new_poly;
		
		return true;
	}
}

radius_execution_context::radius_execution_context(const std::string& r, radius_settings rs)
	: settings_(rs)
	, radius_str(r)
	, radius(boost::lexical_cast<double>(r))
	, minkowski_triangles_(rs.get(radius_settings::MINKOWSKI_TRIANGLES))
	, no_erosion_(rs.get(radius_settings::NO_EROSION))
	, empty_(false) // no longer relevant, bug fixed
{
	/*
	padding_volume = construct_padding_volume_(radius);
	if (ENSURE_2ND_OP_NARROWER && rs.get(radius_settings::NARROWER)) {
		double r2 = radius + 1e-7;
		padding_volume_2 = construct_padding_volume_(r2);
	}
	else {
		padding_volume_2 = padding_volume;
	}
	*/
}

CGAL::Nef_polyhedron_3<Kernel_> radius_execution_context::construct_padding_volume_(const boost::optional<double>& R) {
	double radius = R.get_value_or(this->radius);

	if (settings_.get(radius_settings::SPHERE)) {
		cgal_shape_t ico;
			
		CGAL::make_icosahedron(ico, cgal_point_t(0, 0, 0), 1.0);

		double ml = std::numeric_limits<double>::infinity();
			
		// Take the edge centers and find minimal distance from origin.
		// Or use vertex position
		for (auto e : edges(ico)) {
			auto v1 = e.halfedge()->vertex();
			auto v2 = e.opposite().halfedge()->vertex();

			double v1x = CGAL::to_double(v1->point().cartesian(0));
			double v1y = CGAL::to_double(v1->point().cartesian(1));
			double v1z = CGAL::to_double(v1->point().cartesian(2));

			double v2x = CGAL::to_double(v2->point().cartesian(0));
			double v2y = CGAL::to_double(v2->point().cartesian(1));
			double v2z = CGAL::to_double(v2->point().cartesian(2));

#ifdef ICO_EDGE_CENTRES
			double vx = (v1x + v2x) / 2.;
			double vy = (v1y + v2y) / 2.;
			double vz = (v1z + v2z) / 2.;

			double l = std::sqrt(vx*vx + vy * vy + vz * vz);
			if (l < ml) {
				ml = l;
			}
#else
			double l = std::sqrt(v1x*v1x + v1y * v1y + v1z * v1z);
			if (l < ml) {
				ml = l;
			}
			l = std::sqrt(v2x*v2x + v2y * v2y + v2z * v2z);
			if (l < ml) {
				ml = l;
			}
#endif
		}

		// Divide the coordinates with the miminal distance
		for (auto& v : vertices(ico)) {
			v->point() = CGAL::ORIGIN + ((v->point() - CGAL::ORIGIN) * (radius / ml));
		}

		ico_edge_length = 10.;
		// Now compute ico edge length, we use it later as a treshold for simplification
		for (auto e : edges(ico)) {
			auto v1 = e.halfedge()->vertex();
			auto v2 = e.opposite().halfedge()->vertex();

			double v1x = CGAL::to_double(v1->point().cartesian(0));
			double v1y = CGAL::to_double(v1->point().cartesian(1));
			double v1z = CGAL::to_double(v1->point().cartesian(2));

			double v2x = CGAL::to_double(v2->point().cartesian(0));
			double v2y = CGAL::to_double(v2->point().cartesian(1));
			double v2z = CGAL::to_double(v2->point().cartesian(2));

			double vx = (v1x - v2x);
			double vy = (v1y - v2y);
			double vz = (v1z - v2z);

			double l = std::sqrt(vx*vx + vy * vy + vz * vz);
			if (l < ico_edge_length) {
				ico_edge_length = l;
			}
		}

		return ifcopenshell::geometry::utils::create_nef_polyhedron(ico);
	}
	else {
		auto polycube = ifcopenshell::geometry::utils::create_cube(radius);
		return ifcopenshell::geometry::utils::create_nef_polyhedron(polycube);
	}
}

#include <thread>
#include <future>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

namespace SMS = CGAL::Surface_mesh_simplification;

#ifdef _MSC_VER
#include "windows.h"
#include "psapi.h"
void print_mem_usage() {
	PROCESS_MEMORY_COUNTERS_EX pmc;
	GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
	std::cout << "memory: " << pmc.PrivateUsage << std::endl;
}
#else
void print_mem_usage() {}
#endif


namespace {
	template <size_t N, typename T>
	class queue_shortener {
		T t;
		size_t x;

	public:
		queue_shortener() : x(0) {}

		void add_polyhedron(const CGAL::Nef_polyhedron_3<Kernel_>& u) {
			t.add_polyhedron(u);
			++x;
// does not seem to reduce memory footprint
#if 0
			
			if (x == N) {
				std::cout << "before" << std::endl;
				print_mem_usage();
				auto v = t.get_union();
				v.extract_regularization();
				t = T();
				t.add_polyhedron(v);
				x = 0;
				std::cout << "after" << std::endl;
				print_mem_usage();
			}
#endif
		}

		CGAL::Nef_polyhedron_3<Kernel_> get_union() {
			return t.get_union();
		}
	};

	template <typename Kernel>
	void minkowski_sum_triangles_array_impl(std::vector<std::array<CGAL::Point_3<Kernel>, 3>>& triangles, CGAL::Nef_polyhedron_3<Kernel_>& padding_volume, CGAL::Nef_polyhedron_3<Kernel_>& result) {
		// queue_shortener<30, CGAL::Nef_nary_union_3< CGAL::Nef_polyhedron_3<Kernel_> > > accum;
		CGAL::Nef_nary_union_3< CGAL::Nef_polyhedron_3<Kernel_> > accum;
		for (auto& points : triangles) {
			double A = std::sqrt(CGAL::to_double(CGAL::Triangle_3<Kernel>(points[0], points[1], points[2]).squared_area()));
			if (A < (1.e-5 * 1.e-5 * 0.5)) {
				continue;
			}

			cgal_shape_t T;
			CGAL::Cartesian_converter<Kernel, CGAL::Epeck> C;
			T.make_triangle(C(points[0]), C(points[1]), C(points[2]));

			CGAL::Nef_polyhedron_3<Kernel_> Tnef(T);

			CGAL::Nef_polyhedron_3<Kernel_> padded = CGAL::minkowski_sum_3(Tnef, padding_volume);
			accum.add_polyhedron(padded);
		}
		result = accum.get_union();
	}

	template <typename Poly>
	void minkowski_sum_triangles_single_threaded(typename Poly::Facet_const_iterator begin, typename Poly::Facet_const_iterator end, CGAL::Nef_polyhedron_3<Kernel_>& padding_volume, CGAL::Nef_polyhedron_3<Kernel_>& result) {
		// queue_shortener<30, CGAL::Nef_nary_union_3< CGAL::Nef_polyhedron_3<Kernel_> > > accum;
		CGAL::Nef_nary_union_3< CGAL::Nef_polyhedron_3<Kernel_> > accum;

		size_t num = std::distance(begin, end);

		std::cout << "\n";

		for (auto face = begin; face != end; ++face) {

			if (!face->is_triangle()) {
				std::cout << "Warning: non-triangular face!" << std::endl;
				continue;
			}

			typename Poly::Halfedge_around_facet_const_circulator current_halfedge = face->facet_begin();
			CGAL::Point_3<typename Poly::Traits::Kernel> points[3];

			int i = 0;
			do {
				points[i] = current_halfedge->vertex()->point();
				++i;
				++current_halfedge;
			} while (current_halfedge != face->facet_begin());

			double A = std::sqrt(CGAL::to_double(CGAL::Triangle_3<typename Poly::Traits::Kernel>(points[0], points[1], points[2]).squared_area()));
			if (A < (1.e-5 * 1.e-5 * 0.5)) {
				std::cout << "Skipping triangle with area " << A << std::endl;
				continue;
			}

			cgal_shape_t T;
			CGAL::Cartesian_converter<typename Poly::Traits::Kernel, CGAL::Epeck> C;
			T.make_triangle(C(points[0]), C(points[1]), C(points[2]));

			CGAL::Nef_polyhedron_3<Kernel_> Tnef(T);

			CGAL::Nef_polyhedron_3<Kernel_> padded = CGAL::minkowski_sum_3(Tnef, padding_volume);
			accum.add_polyhedron(padded);

			auto n = std::distance(begin, face);
			if (n % 100) {
				std::cout << "\r" << (n * 100 / num) << "%";
				std::cout.flush();
			}
		}

		std::cout << "\n";

		result = accum.get_union();
	}
	
	void minkowski_sum_triangles_double_multithreaded(const CGAL::Polyhedron_3<CGAL::Epick>& poly_triangulated_epick, CGAL::Nef_polyhedron_3<Kernel_>& padding_volume, CGAL::Nef_polyhedron_3<Kernel_>& result) {
		// We need to copy to non-filtered kernel for multi threading
		// We're using double so that we can sneak in SMS, it would fail otherwise on missing sqrt()
		typedef CGAL::Simple_cartesian<double> TriangleKernel;
		CGAL::Polyhedron_3<TriangleKernel> poly_triangulated;
		util::copy::polyhedron(poly_triangulated, poly_triangulated_epick);

		double stop_ratio = 0.1;
		CGAL::Surface_mesh_simplification::Edge_length_stop_predicate<double> stop(1.e-3);
		int r = SMS::edge_collapse(poly_triangulated, stop, 
			CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index, poly_triangulated))
			.halfedge_index_map(get(CGAL::halfedge_external_index, poly_triangulated))
			.get_cost(SMS::Edge_length_cost<TriangleKernel>()));

		std::cout << "Removed " << r << " edges" << std::endl;

		size_t n_threads = std::thread::hardware_concurrency();
		size_t n_facets = poly_triangulated.size_of_facets();
		size_t facets_per_thread = n_facets / n_threads;

		std::vector< std::future<void> > threadpool;
		std::vector< CGAL::Nef_polyhedron_3<Kernel_> > results(n_threads);
		std::vector< CGAL::Nef_polyhedron_3<Kernel_> > cubes(n_threads);

		for (size_t i = 0; i < n_threads; ++i) {
			// even the cubes have to be reconstructed to avoid race conditions
			auto polycube = ifcopenshell::geometry::utils::create_cube(0.05);
			cubes[i] = ifcopenshell::geometry::utils::create_nef_polyhedron(polycube);
		}

		std::vector< std::vector<std::array<CGAL::Point_3<TriangleKernel>, 3>> > triangles(n_threads);
		threadpool.reserve(n_threads);

		CGAL::Polyhedron_3<TriangleKernel>::Facet_const_iterator begin, end;
		begin = poly_triangulated.facets_begin();

		for (size_t i = 0; i < n_threads; ++i) {
			if (i == (n_threads - 1)) {
				end = poly_triangulated.facets_end();
			} else {
				end = begin;
				std::advance(end, facets_per_thread);
			}

			for (auto it = begin; it != end; ++it) {
				triangles[i].emplace_back();
				CGAL::Polyhedron_3<TriangleKernel>::Halfedge_around_facet_const_circulator current_halfedge = it->facet_begin();
				int ii = 0;
				do {
					const auto& p = current_halfedge->vertex()->point();
					/*
					CGAL::Gmpq px = p.cartesian(0).mpq();
					CGAL::Gmpq py = p.cartesian(1).mpq();
					CGAL::Gmpq pz = p.cartesian(2).mpq();
					*/
					
					double px = p.cartesian(0);
					double py = p.cartesian(1);
					double pz = p.cartesian(2);
					
					triangles[i].back()[ii] = CGAL::Point_3<TriangleKernel>(px, py, pz);
					++ii;
					++current_halfedge;
				} while (current_halfedge != it->facet_begin());
			}
			
			std::future<void> fu = std::async(
				std::launch::async, 
				minkowski_sum_triangles_array_impl<TriangleKernel>,
				std::ref(triangles[i]), 
				std::ref(cubes[i]), 
				std::ref(results[i])
			);
			threadpool.emplace_back(std::move(fu));
			end = begin;
		}

		for (std::future<void>& fu : threadpool) {
			fu.get();
		}

		CGAL::Nef_nary_union_3< CGAL::Nef_polyhedron_3<Kernel_> > thread_join;
		for (auto& r : results) {
			thread_join.add_polyhedron(r);
		}

		result = thread_join.get_union();
	}
}


// Create a bounding box (six-faced Nef poly) around a CGAL Polyhedron
CGAL::Nef_polyhedron_3<Kernel_> create_bounding_box(const cgal_shape_t & input, double radius) {
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


class process_shape_item {
	double radius;
	bool minkowski_triangles_, threaded_;
public:

	process_shape_item(double r, bool mintri, bool threaded)
		: radius(r)
		, minkowski_triangles_(mintri)
		, threaded_(threaded)
	{}
	
#if 1
	void operator()(shape_callback_item* item_ptr, CGAL::Nef_polyhedron_3<Kernel_>* result_ptr, CGAL::Nef_polyhedron_3<Kernel_>* padding_volume_, CGAL::Nef_polyhedron_3<Kernel_>* padding_volume_2_) {
		auto& item = *item_ptr;
		auto& result = *result_ptr;
		auto& padding_volume = *padding_volume_;
		auto& padding_volume_2 = *padding_volume_2_;
#else
	void operator()(shape_callback_item item, CGAL::Nef_polyhedron_3<Kernel_>& result) {
#endif

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

		CGAL::Nef_polyhedron_3<Kernel_> item_nef;
		bool item_nef_succeeded = false;
		if (self_intersections.empty()) {
			if (!(item_nef_succeeded = item.to_nef_polyhedron(item_nef, threaded_))) {
				std::cerr << "no nef for product" << std::endl;
			}
		}
		else {
			std::cerr << "self intersections, not trying to convert to Nef" << std::endl;
		}

		bool result_set = false;
		bool failed = false;

		if (!(minkowski_triangles_ || !item_nef_succeeded || !self_intersections.empty())) {
			item_nef.transform(item.transformation);

			auto T0 = timer::measure("minkowski_sum");
			try {
				CGAL::Nef_polyhedron_3<Kernel_>* item_nef_copy = new CGAL::Nef_polyhedron_3<Kernel_>(item_nef);
				result = CGAL::minkowski_sum_3(*item_nef_copy, padding_volume);
				// So this is funky, we got segfaults in the destructor when exceptions were
				// raised, so we only delete when minkowski (actually the convex_decomposition)
				// succeed. Otherwise, we just have to incur some memory leak.
				// @todo report this to cgal.
				// Still an issue on 5.2. valgrind reports an error as well.
				
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
			}
			else {
				std::cerr << "Max triangle area is " << max_triangle_area << ", using bounding box" << std::endl;
			}
			auto bb = CGAL::Polygon_mesh_processing::bbox(item.polyhedron);
			cgal_point_t lower(bb.min(0) - radius, bb.min(1) - radius, bb.min(2) - radius);
			cgal_point_t upper(bb.max(0) + radius, bb.max(1) + radius, bb.max(2) + radius);
			auto bbpl = ifcopenshell::geometry::utils::create_cube(lower, upper);
			result = ifcopenshell::geometry::utils::create_nef_polyhedron(bbpl);

		}
		else if (!result_set) {

			auto T2 = timer::measure("self_intersection_handling");

			if (self_intersections.size()) {
				std::cerr << self_intersections.size() << " self-intersections for product" << std::endl;
			}

			minkowski_sum_triangles_single_threaded<CGAL::Polyhedron_3<CGAL::Epick>>(
				poly_triangulated.facets_begin(),
				poly_triangulated.facets_end(),
				padding_volume_2, result
				);
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

				auto bounds = create_bounding_box(op->polyhedron, radius);
				CGAL::Nef_polyhedron_3<Kernel_> opening_nef;
				if (!op->to_nef_polyhedron(opening_nef, threaded_)) {
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

		result.extract_regularization();
		
		T1.stop();
	}
};

void radius_execution_context::operator()(shape_callback_item* item) {
	auto it = first_product_for_geom_id.find(item->geom_reference);
	if (it != first_product_for_geom_id.end()) {
		if (it->second != item->src) {
			if (reused_products.find(item->src) == reused_products.end()) {
				std::cout << "Reused " << it->second << std::endl;
				reused_products.insert({ item->src, {
					it->second,
					placements.find(it->second)->second,
					item->transformation}
					});
			}
			return;
		}
	}
	else {
		first_product_for_geom_id.insert(it, { item->geom_reference, item->src });
		// placements only need to be inserted once.
		placements[item->src] = item->transformation.inverse();
	}

	product_geometries[item->src].emplace_back();
	auto result_nef = &product_geometries[item->src].back();
	process_shape_item* task = new process_shape_item(radius, minkowski_triangles_, (bool) threads_);
	
	padding_volumes_.push_back({construct_padding_volume_(), construct_padding_volume_()});
	auto& pp = padding_volumes_.back();
		
	if (!threads_) {
		(*task)(item, result_nef, &pp.first, &pp.second);
	} else {
		/*
		bool placed = false;
		while (!placed) {
			for (auto& fu : threadpool_) {
				if (!fu.valid()) {					
					fu = std::async(std::launch::async, std::ref(*task), item, result_nef);
					placed = true;
					break;
				}
			}
			if (!placed) {
				for (auto& fu : threadpool_) {
					if (fu.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
						try {
							fu.get();
						}
						catch (std::exception& e) {
							std::cerr << e.what() << std::endl;
						}
						catch (...) {
							std::cerr << "unkown error" << std::endl;
						}
						fu = std::async(std::launch::async, std::ref(*task), item, result_nef);
						placed = true;
						break;
					}
				}
			}
		}
		*/
		
		while (threadpool_.size() == threads_) {
			for (size_t i = 0; i < threadpool_.size(); ++i) {
				std::future<void> &fu = threadpool_[i];
				std::future_status status = fu.wait_for(std::chrono::seconds(0));
				if (status == std::future_status::ready) {
					fu.get();
					
					std::swap(threadpool_[i], threadpool_.back());
					threadpool_.pop_back();
				}
			}
		}
		
		std::future<void> fu = std::async(std::launch::async, *task, item, result_nef, &pp.first, &pp.second);
		threadpool_.emplace_back(std::move(fu));
	}
	
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


#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>

void radius_execution_context::set_threads(size_t n) {
	if (!threads_) {
		threads_ = n;
		// threadpool_.resize(n);
		/*
		padding_volumes_.resize(n);
		padding_vol_ptr_.resize(n);
		for (size_t i = 0; i < n; ++i) {
			padding_volumes_[i] = std::make_pair(construct_padding_volume_(), construct_padding_volume_());
			padding_vol_ptr_[i] = &padding_volumes_[i];
		}
		*/		
	}	
}

// Completes the boolean union, extracts exterior and erodes padding radius
void radius_execution_context::finalize() {

	CGAL::Nef_nary_union_3< CGAL::Nef_polyhedron_3<Kernel_> > union_collector;

	for (auto& fu : threadpool_) {
		if (fu.valid()) {
			try {
				fu.get();
			}
			catch (std::exception& e) {
				std::cerr << e.what() << std::endl;
			}
			catch (...) {
				std::cerr << "unkown error" << std::endl;
			}
		}
	}

	auto T = timer::measure("nef_boolean_union");

	for (auto& p : product_geometries) {
		// @todo This part can still be multithreaded
		if (p.second.size() > 1) {
			CGAL::Nef_nary_union_3<CGAL::Nef_polyhedron_3<Kernel_>> per_product_collector;
			for (auto& r : p.second) {
				per_product_collector.add_polyhedron(r);
			}
			p.second = { per_product_collector.get_union() };
			// @todo is this necessary when using the n-ary op?
			// p.second.front().extract_regularization();
		}		
		union_collector.add_polyhedron(p.second.front());
	}

	for (auto& p : reused_products) {
		auto copy = product_geometries[p.second.target].front();
		copy.transform(p.second.inverse);
		copy.transform(p.second.own);
		union_collector.add_polyhedron(copy);
	}

	// @todo spatial sorting?
	auto boolean_result = union_collector.get_union();
	// boolean_result.extract_regularization();
	T.stop();

	if (no_erosion_) {
		exterior = boolean_result;

		if (exterior.is_simple()) {
			polyhedron_exterior = ifcopenshell::geometry::utils::create_polyhedron(exterior);
		}
		else {
			CGAL::convert_nef_polyhedron_to_polygon_mesh(exterior, polyhedron_exterior);
		}
	} else {

		auto T2 = timer::measure("result_nef_processing");
		polyhedron_exterior = ifcopenshell::geometry::utils::create_polyhedron(boolean_result);

		{
			simple_obj_writer tmp_debug("debug-after-boolean");
			tmp_debug(nullptr, polyhedron_exterior.facets_begin(), polyhedron_exterior.facets_end());
		}

		// @todo Wasteful: remove interior on Nef?
		extract_in_place(polyhedron_exterior, LARGEST_AREA);

		if (settings_.get(radius_settings::SPHERE)) {
			typedef CGAL::Simple_cartesian<double> TriangleKernel;
			CGAL::Polyhedron_3<TriangleKernel> poly_simple;

			util::copy::polyhedron(poly_simple, polyhedron_exterior);

			CGAL::Polygon_mesh_processing::triangulate_faces(poly_simple);

			poly_simple.normalize_border();
			if (!poly_simple.is_valid(false, 1)) {
				std::cerr << "invalid before clustering" << std::endl;
				return;
			}

			if (!cluster_vertices(poly_simple, ico_edge_length / 2.)) {
				return;
			}

			poly_simple.normalize_border();
			if (!poly_simple.is_valid(false, 1)) {
				std::cerr << "invalid after clustering" << std::endl;
				return;
			}

			{
				simple_obj_writer tmp_debug("debug-after-custering");
				tmp_debug(nullptr, poly_simple.facets_begin(), poly_simple.facets_end());
			}

			util::copy::polyhedron(polyhedron_exterior, poly_simple);

			polyhedron_exterior.normalize_border();
			if (!polyhedron_exterior.is_valid(false, 1)) {
				std::cerr << "invalid after conversion" << std::endl;
				return;
			}
			
			/*
			// SMS (again) does not really work, geometrical constraints not sattisfied
			Stats stats;
			My_visitor<TriangleKernel> vis(&stats);

			// simplify as small detail create large normal artefacts
			CGAL::Surface_mesh_simplification::Edge_length_stop_predicate<double> stop(length_threshold * length_threshold);
			int r = SMS::edge_collapse(poly_simple, stop,
				CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index, poly_simple))
				.halfedge_index_map(get(CGAL::halfedge_external_index, poly_simple))
				.get_cost(SMS::Edge_length_cost<TriangleKernel>())
				.get_placement(SMS::Midpoint_placement<TriangleKernel>())
				.visitor(vis));

			std::cout << "Removed " << r << " edges" << std::endl;

			
			*/

			// calculate normals
			std::map<cgal_vertex_descriptor_t, Kernel_::Vector_3> vertex_normals;
			boost::associative_property_map<std::map<cgal_vertex_descriptor_t, Kernel_::Vector_3>> vertex_normals_map(vertex_normals);
			CGAL::Polygon_mesh_processing::compute_vertex_normals(polyhedron_exterior, vertex_normals_map);
			
			std::map<cgal_face_descriptor_t, Kernel_::Vector_3> face_normals;
			boost::associative_property_map<std::map<cgal_face_descriptor_t, Kernel_::Vector_3>> face_normals_map(face_normals);
			CGAL::Polygon_mesh_processing::compute_face_normals(polyhedron_exterior, face_normals_map);

			for (auto& v : vertices(polyhedron_exterior)) {

				typedef CGAL::Face_around_target_circulator<cgal_shape_t> face_around_target_circulator;
				face_around_target_circulator circ(v->halfedge(), polyhedron_exterior);
				auto done = circ;
				CGAL::Vector_3<Kernel_> accum, vnorm;
				size_t n = 0;
				do {
					
					cgal_face_descriptor_t f = *circ;
					std::list<cgal_face_descriptor_t> fs = { f };
					auto A = CGAL::Polygon_mesh_processing::area(fs, polyhedron_exterior);
					if (A > 1.e-5) {
						accum += face_normals[f];
						++n;
					}
					++circ;
				} while (circ != done);

				// We don't just use vertex normals, but use connected facet normals if their size is above a certain threshold.
				
				// Let's try how it goes with the vertex clustering, we set n to zero.
				n = 0;

				if (n) {
					vnorm = accum / n;
				}
				else {
					vnorm = vertex_normals[v];
				}

				/*
				std::cout << "N " << CGAL::to_double(vnorm.cartesian(0))
					<< " " << CGAL::to_double(vnorm.cartesian(1))
					<< " " << CGAL::to_double(vnorm.cartesian(2)) << std::endl;

				std::cout << "p0 " << CGAL::to_double(v->point().cartesian(0))
					<< " " << CGAL::to_double(v->point().cartesian(1))
					<< " " << CGAL::to_double(v->point().cartesian(2)) << std::endl;
				*/
				
				v->point() -= vnorm * radius;

				/*
				std::cout << "p1 " << CGAL::to_double(v->point().cartesian(0))
					<< " " << CGAL::to_double(v->point().cartesian(1))
					<< " " << CGAL::to_double(v->point().cartesian(2)) << std::endl;
				*/
			}

			{
				simple_obj_writer tmp_debug("debug-after-erosion");
				tmp_debug(nullptr, polyhedron_exterior.facets_begin(), polyhedron_exterior.facets_end());
			}

			if (!cluster_vertices(polyhedron_exterior, ico_edge_length  / 2.)) {
				return;
			}

			{
				simple_obj_writer tmp_debug("debug-after-custering-again");
				tmp_debug(nullptr, polyhedron_exterior.facets_begin(), polyhedron_exterior.facets_end());
			}

			// util::copy::polyhedron(polyhedron_exterior, poly_simple);
		}
		else {
			throw std::runtime_error("todo");
			
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
			complement_padded = complement; // CGAL::minkowski_sum_3(complement, padding_volume_2);
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
			}
			else {
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
			auto padding_volume = construct_padding_volume_();
			minkowski_sum_triangles_single_threaded<CGAL::Polyhedron_3<CGAL::Epick>>(
				poly_triangulated.facets_begin(),
				poly_triangulated.facets_end(),
				padding_volume, exterior
				);
			if (exterior.is_simple()) {
				polyhedron_exterior = ifcopenshell::geometry::utils::create_polyhedron(exterior);
			}
			else {
				CGAL::convert_nef_polyhedron_to_polygon_mesh(exterior, polyhedron_exterior);
			}

			extract_in_place(polyhedron_exterior, SECOND_LARGEST_AREA);
#endif
		}

		T2.stop();
	}

	std::cout << "exterior poly num facets: " << polyhedron_exterior.size_of_facets() << std::endl;
}

