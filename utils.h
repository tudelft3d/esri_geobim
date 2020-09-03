#ifndef UTILS_H
#define UTILS_H

#include <ifcgeom/kernels/cgal/CgalKernel.h>

#include <CGAL/boost/graph/copy_face_graph.h>

namespace util {

	namespace copy {

		template <class Poly_B, class Poly_A>
		typename std::enable_if<std::is_same<Poly_A, Poly_B>::value>::type polyhedron(Poly_B& poly_b, const Poly_A& poly_a) {
			poly_b = poly_a;
		}

		template <class Poly_B, class Poly_A>
		typename std::enable_if<!std::is_same<Poly_A, Poly_B>::value>::type polyhedron(Poly_B& poly_b, const Poly_A& poly_a) {
			poly_b.clear();
			CGAL::copy_face_graph(poly_a, poly_b);
		}

		template <class Kb, class Ka>
		typename std::enable_if<std::is_same<Ka, Kb>::value>::type transformation(CGAL::Aff_transformation_3<Kb>& trsf_b, const CGAL::Aff_transformation_3<Ka>& trsf_a) {
			trsf_b = trsf_a;
		}

		template <class Kb, class Ka>
		typename std::enable_if<!std::is_same<Ka, Kb>::value>::type transformation(CGAL::Aff_transformation_3<Kb>& trsf_b, const CGAL::Aff_transformation_3<Ka>& trsf_a) {
			CGAL::NT_converter<
				typename Ka::RT,
				typename Kb::RT> converter;
			trsf_b = CGAL::Aff_transformation_3<Kb>(
				converter(trsf_a.hm(0, 0)),
				converter(trsf_a.hm(0, 1)),
				converter(trsf_a.hm(0, 2)),
				converter(trsf_a.hm(0, 3)),
				converter(trsf_a.hm(1, 0)),
				converter(trsf_a.hm(1, 1)),
				converter(trsf_a.hm(1, 2)),
				converter(trsf_a.hm(1, 3)),
				converter(trsf_a.hm(2, 0)),
				converter(trsf_a.hm(2, 1)),
				converter(trsf_a.hm(2, 2)),
				converter(trsf_a.hm(2, 3)),
				converter(trsf_a.hm(3, 3))
				);
		}

	}
	
	template <class HDS>
	class PolyFromMesh : public CGAL::Modifier_base<HDS> {
	private:

		std::list<cgal_point_t> points_;
		std::vector<std::vector<int>> indices_;

	public:
		PolyFromMesh(const std::list<cgal_point_t>& points, const std::vector<std::vector<int>>& indices)
			: points_(points)
			, indices_(indices) {}

		void operator()(HDS& hds) {
			CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, false);

			B.begin_surface(points_.size(), indices_.size());

			for (auto& p : points_) {
				B.add_vertex(p);
			}

			for (auto& fs : indices_) {
				B.begin_facet();
				for (auto& i : fs) {
					B.add_vertex_to_facet(i);
				}
				B.end_facet();
			}

			B.end_surface();
		}
	};

}

#endif
