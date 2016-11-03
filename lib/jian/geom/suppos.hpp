#pragma once

#include "rotate.hpp"
#include "../utils/math.hpp"
#include "../utils/string.hpp"

namespace jian {
	namespace geom {

#define INIT_SUPPOS(sp) \
    auto &&sp##C1 = -(sp.c1); auto &&sp##ROT = sp.rot; auto &&sp##C2 = sp.c2;

#define APPLY_SUPPOS(m, sp) \
    do{geom::translate(m, sp##C1); geom::rotate(m, sp##ROT); geom::translate(m, sp##C2);}while(0)

#define APPLY_SUPPOS_m(m, sp) do{\
    for (int i = 0; i < m.rows(); i++) {\
        for (int j = 0; j < 3; j++) {\
            m(i, j) += sp##C1[j];\
        }\
    }\
    m = m * sp##ROT;\
    for (int i = 0; i < m.rows(); i++) {\
        for (int j = 0; j < 3; j++) {\
            m(i, j) += sp##C2[j];\
        }\
    }\
}while(0)

		class Superposition {
		public:
			Mat rot;
			Vec c1;
			Vec c2;
			val_t rmsd;

			Superposition() = default;

			Superposition(const Mat &m, const Mat &n);

			void init(const Mat &m, const Mat &n);

			template<typename T>
			void apply(T &&t) {
				translate(t, c1);
				rotate(t, rot);
				translate(t, c2);
			}

			template<typename Mat>
			void apply_m(Mat &&m) {
				int i, j;
				for (i = 0; i < m.rows(); i++) {
					for (j = 0; j < 3; j++) {
						m(i, j) += c1[j];
					}
				}
				m = m * rot;
				for (i = 0; i < m.rows(); i++) {
					for (j = 0; j < 3; j++) {
						m(i, j) += c2[j];
					}
				}
			}
		};

		Superposition suppos(const Mat &m, const Mat &n);

		template<typename T, typename U, typename F>
		Superposition suppos(T &t, const U &m, const F &n) {
			auto sp = suppos(m, n);
			for (int i = 0; i < t.rows(); i++) for (int j = 0; j < 3; j++) t(i, j) -= sp.c1[j];
			t = t * sp.rot;
			for (int i = 0; i < t.rows(); i++) for (int j = 0; j < 3; j++) t(i, j) += sp.c2[j];
			return sp;
		}

		template<typename T, typename U, typename F, typename V>
		Superposition suppos(T &src, const U &src_indices, const F &tgt, const V &tgt_indices) {
			Eigen::MatrixXd m(src_indices.size(), 3), n(tgt_indices.size(), 3);
			for (int i = 0; i < src_indices.size(); i++) for (int j = 0; j < 3; j++) {
				m(i, j) = src(src_indices[i], j);
				n(i, j) = tgt(tgt_indices[i], j);
			}

			auto sp = suppos(m, n);

			for (int i = 0; i < src.rows(); i++) for (int j = 0; j < 3; j++) src(i, j) -= sp.c1[j];
			src = src * sp.rot;
			for (int i = 0; i < src.rows(); i++) for (int j = 0; j < 3; j++) src(i, j) += sp.c2[j];
			return sp;
		}

		val_t rmsd(const Mat &x, const Mat &y);

		Mat suppos_axis_polar(double theta_o, double phi_o, double theta_n, double phi_n);

		template<typename Mat, typename O, typename N>
		Mat suppos_axis_xyz(const O &o, const N &n) {
			Mat m = Mat::Identity(3, 3);
			val_t r, x, y, r1, x1, y1, r2, x2, y2, c, s, c1, c2, s1, s2;
			r = std::sqrt(o[0] * o[0] + o[1] * o[1]); x = o[0]; y = o[1];
			if (r != 0) { c = y / r; s = x / r; if (s != 0) m *= z_rot_mat(c, s); }
			r1 = std::sqrt(o[0] * o[0] + o[1] * o[1] + o[2] * o[2]); x1 = std::sqrt(o[0] * o[0] + o[1] * o[1]); y1 = o[2];
			r2 = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]); x2 = std::sqrt(n[0] * n[0] + n[1] * n[1]); y2 = n[2];
			if (r1 != 0 && r2 != 0) {
				c1 = y1 / r1; s1 = x1 / r1; c2 = y2 / r2; s2 = x2 / r2; c = c1*c2 + s1*s2; s = s1*c2 - s2*c1;
				if (s != 0) m *= x_rot_mat(c, s);
			}
			r = std::sqrt(n[0] * n[0] + n[1] * n[1]); x = n[0]; y = n[1];
			if (r != 0) { c = y / r; s = -x / r; if (s != 0) m *= z_rot_mat(c, s); }
			return m;
		}

	} // namespace geom
} // namespace jian

