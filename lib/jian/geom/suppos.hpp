#pragma once

#include "rotate.hpp"
#include "../utils/math.hpp"
#include "../utils/string.hpp"

BEGIN_JN
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

		template<typename NumType>
		class Superposition {
		public:
			using mat_t = MatX<NumType>;
			using vec_t = VecX<NumType>;

			mat_t rot;
			vec_t c1;
			vec_t c2;
			NumType rmsd;

			Superposition() = default;

			Superposition(const mat_t &m, const mat_t &n) {
				init(m, n);
			}

			void init(const mat_t &m, const mat_t &n) {
				int i, j, len;
				mat_t x, y, g, u, v, I, d;
				std::ostringstream stream;

				if (m.rows() != n.rows() || m.cols() != 3 || n.cols() != 3) {
					stream << "jian::geom::suppos error! (" << m.rows() << ' ' << m.cols() << ") (" << n.rows() << ' ' << n.cols() << ")\n";
					throw stream.str();
				}

				len = m.rows();
				x = m;
				y = n;
				c1 = vec_t::Zero(3);
				c2 = vec_t::Zero(3);
				for (i = 0; i < len; i++) for (j = 0; j < 3; j++) { c1[j] += x(i, j); c2[j] += y(i, j); }
				for (i = 0; i < 3; i++) { c1[i] = c1[i] / len; c2[i] = c2[i] / len; }
				for (i = 0; i < len; i++) for (j = 0; j < 3; j++) { x(i, j) -= c1[j]; y(i, j) -= c2[j]; }

				g = x.transpose() * y;
				Eigen::JacobiSVD<mat_t> svd(g, Eigen::ComputeFullU | Eigen::ComputeFullV);
				u = svd.matrixU();
				v = svd.matrixV();

				I = mat_t::Identity(3, 3);
				if (g.determinant() < 0) I(2, 2) = -1;
				rot = u * I * v.transpose();

				d = x * rot - y;
				rmsd = 0;
				for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) rmsd += d(i, j) * d(i, j);
				rmsd = std::sqrt(rmsd / len);

				c1 = -c1;
			}

			template<typename T>
			void apply(T &&t) {
				translate(t, c1);
				rotate(t, rot);
				translate(t, c2);
			}

			template<typename M>
			void apply_m(M &m) {
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

		template<typename NumType>
		Superposition<NumType> suppos(const MatX<NumType> &m, const MatX<NumType> &n) {
			return Superposition<NumType>(m, n);
			//sp.c1 = -(sp.c1);
		}

		template<typename T, typename NumType>
		Superposition<NumType> suppos(T &t, const MatX<NumType> &m, const MatX<NumType> &n) {
			Superposition<NumType> sp(m, n);
            sp.apply_m(t);
//			for (int i = 0; i < t.rows(); i++) for (int j = 0; j < 3; j++) t(i, j) -= sp.c1[j];
//			t = t * sp.rot;
//			for (int i = 0; i < t.rows(); i++) for (int j = 0; j < 3; j++) t(i, j) += sp.c2[j];
			return sp;
		}

		template<typename NumType, typename U, typename F, typename V>
		Superposition<NumType> suppos(MatX<NumType> &src, const U &src_indices, const F &tgt, const V &tgt_indices) {
			MatX<NumType> m(src_indices.size(), 3), n(tgt_indices.size(), 3);
			for (int i = 0; i < src_indices.size(); i++) for (int j = 0; j < 3; j++) {
				m(i, j) = src(src_indices[i], j);
				n(i, j) = tgt(tgt_indices[i], j);
			}

			auto sp = suppos(m, n);

            sp.apply_m(src);

//			for (int i = 0; i < src.rows(); i++) for (int j = 0; j < 3; j++) src(i, j) -= sp.c1[j];
//			src = src * sp.rot;
//			for (int i = 0; i < src.rows(); i++) for (int j = 0; j < 3; j++) src(i, j) += sp.c2[j];
			return sp;
		}

		template<typename NumType1, typename NumType2>
		NumType1 rmsd(const MatX<NumType1> &x, const MatX<NumType2> &y) {
			return Superposition<NumType1>(x, y).rmsd;
		}

		template<typename NumType, typename A, typename B, typename C, typename D>
		MatX<NumType> suppos_axis_polar(A theta_o, B phi_o, C theta_n, D phi_n) {
			MatX<NumType> m = MatX<NumType>::Identity(3, 3);
			double ang = - phi_o;
			if (ang != 0) m *= z_rot_mat(std::cos(ang), std::sin(ang));
			ang = -theta_o + theta_n;
			if (ang != 0) m *= y_rot_mat(std::cos(ang), std::sin(ang));
			ang = phi_n;
			if (ang != 0) m *= z_rot_mat(std::cos(ang), std::sin(ang));
			return m;
		}


		template<typename NumType, typename O, typename N>
		MatX<NumType> suppos_axis_xyz(const O &o, const N &n) {
			MatX<NumType> m = MatX<NumType>::Identity(3, 3);
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
END_JN

