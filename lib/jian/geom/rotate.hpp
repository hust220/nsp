#pragma once

#include <Eigen/Core>
#include <cmath>
#include "../matrix.hpp"

BEGIN_JN
	namespace geom {

		template<typename T, typename U>
		void translate(T &&t, const U &u) {
			for (int i = 0; i < 3; i++) t[i] += u[i];
		}

		// Rotate fixed in the origin
		template<typename T, typename NumType>
		void rotate(T &&t, const MatX<NumType> &mat) {
			NumType x = t[0] * ref(mat, 0, 0) + t[1] * ref(mat, 1, 0) + t[2] * ref(mat, 2, 0);
			NumType y = t[0] * ref(mat, 0, 1) + t[1] * ref(mat, 1, 1) + t[2] * ref(mat, 2, 1);
			NumType z = t[0] * ref(mat, 0, 2) + t[1] * ref(mat, 1, 2) + t[2] * ref(mat, 2, 2);
			t[0] = x; t[1] = y; t[2] = z;
		}

		// Rotate fixed in an specified point
		template<typename T, typename U, typename NumType>
		void rotate(T &&t, const U &origin, const MatX<NumType> &mat) {
			for (int i = 0; i < 3; i++) t[i] -= origin[i];
			rotate(t, mat);
			//NumType x = t[0] * ref(mat, 0, 0) + t[1] * ref(mat, 1, 0) + t[2] * ref(mat, 2, 0);
			//NumType y = t[0] * ref(mat, 0, 1) + t[1] * ref(mat, 1, 1) + t[2] * ref(mat, 2, 1);
			//NumType z = t[0] * ref(mat, 0, 2) + t[1] * ref(mat, 1, 2) + t[2] * ref(mat, 2, 2);
			//t[0] = x; t[1] = y; t[2] = z;
			for (int i = 0; i < 3; i++) t[i] += origin[i];
		}

		template<typename NumType = val_t, class C1 = NumType, class C2 = NumType>
		MatX<NumType> x_rot_mat(C1 c, C2 s) {
			MatX<NumType> rot_mat(3, 3);
			rot_mat <<
				1, 0, 0,
				0, c, s,
				0, -s, c;
			return rot_mat;
		}

		template<typename NumType = val_t, class C1 = NumType, class C2 = NumType>
		MatX<NumType> y_rot_mat(C1 c, C2 s) {
			MatX<NumType> rot_mat(3, 3);
			rot_mat <<
				c, 0, -s,
				0, 1, 0,
				s, 0, c;
			return rot_mat;
		}

		template<typename NumType = val_t, class C1 = NumType, class C2 = NumType>
		MatX<NumType> z_rot_mat(C1 c, C2 s) {
			MatX<NumType> rot_mat(3, 3);
			rot_mat <<
				c, s, 0,
				-s, c, 0,
				0, 0, 1;
			return rot_mat;
		}

		template<typename NumType = val_t, typename T>
		MatX<NumType> rot_mat(int i, T &&v) {
			NumType c = std::cos(v);
			NumType s = std::sin(v);
			assert(i >= 0 && i < 3);
			return (i == 0 ? x_rot_mat(c, s) : (i == 1 ? y_rot_mat(c, s) : z_rot_mat(c, s)));
		}

		// Rotate along with an axis
		template<typename NumType>
		class RotateAlong {
		public:
			using mat_t = MatX<NumType>;
			using vec_t = VecX<NumType>;

			mat_t m_rm;
			vec_t m_beg, m_end;

			RotateAlong() = default;

			template<typename Begin, typename End>
			RotateAlong(const Begin &begin, const End &end, double angle) {
				init(begin, end, angle);
			}

			template<typename Begin, typename End>
			void init(const Begin &begin, const End &end, double t) {
				NumType r1, r2, c1, c2, s1, s2;
				//L l = end - begin;
				vec_t l(3);
				for (int i = 0; i < 3; i++) l[i] = end[i] - begin[i];

				m_beg.resize(3);
				m_end.resize(3);
				for (int i = 0; i < 3; i++) {
					m_beg[i] = begin[i];
					m_end[i] = end[i];
				}

				r1 = std::sqrt(l[0] * l[0] + l[1] * l[1]);
				c1 = l[1] / r1;
				s1 = l[0] / r1;

				r2 = std::sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
				c2 = l[2] / r2;
				s2 = std::sqrt(l[0] * l[0] + l[1] * l[1]) / r2;

				m_rm = mat_t::Identity(3, 3);
				if (r1 != 0) m_rm = m_rm * z_rot_mat(c1, s1);
				if (r2 != 0) m_rm = m_rm * x_rot_mat(c2, s2);
				m_rm = m_rm * z_rot_mat(std::cos(t), std::sin(t));
				if (r2 != 0) m_rm = m_rm * x_rot_mat(c2, -s2);
				if (r1 != 0) m_rm = m_rm * z_rot_mat(c1, -s1);

			}

			template<typename T>
			RotateAlong &operator ()(T &&t) {
				rotate(t, m_beg, m_rm);
				//for (int i = 0; i < 3; i++) {
				//	t[i] -= m_beg[i];
				//}
				//rotate(t, m_rm);
				//for (int i = 0; i < 3; i++) {
				//	t[i] += m_beg[i];
				//}
				return *this;
			}

		};

	} // namespace geom
END_JN

