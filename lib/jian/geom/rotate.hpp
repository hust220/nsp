#pragma once

#include <Eigen/Core>
#include <cmath>
#include "../matrix.hpp"

namespace jian {
	namespace geom {

		using val_t = double;

		template<typename T, typename U>
		void translate(T &&t, const U &u) {
			for (int i = 0; i < 3; i++) t[i] += u[i];
		}

		// Rotate fixed in the origin
		template<typename T, typename Mat>
		void rotate(T &&t, const Mat &mat) {
			val_t x = t[0] * ref(mat, 0, 0) + t[1] * ref(mat, 1, 0) + t[2] * ref(mat, 2, 0);
			val_t y = t[0] * ref(mat, 0, 1) + t[1] * ref(mat, 1, 1) + t[2] * ref(mat, 2, 1);
			val_t z = t[0] * ref(mat, 0, 2) + t[1] * ref(mat, 1, 2) + t[2] * ref(mat, 2, 2);
			t[0] = x; t[1] = y; t[2] = z;
		}

		// Rotate fixed in an specified point
		template<typename T, typename U, typename Mat>
		void rotate(T &&t, const U &origin, const Mat &mat) {
			for (int i = 0; i < 3; i++) t[i] -= origin[i];
			val_t x = t[0] * ref(mat, 0, 0) + t[1] * ref(mat, 1, 0) + t[2] * ref(mat, 2, 0);
			val_t y = t[0] * ref(mat, 0, 1) + t[1] * ref(mat, 1, 1) + t[2] * ref(mat, 2, 1);
			val_t z = t[0] * ref(mat, 0, 2) + t[1] * ref(mat, 1, 2) + t[2] * ref(mat, 2, 2);
			t[0] = x; t[1] = y; t[2] = z;
			for (int i = 0; i < 3; i++) t[i] += origin[i];
		}

		template<typename Mat = Eigen::MatrixXd, typename T>
		Mat rot_mat(int i, T &&v) {
			double c = std::cos(v);
			double s = std::sin(v);
			assert(i >= 0 && i < 3);
			return (i == 0 ? x_rot_mat(c, s) : (i == 1 ? y_rot_mat(c, s) : z_rot_mat(c, s)));
		}

		template<typename Mat = Eigen::MatrixXd, class C1 = val_t, class C2 = val_t>
		Mat x_rot_mat(C1 c, C2 s) {
			Mat rot_mat(3, 3);
			rot_mat <<
				1, 0, 0,
				0, c, s,
				0, -s, c;
			return rot_mat;
		}

		template<typename Mat = Eigen::MatrixXd, class C1 = val_t, class C2 = val_t>
		Mat y_rot_mat(C1 c, C2 s) {
			Mat rot_mat(3, 3);
			rot_mat <<
				c, 0, -s,
				0, 1, 0,
				s, 0, c;
			return rot_mat;
		}

		template<typename Mat = Eigen::MatrixXd, class C1 = val_t, class C2 = val_t>
		Mat z_rot_mat(C1 c, C2 s) {
			Mat rot_mat(3, 3);
			rot_mat <<
				c, s, 0,
				-s, c, 0,
				0, 0, 1;
			return rot_mat;
		}

		// Rotate along with an axis
		class RotateAlong {
		public:
			Mat m_rm;
			Vec m_beg, m_end;

			RotateAlong() = default;

			template<typename L>
			RotateAlong(const L &begin, const L &end, double angle) {
				init(begin, end, angle);
			}

			template<typename L>
			void init(const L &begin, const L &end, double t) {
				val_t r1, r2, c1, c2, s1, s2;
				L l = end - begin;

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

				m_rm = Mat::Identity(3, 3);
				if (r1 != 0) m_rm = m_rm * z_rot_mat(c1, s1);
				if (r2 != 0) m_rm = m_rm * x_rot_mat(c2, s2);
				m_rm = m_rm * z_rot_mat(std::cos(t), std::sin(t));
				if (r2 != 0) m_rm = m_rm * x_rot_mat(c2, -s2);
				if (r1 != 0) m_rm = m_rm * z_rot_mat(c1, -s1);

			}

			template<typename T>
			RotateAlong &operator ()(T &&t) {
				for (int i = 0; i < 3; i++) {
					t[i] -= m_beg[i];
				}
				rotate(t, m_rm);
				for (int i = 0; i < 3; i++) {
					t[i] += m_beg[i];
				}
				return *this;
			}

		};

	} // namespace geom
} // namespace jian

