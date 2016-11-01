#include "suppos.hpp"

namespace jian {
	namespace geom {

		Superposition::Superposition(const Mat &m, const Mat &n) {
			init(m, n);
		}

		Superposition &Superposition::init(const Mat &m, const Mat &n) {
			int i, j, len;
			Mat x, y, g, u, v, I, d;
			std::ostringstream stream;

			if (m.rows() != n.rows() || m.cols() != 3 || n.cols() != 3) {
				stream << "jian::geom::suppos error! (" << m.rows() << ' ' << m.cols() << ") (" << n.rows() << ' ' << n.cols() << ")\n";
				throw stream.str();
			}

			len = m.rows();
			x = m;
			y = n;
			c1 = Vec::Zero(3);
			c2 = Vec::Zero(3);
			for (i = 0; i < len; i++) for (j = 0; j < 3; j++) { c1[j] += x(i, j); c2[j] += y(i, j); }
			for (i = 0; i < 3; i++) { c1[i] = c1[i] / len; c2[i] = c2[i] / len; }
			for (i = 0; i < len; i++) for (j = 0; j < 3; j++) { x(i, j) -= c1[j]; y(i, j) -= c2[j]; }

			g = x.transpose() * y;
			Eigen::JacobiSVD<Mat> svd(g, Eigen::ComputeFullU | Eigen::ComputeFullV);
			u = svd.matrixU();
			v = svd.matrixV();

			I = Eigen::MatrixXd::Identity(3, 3);
			if (g.determinant() < 0) I(2, 2) = -1;
			rot = u * I * v.transpose();

			d = x * rot - y;
			rmsd = 0; for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) rmsd += d(i, j) * d(i, j);
			rmsd = std::sqrt(rmsd / len);

			c1 = -c1;
			return *this;
		}

		Superposition suppos(const Mat &m, const Mat &n) {
			Superposition sp(m, n);
			sp.c1 = -(sp.c1);
			return sp;
		}

		val_t rmsd(const Mat &m, const Mat &n) {
			return Superposition(m, n).rmsd;
		}

		Mat suppos_axis_polar(double theta_o, double phi_o, double theta_n, double phi_n) {
			Mat m = Mat::Identity(3, 3);
			double ang = PI / 2 - phi_o;
			if (ang != 0) m *= z_rot_mat(std::cos(ang), std::sin(ang));
			ang = theta_o - theta_n;
			if (ang != 0) m *= x_rot_mat(std::cos(ang), std::sin(ang));
			ang = phi_n - PI / 2;
			if (ang != 0) m *= z_rot_mat(std::cos(ang), std::sin(ang));
			return m;
		}

	} // namespace geom
} // namespace jian

