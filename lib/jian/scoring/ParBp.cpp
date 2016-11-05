#include "ParBp.hpp"

namespace jian {
	ParBp::ParBp() {
		m_cg.reset(CG::fac_t::create("6p"));
	}

	ParBp::ParBp(const Residue &res1, const Residue &res2) : ParBp() {
		anal(res1, res2);
	}

	void ParBp::set_origin(Vec &origin, const Residue &r) {
		origin = Vec::Zero(3);
		for (int i = 3; i < 6; i++) {
			for (int j = 0; j < 3; j++) {
				origin[j] += r[i][j];
			}
		}
		origin /= 3.0;
		//std::cout << "origin: " << origin << std::endl;
	}

	void ParBp::vec_cross(Vec &v, const Vec &a, const Vec &b) {
		v[0] = a[1] * b[2] - a[2] * b[1];
		v[1] = a[2] * b[0] - a[0] * b[2];
		v[2] = a[0] * b[1] - a[1] * b[0];
	}

	void ParBp::set_norm(Vec &norm, const Vec &origin, const Residue &r) {
		Vec a(3), b(3);
		for (int i = 0; i < 3; i++) {
			a[i] = r[3][i] - origin[i];
			b[i] = r[4][i] - origin[i];
		}
		vec_cross(norm, a, b);
		norm.normalize();
	}

	void ParBp::set_axis_x(Vec &x, const Vec &origin, const Residue &r) {
		int i, j, t;

		t = pdb::res_type(r.name);
		j = ((t == 0 || t == 2) ? 4 : 5);
		for (i = 0; i < 3; i++) {
			x[i] = r[j][i] - origin[i];
		}
	}

	bool ParBp::is_paired() const {
		return d < 8 && std::fabs(z) < 2;
	}

	bool ParBp::is_stacked() const {
		return d < 3 && std::fabs(z) < 4 && std::fabs(theta) > 0.90;
	}

	bool ParBp::is_wc() const {
		return is_paired() && d < 6.1 && std::fabs(alpha) < 15;
	}

	bool ParBp::is_nwc() const {
		return is_paired() && !is_wc();
	}

	ParBp &ParBp::anal(const Residue &res1, const Residue &res2) {
		Vec origin1, origin2, norm1(3), norm2(3), x1(3), y1(3);
		Residue r1 = m_cg->to_cg(res1);
		Residue r2 = m_cg->to_cg(res2);
		double origin2_x, origin2_y;

		set_origin(origin1, r1);
		set_origin(origin2, r2);
		set_norm(norm1, origin1, r1);
		set_norm(norm2, origin2, r2);

		set_axis_x(x1, origin1, r1);
		vec_cross(y1, norm1, x1);
		origin2 -= origin1;
		//d = origin2.norm();
		origin2_x = origin2.dot(x1) / x1.norm();
		origin2_y = origin2.dot(y1) / y1.norm();
		d = std::sqrt(origin2_x*origin2_x+ origin2_y*origin2_y);
		z = origin2.dot(norm1);
		alpha = atan(origin2_y / origin2_x) / 3.1415927 * 180;
		//theta = acos(norm1.dot(norm2)) / 3.1415927 * 180;
		theta = norm1.dot(norm2);
		return *this;
	}

}