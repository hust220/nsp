#include "ParBp.hpp"

BEGIN_JN

ParBp::ParBp() {
	m_cg.reset(CG::fac_t::create("6p"));
}

ParBp::ParBp(const Residue &res1, const Residue &res2) {
	m_cg.reset(CG::fac_t::create("6p"));
	anal(res1, res2);
}

inline void ParBp::set_origin(vec_t &origin, const Residue &r) {
	origin = vec_t::Zero();
	for (int i = 3; i < 6; i++) {
		for (int j = 0; j < 3; j++) {
			origin[j] += r[i][j];
		}
	}
	origin /= 3.0;
}

inline void ParBp::set_bb(vec_t &bb, const Residue &r) {
	int i;

	for (i = 0; i < 3; i++) {
		bb[i] = r[0][i];
	}
}

inline void ParBp::vec_cross(vec_t &v, const vec_t &a, const vec_t &b) {
	v[0] = a[1] * b[2] - a[2] * b[1];
	v[1] = a[2] * b[0] - a[0] * b[2];
	v[2] = a[0] * b[1] - a[1] * b[0];
}

inline void ParBp::set_axis_z(vec_t &z, const vec_t &o, const Residue &r) {
	vec_t a(3), b(3);
	int i, j, t;

	t = pdb::res_type(r.name);
	j = ((t == 0 || t == 2) ? 4 : 5);
	for (i = 0; i < 3; i++) {
		a[i] = r[3][i] - o[i];
		b[i] = r[j][i] - o[i];
	}
	vec_cross(z, a, b);
	z.normalize();
}

inline void ParBp::set_axis_x(vec_t &x, const vec_t &o, const Residue &r) {
	int i, j, t;

	t = pdb::res_type(r.name);
	j = ((t == 0 || t == 2) ? 4 : 5);
	for (i = 0; i < 3; i++) {
		x[i] = r[j][i] - o[i];
	}
	x.normalize();
}

inline void ParBp::set_axis_y(vec_t &y, const vec_t &z, const vec_t &x) {
	vec_cross(y, z, x);
	y.normalize();
}

bool ParBp::is_paired() const {
	//std::cout << d << ' ' << std::fabs(o21_[2]) << ' ' << std::fabs(o12_[2]) << std::endl;
	return d < 8 && std::fabs(o21_[2]) < 2 && std::fabs(o12_[2]) < 2;
	//return d < 8 && std::fabs(o21_[2]) < 2 && std::fabs(o12_[2]) < 2;
}

bool ParBp::is_stacked() const {
	return d < 4 && std::fabs(o21_[2]) < 4 && std::fabs(o21_[2]) < 4 && std::fabs(theta) > 0.90;
}

bool ParBp::is_wc() const {
	return is_paired() && d < 6.1 && std::fabs(alpha) < 15;
}

bool ParBp::is_nwc() const {
	return is_paired() && !is_wc();
}

ParBp &ParBp::anal(const Residue &res1, const Residue &res2) {
	int i;
	Residue r1 = m_cg->to_cg(res1);
	Residue r2 = m_cg->to_cg(res2);

	set_origin(o1, r1);
	set_origin(o2, r2);
	set_bb(bb1, r1);
	set_bb(bb2, r2);

	set_axis_z(ax1[2], o1, r1);
	set_axis_x(ax1[0], o1, r1);
	set_axis_y(ax1[1], ax1[2], ax1[0]);

	set_axis_z(ax2[2], o2, r2);
	set_axis_x(ax2[0], o2, r2);
	set_axis_y(ax2[1], ax2[2], ax2[0]);

	o12 = o1 - o2;
	o21 = o2 - o1;
	bb12 = bb1 - o2;
	bb21 = bb2 - o1;

	d = o21.norm();

	for (i = 0; i < 3; i++) {
		o12_[i] = o12.dot(ax2[i]);
		o21_[i] = o21.dot(ax1[i]);
		bb12_[i] = bb12.dot(ax2[i]);
		bb21_[i] = bb21.dot(ax1[i]);
	}

	alpha = atan(o21_[1] / o21_[0]) / 3.1415927 * 180;
	theta = ax1[2].dot(ax2[2]);
	return *this;
}

END_JN