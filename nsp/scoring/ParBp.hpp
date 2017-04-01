#pragma once

#include <memory>
#include "jian/matrix.hpp"
#include "../pdb.hpp"
#include "../cg.hpp"

BEGIN_JN

Bool is_bp(const Residue &r1, const Residue &r2);

class ParBp {
public:
	using vec_t = Eigen::Vector3d;

	double d, alpha, theta;
	vec_t o1, o2, o12, o21, o12_, o21_;
	vec_t bb1, bb2, bb12, bb21, bb12_, bb21_;
	std::array<vec_t, 3> ax1, ax2;

	std::shared_ptr<CG> m_cg;

	ParBp();

	ParBp(const Residue &res1, const Residue &res2);

	void set_origin(vec_t &origin, const Residue &r);

	void set_bb(vec_t &bb, const Residue &r);

	void vec_cross(vec_t &v, const vec_t &a, const vec_t &b);

	void set_axis_z(vec_t &norm, const vec_t &origin, const Residue &r);

	void set_axis_x(vec_t &x, const vec_t &origin, const Residue &r);

	void set_axis_y(vec_t &y, const vec_t &z, const vec_t &x);

	bool is_paired() const;

	bool is_stacked() const;

	bool is_wc() const;

	bool is_nwc() const;

	ParBp &anal(const Residue &res1, const Residue &res2);
};

END_JN
