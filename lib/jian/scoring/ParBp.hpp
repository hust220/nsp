#pragma once

#include <memory>
#include "../matrix.hpp"
#include "../pdb.hpp"
#include "../cg.hpp"

namespace jian {
	class ParBp {
	public:
		double d, z, alpha, theta;

		std::shared_ptr<CG> m_cg;

		ParBp();

		ParBp(const Residue &res1, const Residue &res2);

		void set_origin(Vec &origin, const Residue &r);

		void vec_cross(Vec &v, const Vec &a, const Vec &b);

		void set_norm(Vec &norm, const Vec &origin, const Residue &r);

		void set_axis_x(Vec &x, const Vec &origin, const Residue &r);

		bool is_paired() const;

		bool is_stacked() const;

		bool is_wc() const;

		bool is_nwc() const;

		ParBp &anal(const Residue &res1, const Residue &res2);
	};


}