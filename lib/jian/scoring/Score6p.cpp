#include "Score6p.hpp"

namespace jian {
	REG_SCORER("6p", Score6p);

	Score6p::Score6p() {
		m_cg = CG::fac_t::create("6p");
	}

	inline double Score6p::en_len(const Residue &r1, const Residue &r2) {
		double d;

		d = geom::distance(r1[2], r2[0]);
		return square(d - 3.1);
	}

	inline double Score6p::en_ang(const Residue &r1, const Residue &r2, const Residue &r3) {
		double d;

		d = geom::angle(r1[1], r2[1], r3[1]);
		return square(d - m_bond_angle_std);
	}

	inline double Score6p::en_dih(const Residue &r1, const Residue &r2, const Residue &r3, const Residue &r4) {
		return 0;
	}

	inline double Score6p::en_crash(const Residue &r1, const Residue &r2) {
		int i, j;
		double d, e;

		e = 0;
		for (i = 0; i < m_res_size; i++) {
			for (j = 0; j < m_res_size; j++) {
				d = geom::distance(r1[i], r2[j]);
				if (d < 3.5) {
					e += square(d - 3.5);
				}
			}
		}
		return e;
	}

} // namespace jian