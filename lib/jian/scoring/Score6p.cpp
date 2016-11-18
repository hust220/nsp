#include "Score6p.hpp"

namespace jian {
	REG_SCORER("6p", Score6p);

	Score6p::Score6p() {
		m_cg = CG::fac_t::create("6p");
	}

	inline double Score6p::en_len(const Chain &c, int beg) {
		double d, e;

		e = 0;
		if (beg + 1 < c.size()) {
			d = geom::distance(c[beg][0], c[beg + 1][0]);
			e += square(d - m_bond_len_std[0]);
			d = geom::distance(c[beg][1], c[beg + 1][0]);
			e += square(d - m_bond_len_std[1]);
		}
		//if (beg + 2 < c.size()) {
		//	d = geom::distance(c[beg][0], c[beg + 2][0]);
		//	e += square(d - m_bond_len_std[2]);
		//}
		return e;
	}

	inline double Score6p::en_ang(const Chain &c, int beg) {
		double d, e;

		e = 0;
		if (beg + 2 < c.size()) {
			d = geom::angle(c[beg][0], c[beg+1][0], c[beg + 2][0]);
			e += square(d - m_bond_angle_std[0]);
		}
		return e;
	}

	inline double Score6p::en_dih(const Chain &c, int beg) {
		double d, e;

		e = 0;
		if (beg + 3 < c.size()) {
			d = geom::dihedral(c[beg][0], c[beg + 1][0], c[beg + 2][0], c[beg + 3][0]);
			e += square(std::sin(0.5*(d - m_bond_dihedral_std[0])));
		}
		//if (beg + 2 < c.size()) {
		//	d = geom::dihedral(c[beg][1], c[beg + 1][0], c[beg + 1][1], c[beg + 2][0]);
		//	e += square(std::sin(0.5*(d - m_bond_dihedral_std[2])));
		//}
		if (beg + 1 < c.size()) {
			d = geom::dihedral(c[beg][0], c[beg][1], c[beg + 1][0], c[beg + 1][1]);
			e += square(std::sin(0.5*(d - m_bond_dihedral_std[3])));
		}
		return e;
	}

	inline double Score6p::en_crash(const Residue &r1, const Residue &r2) {
		int i, j;
		double d, e;
		const Residue *p1, *p2;
		Residue temp1, temp2;

		if (m_cg->is_cg(r1) && m_cg->is_cg(r2)) {
			p1 = &r1;
			p2 = &r2;
		}
		else {
			temp1 = m_cg->to_cg(r1);
			temp2 = m_cg->to_cg(r2);
			p1 = &temp1;
			p2 = &temp2;
		}


		e = 0;
		for (i = 0; i < m_res_size; i++) {
			for (j = 0; j < m_res_size; j++) {
				d = geom::distance(p1->at(i), p2->at(j));
				if (d < 3) {
					e += square(d - 3);
				}
			}
		}

		return e;
	}

	bool Score6p::in_base(int n) {
		return n > 2;
	}

} // namespace jian