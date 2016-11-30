#pragma once

#include "../geom.hpp"
#include "../pdb.hpp"

namespace jian {

	class ResConf {
	public:
		Residue res;
		RowVec p1;
		RowVec p2;
		Mat pp;

		template<typename _Vec>
		ResConf(const Residue &r, const _Vec &v) {
			res = r;
			const Atom &atom = r["P"];
			p1.resize(3);
			for (int i = 0; i < 3; i++) p1[i] = atom[i];
			p2.resize(3);
			for (int i = 0; i < 3; i++) p2[i] = v[i];
			pp = hstack(p1, p2);
		}

		template<typename _Ls, typename _Residues>
		static void extract(_Ls &&ls, const _Residues &residues) {
			int l = size(residues);
			num_t d;
			for (int i = 0; i < l - 1; i++) {
				if (residues[i].has_atom("P", "C4*") && residues[i+1].has_atom("P")) {
					//LOG << residues[i] << std::endl;
					const Atom &p1 = residues[i]["P"];
					const Atom &c1 = residues[i]["C4*"];
					const Atom &p2 = residues[i + 1]["P"];
					d = geom::distance(c1, p2);
					if (d < 4.0) ls.push_back(ResConf(residues[i], p2));
				}
			}
		}

	};


}