#pragma once

#include "../geom.hpp"
#include "../pdb.hpp"

BEGIN_JN

	class ResConf {
	public:
		using Confs = std::deque<ResConf>;
		using MapConfs = std::map<Str, Confs>;

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

		template<typename _Residues>
		static void extract(MapConfs &map, const _Residues &residues) {
			int l = size(residues);
			Num d;
			for (int i = 0; i < l - 1; i++) {
				if (has_atom(residues[i], "P", "C4*") && has_atom(residues[i+1], "P")) {
					const Residue &r1 = residues[i];
					const Residue &r2 = residues[i+1];
					//LOG << residues[i] << std::endl;
					const Atom &p1 = r1["P"];
					const Atom &c1 = r1["C4*"];
					const Atom &p2 = r2["P"];
					d = geom::distance(c1, p2);
					if (d < 4.0) map[r1.name].push_back(ResConf(r1, p2));
				}
			}
		}

	};


}