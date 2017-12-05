#pragma once

#include "geom.hpp"
#include "pdb.hpp"

BEGIN_JN

template<int N>
class FragConf {
public:
	using Confs = Deque<FragConf>;
	using MapConfs = Map<Str, Confs>;

    const static int frag_len = N;

	Array<Residue, N> frag;
	RowVec p1;
	RowVec p2;
	Mat pp;

	template<typename _Frag, typename _Vec>
	FragConf(const _Frag &f, const _Vec &v) {
        for (Int i = 0; i < N; i++) frag[i] = f[i];
		const Atom &atom = f[0]["P"];
		p1.resize(3);
		for (int i = 0; i < 3; i++) p1[i] = atom[i];
		p2.resize(3);
		for (int i = 0; i < 3; i++) p2[i] = v[i];
		pp = hstack(p1, p2);
	}

    template<typename _Frag>
    static Str frag_name(const _Frag &f) {
        std::stringstream stream;
        for (auto && r : f) stream << r.name;
        return stream.str();
    }

	template<typename _Residues>
	static void extract(MapConfs &map, const _Residues &rs) {
		int l = size(rs);
		Num d;
        Deque<Residue> dq;
		for (int i = 0; i < l; i++) {
            if (!has_atom(rs[i], "P")) {
                dq.clear();
                continue;
            }
            const Atom &p2 = rs[i]["P"];
            if (i > 0) {
                if (!has_atom(rs[i-1], "C4*")) {
                    dq.clear();
                    continue;
                }
                const Atom &c1 = rs[i-1]["C4*"];
                d = geom::distance(c1, p2);
                if (d > 4.0) dq.clear();
            }
            if (size(dq) == N) {
                map[frag_name(dq)].push_back(FragConf(dq, p2));
                dq.clear();
            }
            dq.push_back(rs[i]);
		}
	}

};

END_JN

