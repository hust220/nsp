#pragma once

#include "../utils/Debug.hpp"
#include "../pdb.hpp"
#include "../nuc2d.hpp"
#include "../dg.hpp"
#include "../cg.hpp"
#include "../geom.hpp"
#include "HelixPar.hpp"

BEGIN_JN

Chain *build_chain_dg(Str seq, Str ss);

void dg_dist_init(Mat &dist, int l);

template<typename _Ls>
void dg_dist_read_brokens(Mat &dist, _Ls &&brokens) {
	Int l = dist.rows();
	for (auto && i : brokens) {
		if (i + 1 < l) {
			dist(i, i + 1) = 999;
		}
	}
}

void dg_dist_read_loop(Mat &dist, const SSE &sse);

void dg_dist_read_helix(Mat &dist, const SSE &sse);

void dg_dist_read_ss(Mat &dist, Str seq, Str ss);

void dg_dist_read_chain(Mat &dist, const Chain &chain);

template<typename _Ls>
void dg_dist_read_chain(Mat &dist, const Chain &c, _Ls &&indices) {
	Num d;
	int i, j, l, a, b;
	l = size(c);
	for (i = 0; i < l; i++) {
		a = indices[i];
		if (a == -1) continue;
		for (j = i + 1; j < l; j++) {
			b = indices[j];
			if (b == -1) continue;
			if (c[i].empty() || c[j].empty()) continue;
			d = geom::distance(c[i]["C4*"], c[j]["C4*"]);
			dist(a, b) = d;
			dist(b, a) = d;
		}
	}
}

class BuildLoopDG {
public:
	Eigen::MatrixXd _dist_bound;
	DihBound _dih_bound;

	DG m_dg;
	HelixPar helix_par;
	std::shared_ptr<CG> m_cg;

	BuildLoopDG();

	BuildLoopDG &init(const Str &seq, const Str &ss);

	BuildLoopDG &init(const Chain &c, const std::vector<int> &brokens);

	Chain operator ()();

};

END_JN

