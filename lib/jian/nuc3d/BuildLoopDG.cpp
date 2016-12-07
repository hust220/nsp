#include <memory>
#include "../cg.hpp"
#include "../geom.hpp"
#include "BuildLoopDG.hpp"

BEGIN_JN

	BuildLoopDG::BuildLoopDG() {
		m_cg.reset(CG::fac_t::create("1p"));
	}

	Chain *build_chain_dg(S seq, S ss) {
		static BuildLoopDG builder;
		builder.init(seq, ss);
		Chain *chain = new Chain(std::move(builder()));
		return chain;
	}

	void BuildLoopDG::init_dist_bound(Mat &b, int l) {
		b.resize(l, l);
		for (int i = 0; i < l; i++) for (int j = i; j < l; j++) {
			if (j - i == 1) {
				b(i, j) = helix_par.dist_bond;
				b(j, i) = helix_par.dist_bond;
			}
			else if (i != j) {
				b(j, i) = 5;
				b(i, j) = 999;
			}
			else {
				b(i, j) = 0;
			}
		}
	}


	BuildLoopDG &BuildLoopDG::init(const S &seq, const S &ss) {
		int len = seq.size(); 
		init_dist_bound(_dist_bound, len);
		SSTree ss_tree;
		ss_tree.make(seq, ss);
		for (auto &&sse : ss_tree) {
			set_bound_loop(_dist_bound, _dih_bound, &sse);
			set_bound_helix(_dist_bound, _dih_bound, sse.helix);
		}
		return *this;
	}

	BuildLoopDG &BuildLoopDG::init(const Chain &c, const std::vector<int> &brokens) {
		int i, j, l;
		Num d;

		l = c.size();
		init_dist_bound(_dist_bound, l);
		for (auto && i : brokens) {
			if (i + 1 < l) {
				_dist_bound(i, i + 1) = 999;
			}
		}
		for (i = 0; i < l; i++) {
			for (j = i + 1; j < l; j++) {
				if (c[i].empty() || c[j].empty()) continue;
				d = geom::distance(c[i]["C4*"], c[j]["C4*"]);
				_dist_bound(i, j) = d;
				_dist_bound(j, i) = d;
			}
		}
		return *this;
	}

	Chain BuildLoopDG::operator ()() {
		LOG << "## Build Loop By DG" << std::endl;
		Mat &&c = m_dg(_dist_bound);
		return m_cg->to_aa(c, 0, c.rows() - 1);
	}

	void BuildLoopDG::set_bound_loop(Mat &b, DihBound &d, SSE *l) {
		if (!l->loop.empty()) {
			auto it1 = l->loop.begin();
			auto it2 = STD_ next(it1);
			for (; it2 != l->loop.end(); it1++, it2++) {
				if (it1->type == '(' && it2->type == ')') {
					b(it1->num - 1, it2->num - 1) = b(it2->num - 1, it1->num - 1) = helix_par.dist_bp;
				}
				else {
					b(it1->num - 1, it2->num - 1) = b(it2->num - 1, it1->num - 1) = helix_par.dist_bond;
				}
			}
		}
	}

	void BuildLoopDG::set_bound_helix(Mat &b, DihBound &d, const Helix &h) {
		int len = 0;
		std::deque<int> s1, s2;
		for (auto && bp : h) {
			len++;
			s1.push_back(bp.res1.num - 1);
			s2.push_back(bp.res2.num - 1);
		}
		for (int i = 0; i < len; i++) for (int j = i + 1; j < len; j++) {
			b(s1[i], s1[j]) = b(s1[j], s1[i]) = helix_par.dist_a(j - i);
			b(s2[i], s2[j]) = b(s2[j], s2[i]) = helix_par.dist_b(j - i);
			b(s1[i], s2[j]) = b(s2[j], s1[i]) = helix_par.dist_c(j - i);
			b(s1[j], s2[i]) = b(s2[i], s1[j]) = helix_par.dist_d(j - i);
		}
		for (int i = 0; i < len; i++) b(s1[i], s2[i]) = b(s2[i], s1[i]) = helix_par.dist_bp;
	}

END_JN

