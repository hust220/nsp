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
		BEGIN_LOOP_TRAVERSE(ss_tree.head()) {
			set_bound_loop(_dist_bound, _dih_bound, _l);
			set_bound_helix(_dist_bound, _dih_bound, _l->s);
		} END_LOOP_TRAVERSE;
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

	void BuildLoopDG::set_bound_loop(Mat &b, DihBound &d, Hairpin *l) {
		BEGIN_LOOP_EACH(l) {
			if (RES->next != NULL) {
				if (RES->type == '(' && RES->next->type == ')') {
					b(RES->num - 1, RES->next->num - 1) = b(RES->next->num - 1, RES->num - 1) = helix_par.dist_bp;
				}
				else {
					b(RES->num - 1, RES->next->num - 1) = b(RES->next->num - 1, RES->num - 1) = helix_par.dist_bond;
				}
			}
		} END_LOOP_EACH;
	}

	void BuildLoopDG::set_bound_helix(Mat &b, DihBound &d, const helix &h) {
		int len = 0;
		std::deque<int> s1, s2;
		BEGIN_HELIX_EACH(h) {
			len++;
			s1.push_back(BP->res1.num - 1);
			s2.push_back(BP->res2.num - 1);
		} END_HELIX_EACH;
		for (int i = 0; i < len; i++) for (int j = i + 1; j < len; j++) {
			b(s1[i], s1[j]) = b(s1[j], s1[i]) = helix_par.dist_a(j - i);
			b(s2[i], s2[j]) = b(s2[j], s2[i]) = helix_par.dist_b(j - i);
			b(s1[i], s2[j]) = b(s2[j], s1[i]) = helix_par.dist_c(j - i);
			b(s1[j], s2[i]) = b(s2[i], s1[j]) = helix_par.dist_d(j - i);
		}
		for (int i = 0; i < len; i++) b(s1[i], s2[i]) = b(s2[i], s1[i]) = helix_par.dist_bp;
	}

END_JN

