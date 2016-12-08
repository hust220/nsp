#include <iostream>
#include <string>
#include <regex>
#include <vector>
#include <algorithm>
#include "SSTree.hpp"
#include "SSE.hpp"
#include "NASS.hpp"
#include "../pdb.hpp"
#include "../utils/log.hpp"

#define TRUE_LENGTH(a) \
	std::count_if(a.begin(), a.end(), [](auto && c){return c != '&';})

BEGIN_JN

namespace sstree_detail {

	void set_tree_relation(V<SSTree::El *> &s, SSTree::El *l) {
		int num = l->data.num_sons();
		//LOG << num << Endl;
		int i = size(s) - num;
		//LOG << i << Endl;
		if (size(s) != 0 && num != 0) { l->son = s[i]; }
		for (; i < (int)s.size(); i++) {
			if (i + 1 == s.size()) s[i]->brother = NULL; else s[i]->brother = s[i + 1];
		}
		for (int k = 0; k < num; k++) s.pop_back();
		s.push_back(l);
	}

	int char_index(const char &c) {
		if (c == '(') {
			return 1;
		}
		else if (c == ')') {
			return 2;
		}
		else {
			return 0;
		}
	}

	bool find_hairpin_position(const std::deque<res> &v, int &left, int &right) {
		std::vector<std::vector<int>> t{ {0, 1, -1}, {2, 1, 3}, {2, 1, 3}, {0, 1, 3} };
		int state = 0, new_state, num_hinge = 0;
		int i = 0;
		for (auto && r : v) {
			//EACH((r, i), v,
			new_state = t[state][char_index(r.type)];
			if (new_state == -1) throw "Wrong secondary structure!";
			if (state == 1 && (new_state == 2 || new_state == 3)) left = i - 1;
			if ((state == 1 || state == 2) && new_state == 3) right = i;
			if (new_state == 3) return true;
			state = new_state;
			i++;
		}
		return false;
	}

	int len_helix(const std::deque<res> &v, int left, int right) {
		int len = 1;
		while (left >= len && v[left - len].type == '(' &&
			right + len < v.size() && v[right + len].type == ')') {
			len++;
		}
		return len;
	}

	SSE *dig_hairpin(std::deque<res> &v, int left, int right, int len, int hinge) {
		SSE *l = new SSE;
		if (right - left != 1) for (int i = left - hinge + 1; i < right + hinge; i++) {
			l->loop.push_back(res(v[i]));
		}
		for (int i = 0; i < len; i++) l->helix.push_back(bp(v[left - len + 1 + i], v[right + len - 1 - i]));
		v.erase(v.begin() + left - len + 1 + hinge, v.begin() + right + len - hinge);
		for (int i = 0; i < hinge; i++) v[left - len + 1 + i].type = 'Z';
		for (int i = hinge; i < 2 * hinge; i++) v[left - len + 1 + i].type = 'z';
		return l;
	}

	bool find_hairpin(Q<res> &v, V<SSTree::El *> &ls, int hinge) {
		int left, right;
		SSE *l;

		if (v.empty()) return false;

		if (find_hairpin_position(v, left, right)) {
			int len = len_helix(v, left, right);
			if (len < hinge) {
				for (int i = left - len + 1; i < right + len; i++) if (v[i].type == '(' || v[i].type == ')') {
					v[i].type = '.';
				}
				return true;
			}
			else {
				l = dig_hairpin(v, left, right, len, hinge);
			}
		}
		else {
			l = new SSE;
			for (auto && i : v) l->loop.push_back(res(i));
			v.clear();
			if (std::regex_match(l->loop.ss(), std::regex("Z+z+"))) return false;
		}

		int n_res = 0;
		for (auto it = l->loop.begin(); it != l->loop.end() /*&&*/ /*STD_ next(it) != l->loop.end()*/; it++) {
			if (it->type == 'Z' && STD_ next(it)->type == 'z') {
				l->hinges.push_back({ n_res, n_res + 1 });
			}
			n_res++;
		}
		//for (auto &&res : l->loop) {
		//	if (res.type == 'Z' && res.next->type == 'z') {
		//		l->hinges.push_back({ n_res, n_res + 1 });
		//	}
		//	n_res++;
		//}

		for (auto && res : l->loop) {
			if (res.type == 'Z') res.type = '('; else if (res.type == 'z') res.type = ')';
		}

		set_tree_relation(ls, SSTree::El::make(STD_ move(*l)));
		return true;
	}

	SSTree::El *set_tree(const Str &ss, int hinge) {
		std::deque<res> v;
		int i = 0;
		for (auto && c : ss) {
			v.push_back(res(c, i + 1));
			i++;
		}
		std::vector<SSTree::El *> ls;
		while (find_hairpin(v, ls, hinge));
		//STD_ cout << size(ls) << STD_ endl;
		return ls.back();
	}

	void read_seq(SSTree::El *l, const Str &seq, const Str &ss) {
		std::vector<int> v(ss.size());
		int f = 1;
		int i = 0;
		for (auto && c : ss) {
			if (c != '&') {
				v[i] = f;
				f++;
			}
			else {
				v[i] = -1;
			}
			i++;
		}
		for (auto && sse : SSTree::make_range(l)) {
			if (sse.has_loop()) {
				for (auto && res : sse.loop) {
					res.num = v[res.num - 1];
					if (res.num != -1) {
						res.name = seq[res.num - 1];
					}
				}
			}
			if (sse.has_helix()) {
				for (auto && bp : sse.helix) {
					bp.res1.num = v[bp.res1.num - 1];
					bp.res1.name = seq[bp.res1.num - 1];
					bp.res2.num = v[bp.res2.num - 1];
					bp.res2.name = seq[bp.res2.num - 1];
				}
			}
		}
	}

} // namespace sstree_detail

SSTree::SSTree(const Str &seq, const Str &ss, int hinge) {
	make(seq, ss, hinge);
}

void SSTree::make(const Str &seq, const Str &ss, int hinge) {
	free();
	//LOG << "## Make secondary structure tree with no broken tag" << std::endl;
	//LOG << seq << std::endl;
	//LOG << ss << std::endl;
	if (NASS::seq_match_ss(seq, ss)) {
		Str ss_nbt;
		std::copy_if(ss.begin(), ss.end(), std::back_inserter(ss_nbt), [](auto &&c) {return c != '&'; });
		m_beg.el = sstree_detail::set_tree(ss_nbt, hinge);
		sstree_detail::read_seq(m_beg.el, seq, ss_nbt);
	}
	else {
		LOG << "Error:" << std::endl;
		LOG << "Sequence: " << seq << " (" << TRUE_LENGTH(seq) << " nt)" << std::endl;
		LOG << "SS: " << ss << " (" << TRUE_LENGTH(ss) << " nt)" << std::endl;
		LOG << "The lengh of the sequence should equal to the length of the secondary structure!" << std::endl;
		throw "SSTree::make error!";
	}
}

void SSTree::make_b(const Str &seq, const Str &ss, int hinge) {
	free();
	//LOG << "## Make secondary structure tree with broken tag" << std::endl;
	//LOG << seq << std::endl;
	//LOG << ss << std::endl;
	if (NASS::seq_match_ss(seq, ss)) {
		m_beg.el = sstree_detail::set_tree(ss, hinge);
		sstree_detail::read_seq(m_beg.el, seq, ss);
	}
	else {
		LOG << "Error:" << std::endl;
		LOG << "Sequence: " << seq << std::endl;
		LOG << "SS: " << ss << std::endl;
		LOG << "The sequence and the secondary structure don't match!" << std::endl;
		throw "SSTree::make_b error!";
	}
}

END_JN

