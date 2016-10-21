#include <iostream>
#include <string>
#include <regex>
#include <vector>
#include "SSTree.hpp"
#include "loop.hpp"
#include "NucSS.hpp"
#include "../pdb.hpp"
#include "../utils/log.hpp"

namespace jian {

std::pair<int, int> loop_head_tail(loop *l) {
    int left, right;
    if (l != NULL) {
        if (l->has_helix()) {
            left = l->s.head->res1.num-1;
            right = l->s.head->res2.num-1;
            return {left, right};
        } else if (l->has_loop()) {
            left = l->head->num-1;
            LOOP_EACH(l,
                if (RES->next == NULL) {
                    right = RES->num-1;
                }
            );
            return {left, right};
        }
    } else {
        throw "loop_head_tail error!";
    }
}

void ss_read_tree(std::string &ss, loop *l) {
    auto p = loop_head_tail(l);
    ss.resize(p.second - p.first + 1);
    LOOP_TRAVERSE(l,
        HELIX_EACH(L->s,
            ss[BP->res1.num-1] = BP->res1.type;
            ss[BP->res2.num-1] = BP->res2.type;
        );
        LOOP_EACH(L,
            ss[RES->num-1] = RES->type;
        );
    );
}

void seq_read_tree(std::string &seq, loop *l) {
    auto p = loop_head_tail(l);
    seq.resize(p.second - p.first + 1);
    LOOP_TRAVERSE(l,
        HELIX_EACH(L->s,
            seq[BP->res1.num-1] = BP->res1.name;
            seq[BP->res2.num-1] = BP->res2.name;
        );
        LOOP_EACH(L,
            seq[RES->num-1] = RES->name;
        );
    );
}

void print_ss_tree(loop *l) {
    LOOP_TRAVERSE(l, L->print());
}

void free_ss_tree(loop *l) {
    if (l != NULL) {
        free_ss_tree(l->son);
        free_ss_tree(l->brother);
        delete l;
    }
}

namespace sstree_detail {

void set_tree_relation(std::vector<loop *> &s, loop *l) {
    int num = l->num_sons();
    int i = s.size() - num;
    if (s.size() != 0 && num != 0) { l->son = s[i]; }
    for (; i < (int) s.size(); i++) {
        if (i + 1 == s.size()) s[i]->brother = NULL; else s[i]->brother = s[i + 1];
    }
	for (int k = 0; k < num; k++) s.pop_back();
    //FOR((k, num), s.pop_back());
    s.push_back(l);
}

int char_index(const char &c) {
    if (c == '(') {
        return 1;
    } else if (c == ')') {
        return 2;
    } else {
        return 0;
    }
}

bool find_hairpin_position(const std::deque<res> &v, int &left, int &right) {
    std::vector<std::vector<int>> t {{0, 1, -1}, {2, 1, 3}, {2, 1, 3}, {0, 1, 3}};
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
    while (left-len>=0 && v[left-len].type == '(' && 
           right+len < v.size() && v[right+len].type == ')') len++; 
    return len;
}

loop *dig_hairpin(std::deque<res> &v, int left, int right, int len, int hinge) {
    loop *l = new loop;
	if (right - left != 1) for (int i = left - hinge + 1; i < right + hinge; i++) {
		l->push_back(new res(v[i]));
	}
    //if (right - left != 1) FOR((i, left-hinge+1, right+hinge), l->push_back(new res(v[i])));
	for (int i = 0; i < len; i++) l->s.push_back(new bp(v[left - len + 1 + i], v[right + len - 1 - i]));
    //FOR((i, len), l->s.push_back(new bp(v[left-len+1+i], v[right+len-1-i])));
    v.erase(v.begin()+left-len+1+hinge, v.begin()+right+len-hinge);
	for (int i = 0; i < hinge; i++) v[left - len + 1 + i].type = 'Z';
	for (int i = hinge; i < 2 * hinge; i++) v[left - len + 1 + i].type = 'z';
    //FOR((i, hinge), v[left-len+1+i].type = 'Z'); FOR((i, hinge, 2*hinge), v[left-len+1+i].type = 'z');
    return l;
}

bool find_hairpin(std::deque<res> &v, std::vector<loop *> &ls, int hinge) {
    int left, right;
    loop *l;

    if (v.empty()) return false;

    if (find_hairpin_position(v, left, right)) {
        int len = len_helix(v, left, right);
        if (len < hinge) {
			for (int i = left - len + 1; i < right + len; i++) if (v[i].type == '(' || v[i].type == ')') {
				v[i].type = '.';
			}
            //FOR((i, left-len+1, right+len), if (HAS(('(', ')'), v[i].type)) v[i].type = '.';);
            return true;
        } else {
            l = dig_hairpin(v, left, right, len, hinge);
        }
    } else {
        l = new loop;
		for (auto && i : v) l->push_back(new res(i));
        //EACH(i, v, l->push_back(new res(i)));
        v.clear();
        if (std::regex_match(l->ss(), std::regex("Z+z+"))) return false; 
    }
    LOOP_EACH(l, if (RES->type == 'Z' && RES->next->type == 'z') l->hinges.push_back({N_RES, N_RES+1}));
    LOOP_EACH(l, if (RES->type == 'Z') RES->type = '('; else if (RES->type == 'z') RES->type = ')');
    set_tree_relation(ls, l);
    return true;
}

loop *set_tree(const std::string &ss, int hinge) {
    std::deque<res> v; 
    int i = 0;
    for (auto && c : ss) {
        v.push_back(res(c, i+1));
        i++;
    }
    std::vector<loop *> ls; 
    while (find_hairpin(v, ls, hinge));
    return ls.back();
}

void read_seq(loop *l, const std::string &seq, const std::string &ss) {
    std::vector<int> v(ss.size());
    int f = 1;
	int i = 0;
	for (auto && c : ss) {
		if (c != '&') {
			v[i] = f;
			f++;
		} else {
			v[i] = -1;
		}
		i++;
	}
    //EACH((c, i), ss, IF(c != '&', v[i] = f; f++, v[i] = -1));
    LOOP_TRAVERSE(l, 
        if (L->has_loop()) {
            LOOP_EACH(L,
                RES->num = v[RES->num - 1];
                if (RES->num != -1) {
                    RES->name = seq[RES->num - 1];
                }
            );
        }
        if (L->has_helix()) {
            HELIX_EACH(L->s,
                BP->res1.num = v[BP->res1.num - 1];
                BP->res1.name = seq[BP->res1.num - 1];
                BP->res2.num = v[BP->res2.num - 1];
                BP->res2.name = seq[BP->res2.num - 1]
            );
        }
    );
}

} // namespace sstree_detail

struct SSTreeImpl {
    loop *head = NULL;

    ~SSTreeImpl() {
        if (head != NULL) free_ss_tree(head);
    }
};

SSTree::SSTree() : _impl(new SSTreeImpl) {}

SSTree::~SSTree() {
    delete _impl;
}

loop *&SSTree::head() {
    return _impl->head;
}

bool SSTree::empty() const {
    return _impl->head == NULL;
}

void SSTree::make(const std::string &seq, const std::string &ss, int hinge) {
    LOG << "## Make secondary structure tree with no broken tag" << std::endl;
    LOG << seq << std::endl;
    LOG << ss << std::endl;
    if (NucSS::seq_match_ss(seq, ss)) {
        std::string ss_nbt;
        std::copy_if(ss.begin(), ss.end(), std::back_inserter(ss_nbt), [](auto &&c){return c!='&';});
        _impl->head = sstree_detail::set_tree(ss_nbt, hinge);
        sstree_detail::read_seq(_impl->head, seq, ss_nbt);
    } else {
        LOG << "Error:" << std::endl;
        LOG << "Sequence: " << seq << std::endl;
        LOG << "SS: " << ss << std::endl;
        LOG << "The lengh of the sequence should equal to the length of the secondary structure!" << std::endl;
        throw "SSTree::make error!";
    }
}

void SSTree::make_b(const std::string &seq, const std::string &ss, int hinge) {
    LOG << "## Make secondary structure tree with broken tag" << std::endl;
    LOG << seq << std::endl;
    LOG << ss << std::endl;
    if (NucSS::seq_match_ss(seq, ss)) {
        _impl->head = sstree_detail::set_tree(ss, hinge);
        sstree_detail::read_seq(_impl->head, seq, ss);
    } else {
        LOG << "Error:" << std::endl;
        LOG << "Sequence: " << seq << std::endl;
        LOG << "SS: " << ss << std::endl;
        LOG << "The sequence and the secondary structure don't match!" << std::endl;
        throw "SSTree::make_b error!";
    }
}

loop *ss_tree(std::string seq, std::string ss, int hinge) {
    loop *tree = sstree_detail::set_tree(ss, hinge);
    sstree_detail::read_seq(tree, seq, ss);
    return tree;
}

} // namespace jian

