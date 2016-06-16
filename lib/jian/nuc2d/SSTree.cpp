#include <iostream>
#include <string>
#include <regex>
#include <vector>
#include "SSTree.hpp"
#include "loop.hpp"
#include "NucSS.hpp"
#include "../pdb.hpp"

namespace jian {

void free_ss_tree(loop *l) {
    if (l != NULL) {
        free_ss_tree(l->son);
        free_ss_tree(l->brother);
        delete l;
    }
}

struct SSTreeImpl {
    loop *head = NULL;

    ~SSTreeImpl() {
        if (head != NULL) free_ss_tree(head);
    }
};

void set_tree_relation(std::vector<loop *> &s, loop *l) {
    int num = l->num_sons();
    int i = s.size() - num;
    if (s.size() != 0 && num != 0) { l->son = s[i]; }
    for (; i < (int) s.size(); i++) {
        if (i + 1 == s.size()) s[i]->brother = NULL; else s[i]->brother = s[i + 1];
    }
    FOR((k, num), s.pop_back());
    s.push_back(l);
}

int char_index(const char &c) {
    static std::map<char, int> m {{'.', 0}, {'(', 1}, {')', 2}};
    if (m.find(c) != m.end()) return m[c];
    else return 0;
}

bool find_hairpin_position(const std::deque<res> &v, int &left, int &right) {
    std::vector<std::vector<int>> t {{0, 1, -1}, {2, 1, 3}, {2, 1, 3}, {0, 1, 3}};
    int state = 0, new_state, num_hinge = 0;
    EACH((r, i), v,
        new_state = t[state][char_index(r.type)];
        if (new_state == -1) throw "Wrong secondary structure!";
        if (state == 1 && (new_state == 2 || new_state == 3)) left = i - 1;
        if ((state == 1 || state == 2) && new_state == 3) right = i;
        if (new_state == 3) return true;
        state = new_state;
    );
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
    if (right - left != 1) FOR((i, left-hinge+1, right+hinge), l->push_back(new res(v[i])));
    FOR((i, len), l->s.push_back(new bp(v[left-len+1+i], v[right+len-1-i])));
    v.erase(v.begin()+left-len+1+hinge, v.begin()+right+len-hinge);
    FOR((i, hinge), v[left-len+1+i].type = 'Z'); FOR((i, hinge, 2*hinge), v[left-len+1+i].type = 'z');
    return l;
}

bool find_hairpin(std::deque<res> &v, std::vector<loop *> &ls, int hinge) {
    if (v.empty()) return false;
    int left, right; loop *l;
    if (find_hairpin_position(v, left, right)) {
        int len = len_helix(v, left, right);
        if (len < hinge) {
            FOR((i, left-len+1, right+len), if (HAS(('(', ')'), v[i].type)) v[i].type = '.';);
            return true;
        } else {
            l = dig_hairpin(v, left, right, len, hinge);
        }
    } else {
        l = new loop; EACH(i, v, l->push_back(new res(i)));
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
    std::vector<int> v(ss.size()); int f = 1; EACH((c, i), ss, IF(c != '&', v[i] = f; f++, v[i] = -1));
    LOOP_TRAVERSE(l, 
        if (L->has_loop())
            LOOP_EACH(L, RES->num = v[RES->num - 1]; if (RES->num != -1) RES->name = seq[RES->num - 1]);
        if (L->has_helix()) 
            HELIX_EACH(L->s, BP->res1.num = v[BP->res1.num - 1]; BP->res1.name = seq[BP->res1.num - 1];
                             BP->res2.num = v[BP->res2.num - 1]; BP->res2.name = seq[BP->res2.num - 1]);
    );
}

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
    Debug::print("## Make Secondary Structure Tree\n");
    std::cout << seq << std::endl;
    std::cout << ss << std::endl;
    if (NucSS::seq_match_ss(seq, ss)) {
        _impl->head = set_tree(ss, hinge);
        read_seq(_impl->head, seq, ss);
    } else {
        throw "The sequence and the secondary structure don't match!";
    }
}

} // namespace jian

