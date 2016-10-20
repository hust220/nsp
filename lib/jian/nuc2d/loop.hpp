#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <list>
#include <type_traits>
#include "../utils/log.hpp"
#include "helix.hpp"
#include "../pp.hpp"

#define LOOP_EACH(l, c) do{\
    int N_RES = 0;\
    res *RES = (l)->head;\
    for (; RES != NULL; RES = RES->next) {\
        c;\
        N_RES++;\
    }\
    N_RES;\
}while(0)

#define LOOP_TRAVERSE(l, c) do{\
    using type = std::remove_reference_t<decltype(l)>;\
    std::list<type> ls;\
    type L = l;\
    while (true) {\
        if (L == NULL) break;\
        ls.push_back(L);\
        c; \
        if (L->son != NULL) { L = L->son;\
        } else if (L->brother != NULL) {\
            L = L->brother;\
            ls.pop_back(); \
        } else {\
            while (true) {\
                ls.pop_back();\
                if (!ls.empty()) {\
                    L = ls.back()->brother;\
                    if (L == NULL) continue;\
                    else { ls.pop_back(); break; }\
                } else { L = NULL; break; }\
            }\
        }\
    }\
}while(0)

namespace jian {

class loop {
public:
    res *head = NULL;
    helix s;
    std::vector<std::pair<int, int>> hinges;
    loop *son = NULL;
    loop *brother = NULL;

    loop() = default;

    loop(const loop &l) {
        head = res::copy(l.head);
        s = l.s;
        hinges = l.hinges;
    }

    ~loop() {
        res::del(head);
    }

    bool has_helix() const { return s.head != NULL; }
    bool has_loop() const { return head != NULL; }
    bool has_son() const { return son != NULL; }
    bool has_brother() const { return brother != NULL; }
    bool empty() const { return head == NULL; }

    res &at(int n) {
        int index = 0;
        for (auto r = head; r != NULL; r = r->next) {
            if (r->type == '&') continue;
            if (index == n) return *r;
            index++;
        }
        throw "jian::at error! out of range.";
    }

    const res &at(int n) const {
        int index = 0;
        for (auto r = head; r != NULL; r = r->next) {
            if (r->type == '&') continue;
            if (index == n) return *r;
            index++;
        }
        throw "jian::at error! out of range.";
    }

    void push_front(res *r) {
        if (head == NULL) head = r;
        else { head->prev = r; r->next = head; head = r; }
    }

    void push_back(res *r) {
        if (head == NULL) head = r;
        else LOOP_EACH(this, if (RES->next == NULL) { RES->next = r; r->prev = RES; break;});
    }

    void del(res *r) {
        if (r->prev != NULL) r->prev->next = r->next; 
        if (r->next != NULL) r->next->prev = r->prev;
        delete r;
    }

    void insert_after(res *p, res *r) {
        r->prev = p; r->next = p->next; p->next = r; 
        if (p->next != NULL) r->next->prev = r;
    }

    void insert_before(res *p, res *r) {
        r->prev = p->prev; r->next = p; p->prev = r; 
        if (p->prev != NULL) r->prev->next = r;
    }

    void del_head() {
        if (head != NULL) {
            res *r = head; head = head->next; delete r;
            if (head != NULL) head->prev = NULL;
        }
    }

    void del_tail() {
        LOOP_EACH(this, if (RES->next == NULL) {RES->prev->next = NULL; delete RES; break;});
    }

    loop *copy(loop *l) {
        if (l == NULL) { return NULL; }

        loop *p = new loop(*l);
        p->son = copy(l->son);
        p->brother = copy(l->brother);
        return p;
    }

    void del(loop *l) {
        if (l == NULL) { return; }

        del(l->son);
        del(l->brother);
        delete l;
        return;
    }

    int len() const {
        int l = 0;
		LOOP_EACH(this, if(RES->type != '&')l++);
		return l;
    }

    int num_branches() const {
        int len = 0; LOOP_EACH(this, if (RES->type == ')') len++); return len / 2;
    }

    int num_sons() const {
        return hinges.size();
    }

    bool is_open() const {
        return num_branches() == num_sons();
    }

    std::string ss() const {
        std::string ss; LOOP_EACH(this, ss += RES->type); return ss;
    }

    std::string seq() const {
		std::string seq;
		LOOP_EACH(this, if (RES->type != '&') seq += RES->name);
		return seq;
    }

//    friend std::ostream &operator <<(std::ostream &out, const loop &l) {
//        out << "Loop: " << l.seq() << ' ' << l.ss() << ' ';
//        LOOP_EACH(&l, out << RES->num << ' ');
//    }

    operator std::string() const {
        std::ostringstream stream;
        stream << seq() << ' ' << ss() << ' ';
        LOOP_EACH(this, stream << RES->num << ' ');
        return stream.str();
    }

    void print() const {
        LOGI << "Hairpin (" << this << ", son: " << son << ", brother: " << brother << ")" << std::endl;
        if (has_loop()) {
            LOGI << "  Loop: " << seq() << ' ' << ss() << ' ';
            LOOP_EACH(this,
                LOGI << RES->num << ' ';
            );
            LOGI << std::endl;
        }
        if (has_helix()) {
            LOGI << "  Helix: " << s.seq() << ' ' << s.ss() << ' ';
            std::list<int> d; 
            HELIX_EACH(s,
                d.insert(std::next(d.begin(), N_BP), {BP->res1.num, BP->res2.num});
            );
			for (auto && i : d) LOGI << i << ' ';
            //EACH(i, d, LOGI << i << ' ');
            LOGI << std::endl;
        }
        LOGI << std::endl;
    }

    void print_tree() const {
        LOOP_TRAVERSE(this, L->print());
    }

};

} // namespace jian

