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

#define BEGIN_LOOP_EACH(l) \
	do {\
		int N_RES = 0;\
		res *RES = (l)->head;\
		for (; RES != NULL; RES = RES->next, N_RES++) 

#define END_LOOP_EACH \
	} while (0)

#define LOOP_EACH(l, c) \
	BEGIN_LOOP_EACH(l) { \
		c; \
	} END_LOOP_EACH

#define BEGIN_LOOP_TRAVERSE(l) \
	do {\
		using type = std::remove_reference_t<decltype(l)>;\
		std::list<type> _path;\
		type _l = l;\
		JN_BREAK_INIT;\
		while (true) {\
			if (_l == NULL) break;\
			_path.push_back(_l);\
			do

#define END_LOOP_TRAVERSE \
			while (0);\
			if (is_break) break;\
			if (_l->son != NULL) { _l = _l->son;\
			} else if (_l->brother != NULL) {\
				_l = _l->brother;\
				_path.pop_back(); \
			} else {\
				while (true) {\
					_path.pop_back();\
					if (!_path.empty()) {\
						_l = _path.back()->brother;\
						if (_l == NULL) continue;\
						else { _path.pop_back(); break; }\
					} else { _l = NULL; break; }\
				}\
			}\
		}\
	} while(0)

#define LOOP_TRAVERSE(l, c) do{\
    using type = std::remove_reference_t<decltype(l)>;\
    L<type> _path;\
    type _l = l;\
    while (true) {\
        if (_l == NULL) break;\
        _path.push_back(_l);\
        c; \
        if (_l->son != NULL) { _l = _l->son;\
        } else if (_l->brother != NULL) {\
            _l = _l->brother;\
            _path.pop_back(); \
        } else {\
            while (true) {\
                _path.pop_back();\
                if (!_path.empty()) {\
                    _l = _path.back()->brother;\
                    if (_l == NULL) continue;\
                    else { _path.pop_back(); break; }\
                } else { _l = NULL; break; }\
            }\
        }\
    }\
}while(0)

BEGIN_JN

class Hairpin {
public:
	using Path = L<Hairpin *>;

	struct TreePathNode {
		Hairpin *val;
		Path path;
	};

	using TreePath = L<TreePathNode>;

    res *head = NULL;
    helix s;
    V<Pair<int, int>> hinges;
    Hairpin *son = NULL;
    Hairpin *brother = NULL;

    Hairpin() = default;

    Hairpin(const Hairpin &l) {
        head = res::copy(l.head);
        s = l.s;
        hinges = l.hinges;
    }

    ~Hairpin() {
        res::del(head);
    }

	template<typename _Fn>
	void traverse(_Fn &&fn) const {
		Path ls;
		Hairpin *L = m_head;
		while (true) {
			if (L == NULL) break;
			ls.push_back(L);
			if (!fn(L, ls)) return;
			if (L->son != NULL) {
				L = L->son;
			}
			else if (L->brother != NULL) {
				L = L->brother;
				ls.pop_back();
			}
			else {
				while (true) {
					ls.pop_back();
					if (!ls.empty()) {
						L = ls.back()->brother;
						if (L == NULL) continue;
						else { ls.pop_back(); break; }
					}
					else { L = NULL; break; }
				}
			}
		}
	}

	TreePath tree_path() const {
		TreePath path;
		traverse([&path](Hairpin *_l, const Path &_path) {
			path.push_back(TreePathNode{_l, _path});
			return false;
		});
		return path;
	}

	bool has(int i) const {
		BEGIN_LOOP_EACH(this) {
			if (RES->type != '(' && RES->type != ')' && RES->num == i) return true;
		} END_LOOP_EACH;
		BEGIN_HELIX_EACH(s) {
			if (BP->res1.num == i || BP->res2.num == i) return true;
		} END_HELIX_EACH;
		return false;
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
		else {
			BEGIN_LOOP_EACH(this) {
				if (RES->next == NULL) {
					RES->next = r; 
					r->prev = RES; 
					break; 
				}
			} END_LOOP_EACH;
		}
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

    Hairpin *copy(Hairpin *l) {
        if (l == NULL) { return NULL; }

        Hairpin *p = new Hairpin(*l);
        p->son = copy(l->son);
        p->brother = copy(l->brother);
        return p;
    }

    void del(Hairpin *l) {
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

	bool is_ml() const {
		return num_sons() > 1;
	}

    bool is_open() const {
        return num_branches() == num_sons();
    }

    Str ss() const {
        Str ss; LOOP_EACH(this, ss += RES->type); return ss;
    }

    Str seq() const {
		Str seq;
		LOOP_EACH(this, if (RES->type != '&') seq += RES->name);
		return seq;
    }

//    friend std::ostream &operator <<(std::ostream &out, const Hairpin &l) {
//        out << "Loop: " << l.seq() << ' ' << l.ss() << ' ';
//        LOOP_EACH(&l, out << RES->num << ' ');
//    }

    operator Str() const {
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
        LOOP_TRAVERSE(this, _l->print());
    }

};

END_JN

