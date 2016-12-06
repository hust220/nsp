#pragma once

#include <numeric>
#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <list>
#include <type_traits>
#include "../utils/log.hpp"
#include "Loop.hpp"
#include "helix.hpp"
#include "../pp.hpp"
#include "../utils/Tree.hpp"

BEGIN_JN

class SSE {
public:
	Loop loop;
    Helix helix;
    V<Pair<int, int>> hinges;

    SSE *son = NULL;
    SSE *brother = NULL;

    //SSE() = default;

    //SSE(const SSE &l) {
    //    head = res::copy(l.head);
    //    s = l.s;
    //    hinges = l.hinges;
    //}

    //~SSE() {
    //    res::del(head);
    //}

	Pair<int, int> head_tail() const {
		int left, right;
		if (has_helix()) {
			left = helix.head->res1.num - 1;
			right = helix.head->res2.num - 1;
			return { left, right };
		}
		else if (has_loop()) {
			left = loop.head->num - 1;
			for (auto &&res : loop) {
				if (res.next == NULL) {
					right = res.num - 1;
				}
			}
			return { left, right };
		}
		else {
			throw "SSE head_tail error";
		}

	}
	
	bool has(int i) const {
		return STD_ find_if(loop.begin(), loop.end(), [&i](auto &&res) {
			return res.type != '(' && res.type != ')' && res.num == i;
		}) != loop.end() ||
			STD_ find_if(helix.begin(), helix.end(), [&i](auto &&bp) {
			return bp.res1.num == i || bp.res2.num == i;
		}) != helix.end();
	}

    bool has_helix() const {
		return !helix.empty();
	}

    bool has_loop() const {
		return !loop.empty();
	}

    bool has_son() const {
		return son != NULL;
	}

    bool has_brother() const {
		return brother != NULL; 
	}

    //res &at(int n) {
    //    int index = 0;
    //    for (auto r = head; r != NULL; r = r->next) {
    //        if (r->type == '&') continue;
    //        if (index == n) return *r;
    //        index++;
    //    }
    //    throw "jian::at error! out of range.";
    //}

    //const res &at(int n) const {
    //    int index = 0;
    //    for (auto r = head; r != NULL; r = r->next) {
    //        if (r->type == '&') continue;
    //        if (index == n) return *r;
    //        index++;
    //    }
    //    throw "jian::at error! out of range.";
    //}

  //  void push_front(res *r) {
  //      if (head == NULL) head = r;
  //      else { head->prev = r; r->next = head; head = r; }
  //  }

  //  void push_back(res *r) {
  //      if (head == NULL) head = r;
		//else {
		//	BEGIN_LOOP_EACH(this) {
		//		if (RES->next == NULL) {
		//			RES->next = r; 
		//			r->prev = RES; 
		//			break; 
		//		}
		//	} END_LOOP_EACH;
		//}
  //  }

  //  void del(res *r) {
  //      if (r->prev != NULL) r->prev->next = r->next; 
  //      if (r->next != NULL) r->next->prev = r->prev;
  //      delete r;
  //  }

  //  void insert_after(res *p, res *r) {
  //      r->prev = p; r->next = p->next; p->next = r; 
  //      if (p->next != NULL) r->next->prev = r;
  //  }

  //  void insert_before(res *p, res *r) {
  //      r->prev = p->prev; r->next = p; p->prev = r; 
  //      if (p->prev != NULL) r->prev->next = r;
  //  }

  //  void del_head() {
  //      if (head != NULL) {
  //          res *r = head; head = head->next; delete r;
  //          if (head != NULL) head->prev = NULL;
  //      }
  //  }

  //  void del_tail() {
  //      LOOP_EACH(this, if (RES->next == NULL) {RES->prev->next = NULL; delete RES; break;});
  //  }

  //  SSE *copy(SSE *l) {
  //      if (l == NULL) { return NULL; }

  //      SSE *p = new SSE(*l);
  //      p->son = copy(l->son);
  //      p->brother = copy(l->brother);
  //      return p;
  //  }

  //  void del(SSE *l) {
  //      if (l == NULL) { return; }

  //      del(l->son);
  //      del(l->brother);
  //      delete l;
  //      return;
  //  }

    int num_branches() const {
		return STD_ count_if(loop.begin(), loop.end(), [](auto &&res) {return res.type == ')'; }) / 2;
    }

    int num_sons() const {
        return hinges.size();
    }

    bool is_open() const {
        return num_branches() == num_sons();
    }

	bool is_hp() const {
		return !is_open() && num_sons() == 0;
	}

	bool is_il() const {
		return !is_open() && num_sons() == 1;
	}

	bool is_ml() const {
		return !is_open() && num_sons() > 1;
	}

//    friend std::ostream &operator <<(std::ostream &out, const SSE &l) {
//        out << "Loop: " << l.seq() << ' ' << l.ss() << ' ';
//        LOOP_EACH(&l, out << RES->num << ' ');
//    }

   // operator Str() const {
   //     std::ostringstream stream;
   //     stream << seq() << ' ' << ss() << ' ';
   //     LOOP_EACH(this, stream << RES->num << ' ');
   //     return stream.str();
   // }

   // void print() const {
   //     LOGI << "SSE (" << this << ", son: " << son << ", brother: " << brother << ")" << std::endl;
   //     if (has_loop()) {
   //         LOGI << "  Loop: " << seq() << ' ' << ss() << ' ';
   //         LOOP_EACH(this,
   //             LOGI << RES->num << ' ';
   //         );
   //         LOGI << std::endl;
   //     }
   //     if (has_helix()) {
   //         LOGI << "  Helix: " << s.seq() << ' ' << s.ss() << ' ';
   //         std::list<int> d; 
   //         HELIX_EACH(s,
   //             d.insert(std::next(d.begin(), N_BP), {BP->res1.num, BP->res2.num});
   //         );
			//for (auto && i : d) LOGI << i << ' ';
   //         //EACH(i, d, LOGI << i << ' ');
   //         LOGI << std::endl;
   //     }
   //     LOGI << std::endl;
   // }

};

END_JN

