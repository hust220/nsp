#ifndef LOOP_H
#define LOOP_H

#include "helix.h"

namespace jian {

namespace nuc2d {

class loop {
public:
	loop() {}
	loop(loop *l) : type(l->type), len(l->len), isPartial(l->isPartial), 
	                head(res::copy(l->head)), s(l->s), hinges(l->hinges), 
	                sons(l->sons), son(l->son), brother(l->brother) {
	}
	loop(const loop &l) : type(l.type), len(l.len), isPartial(l.isPartial), 
	                      head(res::copy(l.head)), s(l.s), hinges(l.hinges), 
	                      sons(l.sons), son(l.son), brother(l.brother) {
	}
	loop &operator =(const loop &l) {
		type = l.type;
		len = l.len;
		isPartial = l.isPartial;
		res::del(head);
		head = res::copy(l.head);
		s = l.s;
		hinges = l.hinges;
		sons = l.sons;
		son = l.son;
		brother = l.brother;
		return *this;
	}
	static loop *copy(loop *l) {
		if (l == NULL) {
			return NULL;
		}

		loop *l_ = new loop(l);
		l_->son = copy(l->son);
		l_->brother = copy(l->brother);
		return l_;
	}
	static void del(loop *l) {
		if (l == NULL) {
			return;
		}

		del(l->son);
		del(l->brother);
		delete l;
		return;
	}
	~loop() {
		res::del(head);
	}

	int getLen();
	int size();
	int getType();
	int getLoopCounts();
	int isVirtual();
	string getFlag();
	string getSS() const;
	string getSeq() const;
	string ss() const;
	string seq() const;

	int type = 0; // -2:pseudoknot
	int len = 0;
	int isPartial = 0;
	res *head = NULL;
	helix s;
	vector<pair<int, int>> hinges;
	deque<loop> sons;
	loop *son = NULL;
	loop *brother = NULL;
};

} // namespace nuc2d

} // namespace jian

#endif

