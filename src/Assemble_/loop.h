#ifndef LOOP_H
#define LOOP_H

#include "helix.h"
#include "strand.h"

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

	static RNA *create(string, string);
	int getLen();
	int size();
	int getType();
	int getLoopCounts();
	int isVirtual();
	string getFlag();
	string getSS();
	string getSeq();
	void getStrand(vector<strand> &);

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

class loopInfo {
public:
	loopInfo();
	loopInfo(string, string, string, loop *);

	string name;
	int n;
	int len;
	string flag;
	string seq;
	string ss;
	string src;
	string family;
	double score;
	loopInfo *next;
};

/*
class loopInfoEx {
public:
	loopInfoEx();
	loopInfoEx(loopInfo *li, double score);
	
	loopInfo *li;
	loopInfoEx *next;
	double score;
};
*/

class loopInfoList {
public:
	loopInfoList();
	~loopInfoList();

	int getLen();
	void add(loopInfo *li);

	loopInfo *head;

private:
	int len;
};

#endif

