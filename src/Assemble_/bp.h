#ifndef BP_H
#define BP_H

#include "res.h"

class bp {
public:
	bp() {}
	bp(bp *b) : res1(b->res1), res2(b->res2), next(b->next) {}
	bp(const bp &b) : res1(b.res1), res2(b.res2), next(b.next) {}
	bp &operator =(const bp &b) {
		res1 = b.res1;
		res2 = b.res2;
		next = b.next;
		return *this;
	}
	static bp *copy(bp *b) {
		if (b == NULL) {
			return NULL;
		}

		bp *bp_ = new bp(b);
		bp_->next = copy(b->next);
		return bp_;
	}
	static void del(bp *b) {
		if (b == NULL) {
			return;
		}

		del(b->next);
		delete b;
		return;
	}

	res res1, res2;
	bp *next = NULL;
};

#endif

