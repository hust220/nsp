#ifndef R2D_H
#define R2D_H

#include "../Utils.h"
#include "Frag.h"
#include "Nuc.h"

namespace jian {

class R2D {
public:
	R2D() {}
	R2D(std::string seq, std::string ss, int view = 0) {
		(*this)(seq, ss, view);
	}
	~R2D() {
		delTree(_head);
		delTree(_pseudo_head);
	}

	void operator ()(std::string, std::string, int = 0);
	static Frag *setTree(std::string, std::string);
	void print();
	void delTree(Frag *);
	void printTree(Frag *);

	int _len = 0;
	std::string _seq;
	std::string _ss;
	int _view = 0;
	Frag *_head = NULL;
	Frag *_pseudo_head = NULL;

};

} /// namespace jian

#endif

