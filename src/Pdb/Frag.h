#ifndef FRAG_H
#define FRAG_H

#include "../Utils.h"

namespace jian {

class Frag {
public:
	Frag() {}
	Frag(int type, int len, std::string seq, std::string ss, int *num, int flag = -1) : _type(type), _len(len), _seq(seq), _ss(ss), _flag(flag) {
		_num = new int[len];
		for (int i = 0; i < len; i++) {
			_num[i] = num[i];
		}
	}
	~Frag() {
		delete [] _num;
	}

	int _type = -1;
	int _len = 0;
	std::string _seq;
	std::string _ss;
	int *_num = NULL;
	int _flag = -1;
	Frag *_son = NULL;
	Frag *_brother = NULL;
};

} /// namespace jian

#endif


