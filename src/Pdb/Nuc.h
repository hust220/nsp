#ifndef  NUC_H
#define  NUC_H

#include "../Utils.h"

struct Nuc {
	Nuc() {}
	Nuc(char seq, char ss, int num) : _seq(seq), _ss(ss), _num(num) {}
	char  _seq = 'X';
	char  _ss = '.';
	int  _num = -1;
};

#endif



