#ifndef LOOPMODELLING_H
#define LOOPMODELLING_H

#include "../Pdb"

class LoopModelling {
public:
	LoopModelling(string, string);
	void init();
	RNA *run();

	void mc(int);
	void ss2ct();
	double totEn();
	double energy(int, int, int);
	double at(Matr_ *, int, int, int);
	void assign(Matr_ *, int, int, double, int);
	RNA *newRNA();
	Residue *buildNuc(Matr_ *, int, Point *, int);
	Point **cmpltBase();

	string seq;
	string ss;
	string lib;
	int len;
	int *lens;
	int resLen;
	int *type;
	vector<int> brk;
	Matr_ *ct;
	Matr_ *bound;
	Matr_ *coord;
	Matr_ **nucPrmt;
	Matr_ **nucPrmt3;
	Matr_ *nucPrmt4;
	Matr_ **bpPrmt;
	Matr_ **adjPrmt;
	Matr_ **nxtPrmt;
	Matr_ **stackPrmt;
	DG *dg;
};

#endif

