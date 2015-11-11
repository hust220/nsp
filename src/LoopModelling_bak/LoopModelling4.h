#ifndef LOOPMODELLING4_H
#define LOOPMODELLING4_H

#include "../Pdb"

class LoopModelling4 {
public:
	LoopModelling4(string, string);
	void init();
	RNA *run();
	~LoopModelling4();

private:
	void ss2ct();
	void setChir();
	double at(Matr_ *, int, int, int);
	void assign(Matr_ *, int, int, double, int);
	Residue *buildNuc(Point *, double, int);

	string seq;
	string ss;
	int len;
	int *type;
	int resLen;
	Matr_ *ct;
	Matr_ *chir;
	vector<int> brk;
	Matr_ *bound;
	Matr_ *backbone;
	Matr_ **nucPrmt;
	Matr_ **adjPrmt;
	Matr_ **stackPrmt;
	Matr_ **bpPrmt;
	DG *dg;
};

#endif

