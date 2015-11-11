#ifndef LOOPMODELLING5_H
#define LOOPMODELLING5_H

#include "../Pdb"

class LoopModelling5 {
public:
	LoopModelling5(string, string);
	LoopModelling5(string);
	void init();
	Obj<RNA> run();

private:
	void setHelixPar(Frag *);
	void ss2ct();
	double at(Matr_ *, int, int, int);
	void assign(Matr_ *, int, int, double, int);

	void to3pt();
	void mc();
	Obj<RNA> c2a();

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
	Obj<Matr_> coord; // coordinates when using one atom to represent every residue
	Obj<Matr_> constraints;
	DG *dg;
};

#endif

