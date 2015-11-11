#ifndef LOOPMODELLING2_H
#define LOOPMODELLING2_H

#include "../Pdb"

class LoopModelling2 {
public:
	LoopModelling2(string, string);
	RNA *run();
	double getDist();
	Residue *buildNuc(double *, double *, double *, double, Point *, Point *, int);

private:
	string seq;
	string ss;
	int num;
	int len;
	int bins;
	Point *c;
	double *a;
	double *t;
	double *d;
	int dinucQuant[16];
	double **dinuctor[16];
};

#endif

