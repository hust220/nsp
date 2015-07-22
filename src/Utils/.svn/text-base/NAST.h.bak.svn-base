#ifndef NAST_H
#define NAST_H

#include "std.h"
#include "Matr_.h"
#include "Point.h"
#include "MRand.h"

class NAST {
public:
	NAST(Point *, int, string, Matr_ * = NULL);
	Point *run();
	void loadPrmt();
	~NAST();

private:
	double energy(Point *);

	Point *coord;
	int len;
	int cycles;
	int *bp;
	string ss;
	Matr_ *tertiary_contacts;
	Matr_ *bounds;
	Matr_ *basepairs;
	vector<double> dList;
	double kb_all, rb_all, kb_hel, rb_hel, ka_nhel, ra_nhel, ka_hel1, ra_hel1, ka_hel2, ra_hel2;
	double k1d_nhel, k2d_nhel, k3d_nhel, rd_nhel, k1d_hel1, k2d_hel1, k3d_hel1, rd_hel1, k1d_hel2, k2d_hel2, k3d_hel2, rd_hel2, k1d_hel3, k2d_hel3, k3d_hel3, rd_hel3;
	double epsilon, sigma1, sigma2, kpt;
};






#endif











