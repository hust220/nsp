#ifndef ANGANAL_H
#define ANGANAL_H

#include "Analysis.h"

namespace jian {

class AngAnal: public Analysis {
public:
	AngAnal(double = 1.0);
	~AngAnal();

	void readRNA(Obj<RNA>);
	void train();
	double scoring();
	void readParm(char *);
	void readParm(string);
	void initProb();
	void printParm();
	void printProb();
	double getScore();

private:
	void delPoints();
	void initPoints(int);

	double interval;
	int bins;

	int *obsParm;
	double *obsProb;
	int *refParm;
	double *refProb;
	int len;
	Point **p, **o5_, **c5_, **c4_, **c3_, **o3_;
	double score;
};

} /// namespace jian

#endif

