#ifndef LENANAL_H
#define LENANAL_H

#include "Analysis.h"
#include <pdb/util.h>

namespace jian {
namespace scoring {

class LenAnal: public Analysis {
public:
	LenAnal(double = 0.01);
	~LenAnal();

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

} /// namespace scoring
} /// namespace jian

#endif

