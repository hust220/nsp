#ifndef DIHANAL_H
#define DIHANAL_H

#include "Analysis.h"
#include <pdb/util.h>

namespace jian {
namespace scoring {

class DihAnal: public Analysis {
public:
	DihAnal(double = 1.0);
	~DihAnal();

	void readRNA(const RNA &);
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
	Point **p, **o5_, **c5_, **c4_, **c3_, **o3_, **c2_, **c1_, **b1, **b2;
	double score;
};

} /// namespace scoring
} /// namespace jian

#endif

