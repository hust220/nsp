#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <util/util.h>

namespace jian {
namespace scoring {

class Analysis {
public:
	virtual void train() = 0;
	virtual double scoring() = 0;
	virtual void readParm(char *) = 0;
	virtual void readParm(string) = 0;
	virtual void initProb() = 0;
	virtual void printParm() = 0;
	virtual void printProb() = 0;
	virtual double getScore() = 0;
};

} /// namespace scoring
} /// namespace jian

#endif

