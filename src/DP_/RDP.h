#ifndef RDP_H
#define RDP_H

#include "DP.h"

class RDP: public DP {
public:
	RDP();
	RDP(string seq): DP(seq) {
	}
protected:
	int *types_;
	double *score_matrix_;
};





#endif //RDP_H


