#ifndef MRAND_H
#define MRAND_H

#include <stdio.h>
#include <time.h>

#define MATRIX_A              0x9908b0df   
#define UPPER_MASK            0x80000000 
#define LOWER_MASK            0x7fffffff 
#define TEMPERING_MASK_B      0x9d2c5680
#define TEMPERING_MASK_C      0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

namespace jian {

class MRand {
public:
	MRand();
	double run();

private:
	int sstmm();
	double ran_uniform();
	void genrand(double *);

	int N;
	int M;
	unsigned long *mt;
	int mti;
};

} /// namespace jian

#endif

