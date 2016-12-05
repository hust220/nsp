#include <random>

#include "rand.hpp"

BEGIN_JN

	thread_local static std::mt19937 _rand_engine{ 11 };
	thread_local static std::uniform_real_distribution<double> _unif_distr{ 0, 1 };

	double rand() {
		return _unif_distr(_rand_engine);
	}

	void seed(int t) {
		_rand_engine.seed(t);
	}

END_JN

