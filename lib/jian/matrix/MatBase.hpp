#pragma

#include "../utils/traits.hpp"

BEGIN_JN


template<int N, typename T>
class MatBase {
public:
	using Num = T;

	enum {
		DIM = N
	};
};

template<int N, typename T>
class Matrix : public MatBase<N, T> {};

END_JN