#pragma

#include "jian.hpp"

namespace jian {


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

}
