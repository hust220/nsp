#pragma

namespace jian {


	template<int N, typename T>
	class MatBase {
	public:
		using num_t = T;

		enum {
			DIM = N
		};
	};

	template<int N, typename T>
	class Matrix : public MatBase<N, T> {};

}