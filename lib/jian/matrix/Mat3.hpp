#pragma once

#include <array>
#include "MatBase.hpp"

namespace jian {

	template<typename T>
	class Matrix<3, T> : public MatBase<3, T> {
	public:
		std::array<int, 3> size;
		num_t ***p;

		Matrix() {
			init();
		}

		Matrix(int a, int b, int c) {
			init();
			resize(a, b, c);
		}

		Matrix(Matrix<3, num_t> &&m) {
			init();
			std::swap(size, m.size);
			std::swap(p, m.p);
		}

		Matrix<3, num_t> &operator=(Matrix<3, num_t> &&m) {
			std::swap(size, m.size);
			std::swap(p, m.p);
			m.init();
			return *this;
		}

		void init() {
			size = { 0, 0, 0 };
			p = NULL;
		}

		template<typename F>
		void each(F && f) {
			for (int i = 0; i < size[0]; i++) {
				for (int j = 0; j < size[1]; j++) {
					for (int k = 0; k < size[2]; k++) {
						f(p[i][j][k], i, j, k);
					}
				}
			}
		}

		num_t sum() const {
			num_t n = 0;
			for (int i = 0; i < size[0]; i++) {
				for (int j = 0; j < size[1]; j++) {
					for (int k = 0; k < size[2]; k++) {
						n += p[i][j][k];
					}
				}
			}
			return n;
		}

		void resize(int a, int b, int c) {
			int i, j;

			if (p != NULL) clear();
			size[0] = a;
			size[1] = b;
			size[2] = c;
			p = new num_t**[a];
			for (i = 0; i < a; i++) {
				p[i] = new num_t*[b];
				for (j = 0; j < b; j++) {
					p[i][j] = new num_t[c];
				}
			}
		}

		void clear() {
			if (p != NULL) {
				int i, j;
				for (i = 0; i < size[0]; i++) {
					for (j = 0; j < size[1]; j++) {
						delete[] p[i][j];
					}
					delete[] p[i];
				}
				delete[] p;
				p = NULL;
				size = { 0, 0, 0 };
			}
		}

		static Matrix<3, num_t> Zero(int a, int b, int c) {
			Matrix<3, num_t> m(a, b, c);
			m.each([](num_t &n, int a, int b, int c) {
				n = 0;
			});
			return m;
		}

		static Matrix<3, num_t> Ones(int a, int b, int c) {
			Matrix<3, num_t> m(a, b, c);
			m.each([](num_t &n, int a, int b, int c) {
				n = 1;
			});
			return m;
		}

		static Matrix<3, num_t> Constant(int a, int b, int c, T v) {
			Matrix<3, num_t> m(a, b, c);
			m.each([&v](num_t &n, int a, int b, int c) {
				n = v;
			});
			return m;
		}

		~Matrix() {
			clear();
		}

		num_t &operator()(int a, int b, int c) {
			return p[a][b][c];
		}

		const num_t &operator()(int a, int b, int c) const {
			return p[a][b][c];
		}

		void print() const {
			std::cout << *this;
		}

		friend std::ostream &operator <<(std::ostream &stream, const Matrix<3, num_t> &m) {
			for (int i = 0; i < m.size[0]; i++) {
				for (int j = 0; j < m.size[1]; j++) {
					for (int k = 0; k < m.size[2]; k++) {
						stream << m.p[i][j][k] << ' ';
					}
					stream << std::endl;
				}
			}
			return stream;
		}

	};

	using Mat3 = Matrix<3, float>;
	using Mat3f = Matrix<3, float>;
	using Mat3d = Matrix<3, double>;
	using Mat3i = Matrix<3, int>;
	using Mat3u = Matrix<3, unsigned>;

} // namespace jian


