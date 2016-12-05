#pragma once

#include <array>
#include "MatBase.hpp"

BEGIN_JN

	template<typename NumType>
	class Matrix<3, NumType> : public MatBase<3, NumType> {
	public:
		std::array<int, 3> size;
		NumType ***p;

		Matrix() {
			init();
		}

		Matrix(int a, int b, int c) {
			init();
			resize(a, b, c);
		}

		Matrix(Matrix<3, NumType> &&m) {
			init();
			std::swap(size, m.size);
			std::swap(p, m.p);
		}

		Matrix<3, NumType> &operator=(Matrix<3, NumType> &&m) {
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

		NumType sum() const {
			NumType n = 0;
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
			p = new NumType**[a];
			for (i = 0; i < a; i++) {
				p[i] = new NumType*[b];
				for (j = 0; j < b; j++) {
					p[i][j] = new NumType[c];
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

		static Matrix<3, NumType> Zero(int a, int b, int c) {
			Matrix<3, NumType> m(a, b, c);
			m.each([](NumType &n, int a, int b, int c) {
				n = 0;
			});
			return m;
		}

		static Matrix<3, NumType> Ones(int a, int b, int c) {
			Matrix<3, NumType> m(a, b, c);
			m.each([](NumType &n, int a, int b, int c) {
				n = 1;
			});
			return m;
		}

		static Matrix<3, NumType> Constant(int a, int b, int c, NumType v) {
			Matrix<3, NumType> m(a, b, c);
			m.each([&v](NumType &n, int a, int b, int c) {
				n = v;
			});
			return m;
		}

		~Matrix() {
			clear();
		}

		NumType &operator()(int a, int b, int c) {
			return p[a][b][c];
		}

		const NumType &operator()(int a, int b, int c) const {
			return p[a][b][c];
		}

		void print() const {
			std::cout << *this;
		}

		friend std::ostream &operator <<(std::ostream &stream, const Matrix<3, NumType> &m) {
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

END_JN


