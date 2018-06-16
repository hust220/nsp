#pragma once

#include <array>
#include "mat_base.hpp"
#include "fftw3.h"

namespace jian {

#define JN_MAT3_EACH(mat, i, j, k) for (int i = 0; i < (mat).size[0]; i++) for (int j = 0; j < (mat).size[1]; j++) for (int k = 0; k < (mat).size[2]; k++)

template<typename NumType>
class Matrix<3, NumType> : public MatBase<3, NumType> {
    public:
        using value_type = NumType;
        using type = Matrix<3, value_type>;

        std::array<int, 3> size;
        value_type *p;

        Matrix() {
            init();
        }

        Matrix(int a, int b, int c) {
            init();
            resize(a, b, c);
        }

        Matrix(const Matrix<3, value_type> &m) { init(); *this = m; }

        template<typename T> Matrix(const T &m) { init(); *this = m; }

        Matrix(Matrix<3, value_type> &&m) { init(); std::swap(size, m.size); std::swap(p, m.p); m.init(); }

        type &operator=(const Matrix<3, value_type> &m) {
            resize(m.size[0], m.size[1], m.size[2]);
            for (int i = 0; i < size[0]; i++) for (int j = 0; j < size[1]; j++) for (int k = 0; k < size[2]; k++) {
                p[index(i, j, k)] = m(i, j, k);
            }
            return *this;
        }

        template<typename T>
        type &operator=(const T &m) {
            resize(m.size[0], m.size[1], m.size[2]);
            for (int i = 0; i < size[0]; i++) for (int j = 0; j < size[1]; j++) for (int k = 0; k < size[2]; k++) {
                p[index(i, j, k)] = m(i, j, k);
            }
            return *this;
        }

        type &operator=(type &&m) {
            std::swap(size, m.size);
            std::swap(p, m.p);
            m.init();
            return *this;
        }

        void init() {
            size = { 0, 0, 0 };
            p = NULL;
        }

        template<typename T> type &set_all(T && v) {
            JN_MAT3_EACH(*this, i, j, k) p[index(i, j, k)] = v;
            return *this;
        }

        value_type sum() const {
            value_type n = 0;
            for (int i = 0; i < size[0]; i++) {
                for (int j = 0; j < size[1]; j++) {
                    for (int k = 0; k < size[2]; k++) {
                        n += p[index(i,j,k)];
                    }
                }
            }
            return n;
        }

        template<typename T, JN_ENABLE(!JN_IS_SAME(value_type, std::complex<T>))>
        value_type max() const {
            value_type max = p[0];
            for (int i = 0; i < size[0]; i++) for (int j = 0; j < size[1]; j++) for (int k = 0; k < size[2]; k++) {
                if (max < p[index(i, j, k)]) max = p[index(i, j, k)];
            }
            return max;
        }

        type &resize(int a, int b, int c) {
            int i, j;

            if (p != NULL) clear();
            size[0] = a;
            size[1] = b;
            size[2] = c;
            p = new value_type[a*b*c];
            return *this;
        }

        type &clear() {
            if (p != NULL) {
                delete[] p;
                p = NULL;
                size = { 0, 0, 0 };
            }
        }

        static type Zero(int a, int b, int c) {
            type m(a, b, c);
            return m.set_all(0);
        }

        static type Ones(int a, int b, int c) {
            type m(a, b, c);
            return m.set_all(1);
        }

        template<typename T>
        static type Constant(int a, int b, int c, T && v) {
            type m(a, b, c);
            return m.set_all(v);
        }

        ~Matrix() {
            clear();
        }

        Int index(Int a, Int b, Int c) const {
            return a * (size[1] * size[2]) + b * size[2] + c;
        }

        value_type &operator()(int a, int b, int c) {
            return p[index(a,b,c)];
        }

        const value_type &operator()(int a, int b, int c) const {
            return p[index(a, b, c)];
        }

        void print() const {
            std::cout << *this;
        }

        template<typename Mat_>
        const type &fft(Mat_ &&mat) const {
            fftw_plan plan = fftw_plan_dft_3d(size[0], size[1], size[2], reinterpret_cast<fftw_complex*>(p), reinterpret_cast<fftw_complex*>(mat.p), FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_execute(plan); /* repeat as needed */
            fftw_destroy_plan(plan);
            return *this;
        }

        template<typename Mat_>
        const type &ifft(Mat_ &&mat) const {
            fftw_plan plan = fftw_plan_dft_3d(size[0], size[1], size[2], reinterpret_cast<fftw_complex*>(p), reinterpret_cast<fftw_complex*>(mat.p), FFTW_BACKWARD, FFTW_ESTIMATE);
            fftw_execute(plan); /* repeat as needed */
            fftw_destroy_plan(plan);
            return *this;
        }

        template<typename _N>
        type &divide(_N &&n) {
            JN_MAT3_EACH(*this, i, j, k) p[index(i,j,k)] /= n;
            return *this;
        }

        friend std::ostream &operator <<(std::ostream &stream, const type &m) {
            JN_MAT3_EACH(m, i, j, k) stream << m(i, j, k) << ' ';
            return stream;
        }


};

using Mat3 = Matrix<3, Num>;
using Mat3f = Matrix<3, float>;
using Mat3d = Matrix<3, double>;
using Mat3i = Matrix<3, int>;
using Mat3u = Matrix<3, unsigned>;
using Mat3c = Matrix<3, std::complex<Num>>;

}


