#ifndef MATRIX_H
#define MATRIX_H

#include <cstring>
#include "Point.h"
#include "Obj.h"
#include "MLib.h"
#include "../Eigen/Dense"
#include "../Eigen/SVD"
#include "../Eigen/Geometry"
using namespace std;
using namespace Eigen;

namespace jian {

MatrixXf mat_from_file(std::string file);

namespace mat {

template<typename F, typename T1>
auto sum(F &&f, T1 &&t1) ->decltype(f(t1)) {
    return f(std::forward<T1>(t1));
}

template<typename F, typename T1, typename... T2>
auto sum(F &&f, T1 &&t1, T2 && ...t2) ->decltype(f(t1)) {
    return f(std::forward<T1>(t1)) + sum(f, std::forward<T2>(t2)...);
}

template<typename T1>
MatrixXf hstack(T1 &&mat1) {
    return std::forward<T1>(mat1);
}

template<typename T1, typename T2, typename... T3>
MatrixXf hstack(T1 &&mat1, T2 &&mat2, T3 && ...mat3) {
    if (mat1.rows() == 0) {
        return mat2;
    } else if (mat2.rows() == 0) {
        return mat1;
    }
    assert(mat1.cols() == mat2.cols());
    MatrixXf mat(mat1.rows() + mat2.rows(), mat1.cols());
    int i = 0;
    for (; i < mat1.rows(); i++) {
        for (int j = 0; j < mat1.cols(); j++) {
            mat(i, j) = mat1(i, j);
        }
    }
    for (int ii = 0; i < mat.rows(); i++, ii++) {
        for (int j = 0; j < mat2.cols(); j++) {
            mat(i, j) = mat2(ii, j);
        }
    }
    return hstack(mat, std::forward<T3>(mat3)...);
}

} /// namespace mat

class Matr_ {
public:
	Matr_(int, int);
	Matr_(int, int, int);
	Matr_(Matr_ *);
	Matr_(Obj<Matr_>);
	Matr_(string);
	~Matr_();
	Matr_ *copy();

	double *operator [](int);
	
	void read(string);
	void print(int = 8);
	
	void identity();
	void changeRow(int, int);
	Matr_ *colVec(int);
	Matr_ *subMatr_(int, int);

	double det();
	void assign(int, int);
	void assign(Matr_ *);
	Matr_ *add(double);
	Matr_ *add(Matr_ *);
	Matr_ *minus(double);
	Matr_ *minus(Matr_ *);
	Matr_ *multiply(double);
	Matr_ *multiply(Matr_ *);
	Obj<Matr_> dot(Obj<Matr_>);
	Matr_ *sqrt();
	Matr_ *inverse();
	Matr_ *transpose();
	
	void qr(Matr_ *q, Matr_ *r);
	void ed(Matr_ *p1, Matr_ *e, Matr_ *p2);
	void svd(Matr_ *q1, Matr_ *s, Matr_ *q2);

	Obj<Point> point(int);

	double **data;
	int row, col;
};

} /// namespace jian

#endif


