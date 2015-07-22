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


