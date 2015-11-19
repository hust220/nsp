#ifndef SUPPOS_H
#define SUPPOS_H

#include "std.h"
#include "Point.h"
#include "Matr_.h"

namespace jian {

class SupPos {
public:
	tuple<Matrix3f, Point, Point> operator ()(const MatrixXf &, const MatrixXf &);
	tuple<Matrix3f, Point, Point> operator ()(MatrixXf &, const MatrixXf &, const MatrixXf &);
	Matrix3f operator ()(Matr_ *, Matr_ *);

	MatrixXf x;
	MatrixXf y;
	int len = 0;
	Matrix3f rot;
	Point c1;
	Point c2;
	double rmsd = 0;
};

} /// namespace jian

#endif

