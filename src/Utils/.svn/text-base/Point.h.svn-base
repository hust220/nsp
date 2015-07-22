#ifndef POINT_H
#define POINT_H

#include "Obj.h"

namespace jian {

class Point {
public:
	Point(double = 0, double = 0, double = 0, int = 0);
	Point(const Point &);
	Point(Point *);
	double dist(Point *);
	double dist(const Point &);
	double dist(Obj<Point>);
	Point *rotate(Point *, Point *, double);
	Point *rotate(Point &, Point &, double);
	Point *rotate(Point *, double);
	static void rotate(Point *, int, Point, Point, double);
	static void coincide(Point *, int, Point, Point);
	static Point *normalVector(Point *, Point *, Point *);
	static Point *normalVector(Point &, Point &, Point &);
	static Point *rotate(Point *, Point *, Point *, double, double, double);
	static Point *grow(Point *, Point *, double, double, double);
	static double angle(Point *, Point *, Point *);
	static double angle(Point &, Point &, Point &);
	static double angle(Obj<Point>, Obj<Point>, Obj<Point>);
	static double dihedral(Point &, Point &, Point &, Point &);
	static double dihedral(Point *, Point *, Point *, Point *);
	static double dihedral(Obj<Point>, Obj<Point>, Obj<Point>, Obj<Point>);
	static void minmax(double a, double b, double c, double d, double e, double &min, double &max);
	static double gmMax(double a, double b, double c, double d, double e);
	static double gmMin(double a, double b, double c, double d, double e);
	static double gmMin(double a, double b, double c, double d1, double d2, double e1, double e2);
	static Point *t2c(double *, double *, double *, int);
	static Point *t2c(double *, double *, double *, Point *, Point *, int);
	static double chirality(Point *, Point *, Point *, Point *);
	static double chirality(Point &, Point &, Point &, Point &);
	double &operator[](int);
	const double &operator[](int) const;
	friend Point operator -(const Point &);

	double x;
	double y;
	double z;
	int type;
};

} /// namespace jian

#endif
