#ifndef RMSD_H
#define RMSD_H

#include "Model.h"

namespace jian {

namespace pdb {

class RMSD {
public:
	RMSD() {}

	double operator ()(const Model &model1, const Model &model2) {
		return run(model1, model2);
	}
	double run(const Model &, const Model &);
	void setXY();

	double rmsd;

	Model _model1;
	Model _model2;
	int _len;
	MatrixXf x;
	MatrixXf y;
	Matrix3f r;
};

} /// namespace pdb

} /// namespace jian

#endif






