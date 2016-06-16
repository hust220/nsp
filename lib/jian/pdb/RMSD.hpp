#pragma once

#include <Eigen/Dense>
#include "Model.hpp"

namespace jian {

class RMSD {
public:
    double operator ()(const Model &model1, const Model &model2);
    void setXY();
    double run(const Model &model1, const Model &model2);

private:
    double rmsd;
    Model _model1;
    Model _model2;
    int _len;
    Eigen::MatrixXf x;
    Eigen::MatrixXf y;
    Eigen::Matrix3f r;
};

} // namespace jian

