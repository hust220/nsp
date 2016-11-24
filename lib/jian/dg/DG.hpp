#pragma once

#include <vector>
#include <Eigen/Dense>
#include "../matrix.hpp"
#include "DihBound.hpp"

namespace jian {

class DGImpl;

class DG {
public:
    using DihEntry = std::vector<int>;
    using DistBoundType = Mat;
    using DihBoundType = DihBound;

    DG();
    ~DG();
    DG(const Mat &);
    DG(const DistBoundType &, const DihBoundType &);
    Mat operator ()();
    Mat operator ()(const Mat &);
    Mat operator ()(const Mat &, const DihBound &);

private:
    DGImpl *_impl;
};

} // namespace jian

