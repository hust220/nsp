#pragma once

#include <memory>
#include <set>
#include "../utils/Debug.hpp"
#include "../pdb.hpp"
#include "../dg/DG.hpp"
#include "../geom.hpp"
#include "../cg.hpp"

BEGIN_JN

class BuildStrand {
    using Mat = Eigen::MatrixXd;

public:
    std::set<std::string> _coarse_atoms {"C4*"};
    DG dg;

    Chain operator ()(int n, const Mat &a, const Mat &b) {
        return build_strand(n, a, b);
    }

    Chain build_strand(int n, const Mat &a, const Mat &b) {
        if ((a.rows() != _coarse_atoms.size() * 2 && a.rows() != 0) || 
            (b.rows() != _coarse_atoms.size() * 2 && b.rows() != 0)) {
            throw "jian::BuildStrand::build_strand(int, const Mat &, const Mat &) error!";
        }
        Debug::print("build strand:\n");
        auto bound = make_bound(n, a, b);
        auto scaffold = dg(bound);
        superpose_scaffold(scaffold, a, b);
		std::shared_ptr<CG> cg(CG::fac_t::create("1p"));
        auto residues = cg->to_aa(scaffold, 0, int(scaffold.rows())-1);
        auto new_residues = slice(residues, a.rows(), a.rows() + n);
        return new_residues;
    }

    Mat make_bound(int n, const Mat &a, const Mat &b) {
        int size = n * _coarse_atoms.size() + a.rows() + b.rows();
        Mat bound(size, size);
        bound = Mat::Zero(size, size);
        FOR((i, size), FOR((j, i, size), IF(i == j, bound(i, j) = 0, bound(i, j) = 999; bound(j, i) = 2)));
        FOR((i, 1, size), bound(i-1, i) = bound(i, i-1) = 6.1);
        auto mat = hstack(a, b);
        std::vector<int> vec;
        if (a.rows() != 0 && b.rows() != 0) vec = std::vector<int> {0, 1, size - 2, size - 1};
        else if (a.rows() == 0) vec = std::vector<int> {size - 2, size - 1};
        else if (b.rows() == 0) vec = std::vector<int> {0, 1};
        else return bound;
        for (int i = 0; i < mat.rows(); i++) for (int j = i + 1; j < mat.rows(); j++) {
            bound(vec[i], vec[j]) = bound(vec[j], vec[i]) = geom::distance(mat.row(i), mat.row(j));
        }
        return bound;
    }

    template<typename MatType>
    void superpose_scaffold(MatType &scaffold, const Mat &a, const Mat &b) {
        if (a.rows() != 0 || b.rows() != 0) {
            int temp_len = a.rows() + b.rows();
            Mat x(temp_len, 3), y(temp_len, 3);
            if (a.rows() != 0) {
                for (int i = 0; i < 3; i++) {
                    x(0, i) = scaffold(0, i);
                    x(1, i) = scaffold(1, i);
                    y(0, i) = a(0, i);
                    y(1, i) = a(1, i);
                }
            }
            if (b.rows() != 0) {
                for (int i = 0; i < 3; i++) {
                    x(x.rows() - 2, i) = scaffold(scaffold.rows() - 2, i);
                    x(x.rows() - 1, i) = scaffold(scaffold.rows() - 1, i);
                    y(y.rows() - 2, i) = b(0, i);
                    y(y.rows() - 1, i) = b(1, i);
                }
            }
            geom::suppos(scaffold, x, y);
        }
    }

};

END_JN

