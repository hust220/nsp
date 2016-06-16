#include <iostream>
#include <map>
#include <cmath>
#include <Eigen/Dense>
#include "ParseHelix.hpp"
#include "../geom.hpp"

namespace jian {

class ParseHelix {
public:
    using val_t = double;
    using Vec = parse_helix_t::Vec;
    using Mat = Eigen::Matrix<val_t, -1, -1>;
    using Par = struct {val_t theta, phi, h, d, R;};
    using Result = parse_helix_t;

    Par _par {0.562, 1.57, 2.84, -4, 9.7};
    std::map<int, Mat> _cache;

    Mat make_standard_helix(int n) {
        if (_cache.count(n)) return _cache[n];
        else {
            Mat helix = Mat::Zero(n * 2 + 2, 3);
            for (int i = 0; i < n + 1; i++) {
                helix(i, 0) = _par.R*std::cos(i*_par.theta);
                helix(i, 1) = _par.R*std::sin(i*_par.theta);
                helix(i, 2) = i*_par.h;
                helix(2*n+1-i, 0) = _par.R*std::cos(i*_par.theta+_par.phi);
                helix(2*n+1-i, 1) = _par.R*std::sin(i*_par.theta+_par.phi);
                helix(2*n+1-i, 2) = i*_par.h+_par.d;
            }
            return helix;
        }
    }

    Mat get_mat_helix(const Model &helix) {
        int len = num_residues(helix); 
        Mat mat(len, 3); 
        int index = 0; 
        for (auto &chain : helix) for (auto &res : chain) {
            auto atm = atom(res, "C4*");
            for (int i = 0; i < 3; i++) mat(index, i) = atm[i];
            index++;
        }
//        std::cout << mat << std::endl;
        return mat;
    }

    Result parse(const Model &helix) {
        auto mat = get_mat_helix(helix); 
        int len = mat.rows();
        auto sp = geom::suppos(make_standard_helix(len/2-1), mat);
        Vec origin, x, y, z; origin << 0, 0, 0; x << 1, 0, 0; y << 0, 1, 0; z << 0, 0, 1;
        geom::apply_suppos(origin, sp);
        geom::apply_suppos(x, sp); geom::apply_suppos(y, sp); geom::apply_suppos(z, sp);
        val_t theta = geom::angle(std::vector<val_t>{0, 0, 1}, Vec::Zero(), z);
        val_t phi = geom::angle(std::vector<val_t>{1, 0, 0}, Vec::Zero(), std::vector<val_t>{z[0], z[1], 0});
        return Result{origin, x-origin, y-origin, z-origin, theta, phi};
    }

};

static ParseHelix parser;

void dihs_std_helix() {
    auto m = parser.make_standard_helix(20);
    for (int i = 0; i < 15; i++) {
        std::cout << geom::dihedral(m.row(i), m.row(i+1), m.row(i+2), m.row(i+3)) << std::endl;
    }
}

parse_helix_t parse_helix(const Model &helix) {
    return parser.parse(helix);
}

Eigen::MatrixXd make_standard_helix(int n) {
    return parser.make_standard_helix(n);
}
} // namespace jian

