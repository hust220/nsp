#include <iostream>
#include <map>
#include <cmath>
#include <Eigen/Dense>
#include <mutex>
#include "rtsp_parse_helix.hpp"
#include "geom.hpp"
#include "log.hpp"

namespace jian {

namespace parse_helix_detail {

std::mutex mt;

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
        std::lock_guard<std::mutex> gd(mt);
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
//        LOG << mat << std::endl;
        return mat;
    }

    Result parse(const Model &m) {
        return parse(get_mat_helix(m));
    }

    Result parse(const Mat &mat) {
        Vec origin, x, y, z;
		int len;

        len = mat.rows();
		geom::Superposition<double> sp(make_standard_helix(len / 2 - 1), mat);
		origin << 0, 0, 0; x << 1, 0, 0; y << 0, 1, 0; z << 0, 0, 1;
		sp.apply(origin);
		sp.apply(x);
		sp.apply(y);
		sp.apply(z);
        val_t theta = geom::angle(std::vector<val_t>{0, 0, 1}, Vec::Zero(), std::vector<val_t>{z[0]-origin[0],z[1]-origin[1],z[2]-origin[2]});
        val_t phi = geom::angle(std::vector<val_t>{1, 0, 0}, Vec::Zero(), std::vector<val_t>{z[0]-origin[0], z[1]-origin[1], 0});
		if (z[1] < 0) phi = 2*PI-phi;
        return Result{origin, x-origin, y-origin, z-origin, theta, phi};
    }

};

ParseHelix parser;

} // namespace parse_helix_detail

void dihs_std_helix() {
    auto m = parse_helix_detail::parser.make_standard_helix(20);
    for (int i = 0; i < 15; i++) {
        LOG << geom::dihedral(m.row(i), m.row(i+1), m.row(i+2), m.row(i+3)) << std::endl;
    }
}

parse_helix_t parse_helix(const Mat &mat) {
    return parse_helix_detail::parser.parse(mat);
}

parse_helix_t parse_helix(const Model &helix) {
    return parse_helix_detail::parser.parse(helix);
}

Eigen::MatrixXd make_standard_helix(int n) {
    return parse_helix_detail::parser.make_standard_helix(n);
}

}

