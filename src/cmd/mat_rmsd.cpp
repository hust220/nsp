#include <list>
#include <array>
#include <iostream>
#include "nsp.hpp"
#include "pdb.hpp"
#include "geom.hpp"
#include "string.hpp"

namespace jian {

void load_mat(Mat &m, std::string filename) {
    int rows, cols;
    std::ifstream ifile(filename.c_str());
    ifile >> rows >> cols;
    m.resize(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            ifile >> m(i, j);
        }
    }
    ifile.close();
}

REGISTER_NSP_COMPONENT(mat_rmsd) {
    auto g = par.getv("global");
    Mat m1, m2;
    load_mat(m1, g[1]);
    load_mat(m2, g[2]);
    std::cout << geom::rmsd(m1, m2) << std::endl;
}

}

