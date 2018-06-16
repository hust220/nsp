#include "nsp.hpp"
#include "matrix.hpp"

namespace jian {

REGISTER_NSP_COMPONENT(dist_weighted_frequencies) {
    S weight_file = par[std::vector<std::string>{"weights", "w"}][0];
    S freq_file = par[std::vector<std::string>{"frequencies", "f"}][0];
    Eigen::MatrixXd w(85, 85);
    std::ifstream ifile(weight_file);
    for (int i = 0; i < 85; i++) {
        for (int j = 0; j < 85; j++) {
            ifile >> w(i, j);
        }
    }
    ifile.close();
    Eigen::MatrixXd *f[67];
    for (int i = 0; i < 67; i++) {
        f[i] = new Eigen::MatrixXd(85, 85);
    }
    ifile.open(freq_file);
    for (int i = 0; i < 85; i++) {
        for (int j = 0; j < 85; j++) {
            for (int k = 0; k < 67; k++) {
                ifile >> (*(f[k]))(i, j);
                (*(f[k]))(i, j) *= w(i, j);
            }
        }
    }
    ifile.close();
    for (int i = 0; i < 85; i++) {
        for (int j = 0; j < 85; j++) {
            for (int k = 0; k < 67; k++) {
                std::cout << (*(f[k]))(i, j) << ' ';
            }
            std::cout << std::endl;
        }
    }
    for (int i = 0; i < 67; i++) {
        delete f[i];
    }
}

}

