#include "nsp.hpp"
#include "dg.hpp"

namespace jian {

static Mat read_mat(const std::string &fname) {
    int m, n;
    std::ifstream ifile(fname);
    ifile >> m >> n;
    Mat mat(m, n);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            ifile >> mat(i, j);
        }
    }
    ifile.close();
    return mat;
}

REGISTER_NSP_COMPONENT(dg_test) {
    auto m = read_mat(par.get("m"));

    int n = 10;
    par.set(n, "n");

    std::cout << "Bound:" << std::endl;
    std::cout << m << std::endl;

    auto dg = std::make_unique<Dg>(m, 1);

    for (int i = 0; i < n; i++) {
        auto c = dg->sample_it();

//        std::cout << "Coordinates:" << std::endl;
//        std::cout << c << std::endl;

        std::cout << "Distances:" << std::endl;
        std::cout << dg->c2d(c) << std::endl;
    }
}

}

