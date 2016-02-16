#include "Pred2D.h"
#include "../etl.h"

int main(int argc, char **argv) {
    jian::Par par(argc, argv);
    if (par.count("lib")) {
        jian::nuc2d::Pred2D pred(par["lib"][0]);
        std::cout << "min score: " << pred(par["seq"][0]) << std::endl;
    } else {
        jian::nuc2d::Pred2D pred;
        std::cout << "min score: " << pred(par["seq"][0]) << std::endl;
    }
    return 0;
}
