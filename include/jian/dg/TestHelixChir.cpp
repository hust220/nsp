#include "TestHelixChir.h"

int main(int argc, char **argv) {
    jian::Par par(argc, argv);
    jian::dg::TestHelixChir test_helix_chir;
    test_helix_chir(jian::Model(par["model"][0]));
    return 0;
}








