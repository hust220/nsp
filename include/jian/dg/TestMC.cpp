#include "TestMC.h"
#include "../etl.h"

int main(int argc, char **argv) {
    jian::dg::TestMC test_mc;
    jian::Par par(argc, argv);
    test_mc(boost::lexical_cast<int>(par["num"][0]));
    return 0;
}
