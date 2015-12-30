#include "TestTriplex.h"

int main(int argc, char **argv) {
    jian::test::TestTriplex test_triplex;
    jian::Par par(argc, argv);
    std::string type = par["global"][0];
    auto vec = jian::map<std::vector<int>>([](std::string s){return std::stoi(s);}, par["num"]);
    auto model = jian::Model(par["model"][0]);
    if (type == "chir") {
        test_triplex.chir(jian::Model(par["model"][0]), vec);
    } else if (type == "dist") {
        test_triplex.dist(jian::Model(par["model"][0]), vec);
    }
    return 0;
}
