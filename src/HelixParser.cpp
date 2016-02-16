#include <jian/etl.h>
#include <jian/nuc3d/HelixParser.h>

int main(int argc,char **argv) {
    jian::Par par(argc, argv);
    jian::nuc3d::HelixParser helix_parser;
    auto result = helix_parser(jian::pdb::RNA(par["model"][0]));
    auto print = [](auto &f) {
        for (int i = 0; i < 3; i++) std::cout << f[i] << ' '; std::cout << std::endl;
    };
    print(result.origin); print(result.x); print(result.y); print(result.z);
    std::cout << result.theta << std::endl;
    std::cout << result.phi << std::endl;
    return 0;
}

