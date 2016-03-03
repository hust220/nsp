#include <jian/pdb.h>
#include <jian/nuc3d/ParseHelix.h>

int main(int argc,char **argv) {
    jian::Par par(argc, argv);
    jian::nuc3d::ParseHelix parse_helix;
    auto result = parse_helix(jian::RNA(par["model"][0]));
    auto print = [](auto &f) {
        for (int i = 0; i < 3; i++) std::cout << f[i] << ' '; std::cout << std::endl;
    };
    print(result.origin); print(result.x); print(result.y); print(result.z);
    std::cout << result.theta << std::endl;
    std::cout << result.phi << std::endl;
    return 0;
}

