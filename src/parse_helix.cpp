#include "nsp.hpp"
#include <jian/nuc3d/ParseHelix.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(parse_helix) {
    auto result = parse_helix(mol_read_to<Model>(par.get("s")));
    auto print = [](auto &f) {
        for (int i = 0; i < 3; i++) std::cout << f[i] << ' '; std::cout << std::endl;
    };
    print(result.origin); print(result.x); print(result.y); print(result.z);
    std::cout << result.theta << std::endl;
    std::cout << result.phi << std::endl;
}

} // namespace jian

