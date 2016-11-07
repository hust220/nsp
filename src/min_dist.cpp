#include "nsp.hpp"
#include <jian/utils/file.hpp>
#include <jian/pdb.hpp>
#include <jian/geom.hpp>

namespace jian {

static double min_distance(const Residue &r1, const Residue &r2) {
    double d, min = 99999;
    for (auto && a1 : r1) {
        for (auto && a2 : r2) {
            d = geom::distance(a1, a2);
            if (d < min) {
                min = d;
            }
        }
    }
    return min;
}

REGISTER_NSP_COMPONENT(min_dist) {
    Chain chain = read_model_to_chain(par[std::vector<std::string>{"s", "pdb", "p"}][0]);
    int len = chain.size();
    std::deque<int> ls; 
    if (par.has("num")) {
        for (auto && n : par["num"]) {
            ls.push_back(JN_INT(n) - 1);
        }
        std::cout << min_distance(chain[ls[0]], chain[ls[1]]) << std::endl;
    } else if (par.has("list")) {
		BEGIN_READ_FILE(par["list"][0], " ") {
			int a = std::stoi(F[0]) - 1;
			int b = std::stoi(F[1]) - 1;
			std::cout << a + 1 << ' ' << b + 1 << ' ' << min_distance(chain[a], chain[b]) << std::endl;
		} END_READ_FILE;
    }
}

} // namespace jian
















