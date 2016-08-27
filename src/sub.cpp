#include "nsp.hpp"
#include <jian/pdb.hpp>
#include <jian/geom.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(sub) {
    std::deque<int> ls;
    std::vector<std::string> v;
    int beg, end;
    for (auto && s : par[std::vector<std::string>{"num", "n"}]) {
        tokenize(s, v, "-");
        if (v.size() == 1) {
            ls.push_back(JN_INT(s)-1);
        } else if (v.size() == 2) {
            beg = JN_INT(v[0])-1;
            end = JN_INT(v[1])-1;
            for (int i = beg; i <= end; i++) {
                ls.push_back(i);
            }
        }
    }
    std::cout << sub(mol_read_to<Model>(par[std::vector<std::string>{"s", "pdb", "p"}][0]), ls) << std::endl;
}

} // namespace jian
















