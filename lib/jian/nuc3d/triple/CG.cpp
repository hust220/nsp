#include "CG.hpp"

namespace jian {
namespace triple {

CG::res_atoms_map_t CG::res_atoms_map {
    {"DA", {"C5*", "O3*", "C1*", "C2", "C6"}},
    { "A", {"C5*", "O3*", "C1*", "C2", "C6"}},
    {"DT", {"C5*", "O3*", "C1*", "C2", "C4"}},
    { "U", {"C5*", "O3*", "C1*", "C2", "C4"}},
    {"DG", {"C5*", "O3*", "C1*", "C2", "C6"}},
    { "G", {"C5*", "O3*", "C1*", "C2", "C6"}},
    {"DC", {"C5*", "O3*", "C1*", "C2", "C4"}},
    { "C", {"C5*", "O3*", "C1*", "C2", "C4"}}
};


} // namespace triple
} // namespace jian

