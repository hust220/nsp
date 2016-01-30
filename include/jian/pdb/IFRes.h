#ifndef JIAN_PDB_IFRES_H
#define JIAN_PDB_IFRES_H

#include "Residue.h"
#include "../geom.h"

namespace jian {
namespace pdb {

class IFRes : public Residue {
public:
    IFRes() = default;

    IFRes(MolFile &mol_file) : Residue(mol_file) {}

    bool is_paired(const IFRes &if_res) const {
        static std::unordered_map<std::string, unsigned int> cvt {{"A", 1}, {"U", 2}, {"G", 4}, {"C", 8}};

        Atom atoms[4] {(*this)["C"], (*this)["O"], if_res["C"], if_res["O"]};
        std::array<double, 4> dists {geom::distance(atoms[0], atoms[2]), geom::distance(atoms[1], atoms[3]),
                                     geom::distance(atoms[0], atoms[3]), geom::distance(atoms[1], atoms[2])};

        std::cout << num - 1 << '-' << if_res.num - 1 << ' ';
        for (auto i : dists) std::cout << i << ' ';
        std::cout << std::endl;

        return std::set<int>{3, 6, 12}.count(cvt[name] bitor cvt[if_res.name]) and
               (dists[0] > 12 and dists[0] < 13) and (dists[1] > 5  and dists[1] < 7)  and
               (dists[2] > 8  and dists[2] < 11) and (dists[3] > 8  and dists[3] < 11);
    }
};


} // namespace pdb
} // namespace jian

#endif

