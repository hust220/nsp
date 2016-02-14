#ifndef JIAN_PDB_IFMODEL_H
#define JIAN_PDB_IFMODEL_H

#include "Molecule.h"

namespace jian {
namespace pdb {

struct _IF {};

using IFRNA = Model<_RNA, _IF>;
using IFRNAResidue = Residue<_RNA, _IF>;

inline bool is_residue_paired(const IFRNAResidue &res1, const IFRNAResidue &res2) {
    static std::unordered_map<std::string, unsigned int> cvt {{"A", 1}, {"U", 2}, {"G", 4}, {"C", 8}};

    Atom atoms[4] {res1.atom("C"), res1.atom("O"), res2.atom("C"), res2.atom("O")};
    std::array<double, 4> dists {geom::distance(atoms[0], atoms[2]), geom::distance(atoms[1], atoms[3]),
                                 geom::distance(atoms[0], atoms[3]), geom::distance(atoms[1], atoms[2])};

    return std::set<int>{3, 6, 12}.count(cvt[res1._name] | cvt[res2._name]) && 
           (dists[0] > 12 && dists[0] < 13) && (dists[1] > 5  && dists[1] < 7)  &&
           (dists[2] > 8  && dists[2] < 11) && (dists[3] > 8  && dists[3] < 11);
}


} // namespace pdb
} // namespace jian

#endif

