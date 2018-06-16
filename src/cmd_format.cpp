#include "nsp.hpp"
#include "pdb.hpp"
#include "rtsp_format.hpp"

namespace jian {

namespace {

void m_format(const Par &par, S mol_type) {
    S in = par[2];

    Format format;
    Molecule mol;

    mol_read(mol, in, mol_type);
    if (par.has("format")) {
        mol = format(mol);
    }
    else {
        format.sort(mol);
    }
    JN_OUT << mol << std::endl;
}

REGISTER_NSP_COMPONENT(format) {
    m_format(par, "");
}

REGISTER_NSP_COMPONENT(rna) {
    m_format(par, "RNA");
}

REGISTER_NSP_COMPONENT(dna) {
    m_format(par, "DNA");
}

REGISTER_NSP_COMPONENT(protein) {
    m_format(par, "protein");
}

} // namespace

}

