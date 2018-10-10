#include "nsp.hpp"
#include "pdb.hpp"
#include "rtsp_format.hpp"
#include "pdb_reader.hpp"

namespace jian {

namespace {

template<typename Residue_>
void remove_phos_group(Residue_ &residue) {
    Residue_ r = residue;
    r.clear();
    for (auto && atom : residue) {
        if (atom.name != "P" && atom.name != "O1P" && atom.name != "O2P" && atom.name != "OP1" && atom.name != "OP2") {
            r.push_back(atom);
        }
    }
    residue = r;
}

void m_format(const Par &par, S mol_type) {
    S in = par[2];

    Format format;
    Molecule mol;

    if (par.has("format")) {
        mol_read(mol, in, mol_type);
        mol = format(mol);
    }
    else if (par.has("sort")) {
        mol_read(mol, in, mol_type);
        format.sort(mol);
    }
    else if (par.has("rna")) {
        PdbReader reader(mol);
        reader.read(in);
        static std::vector<std::string> chain_names = {
            "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N",
            "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z",
            "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n",
            "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z",
             "0",  "1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9",
            "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
            "20", "21", "22", "23", "24", "25", "26", "27", "28", "29",
            "30", "31", "32", "33", "34", "35", "36", "37", "38", "39",
            "40", "41", "42", "43", "44", "45", "46", "47", "48", "49",
            "50", "51", "52", "53", "54", "55", "56", "57", "58", "59",
            "60", "61", "62", "63", "64", "65", "66", "67", "68", "69",
            "70", "71", "72", "73", "74", "75", "76", "77", "78", "79",
            "80", "81", "82", "83", "84", "85", "86", "87", "88", "89",
            "90", "91", "92", "93", "94", "95", "96", "97", "98", "99",
           "100","101","102","103","104","105","106","107","108","109",
           "110","111","112","113","114","115","116","117","118","119",
           "120","121","122","123","124","125","126","127","128","129",
           "130","131","132","133","134","135","136","137","138","139",
           "140","141","142","143","144","145","146","147","148","149",
           "150","151","152","153","154","155","156","157","158","159",
           "160","161","162","163","164","165","166","167","168","169",
           "170","171","172","173","174","175","176","177","178","179",
           "180","181","182","183","184","185","186","187","188","189",
           "190","191","192","193","194","195","196","197","198","199"
        };

        for (auto && model : mol) {
            int chain_index = 0;
            for (auto && chain : model) {
                chain.name = chain_names[chain_index];
                for (auto && res : chain) {
                    if      (res.name == "A5" || res.name == "A3") res.name = "A";
                    else if (res.name == "U5" || res.name == "U3") res.name = "U";
                    else if (res.name == "G5" || res.name == "G3") res.name = "G";
                    else if (res.name == "C5" || res.name == "C3") res.name = "C";

                    if (res.name == "A" || res.name == "U" || res.name == "G" || res.name == "C") {
                        res = format(res);
                    }
                }
                chain_index++;
            }
        }
    }
    else if (par.has("amber")) {
        mol_read(mol, in, mol_type);
        for (auto && model : mol) {
            for (auto && chain : model) {
                remove_phos_group(chain[0]);
            }
        }
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

