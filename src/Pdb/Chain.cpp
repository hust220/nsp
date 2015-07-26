#include "Chain.h"

using namespace jian;

Chain::Chain(vector<string> &lines, string rnaName, string type) {
    /// Set name
    name = "";
    name += lines[0][21];

    /// set rnaName
    this->rnaName = rnaName;

    /// set residues
    vector<string> strings;
    for (int i = 0; i < (int) lines.size(); i++) {
        if (!strings.empty() && lines[i].substr(22, 5) != strings.back().substr(22, 5)) {
            Residue res(strings, rnaName);
            if (boost::to_upper_copy(type) == "RNA") {
                if (set<string>{"A", "U", "G", "C"}.count(res.name)) {
                    residues.push_back(Residue(strings, rnaName));
                }
            } else if (boost::to_upper_copy(type) == "DNA") {
                if (set<string>{"DA", "DT", "DG", "DC"}.count(res.name)) {
                    residues.push_back(Residue(strings, rnaName));
                }
            }
            strings.clear();
        }
        strings.push_back(lines[i]);
    }
    Residue res(strings, rnaName);
    if (boost::to_upper_copy(type) == "RNA") {
        if (set<string>{"A", "U", "G", "C"}.count(res.name)) {
            residues.push_back(Residue(strings, rnaName));
        }
    } else if (boost::to_upper_copy(type) == "DNA") {
        if (set<string>{"DA", "DT", "DG", "DC"}.count(res.name)) {
            residues.push_back(Residue(strings, rnaName));
        }
    }
    strings.clear();
}

int Chain::atom_nums() {
    int num = 0;
    for (auto &&residue: residues) {
        for (auto &&atom: residue) {
            num++;
        }
    }
    return num;
}

void Chain::push(Residue *residue) {
    residues.push_back(*residue);
}

void Chain::push(Residue &residue) {
    residues.push_back(residue);
}

Residue &Chain::operator [](int n) {
    return residues[n];
}

const Residue &Chain::operator [](int n) const {
    return residues[n];
}

vector<Residue>::iterator Chain::begin() {
    return residues.begin();
}

vector<Residue>::iterator Chain::end() {
    return residues.end();
}



