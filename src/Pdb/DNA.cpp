#include "DNA.h"

namespace jian {

DNA::DNA(char *pdbfile) {
    string str(pdbfile);
    readPDB(pdbfile);
}

DNA::DNA(string pdbfile) {
    readPDB(pdbfile);
}

DNA *DNA::copy() {
    DNA *dna = new DNA;
    dna->name = name;
    dna->len = len;
    dna->chains = chains;
    return dna;
}

void DNA::readPDB(string pdbfile) {
    /* set name */
    if (pdbfile.size() <= 4 || pdbfile.substr(pdbfile.size() - 4, 4).compare(".pdb")) {
        cerr << "Please give me a file ended with '.pdb'\n";
        exit(1);
    }
    name = pdbfile.substr(0, pdbfile.size() - 4);

    /* set chains */
    ifstream ifile(pdbfile.c_str());
    if (!ifile) {
        cerr << "DNA::readPDB error! Open file \"" << pdbfile << "\" failed!" << endl;
        exit(1);
    }
    string line;
    vector<string> lines;
    int n = 0;
    while (getline(ifile, line, '\n')) {
        if (!line.compare(0, 4, "ATOM")) {
            n++;
            if (lines.size() != 0 && line[21] != lines.back()[21]) {
                Chain chain(lines, name, "DNA");
                lines.clear();
                if (!chain.residues.empty()) {
                    chains.push_back(chain);
                }
            }
            lines.push_back(line);
        }
    }

    Chain chain(lines, name, "DNA");
    lines.clear();
    if (!chain.residues.empty()) {
        chains.push_back(chain);
    }

    if (chains.empty()) {
        cerr << "The file '" << pdbfile << "' has nothing!" << endl;
        exit(1);
    }

    ifile.close();
}

Chain &DNA::operator [](int n) {
    return chains[n];
}

void DNA::push(Chain *chain) {
    chains.push_back(*chain);
}

void DNA::push(Chain &chain) {
    chains.push_back(chain);
}

string DNA::getSeq() {
    string str;
    for (int i = 0; i < (int) chains.size(); i++) {
        for (int j = 0; j < (int) chains[i].residues.size(); j++) {
            string temp = chains[i].residues[j].name;
            for (int k = 0; k < (int) temp.size(); k++) {
                if (temp[k] == 'A' || temp[k] == 'U' || temp[k] == 'C' || temp[k] == 'G') {
                    str += temp[k];
                    break;
                }
            }
        }
    }
    return str;
}

string DNA::getChain() {
    string str;
    for (int i = 0; i < (int) chains.size(); i++) {
        str += chains[i].name;
    }
    return str;
}

inline int DNA::getLen() {
    int i, j, flag;

    for (flag = 0, i = 0; i < (int) chains.size(); i++) {
        for (j = 0; j < (int) chains[i].residues.size(); j++) {
            flag++;
        }
    }

    return flag;
    
}

void DNA::updateChains(string ss) {
    // compute counts of chain
    int temp = 0;
    for (int i = 0; i < (int) ss.size(); i++) {
        if (ss[i] == '&' || (ss[i] == ')' && i - 1 >= 0 && ss[i - 1] == '(')) {
            temp++;
        }
    }
    temp++;

    // construct new chains
    Chain *chain = new Chain[temp];
    int i1 = 0;
    for (int j = 0, k = 0; i1 < (int) chains.size(); i1++) {
        for (int i2 = 0; i2 < (int) chains[i1].residues.size();) {
            if (ss[k] == '&' || (ss[k] == ')' && k - 1 >= 0 && ss[k - 1] == '(')) {
                k++;
                j++;
                continue;
            }
            chain[j].residues.push_back(chains[i1].residues[i2]);
            k++;
            i2++;
        }
    }
    chains.clear();
    for (int i1 = 0; i1 < temp; i1++) {
        chains.push_back(chain[i1]);
    }
    delete [] chain;
}

int DNA::totalAtoms() {
    int i, j, k, flag;

    for (flag = 0, i = 0; i < (int) chains.size(); i++) {
        for (j = 0; j < (int) chains[i].residues.size(); j++) {
            for (k = 0; k < (int) chains[i].residues[j].atoms.size(); k++) {
                flag++;
            }
        }
    }

    return flag;
}

double DNA::getDist(int a, int b) {
    int i, j, k, flag;
    double x1, y1, z1, x2, y2, z2;

    if (a > b || b > getLen() || a > getLen()) {
        return -1;
    }

    for (flag = 0, i = 0; i < (int) chains.size(); i++) {
        for (j = 0; j < (int) chains[i].residues.size(); j++) {
            flag++;
            if (flag == a) {
                for (k = 0; k < (int) chains[i].residues[j].atoms.size(); k++) {
                    if (chains[i].residues[j].atoms[k].name == "O5*") {
                        x1 = chains[i].residues[j].atoms[k].x;
                        y1 = chains[i].residues[j].atoms[k].y;
                        z1 = chains[i].residues[j].atoms[k].z;
                    }
                }
            }
            if (flag == b) {
                for (k = 0; k < (int) chains[i].residues[j].atoms.size(); k++) {
                    if (chains[i].residues[j].atoms[k].name == "O3*") {
                        x2 = chains[i].residues[j].atoms[k].x;
                        y2 = chains[i].residues[j].atoms[k].y;
                        z2 = chains[i].residues[j].atoms[k].z;
                    }
                }
            }
        }
    }

    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
}

void DNA::setResNum() {
    int i, j;
    int num;

    for (num = 0, i = 0; i < (int) chains.size(); i++) {
        for (j = 0; j < (int) chains[i].residues.size(); j++) {
            num++;
            chains[i].residues[j].num = num;
        }
    }
}

ostream &operator <<(ostream &output, const DNA &model) {
    int atom_num = 1;
    int residue_num = 1;
    output << fixed << setprecision(3);
    for (auto &chain: model.chains) {
        int flag = 0;
        for (auto &residue: chain.residues) {
            for (auto &atom: residue.atoms) {
                if (flag == 0 && (atom.name == "P" || atom.name == "O1P" || atom.name == "O2P")) 
                    continue;
                output << "ATOM" 
                       << setw(7)  << atom_num << "  "
                       << left << setw(4)  << atom.name
                       << right << setw(3) << residue.name
                       << setw(2)  << chain.name 
                       << setw(4)  << residue_num 
                       << setw(12) << atom.x 
                       << setw(8)  << atom.y 
                       << setw(8)  << atom.z 
                       << "\n";
                atom_num++;
            }
            residue_num++;
            flag++;
        }
        output << "TER" << endl;
    }
    return output;
}

void DNA::print() {
    cout << *this;
}

void DNA::write(string pdbname) {
    ofstream ofile(pdbname.c_str());
    ofile << *this;
    ofile.close();
}

} /// namespace jian


