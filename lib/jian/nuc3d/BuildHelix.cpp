#include <iostream>
#include <fstream>
#include <string>
#include "Connect.hpp"
#include "BuildHelix.hpp"
#include "../utils/Env.hpp"

namespace jian {

Model build_helix(std::string seq) {
    std::string lib = Env::lib() + "/RNA";
    if (seq.size() < 2 || seq.size() % 2 == 1) {
        throw std::string("jian::BuildHelix::operator (std::string) error! Unreasonable length.\nSequence: ") + seq;
    } else if (seq.size() == 2) {
        std::string file_name = lib + "/basepair/" + seq + ".pdb";
        std::ifstream ifile(file_name.c_str());
        if (!ifile) {
            file_name = lib + "/basepair/XX.pdb";
        }
        ifile.close();
        return Model(file_name);
    } else if (seq.size() == 4) {
        std::string file_name = lib + "/basepair/" + seq + ".pdb";
        std::ifstream ifile(file_name.c_str());
        if (!ifile) {
            file_name = lib + "/basepair/XXXX.pdb";
        }
        ifile.close();
        return Model(file_name);
    } else {
        std::string file_name = lib + "/basepair/" + seq.substr(0, 2) + seq.substr(seq.size() - 2, 2) + ".pdb";
        std::ifstream ifile(file_name.c_str());
        if (!ifile) {
            file_name = lib + "/basepair/XXXX.pdb";
        }
        ifile.close();
        Connect connect;
        connect._hinge_size = 1;
        return connect(Model(file_name), build_helix(seq.substr(1, seq.size() - 2)), 1, 2);
    }
}

} // namespace jian

