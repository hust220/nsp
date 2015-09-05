#include "Split.h"

namespace jian {

namespace nuc3d {

using namespace nuc2d;

Split::Split(Par par) {
    /// set default library
    char *env = getenv("NSP");
    env != NULL or die("Please set environment 'NSP'");
    lib = env;

    if (par.count("pdb_file")) name = par["pdb_file"][0];
    if (par.count("pdb")) name = par["pdb"][0];
    if (par.count("name")) name = par["name"][0];
    if (par.count("secondary_structure")) ss = par["secondary_structure"][0];
    if (par.count("ss")) ss = par["ss"][0];
    if (par.count("library_path")) lib = par["library_path"][0];
    if (par.count("library")) lib = par["library"][0];
    if (par.count("lib")) lib = par["lib"][0];
    if (par.count("family")) lib = par["family"][0];
    if (par.count("type")) type = boost::to_upper_copy(par["type"][0]);
    if (par.count("hinge")) _hinge_size = std::stoi(par["hinge"][0]);

    lib += "/";
    lib += type;

    ss != "" || die("Split::Split error! Please provide the 2D structure");
    name != "" || die("Split::Split error! Please provide the pdb file");
    mol = Model(name);
    seq = mol.seq();
}

void Split::operator ()() {
    N2D n2d;
    n2d.hinge_base_pair_num = _hinge_size;
    n2d(seq, ss);
    splitMol(n2d.pseudo_head);
}

void Split::splitMol(loop *l) {
    /// helix
    if (l->s.head != NULL) {
        _helix_num++;

        /// write helix pdb file
        std::string helix_file_name = lib + "/helix-" + std::to_string(_hinge_size) + "/helix_" + mol.name + "_" + std::to_string(_helix_num) + ".pdb";
        std::vector<int> nums;
        for (auto b = l->s.head; b != NULL; b = b->next) {
            nums.push_back(b->res1.num - 1);
            nums.push_back(b->res2.num - 1);
        }
        Model model = mol.sub(nums);
        model.write(helix_file_name);

        /// write to information file
        std::string helix_info_file = lib + "/info-" + std::to_string(_hinge_size) + "/helix";
        std::ofstream ofile(helix_info_file.c_str(), ios::app);
        ofile << "helix_" + mol.name + "_" + std::to_string(_helix_num) << ' ' << l->s.getLen() << ' ' << l->s.seq() << ' ' << l->s.ss() << ' ' << mol.name << endl;
        ofile.close();
    }

    /// loop
    if (l->head != NULL) {
        _loop_num++;

        std::string loop_name = "loop_" + mol.name + "_" + std::to_string(_loop_num);
        std::string loop_file_name = lib + "/loop-" + std::to_string(_hinge_size) + "/" + loop_name + ".pdb";
        std::vector<int> nums;
        for (auto r = l->head; r != NULL; r = r->next) {
            nums.push_back(r->num - 1);
        }
        Model model = mol.sub(nums);
        model.write(loop_file_name);

        /// write to info file
        string l_ss = l->ss();
        std::string loop_info_file = lib + "/info-" + std::to_string(_hinge_size);
        if (std::count_if(l_ss.begin(), l_ss.end(), [](char c) {
            return c == '[' || c == ']';
        })) {
            loop_info_file += "/pseudo_knot";
        } else {
            loop_info_file += "/loop-" + std::to_string(l->hinges.size());
        }

        std::ofstream ofile(loop_info_file.c_str(), ios::app);
        ofile << loop_name << ' ' << l->seq() << ' ' << l->ss() << ' ' << family << endl;
        ofile.close();

    }

    /// son
    if (l->son != NULL) {
        splitMol(l->son);
    }

    /// brother
    if (l->brother != NULL) {
        splitMol(l->brother);
    }
}

} /// namespace nuc3d

} /// namespace jian


