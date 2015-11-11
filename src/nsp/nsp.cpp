#include <nuc3d/util.h>
#include <cluster/Cluster.h>
#include <scoring/util.h>

int main(int argc, char **argv) {
    jian::Par par(argc, argv);
    if (boost::to_lower_copy(par["global"][0]) == "3drna") {
        if (par.count("par")) {
            jian::Par pars(par["par"][0]);
            jian::nuc3d::Assemble ass(pars);
            ass();
        } else {
            jian::nuc3d::Assemble ass(par);
            ass();
        }
    } else if (boost::to_lower_copy(par["global"][0]) == "lm2") {
        jian::nuc3d::LM2 lm2;
        lm2(par["seq"][0], par["ss"][0]);
    } else if (boost::to_lower_copy(par["global"][0]) == "seq2ss") {
        jian::nuc2d::Seq2Ss seq2ss;
        if (par.count("cutoff")) seq2ss._cutoff = std::stof(par["cutoff"][0]);
        seq2ss(par["seq"][0]);
    } else if (boost::to_lower_copy(par["global"][0]) == "seq2tri") {
        jian::nuc2d::Seq2Tri seq2tri;
        if (par.count("cutoff")) seq2tri._cutoff = std::stoi(par["cutoff"][0]);
        seq2tri(par["seq"][0]);
    } else if (boost::to_lower_copy(par["global"][0]) == "seq2qua") {
        jian::nuc2d::Seq2Qua seq2qua;
        if (par.count("cutoff")) seq2qua._cutoff = std::stoi(par["cutoff"][0]);
        seq2qua(par["seq"][0]);
    } else if (boost::to_lower_copy(par["global"][0]) == "dg") {
        auto mat = jian::mat_from_file(par["global"][1]);
        jian::DG dg(mat);
        std::cout << dg() << std::endl;
    } else if (boost::to_lower_copy(par["global"][0]) == "n2d") {
        jian::nuc2d::N2D n2d;
        if (par.count("h")) n2d.hinge_base_pair_num = std::stoi(par["h"][0]);
        n2d(par["ss"][0]);
        n2d.print();
    } else if (boost::to_lower_copy(par["global"][0]) == "sub") {
        jian::Model model(par["global"][1]);
        std::vector<int> nums;
        for (int i = 2; i < par["global"].size(); i++) {
            std::vector<std::string> array;
            jian::tokenize(par["global"][i], array, "-");
            if (array.size() == 1) {
                nums.push_back(std::stoi(array[0]) - 1);
            } else if (array.size() == 2) {
                for (int j = std::stoi(array[0]); j <= std::stoi(array[1]); j++) {
                    nums.push_back(j - 1);
                }
            }
        }
        std::cout << model.sub(nums) << std::endl;
    } else if (boost::to_lower_copy(par["global"][0]) == "split") {
        if (par.count("par")) {
            jian::Par pars(par["par"][0]);
            jian::nuc3d::Split split(pars);
            split();
        } else {
            jian::nuc3d::Split split(par);
            split();
        }
    } else if (boost::to_lower_copy(par["global"][0]) == "rmsd") {
        std::cout << jian::pdb::RMSD()(jian::Model(argv[2]), jian::Model(argv[3])) << std::endl;
    } else if (boost::to_lower_copy(par["global"][0]) == "lm") {
        std::string type = par["type"][0];
        std::string seq = par["seq"][0];
        std::string ss = par["ss"][0];
        boost::to_upper(type);
        int view = (par.count("view") ? 1 : 0);

        jian::nuc3d::LoopModelling2 lm(type);
        lm._view = view;
        std::cout << lm(seq, ss) << std::endl;
    } else if (par["global"][0] == "rebuild" && par["global"][1] == "chain") {
        jian::RNA rna(par["global"][2]);
        std::set<int> break_points;
        jian::Residue old_res;
        int res_num = 0;
        for (auto &&chain: rna.chains) {
            for (auto &&residue: chain.residues) {
                if (! old_res.atoms.empty()) {
                    if (residue["O5*"].dist(old_res["O3*"]) > 3.5) {
                        break_points.insert(res_num - 1);
                    }
                }
                old_res = residue;
                res_num++;
            }
        }
        std::string ss;
        res_num = 0;
        for (int i = 0; i < par["global"][3].size(); i++) {
            if (par["global"][3][i] == '&') {
                continue;    
            }
            ss += par["global"][3][i];
            if (break_points.count(res_num)) {
                ss += '&';
            }
            res_num++;
        }
        std::cout << ss;
    } else if (!strcmp(argv[1], "cluster")) {
        std::ifstream ifile(argv[2]);
        std::string line;
        std::vector<jian::Model> models;
        while (getline(ifile, line)) {
            models.push_back(jian::Model(boost::trim_copy(line)));
        }
        ifile.close();
        jian::Cluster cluster(std::stoi(argv[3]));
        jian::pdb::RMSD rmsd;
        cluster(models.begin(), models.end(), rmsd);
        for (auto &&clu: cluster._clusters) {
            for (auto &&i: clu) {
                std::cout << models[i].name << ' ';
            }
            std::cout << std::endl;
        }
    } else if (par["global"][0] == "train") {
        jian::scoring::Train train(par);
        train();
    } else if (par["global"][0] == "score") {
        jian::scoring::Score score(par);
        score();
//        if (par.count("list")) {
//            std::ifstream ifile(par["list"][0].c_str());
//            std::string line;
//            while (std::getline(ifile, line)) {
//                std::cout << score(jian::RNA(boost::trim_copy(line))) << std::endl;
//            }
//            ifile.close();
//        } else {
//            std::cout << score(par["global"][1]) << std::endl;
//        }
    } else if (par["global"][0] == "train_junction") {
        jian::nuc3d::JunctBuild jb;
        jb.train(jian::Model(par["pdb"][0]), par["ss"][0]);
    } else if (par["global"][0] == "cif") {
        jian::Cif cif(par["global"][1]);
        for (int i = 0; i < cif._loop["_atom_site.group_PDB"].size(); i++) {
            std::cout << cif._loop["_atom_site.label_seq_id"][i] << ' ' << cif._loop["_atom_site.label_atom_id"][i] << ' ' << cif._loop["_atom_site.label_comp_id"][i] << ' ' << cif._loop["_atom_site.label_asym_id"][i] << ' ' << cif._loop["_atom_site.Cartn_x"][i] << ' ' << cif._loop["_atom_site.Cartn_y"][i] << ' ' << cif._loop["_atom_site.Cartn_z"][i] << std::endl;
        }
    } else if (par["global"][0] == "tokenize") {
        std::vector<std::string> frags;
        jian::tokenize(par["global"][1], frags, " ", "''\"\"");
        std::cout << frags.size() << ' ';
        std::copy(frags.begin(), frags.end(), std::ostream_iterator<std::string>(std::cout, ":"));
        std::cout << std::endl;
    } else if (!strcmp(argv[1], "-mol2d")) {
        jian::nuc2d::N2D mol2d(argv[2], 1);
        mol2d.print();
    } else if (!strcmp(argv[1], "-d5p")) {
        jian::D5P model(argv[2]);
        cout << model << endl;
    } else if (boost::to_lower_copy(par["global"][0]) == "dna") {
        jian::DNA mol(par["global"][1]);
        jian::pdb::Format format;
        cout << format(mol);
    } else if (boost::to_lower_copy(par["global"][0]) == "rna") {
        jian::RNA mol(par["global"][1]);
        jian::pdb::Format format;
        cout << format(mol);
    } else if (!strcmp(argv[1], "-test")) {
        jian::D5P model(argv[2]);
        for (auto &chain: model.chains) {
            for (auto &residue: chain.residues) {
                for (auto &atom: residue.atoms) {
                    for (auto &chain2: model.chains) {
                        for (auto &residue2: chain2.residues) {
                            for (auto &atom2: residue2.atoms) {
                                std::cout << atom.pos().dist(atom2.pos()) << ' ';
                            }
                        }
                    }
                    std::cout << std::endl;
                }
            }
        }
    } else if (!strcmp(argv[1], "test:aa")) {
        jian::Model model(argv[2]);
        for (auto &chain: model.chains) {
            for (auto &residue: chain.residues) {
                if (residue.name == argv[3]) {
                    for (auto &atom: residue.atoms) {
                        if (atom.name == argv[4]) {
                            for (auto &atom2: residue.atoms) {
                                if (atom2.name == argv[5]) {
                                    std::cout << atom.pos().dist(atom2.pos()) << std::endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    } else if (!strcmp(argv[1], "-temp")) {
        jian::D5P mol(argv[2]);
        jian::D5P new_mol;
        new_mol.chains.push_back(jian::Chain());
        for (int i = 0; i < 7; i++) {
            new_mol.chains[0].residues.push_back(mol.chains[0].residues[i]);
            new_mol.chains[0].residues.push_back(mol.chains[0].residues[13 - i]);
            new_mol.chains[0].residues.push_back(mol.chains[0].residues[14 + i]);
        }
        std::cout << new_mol;
    } else if (!strcmp(argv[1], "sort")) {
        jian::Model mol(argv[2]);
        jian::pdb::Format format;
        std::cout << format(mol);
    } else if (!strcmp(argv[1], "convert")) {
        jian::Convert cvt;
        jian::Model mol(argv[2]);
        if (!strcmp(argv[3], "-type")) {
            if (!strcmp(argv[4], "DNA")) {
                if (!strcmp(argv[5], "-seq")) {
                    std::string seq = argv[6];
                    int i = 0;
                    for (auto &&chain: mol.chains) {
                        for (auto &&residue: chain.residues) {
                            cvt(residue, std::string("D") + seq[i]);
                            i++;
                        }
                    }
                }
            }
        }
        std::cout << mol;
    } else if (boost::to_lower_copy(par["global"][0]) == "test") {
        jian::RNA rna(par["global"][1]);
        int a = std::stoi(par["global"][2]);
        int b = std::stoi(par["global"][3]);
        int n = 0;
        jian::Point p1, p2;
        for (auto &&chain: rna) {
            for (auto &&res: chain) {
                n++;
                if (n == a) {
                    p1 = res["C4*"].pos();
                } else if (n == b) {
                    p2 = res["C4*"].pos();
                }
            }
        }
        std::cout << p1.dist(p2) << std::endl;
    } else if (boost::to_lower_copy(par["global"][0]) == "test2") {
        jian::RNA rna(par["global"][1]);
        int a = std::stoi(par["global"][2]);
        int b = std::stoi(par["global"][3]);
        int c = std::stoi(par["global"][4]);
        int d = std::stoi(par["global"][5]);
        int n = 0;
        jian::Point p1, p2, p3, p4;
        for (auto &&chain: rna) {
            for (auto &&res: chain) {
                n++;
                if (n == a) {
                    p1 = res["C4*"].pos();
                } else if (n == b) {
                    p2 = res["C4*"].pos();
                } else if (n == c) {
                    p3 = res["C4*"].pos();
                } else if (n == d) {
                    p4 = res["C4*"].pos();
                }
            }
        }
        std::cout << jian::Point::chirality(p1, p2, p3, p4) << std::endl;
        std::cout << jian::Point::dihedral(p1, p2, p3, p4) << std::endl;
    } else if (boost::to_lower_copy(par["global"][0]) == "len") {
        std::cout << jian::Model(par["global"][1]).res_nums() << std::endl;
    } else if (!strcmp(argv[1], "seq")) {
        jian::Model model(argv[2]);
        std::string delimiter = "";
        if (par["global"].size() == 3)
            delimiter = par["global"][2];
        for (auto &&chain: model.chains) {
            for (auto &&residue: chain.residues) {
                std::cout << residue.name << delimiter;
            }
        }
        std::cout << std::endl;
    }
    return 0;
}



