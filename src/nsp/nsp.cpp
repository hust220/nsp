#include "../nsp.h"

int main(int argc, char **argv) {
    jian::Par par(argc, argv);
    std::string par_type = boost::to_lower_copy(par["global"][0]);
    if (par_type == "3drna") {
        if (par.count("par")) {
            jian::Par pars(par["par"][0]);
            jian::nuc3d::Assemble ass(pars);
            ass();
        } else {
            jian::nuc3d::Assemble ass(par);
            ass();
        }
    } else if (par_type == "rmsd") {
        std::cout << jian::pdb::RMSD()(jian::Model(argv[2]), jian::Model(argv[3])) << std::endl;
    } else if (!strcmp(argv[1], "lm")) {
        std::string type = par["type"][0];
        std::string seq = par["seq"][0];
        std::string ss = par["ss"][0];
        boost::to_upper(type);
        int view = (par.count("view") ? 1 : 0);

        jian::nuc3d::LoopModelling2 lm(type);
        lm._view = view;
        std::cout << lm(seq, ss) << std::endl;
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
    //} else if (!strcmp(argv[1], "score")) {
    } else if (par["global"][0] == "score") {
        jian::Score score;
        if (par.count("list")) {
            std::ifstream ifile(par["list"][0].c_str());
            std::string line;
            while (std::getline(ifile, line)) {
                std::cout << score(jian::RNA(boost::trim_copy(line))) << std::endl;
            }
            ifile.close();
        } else {
            std::cout << score(par["global"][1]) << std::endl;
        }
    } else if (par["global"][0] == "train_junction") {
        jian::nuc3d::JunctBuild jb;
        jb.train(jian::Model(par["pdb"][0]), par["ss"][0]);
    } else if (!strcmp(argv[1], "-mol2d")) {
        jian::nuc2d::N2D mol2d(argv[2], 1);
        mol2d.print();
    } else if (!strcmp(argv[1], "-d5p")) {
        jian::D5P model(argv[2]);
        cout << model << endl;
    } else if (!strcmp(argv[1], "-dna")) {
        jian::DNA model(argv[2]);
        cout << model << endl;
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
        jian::SortAtoms()(mol);
        std::cout << mol;
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
    } else if (!strcmp(argv[1], "test")) {
//        MatrixXf a(3, 3);
//        a << 0, 0, 0,
//             0, 1, 0,
//             0, 0, 1;
//
//        std::vector<std::tuple<int, int, double, double>> dists;
//        for (int i = 0; i < a.rows(); i++) {
//            for (int j = 0; j < a.rows(); j++) {
//                double d = (a.row(i) - a.row(j)).squaredNorm();
//                dists.push_back(std::make_tuple(i, j, d, d));
//            }
//        }
//
//        jian::dg::Dgsol dgsol;
//        std::cout << dgsol(dists, 3) << std::endl;
    } else if (!strcmp(argv[1], "seq")) {
        jian::Model model(argv[2]);
        for (auto &&chain: model.chains) {
            for (auto &&residue: chain.residues) {
                std::cout << residue.name.substr(residue.name.size() - 1, 1);
            }
        }
        std::cout << std::endl;
    }
    return 0;
}



