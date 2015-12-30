#include "../nuc3d/Predict.h"
//#include "../nuc3d/LM2.h"
//#include "../nuc3d/LoopModelling2.h"
#include "../nuc3d/BuildJunction.h"
//#include "../nuc3d/BuildTriplex.h"
#include "../nuc3d/Split.h"
#include "../nuc2d/N2D.h"
#include "../nuc2d/Seq2Ss.h"
#include "../nuc2d/Seq2Tri.h"
#include "../nuc2d/Seq2Qua.h"
#include "../cluster/Cluster.h"
//#include "../scoring/util.h"
#include "extract_fragment.h"
#include "../nuc2d/BuildSST.h"
#include "../nuc2d/MergeRings.h"
#include "../nuc2d/GetSS.h"
#include "../pdb/IFModel.h"
#include "../dg/TestMC.h"

namespace jian {

class NSP {
public:
    Par par;
    std::vector<std::string> global;
    std::string type;

    NSP(int argc, char **argv) : par(argc, argv) {
        global = par["global"];
        type = boost::to_lower_copy(global[0]);
    }

    void operator ()() { 
        try {
            if (type == "3drna") {
                if (par.count("par")) {
                    jian::Par pars(par["par"][0]);
                    jian::nuc3d::Predict pred(pars);
                    pred();
                } else {
                    jian::nuc3d::Predict pred(par);
                    pred();
                }
            } else if (type == "extract_fragment") {
                jian::extract_fragment(par["mol"][0], std::stoi(par["len"][0]));
//            } else if (type == "lm2") {
//                jian::nuc3d::LM2 lm2(par);
//                for (int i = 0; i < lm2._num; i++) {
//                    lm2().write(lm2._name + "-" + std::to_string(i + 1) + ".pdb");
//                }
            } else if (type == "seq2ss") {
                jian::nuc2d::Seq2Ss seq2ss;
                if (par.count("cutoff")) seq2ss._cutoff = std::stof(par["cutoff"][0]);
                seq2ss(par["seq"][0]);
            } else if (type == "seq2tri") {
                jian::nuc2d::Seq2Tri seq2tri;
                if (par.count("cutoff")) seq2tri._cutoff = std::stof(par["cutoff"][0]);
                seq2tri(par["seq"][0]);
            } else if (type == "seq2qua") {
                jian::nuc2d::Seq2Qua seq2qua;
                if (par.count("cutoff")) seq2qua._cutoff = std::stoi(par["cutoff"][0]);
                seq2qua(par["seq"][0]);
            } else if (type == "dg") {
                auto mat = jian::mat::mat_from_file(global[1]);
                jian::DG dg(mat);
                std::cout << dg() << std::endl;
            } else if (type == "n2d") {
                jian::nuc2d::N2D n2d;
                if (par.count("view")) n2d.view = 1;
                if (par.count("h")) n2d.hinge_base_pair_num = std::stoi(par["h"][0]);
                n2d(par["ss"][0]);
                n2d.print();
            } else if (type == "sub") {
                jian::Model model(global[1]);
                std::vector<int> nums;
                for (int i = 2; i < global.size(); i++) {
                    std::vector<std::string> array;
                    jian::tokenize(global[i], array, "-");
                    if (array.size() == 1) {
                        nums.push_back(std::stoi(array[0]) - 1);
                    } else if (array.size() == 2) {
                        for (int j = std::stoi(array[0]); j <= std::stoi(array[1]); j++) {
                            nums.push_back(j - 1);
                        }
                    }
                }
                std::cout << model.sub(nums) << std::endl;
            } else if (type == "split") {
                if (par.count("par")) {
                    jian::Par pars(par["par"][0]);
                    jian::nuc3d::Split split(pars);
                    split();
                } else {
                    jian::nuc3d::Split split(par);
                    split();
                }
            } else if (type == "rmsd") {
                std::cout << jian::pdb::RMSD()(jian::Model(global[1]), jian::Model(global[2])) << std::endl;
//            } else if (type == "lm") {
//                std::string type = par["type"][0];
//                std::string seq = par["seq"][0];
//                std::string ss = par["ss"][0];
//                boost::to_upper(type);
//                int view = (par.count("view") ? 1 : 0);
//
//                jian::nuc3d::LoopModelling2 lm(type);
//                lm._view = view;
//                std::cout << lm(seq, ss) << std::endl;
            } else if (type == "rebuild" && global[1] == "chain") {
                jian::RNA rna(global[2]);
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
                for (int i = 0; i < global[3].size(); i++) {
                    if (global[3][i] == '&') {
                        continue;    
                    }
                    ss += global[3][i];
                    if (break_points.count(res_num)) {
                        ss += '&';
                    }
                    res_num++;
                }
                std::cout << ss;
            } else if (type == "cluster") {
                std::ifstream ifile(global[1].c_str());
                std::string line;
                std::vector<jian::Model> models;
                while (getline(ifile, line)) {
                    models.push_back(jian::Model(boost::trim_copy(line)));
                }
                ifile.close();
                jian::Cluster cluster(std::stoi(global[2]));
                jian::pdb::RMSD rmsd;
                cluster(models.begin(), models.end(), rmsd);
                for (auto &&clu: cluster._clusters) {
                    for (auto &&i: clu) {
                        std::cout << models[i].name << ' ';
                    }
                    std::cout << std::endl;
                }
//            } else if (type == "train") {
//                jian::scoring::Train train(par);
//                train();
//            } else if (type == "score") {
//                jian::scoring::Score score(par);
//                score();
            } else if (type == "train_junction") {
                jian::nuc3d::JunctBuild jb;
                jb.train(jian::Model(par["pdb"][0]), par["ss"][0]);
            } else if (type == "cif") {
                jian::Cif cif(global[1]);
                for (int i = 0; i < cif._loop["_atom_site.group_PDB"].size(); i++) {
                    std::cout << cif._loop["_atom_site.label_seq_id"][i] << ' ' << cif._loop["_atom_site.label_atom_id"][i] << ' ' << cif._loop["_atom_site.label_comp_id"][i] << ' ' << cif._loop["_atom_site.label_asym_id"][i] << ' ' << cif._loop["_atom_site.Cartn_x"][i] << ' ' << cif._loop["_atom_site.Cartn_y"][i] << ' ' << cif._loop["_atom_site.Cartn_z"][i] << std::endl;
                }
            } else if (type == "tokenize") {
                std::vector<std::string> frags;
                jian::tokenize(global[1], frags, " ", "''\"\"");
                std::cout << frags.size() << ' ';
                std::copy(frags.begin(), frags.end(), std::ostream_iterator<std::string>(std::cout, ":"));
                std::cout << std::endl;
            } else if (type == "split-models") {
                jian::Pdb pdb(global[1]);
                pdb.print_models();
            } else if (type == "dna") {
                jian::DNA mol(global[1]);
                jian::pdb::Format format;
                cout << format(mol);
            } else if (type == "rna") {
                jian::RNA mol(global[1]);
                jian::pdb::Format format;
                cout << format(mol);
//            } else if (type == "triplex") {
//                auto seq = par["seq"][0];
//                jian::nuc2d::Seq2Tri seq2tri;
//                auto info_list = seq2tri(seq);
//                jian::nuc3d::BuildTriplex build_triplex;
//                build_triplex(seq, info_list[0], 1);
            } else if (type == "sort") {
                jian::Model mol(global[1]);
                jian::pdb::Format format;
                std::cout << format(mol);
            } else if (type == "convert") {
                jian::Convert cvt;
                jian::Model mol(global[1]);
                if (par["type"][0] == "DNA") {
                    std::string seq = par["seq"][0];
                    int i = 0;
                    for (auto &&chain: mol.chains) {
                        for (auto &&residue: chain.residues) {
                            cvt(residue, std::string("D") + seq[i]);
                            i++;
                        }
                    }
                }
                std::cout << mol;
            } else if (type == "test") {
                jian::RNA rna(global[1]);
                int a = std::stoi(global[2]);
                int b = std::stoi(global[3]);
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
            } else if (type == "test2") {
                MatrixXf a = MatrixXf::Zero(3, 3);
                MatrixXf b = MatrixXf::Zero(3, 3);
                MatrixXf c = MatrixXf::Zero(3, 3);
                std::cout << jian::mat::sum([](const MatrixXf &a){return a.rows();}, a, b, c) << std::endl;
            } else if (type == "tree") {
                jian::nuc2d::BuildSST build_sst;
                auto sst = build_sst(par["seq"][0], par["ss"][0]);
                jian::nuc2d::MergeRings merge_rings;
                merge_rings(sst);
            } else if (type == "test-mc") {
                jian::dg::TestMC test_mc;
                test_mc(3);
            } else if (type == "len") {
                std::cout << jian::Model(global[1]).res_nums() << std::endl;
            } else if (type == "ss") {
                jian::nuc2d::GetSS get_ss;
                std::cout << get_ss(jian::pdb::IFModel(global[1])) << std::endl;
            } else if (type == "seq") {
                jian::Model model(global[1]);
                std::string delimiter = "";
                if (global.size() == 3)
                    delimiter = global[2];
                for (auto &&chain: model.chains) {
                    for (auto &&residue: chain.residues) {
                        std::cout << residue.name << delimiter;
                    }
                }
                std::cout << std::endl;
            }
        } catch (const char *inf) {
            std::cout << inf << std::endl;
        }
    }

}; // class NSP


} // namespace jian


