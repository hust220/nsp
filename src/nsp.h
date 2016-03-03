#include <jian/etl.h>
#include <jian/nuc3d/Predict3D.h>
#include <jian/nuc3d/BuildLoop.h>
#include <jian/nuc3d/LM.h>
#include <jian/nuc3d/Split.h>
#include <jian/nuc2d/N2D.h>
#include <jian/nuc2d/Seq2Ss.h>
#include <jian/nuc2d/Seq2Tri.h>
#include <jian/nuc2d/Seq2Qua.h>
#include <jian/cluster/Cluster.h>
#include "extract_fragment.h"
#include <jian/nuc2d/BuildSST.h>
#include <jian/nuc2d/MergeRings.h>
#include <jian/nuc2d/GetSS.h>
#include <jian/pdb/IFModel.h>
#include <jian/dg/TestMC.h>
#include <jian/scoring/AssessPSB.h>

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
                nuc3d::Predict3D pred;
                if (par.count("par")) pred(Par(par["par"][0])); else pred(par);
            } else if (type == "anal-psb") {
                scoring::AssessPSB assess_psb;
                if (par.count("mol")) for (auto && file : par["mol"]) assess_psb.analyze(pdb::PSB(file));
            } else if (type == "extract_fragment") {
                extract_fragment(par["mol"][0], boost::lexical_cast<int>(par["len"][0]));
            } else if (type == "seq2ss") {
                nuc2d::Seq2Ss seq2ss;
                if (par.count("cutoff")) seq2ss._cutoff = boost::lexical_cast<double>(par["cutoff"][0]);
                seq2ss(par["seq"][0]);
            } else if (type == "seq2tri") {
                nuc2d::Seq2Tri seq2tri;
                if (par.count("cutoff")) seq2tri._cutoff = boost::lexical_cast<double>(par["cutoff"][0]);
                seq2tri(par["seq"][0]);
            } else if (type == "seq2qua") {
                nuc2d::Seq2Qua seq2qua;
                if (par.count("cutoff")) seq2qua._cutoff = boost::lexical_cast<int>(par["cutoff"][0]);
                seq2qua(par["seq"][0]);
            } else if (type == "dg") {
                auto mat = mat_from_file(global[1]);
                DG dg(mat);
                std::cout << dg() << std::endl;
            } else if (type == "n2d") {
                nuc2d::N2D n2d;
                if (par.count("view")) n2d.view = 1;
                if (par.count("h")) n2d.hinge_base_pair_num = boost::lexical_cast<int>(par["h"][0]);
                n2d(par["ss"][0]);
                n2d.print();
            } else if (type == "sub") {
                Model model(global[1]);
                std::vector<int> nums;
                for (int i = 2; i < global.size(); i++) {
                    std::vector<std::string> array;
                    tokenize(global[i], array, "-");
                    if (array.size() == 1) {
                        nums.push_back(boost::lexical_cast<int>(array[0]) - 1);
                    } else if (array.size() == 2) {
                        for (int j = boost::lexical_cast<int>(array[0]); j <= boost::lexical_cast<int>(array[1]); j++) {
                            nums.push_back(j - 1);
                        }
                    }
                }
                std::cout << sub(model, nums) << std::endl;
            } else if (type == "split") {
                if (par.count("par")) {
                    Par pars(par["par"][0]);
                    nuc3d::Split split(pars);
                    split();
                } else {
                    nuc3d::Split split(par);
                    split();
                }
            } else if (type == "rmsd") {
                std::cout << RMSD()(Model(global[1]), Model(global[2])) << std::endl;
            } else if (type == "lm") {
                nuc3d::LM lm; std::cout << lm(par["seq"][0], par["ss"][0])[0] << std::endl;
            } else if (type == "rebuild" && global[1] == "chain") {
                auto rna = RNA(global[2]);
                std::set<int> break_points;
                Residue old_res;
                int res_num = 0;
                for (auto &&chain: rna) {
                    for (auto &&residue: chain) {
                        if (! old_res.empty()) {
                            if (geom::distance(atom(residue, "O5*"), atom(old_res, "O3*")) > 3.5) {
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
                std::vector<Model> models;
                while (getline(ifile, line)) {
                    models.push_back(Model(boost::trim_copy(line)));
                }
                ifile.close();
                Cluster cluster(boost::lexical_cast<int>(global[2]));
                RMSD rmsd;
                cluster(models.begin(), models.end(), rmsd);
                for (auto &&clu: cluster._clusters) {
                    for (auto &&i: clu) {
                        std::cout << models[i].name << ' ';
                    }
                    std::cout << std::endl;
                }
            } else if (type == "cif") {
                Cif cif(global[1]);
                for (int i = 0; i < cif._loop["_atom_site.group_PDB"].size(); i++) {
                    std::cout << cif._loop["_atom_site.label_seq_id"][i] << ' ' << cif._loop["_atom_site.label_atom_id"][i] << ' ' << cif._loop["_atom_site.label_comp_id"][i] << ' ' << cif._loop["_atom_site.label_asym_id"][i] << ' ' << cif._loop["_atom_site.Cartn_x"][i] << ' ' << cif._loop["_atom_site.Cartn_y"][i] << ' ' << cif._loop["_atom_site.Cartn_z"][i] << std::endl;
                }
            } else if (type == "tokenize") {
                std::vector<std::string> frags;
                tokenize(global[1], frags, " ", "''\"\"");
                std::cout << frags.size() << ' ';
                std::copy(frags.begin(), frags.end(), std::ostream_iterator<std::string>(std::cout, ":"));
                std::cout << std::endl;
            } else if (type == "split-models") {
                Pdb pdb(global[1]);
                pdb.print_models();
            } else if (type == "dna") {
                auto mol = DNA(global[1]);
                Format format;
                cout << format(mol);
            } else if (type == "rna") {
                auto mol = RNA(global[1]);
                Format format;
                cout << format(mol);
            } else if (type == "sort") {
                Model mol(global[1]);
                Format format;
                std::cout << format(mol);
            } else if (type == "convert") {
                Convert cvt;
                Model mol(global[1]);
                if (par["type"][0] == "DNA") {
                    std::string seq = par["seq"][0];
                    int i = 0;
                    for (auto &&chain: mol) {
                        for (auto &&residue: chain) {
                            cvt(residue, std::string("D") + seq[i]);
                            i++;
                        }
                    }
                }
                std::cout << mol;
            } else if (type == "test") {
                auto rna = RNA(global[1]);
                int a = boost::lexical_cast<int>(global[2]);
                int b = boost::lexical_cast<int>(global[3]);
                int n = 0;
                Point p1, p2;
                for (auto &&chain: rna) for (auto &&res: chain) {
                    n++;
                    if (n == a) p1 = pos(atom(res, "C4*"));
                    else if (n == b) p2 = pos(atom(res, "C4*"));
                }
                std::cout << geom::distance(p1, p2) << std::endl;
            } else if (type == "test2") {
                MatrixXf a = MatrixXf::Zero(3, 3);
                MatrixXf b = MatrixXf::Zero(3, 3);
                MatrixXf c = MatrixXf::Zero(3, 3);
                std::cout << sum([](const MatrixXf &a){return a.rows();}, a, b, c) << std::endl;
            } else if (type == "tree") {
                nuc2d::BuildSST build_sst;
                auto sst = build_sst(par["seq"][0], par["ss"][0]);
                nuc2d::MergeRings merge_rings;
                merge_rings(sst);
            } else if (type == "test-mc") {
                dg::TestMC test_mc;
                test_mc(3);
            } else if (type == "len") {
                std::cout << num_residues(Model(global[1])) << std::endl;
            } else if (type == "ss") {
                nuc2d::GetSS get_ss;
                std::cout << get_ss(pdb::IFRNA(global[1])) << std::endl;
            } else if (type == "seq") {
                Model model(global[1]);
                std::string delimiter = "";
                if (global.size() == 3)
                    delimiter = global[2];
                for (auto &&chain: model) {
                    for (auto &&residue: chain) {
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


