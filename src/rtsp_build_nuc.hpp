#pragma once

#include "dg.hpp"
#include "pdb.hpp"

BEGIN_JN

class BuildNuc {
public:
//    BuildNuc(string type = "RNA");
//    Residue operator() (const string &name, const Mat &scaffold);
//    Residue make_residue(const string &name, const Mat &coords);
    using Mat = MatrixXd;

    int _view = 0;

    S _type = "RNA";
    S _lib;
    std::map<string, Mat> _base_aa_par;
    std::map<string, Mat> _base_cg_par;
    std::map<string, Mat> _phos_sugar_par;
    std::map<string, vector<string>> _atom_names;

    BuildNuc(string type = "RNA"): _type(type) {
        /// Set _lib
        char *lib = getenv("NSP");
        if (lib == NULL) throw "Please designate environment variable 'NSP'";
        _lib = string() + lib + "/" + _type;

        /// Set _base_aa_par and _base_cg_par
        std::map<string, int> temp_map;
        std::map<string, vector<int>> temp_map2;
        if (_type == "RNA") {
            temp_map = {{"A", 11}, {"U", 9}, {"G", 12}, {"C", 9}};
            temp_map2 = {{"A", {0, 6, 8}}, {"U", {0, 3, 6}}, {"G", {0, 6, 9}}, {"C", {0, 3, 6}}};
        } else if (_type == "DNA") {
            temp_map = {{"DA", 11}, {"DT", 10}, {"DG", 12}, {"DC", 9}};
            temp_map2 = {{"DA", {0, 6, 8}}, {"DT", {0, 3, 6}}, {"DG", {0, 6, 9}}, {"DC", {0, 3, 6}}};
        }
        ifstream ifile;
        string par_file;
        for (auto &temp_pair: temp_map) {
            _base_aa_par[temp_pair.first].resize(temp_pair.second, 3);
            par_file = _lib + "/base/" + temp_pair.first;
            ifile.open(par_file.c_str());
            for (int i = 0; i < temp_pair.second; i++) {
                for (int j = 0; j < 3; j++) {
                    ifile >> _base_aa_par[temp_pair.first](i, j);
                }
            }
            ifile.close();

            /// Set _base_cg_par
            _base_cg_par[temp_pair.first].resize(3, 3);
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    _base_cg_par[temp_pair.first](i, j) =
                        _base_aa_par[temp_pair.first](temp_map2[temp_pair.first][i], j);
                }
            }
        }

        /// Set _phos_sugar_par
        int temp_int;
        vector<string> name_list;
        if (_type == "RNA") {
            name_list = {"A", "U", "G", "C"};
            temp_int = 12;
        } else if (_type == "DNA") {
            name_list = {"DA", "DT", "DG", "DC"};
            temp_int = 11;
        }
        for (auto &name: name_list) {
            _phos_sugar_par[name].resize(temp_int, temp_int);
            par_file = _lib + "/pars/" + name + ".avg";
            ifile.open(par_file.c_str());
            for (int i = 0; i < temp_int; i++) {
                for (int j = 0; j < temp_int; j++) {
                    ifile >> _phos_sugar_par[name](i, j);
                }
            }
            ifile.close();

            par_file = _lib + "/pars/" + name + ".stdev";
            ifile.open(par_file.c_str());
            for (int i = 0; i < temp_int; i++) {
                for (int j = 0; j < temp_int; j++) {
                    double temp_double;
                    ifile >> temp_double;
                    _phos_sugar_par[name](i, j) += (i > j ? -3 : 3) * temp_double;
                }
            }
            ifile.close();
        }

        /// Set _atom_names
        if (_type == "RNA") {
            _atom_names = {
                {"A", {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "O2*", "C1*", "N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"}},
                {"U", {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "O2*", "C1*", "N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"}},
                {"G", {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "O2*", "C1*", "N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"}},
                {"C", {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "O2*", "C1*", "N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"}}
            };
        } else if (_type == "DNA") {
            _atom_names = {
                {"DA", {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "C1*", "N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"}},
                {"DT", {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "C1*", "N1", "C2", "O2", "N3", "C4", "O4", "C5", "C7", "C6"}},
                {"DG", {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "C1*", "N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"}},
                {"DC", {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "C1*", "N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"}}
            };
        }
    }

    Residue operator ()(const string &name, const Mat &scaffold) {
        /// Base coordinates
        auto base_coords = _base_aa_par[name];
        auto temp_coords = _base_cg_par[name];
        Mat temp_scaffold(3, 3);
        temp_scaffold = scaffold.block<3, 3>(2, 0);
        geom::suppos(base_coords, temp_coords, temp_scaffold);

        /// Phosphate and sugar coordinates
        auto par = _phos_sugar_par[name];
        vector<int> vec{4, 8, 10};
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double dist = (scaffold.row(i) - scaffold.row(j)).norm();
                par(vec[i], vec[j]) = par(vec[j], vec[i]) = dist;
            }
        }
        DG dg(par);
        dg.log.set_display(_view);
        auto phos_sugar_coords = dg();
        temp_coords.resize(3, 3);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                temp_coords(i, j) = phos_sugar_coords(vec[i], j);
            }
        }
        temp_scaffold = scaffold.block<3, 3>(0, 0);
        geom::suppos(phos_sugar_coords, temp_coords, temp_scaffold);

        /// Merge base coordinates and phosphate-sugar coordinates
        Mat coords(phos_sugar_coords.rows() + base_coords.rows() - 1, 3);
        coords << phos_sugar_coords.block(0, 0, phos_sugar_coords.rows() - 1, 3), base_coords;

        return make_residue(name, coords);
    }

    Residue make_residue(const string &name, const Mat &coords) {
        Residue residue;
        residue.name = name;
        for (int i = 0; i < coords.rows(); i++) {
            residue.push_back(Atom(_atom_names[name][i], coords(i, 0), coords(i, 1), coords(i, 2)));
        }
        return residue;
    }

};

END_JN

