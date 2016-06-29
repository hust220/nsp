#include <iostream>
#include <map>
#include <string>
#include <functional>
#include <deque>
#include <Eigen/Dense>
#include "nsp.hpp"
#include <jian/pdb/Model.hpp>
#include <jian/geom/suppos.hpp>
#include <jian/utils/file.hpp>
#include <jian/utils/Cluster.hpp>

namespace jian {

template<typename T, typename U>
static void print_clusters(T &&cluster, U &&names) {
    int n = 0;
    for (auto && clu : cluster._clusters) {
        std::cout << "Cluster " << n+1 << " (size: " << clu.size() << "): ";
        for (auto && i : clu) {
            std::cout << names[i] << ' ';
        }
        std::cout << std::endl;
        n++;
    }
}

template<typename T>
static void print_clusters(T &&cluster) {
    int n = 0;
    for (auto && clu : cluster._clusters) {
        std::cout << "Cluster " << n+1 << " (size: " << clu.size() << "): ";
        for (auto && i : clu) {
            std::cout << i << ' ';
        }
        std::cout << std::endl;
        n++;
    }
}

static Eigen::MatrixXd * model_to_matrix(const Model &model) {
    int len = num_residues(model);
    Eigen::MatrixXd *mat = new Eigen::MatrixXd(len, 3);
    EACH_RES(model, 
        auto &&atom = RES["C4*"]; 
        for (int i = 0; i < 3; i++) {
            (*mat)(N_RES, i) = atom[i];
        }
    );
    return mat;
}

static Eigen::MatrixXd * model_to_matrix_aa(const Model &model) {
    int len = 0;
    for (auto && chain : model) {
        for (auto && res : chain) {
            for (auto && atom : res) {
                if (std::find(atom.name.begin(), atom.name.end(), 'P') == atom.name.end()) {
                    len++;
                }
            }
        }
    }
    Eigen::MatrixXd *mat = new Eigen::MatrixXd(len, 3);
    int i = 0;
    for (auto && chain : model) {
        for (auto && res : chain) {
            for (auto && atom : res) {
                if (std::find(atom.name.begin(), atom.name.end(), 'P') == atom.name.end()) {
                    for (int k = 0; k < 3; k++) {
                        (*mat)(i, k) = atom[k];
                    }
                    i++;
                }
            }
        }
    }
    return mat;
}

static double dist(Eigen::MatrixXd *m1, Eigen::MatrixXd *m2) {
    auto sp = geom::suppos(*m1, *m2);
    return sp.rmsd;
}

REGISTER_NSP_COMPONENT(cluster) {
    using method_t = std::function<Eigen::MatrixXd*(const Model &)>;
    static std::map<std::string, method_t> methods {
        {"c4", model_to_matrix},
        {"aa", model_to_matrix_aa}
    };
    Cluster cluster(JN_INT(par["k"][0]));
    if (par.has("list")) {
        std::string method_name = "aa";
        if (par.has("method")) method_name = par["method"][0];
        method_t &method = methods[method_name];
        std::deque<Eigen::MatrixXd *> mats; 
        std::deque<std::string> names;
        std::cout << "Reading molecules..." << std::endl;
        EACH_SPLIT_LINE(par["list"][0].c_str(), " ",
            std::cout << "read " << F[0] << std::endl;
            Model model(F[0]);
            mats.push_back(method(model));
            names.push_back(model.name);
        );
        std::cout << "Clustering..." << std::endl;
        cluster(mats.begin(), mats.end(), dist);
        for (auto && i : mats) delete i;
        print_clusters(cluster, names);
    } else if (par.has("matrix")) {
        auto &&mat = mat_from_file(par["matrix"][0]);
        cluster(mat);
        print_clusters(cluster);
    }
}

} // namespace jian

