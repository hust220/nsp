#include "nsp.hpp"
#include <jian/pdb/Model.hpp>
#include <jian/cluster/Cluster.hpp>
#include <jian/geom.hpp>

namespace jian {

MatrixXd * model_to_matrix(const Model &model) {
    int len = num_residues(model);
    MatrixXd *mat = new MatrixXd(len, 3);
    EACH_RES(model, auto &&atom = RES["C4*"]; FOR((i, 3), (*mat)(N_RES, i) = atom[i]));
    return mat;
}

double dist(MatrixXd *m1, MatrixXd *m2) {
    auto sp = geom::suppos(*m1, *m2);
    return sp.rmsd;
}

REGISTER_NSP_COMPONENT(cluster) {
    std::deque<MatrixXd *> mats; std::deque<std::string> names;
    SEELN("Reading molecules...");
    EACH_SPLIT_LINE(par["list"][0].c_str(), " ",
        Model model(F[0]);
        mats.push_back(model_to_matrix(model));
        names.push_back(model.name);
    );
    Cluster cluster(JN_INT(par["k"][0]));
    cluster(mats.begin(), mats.end(), dist);
    EACH(i, mats, delete i);
    EACH((clu, n), cluster._clusters, SEE("Cluster ", n+1, " (", clu.size(), "): "); EACH(i, clu, SEE(names[i], ' ')); SEE('\n'));
}

} // namespace jian

