#ifndef JUNCTBUILD_H
#define JUNCTBUILD_H

#include "../MC.h"
#include "../Pdb.h"

namespace jian {

namespace nuc3d {

class JunctBuild: public mc::System {
public:
    JunctBuild() {}
    Model operator ()(std::string seq, std::string ss);
    void train(const Model &model, std::string ss);
    std::pair<std::vector<Model>, std::vector<int>> get_helices(const Model &model, std::string ss);
    std::tuple<Vector3f, Vector3f, Vector3f> helix_par(const Model &model);

    void move();
    void rollback();
    double energy();
    double min_energy();
    void set_min_state();

    int _helix_len = 6;
    std::vector<std::tuple<Vector3f, Vector3f, Vector3f>> _helices;
};

} /// namespace nuc3d

} /// namespace jian



#endif 

