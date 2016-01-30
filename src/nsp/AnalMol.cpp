#include <jian/util/std.h>
#include <jian/fpl.h>
#include <jian/pdb/Model.h>
#include <jian/geom.h>

int main(int argc, char **argv) {
    jian::Par par(argc, argv);
    jian::Model model(par["pdb"][0]);
    auto ls = jian::fpl::map<std::vector>([](const auto &s){return std::stoi(s);}, par["num"]);
    std::vector<jian::Atom> vec(ls.size());
    model.each_res([&](const auto &res, int index){
        auto result = std::find(ls.begin(), ls.end(), index);
        if (result != ls.end()) {
            vec[std::distance(ls.begin(), result)] = res["C4*"];
        }
    });
    if (par["global"][0] == "dist") {
        std::cout << jian::geom::distance(vec[0], vec[1]) << std::endl;
    } else if (par["global"][0] == "ang") {
        std::cout << jian::geom::angle(vec[0], vec[1], vec[2]) << std::endl;
    } else if (par["global"][0] == "dih") {
        std::cout << jian::geom::dihedral(vec[0], vec[1], vec[2], vec[3]) << std::endl;
    } else if (par["global"][0] == "chir") {
        std::cout << jian::geom::chirality(vec[0], vec[1], vec[2], vec[3]) << std::endl;
    } else {
    }
    return 0;
}











