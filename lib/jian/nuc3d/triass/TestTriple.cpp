#include "../../pdb/Model.h"
#include "../../util/Par.h"
#include "../../geom/geometry.h"
#include "../../fpl.h"

int main(int argc, char **argv) {
    jian::Par par(argc, argv);
    jian::Model model(par["model"][0]);
    auto type = par["global"][0];
    if (type == "dist") {
        int len = model.res_nums(); auto mat = jian::mat::make_mat(len, len);
        model.each_res([&](const auto &res1, int n1){
            model.each_res([&](const auto &res2, int n2){
                jian::mat::ref(mat, n1, n2) = jian::geometry::distance(res1["C4*"], res2["C4*"]);
            });
        });
        std::cout << mat << std::endl;
    } else {
        std::array<jian::Atom, 4> arr;
        auto indices = jian::fpl::map<std::vector>([](const auto &s){return std::stoi(s);}, par["num"]);
        auto flag = jian::fpl::let([&]{std::map<int, int> m; for (int i = 0; i < 4; i++) m[indices[i]] = i; return m;});
        model.each_res([&](const auto &res, int n){
            if (jian::fpl::exists([&](auto i){return i == n;}, indices)) {
                arr[flag[n]] = res["C4*"];
            }
        });
        std::cout << jian::geometry::dihedral(arr[0], arr[1], arr[2], arr[3]) << std::endl;
    }
}










