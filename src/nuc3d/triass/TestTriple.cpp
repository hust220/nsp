#include "../../pdb/Model.h"
#include "../../util/Par.h"
#include "../../geom/geometry.h"

int main(int argc, char **argv) {
    jian::Par par(argc, argv);
    jian::Model model(par["model"][0]);
    int len = model.res_nums(); auto mat = jian::mat::make_mat(len, len);
    model.each_res([&](const auto &res1, int n1){
        model.each_res([&](const auto &res2, int n2){
            jian::mat::ref(mat, n1, n2) = jian::geometry::distance(res1["C4*"], res2["C4*"]);
        });
    });
    std::cout << mat << std::endl;
}










