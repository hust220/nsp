#ifndef JIAN_DG_TESTHELIXCHIR
#define JIAN_DG_TESTHELIXCHIR

#include "../pdb.h"
#include "../util.h"
#include "../geom.h"

namespace jian {
namespace dg {

class TestHelixChir {
public:
    Log log;

    template<typename ModelType> void operator ()(ModelType &&model) {
        test_helix_chir(std::forward<ModelType>(model));
    }

    template<typename ModelType> void test_helix_chir(ModelType &&model) {
        std::deque<Point> pts;
        model.each_res([&pts](const Residue &res, int res_num){
            pts.push_back(res["C4*"].pos());
            if (res_num > 2) {
                log(geometry::chirality(pts[0], pts[1], pts[2], pts[3]), '\n');
                pts.pop_front();
            }
        });
    }
};





} // namespace dg
} // namespace jian

#endif

