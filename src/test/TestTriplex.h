#ifndef JIAN_TEST_TESTTRIPLEX
#define JIAN_TEST_TESTTRIPLEX

#include "../pdb/Model.h"
#include "../geom/geometry.h"

namespace jian {
namespace test {

class TestTriplex {
public:
    void operator ()(std::string pdb_name) {
        test_triplex(Model(pdb_name));
    }

    void test_triplex(const Model &model) {
        std::deque<Residue> deque;
        model.each_res([&](const Residue &residue, int res_num) {
            deque.push_back(residue);
            if (deque.size() == 3) {
                std::cout << geometry::distance(deque[0]["C4*"], deque[2]["C4*"]) << std::endl;
                deque.pop_front();
            }
            res_num++;
        });
    }

    void dist(const Model &model, const std::vector<int> &vec) {
        std::vector<Atom> atoms(2);
        model.each_res([&atoms, &vec](const Residue &residue, int res_num) {
            if (res_num == vec[0]) atoms[0] = residue["C4*"];
            if (res_num == vec[1]) atoms[1] = residue["C4*"];
        });
        std::cout << geometry::distance(atoms[0], atoms[1]) << std::endl;
    }

    void chir(const Model &model, const std::vector<int> &vec) {
        std::deque<Atom> deque(4);
        model.each_res([&](const Residue &residue, int res_num) {
            auto result = std::find(vec.begin(), vec.end(), res_num);
            if (result != vec.end()) {
                deque[std::distance(vec.begin(), result)] = residue["C4*"];
            }
        });
        std::cout << geometry::chirality(deque[0], deque[1], deque[2], deque[3]) << std::endl;
    }
};




} // namespace test
} // namespace jian

#endif

