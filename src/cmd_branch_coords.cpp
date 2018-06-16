#include "nsp.hpp"
#include "rtsp_parse_helix.hpp"
#include "rss.hpp"
#include "rss_get_ss.hpp"
#include "rss_sst.hpp"

namespace jian {

struct BranchCoords {
    Vector<Int> get_helix_nums(SSE *l) {
        Int len = size(l->helix);
        Vector<Int> nums(len * 2);
        Int index = 0;
        for (auto && b : l->helix) {
            nums[index] = b.res1.num - 1;
            nums[2 * len - 1 - index] = b.res2.num - 1;
            index++;
        }
        return nums;
    }

    template<typename T, typename U>
    Mat mat_nums(T && rs, U && nums) {
        Int l = size(nums);
        Mat m(l, 3);
        Int i = 0;
        for (auto && n : nums) {
            auto && atom = rs[n]["C4*"];
            for (Int j = 0; j < 3; j++) m(i, j) = atom[j];
            i++;
        }
        return m;
    }

    void operator ()(const Par &par) {
        auto g = par.getv("global");
        Model m = mol_read_to<Model>(g[1]);
        auto rs = m.residues();
        Str seq = jian::seq(m);
        Str ss = get_ss(rs);
        SST sst = sst_new(seq, ss, 2);
        for_ (sse, tree_nodes(sst.head)) {
            if (sse->has_helix()) {
                auto nums = get_helix_nums(sse);
                auto r = parse_helix(mat_nums(rs, nums));
                JN_OUT << r.theta << ':' << r.phi << ' ';
                for (auto son = sse->son; son != NULL; son = son->bro) {
                    auto nums = get_helix_nums(son);
                    auto r = parse_helix(mat_nums(rs, nums));
                    JN_OUT << r.theta << ':' << r.phi << ' ';
                }
                JN_OUT << std::endl;
            }
        } end_;
        tree_free(sst.head);

    }
};

REGISTER_NSP_COMPONENT(branch_coords) {
    BranchCoords b;
    b(par);
}

}

