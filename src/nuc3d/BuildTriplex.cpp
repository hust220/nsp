#include "BuildTriplex.h"
#include "BuildHelix.h"

namespace jian {
namespace nuc3d {

BuildTriplex::BuildTriplex() {
}

void BuildTriplex::operator ()(std::string seq, nuc2d::Seq2Tri::TupleInfo info, int num) {
    /// Read sequence
    _seq = seq;

    /// Transform information to ct
    _ct = info_to_ct(info);

    /// Calculate scaffold
    get_scaffold();

}

std::vector<std::vector<int>> BuildTriplex::info_to_ct(const nuc2d::Seq2Tri::TupleInfo &info) {
    std::vector<std::vector<int>> ct;
    for (int i = 0; i < info.first.size(); i++) {
        if (!info.first[i].empty()) {
            if (std::count_if(ct.begin(), ct.end(), [&](const std::vector<int> &v){
                    return std::count(v.begin(), v.end(), i);}) == 0) {
                std::vector<int> vec(info.first[i].size() + 1);
                std::copy(info.first[i].begin(), info.first[i].end(), std::next(vec.begin(), 1));
                vec[0] = i;
                ct.push_back(std::move(vec));
            }
        }
    }
    return ct;
}

void BuildTriplex::get_scaffold() {
    /// Get length of bound matrix.
    int len = std::accumulate(_ct.begin(), _ct.end(), 0, [&](int i, const std::vector<int> &vec){
            return i + vec.size();});

    /// Initialize anchor list.
    std::vector<int> anchors(len);
    int index = 0;
    for (auto &&vec: _ct) {
        for (auto &&i: vec) {
            anchors[index] = i;
            index++;
        }
    }
    std::sort(anchors.begin(), anchors.end(), std::less<int>());
    std::map<int, int> anchor_indices;
    for (int i = 0; i < len; i++) anchor_indices[anchors[i]] = i;

    /// Initialize bound matrix.
    _bound = MatrixXf::Zero(len, len);
    for (int i = 0; i < len; i++) {
        for (int j = i; j < len; j++) {
            if (i == j) {
                _bound(i, j) = 0;
            } else {
                _bound(i, j) = 6.1 * (anchors[j] - anchors[i]);
                _bound(j, i) = 6.1;
            }
        }
    }

    /// Add helix parameters into the bound matrix
    for (auto &&vec: _ct) {
        if (vec.size() == 2) {
            int a = anchor_indices[vec[0]];
            int b = anchor_indices[vec[1]];
            _bound(a, b) = _bound(b, a) = 15.1;
        } else if (vec.size() == 3) {
            int a = anchor_indices[vec[0]];
            int b = anchor_indices[vec[1]];
            int c = anchor_indices[vec[2]];
            _bound(a, b) = _bound(b, a) = 15.1;
            _bound(b, c) = _bound(c, b) = 15.1;
        } else if (vec.size() == 4) {
            int a = anchor_indices[vec[0]];
            int b = anchor_indices[vec[1]];
            int c = anchor_indices[vec[2]];
            int d = anchor_indices[vec[3]];
            _bound(a, b) = _bound(b, a) = 15.1;
            _bound(b, c) = _bound(c, b) = 15.1;
            _bound(c, d) = _bound(d, c) = 15.1;
            _bound(a, d) = _bound(d, a) = 15.1;
        } else {
            die("jian::nuc3d::BuildTriplex::get_scaffold error!");
        }
    }

    /// Calculate scaffold coordinates
    DG dg(_bound);
    std::cout << "\nscaffold bound matrix: " << std::endl;
    std::cout << _bound << "\n" << std::endl;
    _scaffold = dg();
    std::cout << "DG...\n";
    std::cout << "scaffold energy: " << dg.E << std::endl;
    std::cout << "scaffold: " << std::endl;
    std::cout << _scaffold << std::endl;

}

} /// namespace nuc3d
} /// namespace jian

