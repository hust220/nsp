#ifndef JIAN_NUC3D_PARSEHELIX
#define JIAN_NUC3D_PARSEHELIX

#include "../geom.h"
#include "../pdb/Molecule.h"

namespace jian {
namespace nuc3d {

class ParseHelix {
public:
    using val_t = double;
    using Vec = Matrix<val_t, 3, 1>;
    using Mat = Matrix<val_t, -1, -1>;
    using Par = struct {val_t theta, phi, h, d, R;};
    using Result = struct {Vec origin, x, y, z; val_t theta, phi;};

    Par _par {0.562, 1.57, 2.84, -4, 9.7};
    std::map<int, Mat> _cache;

    template<typename Helix>
    auto operator ()(const Helix &helix) {
        return parse(helix);
    }

    auto make_standard_helix(int n) {
        if (_cache.count(n)) return _cache[n];
        else {
            Mat helix = Mat::Zero(n * 2 + 2, 3);
            for (int i = 0; i < n + 1; i++) {
                helix(i, 0) = _par.R*std::cos(i*_par.theta);
                helix(i, 1) = _par.R*std::sin(i*_par.theta);
                helix(i, 2) = i*_par.h;
                helix(2*n+1-i, 0) = _par.R*std::cos(i*_par.theta+_par.phi);
                helix(2*n+1-i, 1) = _par.R*std::sin(i*_par.theta+_par.phi);
                helix(2*n+1-i, 2) = i*_par.h+_par.d;
            }
            return helix;
        }
    }

    template<typename Helix>
    auto get_mat_helix(const Helix &helix) {
        int len = pdb::num_residues(helix); Mat mat(len, 3); 
        int index = 0; for (auto &chain : helix) for (auto &res : chain) {
            auto atm = atom(res, "C4*");
            for (int i = 0; i < 3; i++) mat(index, i) = atm[i];
            index++;
        }
        return mat;
    }

    template<typename Helix>
    auto parse(const Helix &helix) {
        auto mat = get_mat_helix(helix); int len = mat.rows();
        auto sp = geom::suppos(make_standard_helix(len/2-1), mat);
        Vec origin, x, y, z; origin << 0, 0, 0; x << 1, 0, 0; y << 0, 1, 0; z << 0, 0, 1;
        geom::translate(origin, -(sp.c1)); geom::rotate(origin, sp.rot); geom::translate(origin, sp.c2);
        geom::translate(x, -(sp.c1)); geom::rotate(x, sp.rot); geom::translate(x, sp.c2);
        geom::translate(y, -(sp.c1)); geom::rotate(y, sp.rot); geom::translate(y, sp.c2);
        geom::translate(z, -(sp.c1)); geom::rotate(z, sp.rot); geom::translate(z, sp.c2);
        val_t theta = geom::angle(std::vector<val_t>{0, 0, 1}, Vec::Zero(), z);
        val_t phi = geom::angle(std::vector<val_t>{1, 0, 0}, Vec::Zero(), std::vector<val_t>{z[0], z[1], 0});
        return Result{origin, x-origin, y-origin, z-origin, theta, phi};
    }

};

} // namespace nuc3d
} // namespace jian


#endif

