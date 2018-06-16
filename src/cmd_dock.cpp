#include "nsp.hpp"
#include "pdb.hpp"
#include "geom.hpp"
#include "matrix.hpp"
#include "fftw3.h"
//#include "simple_fft/fft.h"

namespace jian {

struct Dock {
    // parameters
    Int n_max = 9999;
    Num box;
    Num bin = 0.7;
    Num bin_angle = PI/6;
    Int bins;
    Str prefix = "dock";
    Num interior_rec = -15;
    Num interior_lig = 1;
    Num surface_thickness = 1;
    Num atom_radius = 1;
    Bool b_rotate = true;

    // temporary variables
    Str name1;
    Str name2;
    Chain chain1;
    Chain chain2;
    Chain *rec; // receptor
    Chain *lig; // ligand
    Vec center_r;
    Vec center_l;
    Num rad_r;
    Num rad_l;
    Num pos_r;
    Num pos_l;
    Mat3c grid_r1;
    Mat3c grid_r2;
    Mat3c grid_l1;
    Mat3c grid_l2;
    Mat3c grid_c1;
    Mat3c grid_c2;
    Array<Int, 3> max_pos;

    Vec chain_center(const Chain & rs) {
        Vec v = Vec::Zero(3);
        Int n = 0;
        for (auto && r : rs) {
            for (auto && a : r) {
                for (Int i = 0; i < 3; i++) {
                    v[i] += a[i];
                }
                n++;
            }
        }
        for (Int i = 0; i < 3; i++) {
            v[i] /= n;
        }
        return v;
    }

    template<typename T>
        void chain_translate(Chain& chain, T && v) {
            for (auto && res : chain) {
                for (auto && atom : res) {
                    for (Int i = 0; i < 3; i++) {
                        atom[i] += v[i];
                    }
                }
            }
        }

    Num chain_max_dist(const Chain & rs, const Vec &v) {
        Num max = 0;
        for (auto && r : rs) {
            for (auto && a : r) {
                Num d = geom::distance(a, v);
                if (max < d) max = d;
            }
        }
        return max;
    }

    void set_grid(const Chain &chain, Mat3c &mat, Num v) {
        Int n = std::ceil((atom_radius+surface_thickness)/bin);
        mat.set_all(0);
        for (auto && res : chain) {
            for (auto && atom : res) {
                Array<Int, 3> pos;
                for (Int i = 0; i < 3; i++) {
                    pos[i] = Int(atom[i]/bin);
                }
                for (Int a = -n; a <= n; a++) for (Int b = -n; b <= n; b++) for (Int c = -n; c <= n; c++) {
                    Int i = pos[0]+a;
                    Int j = pos[1]+b;
                    Int k = pos[2]+c;
                    if (i>=0 && i < bins && j >=0 && j < bins && k >=0 && k < bins) {
                        Num d = geom::distance(Vector<Num>{i*bin, j*bin, k*bin}, atom);
                        if (d <= atom_radius) {
                            mat(i, j, k) = v;
                        } else if (mat(i, j, k) == 0.0 && d <= atom_radius+surface_thickness) {
                            mat(i, j, k) = 1;
                        }
                    }
                }
            }
        }
    }

    void print_surface(const Mat3 &mat) {
        for (Int i = 0; i < bins; i++) {
            for (Int j = 0; j < bins; j++) {
                for (Int k = 0; k < bins; k++) {
                    if (mat(i, j, k) == 1) {
                        Out << i*bin << ' ' << j*bin << ' ' << k*bin << Endl;
                    }
                }
            }
        }
    }

    void write_complex(Str filename) {
        Model m;
        m.push_back(*rec);
        m.push_back(*lig);
        m[0].name = "A";
        m[1].name = "B";
        mol_write(m, filename);
    }

    void cycle(Int & v, Num p) {
        Num d = p-v*bin;
        if (d < 0) {
            v-=bins;
            cycle(v, p);
        }
        else if (d > bins*bin) {
            v+=bins;
            cycle(v, p);
        }
    }

    void cal_ligand(Str filename) {
        Out << "[*] Set the grid of ligand..." << Endl;
        set_grid(*lig, grid_l1, interior_lig);
        Out << "[*] FFT of ligand..." << Endl;
        //grid_l1 = grid_l;
        grid_l1.fft(grid_l2);

        for (Int i = 0; i < bins; i++) for (Int j = 0; j < bins; j++) for (Int k = 0; k < bins; k++) {
            grid_c2(i, j, k) = std::conj(grid_r2(i, j, k)) * grid_l2(i, j, k);
        }

        Out << "[*] IFFT of correlation function..." << Endl;
        grid_c2.ifft(grid_c1);
        grid_c1.divide(bins*bins*bins);

        Out << "[*] Find pick..." << Endl;
        Num max = 0;
        for (Int i = 0; i < bins; i++) for (Int j = 0; j < bins; j++) for (Int k = 0; k < bins; k++) {
            Num d = grid_c1(i, j, k).real();
            if (max < d) {
                max = d; max_pos[0] = i; max_pos[1] = j; max_pos[2] = k;
            }
        }
        Out << "[*] Pick: " << max_pos[0] << ' ' << max_pos[1] << ' ' << max_pos[2] << " " << max << Endl;

        for (Int i = 0; i < 3; i++) cycle(max_pos[i], pos_l);
        chain_translate(*lig, Vector<Num>{-max_pos[0]*bin, -max_pos[1]*bin, -max_pos[2]*bin});
        Out << "[*] Writing to '" << filename << "'" << Endl;
        write_complex(filename);
        chain_translate(*lig, Vector<Num>{max_pos[0]*bin, max_pos[1]*bin, max_pos[2]*bin});

    }

    void run() {
        Int l1 = num_atoms(chain1);
        Int l2 = num_atoms(chain2);
        if (l1 > l2) {
            Out << "Receptor:" << name1 << " Ligand:" << name2 << Endl;
            rec = &chain1;
            lig = &chain2;
        }
        else {
            Out << "Receptor:" << name2 << " Ligand:" << name1 << Endl;
            rec = &chain2;
            lig = &chain1;
        }

        center_r = chain_center(*rec);
        center_l = chain_center(*lig);
        Out << "Center of Receptor: " << center_r[0] << ' ' << center_r[1] << ' ' << center_r[2] << Endl;
        Out << "Center of Ligand: " << center_l[0] << ' ' << center_l[1] << ' ' << center_l[2] << Endl;

        rad_r = chain_max_dist(*rec, center_r);
        rad_l = chain_max_dist(*lig, center_l);
        box = 2*rad_r;
        bins = std::ceil(box / bin);
        pos_r = rad_r;
        pos_l = rad_r;
        Out << "bin: " << bin << ", box: " << box << ", bins: " << bins << Endl;

        Out << "[*] Allocating..." << Endl;
        grid_r1.resize(bins, bins, bins);
        grid_r2.resize(bins, bins, bins);
        grid_l1.resize(bins, bins, bins);
        grid_l2.resize(bins, bins, bins);
        grid_c1.resize(bins, bins, bins);
        grid_c2.resize(bins, bins, bins);

        chain_translate(*rec, Vector<Num>{pos_r-center_r[0], pos_r-center_r[1], pos_r-center_r[2]});
        chain_translate(*lig, Vector<Num>{pos_l-center_l[0], pos_l-center_l[1], pos_l-center_l[2]});
        //write_complex(to_str(prefix, "-0.pdb"));

        Out << "[*] Set the grid of receptor..." << Endl;
        set_grid(*rec, grid_r1, interior_rec);
        Out << "[*] FFT of receptor..." << Endl;
        //grid_r1 = grid_r;
        grid_r1.fft(grid_r2);

        if (b_rotate) {
            Int n = 0;
            for (Num dx = 0; dx < PI && n < n_max; dx += bin_angle) {
                for (Num dy = 0; dy < 2*PI && n < n_max; dy += bin_angle) {
                    for (Num dz = 0; dz < 2*PI && n < n_max; dz += bin_angle) {
                        Out << "[*] ========================" << Endl;
                        Out << "[*] Step " << n+1 << " dx:" << dx << " dy:" << dy << " dz:" << dz << Endl;
                        Mat mat = geom::rot_mat(0, dx) * geom::rot_mat(1, dy) * geom::rot_mat(2, dz);
                        for (auto && res : *lig) for (auto && atom : res) {
                            for (Int i = 0; i < 3; i++) atom[i] -= pos_l;
                            geom::rotate(atom, mat);
                            for (Int i = 0; i < 3; i++) atom[i] += pos_l;
                        }

                        cal_ligand(to_str(prefix, '-', n+1, ".pdb"));

                        n++;
                    }
                }
            }
        }
        else {
            cal_ligand(to_str(prefix, ".pdb"));
        }

    }

    void operator ()(const Par &par) {
        auto g = par.getv("global");
        par.set(bin, "b", "bin");
        par.set(bin_angle, "ba", "bin_angle");
        par.set(prefix, "p", "prefix");
        par.set(n_max, "nm", "n_max");
        par.set(atom_radius, "ar", "atom_radius");
        par.set(surface_thickness, "st", "surface_thickness");
        b_rotate = !par.has("nr", "no_rotate");
        name1 = g[1];
        name2 = g[2];
        chain_read_model(chain1, name1);
        chain_read_model(chain2, name2);
        run();
    }
};

REGISTER_NSP_COMPONENT(dock) {
    Dock dock;
    dock(par);
}

}

