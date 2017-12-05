#include <list>
#include <numeric>
#include <array>
#include <iostream>
#include "nsp.hpp"
#include "cmd_en6p.hpp"
#include "pdb.hpp"
#include "geom.hpp"
#include "string.hpp"
#include "cluster.hpp"
#include "cg.hpp"

BEGIN_JN

    namespace {

        class TrajComponent {
        public:
            Str m_func;
            Str m_traj;
            Str m_ref_file;
            Str m_align;
            Chain m_ref;
            std::deque<int> m_common_ref, m_common_tgt;
            int m_bin = 1;
            bool m_loose;
            Mat m_mat_ref;
            Mat m_mat_tgt;
            Par m_par;

            TrajComponent(const Par &par) {
                m_par = par;

                auto v = par.getv("global");
                assert(v.size() >= 3);
                m_func = v[1];
                to_lower(m_func);

                m_traj = v[2];
                if (v.size() == 4) {
                    m_ref_file = v[3];
                    chain_read_model(m_ref, m_ref_file);
                }

                par.set(m_bin, "b", "bin");

                m_loose = par.has("loose");

                par.set(m_align, "align");
                set_common();

            }

            template<typename _Chain>
            Num rmsd(const _Chain &tgt) {
                std::deque<std::array<double, 3>> l1, l2;
                int i, l;

                auto it_ref = m_common_ref.begin();
                auto it_tgt = m_common_tgt.begin();
                while (it_ref != m_common_ref.end() && it_tgt != m_common_tgt.end()) {
                    const Residue & r1 = m_ref[*it_ref];
                    const Residue & r2 = tgt[*it_tgt];
                    if (!m_loose) assert(r1.name == r2.name);
                    for (auto && a1 : r1) {
                        for (auto && a2 : r2) {
                            if (a1.name == a2.name) {
                                l1.push_back({ a1[0],a1[1],a1[2] });
                                l2.push_back({ a2[0],a2[1],a2[2] });
                            }
                        }
                    }
                    it_ref++;
                    it_tgt++;
                }
                

                i = 0;
                l = l1.size();
                Mat mat1(l, 3);
                Mat mat2(l, 3);
                for (i = 0; i < l; i++) {
                    for (int k = 0; k < 3; k++) {
                        mat1(i, k) = l1[i][k];
                        mat2(i, k) = l2[i][k];
                    }
                }
                return geom::rmsd(mat1, mat2);
            }

            void set_common() {
                std::ifstream ifile;
                Str seq1, seq2;
                std::deque<int> dq1, dq2;
                int i, l, n1, n2;

                if (m_align.empty()) {
                    l = size(m_ref);
                    m_common_ref.resize(l);
                    m_common_tgt.resize(l);
                    std::iota(m_common_ref.begin(), m_common_ref.end(), 0);
                    std::iota(m_common_tgt.begin(), m_common_tgt.end(), 0);
                }
                else {
                    FOPEN(ifile, m_align);
                    ifile >> seq1 >> seq2;
                    assert(size(seq1) == size(seq2));
                    l = size(seq1);
                    n1 = 0;
                    n2 = 0;
                    for (i = 0; i < l; i++) {
                        if (seq1[i] != '-' && seq2[i] != '-') {
                            dq1.push_back(n1);
                            dq2.push_back(n2);
                        }
                        if (seq1[i] != '-') n1++;
                        if (seq2[i] != '-') n2++;
                    }
                    FCLOSE(ifile);
                    m_common_ref = std::move(dq1);
                    m_common_tgt = std::move(dq2);
                }
            }

            void rmsd() {
                for_each_model(m_traj, [this](const Model &model, int i) {
                    if (i % m_bin == 0) {
                        JN_OUT << i+1 << ' ' << rmsd(model.residues()) << std::endl;
                    }
                });
            }

            void compress() {
                MolWriter writer(JN_OUT);
                for_each_model(m_traj, [this, &writer](const Model &model, int i) {
                    if (i % m_bin == 0) {
                        LOG << "Reading: model-" << i + 1 << std::endl;
                        writer.write(model);
                    }
                });
            }

            void extract() {
                Num rmsd_target = -1;
                Num num_target = -1;
                m_par.set(rmsd_target, "rmsd");
                m_par.set(num_target, "n");
                MolReader mol_reader(m_traj);
                auto it = mol_reader.model_begin();
                auto it_target = it;
                Num min = 999;
                for (; it != mol_reader.model_end(); it++) {
                    if (num_target != -1) {
                        if (it.n + 1 == num_target) {
                            it_target = it;
                            JN_OUT << *it;
                            break;
                        }
                    }
                    else {
                        Num d = STD_ fabs(rmsd_target - rmsd(it->residues()));
                        if (d < min) {
                            it_target = it;
                            min = d;
                        }
                    }
                }
                JN_OUT << *it_target;
            }

            void split() {
                Str prefix = file::name(m_traj);
                m_par.set(prefix, "p", "prefix");
                Bool aa = m_par.has("aa");
                for_each_model(m_traj, [this, &prefix, &aa](const Model &model, int i) {
                    if (i % m_bin == 0) {
                        LOG << "Writing: model-" << i + 1 << std::endl;
                        SP<CG> cg = CG::fac_t::make("6p");
                        if (aa) {
                            mol_write(cg->to_aa(model), to_str(prefix, '.', i + 1, ".pdb"));
                        }
                        else {
                            mol_write(model, to_str(prefix, '.', i + 1, ".pdb"));
                        }
                    }
                });
            }

            template<typename _Chain>
            Mat *chain_to_mat(const _Chain &rs) {
                Int l = size(rs);
                Mat *mat = new Mat(l * 6, 3);
                for (Int i = 0; i < l; i++) {
                    for (Int j = 0; j < 6; j++) {
                        for (Int k = 0; k < 3; k++) {
                            (*mat)(i * 6 + j, k) = rs[i][j][k];
                        }
                    }
                }
                return mat;
            }

            template<typename T>
            void write_model(const T &m, Str filename, Bool aa) {
                static SP<CG> cg = CG::fac_t::make("6p");
                if (aa) {
                    mol_write(cg->to_aa(m), filename);
                }
                else {
                    mol_write(m, filename);
                }
            }

            void cluster() {
                Str out_dir = ".";
                m_par.set(out_dir, "out_dir");

                Str prefix = file::name(m_traj);
                m_par.set(prefix, "p", "prefix");

                Int k = 5;
                m_par.set(k, "k");

                Bool aa = m_par.has("aa");

                auto cluster = Cluster::fac_t::make("kmeans", Par("k", k));
                Deque<Mat *> mats;
                Map<Mat *, Num> scores;
                Map<Mat *, Int> inds;
                //Map<Mat *, Model> models;
                //auto dist = [](Mat *m1, Mat *m2) {return geom::rmsd(*m1, *m2); };

                for_each_model(m_traj, [this, &mats, &scores, &inds](const Model &model, int i) {
                    if (i % m_bin == 0) {
                        LOG << "Reading: model-" << i + 1 << std::endl;
                        auto rs = model.residues();
                        Mat *m = chain_to_mat(rs);
                        mats.push_back(m);
                        scores[m] = total(en6p_chain(rs));
                        inds[m] = i;
                    }
                });

                if (size(mats) < k) {
                   Int n = 0;
                   for_each_model(m_traj, [this, &aa, &mats, &k, &out_dir, &prefix, &n](const Model &model, int i) {
                      if (i % m_bin == 0) {
                          write_model(model, to_str(out_dir, '/', prefix, '.', n+1, ".center.pdb"), aa);
                          write_model(model, to_str(out_dir, '/', prefix, '.', n+1, ".lowest.pdb"), aa);
                          n++;
                          if (n == size(mats)) {
                            for (Int j = n; j < k; j++) {
                               write_model(model, to_str(out_dir, '/', prefix, '.', j+1, ".center.pdb"), aa);
                               write_model(model, to_str(out_dir, '/', prefix, '.', j+1, ".lowest.pdb"), aa);
                            }
                          }
                      }
                   });
                }
                else {

                   // fraction of lowest scores
                   Num fl = 0.3;
                   m_par.set(fl, "fl");

                   LOG << "Clustering..." << std::endl;
                   std::sort(mats.begin(), mats.end(), [&scores](auto &&m1, auto &&m2){return scores[m1] < scores[m2];});
                   Mat *mat = Cluster::to_mat(mats.begin(), std::next(mats.begin(), std::max(Int(size(mats)*fl), k)), [](Mat *m1, Mat *m2) {return geom::rmsd(*m1, *m2); });
                   (*cluster)(*mat);

                   Map<Int, Str> names_center, names_lowest;

                   Int n = 0;
                   for (auto && c : cluster->m_clusters) {

                       Int lowest_ind = std::distance(c.begin(), std::min_element(c.begin(), c.end(), [&scores, &mats](Int n1, Int n2){return scores[mats[n1]]<scores[mats[n2]];}));
                       Num center_score = scores[mats[c[0]]];
                       Num lowest_score = scores[mats[c[lowest_ind]]];

                       LOG << "Writing Cluster " << n+1 << "..." << std::endl;
                       JN_OUT << "Cluster " << n+1 << ", size " << size(c) << std::endl;
                       for (auto && ind : c) JN_OUT << ind << ' '; JN_OUT << std::endl;
                       JN_OUT << "center score: " << center_score << ", lowest score: " << lowest_score << std::endl;
                       JN_OUT << std::endl;

                       names_center[inds[mats[c[0]]]] = to_str(out_dir, '/', prefix, '.', n+1, ".center.pdb");
                       names_lowest[inds[mats[c[lowest_ind]]]] = to_str(out_dir, '/', prefix, '.', n+1, ".lowest.pdb");
                       n++;
                   }

                   for_each_model(m_traj, [this, &aa, &names_center, &names_lowest](const Model &model, int i) {
                       if (names_center.find(i) != names_center.end()) {
                           write_model(model, names_center[i], aa);
                       }
                       if (names_lowest.find(i) != names_lowest.end()) {
                           write_model(model, names_lowest[i], aa);
                       }
                   });

                   delete mat;
                }

                for (auto && m : mats) delete m;
            }

            void run() {
                if (m_func == "rmsd") {
                    rmsd();
                }
                else if (m_func == "compress") {
                    compress();
                }
                else if (m_func == "extract") {
                    extract();
                }
                else if (m_func == "split") {
                    split();
                }
                else if (m_func == "cluster") {
                    cluster();
                }
            }
        };

        REGISTER_NSP_COMPONENT(traj) {
            TrajComponent traj(par);
            traj.run();
        }
    }

END_JN

