#include "McsmBase.hpp"

#define JN_MCXP_TEMPPAR_SET(a) temp_par.set(PP_CAT(_mc_, a), PP_STRING3(PP_CAT(mc_, a)));
#define JN_MCXP_PAR_SET(a) par.set(PP_CAT(_mc_, a), PP_STRING3(PP_CAT(mc_, a)));
#define JN_MCXP_TEMP(a) LOG << PP_STRING3(PP_CAT(mc_, a)) << ' ' << PP_CAT(_mc_, a) << std::endl;

namespace jian {
	namespace nuc3d {
		namespace mc {

			void MCBase::init(const Par &par) {
				TSP::init(par);

				LOG << "# Set the file of trajectory..." << std::endl;
				par.set(m_traj, "traj");

				LOG << "# Set continuous points..." << std::endl;
				set_continuous_pts();

				LOG << "# Print continuous, angel, dihedral points..." << std::endl;
				for (auto && i : m_continuous_pts) LOG << i << ' '; LOG << std::endl;
				for (auto && i : m_ang_pts) LOG << i << ' '; LOG << std::endl;
				for (auto && i : m_dih_pts) LOG << i << ' '; LOG << std::endl;

				LOG << "# Read initial structure" << std::endl;
				if (par.has("pdb")) chain_read_model(_pred_chain, par.get("pdb"));

				LOG << "# Check constraints" << std::endl;
				validate_constraints();

				//_mc_queue = "heat:30000:20+cool:1000000";
				//par.set(_mc_queue, "mc_queue");
			}

			void MCBase::validate_constraints() {
				int l = _seq.size();
				int i, j;
				for (auto && ct : _constraints.distances) {
					i = ct.key[0];
					j = ct.key[1];
					if (i < 0 || i >= l || j < 0 || j >= l) {
						std::cerr << "Illegal constraints pair (" << i << ", " << j << ") in your constraints file!" << std::endl;
						std::exit(1);
					}
				}
			}

			void MCBase::read_parameters() {
				LOG << "Reading parameters of " << m_par_file << "..." << std::endl;
				Par temp_par(Env::lib() + "/RNA/pars/nuc3d/mc/" + m_par_file + ".par");
				JN_MAP(JN_MCXP_TEMPPAR_SET, JN_MCXP_PARS1, JN_MCXP_PARS2)
			}

			void MCBase::set_parameters(const Par &par) {
				JN_MAP(JN_MCXP_PAR_SET, JN_MCXP_PARS1, JN_MCXP_PARS2)
			}

			void MCBase::print_parameters() {
				JN_MAP(JN_MCXP_TEMP, JN_MCXP_PARS1, JN_MCXP_PARS2)
			}

			void MCBase::set_continuous_pts() {
				int i, n;

				n = 0;
				for (i = 0; i < m_chain_lens.size(); i++) {
					n += m_chain_lens[i];
					m_brk_pts.push_back(n - 1);
				}
				for (int i = 0; i < _seq.size() - 1; i++) {
					if (std::find(m_brk_pts.begin(), m_brk_pts.end(), i) == m_brk_pts.end()) {
						m_continuous_pts.push_back(i);
					}
				}
				for (int i = 0; i < _seq.size() - 2; i++) {
					if (std::find(m_brk_pts.begin(), m_brk_pts.end(), i) == m_brk_pts.end() &&
						std::find(m_brk_pts.begin(), m_brk_pts.end(), i + 1) == m_brk_pts.end()) {
						m_ang_pts.push_back(i);
					}
				}
				for (int i = 0; i < _seq.size() - 3; i++) {
					if (std::find(m_brk_pts.begin(), m_brk_pts.end(), i) == m_brk_pts.end() &&
						std::find(m_brk_pts.begin(), m_brk_pts.end(), i + 1) == m_brk_pts.end() &&
						std::find(m_brk_pts.begin(), m_brk_pts.end(), i + 2) == m_brk_pts.end()) {
						m_dih_pts.push_back(i);
					}
				}
			}

			void MCBase::mc_write() {
				write_traj();
				write_en();
				mc_num_writing++;
			}

			void MCBase::write_traj() {
				if (!(m_traj.empty())) {
					if (mc_num_writing == 1) {
						file::clean(m_traj);
					}
					std::ofstream output(m_traj.c_str(), std::ios::app);
					m_writer.bind_stream(output);
					m_writer.write_model([&]() {
						this->m_writer.write(this->_pred_chain);
					});
					output.close();
				}
			}

			void MCBase::mc_sample() {
				mc_sample_res();
			}

			void MCBase::mc_sample_res() {
				backup();
				if (rand() < 0.5) {
					// translate
					int index = int(rand() * 3);
					double dist = (rand() - 0.5) * 2 * _mc_max_shift;
					for (int i = 0; i < _seq.size(); i++) {
						if (is_selected(i)) {
							for (auto && atom : _pred_chain[i]) {
								atom[index] += dist;
							}
							space_update_item(i);
						}
					}
				}
				else {
					// rotate
					int index = int(rand() * 3);
					double dih = (rand() - 0.5) * PI / 6;
					auto &&rot = geom::rot_mat(index, dih);
					auto &&origin = rotating_center();
					for (int i = 0; i < _seq.size(); i++) {
						if (is_selected(i)) {
							for (auto && atom : _pred_chain[i]) {
								geom::rotate(atom, origin, rot);
							}
							space_update_item(i);
						}
					}
				}
			}

			void MCBase::mc_back() {
				for (int i = 0; i < _seq.size(); i++) {
					if (is_selected(i)) {
						for (auto && atom : _pred_chain[i]) {
							atom = _moved_atoms.front();
							_moved_atoms.pop_front();
						}
						space_update_item(i);
					}
				}
			}

			void MCBase::backup() {
				_moved_atoms.clear();
				for (int i = 0; i < _seq.size(); i++) {
					if (is_selected(i)) {
						for (auto && atom : _pred_chain[i]) {
							_moved_atoms.push_back(atom);
						}
					}
				}
			}

			void MCBase::init_space() {
				m_item_space.resize(_seq.size());
				for (int i = 0; i < _seq.size(); i++) {
					space_val_t &s = space_val(i);
					s.push_back(i);
					m_item_space[i] = &s;
				}
			}

			int MCBase::space_index(double n) const {
				return int((n + 1000) / m_box_size);
			}

			MCBase::item_t &MCBase::item(int i) {
				return _pred_chain[i][0];
			}

			MCBase::space_val_t &MCBase::space_val(int i) {
				item_t &a = item(i);
				return m_space[space_index(a[0])][space_index(a[1])][space_index(a[2])];
			}

			void MCBase::space_update_item(int i) {
				space_val_t &n = space_val(i);
				space_val_t &o = *(m_item_space[i]);
				if (&o != &n) {
					o.erase(std::find(o.begin(), o.end(), i));
					m_item_space[i] = &n;
					n.push_back(i);
				}
			}

			void MCBase::run() {
				LOG << "# Check initial structure..." << std::endl;
				if (num_residues(_pred_chain) == 0) {
					throw "Please give an initial structure before the optimization procedure!";
				}

				LOG << "# Read parameters..." << std::endl;
				m_par_file = m_cg_type;
				_par->set(m_par_file, "par_file");
				read_parameters();

				LOG << "# Set parameters..." << std::endl;
				set_parameters(*_par);

				LOG << "# Print parameters..." << std::endl;
				print_parameters();

				LOG << "# Initializing running..." << std::endl;
				before_run();

				LOG << "# Save helices" << std::endl;
				save_fixed_ranges();

				LOG << "# Carrying on CG processing with the Chain..." << std::endl;
				_pred_chain = m_cg->to_cg(_pred_chain);

				LOG << "# Init space..." << std::endl;
				init_space();

				LOG << "# MC..." << std::endl;
				mc_run();

				LOG << "# Print Constraints and Distances..." << std::endl;
				print_final_constraints();

				LOG << "# Finishing running..." << std::endl;
				finish_run();

				LOG << "# Coarsed Grained To All Atom..." << std::endl;
				cg_to_aa(chain_to_coords());

				LOG << "# Restore helix..." << std::endl;
				restore_fixed_ranges();

				LOG << "# Transform..." << std::endl;
				this->transform();

			}

			void MCBase::transform() {
				Model m;
				m.push_back(_pred_chain);
				_pred_chain = std::move(jian::transform(m, _seq, _type)[0]);
			}

			Mat MCBase::chain_to_coords() {
				int n_atoms = num_atoms(_pred_chain);
				Mat c(n_atoms, 3);
				int n_atom = 0;
				for (int i = 0; i < _pred_chain.size(); i++) {
					for (auto && atom : _pred_chain[i]) {
						for (int j = 0; j < 3; j++) {
							c(n_atom, j) = atom[j];
						}
						n_atom++;
					}
				}
				return c;
			}

			void MCBase::print_final_constraints() {
				double d;
				int i, j;
				LOG << "Print Contacts:" << std::endl;
				for (auto && c : _constraints.contacts) {
					i = c.key[0];
					j = c.key[1];
					d = dist_two_res(_pred_chain[i], _pred_chain[j]);
					LOG << i << ' ' << j << " weight:" << c.weight << " dist:" << d << std::endl;
				}
				LOG << "Print Distances:" << std::endl;
				for (auto && c : _constraints.distances) {
					i = c.key[0];
					j = c.key[1];
					d = dist_two_res(_pred_chain[i], _pred_chain[j]);
					LOG << i << ' ' << j << " value:" << c.value << " weight:" << c.weight << " dist:" << d << std::endl;
				}
			}

			void MCBase::cg_to_aa(const Mat &c) {
				_pred_chain = m_cg->to_aa(c, 0, c.rows() - 1);
			}

			void MCBase::before_run() {}
			void MCBase::finish_run() {}

			std::string MCBase::file_parameters() const {
				return "3drna";
			}

			void MCBase::save_fixed_ranges() {}
			void MCBase::restore_fixed_ranges() {}

		}
	}
}
