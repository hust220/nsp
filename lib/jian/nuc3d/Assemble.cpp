#include <algorithm>
#include "../utils/Par.hpp"
#include "../utils/Env.hpp"
#include "../utils/file.hpp"
#include "../jian.hpp"
#include "Assemble.hpp"
#include "BuildHelix.hpp"

BEGIN_JN
namespace nuc3d {

void chain_read_record(Chain &chain, const record_t &templ_res) {
    STD_ ostringstream stream;
    stream << Env::lib() << "/RNA/templates/" << templ_res._name << ".pdb";
    chain_read_model(chain, stream.str());
}

void find_loop_records(SSE *l, records_t &records, S name, 
	const pdbs_t &used_pdbs, const pdbs_t &disused_pdbs, family_t family)
{
    if (l->loop.empty()) return;

    S seq = l->loop.seq(),
    ss = l->loop.ss(),
    p_ss = NASS::pure_ss(ss),
    lower_ss = NASS::lower_ss(p_ss, 1);

    int num_sons = l->num_sons();

	STD_ ostringstream stream;
    stream << Env::lib() << "/RNA/records/" << (l->is_open() ? "open_" : "") << "loop";
//LOGI << "FIND LOOP RECORDS OF: " << l << ' ' << stream.str() << std::endl;
	STD_ ifstream ifile;
	FOPEN(ifile, stream.str());

    record_t templ_rec;
    S line;
    int num = 0;
    while (std::getline(ifile, line) && set_loop_rec(templ_rec, line)) {
        if (std::find(disused_pdbs.begin(), disused_pdbs.end(), templ_rec._src) != disused_pdbs.end()) {
            continue;
        } else if ((!used_pdbs.empty()) &&
                   std::find(used_pdbs.begin(), used_pdbs.end(), templ_rec._src) == used_pdbs.end()) {
            continue;
        } else if (NASS::pure_ss(NASS::lower_ss(templ_rec._ss, 1)) == lower_ss) {
            templ_rec._score = (templ_rec._ss == ss ? 5 : 0);
            if (templ_rec._src == jian::upper(name)) {
                templ_rec._score += 5;
            }
            if (templ_rec._family == family && family != "other") {
                templ_rec._score += 2;
            }
            for (uint i = 0; i < templ_rec._seq.size(); i++) {
                if (seq[i] == templ_rec._seq[i]) {
                    if (ss[i] == templ_rec._ss[i] && ss[i] != '.' && ss[i] != '(' && ss[i] != ')') {
                        templ_rec._score += 2;
                    } else if (ss[i] == '(' || ss[i] == ')') {
                        templ_rec._score += 0.2;
                    } else {
                        templ_rec._score += 1;
                    }
                }
            }
            records.push_back(templ_rec);
            num++;
        }
    }
    ifile.close();
    std::sort(records.begin(), records.end(), [](auto &&loop1, auto &&loop2) {
        return loop1._score > loop2._score;
    });
}

void find_helix_records(SSE *l, records_t &records, S name, S family) {
    if (l->helix.empty()) return;

	S seq = l->helix.seq();
	S ss = l->helix.ss();

    int len = size(l->helix);

	STD_ ostringstream stream;
    stream << Env::lib() << "/RNA/records/helix";
	STD_ ifstream ifile;
    FOPEN(ifile, stream.str());

    record_t templ_rec;
    S line;
    int num = 0;
    while (std::getline(ifile, line) && set_helix_rec(templ_rec, line)) {
        if (templ_rec._len == len) {
            templ_rec._score = 0;
            if (templ_rec._src == jian::upper(name)) {
                templ_rec._score += 5;
            }
            if (templ_rec._family == family && family != "other") {
                templ_rec._score += 2;
            }
            for (uint i = 0; i < templ_rec._seq.size(); i++) {
                if (seq[i] == templ_rec._seq[i]) {
                    templ_rec._score++;
                }
            }
            records.push_back(std::move(templ_rec));
            num++;
        }
    }
    ifile.close();

    std::sort(records.begin(), records.end(), [](auto &&loop1, auto &&loop2) {
        return loop1._score > loop2._score;
    });
}

    Assemble::Assemble(const Par &par) {
        TSP::init(par);

        log << "# Set disused pdbs..." << STD_ endl;
        par.setv(m_disused_pdbs, "disused_pdbs");
        //for (Str & pdb : m_disused_pdbs) STD_ cout << pdb << STD_ endl;

        m_sample = par.has("sample");
        par.set(m_sample_mode, "sample_mode");

		m_loop_building_method = "partial_dg";
		par.set(m_loop_building_method, "loop_building", "lb");
		to_lower(m_loop_building_method);

        log << "# Check input" << std::endl;
		S info_errors;
		auto b = NASS::check_ss(_ss, info_errors);
        if (!b) throw info_errors;

        log << "# Constructing 2D Structure Tree" << std::endl;
        _ss_tree.make(_seq, _ss, _hinge);

		log << "# Searching records..." << std::endl;
		find_records();

		//log << "# Completing records..." << std::endl;
		//complete_records();

		log << "# Printing records..." << std::endl;
        print_records();

    }

    void Assemble::run() {
        std::ostringstream stream;
        int n;

        predict_one();
        if (m_sample) {
            n = 1;
            log << "# Writing sampling structure " << n << std::endl;
            stream.str("");
            stream << _name << ".assemble." << n << ".pdb";
            mol_write(_pred_chain, stream.str());
            for (n = 2; n <= _num; n++) {
                if (m_sample_mode == "sample_one") {
                    sample_one_template();
                } else if (m_sample_mode == "sample_all") {
                    sample_all_templates();
                } else {
                    throw "Illegal sampling mode!";
                }
                assemble();
                log << "# Writing sampling structure " << n << std::endl;
                stream.str("");
                stream << _name << ".assemble." << n << ".pdb";
                mol_write(_pred_chain, stream.str());
            }
        }
    }

    void Assemble::predict_one() {
        log << "# Select Templates" << std::endl;
        select_templates(); 

        log << "# Print Templates" << std::endl;
        print_templates();

        log << "# Assemble" << std::endl;
        assemble();
    }

    double Assemble::score_templates() {
        double e = 0;
        for (auto && pair : m_records) {
            auto & l = pair.first;
            if (l->has_loop()) {
                e += m_selected_record[l].first._score;
            }
            if (l->has_helix()) {
                e += m_selected_record[l].second._score;
            }
        }
        return e;
    }

    bool Assemble::lack_templates() {
        SSE *l;
        for (auto && pair : m_templates) {
            l = pair.first;
            if (l->has_loop()) {
                if (m_records[l].first.empty()) return true;
            }
        }
        return false;
    }

	void Assemble::print_templates() {
		for (auto &&sse : _ss_tree) {
			log << &sse << " : " << m_templates[&sse].first.model_name << ' ' << m_templates[&sse].second.model_name << std::endl;
		}
	}

    void Assemble::assemble() {
        log << "## Score of Templates: " << score_templates() << std::endl;

        log << "## Print Templates." << std::endl;
        print_templates();

        log << "## Position Templates." << std::endl;
        position_templates();

        log << "## Assemble Templates." << std::endl;
        assemble_templates(&_ss_tree.front());

        log << "## Transform." << std::endl;
        this->transform();
    }

    void Assemble::transform() {
        Model m;
        m.push_back(_pred_chain);
        _pred_chain = std::move(jian::transform(m, _seq, _type)[0]);
    }

	void Assemble::set_loop_template(SSE *l, Bool is_first) {
		if (l != NULL && l->has_loop()) {
			if (m_loop_building_method == "all_dg" ||
				(m_loop_building_method == "partial_dg" && m_records[l].first.empty())) {
				build_loop_dg.init(l->loop.seq(), NASS::lower_ss(l->loop.ss()));
				m_templates[l].first = build_loop_dg();
				m_selected_record[l].first = TemplRec{};
			}
			else if (m_loop_building_method == "all_raw" ||
				(m_loop_building_method == "partial_raw" && m_records[l].first.empty())) {
				build_loop_raw.init(l->loop.seq(), NASS::lower_ss(l->loop.ss()));
				m_templates[l].first = build_loop_raw();
				m_selected_record[l].first = TemplRec{};
			}
			else if (m_loop_building_method == "all_mc" ||
				(m_loop_building_method == "partial_mc" && m_records[l].first.empty())) {
				build_loop_raw.init(l->loop.seq(), NASS::lower_ss(l->loop.ss()));
				sample_loop.init(build_loop_raw(), l->loop.ss());
				m_templates[l].first = sample_loop();
				m_selected_record[l].first = TemplRec{};
			}
			else {
				int n = (is_first ? 0 : int(rand() * m_records[l].first.size()));
				m_templates[l].first = load_pdb(m_records[l].first[n]);
				m_selected_record[l].first = m_records[l].first[n];
			}
		}
	}

	void Assemble::set_helix_template(SSE *sse, Bool is_first) {
		if (sse->has_helix()) {
			auto &records = m_records[sse].second;
			if (records.empty()) {
				chain_append(m_templates[sse].second, build_helix(sse->helix.seq()));
				m_selected_record[sse].second = TemplRec{};
			}
			else {
				m_templates[sse].second = load_pdb(m_records[sse].second[0]);
				m_selected_record[sse].second = m_records[sse].second[0];
			}
		}
	}

    void Assemble::select_templates() {
		log << _ss_tree << STD_ endl;
		for (auto && sse : _ss_tree) {
			log << "# SSE: " << &sse << STD_ endl;
			log << sse << STD_ endl;
			log << "# Select template of loop of SSE: " << &sse << STD_ endl;
			set_loop_template(&sse, true);
			log << "# Select template of helix of SSE: " << &sse << STD_ endl;
			set_helix_template(&sse, true);
		}
    }

	void Assemble::sample_helix_template() {
		assert(m_records.size() == 1);
		records_t & r = m_records.begin()->second.second;
		uint n = static_cast<uint>(rand() * r.size());
		m_templates[0].second = load_pdb(r[n]);
		m_selected_record[0].second = r[n];
	}

    void Assemble::sample_one_template() {
		SSE * l = select_loop();
		if (l == NULL) {
			sample_helix_template();
		}
		else {
			sample_loop_template(select_loop());
		}
    }

    void Assemble::sample_all_templates() {
        for (auto && pair : m_templates) {
            SSE *l = pair.first;
            if (l->has_loop()) {
                sample_loop_template(l);
            }
        }
    }

    void Assemble::sample_loop_template(SSE *l) {
		set_loop_template(l, false);
    }

    SSE *Assemble::select_loop() {
        std::deque<SSE *> ls;
        for (auto && pair : m_records) {
            auto &l = pair.first;
            if (l->has_loop()) {
                ls.push_back(l);
            }
        }
		return (ls.empty() ? NULL : ls[static_cast<uint>(ls.size() * jian::rand())]);
    }

    Chain Assemble::load_pdb(const TemplRec &templ_res, S type) {
        Chain chain;
        chain_read_record(chain, templ_res);
        return chain;
    }

    void Assemble::assemble_templates(SSE *l) {
        _pred_chain.resize(_seq.size());
		for (auto &&sse : _ss_tree) {
			if (sse.has_loop()) {
				auto &temp_residues = m_templates[&sse].first;
				int i = 0;
				for (auto && res : sse.loop) {
					if (res.type != '&') {
						_pred_chain[res.num - 1] = temp_residues[i];
						i++;
					}
				}
			}
			if (sse.has_helix()) {
				int len = size(sse.helix);
				auto &temp_residues = m_templates[&sse].second;
				int n_bp = 0;
				for (auto && bp : sse.helix) {
					_pred_chain[bp.res1.num - 1] = temp_residues[n_bp];
					_pred_chain[bp.res2.num - 1] = temp_residues[2 * len - 1 - n_bp];
					n_bp++;
				}
			}
		}
    }

    void Assemble::position_templates() {
        position_templates(_ss_tree.root(), L<Mat>());
    }

    void Assemble::position_templates(SSTree::El *l, L<Mat> mats) {
        if (l == NULL) return;
        auto &loop = m_templates[&l->data].first;
        auto &helix = m_templates[&l->data].second;
        if (l->data.has_helix()) {
            // position helix
            if (!mats.empty()) {
                auto mat1 = mats.front();
                position_model(helix, mat1);
            }
            if (l->data.has_loop()) {
                // position SSE
                int len = helix.size();
                auto mat2 = model_mat(helix, std::list<int>{len/2-2, len/2-1, len/2, len/2+1});
                position_model(loop, mat2);
                // position son
                if (l->has_son()) {
                    auto son_mats = loop_mats(loop, &l->data);
                    position_templates(l->son, son_mats);
                }
                // position brother
                if (l->has_brother()) {
                    mats.pop_front();
                    position_templates(l->brother, mats);
                }
            } else {
                // position brother
                if (l->has_brother()) {
                    mats.pop_front();
                    position_templates(l->brother, mats);
                }
            }
        } else {
            // position son
            if (l->has_son()) {
                auto son_mats = loop_mats(loop, &l->data);
                position_templates(l->son, son_mats);
            }
        }
    }

    void Assemble::position_model(Chain &chain, const Mat &mat) {
        int len = chain.size();
        auto mat2 = model_mat(chain, std::list<int>{0, 1, len - 2, len - 1});
        auto sp = geom::suppos(mat2, mat);
        INIT_SUPPOS(sp);
        for (auto && res : chain) {
            for (auto && atom : res) {
                APPLY_SUPPOS(atom, sp);
            }
        }
//        EACH_ATOM(model, 
//            APPLY_SUPPOS(ATOM, sp);
//        );
    }

    Mat Assemble::model_mat(const Chain &chain, const Li &list) {
        Mat mat(_suppos_atoms.size() * list.size(), 3);
        int index = 0;
        for (uint n = 0; n < chain.size(); n++) {
            if (std::find(list.begin(), list.end(), n) != list.end()) {
                for (auto && atom : chain[n]) {
                    if (_suppos_atoms.find(atom.name) != _suppos_atoms.end()) {
                        for (uint i = 0; i < 3; i++) {
                            mat(index, i) = atom[i];
                        }
                        index++;
                    }
                }
            }
        }
        return mat;
    }

    std::list<Mat> Assemble::loop_mats(const Chain &chain, SSE *l) {
        std::list<Mat> mats;
        if (l->num_sons() == 0) {
            return mats;
        } else {
            for (auto && hinge : l->hinges) {
                int a = hinge.first,
                    b = hinge.second;
                mats.push_back(model_mat(chain, std::list<int>{a-1,a,b,b+1}));
            }
            return mats;
        }
    }

    void Assemble::find_records() {
		for (auto &&sse : _ss_tree) {
			_find_self[&sse] = false;
			find_loop_records(&sse);
			find_helix_records(&sse);
		}
    }

    void Assemble::print_records() {
        log << "Records searching results:" << std::endl;
		for (auto &&sse : _ss_tree) {
			log << "SSE(" << &sse << "):" << std::endl;
			if (sse.has_helix()) {
				log << "  Helix " << Str(sse.helix) << " (" << m_records[&sse].second.size() << ")" << std::endl;
			}
			if (sse.has_loop()) {
				log << "  Loop " << Str(sse.loop) << " (" << m_records[&sse].first.size() << ")" << std::endl;
			}
		}
    }

    void Assemble::find_loop_records(SSE *l) {
        if (l->loop.empty()) return;

        S seq = l->loop.seq(), ss = l->loop.ss(), p_ss = NASS::pure_ss(ss), lower_ss = NASS::lower_ss(p_ss, 1), family = _family;
        int num_sons = l->num_sons();

        S info_file = _lib + "/RNA/records/" + (l->is_open() ? "open_" : "") + "loop";
		std::ifstream ifile;
        FOPEN(ifile, info_file);

        TemplRec templ_rec;
        S line;
        int num = 0;
        while (std::getline(ifile, line) && set_loop_rec(templ_rec, line)) {
            if (std::find(m_disused_pdbs.begin(), m_disused_pdbs.end(), templ_rec._src) != m_disused_pdbs.end()) {
                continue;
            } else if (NASS::pure_ss(NASS::lower_ss(templ_rec._ss, 1)) == lower_ss) {
                templ_rec._score = (templ_rec._ss == ss ? 5 : 0);
                if (templ_rec._src == _name) {
                    templ_rec._score += 5;
                }
                if (templ_rec._family == family && family != "other") {
                    templ_rec._score += 2;
                }
                for (uint i = 0; i < templ_rec._seq.size(); i++) {
                    if (seq[i] == templ_rec._seq[i]) {
                        if (ss[i] == templ_rec._ss[i] && ss[i] != '.' && ss[i] != '(' && ss[i] != ')') {
                            templ_rec._score += 2;
                        } else if (ss[i] == '(' || ss[i] == ')') {
                            templ_rec._score += 0.2;
                        } else {
                            templ_rec._score += 1;
                        }
                    }
                }
                m_records[l].first.push_back(templ_rec);
                num++;
                if ((!(_source_pdb.empty())) && (templ_rec._src == _source_pdb)) {
                    m_records[l].first.clear(); 
                    m_records[l].first.push_back(templ_rec);
                    num = 1;
                    _find_self[l] = true;
                    break;
                }
            }
        }
        ifile.close();
        std::sort(m_records[l].first.begin(), m_records[l].first.end(), []( const TemplRec &loop1, const TemplRec &loop2) {
            return loop1._score > loop2._score; });
    }

    void Assemble::find_helix_records(SSE *l) {
        if (l->helix.empty()) return;

        S seq = l->helix.seq(), ss = l->helix.ss(), family = _family;
        int len = size(l->helix);

        S info_file = _lib + "/RNA/records/helix";
        std::ifstream ifile(info_file.c_str());
        if (!ifile) throw "jian::FindTemplates error! Can't find '" + info_file + "'!";

        TemplRec templ_rec;
        S line;
        int num = 0;
        while (std::getline(ifile, line) && set_helix_rec(templ_rec, line)) {
            /* if (std::find(_disused_pdbs.begin(), _disused_pdbs.end(), templ_rec._src) != _disused_pdbs.end()) {
                continue;
			} else */ 
			if (templ_rec._len == len) {
				templ_rec._score = 0;
				if (templ_rec._src == jian::upper(_name)) {
					templ_rec._score += 5;
				}
				if (templ_rec._family == family && family != "other") {
					templ_rec._score += 2;
				}
				for (uint i = 0; i < templ_rec._seq.size(); i++) {
					if (seq[i] == templ_rec._seq[i]) {
						templ_rec._score++;
					}
				}
				m_records[l].second.push_back(std::move(templ_rec));
				num++;
			}
        }
        ifile.close();
        std::sort(m_records[l].second.begin(), m_records[l].second.end(), []( const TemplRec &loop1, const TemplRec &loop2) {
            return loop1._score > loop2._score; });
    }

} // namespace nuc3d
END_JN

