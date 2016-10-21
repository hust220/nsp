#include "../utils/Par.hpp"
#include "../utils/Env.hpp"
#include "../utils/file.hpp"
#include "../jian.hpp"
#include "Ass.hpp"

namespace jian {
namespace nuc3d {

void chain_read_record(Chain &chain, const record_t &templ_res) {
    std::ostringstream stream;
    stream << Env::lib() << "/RNA/templates/" << templ_res._name << ".pdb";
    chain_read_model(chain, stream.str());
}

void find_loop_records(loop *l, records_t &records, std::string name,
                       const pdbs_t &used_pdbs, const pdbs_t &disused_pdbs, family_t family) {
    if (l->empty()) return;

    std::string seq = l->seq(),
                ss = l->ss(),
                p_ss = NucSS::pure_ss(ss),
                lower_ss = NucSS::lower_ss(p_ss, 1);

    int num_sons = l->num_sons();

    std::ostringstream stream;
    stream << Env::lib() << "/RNA/records/" << (l->is_open() ? "open_" : "") << "loop";
//LOGI << "FIND LOOP RECORDS OF: " << l << ' ' << stream.str() << std::endl;
	std::ifstream ifile;
	FOPEN(ifile, stream.str());

    record_t templ_rec;
    std::string line;
    int num = 0;
    while (std::getline(ifile, line) && set_loop_rec(templ_rec, line)) {
        if (std::find(disused_pdbs.begin(), disused_pdbs.end(), templ_rec._src) != disused_pdbs.end()) {
            continue;
        } else if ((!used_pdbs.empty()) &&
                   std::find(used_pdbs.begin(), used_pdbs.end(), templ_rec._src) == used_pdbs.end()) {
            continue;
        } else if (NucSS::pure_ss(NucSS::lower_ss(templ_rec._ss, 1)) == lower_ss) {
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

void find_helix_records(loop *l, records_t &records, std::string name, std::string family) {
    if (l->s.empty()) return;

    std::string seq = l->s.seq(),
                ss = l->s.ss();

    int len = l->s.len();

    std::ostringstream stream;
    stream << Env::lib() << "/RNA/records/helix";
	std::ifstream ifile;
    FOPEN(ifile, stream.str());

    record_t templ_rec;
    std::string line;
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
        JobPredict3D::init(par);

        m_sample = par.has("sample");
        par.set(m_sample_mode, "sample_mode");

        LOG << "# Check input" << std::endl;
        if (!NucSS::check_ss(_ss)) throw "The secondary structure includes illegal characters!";

        LOG << "# Construct 2D Structure Tree" << std::endl;
        _ss_tree.make(_seq, _ss, _hinge);

        LOG << "# Set Virtual Loops" << std::endl;
        set_virtual_loops();

        LOG << "# Searching Templates..." << std::endl;
        find_templates();

        LOG << "# Printing records..." << std::endl;
        print_records();

    }

    void Assemble::set_virtual_loops() {
        LOOP_TRAVERSE(_ss_tree.head(), 
            if (L->has_loop()) {
                if (_method == "FADG") {
                    _is_virtual[L] = false;
//                } else if (_strategy == "loose" && (L->is_open() || L->num_sons() >= 2)) {
//                    _is_virtual[L] = true;
                } else {
                    _is_virtual[L] = false;
                }
            }
        );
    }

    void Assemble::run() {
        std::ostringstream stream;
        int n;

        predict_one();
        LOG << "# Writing assembly structure." << std::endl;
        if (_par->has("out")) {
            mol_write(_pred_chain, _par->get("out"));
        } else {
            stream << _name << ".assemble.pdb";
            mol_write(_pred_chain, stream.str());
        }
        if (m_sample) {
            n = 1;
            LOG << "# Writing sampling structure " << n << std::endl;
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
                LOG << "# Writing sampling structure " << n << std::endl;
                stream.str("");
                stream << _name << ".assemble." << n << ".pdb";
                mol_write(_pred_chain, stream.str());
            }
        }
    }

    void Assemble::predict_one() {
        LOG << "# Select Templates" << std::endl;
        select_templates(); 

        LOG << "# Print Templates" << std::endl;
        print_templates();

        LOG << "# Assemble" << std::endl;
        assemble();
    }

    double Assemble::score_templates() {
        double e = 0;
        for (auto && pair : _records) {
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
        loop *l;
        for (auto && pair : _templates) {
            l = pair.first;
            if (l->has_loop()) {
                if (_records[l].first.empty()) return true;
            }
        }
        return false;
    }

    void Assemble::print_templates() {
        LOOP_TRAVERSE(_ss_tree.head(),
			LOG << L << " : " << _templates[L].first.model_name << ' ' << _templates[L].second.model_name << std::endl;
		);
    }

    void Assemble::assemble() {
        LOG << "## Score of Templates: " << score_templates() << std::endl;

        LOG << "## Print Templates." << std::endl;
        print_templates();

        LOG << "## Position Templates." << std::endl;
        position_templates();

        LOG << "## Assemble Templates." << std::endl;
        assemble_templates(_ss_tree.head());

        LOG << "## Transform." << std::endl;
        this->transform();
    }

    void Assemble::transform() {
        Model m;
        m.push_back(_pred_chain);
        _pred_chain = std::move(jian::transform(m, _seq, _type)[0]);
    }

    void Assemble::select_templates() {
        LOOP_TRAVERSE(_ss_tree.head(),
            if (L->has_loop()) {
                if (_records[L].first.empty()) {
                    build_loop_dg.init(L->seq(), NucSS::lower_ss(L->ss())); 
                    _templates[L].first = build_loop_dg();
                    m_selected_record[L].first = TemplRec{};
                } else {
                    _templates[L].first = load_pdb(_records[L].first[0]);    
                    m_selected_record[L].first = _records[L].first[0];
                }
            }
            if (L->has_helix()) {
                _templates[L].second = load_pdb(_records[L].second[0]);
                m_selected_record[L].second = _records[L].second[0];
            }
        );
    }

	void Assemble::sample_helix_template() {
		assert(_records.size() == 1);
		records_t & r = _records.begin()->second.second;
		uint n = static_cast<uint>(rand() * r.size());
		_templates[0].second = load_pdb(r[n]);
		m_selected_record[0].second = r[n];
	}

    void Assemble::sample_one_template() {
		loop * l = select_loop();
		if (l == NULL) {
			sample_helix_template();
		}
		else {
			sample_loop_template(select_loop());
		}
    }

    void Assemble::sample_all_templates() {
        for (auto && pair : _templates) {
            loop *l = pair.first;
            if (l->has_loop()) {
                sample_loop_template(l);
            }
        }
    }

    void Assemble::sample_loop_template(loop *l) {
        if (_records[l].first.empty()) {
            build_loop_dg.init(l->seq(), NucSS::lower_ss(l->ss())); 
            _templates[l].first = build_loop_dg();
            m_selected_record[l].first = TemplRec{};
        } else {
            int n = int(rand() * _records[l].first.size());
            _templates[l].first = load_pdb(_records[l].first[n]);    
            m_selected_record[l].first = _records[l].first[n];
        }
    }

    loop *Assemble::select_loop() {
        std::deque<loop *> ls;
        for (auto && pair : _records) {
            auto &l = pair.first;
            if (l->has_loop()) {
                ls.push_back(l);
            }
        }
		return (ls.empty() ? NULL : ls[static_cast<uint>(ls.size() * jian::rand())]);
    }

    Chain Assemble::load_pdb(const TemplRec &templ_res, std::string type) {
        Chain chain;
        chain_read_record(chain, templ_res);
        return chain;
    }

    void Assemble::assemble_templates(loop *l) {
        _pred_chain.resize(_seq.size());
        LOOP_TRAVERSE(l,
//L->print();
            if (L->has_loop()) {
                auto &temp_residues = _templates[L].first;
//                LOGI << "loop: \n" << temp_residues << std::endl;
                int i = 0;
                LOOP_EACH(L,
                    if (RES->type != '&') {
                        _pred_chain[RES->num - 1] = temp_residues[i];
                        i++;
                    }
                );
            }
            if (L->has_helix()) {
                int len = L->s.len();
                auto &temp_residues = _templates[L].second;
//                LOGI << "helix: \n" << temp_residues << std::endl;
                HELIX_EACH(L->s,
                    _pred_chain[BP->res1.num - 1] = temp_residues[N_BP];
                    _pred_chain[BP->res2.num - 1] = temp_residues[2*len-1-N_BP];
                );
            }
        );
    }

    void Assemble::position_templates() {
        position_templates(_ss_tree.head(), std::list<Mat>());
    }

    void Assemble::position_templates(loop *l, std::list<Mat> mats) {
        if (l == NULL) return;
        auto &loop = _templates[l].first;
        auto &helix = _templates[l].second;
        if (l->has_helix()) {
            // position helix
            if (!mats.empty()) {
                auto mat1 = mats.front();
                position_model(helix, mat1);
            }
            if (l->has_loop()) {
                // position loop
                int len = helix.size();
                auto mat2 = model_mat(helix, std::list<int>{len/2-2, len/2-1, len/2, len/2+1});
                position_model(loop, mat2);
                // position son
                if (l->has_son()) {
                    auto son_mats = loop_mats(loop, l);
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
                auto son_mats = loop_mats(loop, l);
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

    Mat Assemble::model_mat(const Chain &chain, const std::list<int> &list) {
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

    std::list<Mat> Assemble::loop_mats(const Chain &chain, loop *l) {
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

    void Assemble::find_templates() {
        LOOP_TRAVERSE(_ss_tree.head(), 
            _find_self[L] = false;
            find_loop_records(L); 
            find_helix_records(L)
        );
    }

    void Assemble::print_records() {
        LOG << "Records searching results:" << std::endl;
        LOOP_TRAVERSE(_ss_tree.head(),
            LOG << "Hairpin(" << L << "):" << std::endl;
            if (L->has_helix()) {
                LOG << "  Helix " << std::string(L->s) << " (" << _records[L].second.size() << ")" << std::endl;
            }
            if (L->has_loop()) {
                LOG << "  Loop " << std::string(*L) << " (" << _records[L].first.size() << ")" << std::endl;
            }
        );
    }

    void Assemble::find_loop_records(loop *l) {
        if (l->empty()) return;

        std::string seq = l->seq(), ss = l->ss(), p_ss = NucSS::pure_ss(ss), lower_ss = NucSS::lower_ss(p_ss, 1), family = _family;
        int num_sons = l->num_sons();

        std::string info_file = _lib + "/RNA/records/" + (l->is_open() ? "open_" : "") + "loop";
		std::ifstream ifile;
        FOPEN(ifile, info_file);

        TemplRec templ_rec;
        std::string line;
        int num = 0;
        while (std::getline(ifile, line) && set_loop_rec(templ_rec, line)) {
            if (std::find(m_disused_pdbs.begin(), m_disused_pdbs.end(), templ_rec._src) != m_disused_pdbs.end()) {
                continue;
            } else if (NucSS::pure_ss(NucSS::lower_ss(templ_rec._ss, 1)) == lower_ss) {
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
                _records[l].first.push_back(templ_rec);
                num++;
                if ((!(_source_pdb.empty())) && (templ_rec._src == _source_pdb)) {
                    _records[l].first.clear(); 
                    _records[l].first.push_back(templ_rec);
                    num = 1;
                    _find_self[l] = true;
                    break;
                }
            }
        }
        ifile.close();
        std::sort(_records[l].first.begin(), _records[l].first.end(), []( const TemplRec &loop1, const TemplRec &loop2) {
            return loop1._score > loop2._score; });
    }

    void Assemble::find_helix_records(loop *l) {
        if (l->s.empty()) return;

        std::string seq = l->s.seq(), ss = l->s.ss(), family = _family;
        int len = l->s.len();

        std::string info_file = _lib + "/RNA/records/helix";
        std::ifstream ifile(info_file.c_str());
        if (!ifile) throw "jian::FindTemplates error! Can't find '" + info_file + "'!";

        TemplRec templ_rec;
        std::string line;
        int num = 0;
        while (std::getline(ifile, line) && set_helix_rec(templ_rec, line)) {
            /* if (std::find(m_disused_pdbs.begin(), m_disused_pdbs.end(), templ_rec._src) != m_disused_pdbs.end()) {
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
				_records[l].second.push_back(std::move(templ_rec));
				num++;
			}
        }
        ifile.close();
        std::sort(_records[l].second.begin(), _records[l].second.end(), []( const TemplRec &loop1, const TemplRec &loop2) {
            return loop1._score > loop2._score; });
    }

} // namespace nuc3d
} // namespace jian

