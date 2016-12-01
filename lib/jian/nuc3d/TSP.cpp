#include <sstream>
#include "TSP.hpp"
#include "../utils/Debug.hpp"
#include "../utils/Par.hpp"
#include "../utils/Env.hpp"
#include "../utils/rand.hpp"
#include "../utils/log.hpp"
#include "../utils/string.hpp"

#ifdef JN_PARA
#include "../mpi.hpp"
#endif

namespace jian {

void TSP::init(const Par &pars) {
    _par = new Par(pars);
    _lib = Env::lib();
	_name = pars.get("job_name", "job", "name");
#ifdef JN_PARA
	if (g_mpi->m_size == 1) log_file(to_str(_name, ".tsp.log"));
	else log_file(to_str(_name, ".", g_mpi->m_rank + 1, ".tsp.log"));
#else
	log_file(to_str(_name, ".tsp.log"));
#endif
    //m_cmd = (*_par)[1];
	LOG << "# Reading sequence" << std::endl;
	read_seq();
    LOG << "# Reading secondary structure" << std::endl;
    read_ss();
	read_cg();
    pars.set(_lib, "lib", "library_path");
    pars.set(_num, "n", "num", "number");
    pars.set(_hinge, "h", "hinge");
    pars.set(_strategy, "strategy");
    pars.set(_family, "family");
    pars.set(_type, "t", "type", "mol_type");
    pars.set(_seed, "seed");
    seed(int(_seed));
    pars.set(m_mode, "mode");
    set_constraints();
	set_disused_pdbs();
    pars.set(_method, "method");
    pars.set(_is_test, "test");
    pars.set(_native, "native");
    pars.set(_source_pdb, "source_pdb");
    jian::to_upper(_source_pdb);
}

TSP::~TSP() {
    delete _par;
}

void TSP::set_disused_pdbs() {
	if (_par->has("disused_pdbs")) {
		_disused_pdbs = (*_par)["disused_pdbs"];
	}
	for (auto && s : _disused_pdbs) jian::to_upper(s);
}

void TSP::set_constraints() {
	_par->set(_file_distances, "distances");
	_par->set(_file_dca, "dca");
	_par->set(_file_contacts, "contacts");
	_constraints.read_dca(_file_dca, int(_seq.size() * 0.5));
	_constraints.read_contacts(_file_contacts);
	_constraints.read_distances(_file_distances);
}

void TSP::predict() {
	LOG << "# Displaying starting information..." << std::endl;
	display_start_information();

	run();

	LOG << "# Writing final model..." << std::endl;
	write_final_model();

	LOG << "# Displaying ending information..." << std::endl;
	display_end_information();
}

void TSP::run() {
}

void TSP::write_final_model() {
#ifdef JN_PARA
	if (g_mpi->m_size == 1) mol_write(_pred_chain, to_str(_name, ".pred.pdb"));
	else mol_write(_pred_chain, to_str(_name, ".", g_mpi->m_rank + 1, ".pred.pdb"));
#else
	mol_write(_pred_chain, to_str(_name, ".pred.pdb"));
#endif
}


void TSP::display_start_information() {
    _start_time = std::time(nullptr);
    LOG << "=========================================================\n"
              << "New Job: " << _name << '\n'
              << "Time: " << std::asctime(std::localtime(&(_start_time))) << '\n'
              << "Sequence: " << _seq << '\n'
              << "2D Structure: " << _ss << '\n'
              << "Molecular Type: " << _type << '\n'
              << "Number: " << _num << '\n'
              << "Method: " << _method << '\n'
              << "----------------------------------------\n";
}

void TSP::display_end_information() {
    _end_time = std::time(nullptr);
    LOG << "----------------------------------------\n"
              << "Finish Time: " << std::asctime(std::localtime(&(_end_time))) << '\n'
              << "Time elapsed: " << _end_time - _start_time << "s\n"
              << "=========================================================\n\n";
}

void TSP::read_seq() {
	str_t seq;
	jian::tokenize_v v;
	jian::tokenize_v w;

	_seq = "";
	_par->set(seq, "seq", "sequence");
	tokenize(seq, v, "&");
	if (v.size() > 0) {
		for (auto && s : v) {
			tokenize(s, w, ":");
			if (w.size() == 1 && s.back() != ':') {
				_seq += w[0];
				m_chain_lens.push_back(w[0].size());
				m_chain_names.push_back("A");
			}
			else if (w.size() == 2) {
				_seq += w[1];
				m_chain_lens.push_back(w[1].size());
				m_chain_names.push_back(w[0]);
			}
			else {
				throw "Illegal sequence!";
			}
		}
		jian::to_upper(_seq);
	}
	else {
		throw "Illegal sequence!";
	}
}

void TSP::read_cg() {
	m_cg_type = "6p";
	_par->set(m_cg_type, "cg");
	m_cg.reset(CG::fac_t::create(m_cg_type));
}

void TSP::read_ss() {
    _par->set(_ss, "ss", "secondary_structure");
}

} // namespace jian

