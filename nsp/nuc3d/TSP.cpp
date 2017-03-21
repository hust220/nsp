#include <sstream>
#include "TSP.hpp"
#include "jian/utils/Debug.hpp"
#include "jian/utils/Par.hpp"
#include "jian/utils/Env.hpp"
#include "jian/utils/rand.hpp"
#include "jian/utils/log.hpp"
#include "jian/utils/string.hpp"

#ifdef JN_PARA
#  include "jian/mpi.hpp"
#endif

BEGIN_JN

void TSP::set_logfile_name() {
	if (_par->has("log:off")) {
		log.file("");
	}
	else if (_par->has("log")) {
		log.file(_par->get("log"));
	}
	else {
#ifdef JN_PARA
		log.file(to_str(m_out_dir, '/', _name, ".p", mpi_rank() + 1, ".log"));
#else
		log.file(to_str(m_out_dir, '/', _name, ".p1.log"));
#endif
	}
}

void TSP::init(const Par &pars) {
	_par = new Par(pars);
	_lib = Env::lib();
	_name = pars.get("name", "job_name", "job");
    pars.set(m_out_dir, "out_dir");
	m_will_log = !_name.empty();
	if (m_will_log) set_logfile_name();
	log << "# Reading sequence" << std::endl;
	read_seq();
	log << "# Reading secondary structure" << std::endl;
	read_ss();
	read_cg();
	pars.set(_lib, "lib", "library_path");
	pars.set(_num, "n", "num", "number");
	pars.set(_hinge, "h", "hinge");
	pars.set(_strategy, "strategy");
	pars.set(_family, "family");
	pars.set(_type, "t", "type", "mol_type");
	pars.set(_seed, "seed");
#ifdef JN_PARA
	_seed += mpi_rank() * 123;
#endif
	seed(int(_seed));
	pars.set(m_mode, "mode");
	set_constraints();
	pars.set(_method, "method");
	pars.set(_is_test, "test");
	pars.set(_native, "native");
	pars.set(_source_pdb, "source_pdb");
	to_upper(_source_pdb);
}

TSP::~TSP() {
	delete _par;
	//log_pop();
}

void TSP::set_constraints() {
	_par->set(_file_distances, "distances");
	_par->set(_file_dca, "dca");
	_par->set(_file_contacts, "contacts");
   if (_par->has("constraints", "c")) {
      Str file = _par->get("constraints", "c");
      Str type = "dca";
      _par->set(type, "constraints_type", "ct");
      if (type == "dca") _file_dca = file;
      else if (type == "distances") _file_distances = file;
      else if (type == "contacts") _file_contacts = file;
   }
	_constraints.read_dca(_file_dca, int(_seq.size() * 0.5));
	_constraints.read_contacts(_file_contacts);
	_constraints.read_distances(_file_distances);
}

void TSP::predict() {
	log << "# Displaying starting information..." << std::endl;
	display_start_information();

	run();

	log << "# Writing final model..." << std::endl;
	write_final_model();

	log << "# Displaying ending information..." << std::endl;
	display_end_information();
}

void TSP::run() {
}

void TSP::write_final_model() {
#ifdef JN_PARA
	mol_write(_pred_chain, to_str(m_out_dir, '/', _name, ".pred.p", mpi_rank() + 1, ".pdb"));
#else
	mol_write(_pred_chain, to_str(m_out_dir, '/', _name, ".pred.p1.pdb"));
#endif
}


void TSP::display_start_information() {
	_start_time = std::time(nullptr);
	log << "=========================================================\n"
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
	log << "----------------------------------------\n"
		<< "Finish Time: " << std::asctime(std::localtime(&(_end_time))) << '\n'
		<< "Time elapsed: " << _end_time - _start_time << "s\n"
		<< "=========================================================\n\n";
}

void TSP::read_seq() {
	Str seq;
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

END_JN

