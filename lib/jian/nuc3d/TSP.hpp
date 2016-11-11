#pragma once

#include <string>
#include <ctime>
#include <fstream>
#include <memory>
#include "Constraints.hpp"
#include "../nuc2d/SSTree.hpp"
#include "../utils/log.hpp"
#include "../utils/Par.hpp"
#include "../dca.hpp"
#include "../pdb.hpp"
#include "../cg.hpp"

namespace jian {

class TSP {
public:
    std::string _name = "assemble";
    std::string _seq;
    std::string _ss;
    std::string _lib;
    std::string _family = "other";
    std::string _type = "RNA";
    std::string _file_distances;
	std::string _file_dca;
	std::string _file_contacts;
	std::string m_cg_type;
	std::shared_ptr<CG> m_cg;
    Constraints _constraints;
    std::deque<std::string> _disused_pdbs;
	std::string m_out_pdb;
	Chain _pred_chain;
    std::time_t _start_time, _end_time;
    int _hinge = 2;
    int _num = 1;
    int _num_sampling = 100;
    bool _is_test = false; // Is this a test case or not?
    std::string _method = "FA";
    std::string _native;
    SSTree _ss_tree;
    std::string _strategy = "loose";
    std::string _source_pdb;
    double _seed = 11;
    std::string m_mode = "auto";
    const Par *_par = NULL;
    std::string m_cmd;
	std::vector<int> m_chain_lens{};
	std::vector<std::string> m_chain_names{};

    TSP() = default;
    ~TSP();
    void init(const Par &pars);
    void set_constraints();
	void set_disused_pdbs();
	void predict();
	void write_final_model();
    void display_start_information();
    void display_end_information();
	void read_seq();
	void read_cg();
	virtual void run();
    virtual void read_ss();

};

} // namespace jian

