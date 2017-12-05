#pragma once

#include <string>
#include <ctime>
#include <fstream>
#include <memory>
#include "rtsp_constraints.hpp"
#include "rss_ss_tree.hpp"
#include "log.hpp"
#include "par.hpp"
#include "dca.hpp"
#include "pdb.hpp"
#include "cg.hpp"

BEGIN_JN

class TSP {
public:
	Log log;
    S _name = "assemble";
    S _seq;
    S _ss;
    S _lib;
    S _family = "other";
    S _type = "RNA";
    S _file_distances;
	S _file_dca;
	S _file_contacts;
	S m_cg_type;
	SP<CG> m_cg;
    Constraints _constraints;
	S m_out_pdb;
	Chain _pred_chain;
    std::time_t _start_time, _end_time;
    I _hinge = 2;
    I _num = 1;
    I _num_sampling = 100;
    SSTree _ss_tree;
    S _source_pdb;
    D _seed = 11;
    const Par *_par = NULL;
    S m_cmd;
	Vi m_chain_lens{};
	Vs m_chain_names{};
	B m_will_log;
    Str m_out_dir = ".";

    TSP() = default;
    ~TSP();
	void set_logfile_name();
    void init(const Par &pars);
    void set_constraints();
	void predict();
	void write_final_model();
    void display_start_information();
    void display_end_information();
	void read_seq();
	void read_cg();
	virtual void run();
    virtual void read_ss();

};

END_JN

