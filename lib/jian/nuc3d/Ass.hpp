#pragma once

#include <iostream>
#include <set>
#include <sstream>
#include "../pdb.hpp"
#include "../geom.hpp"
#include "../nuc2d.hpp"
#include "../utils/rand.hpp"
#include "TemplRec.hpp"
#include "transform.hpp"
#include "BuildLoopDG.hpp"
#include "JobPredict3D.hpp"

namespace jian {
namespace nuc3d {

using record_t = TemplRec;
using records_t = std::deque<record_t>;
using pdbs_t = std::deque<std::string>;
using family_t = std::string;

void chain_read_record(Chain &chain, const record_t &templ_res);

void find_loop_records(loop *l, records_t &records, std::string name = "",
                       const pdbs_t &used_pdbs = {}, const pdbs_t &disused_pdbs = {}, family_t family = "");

void find_helix_records(loop *l, records_t &records, std::string name = "", std::string family = "");

class Assemble : public JobPredict3D {
public:    
    using Mat = Eigen::MatrixXd;
	using RecordsPair = std::pair<records_t, records_t>;

    std::map<loop *, RecordsPair> _records;
    std::map<loop *, std::pair<Chain, Chain>> _templates;
    std::map<loop *, std::pair<TemplRec, TemplRec>> m_selected_record;
    int _it_num = 0;
    std::map<loop *, bool> _is_virtual;
    std::map<loop *, bool> _is_sampled;
    std::map<loop *, bool> _find_self;
    Chain _pred_chain;
    std::set<std::string> _suppos_atoms{"C5*", "O3*", "C1*"};
    bool m_sample = false;
    std::string m_sample_mode = "sample_one";

    BuildLoopDG build_loop_dg;

    Assemble(const Par &par);

    void set_virtual_loops();

    void run();

    void predict_one();

    double score_templates();

    bool lack_templates();

    void print_templates();

    void assemble();

    void transform();

    void select_templates();

    void sample_one_template();

    void sample_all_templates();

    void sample_loop_template(loop *l);

	void sample_helix_template();

    loop *select_loop();

    Chain load_pdb(const TemplRec &templ_res, std::string type = "");

    void assemble_templates(loop *l);

    void position_templates();

    void position_templates(loop *l, std::list<Mat> mats);

    void position_model(Chain &chain, const Mat &mat);

    Mat model_mat(const Chain &chain, const std::list<int> &list);

    std::list<Mat> loop_mats(const Chain &chain, loop *l);

    void find_templates();

    void print_records();

    void find_loop_records(loop *l);

    void find_helix_records(loop *l);

};

} // namespace nuc3d
} // namespace jian

