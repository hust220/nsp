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
#include "BuildLoopRaw.hpp"
#include "SampleLoop.hpp"
#include "TSP.hpp"

BEGIN_JN

namespace nuc3d {

	using record_t = TemplRec;
	using records_t = Deque<record_t>;
	using pdbs_t = Deque<Str>;
	using family_t = Str;

	void chain_read_record(Chain &chain, const record_t &templ_res);

	void find_loop_records(
		SSE *l, records_t &records, S name = "",
		const pdbs_t &used_pdbs = {}, const pdbs_t &disused_pdbs = {}, family_t family = ""
	);

	void find_helix_records(
		SSE *l, records_t &records, Str name = "", Str family = ""
	);

	class Assemble : public TSP {
	public:
		//using Mat = Eigen::MatrixXd;
		using RecordsPair = Pair<records_t, records_t>;

		M<SSE *, RecordsPair> m_records;
		M<SSE *, Pair<Chain, Chain>> m_templates;
		M<SSE *, Pair<Q<record_t *>, Q<record_t *>>> m_records_cache;
		M<SSE *, Pair<Q<Chain *>, Q<Chain *>>> m_templates_cache;
		M<SSE *, Pair<TemplRec, TemplRec>> m_selected_record;
		I _it_num = 0;
		M<SSE *, B> _is_virtual;
		M<SSE *, B> _is_sampled;
		M<SSE *, B> _find_self;
		Set<Str> _suppos_atoms{ "C5*", "O3*", "C1*" };
		B m_sample = false;
		S m_sample_mode = "sample_one";
		S m_loop_building_method;
		I m_num_sse_templs;

		BuildLoopDG build_loop_dg;
		BuildLoopRaw build_loop_raw;
		SampleLoop sample_loop;
        Par::pars_t m_disused_pdbs;

		Assemble(const Par &par);

		virtual void run();

		void predict_one();

		double score_templates();

		bool lack_templates();

		void print_templates();

		void assemble();

		void transform();

		void select_templates();

		void set_loop_template(SSE *l, Bool is_first);

		void set_helix_template(SSE *l, Bool is_first);

		void sample_one_template();

		void sample_all_templates();

		void sample_loop_template(SSE *l);

		void sample_helix_template();

		SSE *select_loop();

		Chain load_pdb(const TemplRec &templ_res, S type = "");

		void assemble_templates(SSE *l);

		void position_templates();

		void position_templates(SSTree::El *l, L<Mat> mats);

		void position_model(Chain &chain, const Mat &mat);

		Mat model_mat(const Chain &chain, const Li &list);

		std::list<Mat> loop_mats(const Chain &chain, SSE *l);

		void find_records();

		void print_records();

		void find_loop_records(SSE *l);

		void find_helix_records(SSE *l);

	};

} // namespace nuc3d

END_JN

