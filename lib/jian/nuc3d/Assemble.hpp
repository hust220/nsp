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
		Hairpin *l, records_t &records, S name = "",
		const pdbs_t &used_pdbs = {}, const pdbs_t &disused_pdbs = {}, family_t family = ""
	);

	void find_helix_records(
		Hairpin *l, records_t &records, Str name = "", Str family = ""
	);

	class Assemble : public TSP {
	public:
		//using Mat = Eigen::MatrixXd;
		using RecordsPair = Pair<records_t, records_t>;

		M<Hairpin *, RecordsPair> m_records;
		M<Hairpin *, Pair<Chain *, Chain *>> m_templates;
		M<Hairpin *, Pair<Q<record_t *>, Q<record_t *>>> m_records_cache;
		M<Hairpin *, Pair<Q<Chain *>, Q<Chain *>>> m_templates_cache;
		M<Hairpin *, Pair<TemplRec, TemplRec>> m_selected_record;
		I _it_num = 0;
		M<Hairpin *, B> _is_virtual;
		M<Hairpin *, B> _is_sampled;
		M<Hairpin *, B> _find_self;
		Set<Str> _suppos_atoms{ "C5*", "O3*", "C1*" };
		B m_sample = false;
		S m_sample_mode = "sample_one";
		S m_loop_building_method;
		I m_num_sse_templs;

		BuildLoopDG build_loop_dg;
		BuildLoopRaw build_loop_raw;
		SampleLoop sample_loop;

		Assemble(const Par &par);

		void set_virtual_loops();

		virtual void run();

		void predict_one();

		double score_templates();

		bool lack_templates();

		void print_templates();

		void assemble();

		void transform();

		void select_templates();

		void set_loop_template(Hairpin *l, bool is_first);

		void sample_one_template();

		void sample_all_templates();

		void sample_loop_template(Hairpin *l);

		void sample_helix_template();

		Hairpin *select_loop();

		Chain load_pdb(const TemplRec &templ_res, S type = "");

		void assemble_templates(Hairpin *l);

		void position_templates();

		void position_templates(Hairpin *l, L<Mat> mats);

		void position_model(Chain &chain, const Mat &mat);

		Mat model_mat(const Chain &chain, const Li &list);

		std::list<Mat> loop_mats(const Chain &chain, Hairpin *l);

		void find_records();

		void complete_records();

		void print_records();

		void find_loop_records(Hairpin *l);

		void find_helix_records(Hairpin *l);

	};

} // namespace nuc3d

END_JN

