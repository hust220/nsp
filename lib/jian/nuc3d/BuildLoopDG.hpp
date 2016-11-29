#pragma once

#include "../utils/Debug.hpp"
#include "../pdb.hpp"
#include "../nuc2d.hpp"
#include "../dg.hpp"
#include "../cg.hpp"
#include "HelixPar.hpp"

namespace jian {

	Chain *build_chain_dg(str_t seq, str_t ss);

	class BuildLoopDG {
	public:
		Eigen::MatrixXd _dist_bound;
		DihBound _dih_bound;

		DG m_dg;
		HelixPar helix_par;
		std::shared_ptr<CG> m_cg;

		BuildLoopDG();

		void init_dist_bound(Mat &b, int l);

		BuildLoopDG &init(const str_t &seq, const str_t &ss);

		BuildLoopDG &init(const Chain &c, const std::vector<int> &brokens);

		Chain operator ()();

		void set_bound_loop(Mat &b, DihBound &d, loop *l);

		void set_bound_helix(Mat &b, DihBound &d, const helix &h);
	};

} // namespace jian

