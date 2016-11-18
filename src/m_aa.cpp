#include <memory>
#include <jian/pdb.hpp>
#include <jian/cg.hpp>
#include <jian/scoring/ParBp.hpp>
#include <jian/geom.hpp>
//#include <jian/mpi.hpp>
#include "nsp.hpp"

namespace jian {
	namespace {

		REGISTER_NSP_COMPONENT(aa) {
			std::string list_file, prefix;
			Chain chain;
			int l, n, i, j;
			std::shared_ptr<CG> cg;
			std::string content = "hi";

			list_file = par.get("l", "list");
			//prefix = par.get("prefix");
			cg.reset(CG::fac_t::create("6p"));

			n = 0;
			ParBp par_bp;
			BEGIN_READ_FILE(list_file, " ") {
				chain_read_model(chain, F[0]);
				chain = cg->to_cg(chain);
				l = chain.size();
				for (i = 0; i < l; i++) {
					for (j = i + 1; j < l; j++) {
						par_bp.anal(chain[i], chain[j]);

						if (par_bp.is_paired()) {
							std::cout
								<< i + 1 << '-' << chain[i].name << ' '
								<< j + 1 << '-' << chain[j].name << '\n'
								<< "theta: " << par_bp.theta << '\n'
								<< "o21_: \n"
								<< par_bp.o21_ << '\n'
								<< "o12_: \n"
								<< par_bp.o12_ << std::endl;
						}

						//int m, n;
						//for (m = 0; m < 6; m++) {
						//	for (n = 0; n < 6; n++) {
						//		d = geom::distance(chain[i][m], chain[j][n]);
						//		if (d < 3) {
						//			std::cout
						//				<< d << ' ' << chain[i].name << i + 1 << ' ' << chain[j].name << j + 1 << ' '
						//				<< chain[i][m].name << ' ' << m << ' ' << chain[j][n].name << ' ' << n
						//				<< std::endl;
						//		}
						//	}
						//}

						/*
						Chain c;
						c.push_back(chain[i]);
						c.push_back(chain[j]);
						std::ostringstream stream;
						stream << prefix << "." << n + 1 << ".pdb";
						mol_write(c, stream.str());
						std::cout
							<< n + 1 << ' ' << chain[i].name << i + 1 << ' ' << chain[j].name << j + 1 << ' '
							<< par_bp.d << ' ' << par_bp.o12_[2] << ' ' << par_bp.o21_[2] << ' ' << par_bp.alpha << ' ' << par_bp.theta
							<< par_bp.bb12_[0] << ' ' << par_bp.bb12_[1] << ' ' << par_bp.bb12_[2]
							<< std::endl;
						*/
						n++;
					}
				}
			} END_READ_FILE;
		}
	}
} // namespace jian

