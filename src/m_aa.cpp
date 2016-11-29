#include <memory>
#include <jian/pdb.hpp>
#include <jian/cg.hpp>
#include <jian/scoring/ParBp.hpp>
#include <jian/geom.hpp>
#include <jian/utils/Serial.hpp>
#include <jian/nuc3d/BuildChain.hpp>
//#include <jian/mpi.hpp>
#include "nsp.hpp"

namespace jian {
	namespace {

		REGISTER_NSP_COMPONENT(test_build_chain) {
			BuildChain builder;
			int n = 10;

			par.set(n, "n");
			builder(n);
			JN_OUT << builder.m_chain << std::endl;
		}


		REGISTER_NSP_COMPONENT(test_count) {
			std::vector<int> v{ 1, 2, 3 };
			std::vector<std::vector<int>> ls{ {1, 2}, {2, 3, 4} };
			each<int>(ls, [](auto i) {
				LOG << i << std::endl;
				return true;
			});
			LOG << count(v) << ' ' << count(ls) << ' ' << count<int>(ls) << std::endl;
		}

		REGISTER_NSP_COMPONENT(test_aa) {
			int n, l, i, j, k;
			Chain chain;
			std::shared_ptr<CG> cg;
			std::string list_file;

			list_file = par.get("l", "list");
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
							std::cout << chain[i].name << ' ' << chain[j].name << ' ';
							for (k = 0; k < 6; k++) {
								std::cout
									<< geom::distance(chain[i][k], chain[j][k]) << ' ';
							}
							std::cout << std::endl;
						}
						n++;
					}
				}
			} END_READ_FILE;
		}

		REGISTER_NSP_COMPONENT(aa) {
			std::string list_file, prefix;
			Chain chain;
			int l, n, i, j;
			std::shared_ptr<CG> cg;
			std::string content = "hi";

			Serial serial;
			int a = 1;
			double b = 3.5;
			std::string c = "hi";
			LOG << "stringify" << std::endl;
			std::string s = serial.stringify(a, b);
			LOG << "parse" << std::endl;
			short e, f;
			serial.parse(s, e, f, b, c);
			LOG << e << ' ' << f << ' ' << c << std::endl;

			return;

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

