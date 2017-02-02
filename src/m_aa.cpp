#include <memory>
#include <nsp/pdb.hpp>
#include <nsp/cg.hpp>
#include <nsp/scoring/ParBp.hpp>
#include <jian/geom.hpp>
#include <jian/utils/Serial.hpp>
#include <nsp/nuc3d/BuildChain.hpp>
#include <jian/utils/file.hpp>
//#include <jian/mpi.hpp>
#include "nsp.hpp"

BEGIN_JN

namespace {

	REGISTER_NSP_COMPONENT(test_file) {
		auto file = FileLines("aa.txt");
		JN_OUT << (file.begin() == file.end()) << STD_ endl;
		for (auto &&it : FileLines("aa.txt")) {
			JN_OUT << it.line << STD_ endl;
		}
	}

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
		S list_file;

		list_file = par.get("l", "list");
		cg.reset(CG::fac_t::create("6p"));

		n = 0;
		ParBp par_bp;
		for (auto &&it : FileLines(list_file)) {
			chain_read_model(chain, it.arr[0]);
			chain = cg->to_cg(chain);
			l = chain.size();
			for (i = 0; i < l; i++) {
				for (j = i + 1; j < l; j++) {
					par_bp.anal(chain[i], chain[j]);
					if (par_bp.is_paired()) {
						JN_OUT << chain[i].name << ' ' << chain[j].name << ' ';
						//std::cout << par_bp.theta;
						for (k = 0; k < 6; k++) {
							JN_OUT
								<< geom::distance(chain[i][k], chain[j][k]) << ' ';
						}
						JN_OUT << std::endl;
					}
					n++;
				}
			}
		}
	}

	REGISTER_NSP_COMPONENT(aa) {
		S list_file, prefix;
		Chain chain;
		int l, n, i, j;
		std::shared_ptr<CG> cg;
		S content = "hi";

		Serial serial;
		int a = 1;
		double b = 3.5;
		S c = "hi";
		LOG << "stringify" << std::endl;
		S s = serial.stringify(a, b);
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
		for (auto &&it : FileLines(list_file)) {
			chain_read_model(chain, it.arr[0]);
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
		}
	}
}

END_JN

