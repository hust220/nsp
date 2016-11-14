#include <memory>
#include <jian/pdb.hpp>
#include <jian/cg.hpp>
#include <jian/scoring/ParBp.hpp>
#include <jian/geom.hpp>
#include "mpi.h"
#include "nsp.hpp"

namespace jian {
	namespace {

		class MPI {
		public:
			MPI() {
				MPI_Init(&NSP::instance().m_argc, &NSP::instance().m_argv);
			}

			int rank() {
				int rank;
				MPI_Comm_rank(MPI_COMM_WORLD, &rank);
				return rank;
			}

			void send(std::string s, int rank) {
				MPI_Send(s.c_str(), s.size(), MPI_CHAR, rank, 0, MPI_COMM_WORLD);
			}

			std::string recv(int rank) {
				char hi[1024];
				MPI_Status status;
				//int n = MPI_Recv(hi, 1024, MPI_CHAR, rank, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
				MPI_Recv(hi, 1024, MPI_CHAR, rank, 0, MPI_COMM_WORLD, &status);
				int n = status.internal[0];
				return std::string(hi, n);
			}

			~MPI() {
				MPI_Finalize();
			}
		};
		REGISTER_NSP_COMPONENT(aa) {
			std::string list_file, prefix;
			Chain chain;
			int l, n, i, j;
			std::shared_ptr<CG> cg;
			std::string content = "hi";

			par.set(content, "content");

			MPI mpi;
			int rank = mpi.rank();
			if (rank == 0) {
				mpi.send(content, 1);
				LOG << "Rank 0 send '" << content << "'" << std::endl;
			}
			else if (rank == 1) {
				LOG << "Rank 1 received string " << mpi.recv(0) << " from Rank 0\n" << std::endl;
			}
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
							std::cout << i << ' ' << j << "\n" << par_bp.o21_ << std::endl;
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

