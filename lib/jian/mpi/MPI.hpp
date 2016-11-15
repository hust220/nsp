#pragma once

#ifdef JN_PARA

#include <string>
#include <memory>
#include "mpi.h"
#include "../utils/Env.hpp"

namespace jian {


	class MPI {
	public:
		int m_rank;
		int m_size;

		MPI() {
			MPI_Init(&g_argc, &g_argv);
			m_rank = rank();
			m_size = size();
		}

		int rank() {
			int rank;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			return rank;
		}

		int size() {
			int size;
			MPI_Comm_size(MPI_COMM_WORLD, &size);
			return size;
		}

		void send(std::string s, int rank) {
			MPI_Send(s.c_str(), s.size(), MPI_CHAR, rank, 0, MPI_COMM_WORLD);
		}

		std::string recv(int rank) {
			char hi[1024];
			MPI_Status status;
			//int n = MPI_Recv(hi, 1024, MPI_CHAR, rank, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			MPI_Recv(hi, 1024, MPI_CHAR, rank, 0, MPI_COMM_WORLD, &status);
			int n = ((int*)&status)[0];
			return std::string(hi, n);
		}

		~MPI() {
			MPI_Finalize();
		}
	};

	extern std::shared_ptr<MPI> g_mpi;

}

#endif
