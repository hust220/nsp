#pragma once

#ifdef JN_PARA

#include <string>
#include <memory>
#include "mpi.h"
#include "../utils/Env.hpp"

BEGIN_JN

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

	void send(S s, int rank) {
		MPI_Send(s.c_str(), s.size(), MPI_CHAR, rank, 0, MPI_COMM_WORLD);
	}

	S recv(int rank) {
		char hi[1024];
		MPI_Status status;
		//int n = MPI_Recv(hi, 1024, MPI_CHAR, rank, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		MPI_Recv(hi, 1024, MPI_CHAR, rank, 0, MPI_COMM_WORLD, &status);
		int n = ((int*)&status)[0];
		return S(hi, n);
	}

	~MPI() {
		MPI_Finalize();
	}
};

extern SP<MPI> g_mpi;

END_JN

#endif
