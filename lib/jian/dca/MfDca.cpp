#include "MfDca.hpp"

BEGIN_JN
	namespace dca {

		REG_DCA_FAC("mf", MfDca);

		MfDca::MfDca() : Dca() {}

		MfDca::MfDca(S mol_type, float pw) : Dca(mol_type, pw) {}

		void MfDca::calculate_C() {
			int i, j, k, l;
			C = Matf::Zero(M * (q - 1), M * (q - 1));
			for (i = 0; i < M; i++) for (j = 0; j < M; j++) for (k = 0; k < q - 1; k++) for (l = 0; l < q - 1; l++) {
				C(i*(q - 1) + k, j*(q - 1) + l) = Pij(i, j, k, l) - Pi(i, k) * Pi(j, l);
			}
		}

		void MfDca::calculate_eij() {
			calculate_C();
			eij = C.inverse();
			for (int i = 0; i < C.rows(); i++) {
				for (int j = 0; j < C.cols(); j++) {
					eij(i, j) = std::exp(-(eij(i, j)));
				}
			}
		}

		void MfDca::set_mu(const Matf &m, const Vecf &pi, const Vecf &pj, Vecf &mu1, Vecf &mu2) {
			float epsilon = 1e-4f;
			float diff = 1.0f;
			Vecf v1, v2;
			float sum1, sum2;
			int i;
			float d;

			while (diff > epsilon) {
				diff = 0;
				v1 = m * mu2;
				v2 = m.transpose() * mu1;
				sum1 = 0;
				sum2 = 0;
				for (i = 0; i < q; i++) {
					v1[i] = pi[i] / v1[i];
					v2[i] = pj[i] / v2[i];
					sum1 += v1[i];
					sum2 += v2[i];
				}
				for (i = 0; i < q; i++) {
					v1[i] /= sum1;
					d = std::fabs(mu1[i] - v1[i]);
					if (d > diff) diff = d;
					v2[i] /= sum2;
					d = std::fabs(mu2[i] - v2[i]);
					if (d > diff) diff = d;
				}
				mu1 = v1;
				mu2 = v2;
			}
		}

		float MfDca::cal_di(int i, int j) {
			// set mu
			Vecf mu1 = Vecf::Constant(q, 1.0f / q);
			Vecf mu2 = Vecf::Constant(q, 1.0f / q);
			Matf Pdir = Matf::Ones(q, q);
			Pdir.block(0, 0, q - 1, q - 1) = eij.block(i*(q - 1), j*(q - 1), q - 1, q - 1);
			Vecf pi = Pi.row(i);
			Vecf pj = Pi.row(j);
			set_mu(Pdir, pi, pj, mu1, mu2);

			// calculate Pdir

			int a, b;
			float sum = 0;
			for (a = 0; a < q; a++) for (b = 0; b < q; b++) {
				Pdir(a, b) *= mu1(a) * mu2(b);
				sum += Pdir(a, b);
			}
			for (a = 0; a < q; a++) for (b = 0; b < q; b++) {
				Pdir(a, b) /= sum;
			}

			// calulate DI
			float DI = 0;
			for (a = 0; a < q; a++) for (b = 0; b < q; b++) {
				DI += Pdir(a, b) * std::log(Pdir(a, b) / pi[a] / pj[b]);
			}
			return DI;
		}

	} // namespace dca
END_JN


