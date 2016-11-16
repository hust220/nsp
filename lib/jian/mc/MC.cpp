#include <algorithm>
#include <numeric>
#include "../utils/log.hpp"
#include "../utils/Serial.hpp"
#include "MC.hpp"

namespace jian {

	MC::MC() {
		Par par(Env::lib() + "/RNA/pars/mc/mc.par");
		par.set(_mc_cycle_steps, "mc_cycle_steps");
		par.set(_mc_write_steps, "mc_write_steps");
		par.set(_mc_heat_steps, "mc_heat_steps");
		par.set(_mc_cool_steps, "mc_cool_steps");
		par.set(_mc_heat_rate, "mc_heat_rate");
		par.set(_mc_init_tempr, "mc_init_tempr");
		par.set(_mc_dec_rate, "mc_dec_rate");
	}

	void MC::mc_next_step() {
		_mc_step++;
	}

	bool MC::mc_is_heating() const {
		return _mc_state == MC_HEATING;
	}

	bool MC::mc_is_cooling() const {
		return _mc_state == MC_COOLING;
	}

	void MC::mc_write() {
		LOG << __FUNCTION__ << "..." << std::endl;
	}

	void MC::mc_run() {
		tokenize_v v, w;
		int steps;

		LOG << __FUNCTION__ << "..." << std::endl;
		tokenize(_mc_queue, v, "+");
		_mc_tempr = _mc_init_tempr;
		for (auto && s : v) {
			tokenize(s, w, ":");
			if (w.size() >= 2) {
				if (w[0] == "heat") {
					steps = JN_INT(w[1]);
					if (w.size() == 3) _mc_tempr = JN_DBL(w[2]);
					mc_heat(steps);
				}
				else if (w[0] == "cool") {
					steps = JN_INT(w[1]);
					if (w.size() == 3) _mc_tempr = JN_DBL(w[2]);
					mc_cool(steps);
				}
				else if (w[0] == "warm") {
					steps = JN_INT(w[1]);
					if (w.size() == 3) _mc_tempr = JN_DBL(w[2]);
					mc_warm(steps);
				}
				else if (w[0] == "remc") {
					steps = JN_INT(w[1]);
					_mc_lowest_tempr = JN_DBL(w[2]);
					_mc_highest_tempr = JN_DBL(w[3]);
					mc_remc(steps);
				}
				else {
					throw w[0] + ": illegal mc mode!";
				}
			}
			else {
				throw "jian::MC::mc_run error!";
			}
		}
	}

	void MC::mc_heat(int steps) {
		LOG << __FUNCTION__ << "..." << std::endl;
		_mc_state = MC_HEATING;
		mc_base(steps, [this]() {
			if (_mc_local_succ_rate < _mc_heat_rate) {
				_mc_tempr *= 1.1;
				return true;
			}
			else if (_mc_local_succ_rate > _mc_heat_rate) {
				_mc_tempr *= 0.9;
				return true;
			}
			else {
				return true;
			}
		});
	}

	void MC::mc_cool(int steps) {
		LOG << __FUNCTION__ << "..." << std::endl;
		_mc_state = MC_COOLING;
		int n_rate = 0, n_en = 0;
		double en = 0;
		mc_base(steps, [this, &en, &n_en, &n_rate]() {
			_mc_tempr *= _mc_dec_rate;
			if (std::fabs(en - _mc_en) <= 0.05) n_en++; else n_en = 0;
			en = _mc_en;
			if (_mc_local_succ_rate <= 0.05) n_rate++; else n_rate = 0;
			return n_rate < 5 && n_en < 5;
		});
	}

	void MC::mc_warm(int steps) {
		LOG << __FUNCTION__ << "..." << std::endl;
		_mc_state = MC_WARMING;
		mc_base(steps, [this]() {
			return true;
		});
	}

	void MC::mc_remc(int steps) {
		LOG << __FUNCTION__ << "..." << std::endl;
#ifndef JN_PARA
		throw "REMC is only supported in the parallel version! "
			"You may need to recompile nsp if you really want to use REMC!";
#else
		_mc_state = MC_REMC;
		_mc_tempr = _mc_lowest_tempr + g_mpi->m_rank * (_mc_highest_tempr - _mc_lowest_tempr) / (g_mpi->m_size - 1.0);
		Serial serial;

		auto exchange_tempr = [](auto &&tempr, auto &&en) {
			int l = tempr.size();
			double d;

			for (int i = 0; i < l - 1; i++) {
				d = en[i] - en[i + 1];
				d *= 1.0 / tempr[i] - 1.0 / tempr[i + 1];
				d = std::exp(d);
				if (rand() < std::min(1.0, d)) {
					std::swap(tempr[i], tempr[i+1]);
					i++;
				}
			}
		};

		mc_base(steps, [this, &serial, &exchange_tempr]() {
			if (g_mpi->m_rank == 0) {
				std::vector<double> tempr(g_mpi->m_size);
				std::vector<double> en(g_mpi->m_size);
				tempr[0] = _mc_tempr;
				en[0] = _mc_en;
				for (int i = 1; i < g_mpi->m_size; i++) {
					serial.parse(g_mpi->recv(i), tempr[i], en[i]);
				}
				exchange_tempr(tempr, en);
				_mc_tempr = tempr[0];
				for (int i = 1; i < g_mpi->m_size; i++) {
					g_mpi->send(serial.stringify(tempr[i]), i);
				}
			}
			else {
				g_mpi->send(serial.stringify(_mc_tempr, _mc_en), 0);
				serial.parse(g_mpi->recv(0), _mc_tempr);
			}
			return true;
		});
#endif
	}

	double MC::mc_total_energy() {
		return 0;
	}

	double MC::mc_partial_energy() {
		return rand() * 20;
	}

	void MC::mc_select() {
		LOG << __FUNCTION__ << "..." << std::endl;
	}

	void MC::mc_sample() {
		LOG << __FUNCTION__ << "..." << std::endl;
	}

	void MC::mc_back() {
		LOG << __FUNCTION__ << "..." << std::endl;
	}

} // namespace jian

