#include "../utils/log.hpp"
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
			if (w.size() == 2 || w.size() == 3) {
				if (w.size() == 3) {
					_mc_tempr = JN_DBL(w[2]);
				}
				steps = JN_INT(w[1]);
				if (w[0] == "heat") {
					mc_heat(steps);
				}
				else if (w[0] == "cool") {
					mc_cool(steps);
				}
				else if (w[0] == "warm") {
					mc_warm(steps);
				}
				else {
					throw w[0] + ": illegal mc action!";
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
		mc_base(steps, [&]() {
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
		mc_base(steps, [&]() {
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
		mc_base(steps, [&]() {
			return true;
		});
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

