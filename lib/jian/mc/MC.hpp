#pragma once

#include <cmath>
#include <iostream>
#include "../utils/rand.hpp"
#include "../utils/Par.hpp"
#include "../utils/Env.hpp"

namespace jian {

	class MC {
	public:
		enum mc_state_t {
			MC_READY,
			MC_HEATING,
			MC_COOLING,
			MC_WARMING,
			MC_DONE
		};
		double _mc_init_tempr;
		double _mc_tempr;
		int _mc_cycle_steps;
		int _mc_write_steps;
		int _mc_step = 0;
		double _mc_local_succ_rate;
		double _mc_en;
		int _mc_heat_steps;
		int _mc_cool_steps;
		double _mc_heat_rate;
		mc_state_t _mc_state = MC_READY; // 0: ready, 1: heating, 2: cooling, 3: done
		double _mc_dec_rate;
		std::string _mc_queue;

		MC();

		template<typename Fn>
		void mc_base(int steps, Fn &&ctrl_tempr) {
			int local_succ_num = 0;
			_mc_en = mc_total_energy();
			_mc_step = 0;
			mc_write();
			for (; _mc_step < steps; _mc_step++) {
				mc_select();
				auto &&en_old = mc_partial_energy();
				mc_sample();
				auto &&en_new = mc_partial_energy();
				auto &&en_diff = en_new - en_old;
				if (en_new > en_old && jian::rand() > std::exp(-en_diff / _mc_tempr)) {
					mc_back();
				}
				else {
					_mc_en += en_diff;
					local_succ_num++;
				}

				if (_mc_step % _mc_write_steps == _mc_write_steps - 1) {
					mc_write();
				}
				if (_mc_step % _mc_cycle_steps == _mc_cycle_steps - 1) {
					_mc_local_succ_rate = double(local_succ_num) / _mc_cycle_steps;
					local_succ_num = 0;
					if (!ctrl_tempr()) break;
				}
			}
			mc_write();
		}

		virtual void mc_write();
		void mc_run();
		virtual void mc_heat(int);
		virtual void mc_cool(int);
		virtual void mc_warm(int);
		virtual double mc_total_energy();
		virtual double mc_partial_energy();
		virtual void mc_select();
		virtual void mc_sample();
		virtual void mc_back();
		bool mc_is_heating() const;
		bool mc_is_cooling() const;

	};

} // namespace jian

