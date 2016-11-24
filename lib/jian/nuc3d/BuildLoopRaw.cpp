#include "BuildLoopRaw.hpp"

namespace jian {

	BuildLoopRaw::BuildLoopRaw() {
		m_build_loop_dg = std::make_shared<BuildLoopDG>();
	}

	BuildLoopRaw &BuildLoopRaw::init(const std::string &seq, const std::string &ss) {
		//Hinges &&hinges = ss_to_hinges(NASS::hinge_ss(ss));
		m_seq = seq;
		m_ss = ss;
		m_hinges = ss_to_hinges(m_ss);
		m_helices = make_helices(m_hinges);
		//print_helices(m_helices);
		return *this;
	}

	Chain BuildLoopRaw::operator ()() {
		return helices_to_model(m_helices, m_hinges);
	}

	BuildLoopRaw::Hinges BuildLoopRaw::ss_to_hinges(const std::string &ss) {
		if (_cache.count(ss)) return _cache[ss];

		Hinges hinges;
		std::deque<char> stack;
		std::deque<int> stack2;
		int i, j, l;

		for (i = 0; i < ss.size(); i++) {
			if (ss[i] == '(' || ss[i] == ')') {
				stack.push_back(ss[i]);
				stack2.push_back(i);
				l = stack.size();
				if (l >= 4 && stack[l - 4] == '(' && stack[l - 3] == '(' && stack[l - 2] == ')' && stack[l - 1] == ')') {
					hinges.push_back({ stack2[l - 4], stack2[l - 3], stack2[l - 2], stack2[l - 1] });
					for (j = 0; j < 4; j++) {
						stack.pop_back();
						stack2.pop_back();
					}
				}
			}
			else {
				// pass
			}
		}
		if (!stack.empty()) {
			throw std::string("jian::ss_to_hinges error! ss: ") + ss;
		}
		_cache[ss] = hinges;
		return hinges;
	}

	Model BuildLoopRaw::make_internal_loop(const Hinges &hinges) {
		return Model();
	}

	void BuildLoopRaw::print_helices(const Helices &helices) {
		for (auto &helix : helices) {
			LOG << helix.theta << ':' << helix.phi << ' ';
		}
		LOG << std::endl;
	}

	BuildLoopRaw::Helices BuildLoopRaw::make_helices(const Hinges &hinges) {
		int len = hinges.size();
		Helices helices(len);
		for (int i = 0; i < len; i++) {
			helices[i] = {i*2*PI/len, 0};
		}
		return helices;
	}

	BuildLoopRaw::Helices BuildLoopRaw::make_helices_random(const Hinges &hinges) {
		int len = hinges.size();
		Helices helices(len);
		std::vector<int> v(len);
		std::iota(v.begin(), v.end(), 0);
		std::random_shuffle(v.begin(), v.end());
		helices[v[0]] = { 0, 0 };
		helices[v[1]] = { PI, 0 };
		for (int i = 2; i < len; i++) {
			helices[v[i]] = { int(rand() * 6)*PI / 6, int(rand() * 12)*PI / 6 };
		}
		return helices;
	}

	Chain BuildLoopRaw::helices_to_model(const Helices &helices, const Hinges &hinges) {
		//Model model;
		Chain chain;
		//model.resize(1);
		//model[0].resize(m_ss.size());
		chain.resize(m_ss.size());
		Model &&pairs = read_standard_pairs();
		auto result = parse_helix(pairs);
		num_t radius = 30;

		for (int i = 0; i < helices.size(); i++) {
			Model helix = pairs;
			std::vector<num_t> origin {
				radius*std::sin(helices[i].theta)*std::cos(helices[i].phi), 
				radius*std::sin(helices[i].theta)*std::sin(helices[i].phi),
				radius*std::cos(helices[i].theta)
			};
			adjust_helix(
				helix, result.origin, origin,
				result.theta, result.phi, helices[i].theta, helices[i].phi
			);
			append_helix(chain, helix, hinges[i]);
		}

		complete_chain(chain);
		return chain;
	}

	void BuildLoopRaw::complete_chain(Chain &chain) {
		Mat x, y;
		int i, j, n, l;
		std::vector<int> brokens;

		l = std::accumulate(chain.begin(), chain.end(), 0, [](int n, const Residue &res) {
			return n + (res.empty() ? 0 : 1);
		});
		x.resize(l, 3);
		y.resize(l, 3);

		for (auto && hinge : m_hinges) {
			if (hinge[0] - 1 >= 0) {
				brokens.push_back(hinge[0] - 1);
			}
		}
		m_build_loop_dg->init(chain, brokens);
		Chain &&c = (*m_build_loop_dg)();
		//JN_OUT << c << std::endl;
		i = 0;
		n = 0;
		for (auto && res : chain) {
			if (!res.empty()) {
				Atom &atom1 = c[i]["C4*"];
				Atom &atom2 = chain[i]["C4*"];
				for (j = 0; j < 3; j++) {
					x(n, j) = atom1[j];
					y(n, j) = atom2[j];
				}
				n++;
			}
			i++;
		}

		geom::Superposition<num_t> sp(x, y);
		for (auto && res : c) {
			for (auto && atom : res) {
				sp.apply(atom);
			}
		}

		i = 0;
		for (auto && res : chain) {
			if (res.empty()) {
				res = c[i];
			}
			i++;
		}
	}

	void BuildLoopRaw::append_helix(Chain &c, const Model &helix, const Hinge &hinge) {
		int n_res = 0;
		for (auto && chain : helix) for (auto && res : chain) {
			c[hinge[n_res]] = res;
			n_res++;
		}
	}

	Model BuildLoopRaw::read_standard_pairs() {
		return mol_read_to<Model>(_lib + "/" + type + "/pars/nuc3d/BuildLoop/pairs.pdb");
	}


}