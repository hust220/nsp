#include "BuildLoopRaw.hpp"

namespace jian {

	BuildLoopRaw &BuildLoopRaw::init(const std::string &seq, const std::string &ss) {
		//Hinges &&hinges = ss_to_hinges(NASS::hinge_ss(ss));
		m_seq = seq;
		m_ss = ss;
		set_hinges();
		set_frags();
		print_frags();
		set_pos();
		return *this;
	}

	Chain BuildLoopRaw::operator ()() {
		return make_chain();
	}

	void BuildLoopRaw::set_hinges() {
		if (m_cache_hinges.count(m_ss)) {
			m_hinges = m_cache_hinges[m_ss];
			return;
		}

		Hinges hinges;
		std::deque<char> stack;
		std::deque<int> stack2;
		int i, j, l;

		for (i = 0; i < m_ss.size(); i++) {
			if (m_ss[i] == '(' || m_ss[i] == ')') {
				stack.push_back(m_ss[i]);
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
			throw std::string("jian::ss_to_hinges error! ss: ") + m_ss;
		}
		m_cache_hinges[m_ss] = hinges;
		m_hinges = hinges;
	}

	bool BuildLoopRaw::is_open() {
		return m_hinges.empty() || m_hinges.back()[0] != 0;
	}

	void BuildLoopRaw::set_frags() {
		int j;

		if (m_cache_frags.count(m_ss)) {
			m_frags = m_cache_frags[m_ss];
			return;
		}

		Frags frags;
		Frag frag;
		if (is_open()) {
			for (j = 0; j < size(m_ss) && m_ss[j] != '(' && m_ss[j] != ')'; j++) frag.push_back(j);
			frags.push_back(std::move(frag));
			while (j < size(m_ss)) {
				j += 4;
				for (; j < size(m_ss) && m_ss[j] != '(' && m_ss[j] != ')'; j++) frag.push_back(j);
				frags.push_back(std::move(frag));
			}
			//for (i = 0; i < size(m_hinges) - 1; i++) {
			//	for (j = m_hinges[i][3] + 1; j < size(m_ss) && m_ss[j] != '(' && m_ss[j] != ')'; j++) {
			//		frag.push_back(j);
			//	}
			//	frags.push_back(std::move(frag));
			//}
			//for (j = m_hinges[i][3] + 1; j < size(m_ss); j++) frag.push_back(j);
			//frags.push_back(std::move(frag));
		}
		else {
			for (j = 2; j < size(m_ss) - 2 && m_ss[j] != '(' && m_ss[j] != ')'; j++) frag.push_back(j);
			frags.push_back(std::move(frag));
			while (j < size(m_ss) - 2) {
				j += 4;
				for (; j < size(m_ss) - 2 && m_ss[j] != '(' && m_ss[j] != ')'; j++) frag.push_back(j);
				frags.push_back(std::move(frag));
			}
			//for (i = 0; i + 2 < size(m_hinges); i++) {
			//	for (j = m_hinges[i][3] + 1; j < size(m_ss) && m_ss[j] != '(' && m_ss[j] != ')'; j++) {
			//		frag.push_back(j);
			//	}
			//	frags.push_back(std::move(frag));
			//}
			//for (j = m_hinges[i][3] + 1; j + 2 < size(m_ss); j++) frag.push_back(j);
			//frags.push_back(std::move(frag));
		}

		m_cache_frags[m_ss] = frags;
		m_frags = frags;

	}

	void BuildLoopRaw::print_helices(const Helices &helices) {
		for (auto &helix : helices) {
			LOG << helix.theta << ':' << helix.phi << ' ';
		}
		LOG << std::endl;
	}

	void BuildLoopRaw::print_frags() {
		for (auto &&frag : m_frags) {
			LOG << "Frag: ";
			for (auto && i : frag) LOG << i << ' ';
			LOG << " End." << std::endl;
		}
	}

	void BuildLoopRaw::set_pos() {
		int i1, i2;
		num_t phi, j, n;

		n = std::accumulate(m_frags.begin(), m_frags.end(), 0.0, [](num_t n, auto &&frag) {
			return n + size(frag);
		}) + size(m_hinges) * 5 + 1;
		m_radius = 7.0 * n / (2.0 * PI);
		LOG << "radius: " << m_radius << std::endl;
		phi = 2.0 * PI / n;
		LOG << "phi: " << phi << std::endl;
		m_res_pos.resize(size(m_ss));
		m_helices.resize(size(m_hinges));

		i1 = 0;
		i2 = 0;
		j = 0;
		while (i1 < m_frags.size() || i2 < m_hinges.size()) {
			if (i1 < m_frags.size()) {
				for (auto && i : m_frags[i1]) {
					m_res_pos[i] = { PI / 2.0, phi * j };
					j++;
				}
				i1++;
			}
			if (i2 < m_hinges.size()) {
				j += 2;
				m_helices[i2] = { PI / 2.0, phi * j };
				j += 3;
				i2++;
			}
		}

		//int len = m_hinges.size();
		//Helices helices(len);
		//for (int i = 0; i < len; i++) {
		//	helices[i] = {i*2*PI/len, 0};
		//}
		//m_helices = std::move(helices);
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

	Chain BuildLoopRaw::make_chain() {
		Chain chain;
		chain.resize(m_ss.size());
		Model &&pairs = read_standard_pairs();
		auto result = parse_helix(pairs);
		int i;

		for (i = 0; i < m_helices.size(); i++) {
			Model helix = pairs;
			std::vector<num_t> origin {
				m_radius*std::sin(m_helices[i].theta)*std::cos(m_helices[i].phi),
				m_radius*std::sin(m_helices[i].theta)*std::sin(m_helices[i].phi),
				m_radius*std::cos(m_helices[i].theta)
			};
			if (!is_open() && i == size(m_helices) - 1) {
				adjust_helix(
					helix, result.origin, origin,
					result.theta, result.phi, m_helices[i].theta, m_helices[i].phi
				);
			}
			else {
				adjust_helix(
					helix, result.origin, origin,
					result.theta, result.phi, m_helices[i].theta, 180+m_helices[i].phi
				);
			}
			append_helix(chain, helix, m_hinges[i]);
		}

		for (i = 0; i < chain.size(); i++) {
			if (chain[i].empty()) {
				chain[i].push_back(Atom("C4*", 
					m_radius*std::sin(m_res_pos[i].theta)*std::cos(m_res_pos[i].phi),
					m_radius*std::sin(m_res_pos[i].theta)*std::sin(m_res_pos[i].phi),
					m_radius*std::cos(m_res_pos[i].theta)));
			}
		}

		complete_chain(chain);
		transform(chain);
		return chain;
	}

	void BuildLoopRaw::transform(Chain &chain) {
		Model m;
		m.push_back(chain);
		chain = std::move(jian::transform(m, m_seq, "RNA")[0]);
	}

	void BuildLoopRaw::complete_chain(Chain &chain) {
		Mat x, y;
		int i, j, l;

		l = chain.size();
		x.resize(l, 3);
		y.resize(l, 3);

		for (i = 0; i < l; i++) {
			Atom &atom = chain[i]["C4*"];
			for (j = 0; j < 3; j++) {
				x(i, j) = atom[j];
			}
		}
		Chain c = CG::fac_t::make("1p")->to_aa(x, 0, l - 1);
		for (i = 0; i < l; i++) {
			Atom &atom = c[i]["C4*"];
			for (j = 0; j < 3; j++) {
				y(i, j) = atom[j];
			}
		}

		geom::Superposition<num_t> sp(y, x);
		for (auto && res : c) {
			for (auto && atom : res) {
				sp.apply(atom);
			}
		}

		i = 0;
		for (auto && res : chain) {
			if (res.size() == 1) {
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