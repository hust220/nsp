#pragma once

#include <numeric>
#include "../pdb.hpp"
#include "../geom.hpp"
#include "../nuc2d.hpp"
#include "../utils/log.hpp"
#include "../utils/rand.hpp"
#include "../utils/Env.hpp"
#include "ParseHelix.hpp"
#include "BuildLoopDG.hpp"
#include "transform.hpp"

namespace jian {

class BuildLoopRaw {
public:
	struct Pos {
		num_t theta, phi;
	};
	using Helix = Pos;
    using Helices = std::vector<Helix>;
    using Hinge = std::array<int, 4>;
    using Hinges = std::deque<Hinge>;
	using Frag = std::deque<int>;
	using Frags = std::deque<Frag>;

    std::string _lib = Env::lib();
	std::map<std::string, Hinges> m_cache_hinges;
	std::map<std::string, Frags> m_cache_frags;
    std::string type = "RNA";
	std::string m_seq;
	std::string m_ss;
	Hinges m_hinges;
	Helices m_helices;
	Frags m_frags;
	std::vector<Pos> m_res_pos;
	num_t m_radius;

	BuildLoopRaw &init(const std::string &seq, const std::string &ss);

	Chain operator ()();

	void set_hinges();

	void set_frags();

	bool is_open();

	void print_helices(const Helices &helices);

	void print_frags();

	void set_pos();

	void transform(Chain &chain);

	Helices make_helices_random(const Hinges &hinges);

	Chain make_chain();

	void complete_chain(Chain &chain);

	void append_helix(Chain &chain, const Model &helix, const Hinge &hinge);

	Model read_standard_pairs();

    template<typename O, typename N>
    void adjust_helix(Model &model, const O &o, const N &n, num_t theta_o, num_t phi_o, num_t theta_n, num_t phi_n) {
        auto rot = geom::suppos_axis_polar<num_t>(theta_o, phi_o, theta_n, phi_n);
        std::vector<num_t> v1 {-o[0], -o[1], -o[2]}, v2 {n[0], n[1], n[2]};
        LOG 
			<< o[0] << ',' << o[1] << ',' << o[2] << ':' << theta_o << ',' << phi_o << ' '
			<< n[0] << ':' << n[1] << ':' << n[2] << ':' << theta_n << ',' << phi_n << std::endl;
        for (auto &chain : model) for (auto &residue : chain) for (auto &atom : residue) {
            geom::translate(atom, v1);
			geom::rotate(atom, rot);
			geom::translate(atom, v2);
        }
    }

};

} // namespace jian

