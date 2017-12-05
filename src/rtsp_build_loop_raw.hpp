#pragma once

#include <numeric>
#include "pdb.hpp"
#include "geom.hpp"
#include "rss.hpp"
#include "log.hpp"
#include "rand.hpp"
#include "env.hpp"
#include "rtsp_parse_helix.hpp"
#include "rtsp_build_loop_dg.hpp"
#include "rtsp_mutate.hpp"

BEGIN_JN

class BuildLoopRaw {
public:
	struct Pos {
		Num theta, phi;
	};
	using Helix = Pos;
    using Helices = std::vector<Helix>;
    using Hinge = std::array<int, 4>;
    using Hinges = std::deque<Hinge>;
	using Frag = std::deque<int>;
	using Frags = std::deque<Frag>;

    Str _lib = Env::lib();
	std::map<Str, Hinges> m_cache_hinges;
	std::map<Str, Frags> m_cache_frags;
    Str type = "RNA";
	Str m_seq;
	Str m_ss;
	Hinges m_hinges;
	Helices m_helices;
	Frags m_frags;
	std::vector<Pos> m_res_pos;
	Num m_radius;

	BuildLoopRaw &init(const Str &seq, const Str &ss);

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
    void adjust_helix(Model &model, const O &o, const N &n, Num theta_o, Num phi_o, Num theta_n, Num phi_n) {
        auto rot = geom::suppos_axis_polar<Num>(theta_o, phi_o, theta_n, phi_n);
        std::vector<Num> v1 {-o[0], -o[1], -o[2]}, v2 {n[0], n[1], n[2]};
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

END_JN

