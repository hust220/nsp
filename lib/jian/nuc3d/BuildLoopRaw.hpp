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

namespace jian {

class BuildLoopRaw {
public:
	struct Helix {
		num_t theta, phi;
	};
    using Hinge = std::array<int, 4>;
    using Hinges = std::deque<Hinge>;
    using Helices = std::vector<Helix>;

    std::string _lib = Env::lib();
    std::map<std::string, Hinges> _cache;
    std::string type = "RNA";
	std::shared_ptr<BuildLoopDG> m_build_loop_dg;
	std::string m_seq;
	std::string m_ss;
	Hinges m_hinges;
	Helices m_helices;

	BuildLoopRaw();

	BuildLoopRaw &init(const std::string &seq, const std::string &ss);

	Chain operator ()();

	Hinges ss_to_hinges(const std::string &ss);

	Model make_internal_loop(const Hinges &hinges);

	void print_helices(const Helices &helices);

	Helices make_helices(const Hinges &hinges);

	Helices make_helices_random(const Hinges &hinges);

	Chain helices_to_model(const Helices &helices, const Hinges &hinges);

	void complete_chain(Chain &chain);

	void append_helix(Chain &chain, const Model &helix, const Hinge &hinge);

	Model read_standard_pairs();

    template<typename O, typename N>
    void adjust_helix(Model &model, const O &o, const N &n, num_t theta_o, num_t phi_o, num_t theta_n, num_t phi_n) {
        auto rot = geom::suppos_axis_polar<num_t>(theta_o, phi_o, theta_n, phi_n);
        std::vector<num_t> v1 {-o[0], -o[1], -o[2]}, v2 {n[0], n[1], n[2]};
        LOG 
			<< o[0] << ':' << o[1] << ':' << o[2] << ' '
			<< n[0] << ':' << n[1] << ':' << n[2] << std::endl;
        for (auto &chain : model) for (auto &residue : chain) for (auto &atom : residue) {
            geom::translate(atom, v1);
			geom::rotate(atom, rot);
			geom::translate(atom, v2);
        }
    }

};

} // namespace jian

