#include "../utils/Env.hpp"
#include "../pdb.hpp"
#include "../nuc2d.hpp"
#include "split.hpp"

BEGIN_JN
namespace nuc3d {

class Split {
public:
    using nums_t = std::vector<int>;

    S _name;
    S _seq;
    S _ss;
    S _family;
    S _lib;
    Model _mol;
    int _hinge;
    int _helix_num;
    int _loop_num;

    Split(Par par) {
        _helix_num = 0;
        _loop_num = 0;

        par.set(_name, "name", "pdb", "pdb_file");
        assert(!_name.empty());
        mol_read(_mol, _name);
        _name = _mol.name;
        _seq = jian::seq(_mol);

        par.set(_ss, "ss", "secondary_structure");
        assert(!_ss.empty());

        _lib = Env::lib();
        par.set(_lib, "lib", "library", "library_path");

        _family = "other";
        par.set(_family, "family");

        _hinge = 2;
        par.set(_hinge, "hinge");

        NASS::check_ss(_ss);
    }

    void run() {
        SSTree ss_tree;
        ss_tree.make(_seq, _ss, _hinge);
		for (auto &&sse : ss_tree) {
			parse_helix(&sse);
			parse_loop(&sse);
		}
    }

    void parse_helix(SSE *l) {
        if (l->has_helix()) {
            _helix_num++;

            // write to pdb file
            int len = size(l->helix);
            S pdb_path = _lib + "/RNA/templates/";
            S helix_name = _name + "-" + std::to_string(_helix_num) + "-helix-" + std::to_string(len);
            std::vector<int> nums(len * 2);
            int index = 0;
			for (auto && b : l->helix) {
            //for (auto b = l->helix.head; b != NULL; b = b->next) {
                nums[index] = b.res1.num - 1;
                nums[2 * len - 1 - index] = b.res2.num - 1;
                index++;
            }
            mol_write(sub(_mol, nums), pdb_path + helix_name + ".pdb");

            // write to records file
            S helix_info_file = _lib + "/RNA/records/helix";
            std::ofstream ofile(helix_info_file.c_str(), std::ios::app);
            ofile << helix_name << ' ' << size(l->helix) << ' ' << l->helix.seq() << ' ' << l->helix.ss() << ' ' << _family << std::endl;
            ofile.close();
        } else {
            // pass
        }
    }

    S get_loop_name(SSE *l, S name, int loop_num) {
        auto type = l->hinges.size();
        std::ostringstream stream;
        stream << name << '-' << loop_num << '-' << (l->is_open() ? "open_" : "") << "loop-" << type;
        return stream.str();
    }

    nums_t get_loop_nums(SSE *l) {
        nums_t nums;
		for (auto && res : l->loop) {
			nums.push_back(res.num - 1);
		}
        return nums;
    }

    void parse_loop(SSE *l) {
        if (l->has_loop()) {
            _loop_num++;

            //int len = size(l->loop);
            S pdb_path = _lib + "/RNA/templates/";
            S loop_name = get_loop_name(l, _name, _loop_num);
            mol_write(sub(_mol, get_loop_nums(l)), pdb_path + loop_name + ".pdb");

            // write to info file
            S l_ss = l->loop.ss();
            S loop_info_file = _lib + "/RNA/records/" + (l->is_open() ? "open_" : "") + "loop";

            std::ofstream ofile(loop_info_file.c_str(), std::ios::app);
            ofile << loop_name << ' ' << l->num_sons() << ' ' << l->loop.seq() << ' ' << l->loop.ss() << ' ' << _family << std::endl;
            ofile.close();
        }
    }

};

void split(const Par &par) {
    Split impl(par);
    impl.run();
}

} // namespace nuc3d
END_JN


