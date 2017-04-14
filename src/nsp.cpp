#include "nsp.hpp"
#include "lua.hpp"

BEGIN_JN

NSP &NSP::instance() {
	static NSP nsp;
	return nsp;
}

template<typename _Methods>
static void show_help(Str filename, _Methods && m) {
    for (auto &&it : FileLines(filename)) {
        std::cout << it.line << std::endl;
    }
    for (auto && pair : m) { std::cout << pair.first << ' '; }
    std::cout << std::endl;
    std::cout << std::endl;
}

void NSP::run(int argc, char **argv) {
	Par par(argc, argv);
    auto g = par.getv("global");
	NSP::instance().m_argc = argc;
	NSP::instance().m_argv = argv;
	if (par.has("log_level")) {
		log_level(std::stoi(par["log_level"][0]));
	}
	if (par.has("log")) {
		log_file(par.get("log"));
	}
	if (par.has("out", "o")) {
		out_file(par.get("out", "o"));
	}
	auto &m = instance()._methods;
	S path = Env::lib() + "/RNA/pars/src/";
	if (par._orig_pars.size() <= 1) {
        show_help(path + "nsp.md", m);
	}
	else if (par.has("help") || par.has("h") || par.has("-help")) {
		S name = path + par[1] + ".md";
		for (auto &&it : FileLines(name)) {
			std::cout << it.line << std::endl;
		}
	}
	else {
        if (m.find(par[1]) != m.end()) {
            m[par[1]](par);
        }
        else {
            Str filename = to_str(Env::lib(), "/RNA/scripts/", par[1], ".lua");
            std::ifstream ifile(filename.c_str());
            if (ifile) {
                lua_run(filename, par);
            }
            else {
                show_help(path + "nsp.md", m);
            }
        }
	}
}

END_JN


