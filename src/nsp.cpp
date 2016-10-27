#include "nsp.hpp"

namespace jian {

	NSP &NSP::instance() {
		static NSP nsp;
		return nsp;
	}

	void NSP::run(int argc, char **argv) {
		Par par(argc, argv);
		if (par.has("log_level")) {
			set_log_level(std::stoi(par["log_level"][0]));
		}
		if (par.has("log")) {
			log_file(par.get("log"));
		}
		auto &m = instance()._methods;
		std::string path = Env::lib() + "/RNA/pars/src/";
		if (par._orig_pars.size() <= 1 || m.find(par[1]) == m.end()) {
			std::string name = path + "nsp.md";
			EACH_LINE(name, std::cout << L << std::endl;);
			for (auto && pair : m) { std::cout << pair.first << ' '; }
			std::cout << std::endl;
		}
		else if (par.has("help") || par.has("h") || par.has("-help")) {
			std::string name = path + par[1] + ".md";
			EACH_LINE(name, std::cout << L << std::endl;);
		}
		else {
			m[par[1]](par);
		}
	}

} // namespace jian


