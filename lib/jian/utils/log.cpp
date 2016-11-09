#include "log.hpp"

namespace jian {

	void out_file(std::string filename) {
		OUTMANAGER.set_file(filename);
	}

	std::string out_file() {
		return OUTMANAGER.get_file();
	}

	void log_file(std::string file_name) {
		LOGMANAGER.set_file(file_name);
	}

	std::string log_file() {
		return LOGMANAGER.get_file();
	}

	void log_level(int l) {
		LOGMANAGER.set_level(l);
	}

	int log_level() {
		return LOGMANAGER.get_level();
	}

}

