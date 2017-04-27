#pragma once

#include <string>
#include <iostream>
#include "traits.hpp"
#include "Par.hpp"

BEGIN_JN

extern int g_argc;
extern char **g_argv;
extern Par g_par;

class Env {
private:

	Env() = default;

	Env(const Env &env) = default;

	Env &operator =(const Env &env) = default;

	S m_path_lib;

public:

	static Env &instance() {
		static Env env;
		return env;
	}

	static S lib() {
		return Env::instance().get_lib();
	}

	void set_lib(S s) {
		m_path_lib = s;
	}

	S get_lib() {
		if (m_path_lib.empty()) {
			char *path = std::getenv("NSP");
			if (path == NULL) {
				std::cout << "Please tell me the path of the NSP library: ";
				std::cin >> m_path_lib;
			}
			else {
				m_path_lib = path;
			}
		}
		return m_path_lib;
	}
};

}
