#pragma once

#include <string>
#include <iostream>

namespace jian {

	class Env {
	private:

		Env() = default;

		Env(const Env &env) = default;

		Env &operator =(const Env &env) = default;

		std::string m_path_lib;

	public:

		static Env &instance() {
			static Env env;
			return env;
		}

		static std::string lib() {
			return Env::instance().get_lib();
		}

		void set_lib(std::string s) {
			m_path_lib = s;
		}

		std::string get_lib() {
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
