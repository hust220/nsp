#pragma once

#include "string.hpp"
#include <regex>
#include <fstream>

namespace jian {

	struct file {
		static std::string name(const std::string & file_path);

		static std::string type(const std::string & file_path);

		static void clean(const std::string & file_name);

		template<typename T>
		static void open(T & stream, const char * file_name) {
			stream.open(file_name);
		}

		template<typename T>
		static void open(T & stream, const std::string & file_name) {
			stream.open(file_name.c_str());
		}
	};

#define FOPEN(f, s) \
	jian::file::open(f, s); \
    if (!(f)) {\
        std::ostringstream stream;\
        stream << "Can't open file '" << s << "'!" << std::endl;\
        throw stream.str();\
    }

#define FCLOSE(f) (f).close()

#define BEGIN_READ_FILE(f, t) do {\
	int N = 0;\
	std::ifstream ifile;\
	FOPEN(ifile, f);\
	std::string L;\
	tokenize_v F;\
	while (std::getline(ifile, L)) {\
		::jian::tokenize(L, F, t); \
		do


#define END_READ_FILE while(0);N++;}FCLOSE(ifile);}while(0)



} // namespace jian

