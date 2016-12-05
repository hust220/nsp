#pragma once

#include "string.hpp"
#include <regex>
#include <fstream>

BEGIN_JN

struct file {
	static S name(const S & file_path);

	static S type(const S & file_path);

	static void clean(const S & file_name);

	template<typename T>
	static void open(T & stream, const char * file_name) {
		stream.open(file_name);
	}

	template<typename T>
	static void open(T & stream, const S & file_name) {
		stream.open(file_name.c_str());
	}
};

#define FOPEN(fstream, filename) \
	jian::file::open(fstream, filename); \
    if (!(fstream)) {\
        std::ostringstream stream;\
        stream << "Can't open file '" << filename << "'!" << std::endl;\
        throw stream.str();\
    }

#define FCLOSE(fstream) (fstream).close()

#define BEGIN_READ_FILE(f, t) do {\
	int N = 0;\
	std::ifstream ifile;\
	FOPEN(ifile, f);\
	S L;\
	tokenize_v F;\
	while (std::getline(ifile, L)) {\
		::jian::tokenize(L, F, t); \
		do

#define END_READ_FILE while(0);N++;}FCLOSE(ifile);}while(0)

END_JN

