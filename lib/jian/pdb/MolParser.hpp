#pragma once

#include <map>
#include <string>
#include <functional>
#include <fstream>
#include <list>
#include "../utils/file.hpp"
#include "../utils/Factory.hpp"

namespace jian {

	struct MolParsedLine {
		std::string atom_type;
		std::string atom_name;
		std::string res_name;
		std::string chain_name;
		std::string res_flag;
		int atom_num;
		int res_num;
		int chain_num;
		int model_num;
		double x, y, z;
	};


	class MolParser {
	public:
		using creater_t = MolParser *(const std::string &);

		MolParsedLine *_curr_line = NULL;
		MolParsedLine *_next_line = NULL;
		std::list<MolParsedLine *> gc_line;
		std::string file_name;
		std::string file_type;
		std::string mol_type = "";
		std::ifstream ifile;

		MolParser(const std::string &f);

		~MolParser();

		MolParsedLine *parse_line();

		bool eof();

		virtual MolParsedLine *getline() = 0;

		static MolParser *make(const std::string &file_type, const std::string &file_path, std::string mol_type);

	};

#define REG_MOL_PARSER(name, Type) REGISTER_FACTORY(jian::MolParser::creater_t, name, Type)
	using FacMolParser = Factory<MolParser::creater_t>;

} // namespace jian

