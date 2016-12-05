#pragma once

#include <map>
#include <string>
#include <functional>
#include <fstream>
#include <list>
#include "../utils/file.hpp"
#include "../utils/Factory.hpp"

BEGIN_JN

	struct MolParsedLine {
		S atom_type;
		S atom_name;
		S res_name;
		S chain_name;
		S res_flag;
		int atom_num;
		int res_num;
		int chain_num;
		int model_num;
		double x, y, z;
	};


	class MolParser {
	public:
		using creater_t = MolParser *(const S &);

		MolParsedLine *_curr_line = NULL;
		MolParsedLine *_next_line = NULL;
		std::list<MolParsedLine *> gc_line;
		S file_name;
		S file_type;
		S mol_type = "";
		std::ifstream ifile;

		MolParser(const S &f);

		~MolParser();

		MolParsedLine *parse_line();

		bool eof();

		virtual MolParsedLine *getline() = 0;

		static MolParser *make(const S &file_type, const S &file_path, S mol_type);

	};

#define REG_MOL_PARSER(name, Type) REGISTER_FACTORY(jian::MolParser::creater_t, name, Type)
	using FacMolParser = Factory<MolParser::creater_t>;

END_JN

