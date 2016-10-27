#include <iostream>
#include <iomanip>
#include <string>
#include <set>
#include <list>
#include "../utils/file.hpp"
#include "io.hpp"
#include "molmanip.hpp"
#include "../pp.hpp"


namespace jian {

	struct Line {
		int atom_num = 1, residue_num = 1, model_num = 1;
		std::string atom_name = "X", residue_name = "X", chain_name = "X", atom_name_label = "X";
		double x = 0, y = 0, z = 0, a = 1, b = 0;
		std::ostream stream;

		Line() : stream(std::cout.rdbuf()) {}

		Line(std::ostream &out) : stream(out.rdbuf()) {}

		void bind_stream(std::ostream &s) {
			stream.rdbuf(s.rdbuf());
		}

		void write() {
			std::replace(atom_name.begin(), atom_name.end(), '*', '\'');
			if (atom_name == "O1P") atom_name = "OP1";
			if (atom_name == "O2P") atom_name = "OP2";
			if (std::count_if(chain_name.begin(), chain_name.end(), [](auto && c) {return c != ' '; }) == 0) {
				chain_name = "X";
			}
			stream
				<< std::fixed
				<< "ATOM"
				<< std::setw(7) << atom_num
				<< "  "
				<< std::left
				<< std::setw(4) << atom_name
				<< std::right
				<< std::setw(3) << residue_name
				<< std::setw(2) << chain_name
				<< std::setw(4) << residue_num
				<< std::setprecision(3)
				<< std::setw(12) << x
				<< std::setw(8) << y
				<< std::setw(8) << z
				<< std::setprecision(2)
				<< std::setw(6) << a
				<< std::setw(6) << b
				<< std::setw(12) << atom_name_label
				<< "  "
				<< std::endl;
		}

		void read(const Atom &atom) {
			atom_name = atom.name;
			x = atom[0];
			y = atom[1];
			z = atom[2];
			atom_name_label = atom.name.substr(0, 1);
		}

		void write(const Atom &atom) {
			read(atom);
			write();
			atom_num++;
		}

		void write(const Residue &residue) {
			residue_name = residue.name;
			for (auto && atom : residue) {
				write(atom);
			}
			residue_num++;
			residue_name = "X";
		}

		void write(const Chain &chain) {
			chain_name = chain.name;
			for (auto &&residue : chain) {
				write(residue);
			}
			write_chain_end();
			residue_name = "X";
			chain_name = "X";
		}

		void write(const Model &model) {
			write_model_begin();
			for (auto && chain : model) {
				write(chain);
			}
			write_model_end();
			model_num++;
			atom_num = 1;
			residue_num = 1;
			residue_name = "X";
			chain_name = "X";
		}

		void write(const Molecule &mol) {
			for (auto && model : mol) {
				write(model);
			}
			write_file_end();
		}

		void write_model_begin() {
			stream
				<< std::left
				<< std::setw(13) << "MODEL"
				<< model_num
				<< std::right
				<< std::endl;
		}

		void write_model_end() {
			stream << "ENDMDL" << std::endl;
		}

		void write_file_end() {
			stream << "END" << std::endl;
		}

		void write_chain_end() {
			stream
				<< "TER "
				<< std::setw(7) << atom_num
				<< "  "
				<< std::left
				<< std::setw(4) << " "
				<< std::right
				<< std::setw(3) << residue_name
				<< std::setw(2) << chain_name
				<< std::setw(4) << residue_num - 1
				<< std::endl;
		}

	};

	bool diff_model(const molstream &parser) {
		auto && line = parser._next_line;
		auto && old_line = parser._curr_line;
		return line->model_num != old_line->model_num;
	}

	bool diff_chain(const molstream &parser) {
		auto && line = parser._next_line;
		auto && old_line = parser._curr_line;
		return line->chain_name != old_line->chain_name ||
			diff_model(parser);
	}

	bool diff_residue(const molstream &parser) {
		auto && line = parser._next_line;
		auto && old_line = parser._curr_line;
		return line->res_num != old_line->res_num || line->res_name != old_line->res_name || line->res_flag != old_line->res_flag ||
			diff_chain(parser);
	}

	bool diff_atom(const molstream &parser) {
		auto && line = parser._next_line;
		auto && old_line = parser._curr_line;
		return line->atom_num != old_line->atom_num || line->atom_name != old_line->atom_name ||
			diff_residue(parser);
	}

	void chain_read_model(Chain &s, std::string f, std::string type) {
		molstream *parser = molstream::make(jian::file::type(f), f, type);
		do {
			(*parser) >> s;
		} while (!parser->eof() && !diff_model(*parser));
		delete parser;
	}

	Chain read_model_to_chain(std::string f, std::string type) {
		Chain chain;
		chain_read_model(chain, f, type);
		return chain;
	}

	void append_chain_to_file(const Chain &chain, const std::string &file_name, int n) {
		std::ofstream output(file_name.c_str(), std::ios::app);
		output << "MODEL " << n << std::endl;
		output << chain;
		output << "ENDMDL" << std::endl;
		output.close();
	}

	molstream &operator >> (molstream &parser, Atom &atom) {
		MolParsedLine *line = parser.parse_line();
		if (line != NULL) {
			atom.init(line->atom_name, line->x, line->y, line->z, line->atom_num);
		}
		return parser;
	}

	molstream &operator >> (molstream &parser, Residue &residue) {
		//    std::cout << "molstream >> residue" << std::endl;
		residue.name = parser._next_line->res_name;
		residue.num = parser._next_line->res_num;
		for (int i = 0; !parser.eof(); i++) {
			if (i == 0 || !diff_residue(parser)) {
				Atom atom;
				parser >> atom;
				residue.push_back(atom);
			}
			else {
				break;
			}
		}
		return parser;
	}

	molstream &operator >> (molstream &parser, Chain &chain) {
		//    std::cout << "molstream >> chain" << std::endl;
		chain.name = parser._next_line->chain_name;
		chain.model_name = parser.file_name;
		for (int i = 0; !parser.eof(); i++) {
			if (i == 0 || !diff_chain(parser)) {
				Residue residue;
				parser >> residue;
				if (res_is_type(residue, parser.mol_type)) {
					chain.push_back(residue);
				}
			}
			else {
				break;
			}
		}
		return parser;
	}

	molstream &operator >> (molstream &parser, Model &model) {
		//    std::cout << "molstream >> model" << std::endl;
		model.num = parser._next_line->model_num;
		model.name = parser.file_name;
		for (int i = 0; !parser.eof(); i++) {
			if (i == 0 || !diff_model(parser)) {
				Chain chain;
				parser >> chain;
				model.push_back(chain);
			}
			else {
				break;
			}
		}
		return parser;
	}

	molstream &operator >> (molstream &parser, Molecule &mol) {
		mol.name = parser.file_name;
		for (int i = 0; !parser.eof(); i++) {
			Model model;
			parser >> model;
			mol.push_back(model);
		}
		return parser;
	}


	std::ostream &operator <<(std::ostream &output, const Molecule &mol) {
		Line(output).write(mol);
		return output;
	}

	std::ostream &operator <<(std::ostream &output, const Model &model) {
		Line l(output);
		l.write(model);
		l.write_file_end();
		return output;
	}

	std::ostream &operator <<(std::ostream &output, const Chain &chain) {
		Line l(output);
		l.write_model_begin();
		l.write(chain);
		l.write_model_end();
		l.write_file_end();
		return output;
	}

	std::ostream &operator <<(std::ostream &output, const Residue &residue) {
		Line(output).write(residue);
		return output;
	}

	std::ostream &operator <<(std::ostream &output, const Atom &atom) {
		Line(output).write(atom);
		return output;
	}

}

