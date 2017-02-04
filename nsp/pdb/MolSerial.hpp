#pragma once

#include <fstream>
#include <sstream>
#include "Molecule.hpp"

BEGIN_JN

	class MolSerial {
	public:
		std::fstream m_stream;

		MolSerial(S out, std::ios::openmode mode) {
			m_stream.open(out, std::ios::binary | mode);
		}

		~MolSerial() {
			m_stream.close();
		}

		void write(char c) {
			m_stream.write(reinterpret_cast<const char*>(&c), sizeof(char));
		}

		void read(char &c) {
			m_stream.read(reinterpret_cast<char*>(&c), sizeof(char));
		}

		template<typename T, JN_ENABLE(std::is_integral<T>::value)>
		void write(T n) {
			short m = short(n);
			m_stream.write(reinterpret_cast<const char*>(&m), sizeof(short));
		}

		template<typename T, JN_ENABLE(std::is_integral<T>::value)>
		void read(T &n) {
			short m;
			m_stream.read(reinterpret_cast<char*>(&m), sizeof(short));
			n = T(m);
		}

		template<typename T, JN_ENABLE(std::is_floating_point<T>::value)>
		void write(T n) {
			float m = float(n);
			m_stream.write(reinterpret_cast<const char*>(&m), sizeof(float));
		}

		template<typename T, JN_ENABLE(std::is_floating_point<T>::value)>
		void read(T &n) {
			float m;
			m_stream.read(reinterpret_cast<char*>(&m), sizeof(float));
			n = T(m);
		}

		void write(const S &s) {
			//write(s.size());
			m_stream.write(s.c_str(), s.size() + 1);
		}

		void read(S &s) {
			char c[2];
			std::ostringstream stream;

			while (true) {
				m_stream.read(c, 1);
				if (c[0] == '\0') break;
				stream << c[0];
			}
			s = stream.str();
		}

		template<typename T1, typename T2, typename... Rest>
		void write(T1 && t1, T2 && t2, Rest && ...rest) {
			write(t1);
			write(t2, rest...);
		}

		template<typename T1, typename T2, typename... Rest>
		void read(T1 && t1, T2 && t2, Rest && ...rest) {
			read(t1);
			read(t2, rest...);
		}

		void write(const Atom &atom) {
			write('A', atom.name, atom.num, atom.mass, atom[0], atom[1], atom[2]);
		}

		void read(Atom &atom) {
			read(atom.name, atom.num, atom.mass, atom[0], atom[1], atom[2]);
		}

		void write(const Residue &res) {
			write('R', res.name, res.num, res.m_cg);
			for (auto && atom : res) {
				write(atom);
			}
			write('r');
		}

		void read(Residue &res) {
			char c;

			read(res.name, res.num, res.m_cg);
			while (true) {
				read(c);
				if (c == 'r') {
					return;
				}
				else if (c == 'A') {
					res.push_back(Atom{});
					read(res.back());
				}
				else {
					throw "This is not a standard '.jn' file! Illegal " + c;
				}
			}
		}

		void write(const Chain &chain) {
			write('C', chain.name, chain.type, chain.model_name, chain.m_cg);
			for (auto && res : chain) {
				write(res);
			}
			write('c');
		}

		void read(Chain &chain) {
			char c;

			read(chain.name, chain.type, chain.model_name, chain.m_cg);
			while (true) {
				read(c);
				if (c == 'c') {
					return;
				}
				else if (c == 'R') {
					chain.push_back(Residue{});
					read(chain.back());
				}
				else {
					throw "This is not a standard '.jn' file! Illegal " + c;
				}
			}
		}

		void write(const Model &model) {
			write('M', model.name, model.num, model.type, model.m_cg);
			for (auto && chain : model) {
				write(chain);
			}
			write('m');
		}

		void read(Model &model) {
			char c;

			read(model.name, model.num, model.type, model.m_cg);
			while (true) {
				read(c);
				if (c == 'm') {
					return;
				}
				else if (c == 'C') {
					model.push_back(Chain{});
					read(model.back());
				}
				else {
					throw "This is not a standard '.jn' file! Illegal " + c;
				}
			}
		}

		void write_mol(S filename) {
			for_each_model(filename, [this](const Model &model, int i) {
				write(model);
			});
			write('o');
		}

		void write(const Molecule &mol) {
			for (auto && model : mol) {
				write(model);
			}
			write('o');
		}

		void read(Molecule &mol) {
			char c;
			while (true) {
				read(c);
				if (c == 'o') {
					break;
				}
				else if (c == 'M') {
					mol.push_back(Model{});
					read(mol.back());
				}
				else {
					throw "This is not a standard '.jn' file ! Illegal " + c;
				}
			}

		}

	};

}
