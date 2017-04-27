#pragma once

#include <iostream>
#include <fstream>
#include <regex>
#include "../pdb/Model.hpp"
#include "../geom.hpp"

BEGIN_JN

class Script {
public:
    using Word = struct {S value; int type;};
    using Var = struct {int addr; int type;};
    
    enum word_type {WORD_NUMBER, WORD_STRING, WORD_VARIABLE, WORD_SYMBOL};
    enum var_type {VAR_NUMBER, VAR_STRING, VAR_ATOM, VAR_RES, VAR_MODEL};

    std::map<std::string, Var> _vars;
    std::map<int, int> _num_refs;
    std::string<std::string, Var> _strings;
    std::vector<std::vector<int>> _states {
        {0, 1, 2, 3, 4, 5}, {-1, 1, 1, -1, -1, -1}, {-1, 2, 2, -1, -1, -1}, 
        {3, 3, 3, -1, 3, 3}, {4, 4, 4, 4, -1, 4}, {-1, -1, -1, -1, -1, -1}
    };

    S _line;

    Var var(const Word &word) {
        if (word.type == WORD_NUMBER) return Var{JN_INT(word.value), VAR_NUMBER};
        else if (word.type == WORD_STRING) {
            if (_strings.count(word.value)) {
                return _strings[word.value];
            } else {
                S *s = new std::string(word.value);
                _strings[word.value] = Var{s, VAR_STRING};
            }
        } else if (word.type == WORD_VARIABLE) {
            return _vars[word.value];
        }
    }

    void set_var(const S &name, void * addr,  int type) {
        _vars[name].addr = addr;
        _vars[name].type = type;
    }

    template<typename T>
    void free_var(T &&var) {
    }

    template<typename T, typename U>
    Var make_var(T &&addr, U &&type) {
        return Var{std::forward<T>(addr), std::forward<T>(type)};
    }

    template<typename T, typename U>
    Word make_word(T &&value, U &&type) {
        return Word{std::forward<T>(value), std::forward<T>(type)};
    }

    int char_index(const char &c) {
        if (c == ' ' || c == '\t') return 0;
        else if (isalpha(c)) return 1;
        else if (isdigit(c)) return 2;
        else if (c == '"') return 3;
        else if (c == '\'') return 4;
        else return 5;
    }

    template<typename T>
    Word read_word(T &&ifile, const S &s = "") {
        int state = 0, new_state;
        char c;
        S word;
        while (ifile) {
            c = ifile.get();
            if ((c != '\n' && c != '\r') && (_line.back() == '\n' || _line.back() == '\r')) _line = "";
            _line += c;
            new_state = _states[state][char_index(c)];
            if (new_state == 4) {
                if (state == 1) {
                    ifile.unget();
                    if (s != "" && s != word) { throw "Syntax error!"; }
                    return make_word(word, WORD_VARIABLE);
                } else if (state == 2) {
                    ifile.unget();
                    if (s != "" && s != word) { throw "Syntax error!"; }
                    return make_word(word, WORD_NUMBER);
                } else if (state == 3) {
                    if (s != "" && s != word+c) { throw "Syntax error!"; }
                    return make_word(word.substr(1), WORD_STRING);
                } else if (state == 4) {
                    if (s != "" && s != word+c) { throw "Syntax error!"; }
                    return make_word(word.substr(1), WORD_STRING);
                } else {
                    if (s != "" && s != word) { throw "Syntax error!"; }
                    ifile.unget();
                    return make_word(word, WORD_SYMBOL);
                }
            } else if (new_state != 0) word += c;
            state = new_state;
        }
    }

    template<typename T>
    std::function<Name()> read_each_atom(T &&ifile) {
        auto each_atom = [&](auto &&atom_name, auto &&model_name, auto &&f)->Name{
            EACH_ATOM(Model(model_name),
                this->set_var(atom_name, &ATOM, 0);
                f();
            );
        };
        auto &&atom_name = read_word(ifile).value;
        read_word(ifile, "in");
        auto &&model_name = read_word(ifile).value;
        read_word(ifile, ":");
        auto f = read_statement(ifile);
        return std::bind(each_atom, atom_name, model_name, f);
    }

    template<typename T>
    std::function<Name()> read_each_res(T &&ifile) {
        auto each_atom = [&](auto &&res_name, auto &&model_name, auto &&f)->Name{
            EACH_RES(Model(model_name),
                this->set_var(res_name, &RES, 1);
                f();
            );
        };
        auto &&res_name = read_word(ifile).value;
        read_word(ifile, "in");
        auto &&model_name = read_word(ifile).value;
        read_word(ifile, ":");
        auto f = read_statement(ifile);
        return std::bind(each_atom, res_name, model_name, f);
    }

    template<typename T>
    std::function<Name()> read_print(T &&ifile) {
        auto print = [&](auto &&name)->Name{
            if (_vars[name].type == 0) {
                auto &&atom = *((Atom *)(_vars[name].addr));
                std::cout << atom.name << ' ' << atom[0] << ' ' << atom[1] << ' ' << atom[2] << std::endl;
            } else if (_vars[name].type == 1) {
                auto &&res = *((Residue *)(_vars[name].addr));
                auto &&atom = res["C4*"];
                std::cout << atom.name << ' ' << atom[0] << ' ' << atom[1] << ' ' << atom[2] << std::endl;
            }
        };
        auto &&var = read_word(ifile).value;
        return std::bind(print, var);
    }

    template<typename T>
    std::function<Name()> read_distance(T &&ifile) {
        auto distance = [&](auto &&name1, auto &&name2)->Var{
            auto &&res1 = *((Residue *)(_vars[name1].addr));
            auto &&res2 = *((Residue *)(_vars[name2].addr));
            auto &&atom1 = res1["C4*"];
            auto &&atom2 = res2["C4*"];
            std::cout << geom::distance(atom1, atom2) << std::endl;
        };
        auto &&var1 = read_word(ifile).value;
        auto &&var2 = read_word(ifile).value;
        return std::bind(distance, var1, var2);
    }

    template<typename T, typename U>
    std::function<Name()> read_assign(T &&ifile, U &&word1) {
        auto assign = [&](auto &&word1, auto &&val2) ->Var {
            _vars[word1.value] = var2;
            return var2;
        }
        auto &&val2 = read_statement(ifile)();
        return std::bind(assign, word1, val2);
    }

    template<typename T, typename U, typename V>
    std::function<Name()> read_exec_func(T &&ifile, U &&func, V &&word) {
        auto model = [&](auto &&var) ->Var {
            if (var.type == VAR_STRING) {
                Model *model = new Model(*((S *)(var.addr)));
                return make_var(model, VAR_MODEL);
            }
        }
        if (func.value == "model") {
            return std::bind(model, var(word));
        }
    }

    template<typename T>
    std::function<Name()> read_statement(T &&ifile) {
        while (ifile) {
            auto &&word = read_word(ifile);
            if (word.value == "each") {
                word = read_word(ifile);
                if (word.value == "atom") {
                    return read_each_atom(ifile);
                } else if (word.value == "res") {
                    return read_each_res(ifile);
                }
            } else if (word.value == "print") {
                return read_print(ifile);
            } else if (word.value == "distance") {
                return read_distance(ifile);
            } else if (word.value == "\n" || word.value == "\r") {
                return read_statement(ifile);
            } else {
                auto &&word2 = read_word(ifile);
                if (word.type == WORD_VARIABLE && word2.type == WORD_VARIABLE) {
                    return read_exec_func(ifile, word, word2);
                } else if (word.type == WORD_VARIABLE && word2.value == "=") {
                    return read_assign(ifile, word);
                }
                std::cout << "Unknown symbol: " << word.value << std::endl;
                throw "Syntax error!";
            }
        }
        throw "End of file.";
    }

    void run(const S &file_name) {
        std::ifstream ifile(file_name.c_str());
        while (ifile) {
            try {
                read_statement(ifile)();
            } catch (const char *inf) {
                std::cout << inf << std::endl;
                break;
            }
        }
        ifile.close();
    }

    S anonym() {
        static int n = 0;
        n++;
        return "anonym_"s + std::to_string(n);
    }

    void exec(const S &content) {}

};

}

