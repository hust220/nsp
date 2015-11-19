#include "Cif.h"

namespace jian {
   
Cif::Cif(std::string file_name) {
    std::ifstream ifile(file_name.c_str());
    read_name(ifile);
    while (read_block(ifile));
    ifile.close();
}

void Cif::read_name(std::ifstream &ifile) {
    _name = read_word(ifile);
    _name = _name.substr(_name.size() - 4, 4);
    read_word(ifile);
}

int Cif::read_block(std::ifstream &ifile) {
    std::string key = read_word(ifile);
    if (key == "") {
        return 0;    
    } else if (key == "#") {
        return 1;
    }
    std::string value;
    if (key == "loop_") {
        read_loop(ifile);
        return 1;
    } else {
        while (true) {
            value = read_word(ifile);
            if (value == "" || value == "#") {
                die("Cif::read_block error!");
            }
            (*this)[key] = value;
            key = read_word(ifile);
            if (key == "#") {
                return 1;
            } else if (key == "") {
                return 0;
            }
        }
    }
}

void Cif::read_loop(std::ifstream &ifile) {
    std::vector<std::string> names;
    std::string word;
    while (true) {
        word = read_word(ifile);
        if (word[0] == '_') {
            names.push_back(word);
        } else {
            while (true) {
                for (int i = 0; i < names.size(); i++) {
                    _loop[names[i]].push_back(word);
                    word = read_word(ifile);
                    if (word == "" || word == "#") {
                        BOOST_ASSERT(i + 1 == names.size() && "Cif::read_loop error!");
                        return;
                    }
                }
            }
        }
    }
}

std::string Cif::read_word(ifstream &ifile) {
    std::string word;
    A(ifile, word);
    return word;
}

void Cif::A(ifstream &ifile, std::string &word) {
    char c = ifile.get();
    switch (char_type(c)) {
        case 1: 
        case 2: A(ifile, word); break;
        case 3: C(ifile, word); break;
        case 4: D(ifile, word); break;
        case 5: E(ifile, word); break;
        case 6: word += c; B(ifile, word); break;
        case 7: return;
    }
}

void Cif::B(ifstream &ifile, std::string &word) {
    char c = ifile.get();
    switch (char_type(c)) {
        case 1: 
        case 2: return;
        case 3: C(ifile, word); break;
        case 4: D(ifile, word); break;
        case 5: E(ifile, word); break;
        case 6: word += c; B(ifile, word); break;
        case 7: return;
    }
}

void Cif::C(ifstream &ifile, std::string &word) {
    char c = ifile.get();
    switch (char_type(c)) {
        case 1: word += " "; C(ifile,word); break;
        case 2: C(ifile, word); break;
        case 4: 
        case 5:
        case 6: word += c; C(ifile, word); break;
        case 3: B(ifile, word); break;
        case 7: return;
    }
}

void Cif::D(ifstream &ifile, std::string &word) {
    char c = ifile.get();
    switch (char_type(c)) {
        case 1: word += " "; D(ifile,word); break;
        case 2: D(ifile, word); break;
        case 3: 
        case 5:
        case 6: word += c; D(ifile, word); break;
        case 4: B(ifile, word); break;
        case 7: return;
    }
}

void Cif::E(ifstream &ifile, std::string &word) {
    char c = ifile.get();
    switch (char_type(c)) {
        case 1: word += " "; E(ifile,word); break;
        case 3: 
        case 4: 
        case 5:
        case 6: word += c; E(ifile, word); break;
        case 2: F(ifile, word); break;
        case 7: return;
    }
}

void Cif::F(ifstream &ifile, std::string &word) {
    char c = ifile.get();
    switch (char_type(c)) {
        case 1: word += " "; E(ifile,word); break;
        case 3: 
        case 4: 
        case 6: word += c; E(ifile, word); break;
        case 2: F(ifile, word); break;
        case 5: return;
        case 7: return;
    }
}

int Cif::char_type(char c) {
    if (c == ' ' || c == '\t') { 
        return 1;
    } else if (c == '\r' || c == '\n') {
        return 2;
    } else if (c == '\'') {
        return 3;
    } else if (c == '"') {
        return 4;
    } else if (c == ';') {
        return 5;
    } else if (c == std::char_traits<char>::eof()) {
        return 7;
    } else {
        return 6;
    }
}

} /// namespace jian

