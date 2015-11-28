#include "MLib.h"

namespace jian {

void tokenize(const string &str, vector<string> &tokens, const string &delimiters) {
    tokens.clear();
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);
    while (string::npos != pos || string::npos != lastPos) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

void tokenize(const string &str, vector<string> &tokens, const string &delimiters, const string &temp) {
    tokens.clear();
    std::vector<std::pair<string::size_type, string::size_type>> vec;
    string::size_type first_i, first_j, second_i, second_j;
    int expected = 0;
    for (string::size_type i = 0; i < str.size(); i++) {
        int flag = 0;
        string::size_type j;
        for (j = 0; j < temp.size(); j++) {
            if (str[i] == temp[j]) {
                if (j % 2 == 0 && expected == 0) {
                    flag = 1;
                    break;
                } else if (j % 2 == 1 && expected == 1) {
                    flag = 2;
                    break;
                }
            }
        }
        if (flag == 1) {
            first_i = i;
            first_j = j;
            expected = 1;
        } else if (flag == 2 && j - first_j == 1) {
            second_i = i;
            second_j = j;
            expected = 0;
            vec.push_back(std::make_pair(first_i, second_i));
        }
    }
//    for (int i = 0; i < temp.size(); i += 2) {
//        string::size_type pos1 = str.find_first_of(string() + temp[i], 0);
//        string::size_type pos2 = str.find_first_of(string() + temp[i + 1], pos1 + 1);
//        while (string::npos != pos1 && string::npos != pos2) {
//            vec.push_back(make_pair(pos1, pos2));
//            pos1 = str.find_first_of(string() + temp[i], pos2 + 1);
//            pos2 = str.find_first_of(string() + temp[i + 1], pos1 + 1);
//        }
//    }
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);
    while (any_of(vec.begin(), vec.end(), [&pos](const pair<string::size_type, string::size_type> &p){
        return pos != string::npos && p.first < pos && pos < p.second;
    })) {
        pos = str.find_first_of(delimiters, pos + 1);
    }
    while (string::npos != pos || string::npos != lastPos) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
        while (any_of(vec.begin(), vec.end(), [&pos](const pair<string::size_type, string::size_type> &p){
            return pos != string::npos && p.first < pos && pos < p.second;
        })) {
            pos = str.find_first_of(delimiters, pos + 1);
        }
    }
}

string upper(string str) {
    transform(begin(str), end(str), begin(str), ::toupper);
    return str;
}

string lower(string str) {
    transform(begin(str), end(str), begin(str), ::tolower);
    return str;
}

int die(std::string str) {
    std::cout << str << std::endl;
    exit(1);
}

std::string env(std::string str) {
    char *path = getenv(str.c_str());
    std::string env_path = (path ? path : throw "Please set the 'NSP' environment variable to the path of the library of nsp!");
    return env_path;
}

} /// namespace jian

