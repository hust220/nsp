#include "string.hpp"
#include <algorithm>

namespace jian {

void tokenize(const std::string &str, std::vector<std::string> &tokens, const std::string &delimiters) {
    tokens.clear();
    auto lastPos = str.find_first_not_of(delimiters, 0);
    auto pos = str.find_first_of(delimiters, lastPos);
    while (std::string::npos != pos || std::string::npos != lastPos) {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
    }
}

void tokenize(const std::string &str, std::vector<std::string> &tokens, const std::string &delimiters, const std::string &temp) {
    tokens.clear();
    std::vector<std::pair<std::string::size_type, std::string::size_type>> vec;
    std::string::size_type first_i, first_j, second_i, second_j;
    int expected = 0;
    for (std::string::size_type i = 0; i < str.size(); i++) {
        int flag = 0;
        std::string::size_type j;
        for (j = 0; j < temp.size(); j++) {
            if (str[i] == temp[j]) {
                if (j % 2 == 0 && expected == 0) { flag = 1; break;
                } else if (j % 2 == 1 && expected == 1) { flag = 2; break; }
            }
        }
        if (flag == 1) {
            first_i = i; first_j = j; expected = 1;
        } else if (flag == 2 && j - first_j == 1) {
            second_i = i; second_j = j; expected = 0;
            vec.push_back(std::make_pair(first_i, second_i));
        }
    }
    auto lastPos = str.find_first_not_of(delimiters, 0);
    auto pos = str.find_first_of(delimiters, lastPos);
    while (std::any_of(vec.begin(), vec.end(), [&pos](const auto &p){
        return pos != std::string::npos && p.first < pos && pos < p.second;
    })) {
        pos = str.find_first_of(delimiters, pos + 1);
    }
    while (std::string::npos != pos || std::string::npos != lastPos) {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
        while (std::any_of(vec.begin(), vec.end(), [&pos](const std::pair<std::string::size_type, std::string::size_type> &p){
            return pos != std::string::npos && p.first < pos && pos < p.second;
        })) {
            pos = str.find_first_of(delimiters, pos + 1);
        }
    }
}

std::string upper(const std::string &str) {
    std::string s = str; 
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    return s;
}

std::string lower(const std::string &str) {
    std::string s = str; 
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}

std::string& trim(std::string &s) {  
    if (s.empty()) {  
        return s;  
    }
    s.erase(0, s.find_first_not_of(" "));  
    s.erase(s.find_last_not_of(" ") + 1);  
    return s;  
}  

std::string trim_copy(const std::string &s) {
    if (s.empty()) {
        return s;
    }
    std::string::size_type beg = s.find_first_not_of(" ");
    std::string::size_type end = s.find_last_not_of(" ");
    if (beg == std::string::npos || end == std::string::npos) {
        return "";
    } else {
        return std::string(s.begin() + beg, s.begin() + end + 1);
    }
}

} // namespace jian

