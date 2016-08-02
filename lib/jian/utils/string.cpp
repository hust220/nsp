#include "string.hpp"
#include <algorithm>

namespace jian {

void tokenize(const str_t &str, tokenize_v_t &tokens, const str_t &delimiters) {
    tokens.clear();
    auto lastPos = str.find_first_not_of(delimiters, 0);
    auto pos = str.find_first_of(delimiters, lastPos);
    while (str_t::npos != pos || str_t::npos != lastPos) {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
    }
}

void tokenize(const str_t &str, tokenize_v_t &tokens, const str_t &delimiters, const str_t &temp) {
    tokens.clear();
    std::vector<std::pair<str_t::size_type, str_t::size_type>> vec;
    str_t::size_type first_i, first_j, second_i, second_j;
    int expected = 0;
    for (str_t::size_type i = 0; i < str.size(); i++) {
        int flag = 0;
        str_t::size_type j;
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
        return pos != str_t::npos && p.first < pos && pos < p.second;
    })) {
        pos = str.find_first_of(delimiters, pos + 1);
    }
    while (str_t::npos != pos || str_t::npos != lastPos) {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
        while (std::any_of(vec.begin(), vec.end(), [&pos](const std::pair<str_t::size_type, str_t::size_type> &p){
            return pos != str_t::npos && p.first < pos && pos < p.second;
        })) {
            pos = str.find_first_of(delimiters, pos + 1);
        }
    }
}

str_t upper(const str_t &str) {
    str_t s = str; 
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    return s;
}

str_t lower(const str_t &str) {
    str_t s = str; 
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}

str_t& trim(str_t &s) {  
    if (s.empty()) {  
        return s;  
    }
    s.erase(0, s.find_first_not_of(" "));  
    s.erase(s.find_last_not_of(" ") + 1);  
    return s;  
}  

str_t trim_copy(const str_t &s) {
    if (s.empty()) {
        return s;
    }
    str_t::size_type beg = s.find_first_not_of(" ");
    str_t::size_type end = s.find_last_not_of(" ");
    if (beg == str_t::npos || end == str_t::npos) {
        return "";
    } else {
        return str_t(s.begin() + beg, s.begin() + end + 1);
    }
}

} // namespace jian

