#include "string.hpp"
#include <algorithm>

BEGIN_JN

void tokenize(const Str &str, tokenize_v &tokens, const Str &delimiters) {
    tokens.clear();
    auto lastPos = str.find_first_not_of(delimiters, 0);
    auto pos = str.find_first_of(delimiters, lastPos);
    while (Str::npos != pos || Str::npos != lastPos) {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
    }
}

void tokenize(const Str &str, tokenize_v &tokens, const Str &delimiters, const Str &temp) {
    tokens.clear();
    std::vector<std::pair<Str::size_type, Str::size_type>> vec;
    Str::size_type first_i, first_j, second_i, second_j;
    int expected = 0;
    for (Str::size_type i = 0; i < str.size(); i++) {
        int flag = 0;
        Str::size_type j;
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
        return pos != Str::npos && p.first < pos && pos < p.second;
    })) {
        pos = str.find_first_of(delimiters, pos + 1);
    }
    while (Str::npos != pos || Str::npos != lastPos) {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
        while (std::any_of(vec.begin(), vec.end(), [&pos](const std::pair<Str::size_type, Str::size_type> &p){
            return pos != Str::npos && p.first < pos && pos < p.second;
        })) {
            pos = str.find_first_of(delimiters, pos + 1);
        }
    }
}

Str upper(const Str &str) {
    Str s = str; 
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    return s;
}

Str lower(const Str &str) {
    Str s = str; 
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}

void to_upper(Str &s) {
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
}

Str to_upper_copy(const Str &str) {
    Str s = str; 
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    return s;
}

void to_lower(Str &s) {
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
}

Str to_lower_copy(const Str &str) {
    Str s = str; 
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}

void trim(Str &s) {  
    if (!(s.empty())) {  
        s.erase(0, s.find_first_not_of(" "));  
        s.erase(s.find_last_not_of(" ") + 1);  
    }
}  

Str trim_copy(const Str &s) {
    if (s.empty()) {
        return s;
    }
    Str::size_type beg = s.find_first_not_of(" ");
    Str::size_type end = s.find_last_not_of(" ");
    if (beg == Str::npos || end == Str::npos) {
        return "";
    } else {
        return Str(s.begin() + beg, s.begin() + end + 1);
    }
}

END_JN

