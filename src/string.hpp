#pragma once

#include <string>
#include <vector>
#include <map>
#include <utility>
#include "lexical_cast.hpp"
#include "jian.hpp"

#define JN_INT(n) jian::lexical_cast<int>(n)
#define JN_FLT(n) jian::lexical_cast<float>(n)
#define JN_DBL(n) jian::lexical_cast<double>(n)
#define JN_NUM(n) jian::lexical_cast<Num>(n)
#define JN_STR(n) jian::lexical_cast<S>(n)
#define JINT JN_INT
#define JSTR JN_STR

namespace jian {

using tokenize_v = std::vector<Str>;
void tokenize(const Str &str, tokenize_v &tokens, const Str &delimiters = " ");
void tokenize(const Str &str, tokenize_v &tokens, const Str &delimiters, const Str &temp);

template<typename... _Args>
inline Str to_str(_Args && ...args) {
    std::ostringstream stream;
    stream_push(stream, args...);
    return stream.str();
}

Str upper(const Str &str);
Str lower(const Str &str);

void to_upper(Str &str);
Str to_upper_copy(const Str &str);

void to_lower(Str &str);
Str to_lower_copy(const Str &str);

void trim(Str &);
Str trim_copy(const Str &);

template<typename _Interval, typename _Ls>
static Str join(_Interval && interval, _Ls && ls) {
    std::stringstream stream;
    Int i = 0;
    for (const auto & s : ls) {
        if (i != 0) stream << interval;
        stream << s;
        i++;
    }
    return stream.str();
}

template<typename... Pars_>
::std::string format(std::string fmt, Pars_ && ...pars) {
    int count = snprintf(NULL, 0, fmt.c_str(), pars...);
    ::std::string buf;
    buf.resize(count);
    sprintf(&(buf[0]), fmt.c_str(), pars...);
    return std::move(buf);
}

inline Vs string_tokenize(const Str &str, const Str &delimiters) {
    Vs tokens;
    auto lastPos = str.find_first_not_of(delimiters, 0);
    auto pos = str.find_first_of(delimiters, lastPos);
    while (Str::npos != pos || Str::npos != lastPos) {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
    }
    return std::move(tokens);
}

inline Vs string_tokenize(const Str &str, const Str &delimiters, const Str &temp) {
    Vs tokens;
    using pair_t = ::std::pair<Str::size_type, Str::size_type>;
    ::std::vector<pair_t> vec;
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
            vec.push_back(::std::make_pair(first_i, second_i));
        }
    }
    auto lastPos = str.find_first_not_of(delimiters, 0);
    auto pos = str.find_first_of(delimiters, lastPos);
    while (::std::any_of(vec.begin(), vec.end(), [&pos](const pair_t &p){
        return pos != Str::npos && p.first < pos && pos < p.second;
    })) {
        pos = str.find_first_of(delimiters, pos + 1);
    }
    while (Str::npos != pos || Str::npos != lastPos) {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
        while (::std::any_of(vec.begin(), vec.end(), [&pos](const pair_t &p){
            return pos != Str::npos && p.first < pos && pos < p.second;
        })) {
            pos = str.find_first_of(delimiters, pos + 1);
        }
    }
    return std::move(tokens);
}

inline void string_trim(Str &s) {
    if (!(s.empty())) {  
        s.erase(0, s.find_first_not_of(" "));  
        s.erase(s.find_last_not_of(" ") + 1);  
    }
}  

inline Str string_trim_c(const Str &s) {
    if (s.empty()) {
        return s;
    }
    Str::size_type beg = s.find_first_not_of(" \t\n\r");
    Str::size_type end = s.find_last_not_of(" \t\n\r");
    if (beg == Str::npos || end == Str::npos) {
        return "";
    } else {
        return Str(s.begin() + beg, s.begin() + end + 1);
    }
}

inline bool string_starts_with(const Str &str1, const Str &str2) {
    int sz = str2.size();
    for (int i = 0; i < sz; i++) {
        if (str1[i] != str2[i]) return false;
    }
    return true;
}

inline bool string_ends_with(const Str &str1, const Str &str2) {
    int sz = str2.size();
    auto it1 = str1.rbegin();
    auto it2 = str2.rbegin();
    for (int i = 0; i < sz; i++) {
        if (*it1 != *it2) return false;
        it1++;
        it2++;
    }
    return true;
}

}

