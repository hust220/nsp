#pragma once

#include <string>
#include <vector>
#include <map>
#include <utility>
#include "lexical_cast.hpp"
#include "traits.hpp"

#define JN_INT(n) jian::lexical_cast<int>(n)
#define JN_FLT(n) jian::lexical_cast<float>(n)
#define JN_DBL(n) jian::lexical_cast<double>(n)
#define JN_NUM(n) jian::lexical_cast<Num>(n)
#define JN_STR(n) jian::lexical_cast<S>(n)
#define JINT JN_INT
#define JSTR JN_STR

BEGIN_JN

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

END_JN

