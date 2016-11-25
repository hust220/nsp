#pragma once

#include <string>
#include <vector>
#include <utility>
#include "lexical_cast.hpp"
#include "traits.hpp"

#define JN_INT(n) jian::lexical_cast<int>(n)
#define JN_FLT(n) jian::lexical_cast<float>(n)
#define JN_DBL(n) jian::lexical_cast<double>(n)
#define JN_STR(n) jian::lexical_cast<std::string>(n)
#define JINT JN_INT
#define JSTR JN_STR

namespace jian {

	using str_t = std::string;

	using tokenize_v = std::vector<str_t>;
	void tokenize(const str_t &str, tokenize_v &tokens, const str_t &delimiters = " ");
	void tokenize(const str_t &str, tokenize_v &tokens, const str_t &delimiters, const str_t &temp);

	str_t upper(const str_t &str);
	str_t lower(const str_t &str);

	void to_upper(str_t &str);
	str_t to_upper_copy(const str_t &str);

	void to_lower(str_t &str);
	str_t to_lower_copy(const str_t &str);

	void trim(str_t &);
	str_t trim_copy(const str_t &);

} // namespace jian

