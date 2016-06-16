#pragma once

#include <string>
#include <vector>
#include <utility>
#include <boost/lexical_cast.hpp>

#define JN_INT(n) boost::lexical_cast<int>(n)
#define JN_FLT(n) boost::lexical_cast<float>(n)
#define JN_DBL(n) boost::lexical_cast<double>(n)
#define JN_STR(n) boost::lexical_cast<std::string>(n)
#define JINT JN_INT
#define JSTR JN_STR

namespace jian {

void tokenize(const std::string &str, std::vector<std::string> &tokens, const std::string &delimiters = " ");
void tokenize(const std::string &str, std::vector<std::string> &tokens, const std::string &delimiters, const std::string &temp);
std::string upper(const std::string &str);
std::string lower(const std::string &str);
std::string &trim(std::string &);
std::string trim_copy(const std::string &);

} // namespace jian

