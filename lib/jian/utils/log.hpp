#pragma once

#include <iostream>
#include <string>

namespace jian {

#define LOG *(logger())

std::ostream *logger();

void log_file(const std::string &);

}

