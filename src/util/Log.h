#ifndef LOG_H
#define LOG_H

#include "std.h"

namespace jian {

class Log {
public:
    Log() {}
    Log(std::string log_file): _log_file(log_file) {}

    void print(std::ostream &out, std::string content, int space_nums = 0) {
        for (int i = 0; i < space_nums; i++) { out << ' '; }
        out << content;
    }

    void println(std::ostream &out, std::string content, int space_nums = 0) {
        for (int i = 0; i < space_nums; i++) { out << ' '; }
        out << content << std::endl;
    }

    void operator ()(std::string content, int space_nums = 0) {
        if (_log_file == "") {
            println(std::cout, content, space_nums);
        } else {
            std::ofstream out(_log_file.c_str(), std::ios_base::app);
            out || die("Log::Log error! Couldn't open '" + _log_file + "'!");
            println(out, content, space_nums);
            out.close();
        }
    }

private:
    std::string _log_file;
};



}

#endif

