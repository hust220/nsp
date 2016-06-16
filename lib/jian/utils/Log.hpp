#pragma once

#include <iostream>
#include <string>
#include <fstream>

namespace jian {

class Log {
public:
    std::string _log_file = "";
    bool _display = true;

    Log() {}

    Log(const std::string &log_file, bool display = true) : _log_file(log_file), _display(display) {}

    void set_display(bool b) {
        _display = b;
    }

    std::string bind(const std::string &log_file = "") {
        _log_file = log_file;
        return _log_file;
    }

    void clear() {
        if (_log_file != "") {
            std::ofstream out(_log_file.c_str());
            out.close();
        }
    }

    void operator ()() {}

    template<typename Data, typename... Datum> 
    void operator ()(Data &&data, Datum && ...datum) {
        if (!_display) return;
        if (_log_file == "") {
            std::cout << data;
        } else {
            std::ofstream out(_log_file.c_str(), std::ios_base::app);
            if (!out) throw "Log::Log error! Couldn't open '" + _log_file + "'!";
            out << data;
            out.close();
        }
        (*this)(datum...);
    }

};

}

