#pragma once

#include <iostream>
#include <string>

#define TRACE_IN Trace::enter_function(__FUNCTION__)
#define TRACE_OUT Trace::exit_function(__FUNCTION__)

BEGIN_JN

class Trace {
public:
    static void enter_function(const S &func_name) {
        log("\n<", func_name, ">\n");
    }

    static void exit_function(const S &func_name) {
        log("\n</", func_name, ">\n");
    }

    static void log() {}

    template<typename T, typename... U> 
    static void log(T &&data, U && ...datum) {
        std::cout << data;
        log(datum...);
    }

};

}

