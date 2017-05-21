#pragma once

#include <iostream>
#include <string>

#define DEBUG_IN Debug::enter_function(__FUNCTION__)
#define DEBUG_OUT Debug::exit_function(__FUNCTION__)

BEGIN_JN

class Debug {
public:
    static void enter_function(const S &func_name) {
        print("\n<", func_name, ">\n");
    }

    static void exit_function(const S &func_name) {
        print("\n</", func_name, ">\n");
    }

    static void print() {}

    template<typename T, typename... U> 
    static void print(T &&data, U && ...datum) {
        std::cout << data;
        print(datum...);
    }

    template<typename... T> 
    static void println(T && ...data) {
        print(data...);
        print("\n");
    }

};

}

