#include "log.hpp"
#include <map>
#include <thread>
#include <fstream>
#include <mutex>

namespace jian {

static std::mutex mt;

class Impl {
public:
    std::map<std::thread::id, std::ostream *> loggers;

    Impl(std::string file) {
        log_file(file);
    }

    std::thread::id get_id() {
        return std::this_thread::get_id();
    }

    void log_file(const std::string &file_name) {
        std::lock_guard<std::mutex> gd(mt);
        auto id = get_id();
//        std::cout << "thead " << id << " log_file" << std::endl;
        if (loggers.count(id) && loggers[id] != &std::cout) {
            delete loggers[id];
        }
        if (file_name == "") {
            loggers[id] = &std::cout;
        } else {
            loggers[id] = new std::ofstream(file_name.c_str());
        }
    }

    std::ostream *logger() {
        std::lock_guard<std::mutex> gd(mt);
        auto id = get_id();
//        std::cout << "thead " << id << " logger" << std::endl;
//        std::cout << "thead " << id << " " << loggers.size() << std::endl;
        if (loggers.find(id) != loggers.end()) {
            return loggers[id];
        } else {
            return &std::cout;
        }
    }

    ~Impl() {
        for (auto && pair : loggers) {
            if (pair.second != &std::cout) {
                delete pair.second;
            }
        }
    }
};

static Impl impl("");

void log_file(const std::string &file_name) {
    impl.log_file(file_name);
}

std::ostream *logger() {
    return impl.logger();
}

}

