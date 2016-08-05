#include "log.hpp"
#include <map>
#include <thread>
#include <fstream>
#include <mutex>

namespace jian {

namespace log_detail {

std::mutex mt;

class Impl {
public:
    std::map<std::thread::id, std::ostream *> loggers;
    std::map<std::thread::id, int> levels;
    int level = LOG_LEVEL_INFO;

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

    void set_log_level(int l) {
        std::lock_guard<std::mutex> gd(mt);
        level = l;
    }

    void set_this_log_level(int l) {
        std::lock_guard<std::mutex> gd(mt);
        auto id = get_id();
        levels[id] = l;
    }

    int log_level() {
        std::lock_guard<std::mutex> gd(mt);
        return level;
//        auto id = get_id();
//        if (levels.count(id)) {
//            return levels[id];
//        } else {
//            levels[id] = level;
//            return level;
//        }
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

Impl impl("");

} // log_detail

void log_file(const std::string &file_name) {
    log_detail::impl.log_file(file_name);
}

void set_log_level(int l) {
    log_detail::impl.set_log_level(l);
}

void set_this_log_level(int l) {
    log_detail::impl.set_this_log_level(l);
}

int log_level() {
    return log_detail::impl.log_level();
}

std::ostream *logger() {
    return log_detail::impl.logger();
}

}

