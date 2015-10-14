#ifndef LOG_H
#define LOG_H

#include "std.h"

namespace jian {

class Log {
public:
	Log() {}
	Log(std::string log_file): _log_file(log_file) {
		_out.open(log_file.c_str(), std::ios_base::app);
		_out || die("Log::Log error! Couldn't open '" + log_file + "'!");
	}
	Log(const Log &log);
	Log &operator =(const Log &log);
	~Log() {
		if (!_log_file.empty()) _out.close();
	}

	void operator ()(std::string content, int space_nums = 0) {
		println(content, space_nums);
	}

	void print(std::string content, int space_nums = 0) {
		std::ostream &out = (_log_file.empty() ? std::cout : _out);
		for (int i = 0; i < space_nums; i++) {
			out << ' ';
		}
		out << content;
	}

	void println(std::string content, int space_nums = 0) {
		std::ostream &out = (_log_file.empty() ? std::cout : _out);
		for (int i = 0; i < space_nums; i++) {
			out << ' ';
		}
		out << content << std::endl;
	}

	void set_log_file(std::string log_file) {
		if (!_log_file.empty()) _out.close();
		_log_file = log_file;
		_out.open(log_file.c_str(), std::ios_base::app);
		_out || die("Log::Log error! Couldn't open '" + log_file + "'!");
	}

private:
	std::ofstream _out;
	std::string _log_file;
};



}

#endif

