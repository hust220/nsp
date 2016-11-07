#pragma once

#include <iostream>
#include <string>

namespace jian {

	extern std::ostream jnout;

	const int LOG_LEVEL_NONE = 0;
	const int LOG_LEVEL_FATAL = 1;
	const int LOG_LEVEL_ERROR = 2;
	const int LOG_LEVEL_WARNING = 3;
	const int LOG_LEVEL_INFO = 4;
	const int LOG_LEVEL_DEBUG = 5;
	const int LOG_LEVEL_VERBOSE = 6;

#define LOG  if (log_level() >= LOG_LEVEL_INFO) *(logger())
#define LOGN if (log_level() >= LOG_LEVEL_NONE) *(logger())
#define LOGF if (log_level() >= LOG_LEVEL_FATAL) *(logger())
#define LOGE if (log_level() >= LOG_LEVEL_ERROR) *(logger())
#define LOGW if (log_level() >= LOG_LEVEL_WARNING) *(logger())
#define LOGI if (log_level() >= LOG_LEVEL_INFO) *(logger())
#define LOGD if (log_level() >= LOG_LEVEL_DEBUG) *(logger())
#define LOGV if (log_level() >= LOG_LEVEL_VERBOSE) *(logger())

	std::ostream *logger();

	void set_log_level(int);

	void set_this_log_level(int);

	int log_level();

	void log_file(const std::string &);

}

