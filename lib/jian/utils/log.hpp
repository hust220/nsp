#pragma once

#include <map>
#include <thread>
#include <fstream>
#include <mutex>
#include <iostream>
#include <string>
#include "platform.hpp"

#define OUTMANAGER (Logger::out())
#define OUTSTREAM (*(OUTMANAGER.get_stream()))
#define JN_OUT  OUTSTREAM

#define LOGMANAGER (Logger::log())
#define LOGLEVEL LOGMANAGER.get_level()
#define LOGSTREAM (*(LOGMANAGER.get_stream()))
#define LOG  if (LOGLEVEL >= LOG_LEVEL_INFO) LOGSTREAM
#define LOGN if (LOGLEVEL >= LOG_LEVEL_NONE) LOGSTREAM
#define LOGF if (LOGLEVEL >= LOG_LEVEL_FATAL) LOGSTREAM
#define LOGE if (LOGLEVEL >= LOG_LEVEL_ERROR) LOGSTREAM
#define LOGW if (LOGLEVEL >= LOG_LEVEL_WARNING) LOGSTREAM
#define LOGI if (LOGLEVEL >= LOG_LEVEL_INFO) LOGSTREAM
#define LOGD if (LOGLEVEL >= LOG_LEVEL_DEBUG) LOGSTREAM
#define LOGV if (LOGLEVEL >= LOG_LEVEL_VERBOSE) LOGSTREAM

namespace jian {

	const int LOG_LEVEL_NONE = 0;
	const int LOG_LEVEL_FATAL = 1;
	const int LOG_LEVEL_ERROR = 2;
	const int LOG_LEVEL_WARNING = 3;
	const int LOG_LEVEL_INFO = 4;
	const int LOG_LEVEL_DEBUG = 5;
	const int LOG_LEVEL_VERBOSE = 6;


	template<typename CharType, class CharTraits = ::std::char_traits<CharType> >
	class basic_nullbuf : public ::std::basic_streambuf<CharType, CharTraits> {
		typedef ::std::basic_streambuf<CharType, CharTraits>  base_type;
	public:
		// Types
		typedef typename base_type::char_type    char_type;
		typedef typename base_type::traits_type  traits_type;
		typedef typename base_type::int_type     int_type;
		typedef typename base_type::pos_type     pos_type;
		typedef typename base_type::off_type     off_type;

	protected:
		virtual  ::std::streamsize  xsputn(char_type const* /*s*/, ::std::streamsize n) { return n; } // "s" is unused
		virtual  int_type           overflow(int_type c = traits_type::eof()) { return traits_type::not_eof(c); }
	};

	typedef basic_nullbuf<char>      nullbuf;
	typedef basic_nullbuf<wchar_t>  wnullbuf;

#ifdef JN_OS_WIN
# pragma warning(push)
# pragma warning(disable: 4355) // 'this' : used in base member initializer list
#endif

	template<typename T>
	struct member_from_base {
		T member;
	};

	template< typename _CharType, class _CharTraits = ::std::char_traits<_CharType> >
	class basic_onullstream : 
		protected member_from_base<basic_nullbuf<_CharType, _CharTraits>>,
		public ::std::basic_ostream<_CharType, _CharTraits>
	{
		typedef member_from_base<basic_nullbuf<_CharType, _CharTraits>> pbase_type;
		typedef ::std::basic_ostream<_CharType, _CharTraits> base_type;
	public:
		basic_onullstream() : pbase_type(), base_type(&this->pbase_type::member) {}
	};

#ifdef JN_OS_WIN
# pragma warning(default: 4355)
#endif

	typedef basic_onullstream<char>      onullstream;
	typedef basic_onullstream<wchar_t>  wonullstream;

	// Swallow all types
	class Logger {
	public:
		std::map<std::thread::id, std::ostream *> streams;
		std::map<std::thread::id, int> levels;
		std::map<std::thread::id, std::string> filenames;
		int level = LOG_LEVEL_INFO;
		std::mutex mt;
		onullstream onstream;

		Logger(std::string file) {
			set_file(file);
		}

		std::thread::id get_id() {
			return std::this_thread::get_id();
		}

		void set_file(std::string filename) {
			std::lock_guard<std::mutex> gd(mt);
			auto id = get_id();
			if (streams.count(id) && streams[id] != &onstream && streams[id] != &std::cout) {
				delete streams[id];
			}
			if (filename == "") {
				streams[id] = &onstream;
			}
			else if (filename == "std.out") {
				streams[id] = &std::cout;
			}
			else {
				streams[id] = new std::ofstream(filename.c_str());
			}
			filenames[id] = filename;
		}

		std::string get_file() {
			auto id = get_id();
			return filenames[id];
		}

		void set_level(int l) {
			std::lock_guard<std::mutex> gd(mt);
			auto id = get_id();
			levels[id] = l;
		}

		void set_all_level(int l) {
			std::lock_guard<std::mutex> gd(mt);
			level = l;
		}

		int get_level() {
			return level;
		}

		std::ostream *get_stream() {
			std::lock_guard<std::mutex> gd(mt);
			auto id = get_id();
			if (streams.find(id) != streams.end()) {
				return streams[id];
			}
			else {
				return &std::cout;
			}
		}

		~Logger() {
			for (auto && pair : streams) {
				if (pair.second != &std::cout && pair.second != &onstream) {
					delete pair.second;
				}
			}
		}

		static Logger &log() {
			static Logger g_log("");
			return g_log;
		}

		static Logger &out() {
			static Logger g_out("");
			return g_out;
		}
	};

	void out_file(std::string f);

	std::string out_file();

	void log_file(std::string f);

	std::string log_file();

	void log_level(int l);

	int log_level();

}

