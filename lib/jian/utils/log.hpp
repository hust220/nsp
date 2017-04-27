#pragma once

#include <map>
#include <thread>
#include <fstream>
#include <deque>
#include <mutex>
#include <iostream>
#include <string>
#include "platform.hpp"
#include "string.hpp"
#include "traits.hpp"

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

BEGIN_JN

#ifdef JN_OS_WIN
#  pragma warning(push)
#  pragma warning(disable: 4355) // 'this' : used in base member initializer list
#endif

const int LOG_LEVEL_NONE = 0;
const int LOG_LEVEL_FATAL = 1;
const int LOG_LEVEL_ERROR = 2;
const int LOG_LEVEL_WARNING = 3;
const int LOG_LEVEL_INFO = 4;
const int LOG_LEVEL_DEBUG = 5;
const int LOG_LEVEL_VERBOSE = 6;

template<typename CharType, class CharTraits = STD_ char_traits<CharType> >
class basic_nullbuf : public STD_ basic_streambuf<CharType, CharTraits>
{
	typedef STD_ basic_streambuf<CharType, CharTraits>  base_type;
public:
	// Types
	typedef typename base_type::char_type    char_type;
	typedef typename base_type::traits_type  traits_type;
	typedef typename base_type::int_type     int_type;
	typedef typename base_type::pos_type     pos_type;
	typedef typename base_type::off_type     off_type;

protected:
	virtual  ::std::streamsize  xsputn(char_type const* /*s*/, STD_ streamsize n) { return n; } // "s" is unused
	virtual  int_type           overflow(int_type c = traits_type::eof()) { return traits_type::not_eof(c); }
};

typedef basic_nullbuf<char>      nullbuf;
typedef basic_nullbuf<wchar_t>  wnullbuf;

template< typename _CharType, class _CharTraits = STD_ char_traits<_CharType> >
class basic_onullstream : 
	protected member_from_base<basic_nullbuf<_CharType, _CharTraits>>,
	public STD_ basic_ostream<_CharType, _CharTraits>
{
	typedef member_from_base<basic_nullbuf<_CharType, _CharTraits>> pbase_type;
	typedef STD_ basic_ostream<_CharType, _CharTraits> base_type;
public:
	basic_onullstream() : pbase_type(), base_type(&this->pbase_type::member) {}
};

typedef basic_onullstream<char>      onullstream;
typedef basic_onullstream<wchar_t>  wonullstream;

template<typename _CharType, typename _CharTraits = STD_ char_traits<_CharType>>
class BasicLog : 
	protected member_from_base<basic_nullbuf<_CharType, _CharTraits>>,
	protected member_from_base<STD_ basic_filebuf<_CharType, _CharTraits>>,
	public    STD_ basic_ostream<_CharType, _CharTraits>
{
public:
	using nullbuf = member_from_base<basic_nullbuf<_CharType, _CharTraits>>;
	using filebuf = member_from_base<STD_ basic_filebuf<_CharType, _CharTraits>>;
	using base_type = STD_ basic_ostream<_CharType, _CharTraits>;
	using string_type = BasicStr<_CharType, _CharTraits>;

	basic_onullstream<_CharType> onstream;
	STD_ basic_ofstream<_CharType> ofile;
	string_type m_filename;

	BasicLog() : base_type(&this->nullbuf::member) {}

	BasicLog(string_type filename) : BasicLog() {
		file(filename);
	}

	void file(string_type filename, STD_ ios::openmode openmode = STD_ ios::out) {
		m_filename = filename;
		if (filename == "") {
			this->rdbuf(&this->nullbuf::member);
		}
		else if (filename == "std.out") {
			this->rdbuf(Out.rdbuf());
		}
		else if (filename == "std.err") {
			this->rdbuf(Err.rdbuf());
		}
		else {
			auto &buf = this->filebuf::member;
			if (buf.is_open()) buf.close();
			buf.open(filename, openmode);
			this->rdbuf(&buf);
		}
	}

	string_type file() const {
		return m_filename;
	}

};

using Log = BasicLog<char>;
using Logw = BasicLog<wchar_t>;

// Swallow all types
class Logger
{
public:
	std::map<std::thread::id, std::ostream *> streams;
	std::map<std::thread::id, int> levels;
	std::map<std::thread::id, Str> filenames;
	std::map<std::thread::id, std::deque<Str>> chains;
	int level = LOG_LEVEL_INFO;
	std::mutex mt;
	onullstream onstream;

	Logger(Str file) {
		set_file(file);
	}

	std::thread::id get_id() {
		return std::this_thread::get_id();
	}

	void push() {
		auto id = get_id();
		chains[id].push_back(filenames[id]);
	}

	void pop() {
		auto id = get_id();
		set_file(chains[id].back(), std::ios::app);
		chains[id].pop_back();
	}

	void set_file(Str filename, std::ios::openmode openmode = std::ios::out) {
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
			streams[id] = new std::ofstream(filename.c_str(), openmode);
		}
		filenames[id] = filename;
	}

	Str get_file() {
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
		static Logger g_log("std.out");
		return g_log;
	}

	static Logger &out() {
		static Logger g_out("std.out");
		return g_out;
	}
};

inline void log_push() {
	LOGMANAGER.push();
}

inline void log_pop() {
	LOGMANAGER.pop();
}

inline void out_push() {
	OUTMANAGER.push();
}

inline void out_pop() {
	OUTMANAGER.pop();
}

inline void out_file(Str filename, std::ios::openmode openmode = std::ios::out) {
	OUTMANAGER.set_file(filename, openmode);
}

inline Str out_file() {
	return OUTMANAGER.get_file();
}

inline void log_file(Str file_name, std::ios::openmode openmode = std::ios::out) {
	LOGMANAGER.set_file(file_name, openmode);
}

inline Str log_file() {
	return LOGMANAGER.get_file();
}

inline void log_level(int l) {
	LOGMANAGER.set_level(l);
}

inline int log_level() {
	return LOGMANAGER.get_level();
}

#ifdef JN_OS_WIN
#  pragma warning(default: 4355)
#endif

END_JN

