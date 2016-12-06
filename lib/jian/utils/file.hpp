#pragma once

#include "string.hpp"
#include <regex>
#include <fstream>
#include "Range.hpp"

BEGIN_JN

struct file {
	static S name(const S & file_path);

	static S type(const S & file_path);

	static void clean(const S & file_name);

	template<typename T>
	static void open(T & stream, const char * file_name) {
		stream.open(file_name);
	}

	template<typename T>
	static void open(T & stream, const S & file_name) {
		stream.open(file_name.c_str());
	}
};

#define FOPEN(fstream, filename) \
	jian::file::open(fstream, filename); \
    if (!(fstream)) {\
        std::ostringstream stream;\
        stream << "Can't open file '" << filename << "'!" << std::endl;\
        throw stream.str();\
    }

#define FCLOSE(fstream) (fstream).close()

#define BEGIN_READ_FILE(f, t) do {\
	int N = 0;\
	std::ifstream ifile;\
	FOPEN(ifile, f);\
	S L;\
	tokenize_v F;\
	while (std::getline(ifile, L)) {\
		::jian::tokenize(L, F, t); \
		do

#define END_READ_FILE while(0);N++;}FCLOSE(ifile);}while(0)

struct FileLine {
	Str line;
	tokenize_v arr;
	int n;
};

class FileRangeIt :
	public RangeIt<FileRangeIt, FileLine>
{
public:
	using Val = FileLine;
	using Self = FileRangeIt;
	using Base = RangeIt<Self, Val>;

	STD_ ifstream *stream;
	Val file_line;
	Str delimiters = " ";
	Val *val;

	FileRangeIt() 
	{
		val = NULL;
	}

	void next_line() {
		if (STD_ getline(*stream, file_line.line)) {
			tokenize(file_line.line, file_line.arr, delimiters);
			file_line.n++;
			val = &file_line;
		}
		else {
			val = NULL;
		}
	}

	Self &operator ++() {
		if (val == NULL) throw "ListRangeIt operator ++ error!";
		next_line();
		return *this;
	}
};

class FileRange :
	public Range<FileRangeIt>
{
public:
	using It = FileRangeIt;

	STD_ ifstream *stream;
	Str delimiters;

	FileRange(STD_ ifstream &stream_, Str delimiters_ = " ") :
		stream(&stream_), delimiters(delimiters_)
	{}

	virtual It begin() const {
		It it;
		it.stream = stream;
		it.delimiters = delimiters;
		it.next_line();
		return it;
	}

	virtual It end() const {
		return It{};
	}

};

class File :
	public Entity<FileRange>
{
public:
	Str filename;
	Str delimiters;
	STD_ ifstream stream;

	File() = delete;

	File(Str filename_, Str delimiters_) :
		filename(filename_), delimiters(delimiters_), FileRange(&stream, delimiters)
	{
		stream.open(filename);
	}

	~File() {

	}
};



END_JN

