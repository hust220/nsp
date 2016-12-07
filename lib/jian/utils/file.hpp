#pragma once

#include "string.hpp"
#include <regex>
#include <fstream>
#include "Entity.hpp"

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

struct FileLine
{
	using Data = FileLine;

	Str line;
	tokenize_v arr;
	int n = -2;
};

class FileLinesIt :
	//protected member_from_base<FileLine>,
	public BasicIt<FileLinesIt, FileLine>
{
public:
	using It = FileLinesIt;
	using Data = FileLine;
	using El = FileLine;
	//using ElMember = member_from_base<El>;

	//El *el;
	SP<El> el;
	STD_ ifstream *stream;
	Str delimiters;

	FileLinesIt() :
		stream(NULL), delimiters(" "), el(STD_ make_shared<El>())
	{
		el->n = -2;
	}

	FileLinesIt(STD_ ifstream *stream_, Str delimiters_ = " ") :
		stream(stream_), delimiters(delimiters_), el(STD_ make_shared<El>())
	{
		if (stream != NULL) {
			el->n = -1;
			next_line();
		}
		else {
			el->n = -2;
		}
	}

	virtual Data &operator *() const {
		return *el;
	}

	virtual bool operator ==(It other) const {
		return el->n == other.el->n;
	}

	void next_line()
	{
		if (stream == NULL || el->n == -2) throw "next_line error!";
		STD_ getline(*stream, el->line);
		if (!stream->eof()) {
			tokenize(el->line, el->arr, delimiters);
			el->n++;
		}
		else {
			el->n = -2;
		}
	}

	It &operator ++()
	{
		next_line();
		return *this;
	}
};

class FileLinesRg :
	public BasicRg<FileLinesRg, FileLinesIt>
{
public:
	using Data = Str;
	using El = FileLine;
	using It = FileLinesIt;
	using Rg = FileLinesRg;

	STD_ ifstream *stream = NULL;
	Str delimiters;

protected:
	void free()
	{
		if (stream != NULL) delete stream;
		stream = NULL;
	}
};

class FileLines :
	public Entity<FileLinesRg>
{
public:
	using Data = Str;
	using El = FileLine;
	using It = FileLinesIt;
	using Rg = FileLinesRg;
	using Nt = FileLines;

	Str filename;

	FileLines() :
		filename("")
	{}

	FileLines(Str filename_, Str delimiters_ = " ")
	{
		bind(filename_, delimiters_);
	}

	FileLines &bind(Str filename_, Str delimiters_ = " ")
	{
		free();
		filename = filename_;
		stream = new STD_ ifstream(filename.c_str());
		delimiters = delimiters_;
		m_beg = It(stream, delimiters);
		return *this;
	}
};

END_JN

