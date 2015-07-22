#include "MLib.h"

namespace jian {

void tokenize(const string &str, vector<string> &tokens, const string &delimiters) {
	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	string::size_type pos     = str.find_first_of(delimiters, lastPos);
	while (string::npos != pos || string::npos != lastPos) {
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
}

string upper(string str) {
	transform(begin(str), end(str), begin(str), ::toupper);
	return str;
}

string lower(string str) {
	transform(begin(str), end(str), begin(str), ::tolower);
	return str;
}

int die(std::string str) {
	std::cout << str << std::endl;
	exit(1);
}

std::string env(std::string str) {
    char *path = getenv(str.c_str());
    std::string env_path = (path ? path : ".");
    return env_path;
}

} /// namespace jian

