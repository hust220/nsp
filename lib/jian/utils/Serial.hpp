#pragma once

#include <string>
#include <sstream>

namespace jian {

	class Serial {
	public:

		void stringify(std::ostringstream &stream) {}

		template<typename Head, typename... Tail>
		void stringify(std::ostringstream &stream, Head && head, Tail && ...tail) {
			stream << head;
			stringify(stream, tail...);
		}

		template<typename Head, typename... Tail>
		std::string stringify(Head && head, Tail && ...tail) {
			std::ostringstream stream;
			stringify(stream, head, tail...);
			return stream.str();
		}

		void parse(std::istringstream &stream) {}

		template<typename Head, typename... Tail>
		void parse(std::istringstream &stream, Head && head, Tail && ...tail) {
			stream >> head;
			parse(stream, tail...);
		}

		template<typename Head, typename... Tail>
		void parse(const std::string &s, Head && head, Tail && ...tail) {
			std::istringstream stream(s);
			parse(stream, head, tail...);
		}

	};

}