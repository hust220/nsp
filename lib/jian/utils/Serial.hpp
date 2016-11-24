#pragma once

#include <string>
#include <sstream>
#include "../utils/traits.hpp"

namespace jian {

	template<typename T, JN_ENABLE(std::is_integral<T>::value || std::is_floating_point<T>::value || JN_IS_SAME(T, char))>
	inline std::ostream &operator <(std::ostream &stream, const T &t) {
		stream.write(reinterpret_cast<const char *>(&t), sizeof t);
		return stream;
	}

	inline std::ostream &operator <(std::ostream &stream, const std::string &t) {
		stream.write(t.c_str(), t.size() + 1);
		return stream;
	}

	template<typename T, JN_ENABLE(std::is_integral<T>::value || std::is_floating_point<T>::value || JN_IS_SAME(T, char))>
	inline std::istream &operator >(std::istream &stream, T &t) {
		stream.read(reinterpret_cast<char *>(&t), sizeof t);
		return stream;
	}

	inline std::istream &operator >(std::istream &stream, std::string &t) {
		char c[2];
		std::ostringstream ostream;

		while (stream.read(c, 1)) {
			if (c[0] == '\0') break;
			ostream << c[0];
		}
		t = ostream.str();
		return stream;
	}

	class Serial {
	public:

		//template<typename T, JN_ENABLE(std::is_integral<T>::value || std::is_floating_point<T>::value || JN_IS_SAME(T, char))>
		//void write(std::stringstream &stream, const T &t) {
		//	stream.write(reinterpret_cast<const char *>(&t), sizeof t);
		//}

		//void write(std::stringstream &stream, const std::string &t) {
		//	stream.write(t.c_str(), t.size() + 1);
		//}

		void stringify_(std::stringstream &stream) {}

		template<typename Head, typename... Tail>
		void stringify_(std::stringstream &stream, Head && head, Tail && ...tail) {
			//write(stream, head);
			stream < head;
			stringify_(stream, tail...);
		}

		template<typename Head, typename... Tail>
		std::string stringify(Head && head, Tail && ...tail) {
			std::stringstream stream;
			std::string s;
			std::size_t size;

			stringify_(stream, head, tail...);
			size = static_cast<size_t>(stream.tellp());
			s.resize(size);
			stream.read(&s[0], size);
			return s;
		}

		//template<typename T, JN_ENABLE(std::is_integral<T>::value || std::is_floating_point<T>::value || JN_IS_SAME(T, char))>
		//void read(std::stringstream &stream, T &t) {
		//	stream.read(reinterpret_cast<char *>(&t), sizeof t);
		//}

		//void read(std::stringstream &stream, std::string &t) {
		//	char c[2];
		//	std::ostringstream ostream;

		//	while (stream.read(c, 1)) {
		//		if (c[0] == '\0') break;
		//		ostream << c[0];
		//	}
		//	t = ostream.str();
		//}

		void parse_(std::stringstream &stream) {}

		template<typename Head, typename... Tail>
		void parse_(std::stringstream &stream, Head && head, Tail && ...tail) {
			//read(stream, head);
			stream > head;
			parse_(stream, tail...);
		}

		template<typename Head, typename... Tail>
		void parse(const std::string &s, Head && head, Tail && ...tail) {
			std::stringstream stream;
			stream.write(s.c_str(), s.size());
			parse_(stream, head, tail...);
		}

	};

}