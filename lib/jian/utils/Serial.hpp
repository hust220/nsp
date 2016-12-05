#pragma once

#include <string>
#include <sstream>
#include "string.hpp"
#include "../utils/traits.hpp"

BEGIN_JN

	template<typename T, JN_ENABLE(std::is_integral<T>::value || std::is_floating_point<T>::value || JN_IS_SAME(T, char))>
	inline std::ostream &operator <(std::ostream &stream, const T &t) {
		stream.write(reinterpret_cast<const char *>(&t), sizeof t);
		return stream;
	}

	inline std::ostream &operator <(std::ostream &stream, const S &t) {
		stream.write(t.c_str(), t.size() + 1);
		return stream;
	}

	template<typename T, JN_ENABLE(std::is_integral<T>::value || std::is_floating_point<T>::value || JN_IS_SAME(T, char))>
	inline std::istream &operator >(std::istream &stream, T &t) {
		stream.read(reinterpret_cast<char *>(&t), sizeof t);
		return stream;
	}

	inline std::istream &operator >(std::istream &stream, S &t) {
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

		void stringify_(std::stringstream &stream) {}

		template<typename Head, typename... Tail>
		void stringify_(std::stringstream &stream, Head && head, Tail && ...tail) {
			//write(stream, head);
			stream < head;
			stringify_(stream, tail...);
		}

		template<typename Head, typename... Tail>
		S stringify(Head && head, Tail && ...tail) {
			std::stringstream stream;
			Str s;
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

		//void read(std::stringstream &stream, S &t) {
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
		void parse(const S &s, Head && head, Tail && ...tail) {
			std::stringstream stream;
			stream.write(s.c_str(), s.size());
			parse_(stream, head, tail...);
		}

	};

}