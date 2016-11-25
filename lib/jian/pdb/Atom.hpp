#pragma once

#include <string>
#include <array>
#include <map>
#include "../matrix.hpp"

#define JN_DEFAULT_CONSTRUCTORS(type) \
	type() = default;\
	type(const type &) = default;\
	type(type &&) = default;\
	type &operator =(const type &) = default;\
	type &operator =(type &&) = default

namespace jian {

	class Atom : public std::array<double, 3> {
	public:
		std::string name;
		int num;
		double mass;

		JN_DEFAULT_CONSTRUCTORS(Atom);

		Atom(std::string name, double x, double y, double z, int num = -1) {
			init(name, x, y, z, num);
		}

		void init(std::string name, double x, double y, double z, int num = -1);

		void set_name(const std::string &s);
		void set_mass();

	};

} // namespace jian

