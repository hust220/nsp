#pragma once

#include <numeric>
#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <list>
#include <type_traits>
#include "rss_loop.hpp"
#include "rss_helix.hpp"
#include "pp.hpp"
#include "entity_tree.hpp"
#include "log.hpp"

BEGIN_JN

struct SSE {
	Loop loop;
	Helix helix;
	Vector<Pair<Int, Int>> hinges;
    SSE *son = NULL, *bro = NULL;

	Pair<int, int> head_tail() const {
		int left, right;
		if (has_helix()) {
			left = helix.front().res1.num - 1;
			right = helix.front().res2.num - 1;
			return{ left, right };
		}
		else if (has_loop()) {
			left = loop.front().num - 1;
			right = loop.back().num - 1;
			return{ left, right };
		}
		else {
			throw "SSE head_tail error";
		}

	}

	bool has(int i) const {
		return STD_ find_if(loop.begin(), loop.end(), [&i](auto &&res) {
			return res.type != '(' && res.type != ')' && res.num == i;
		}) != loop.end() ||
			STD_ find_if(helix.begin(), helix.end(), [&i](auto &&bp) {
			return bp.res1.num == i || bp.res2.num == i;
		}) != helix.end();
	}

	bool has_helix() const {
		return !helix.empty();
	}

	bool has_loop() const {
		return !loop.empty();
	}

	bool has_son() const {
		return !hinges.empty();
	}

	int num_branches() const {
		return STD_ count_if(loop.begin(), loop.end(), [](auto &&res) {return res.type == ')'; }) / 2;
	}

	int num_sons() const {
		return hinges.size();
	}

	bool is_open() const {
		return num_branches() == num_sons();
	}

	bool is_hp() const {
		return !is_open() && num_sons() == 0;
	}

	bool is_il() const {
		return !is_open() && num_sons() == 1;
	}

	bool is_ml() const {
		return !is_open() && num_sons() > 1;
	}

	friend STD_ ostream &operator <<(STD_ ostream &stream, const SSE &sse) {
		stream << "SSE (" << &sse << ") : " << STD_ endl;
		stream << sse.helix << STD_ endl;
		stream << sse.loop;
		return stream;
	}


};

END_JN
