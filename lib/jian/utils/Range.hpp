#pragma once

#include "traits.hpp"

BEGIN_JN

template<typename _Derived, typename _Node>
class RangeIt :
	public STD_ iterator<
		STD_ input_iterator_tag,    // iterator_category
	    _Node,					    // value_type
		int                         // diffrence type
	>
{
public:
	using Derived = _Derived;
	using Val = _Node;
	using Self = RangeIt<Derived, Val>;
	using Base = STD_ iterator<STD_ input_iterator_tag, Val, int>;

	Val *val;

	Val &operator *() const {
		return *val;
	}

	Val &operator ->() const {
		return *val;
	}

	bool operator ==(Derived other) const {
		return val == other.val;
	}

	bool operator !=(Derived other) const {
		return !(*this == other);
	}

	Derived operator++(int) {
		Derived retval = *this;
		++(*this);
		return retval;
	}

};

template</*typename _Derived, */typename _It>
class Range
{
public:
	//using Derived = _Derived;

	virtual _It begin() const = 0;

	virtual _It end() const = 0;

	bool empty() const {
		return begin() == end();
	}

	int size() const {
		return STD_ distance(begin(), end());
	}

};

template<typename _Range>
class Entity :
	public _Range
{
public:
	Entity() = default;

	_Range range() const {
		return _Range(*this);
	}
};

END_JN

