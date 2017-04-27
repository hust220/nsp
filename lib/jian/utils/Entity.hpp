#pragma once

#include "traits.hpp"

BEGIN_JN

template<typename _Derived, typename _Data>
class BasicEl 
{
public:
	using Data = _Data;
	using El = _Derived;
	using Self = BasicEl<El, Data>;

	Data data;

	template<typename _Ty>
	static El *make(_Ty &&data) {
		El *el = new El;
		el->data = STD_ forward<_Ty>(data);
		return el;
	}

};

template<typename _Derived, typename _El>
class BasicIt :
	public STD_ iterator<
		STD_ input_iterator_tag,    // iterator_category
	    typename _El::Data,					    // value_type
		int                         // diffrence type
	>
{
public:
	using El = _El;
	using Data = typename El::Data;
	using It = _Derived;
	using Self = STD_ iterator<
		STD_ input_iterator_tag,
		Data,
		int
	>;

	virtual Data &operator *() const = 0;

	Data *operator ->() const
	{
		return &(*(*this));
	}

	virtual bool operator ==(It other) const = 0;

	bool operator !=(It other) const {
		return !(*this == other);
	}

	virtual It &operator ++() = 0;

	friend It operator ++(It &it, int) {
		It retval = it;
		++it;
		return retval;
	}

};

template<typename _Derived, typename _It>
class BasicRg
{
public:
	using It = _It;
	using Data = typename It::Data;
	using El = typename It::El;
	using Rg = _Derived;

	It m_beg;
	It m_end;

	static Rg make_range(It beg, It end) {
		Rg rg;
		rg.m_beg = beg;
		rg.m_end = end;
		return rg;
	}

	static Rg make_range(El *el) {
		Rg rg;
		rg.m_beg.el = el;
		return rg;
	}

	Rg range() const {
		return make_range(m_beg, m_end);
	}

	It begin() const {
		return m_beg;
	}

	It end() const {
		return m_end;
	}

	Data &front() const {
		return *m_beg;
	}

	bool empty() const {
		return m_beg == m_end;
	}

	int size() const {
		return STD_ distance(m_beg, m_end);
	}

};

template<typename _Rg>
class Entity :
	public _Rg
{
public:
	using Rg = _Rg;
	using It = typename Rg::It;

	Entity() = default;

	~Entity() {
		Rg::free();
	}

};

END_JN

