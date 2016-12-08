#pragma once

#include "Entity.hpp"

BEGIN_JN

template<typename _Derived, typename _Data>
class BasicListEl :
	public BasicEl<_Derived, _Data>
{
public:
	using El = _Derived;
	using Data = _Data;

	El *next = NULL;
};

template<typename _Data>
class PlainListEl :
	public BasicListEl<PlainListEl<_Data>, _Data>
{};

template<typename _Derived, typename _El>
class BasicListIt :
	public BasicIt<_Derived, _El>
{
public:
	using It = _Derived;
	using El = _El;
	using Data = typename El::Data;

	El *el;

	BasicListIt() :
		el(NULL)
	{}

	virtual Data &operator *() const
	{
		return el->data;
	}

	virtual Bool operator ==(It other) const
	{
		return el == other.el;
	}

	virtual It &operator ++()
	{
		if (el == NULL) throw "ListRangeIt operator ++ error!";
		el = el->next;
		return *(It *)this;
	}

};

template<typename _El>
class PlainListIt :
	public BasicListIt<PlainListIt<_El>, _El>
{};

template<typename _Derived, typename _It>
class BasicListRg :
	public BasicRg<_Derived, _It>
{
public:
	using Rg = _Derived;
	using It = _It;
	using El = typename It::El;
	using Data = typename It::Data;
	using Base = BasicRg<Rg, It>;

	It last() const {
		auto it = Base::m_beg;
		for (; STD_ next(it) != Base::m_end; ++it);
		return it;
	}

	Data &back() const {
		return *last();
	}

protected:
	void free(El *el) {
		if (el != NULL) {
			free(el->next);
			delete el;
		}
	}
	void free() {
		free(Base::m_beg.el);
		Base::m_beg.el = NULL;
	}

};

template<typename _It>
class PlainListRg :
	public BasicListRg<PlainListRg<_It>, _It>
{};

template<typename _Rg>
class BasicListNt :
	public Entity<_Rg>
{
public:
	using Rg = _Rg;
	using Data = typename Rg::Data;
	using El = typename Rg::El;
	using It = typename Rg::It;
	using Nt = BasicListNt<Rg>;

	BasicListNt() = default;

	BasicListNt(const Nt &nt) {
		this->free();
		for (auto && i : nt) {
			push_back(i);
		}
	}

	BasicListNt(Nt &&nt) {
		this->free();
		this->m_beg = nt.m_beg;
		this->m_end = nt.m_end;
		nt.m_beg = It();
		nt.m_end = It();
	}

	BasicListNt &operator =(const Nt &nt) {
		this->free();
		for (auto && i : nt) {
			push_back(i);
		}
		return *this;
	}

	BasicListNt &operator =(Nt &&nt) {
		this->free();
		this->m_beg = nt.m_beg;
		this->m_end = nt.m_end;
		nt.m_beg = It();
		nt.m_end = It();
		return *this;
	}

	template<typename _Data>
	void push_back(_Data &&data) {
		El *el = El::make(std::forward<_Data>(data));
		if (this->empty()) {
			this->m_beg.el = el;
		}
		else {
			this->last().el->next = el;
		}
	}

	template<typename _First, typename _Second, typename... _Rest>
	void push_back(_First &&first, _Second &&second, _Rest &&...rest) {
		push_back(STD_ forward<_First>(first));
		push_back(STD_ forward<_Second>(second), STD_ forward<_Rest>(rest)...);
	}

	template<typename _Data>
	void push_front(_Data &&data) {
		El *el = make(std::forward<_Data>(data));
		el->next = this->m_beg.el;
		this->m_beg.el = el;
	}

	template<typename _First, typename _Second, typename... _Rest>
	void push_front(_First &&first, _Second &&second, _Rest &&...rest) {
		push_front(STD_ forward<_First>(first));
		push_front(STD_ forward<_Second>(second), STD_ forward<_Rest>(rest)...);
	}

};

template<typename _Data>
using SList = BasicListNt<PlainListRg<PlainListIt<PlainListEl<_Data>>>>;

END_JN