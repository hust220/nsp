#pragma once

#include "Range.hpp"

BEGIN_JN

template<typename _Node>
class ListRangeIt :
	public RangeIt<ListRangeIt<_Node>, _Node>
{
public:
	using Val = _Node;
	using Self = ListRangeIt<Val>;
	using Base = RangeIt<Self, Val>;

	Val *val;

	ListRangeIt(Val *v = NULL) : val(v) {}

	Self &operator ++() {
		if (val == NULL) throw "ListRangeIt operator ++ error!";
		val = val->next;
		return *this;
	}
};

template<typename _Node>
class ListRange :
	public Range<ListRangeIt<_Node>>
{
public:
	using Self = ListRange<_Node>;
	using Base = Range<ListRange<_Node>>;
	using Node = _Node;
	using It = ListRangeIt<_Node>;

	Node *head;

	ListRange(Node *node = NULL) : head(node) {}

	virtual It begin() const {
		return It(head);
	}

	virtual It end() const {
		return It(NULL);
	}

};

template<typename _Range>
class ListEntity : 
	public Entity<_Range>
{
public:
	using Node = typename _Range::Node;
	using Range = _Range;

	~ListEntity() {
		free(Range::head);
	}

	static void free(Node *node) {
		if (node != NULL) {
			free(node->next);
			delete node;
		}
	}

	template<typename _Nd>
	void push_back(_Nd &&node) {
		Node *n = new Node(std::forward<_Nd>(node));
		if (head == NULL) {
			head = n;
			head->next = NULL;
		}
		else {
			Node *it = head;
			for (; it->next != NULL; it = it->next);
			it->next = n;
			n->next = NULL;
		}
	}

	template<typename _First, typename _Second, typename... _Rest>
	void push_back(_First &&first, _Second &&second, _Rest &&...rest) {
		push_back(std::forward<_First>(first));
		push_back(std::forward<_Rest>(second), std::forward<_Rest>(rest)...)
	}

	template<typename _Nd>
	void push_front(_Nd &&node) {
		Node *n = new Node(std::forward<_Nd>(node));
		n->next = head;
		head = n->next;
	}

	template<typename _First, typename _Second, typename... _Rest>
	void push_front(_First &&first, _Second &&second, _Rest &&...rest) {
		push_front(std::forward<_First>(first));
		push_front(std::forward<_Rest>(second), std::forward<_Rest>(rest)...)
	}
};

END_JN