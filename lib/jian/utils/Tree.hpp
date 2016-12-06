#include <iterator>
#include "Range.hpp"
#include "traits.hpp"

BEGIN_JN

template<typename Val>
using TreePath = L<Val *>;

template<typename _Node>
class PreTreeIt : 
	public RangeIt<PreTreeIt<_Node>, _Node>
{
public:
	using Val = _Node;
	using Self = PreTreeIt<Val>;
	using Base = RangeIt<Self, Val>;

	Val *val;
	TreePath<Val> path;

	PreTreeIt(Val *v = NULL) : val(v) {}

	Self &operator ++() {
		if (val->son != NULL) {
			val = val->son;
		}
		else if (val->brother != NULL) {
			val = val->brother;
			path.pop_back();
		}
		else {
			while (true) {
				path.pop_back();
				if (!path.empty()) {
					val = path.back()->brother;
					if (val == NULL) continue;
					else {
						path.pop_back();
						break;
					}
				}
				else {
					val = NULL;
					break;
				}
			}
		}

		if (val == NULL) return *this;
		path.push_back(val);
		return *this;
	}

};

template<typename _Node>
class TreeRange :
	public Range<PreTreeIt<_Node>>
{
public:
	using Node   = _Node;
	using PreIt = PreTreeIt<Node>;

	Node *head = NULL;

	TreeRange(Node *node = NULL) : head(node) {}

	virtual PreIt begin() const {
		PreIt it(head);
		if (head != NULL) {
			it.path.push_back(head);
		}
		return it;
	}

	PreIt pre_begin() const {
		return begin();
	}

	virtual PreIt end() const {
		return PreIt(NULL);
	}

	PreIt pre_end() const {
		return end();
	}

	TreePath<Node> path() const {
		TreePath<Node> path;
		for (auto && node : *this) {
			path.push_back(&node);
		}
		return path;
	}

};

template<typename _Range>
class Tree :
	public Entity<_Range>
{
public:
	using Node = typename _Range::Node;
	using Range = _Range;

	Tree() = default;

	~Tree() {
		free(Range::head);
	}

	static void free(Node *node) {
		if (node != NULL) {
			free(node->son);
			free(node->brother);
			delete node;
		}
	}
};

template<typename _Node>
TreeRange<_Node> pretree(_Node *node) {
	return TreeRange<_Node>{node};
}

END_JN