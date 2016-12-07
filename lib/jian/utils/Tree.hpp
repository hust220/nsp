#include <iterator>
#include "Entity.hpp"
#include "traits.hpp"

BEGIN_JN

template<typename _Data>
using TreePath = L<_Data *>;

template<typename _Derived, typename _Data>
class BasicTreeEl :
	public BasicEl<_Derived, _Data>
{
public:
	using El = _Derived;
	using Data = _Data;
	using Base = BasicEl<El, Data>;

	El *son = NULL;
	El *brother = NULL;

	static El *deep_copy(El *other) {
		if (other == NULL) return NULL;
		El *el = Base::make(other->data);
		el->son = deep_copy(other->son);
		el->brother = deep_copy(other->brother);
		return el;
	}

	bool has_son() const {
		return son != NULL;
	}

	bool has_brother() const {
		return brother != NULL;
	}

	int num_brothers() const {
		return (has_brother() ? brother->num_brothers() + 1 : 0);
	}

	int num_sons() const {
		return (has_son() ? son->num_brothers() + 1 : 0);
	}

	El *set_son(Data data) {
		son = Base::make(data);
		return son;
	}

	El *set_brother(Data data) {
		brother = Base::make(data);
		return brother;
	}

	El *last_brother() const {
		El *it = this;
		for (; it->brother != NULL; it = it->brother);
		return it;
	}

};

template<typename _Data>
class PlainTreeEl :
	public BasicTreeEl<PlainTreeEl<_Data>, _Data>
{};

template<typename _Derived, typename _El>
class BasicTreeIt : 
	public BasicIt<_Derived, _El>
{
public:
	using It = _Derived;
	using El = _El;
	using Data = typename El::Data;

	El *el = NULL;
	TreePath<El> path;
	//L<Rg> path;

	BasicTreeIt() = default;

	virtual Data &operator *() const
	{
		return el->data;
	}

	virtual bool operator ==(It other) const
	{
		return el == other.el;
	}

	virtual It &operator ++() {
		if (el == NULL) throw "BasicTreeIt operator ++ error!";
		path.push_back(el);
		if (el->son != NULL) {
			el = el->son;
		}
		else if (el->brother != NULL) {
			el = el->brother;
			path.pop_back();
		}
		else {
			while (true) {
				path.pop_back();
				if (!path.empty()) {
					if (path.back()->brother == NULL) continue;
					else {
						el = path.back()->brother;
						path.pop_back();
						break;
					}
				}
				else {
					el = NULL;
					break;
				}
			}
		}

		return *(It *)this;
	}

};

template<typename _El>
class PlainTreeIt :
	public BasicTreeIt<PlainTreeIt<_El>, _El>
{};

template<typename _Derived, typename _It>
class BasicTreeRg :
	public BasicRg<_Derived, _It>
{
public:
	using Rg = _Derived;
	using It = _It;
	using El = typename It::El;
	using Data = typename It::Data;

	El *root() const {
		return this->m_beg.el;
	}

	TreePath<El> path() const {
		TreePath<El> path;
		for (auto it = this->begin(); it != this->end(); it++) {
			path.push_back(it.el);
		}
		return path;
	}

protected:
	void free(El *el) {
		if (el != NULL) {
			free(el->son);
			free(el->brother);
			delete el;
		}
	}
	void free() {
		free(this->m_beg.el);
		this->m_beg.el = NULL;
	}

};

template<typename _It>
class PlainTreeRg :
	public BasicTreeRg<PlainTreeRg<_It>, _It>
{};

template<typename _Rg>
class TreeNt :
	public Entity<_Rg>
{
public:
	using Rg = _Rg;
	using It = typename Rg::It;
	using El = typename Rg::El;
	using Data = typename Rg::Data;
	using Nt = TreeNt<Rg>;
	using Base = Entity<Rg>;

	TreeNt() = default;

	TreeNt(const Nt &nt) {
		this->free();
		this->m_beg.el = El::deep_copy(nt.m_beg.el);
	}

	TreeNt(Nt &&nt) {
		this->free();
		this->m_beg = nt.m_beg;
		this->m_end = nt.m_end;
		nt.m_beg = It();
		nt.m_end = It();
	}

	Nt &operator =(const Nt &nt) {
		this->free();
		this->m_beg.el = El::deep_copy(nt.m_beg.el);
		return *this;
	}

	Nt &operator =(Nt &&nt) {
		this->free();
		this->m_beg = nt.m_beg;
		this->m_end = nt.m_end;
		nt.m_beg = It();
		nt.m_end = It();
		return *this;
	}

	template<typename _Data>
	El *set_root(_Data &&data) {
		this->m_beg.el = El::make(STD_ forward<_Data>(data));
		return this->root();
	}
};

template<typename _Data>
using Tree = TreeNt<PlainTreeRg<PlainTreeIt<PlainTreeEl<_Data>>>>;

END_JN