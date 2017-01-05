#pragma once

#include <type_traits>
#include <functional>
#include <vector>
#include <list>
#include <deque>
#include <map>
#include <utility>
#include <memory>
#include <set>
#include <tuple>
#include <fstream>
#include <sstream>
#include <iostream>
#include "platform.hpp"

#ifdef JN_OS_WIN
#  pragma warning(disable: 4018)
#endif

#define STD_ ::std::
#define JN_  ::jian::

#define BEGIN_JN namespace jian {
#define END_JN   }

BEGIN_JN

inline void stream_push(std::ostream &stream) {}

template<typename _First, typename... _Tail>
inline void stream_push(std::ostream &stream, _First &&first, _Tail && ...tail) {
	stream << first;
	stream_push(stream, tail...);
}

template<typename... _Args>
void die(_Args &&...args) {
	stream_push(STD_ cerr, args...);
	exit(1);
}

#ifdef JN_PRECISION
using Num = JN_PRECISION;
#else
using Num = double;
using val_t = double;
#endif

using N = Num;

using Int  = int;
using I = Int;
using Double  = double;
using D = Double;
using Float  = float;
using F = Float;
using Bool = bool;
using B = Bool;
using Char = char;
using C = Char;
using WChar = wchar_t;
using WC = WChar;
using W = WC;

template<typename _CharType, typename _CharTraits = STD_ char_traits<_CharType>>
using BasicStr = STD_ basic_string<_CharType, _CharTraits>;
using Str = STD_ string;
using S = Str;
using WStr = STD_ wstring;
using WS = WStr;

template<typename _Type>
using Ptr = _Type *;
template<typename _Type>
using P = Ptr<_Type>;

template<typename _Type>
using SPtr = STD_ shared_ptr<_Type>;
template<typename _Type>
using SP = SPtr<_Type>;

template<typename _Type>
using UPtr = STD_ unique_ptr<_Type>;
template<typename _Type>
using UP = UPtr<_Type>;

template<typename _Fty>
using Fn = STD_ function<_Fty>;
template<typename _Fty>
using Function = STD_ function<_Fty>;

template<typename _Type, I _N>
using A = STD_ array<_Type, _N>;
template<typename _Type, I _N>
using Array = STD_ array<_Type, _N>;
template<int _N>
using An = A<N, _N>;
template<int _N>
using Ai = A<I, _N>;
template<int _N>
using Ad = A<D, _N>;
template<int _N>
using Af = A<F, _N>;
template<int _N>
using Ac = A<C, _N>;
template<int _N>
using As = A<S, _N>;
template<int _N>
using Ab = A<B, _N>;

template<typename _Type>
using V = STD_ vector<_Type>;
template<typename _Type>
using Vector = STD_ vector<_Type>;
using Vn = V<N>;
using Vi = V<I>;
using Vd = V<D>;
using Vf = V<F>;
using Vc = V<C>;
using Vs = V<S>;
using Vb = V<B>;

template<typename _Type>
using L = STD_ list<_Type>;
template<typename _Type>
using List = STD_ list<_Type>;
using Ln = L<N>;
using Li = L<I>;
using Ld = L<D>;
using Lf = L<F>;
using Lc = L<C>;
using Ls = L<S>;
using Lb = L<B>;

template<typename _Type>
using Q = STD_ deque<_Type>;
template<typename _Type>
using Deque = STD_ deque<_Type>;
using Qn = Q<N>;
using Qi = Q<I>;
using Qd = Q<D>;
using Qf = Q<F>;
using Qc = Q<C>;
using Qs = Q<S>;
using Qb = Q<B>;

template<typename _Type>
using T = STD_ set<_Type>;
template<typename _Type>
using Set = STD_ set<_Type>;
using Tn = T<N>;
using Ti = T<I>;
using Td = T<D>;
using Tf = T<F>;
using Tc = T<C>;
using Ts = T<S>;
using Tb = T<B>;

template<typename _Type1, typename _Type2>
using Pair = STD_ pair<_Type1, _Type2>;

template<typename... _Types>
using Tuple = STD_ tuple<_Types...>;

template<typename _KeyType, typename _Type>
using M = STD_ map<_KeyType, _Type>;
template<typename _KeyType, typename _Type>
using Map = STD_ map<_KeyType, _Type>;
template<typename _Type>
using Mn = M<N, _Type>;
template<typename _Type>
using Mi = M<I, _Type>;
template<typename _Type>
using Md = M<D, _Type>;
template<typename _Type>
using Mf = M<F, _Type>;
template<typename _Type>
using Mc = M<C, _Type>;
template<typename _Type>
using Ms = M<S, _Type>;
template<typename _Type>
using Mb = M<B, _Type>;

#define Out  STD_ cout
#define In   STD_ cin 
#define Err  STD_ cerr
#define Endl STD_ endl

#define JN_ENABLE(t) std::enable_if_t<t, Int> = 42
#define JN_IS_SAME(A, B) std::is_same<A, B>::value

#define JN_DEFAULT_CONSTRUCTORS(type) \
	type() = default;\
	type(const type &) = default;\
	type(type &&) = default;\
	type &operator =(const type &) = default;\
	type &operator =(type &&) = default

template<typename T>
inline int size(T &&t) {
	return int(t.size());
}

template<typename T>
struct member_from_base
{
	using type = T;
	T member;
};

template<typename _ValType>
class refs : public std::vector<std::reference_wrapper<_ValType>> {
public:
	using value_type = _ValType;

	_ValType &operator [](int i) const {
		return std::vector<std::reference_wrapper<_ValType>>::operator[](i);
	}

	_ValType &at(int i) const {
		return std::vector<std::reference_wrapper<_ValType>>::at(i);
	}

	refs<_ValType> &append(_ValType &val) {
		std::vector<std::reference_wrapper<_ValType>>::push_back(std::ref(val));
		return *this;
	}

	template<typename _Ls>
	refs<_ValType> &append(_Ls &&ls) {
		for (auto && item : ls) {
			append(item);
		}
		return *this;
	}

};

template<typename T, typename U>
struct is_decay_same {
	enum { value = std::is_same<std::decay_t<T>, std::decay_t<U>>::value };
};

template<typename... T>
struct is_all_decay_same;

template<typename T>
struct is_all_decay_same<T> {
	enum { value = true };
};

template<typename T, typename U, typename... F>
struct is_all_decay_same<T, U, F...> {
	enum { value = is_decay_same<T, U>::value && is_all_decay_same<U, F...>::value };
};

template<typename T, template<typename...> class U>
struct is_template_same {
private:
	template<template<typename...> class F, typename K, typename... Pars>
	static K check(F<K, Pars...>);
public:
	enum { value = std::is_same<U<decltype(check(std::declval<std::decay_t<T>>()))>, std::decay_t<T>>::value };
};

template<typename T, typename U> struct uniform_const { using type = std::remove_const_t<T>; };
template<typename T, typename U> struct uniform_const<T, const U> { using type = const std::remove_const_t<T>; };
template<typename T, typename U> struct uniform_const<T, const U &> { using type = const std::remove_const_t<T>; };
template<typename T, typename U> struct uniform_const<T, const U &&> { using type = const std::remove_const_t<T>; };
template<typename T, typename U> using uniform_const_t = typename uniform_const<T, U>::type;

template<typename T>
struct template_first_parameter {
private:
	template<typename K> struct check { using type = std::false_type; };
	template<template<typename...> class F, typename K, typename... Pars> struct check<F<K, Pars...>> { using type = K; };
	template<template<typename, int...> class F, typename K, int... ints> struct check<F<K, ints...>> { using type = K; };
	template<typename K, typename U> struct modify { using type = K; };
	template<typename K, typename U> struct modify<K, const U> { using type = const K; };
	template<typename K, typename U> struct modify<K, const U &> { using type = const K &; };
	template<typename K, typename U> struct modify<K, const U &&> { using type = const K &&; };
	template<typename K, typename U> struct modify<K, U&> { using type = K&; };
	template<typename K, typename U> struct modify<K, U&&> { using type = K&&; };
public:
	using type = typename modify<typename check<std::decay_t<T>>::type, T>::type;
};

template<typename T>
using template_first_parameter_t = typename template_first_parameter<T>::type;

template<typename _ValType>
inline int count(const _ValType &val) {
	return 1;
}

template<typename _ValType, typename... _Args, template<typename...> class _Ls>
inline int count(const _Ls<_ValType, _Args...> &ls) {
	return int(ls.size());
}

template<typename _ValType, typename _Ls, JN_ENABLE(JN_IS_SAME(_Ls, _ValType))>
inline int count(const _Ls &ls) {
	return 1;
}

template<typename _ValType, typename _Ls, JN_ENABLE(JN_IS_SAME(typename _Ls::value_type, _ValType))>
inline int count(const _Ls &ls) {
	return size(ls);
}

template<typename _ValType, typename _Ls, JN_ENABLE(!JN_IS_SAME(typename _Ls::value_type, _ValType))>
inline int count(const _Ls &ls) {
	int c = 0;
	for (auto && item : ls) c += count<_ValType>(item);
	return c;
}

template<typename _ValType, typename _Fn>
inline bool each(_ValType &&ls, _Fn &&fn) {
	fn(ls);
	return true;
}

template<typename _ValType, typename... _Args, template<typename...> class _Ls, typename _Fn>
inline bool each(_Ls<_ValType, _Args...> &&ls, _Fn &&fn) {
	for (auto && item : ls) {
		if (!fn(item)) return false;
	}
	return true;
}

template<typename _ValType, typename _Ls, typename _Fn, JN_ENABLE(JN_IS_SAME(_Ls, _ValType))>
inline bool each(_Ls &&ls, _Fn &&fn) {
	fn(ls);
	return true;
}

template<typename _ValType, typename _Ls, typename _Fn, JN_ENABLE(JN_IS_SAME(typename std::decay_t<_Ls>::value_type, _ValType))>
inline bool each(_Ls &&ls, _Fn &&fn) {
	for (auto && item : ls) {
		if (!fn(item)) return false;
	}
	return true;
}

template<typename _ValType, typename _Ls, typename _Fn, JN_ENABLE(!JN_IS_SAME(typename std::decay_t<_Ls>::value_type, _ValType))>
inline bool each(_Ls &&ls, _Fn &&fn) {
	for (auto && item : ls) {
		if (!each<_ValType>(item, fn)) return false;
	}
	return true;
}

END_JN

