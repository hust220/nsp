#ifndef JIAN_MPL_LIST
#define JIAN_MPL_LIST

BEGIN_JN
namespace mpl {

#define Define(A, n) struct A {enum {data = n};}

template<typename T, typename U>
struct Pair {
    typedef T head;
    typedef U tail;
};

struct NullType;

template<typename T>
struct car {
    typedef typename T::head type;
};

template<typename T>
struct cdr {
    typedef typename T::tail type;
};

template<typename T>
struct cadr {
    typedef typename T::tail::head type;
};

template<typename T>
struct caddr {
    typedef typename T::tail::tail::head type;
};

template<typename... U>
struct list;

template<>
struct list<> {
    typedef NullType type;
};

template<typename T, typename... U>
struct list<T, U...> {
    typedef Pair<T, typename list<U...>::type> type;
};

} // namespace mpl
END_JN










#endif

