#ifndef JIAN_FPL_PAIR
#define JIAN_FPL_PAIR

#include "../util/std.h"

namespace jian {
namespace fpl {

template<typename Head, typename Tail>
struct Pair {
    Head head;
    Tail tail;
};

template<typename Val>
class List {
public:
    Val head;
    List *tail;

    List(const List<Val> &ls) {
        head = ls.head;
        tail = new List<Val>(*ls.tail);
    }

    List(List<Val> &&ls) {
        std::swap(head, ls.head);
        std::swap(tail, ls.tail);
    }

    List(const Val &val, List<Val> *ls) {
        head = val;
        tail = ls;
    }

    List(const Val &val, List<Val> &&ls) {
        head = val;
        tail = new List<Val>(std::forward<List<Val>>(ls));
    }

    ~List() {
        if (tail != nullptr) delete tail;
    }
};

template<typename Par>
List<Par> list(const Par &par) {
    return List<Par>(par, nullptr);
}

template<typename Par, typename... Pars>
List<Par> list(const Par &par, const Par &par2, const Pars & ...pars) {
    return List<Par>(par, list<Par>(par2, pars...));
}

template<typename Fn, typename Par>
void each2(Fn &&f, const List<Par> &ls) {
    f(ls.head);
    if (ls.tail != nullptr) each2(f, *ls.tail);
}

template<typename Val>
Val &car(List<Val> &ls) {
    return ls.head;
}

template<typename Val>
const Val &car(const List<Val> &ls) {
    return ls.head;
}

template<typename Val>
List<Val> &cdr(List<Val> &ls) {
    return *(ls.tail);
}

template<typename Val>
const List<Val> &cdr(const List<Val> &ls) {
    return *(ls.tail);
}

template<typename Fn, typename Val, template<typename...> class... Lists>
List<Val> map2(Fn &&f, const Lists<Val> & ...lists) {
    return List<Val>(f(car(lists)...), map2(f, cdr(lists)...));
}

} // namespace fpl
} // namespace jian

#endif





