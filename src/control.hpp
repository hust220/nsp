#ifndef JIAN_ETL_CONTROL
#define JIAN_ETL_CONTROL

BEGIN_JN

namespace static_if_detail {

struct identity {
    template<typename T>
    T operator()(T&& x) const {
        return std::forward<T>(x);
    }
};

template<bool Cond>
struct statement {
    template<typename F>
    void then(const F& f){
        f(identity());
    }

    template<typename F>
    void else_(const F&){}
};

template<>
struct statement<false> {
    template<typename F>
    void then(const F&){}

    template<typename F>
    void else_(const F& f){
        f(identity());
    }
};

} // namespace static_if_detail

template<bool Cond, typename F>
static_if_detail::statement<Cond> static_if(const F &f){
    static_if_detail::statement<Cond> if_;
    if_.then(f);
    return if_;
}

END_JN

#endif


