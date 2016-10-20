#ifndef JIAN_ETL_MATH
#define JIAN_ETL_MATH

#define PI 3.1415927

namespace jian {

template<typename T> 
inline T square(const T &data) {
    return data * data;
}

template<typename F, typename T1>
inline auto sum(F &&f, T1 &&t1) {
    return f(std::forward<T1>(t1));
}

template<typename F, typename T1, typename... T2>
inline auto sum(F &&f, T1 &&t1, T2 && ...t2) {
    return f(std::forward<T1>(t1)) + sum(f, std::forward<T2>(t2)...);
}

} // namespace jian

#endif





