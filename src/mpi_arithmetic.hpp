#ifndef JIAN_MPL_ARITHMETIC
#define JIAN_MPL_ARITHMETIC

BEGIN_JN
namespace mpl {

template<int N>
struct I {
    enum {value = N};
};

template<typename A, typename B>
struct plus;

template<int a, int b>
struct plus<I<a>, I<b>> {
    typedef I<a + b> type;
};

template<typename A, typename B>
struct minus;

template<int a, int b>
struct minus<I<a>, I<b>> {
    typedef I<a - b> type;
};

template<typename A, typename B>
struct multiplies;

template<int a, int b>
struct multiplies<I<a>, I<b>> {
    typedef I<a * b> type;
};

template<typename A, typename B>
struct divides;

template<int a, int b>
struct divides<I<a>, I<b>> {
    typedef I<a / b> type;
};

} // namespace mpl
END_JN










#endif

