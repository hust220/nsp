//template<typename T>
//struct template_template_helper;
//
//template<template<typename...> class T, typename... U>
//struct template_template_helper<T<U...>> {
//template<typename... V>
//using Rebind = T<V...>;
//
//using Original = std::tuple<U...>;
//};

template<typename T>
void foo(T &&val) {
    template<typename... F>
    using Z<F...> = T;
}


int main() {
//using H = template_template_helper<std::vector<int>>;
//static_assert(std::is_same<H::Rebind<char>, std::vector<char>>::value, "");
//static_assert(std::is_same<H::Original, std::tuple<int, std::allocator<int>>>::value, "");

