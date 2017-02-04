#pragma once

BEGIN_JN

template<typename Data>
class Tree {
public:
    Data _data;
    std::list<Tree<Data>> _sons;

    Tree() {}
    Tree(const Tree &tree) : _data(tree._data), _sons(tree._sons) {}
    Tree(Tree &&tree) : _data(std::move(tree._data)), _sons(std::move(tree._sons)) {}
    Tree<Data> &operator =(const Tree &tree) { _data = tree._data; _sons = tree._sons; }
    Tree<Data> &operator =(Tree &&tree) { std::swap(_data, tree._data); std::swap(_sons, tree._sons); }
    Tree(const Data &data) : _data(data) {}
    Tree(Data &&data) : _data(std::move(data)) {}
    Tree<Data> &operator =(const Data &data) { _data = data; }
    Tree<Data> &operator =(Data &&data) { std::swap(_data, data); }

    Tree<Data> * self() {
        return this;
    }

    const Tree<Data> * self() const {
        return this;
    }

    Data &head() {
        return _data;
    }

    const Data &head() const {
        return _data;
    }

    std::list<Tree<Data>> &sons() {
        return _sons;
    }

    const std::list<Tree<Data>> &sons() const {
        return _sons;
    }

    void append() {}

    template<typename T1, typename... T2> void append(T1 &&t1, T2 && ...t2) {
        _sons.push_back(Tree<Data>(std::forward<T1>(t1)));
        append(std::forward<T2>(t2)...);
    }

    template<typename Fun> void apply(Fun &&f) {
        f(*this);
        for (auto &&tree: _sons) tree.apply(std::forward<Fun>(f));
    }

    template<typename Fun> void apply(Fun &&f) const {
        f(*this);
        for (auto &&tree: _sons) tree.apply(std::forward<Fun>(f));
    }

    template<typename Fun, typename Con> auto fold(Fun &&f, Con &&con) const -> decltype(f(_data)) {
        auto result = f(_data);
        for (auto &&tree: _sons) result = con(result, tree.fold(f, std::forward<Con>(con)));
        return result;
    }

};

END_JN

