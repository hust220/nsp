#pragma once

#include <functional>
#include <vector>

BEGIN_JN

template <typename T>
class FSM {
public:
    using table_t = std::vector<std::vector<int>>;
    using input_t = T;
    using map_t = std::function<int(const input_t &)>;
    using state_t = int;

    table_t table;
    map_t map;
    state_t state;

    FSM(const table_t &t, const map_t &m) : table(t), map(m) {
        state = 0;
    }

    void init() {
        state = 0;
    }

    state_t accept(const input_t &input) {
        state = table[state][map(input)];
        return state;
    }
};

}

