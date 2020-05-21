#pragma once

#include "jian.hpp"
#include "entity_list_range.hpp"
#include "rss_res.hpp"

namespace jian {

class Loop : public Deque<res>
{
    public:

        Str ss() const {
            std::ostringstream stream;
            for (auto &&res : *this) stream << res.type;
            return stream.str();
        }

        Str seq() const {
            std::ostringstream stream;
            for (auto &&res : *this) {
                if (res.type != '&') {
                    stream << res.name;
                }
            }
            return stream.str();
        }

        Deque<Int> nums() const {
            Deque<Int> ls;
            for (auto && res : *this) {
                ls.push_back(res.num - 1);
            }
            return ls;
        }

        operator Str() const {
            std::ostringstream stream;
            stream << seq() << ' ' << ss() << ' ';
            for (auto && res : *this) stream << res.num << ' ';
            return stream.str();
        }

        friend std::ostream &operator <<(std::ostream &stream, Loop l) {
            stream << "Loop: " << ' ' << (Str)l;
            return stream;
        }

};

}
