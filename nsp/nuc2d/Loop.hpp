#pragma once

#include "jian/utils/traits.hpp"
#include "jian/utils/ListRange.hpp"
#include "res.hpp"

BEGIN_JN

class Loop : public Deque<res>
{
    public:

        Str ss() const {
            STD_ ostringstream stream;
            for (auto &&res : *this) stream << res.type;
            return stream.str();
        }

        Str seq() const {
            STD_ ostringstream stream;
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

        friend STD_ ostream &operator <<(STD_ ostream &stream, Loop l) {
            stream << "Loop: " << ' ' << (Str)l;
            return stream;
        }

};

END_JN
