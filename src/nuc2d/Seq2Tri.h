#ifndef JIAN_NUC2D_SEQ2TRI_H
#define JIAN_NUC2D_SEQ2TRI_H

#include <util/util.h>

namespace jian {
namespace nuc2d {

namespace seq2tri {
struct MyHash {
    std::size_t operator ()(const std::vector<int> &s) const {
        long hash = 5381;  
        for(int i = 0; i < s.size(); i++)  {  
            hash = ((hash << 5) + hash) + s[i];  
        }  
        return hash;
    }
};
}

class Seq2Tri {
public:
    std::vector<std::tuple<int, int, int>> operator ()(std::string seq);
    double get(const std::vector<int> &);
    double score(char, char, char);
    double score(int, int, int);
    std::string num_to_seq(const std::vector<int> &);
    std::vector<std::vector<std::tuple<int, int, int>>> backtrack(const std::vector<int> &);

    std::unordered_map<std::vector<int>, double, seq2tri::MyHash> _scores;

    std::string _seq;
    std::vector<int> _types;
    std::vector<std::tuple<int, int, int>> _pairs;
    std::map<char, int> _convert{{'A', 0}, {'U', 1}, {'G', 2}, {'C', 3}};
    std::hash<std::string> hash_fn;
    int _min_hairpin_size = 4;
};

} /// namespace nuc2d
} /// namespace jian

#endif
