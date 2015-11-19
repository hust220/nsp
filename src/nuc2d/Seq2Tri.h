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
    typedef std::tuple<int, int, int> Tuple;
    typedef std::vector<Tuple> Tuples;
    typedef std::vector<std::vector<int>> TupleList;
    typedef std::pair<TupleList, double> TupleInfo;
    typedef std::vector<TupleInfo> InfoList;

    Seq2Tri();
    InfoList operator ()(std::string seq);
    InfoList sub_info(const std::vector<int> &sequence, double cutoff);
    InfoList best_info(const std::vector<int> &sequence);
    TupleInfo strand_info(const std::vector<int>  &sequence);
    TupleInfo score(const TupleInfo &info1, const TupleInfo &info2, int a, int b, int c);
    TupleInfo score(const TupleInfo &info, int m);
    std::pair<double, bool> tuple_energy(char i, char j, char k);
    std::pair<double, bool> pair_energy(char m, char n);
    double stack_energy(char i, char j);

    std::unordered_map<std::vector<int>, InfoList, seq2tri::MyHash> _info;

    std::string _seq;
    std::vector<int> _types;
    std::vector<std::tuple<int, int, int>> _pairs;
    std::map<char, int> _convert{{'A', 0}, {'U', 1}, {'T', 1}, {'G', 2}, {'C', 3}};
    std::hash<std::string> hash_fn;
    int _min_hairpin_size;
    double _cutoff;
    Matrix4f _pair_energy;
    Matrix4f _stack_energy;
};

} /// namespace nuc2d
} /// namespace jian

#endif
