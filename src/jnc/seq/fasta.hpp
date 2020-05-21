#include "file.hpp"

namespace jian {

struct FastaItem {
    Str name;
    Str seq;
};

using Fasta = Deque<FastaItem>;

void fasta_read(Fasta &fa, Str fastafile);

void fasta_write(std::ostream &stream, const Fasta &fa, Int n = -1, Int line_width = 60);

}

