#include <jian/utils/file.hpp>

BEGIN_JN

struct FastaItem {
    Str name;
    Str seq;
};

using Fasta = Deque<FastaItem>;

void fasta_read(Fasta &fa, Str fastafile);

void fasta_write(std::ostream &stream, const Fasta &fa, Int n = -1);

END_JN

