#include "fasta.hpp"

BEGIN_JN

void fasta_read(Fasta &fa, Str fastafile) {
    FastaItem x;
    for (auto &&it : FileLines(fastafile)) {
        if (it.line[0] == '>') {
            fa.push_back(x);
            x.name = (JN_ trim_copy(it.line.substr(1)));
            x.seq = "";
        }
        else {
            x.seq += JN_ trim_copy(it.line);
        }
    }
    fa.push_back(x);
    fa.pop_front();
}

static void fasta_write_seq(std::ostream &stream, Str seq) {
    for (Int i = 0; i < size(seq); i++) {
        stream << seq[i];
        if (i % 60 == 59 || i == size(seq) - 1) {
            stream << std::endl;
        }
    }
}

void fasta_write(std::ostream &stream, const Fasta &fa, Int n) {
    for (Int i = 0; i < size(fa); i++) {
        if (n == -1 || i < n) {
            stream << "> " << fa[i].name << std::endl;
            fasta_write_seq(stream, fa[i].seq);
        }
    }
}

END_JN

