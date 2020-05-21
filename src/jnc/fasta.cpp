#include "fasta.hpp"

namespace jian {

void fasta_read(Fasta &fa, Str fastafile) {
    FastaItem x;
    for (auto &&it : FileLines(fastafile)) {
        if (it.line[0] == '>') {
            fa.push_back(x);
            x.name = (jian::trim_copy(it.line.substr(1)));
            x.seq = "";
        }
        else {
            x.seq += jian::trim_copy(it.line);
        }
    }
    fa.push_back(x);
    fa.pop_front();
}

static void fasta_write_seq(std::ostream &stream, Str seq, int max = 60) {
    for (Int i = 0; i < size(seq); i++) {
        stream << seq[i];
        if ((max != -1 && i % max == max-1) || i == size(seq) - 1) {
            stream << std::endl;
        }
    }
}

void fasta_write(std::ostream &stream, const Fasta &fa, Int n, Int line_width) {
    for (Int i = 0; i < size(fa); i++) {
        if (n == -1 || i < n) {
            stream << "> " << fa[i].name << std::endl;
            fasta_write_seq(stream, fa[i].seq, line_width);
        }
    }
}

}

