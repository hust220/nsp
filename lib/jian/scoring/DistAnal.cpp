#include "../pdb.hpp"
#include "../geom.hpp"
#include "DistAnal.hpp"

namespace jian {

void DistAnal::free_freq(double *f) {
    delete f;
}

DistAnal::DistAnal(double interval, int cutoff) {
    this->interval = interval;
    this->cutoff = cutoff;
    bins = int(ceil(cutoff / interval));
}

DistAnal::~DistAnal() {
    delete [] num;
    delete [] type;
    delete [] ntLen;
    for (int i = 0; i < len; i++) {
        delete [] list[i];
    }
    delete [] list;

    free_freq(freq);
}

int DistAnal::res_type(std::string name) {
    return (name == "A" ? 1 : (name == "U" ? 2 : (name == "G" ? 3 : (name == "C" ? 4 : (name == "T" ? 2 : -1)))));
}

void DistAnal::read_mol(const Chain &chain) {
    delete [] num;
    delete [] type;
    delete [] ntLen;
    for (int i = 0; i < len; i++) {
        delete [] list[i];
    }
    delete [] list;

    len = chain.size();
    num = new int[len];
    type = new int[len];
    ntLen = new int[len];
    list = new Point *[len];

    Point p;
    int m = 0, n = 0;
    for (auto && res : chain) {
        m++;
        ntLen[n] = res.size();
        type[n] = res_type(res.name);
        list[n] = new Point[ntLen[n]];
        int k = 0;
        for (auto && atom : res) {
            for (int l = 0; l < 3; l++) list[n][k][l] = atom[l];
            list[n][k][3] = (res[0].name == "P" ? 0 : 3);
            list[n][k][3] += (type[n] == 1 ? k : (type[n] == 2 ? 22+k : (type[n] == 3 ? 42+k : 65+k)));
            if (atom.name == "O5*" && n != 0) {
                double dist = geom::distance(list[n][k], p);
                if (dist > 4) m += 10000; 
            } else if (atom.name == "O3*") {
                for (int l = 0; l < 3; l++) p[l] = list[n][k][l];
            }
            k++;
        }
        num[n] = m;
        n++;
    }
}

void DistAnal::train() {
/*
    for (int i = 0; i < len; i++) {
        for (int j = i + 1; j < len; j++) {
            for (int k = 0; k < ntLen[i]; k++) {
                for (int l = 0; l < ntLen[j]; l++) {
                    int type1 = list[i][k][3];
                    int type2 = list[j][l][3];
                    if (num[j] - num[i] == 1 && (((type1 >= 0 && type1 <= 11) || 
                            (type1 >= 22 && type1 <= 33) || (type1 >= 42 && type1 <= 53) || 
                            (type1 >= 62 && type1 <= 73)) && ((type2 >= 0 && type2 <= 11) || 
                            (type2 >= 22 && type2 <= 33) || (type2 >= 42 && type2 <= 53) || 
                            (type2 >= 62 && type2 <= 73)))) 
                        continue;
                    double temp = geom::distance(list[i][k], list[j][l]);
                    if (temp >= cutoff) continue;
                    obsParm[(list[i][k][3] * 85 + list[j][l][3]) * bins + int(temp / interval)]++;
                    obsParm[(list[j][l][3] * 85 + list[i][k][3]) * bins + int(temp / interval)]++;
                }
            }
        }
    }
*/
}

void DistAnal::read_parm(std::string filename) {
    std::ifstream ifile(filename.c_str());
    freq = new double [85 * 85 * 67];
    for (int i = 0; i < 85; i++) {
        for (int j = 0; j < 85; j++) {
            for (int k = 0; k < 67; k++) {
                ifile >> freq[(i * 85 + j) * 67 + k];
            }
        }
    }
    ifile.close();
}

DistAnal &DistAnal::run(const Chain &chain) {
    read_mol(chain);
    score = 0;
    for (int i = 0; i < len; i++) {
        for (int j = i + 1; j < len; j++) {
            for (int k = 0; k < ntLen[i]; k++) {
                for (int l = 0; l < ntLen[j]; l++) {
                    int type1 = list[i][k][3];
                    int type2 = list[j][l][3];
                    if (num[j] - num[i] == 1 && (((type1 >= 0 && type1 <= 11) || (type1 >= 22 && type1 <= 33) || (type1 >= 42 && type1 <= 53) || (type1 >= 62 && type1 <= 73)) && ((type2 >= 0 && type2 <= 11) || (type2 >= 22 && type2 <= 33) || (type2 >= 42 && type2 <= 53) || (type2 >= 62 && type2 <= 73)))) continue;
                    double temp = geom::distance(list[i][k], list[j][l]);
                    if (temp >= cutoff) continue;
                    double a = freq[int(list[i][k][3] * 85 + list[j][l][3]) * bins + int(temp / interval)];
                    double b = freq[int(list[j][l][3] * 85 + list[i][k][3]) * bins + int(temp / interval)];
                    if (a != 0) {
                        if (((type1 > 11 && type1 < 22) || (type1 > 33 && type1 < 42) || (type1 > 54 && type1 < 62) || (type1 > 74 && type1 < 85)) && ((type2 > 11 && type2 < 22) || (type2 > 33 && type2 < 42) || (type2 > 54 && type2 < 62) || (type2 > 74 && type2 < 85))) {
                            score -= 2.5 * log(a);
                        } else {
                            score -= log(a);
                        }
                    } else {
                        score += penalty;
                    }
                    if (b != 0) {
                        if (((type1 > 11 && type1 < 22) || (type1 > 33 && type1 < 42) || (type1 > 54 && type1 < 62) || (type1 > 74 && type1 < 85)) && ((type2 > 11 && type2 < 22) || (type2 > 33 && type2 < 42) || (type2 > 54 && type2 < 62) || (type2 > 74 && type2 < 85))) {
                            score -= 2.5 * log(b);
                        } else {
                            score -= log(b);
                        }
                    } else {
                        score += penalty;
                    }
                }
            }
        }
    }
    score = score / (len * (len - 1));
    return *this;
}

} // namespace jian

