#ifndef JIAN_NUC3D_PREDICT3D
#define JIAN_NUC3D_PREDICT3D

#include "../nuc2d/N2D.h"
#include "../nuc2d/util.h"
#include "../pdb.h"
#include "Assemble.h"
#include "triass/TriAss.h"

namespace jian {
namespace nuc3d {

class Predict3D {
public:
    using Constraints = std::deque<std::pair<int, int>>;

    std::map<std::string, std::function<BasicPredict3D*(const Par &)>> _methods {
        {"assemble", [](const Par &par){return new Assemble(par);}},
        {"tri-assemble", [](const Par &par){return new triass::TriAss(par);}}
    };

    Log log;

    void operator ()(const Par &par) {
        auto job = JobPredict3D(par);
        log = job.log;
        auto c = read_constraints(job);
        display_start_information_job(job);
        auto method = _methods[job._method](par);
        int max_it_num = 100;
        for (int i = 0; i < job._num; i++) {
            for (int j = 0; j < max_it_num; j++) {
                auto model = method->predict();
                if (satisfy(model, c) || j == max_it_num - 1) {
                    write_pdb(model, job._name + "-" + boost::lexical_cast<std::string>(i + 1) + ".pdb");
                    break;
                }
            }
        }
        delete method;
        display_end_information_job(job);
    }

    Constraints read_constraints(const JobPredict3D &job) {
        Constraints c;
        if (job._constraints == "") return c;
        std::ifstream ifile(job._constraints.c_str());
        std::string line; std::vector<std::string> v;
        while (ifile >> line) {
            tokenize(line, v, " ");
            if (v.size() == 2) c.push_back({boost::lexical_cast<int>(v[0]), boost::lexical_cast<int>(v[1])});
        }
        ifile.close();
        return c;
    }

    void display_start_information_job(JobPredict3D &job) {
        log.clear();
        job._start_time = std::time(nullptr);
        log("=========================================================\n",
            "New Job: ", job._name, '\n',
            "Time: ", std::asctime(std::localtime(&(job._start_time))), '\n',
            "Sequence: ", job._seq, '\n',
            "2D Structure: ", job._ss, '\n',
            "Molecular Type: ", job._type, '\n',
            "Number: ", job._num, '\n',
            "Number of Sampling: ", job._num_sampling, '\n',
            "Method: ", job._method, '\n',
            "----------------------------------------\n");
    }

    void display_end_information_job(JobPredict3D &job) {
        job._end_time = std::time(nullptr);
        log("----------------------------------------\n",
            "Finish Time: ", std::asctime(std::localtime(&(job._end_time))), '\n',
            "Time elapsed: ", job._end_time - job._start_time, "s\n",
            "=========================================================\n\n");
    }

    template<typename T>
    bool satisfy(T &&model, const Constraints &c) {
        double cutoff = 20;
        int i = 0; for (auto &&chain1 : model) for (auto &&res1 : chain1) {
            int j = 0; for (auto &&chain2 : model) for (auto &&res2 : chain2) {
                if (belong(i, j, c) && residue_distance(res1, res2) > cutoff) return false;
                j++;
            }
            i++;
        }
        return true;
    }

    bool belong(int i, int j, const Constraints &c) {
        for (auto &&pair : c) if (pair.first == i || pair.second == j) return true;
        return false;
    }

    template<typename T, typename U>
    double residue_distance(T &&r1, U &&r2) {
        return geom::distance(r1["C4*"], r2["C4*"]);
    }

};

} // namespace nuc3d
} // namespace jian

#endif

