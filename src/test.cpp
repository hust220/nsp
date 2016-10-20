#include "nsp.hpp"
#include <jian/mc/MC.hpp>
#include <jian/utils/string.hpp>
#include <jian/utils/math.hpp>

namespace jian {

namespace test_detail {

std::string cmd(const std::string &s) {
    FILE *fp = popen(s.c_str(), "r");
    if (fp != NULL) {
        char buf[1024];
        std::ostringstream stream;
        while (std::fgets(buf, 1024, fp) != NULL) {
            stream << buf;
        }
        pclose(fp);
        return stream.str();
    } else {
        return "";
    }
}

class Refine : public MC {
public: 
    using coeffs_t = std::array<float, 4>;

    int m_index = 0;
    float m_old_value;
    coeffs_t m_coeffs;
    std::array<std::array<std::array<float, 4>, 10000>, 32> x;
    std::array<std::array<float, 10000>, 32> y_;
    std::array<std::array<float, 10000>, 32> y;
    std::array<float, 32> mins;
    
    Refine(const std::string &f, const coeffs_t &coeffs) : m_coeffs(coeffs) {
        std::ifstream stream;
        stream.open(f.c_str());
        for (int i = 0; i < 32; i++) {
            for (int j = 0; j < 10000; j++) {
                for (int k = 0; k < 4; k++) {
                    stream >> x[i][j][k];
                }
                stream >> y[i][j];
            }
        }
        for (int i = 0; i < 32; i++) {
            mins[i] = *std::min_element(y[i].begin(), y[i].end());
        }
        stream.close();
    }

    virtual void mc_select() {
        m_index = 1+int(jian::rand() * 3);
    }

    virtual void mc_sample() {
        m_old_value = m_coeffs[m_index];

//        static std::vector<double> v {-2, -1.25, -0.8, -0.5, 0.5, 0.8, 1.25, 2};
        static std::vector<double> v {-1.25, -0.8, 0.8, 1.25};
        int d = int(jian::rand() * v.size());
        m_coeffs[m_index] *= v[d];
        if (m_coeffs[m_index] > 10000 || m_coeffs[m_index] < -10000) {
            mc_back();
            mc_sample();
        }
//        double d = (jian::rand() - 0.5) * 200;
//        m_coeffs[m_index] += d;
    }

    virtual void mc_back() {
        m_coeffs[m_index] = m_old_value;
    }

    virtual double mc_partial_energy() {
        double e = 0;
        int n;
        for (int i = 0; i < 32; i++) {
            for (int j = 0; j < 1000; j++) {
                y_[i][j] = 0;
                for (int k = 0; k < 4; k++) {
                    y_[i][j] += x[i][j][k] * m_coeffs[k];
                }
            }
            auto it = std::min_element(y_[i].begin(), y_[i].end());
            n = std::distance(y_[i].begin(), it);
            //e += square(y[i][n]-mins[i]);
            e += y[i][n];
        }
        return e/32.0;
    }

    virtual void mc_write() {
        std::cout << _mc_step << ' ' << m_coeffs[0] << ' ' << m_coeffs[1] << ' ' << m_coeffs[2] << ' ' << m_coeffs[3] << ' ' << mc_partial_energy() << "(en) " << _mc_tempr << "(tempr) " << _mc_local_succ_rate << "(rate)" << std::endl;
    }
};

} // namespace test_detail

REGISTER_NSP_COMPONENT(test) {
//    std::cout << test_detail::cmd(par[2]);
    test_detail::Refine::coeffs_t coeffs = {1, 1, 1, 1};
    if (par.has("init")) {
        int i = 0;
        for (auto && s : par["init"]) {
            coeffs[i] = std::stod(s);
            i++;
        }
    }
    if (par.has("seed")) {
        jian::seed(std::stod(par["seed"][0]));
    }
    test_detail::Refine refine(par["f"][0], coeffs);
    par.set(refine._mc_heat_steps, "heat_steps");
    par.set(refine._mc_init_tempr, "init_tempr");
    refine.mc();
}

} // namespace jian


