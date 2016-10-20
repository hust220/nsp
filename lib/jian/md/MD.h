#ifndef JIAN_MD_MD
#define JIAN_MD_MD

#include "../etl.h"

namespace jian {
namespace md {

template<typename Mat = MatrixXf>
class MD {
public:
    typedef decltype(ref(Mat(), 0, 0)) Val;

    int _rand_seed = 12345;
    std::string _name = "test"
    int _max_steps = 1000;
    int _record_cycles = 10;
    double _atom_radius = 0.5;
    int _num_atoms;
    double _min_time = 0.001;
    double _box_size = 20;

    void operator ()(int num_atoms, const Mat &pos, const Mat &vel) {
        _num_atoms = num_atoms;
        md(pos, vel);
    }

    void md(const Mat &init_pos,const Mat &init_vel) {
        std::list<Mat> pos_list, vel_list;
        pos_list.push_back(init_pos);
        vel_list.push_back(init_vel);
        for (int i = 0; i < _max_steps; i++) {
            auto &pos = pos_list.back();
            auto &vec = vec_list.back();
            auto acc = get_acc(pos);
            pos_list.push_back(get_pos(pos, vec, acc));
            vec_list.push_back(get_vec(pos, vec, acc));
            pos_list.pop_front();
            vec_list.pop_front();
            log("step ", i + 1, '\n', pos_list.back(), '\n');
        }
    }

    Mat get_acc(const Mat &pos) {
        Mat acc = make_mat<Mat>(_num_atoms, 3);
        for (int i = 0; i < _num_atoms; i++) for (int j = 0; j < 3; j++) {
            ref(acc, i, j) = 0;
        }
        return acc;
    }

    Mat get_pos(const Mat &pos, const Mat &vel, const Mat &acc) {
        Mat pos = make_mat<Mat>(_num_atoms, 3);
        for (int i = 0; i < _num_atoms; i++) for (int j = 0; j < 3; j++) {
            set_pos_coeff(pos, i, j, ref(pos, i, j) + ref(vel, i, j) * _min_time);
        }
        return pos;
    }

    void set_pos_coeff(const Mat &pos, int i, int j, const Val &val) {
        if (val > _box_size / 2) ref(pos, i, j) = val - _box_size;
        else if (val < -_box_size / 2) ref(pos, i, j) = val + _box_size;
        else ref(pos, i, j) = val;
    }

    Mat get_vel(const Mat &pos, const Mat &vel, const Mat &acc) {
        Mat vel = make_mat<Mat>(_num_atoms, 3);
        for (int i = 0; i < _num_atoms; i++) for (int j = 0; j < 3; j++) {
            ref(vel, i, j) += ref(acc, i, j) * _min_time;
        }
        return pos;
    }
};

} // namespace md
} // namespace jian

#endif

