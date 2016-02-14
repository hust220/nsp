#ifndef JIAN_DG_SMOOTH_H
#define JIAN_DG_SMOOTH_H

#include "Job.h"

namespace jian {
namespace dg {

class Smooth : public virtual Job {
public:
    Smooth() = default;
    Smooth(const Smooth &) = default;
    Smooth &operator =(const Smooth &) = default;

    void smooth() {
        /// check minimum distance
        for (int i = 0; i < bound.rows(); i++) for (int j = 0; j < bound.cols(); j++) {
            if (i == j) continue;
            if (bound(i, j) < _min_dist) bound(i, j) = _min_dist;
        }

        /// smooth
        int flag = 0, temp = 999;
        int n = 0;
        while (temp && n < 5) {
            flag = 0;
            for (int i = 0; i < len - 2; i++) for (int j = i + 1; j < len - 1; j++) for (int k = j + 1; k < len; k++) {
                /** max **/
                double max;
                /* i, j */
                max = at(bound, j, k, len) + at(bound, i, k, len);
                if (at(bound, i, j, len) > max && at(bound, j, i, len) < max) {
                    assign(bound, i, j, max, len);
                    flag++;
                }
                /* j, k */
                max = at(bound, i, j, len) + at(bound, i, k, len);
                if (at(bound, j, k, len) > max && at(bound, k, j, len) < max) {
                    assign(bound, j, k, max, len);
                    flag++;
                }
                /* i, k */
                max = at(bound, i, j, len) + at(bound, j, k, len);
                if (at(bound, i, k, len) > max && at(bound, k, i, len) < max) {
                    assign(bound, i, k, max, len);
                    flag++;
                }

                /** min **/
                /* i, j */
                double min = 0;
                double min1 = at(bound, k, j, len) - at(bound, i, k, len);
                double min2 = at(bound, k, i, len) - at(bound, j, k, len);
                if (min1 > 0) min = min1; else if (min2 > 0) min = min2;
                if (at(bound, j, i, len) < min && at(bound, i, j, len) > min) {
                    assign(bound, j, i, min, len);
                    flag++;
                }
                /* j, k */
                min = 0;
                min1 = at(bound, j, i, len) - at(bound, i, k, len);
                min2 = at(bound, k, i, len) - at(bound, i, j, len);
                if (min1 > 0) min = min1; else if (min2 > 0) min = min2;
                if (at(bound, k, j, len) < min && at(bound, j, k, len) > min) {
                    assign(bound, k, j, min, len);
                    flag++;
                }
                /* i, k */
                min = 0;
                min1 = at(bound, j, i, len) - at(bound, j, k, len);
                min2 = at(bound, k, j, len) - at(bound, i, j, len);
                if (min1 > 0) min = min1; else if (min2 > 0) min = min2;
                if (at(bound, k, i, len) < min && at(bound, i, k, len) > min) {
                    assign(bound, k, i, min, len);
                    flag++;
                }
            }
            if (flag == temp) n++;
            temp = flag;
        }

    }

    double at(const Mat &bound, int i, int j, int len) const {
        if (j >= len) {
            j -= len;
            return at(bound, j, i, len);
        }
        if (i >= len) {
            i -= len;
            return at(bound, j, i, len);
        }
        return bound(i, j);
    }

    void assign(Mat &bound, int i, int j, double d, int len) {
        if (j >= len) {
            j -= len;
            assign(bound, j, i, d, len);
            return;
        }
        if (i >= len) {
            i -= len;
            assign(bound, j, i, d, len);
            return;
        }
        bound(i, j) = d;
    }

};

} // namespace dg
} // namespace jian

#endif


