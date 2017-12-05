#include "nsp.hpp"
#include "matrix.hpp"

BEGIN_JN

REGISTER_NSP_COMPONENT(dist_frequencies) {
    S file_name = par["f"][0];
    std::vector<Eigen::MatrixXd *> mats(67);
    for (int i = 0; i < 67; i++) {
        mats[i] = new Eigen::MatrixXd(85, 85);
    }
    std::ifstream ifile(file_name.c_str());
    double arr[67];
    double sum;
    for (int i = 0; i < 85; i++) {
        for (int j = 0; j < 85; j++) {
            sum = 0;
            for (int k = 0; k < 67; k++) {
                ifile >> arr[k];
                sum += arr[k];
            }
            for (int k = 0; k < 67; k++) {
                (*(mats[k]))(i, j) = arr[k] / sum;
            }
        }
    }
    sum = 0;
    for (int k = 0; k < 67; k++) {
        ifile >> arr[k];
        sum += arr[k];
    }
    ifile.close();
    for (int k = 0; k < 67; k++) {
        arr[k] /= sum;
    }
    for (int i = 0; i < 85; i++) {
        for (int j = 0; j < 85; j++) {
            for (int k = 0; k < 67; k++) {
                if (arr[k] != 0) {
                    (*(mats[k]))(i, j) /= arr[k];
                }
                std::cout << (*(mats[k]))(i, j) << ' ';
            }
            std::cout << std::endl;
        }
    }
    for (int i = 0; i < 67; i++) {
        delete mats[i];
    }
}

END_JN

