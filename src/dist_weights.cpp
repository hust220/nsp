#include "nsp.hpp"
#include <jian/matrix.hpp>

namespace jian {

static double alpha(double max, double min, double beta) {
    if (max > 0 && min >= 0) {
        return beta / ((1 - beta) * max);
    } else if (min < 0 && max <= 0) {
        return -beta / ((1 + beta) * min);
    } else if (max > 0 && min < 0) {
        return std::min(beta / ((1 - beta) * max), -beta / ((1 + beta) * min));
    } else if (min == 0 && max == 0) {
        return 1;
    }
}

REGISTER_NSP_COMPONENT(dist_weights) {
    double beta = std::stod(par["beta"][0]);
    std::string file_name = par["f"][0];
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(85, 85);
    std::ifstream ifile(file_name.c_str());
    int sum = 0;
    int d;
    for (int i = 0; i < 85; i++) {
        for (int j = 0; j < 85; j++) {
            for (int k = 0; k < 67; k++) {
                ifile >> d;
                mat(i, j) += d;
                sum += d;
            }
        }
    }
    ifile.close();
    mat /= double(sum);
//    std::cout << mat << std::endl;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(mat);
    Eigen::VectorXd val_us = eigen_solver.eigenvalues();
    Eigen::MatrixXd vec_us = eigen_solver.eigenvectors();
    double a = alpha(val_us.maxCoeff(), val_us.minCoeff(), beta);
    Eigen::MatrixXd m = mat * a;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver2(m);
    Eigen::VectorXd val = eigen_solver2.eigenvalues();
    Eigen::MatrixXd vec = eigen_solver2.eigenvectors();
//    std::cout << m << std::endl;
    Eigen::MatrixXd lambdas = Eigen::MatrixXd::Zero(85, 85);
    for (int i = 0; i < 85; i++) {
        lambdas(i, i) = val[i] / (val[i] + 1);
    }
    Eigen::MatrixXd direct = vec * lambdas * vec.inverse();
    std::cout << direct << std::endl;
//    std::cout << direct + direct * direct + direct * direct * direct + direct * direct * direct * direct - m << std::endl;
//    for (int i = 0; i < 67; i++) {
//        (*mats[i]) = (*mats[i]) * (Eigen::MatrixXd::Identity(85, 85) + (*mats[i])).inverse();
//    }
//    std::vector<Eigen::MatrixXd *> mats_ref(67);
//    for (int i = 0; i < 67; i++) {
//        mats_ref[i] = new Eigen::MatrixXd(85, 85);
//    }
//    for (int i = 0; i < 85; i++) {
//        for (int j = 0; j < 85; j++) {
//            for (int k = 0; k < 67; k++) {
//                (*(mats_ref[k]))(i, j) = arr[k] / sum;
//            }
//        }
//    }
//    for (int i = 0; i < 67; i++) {
//        (*mats_ref[i]) = (*mats_ref[i]) * (Eigen::MatrixXd::Identity(85, 85) + (*mats_ref[i])).inverse();
//    }
//    double d;
//    for (int i = 0; i < 85; i++) {
//        for (int j = 0; j < 85; j++) {
//            for (int k = 0; k < 67; k++) {
//                d = (*(mats_ref[k]))(i, j);
//                if (d == 0) {
//                    (*(mats[k]))(i, j) = 0;
//                } else {
//                    (*(mats[k]))(i, j) /= d;
//                }
//                std::cout << (*(mats[k]))(i, j) << ' ';
//            }
//            std::cout << std::endl;
//        }
//    }
//    for (int i = 0; i < 67; i++) {
//        delete mats[i];
//        delete mats_ref[i];
//    }
}

} // namespace jian

