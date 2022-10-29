#include <iostream>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

int main() {
    MatrixXd m = MatrixXd::Random(3,3);
    m = (m + MatrixXd::Constant(3,3,1.0)) * 10;

    std::cout << "m = " << std::endl << m << std::endl;

    VectorXd v(3);

    v << 1,0,0;

    std::cout << "m * v = " << std::endl << m * v << std::endl;
}