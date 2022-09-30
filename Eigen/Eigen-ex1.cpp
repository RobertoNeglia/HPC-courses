#include <Eigen/Eigen>
#include <iostream>

int main(int argc, char** argv) {
    int n;
    
    n = 12;
    
    Eigen::SparseMatrix<double> A(n,n);
    Eigen::VectorXd x;

    

    for(std::size_t i = 0; i < n; i++) {
        A.coeffRef(i,i) = 1.0;
    }

    x = VectorXd::Constant(n, 1);

    



}