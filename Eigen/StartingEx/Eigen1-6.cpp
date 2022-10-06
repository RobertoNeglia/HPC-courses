#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

using namespace std;
using namespace Eigen;

int main(int argc, char** argv) {
    // Load matrix
    SparseMatrix<double> mat;
    loadMarket(mat, "/home/jellyfish/shared-folder/Matrices/gr_30_30.mtx");

    VectorXd xe = VectorXd::Constant(mat.rows(), 1);        // exact solution
    VectorXd b = mat*xe;                                    // compute right-hand side

    SparseLU<SparseMatrix<double>> solver(mat);
    solver.compute(mat);
    if(solver.info() != Success) {
        cout << "cannot factorize matrix." << endl;
        return 0;
    }

    VectorXd x = solver.solve(b);

    double relative_error = (x-xe).norm() / xe.norm();

    cout << relative_error << endl;

    return 0;
}