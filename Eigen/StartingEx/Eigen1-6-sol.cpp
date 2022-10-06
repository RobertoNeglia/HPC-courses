#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

using namespace std;
using namespace Eigen;

int main(int argc, char** argv) {
    // Load matrix
    SparseMatrix<double> mat;
    loadMarket(mat, "/home/jellyfish/shared-folder/Matrices/gr_30_30.mtx");

    VectorXd xe = VectorXd::Constant(mat.rows(), 1); // define exact solution
    VectorXd b = mat*xe; // compute right-hand side

    // Solving
    // SimplicialLDLT<Eigen::SparseMatrix<double> > solver(mat); // performs LDLT factorization of mat
    SparseLU<SparseMatrix<double> > solver; // performs LU factorization of mat
    solver.compute(mat);

    if(solver.info()!=Success) { // first sanity check
        cout << "cannot factorize the matrix" << endl; // decomposition failed
        return 0;
    }
    
    VectorXd x = solver.solve(b); // use the factorization to solve for the
    if(solver.info()!=Success) { // second sanity check
        cout << "cannot solve the linear system" << endl; // solving failed
        return 0;
    }
    double relative_error = (x-xe).norm()/(xe).norm(); // compute and print relative error
    cout << relative_error << endl;
    return 0;
}