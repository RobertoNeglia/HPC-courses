#include <iostream>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

int main(int argc, char** argv) {
    SparseMatrix<double> mat(10,10);                        // define matrix
    for(int i = 0; i < 10; i++) {
        mat.coeffRef(i,i) = 1.0;
    }

    cout << "Matrix: " << endl << mat << endl;

    VectorXd b = VectorXd::Constant(mat.rows(), 1);         // degine right-hand side

    cout << "Right-hand side: " << endl << b << endl;

    //Solving
    ConjugateGradient<SparseMatrix<double>> solver(mat);
    solver.compute(mat);
    if(solver.info() != Success) {
        cout << "cannot factorize the matrix" << endl;
        return 0;
    }

    VectorXd x = solver.solve(b);

    cout << "Solution: " << endl << x << endl;
    return 0;
}