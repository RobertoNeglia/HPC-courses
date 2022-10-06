#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

using namespace std;
using namespace Eigen;

int main(int argc, char** argv) {
    // Load matrix
    SparseMatrix<double> mat;
    loadMarket(mat, "/home/jellyfish/shared-folder/Matrices/mhd416a.mtx");

    VectorXd xe = VectorXd::Constant(mat.rows(), 1);
    VectorXd b = mat * xe;

    cout << b << endl;

    return 0;

}