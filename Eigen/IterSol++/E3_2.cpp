    #include <cstdlib>                      // System includes
    #include <iostream>                      
    #include <Eigen/SparseCore>
    #include <Eigen/IterativeLinearSolvers>
    #include <unsupported/Eigen/SparseExtra>

    using std::endl;
    using std::cout;

    #include "bcgstab.hpp"                          

    int main(int argc, char** argv)
    {
    using namespace LinearAlgebra;
    // Some useful alias
    using SpMat=Eigen::SparseMatrix<double>;
    using SpVec=Eigen::VectorXd;
    

    SpMat M;
    Eigen::loadMarket(M, "bcsstm12.mtx");

    M = M + SpMat(M.transpose());

    int n = M.rows();

        


    double tol = 1.e-10;                 // Convergence tolerance
    int result, maxit = 1000;           // Maximum iterations

    std::cout<<"Matrix size:"<<M.rows()<<"X"<<M.cols()<<std::endl;
    std::cout<<"Non zero entries:"<<M.nonZeros()<<std::endl;

    SpMat B = SpMat(M.transpose()) - M;  // Check symmetry
    std::cout<<"Norm of M-M.t: "<<B.norm()<<std::endl;

    // Create Rhs b
    SpVec e=SpVec::Ones(M.rows());
    SpVec b=M* e;
    SpVec x(M.rows());
    Eigen::DiagonalPreconditioner<double> D(M); // Create diagonal preconditioner

    // First with eigen bicgstab
    Eigen::BiCGSTAB<SpMat> bicgstab;
    bicgstab.setMaxIterations(maxit);
    bicgstab.setTolerance(tol);
    bicgstab.compute(M);
    x = bicgstab.solve(b);
    std::cout <<" Eigen native bicgstab"<<std::endl;
    std::cout << "#iterations:     " << bicgstab.iterations() << std::endl;
    std::cout << "estimated error: " << bicgstab.error()      << std::endl;
    std::cout << "effective error: "<<(x-e).norm()<<std::endl;

    // Now with hand-made bicgstab
    x=0*x;
    result = BiCGSTAB(M, x, b, D, maxit, tol);        // Solve system

    std::cout <<" hand-made bicgstab "<<std::endl;
    cout << "bicgstab flag = " << result << endl;
    cout << "iterations performed: " << maxit << endl;
    cout << "tolerance achieved  : " << tol << endl;
    std::cout << "Error norm: "<<(x-e).norm()<<std::endl;

    return result;
    }
