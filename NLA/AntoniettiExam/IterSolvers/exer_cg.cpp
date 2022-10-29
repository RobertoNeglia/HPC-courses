#include <cstdlib>                      // System includes
#include <iostream>                      
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/SparseExtra>

using std::endl;
using std::cout;

#include "cg.hpp"                          

int main(int argc, char** argv)
{
  using namespace LinearAlgebra;
  // Some useful alias
  using SpMat=Eigen::SparseMatrix<double>;
  using SpVec=Eigen::VectorXd;

  // Load matrix
  SpMat A; 
  Eigen::loadMarket(A, "bcsstm12.mtx");
  A = SpMat(A.transpose()) + A;

  double tol = 1.e-14;                 // Convergence tolerance
  int result, maxit = 5000;           // Maximum iterations

  std::cout<<"Matrix size:"<<A.rows()<<"X"<<A.cols()<<std::endl;
  std::cout<<"Non zero entries:"<<A.nonZeros()<<std::endl;

  // Create Rhs b
  SpVec e = SpVec::Ones(A.rows());
  SpVec b = A*e;
  SpVec x(A.rows());

  // Create preconditioner
  Eigen::DiagonalPreconditioner<double> D(A); 
  Eigen::IncompleteLUT<double> ILU(A);

  // Eigen CG with diagonal precond
  Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper> cg;
  cg.setMaxIterations(maxit);
  cg.setTolerance(tol);
  cg.compute(A);
  x = cg.solve(b);
  std::cout <<" Eigen native CG"<<std::endl;
  std::cout << "#iterations:     " << cg.iterations() << std::endl;
  std::cout << "estimated error: " << cg.error()      << std::endl;
  std::cout << "effective error: "<<(x-e).norm()<<std::endl;

  // Hand-made CG with diagonal precond
  x = 0*x;
  result = CG(A, x, b, D, maxit, tol);        // Solve system

  std::cout <<" hand-made CG "<<std::endl;
  cout << "iterations performed: " << maxit << endl;
  cout << "tolerance achieved  : " << tol << endl;
  std::cout << "Error norm: "<<(x-e).norm()<<std::endl;

  // Eigen CG with ILU precond
  Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper, Eigen::IncompleteLUT<double>> cg_ilu;
  cg_ilu.setMaxIterations(maxit);
  cg_ilu.setTolerance(tol);
  cg_ilu.compute(A);
  x = cg_ilu.solve(b);
  std::cout <<" Eigen native CG"<< std::endl;
  std::cout << "#iterations:     " << cg_ilu.iterations() << std::endl;
  std::cout << "estimated error: " << cg_ilu.error()      << std::endl;
  std::cout << "effective error: "<<(x-e).norm()<<std::endl;

  // Hand-made CG with ILU precond
  x = 0*x;
  result = CG(A, x, b, ILU, maxit, tol);        
  std::cout <<" hand-made CG "<< std::endl;
  cout << "iterations performed: " << maxit << endl;
  cout << "tolerance achieved  : " << tol << endl;
  std::cout << "Error norm: "<<(x-e).norm()<<std::endl;

  return result;
}
