# Domain Decomposition and Algebraic Multigrid preconditioners

In this lab we want to explore the Additive Schwarz and Algebraic Multigrid methods implemented in the LIS library in order to define advanced preconditioners for sparse linear systems.

## Additive Schwarz preconditioner

Domain decomposition (DD) methods are techniques based on a decomposition of the spatial domain of the problem into several subdomains. Such reformulations are usually motivated by the need to create solvers which are easily parallelized on coarse grain parallel computers, though sometimes they can also reduce the complexity of solvers on sequential computers.

In order to obtain a domain decomposition for the linear system $Ax = b$, one needs to decompose the unknowns in the vector $x$ into subsets. 
For the sake of simplicity, we consider only a two subdomain case. Let $\mathcal{N}_i$ , $i = 1,2$ be a partition of the indices corresponding to the vector $x$ and  $R_i$ , $i = 1,2$ denote the matrix that, when applied to the vector $x$, returns only those values associated with the indices in $N_i$. The matrices $R_i$ are often referred to as the restriction operators, while $R^{\rm T}_i$ are the interpolation matrixes.

One can define the restriction of the matrix $A$ to the first and second unknowns using the restriction matrices:
$$
A_j = R_j A R^{\rm T}_j,\quad j = 1, 2.
$$
Thus the matrix $A_j$ is simply the subblock of A associated with the indices in $\mathcal{N}_j$. The preconditioned system obtained by applying the Additive Schwarz method reads
$$
(R^{\rm T}_1 A^{-1}_1 R_1 + R^{\rm T}_2 A^{-1}_2 R_2) A x = (R^{\rm T}_1 A^{-1}_1 R_1 + R^{\rm T}_2 A^{-1}_2 R_2) b.
$$
Since the sizes of the matrices $A_1$ and $A_2$ are smaller than the original matrix $A$, the inverse matrices $A_1^{-1}$ and $A_2^{-1}$ are easier to compute. Using this preconditioner for a stationary iterative method yields
$$
x^{n+1} = x^n + (R^{\rm T}_1 A^{-1}_1 R_1 + R^{\rm T}_2 A^{-1}_2 R_2) (b- Ax^n).
$$

Even if the matrices $A_i$, $i=1,2$, are easier to be inverted, in LIS it is preferred not to compute the inverse matrices exactly. Another preconditioned technique must be used together with the Additive Schwarz method in order to define the matrices $B_i$ approximating $A^{-1}_i$. The corresponding preconditioner is computed as $R^{\rm T}_1 B_1 R_1 + R^{\rm T}_2 B_2 R_2$. We can check this by the following test:

```
mpirun -n 4 ./test1 testmat0.mtx 2 sol.mtx hist.txt -i cg

mpirun -n 4 ./test1 testmat0.mtx 2 sol.mtx hist.txt -i cg -adds true
```
We observe that the option `-adds true` needed to specify the use of the Additive Schwarz method does not provide any preconditioner. In both cases, we get the following result

```
number of processes = 4
matrix size = 100 x 100 (460 nonzero entries)

initial vector x      : all components set to 0
precision             : double
linear solver         : CG
preconditioner        : none
convergence condition : ||b-Ax||_2 <= 1.0e-12 * ||b-Ax_0||_2
matrix storage format : CSR
linear solver status  : normal end
CG: number of iterations = 15
```

If instead we write `mpirun -n 4 ./test1 testmat0.mtx 2 sol.mtx hist.txt -i cg -adds true -p jacobi`, we get
```
initial vector x      : all components set to 0
precision             : double
linear solver         : CG
preconditioner        : Jacobi + Additive Schwarz
convergence condition : ||b-Ax||_2 <= 1.0e-12 * ||b-Ax_0||_2
matrix storage format : CSR
linear solver status  : normal end

CG: number of iterations = 14
```

## Some examples:

Assess the performances of the Additive Schwarz preconditioner implemented in LIS on different sparse linear systems:

```
mpirun -n 2 ./test1 testmat2.mtx 2 sol.mtx hist.txt -i gmres
mpirun -n 2 ./test1 testmat2.mtx 2 sol.mtx hist.txt -i gmres -adds true -p ssor
```
In this case the number of iterations for reaching the desired tolerance reduces from 36 to 14! 

```
mpirun -n 2 ./test1 testmat2.mtx 2 sol.mtx hist.txt -i bicgstab
mpirun -n 2 ./test1 testmat2.mtx 2 sol.mtx hist.txt -i bicgstab -adds true -p ilu -ilu_fill 2
```
In this case the number of iterations reduces from 36 to 6 and also the total elapsed time for the computation is reduced. 

We can also set the number of iterations of the Additive Schwarz method using the `adds_iter` option. The effect of this parameter can be observed in the following example.
```
mpirun -n 4 ./test2 100 100 0 sol.mtx hist.txt -i gmres -p ssor -adds true -adds_iter 2
mpirun -n 4 ./test2 100 100 0 sol.mtx hist.txt -i gmres -p ssor -adds true -adds_iter 3
mpirun -n 4 ./test2 100 100 0 sol.mtx hist.txt -i gmres -p ssor -adds true -adds_iter 4
```
The number of iterations drops from 104 to 81 (using `adds_iter = 3`) and 70 (using `adds_iter = 4`). However we can observe that the total time for the solution of the linear system is similar because the time required for computing the preconditioner increases.


## Exercise 1

Evaluate and compare the effect of the Additive Schwarz preconditioner for solving the linear systems deefined by the matrices [bcsstm12.mtx](https://suitesparse-collection-website.herokuapp.com/MM/HB/bcsstm12.tar.gz) and [nos1.mtx](https://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/lanpro/nos1.mtx.gz).

```
wget https://suitesparse-collection-website.herokuapp.com/MM/HB/bcsstm12.tar.gz
tar -xf bcsstm12.tar.gz
mv bcsstm12/bcsstm12.mtx .
rm -rf bcsstm12 bcsstm12.tar.gz 

wget https://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/lanpro/nos1.mtx.gz
gzip -dk nos1.mtx.gz 

mpirun -n 4 ./test1 bcsstm12.mtx 2 sol.txt hist.txt -i bicgstab -maxiter 5000 -tol 1e-10

mpirun -n 4 ./test1 bcsstm12.mtx 2 sol.txt hist.txt -i bicgstab -tol 1e-10 -p sainv 

mpirun -n 4 ./test1 bcsstm12.mtx 2 sol.txt hist.txt -i bicgstab -tol 1e-10 -p sainv -adds true

mpirun -n 4 ./test1 bcsstm12.mtx 2 sol.txt hist.txt -i bicgstab -tol 1e-10 -p sainv -adds true -adds_iter 2

mpirun -n 4 ./test1 bcsstm12.mtx 2 sol.txt hist.txt -i bicgstab -tol 1e-10 -p sainv -adds true -adds_iter 3
```
In this case we can see that the number of iterations required to convergence does not depend monotonically on the choice of `adds_iter`.

```
wget https://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/lanpro/nos1.mtx.gz
gzip -dk nos1.mtx.gz 

mpirun -n 4 ./test1 nos1.mtx 1 sol.txt hist.txt -i cg
mpirun -n 4 ./test1 nos1.mtx 1 sol.txt hist.txt -i cg -p ssor
mpirun -n 4 ./test1 nos1.mtx 1 sol.txt hist.txt -i cg -p ssor -adds true -adds_iter 2
mpirun -n 4 ./test1 nos1.mtx 1 sol.txt hist.txt -i cg -p jacobi -adds true -adds_iter 5
```

## Algebraic Multigrid preconditioner

Algebraic multigrid (AMG) solves linear systems based on multigrid principles, but in a way that only depends on the coefficients in the underlying matrix. The AMG method determines coarse “grids”, inter-grid transfer operators, and coarse-grid equations based solely on the matrix entries.

The AMG preconditioner implemented in LIS is called `SA-AMG` and can be specified as an option by typing `-p saamg`. Unfortunately this preconditioner is incompatible with the use of mpi for multi-processors computations. Therefore, a serial reconfiguration of the LIS installation is needed in order to test the AMG method. 

- Download the file `lis_2.0.34.tar` from webeep and move it to the `shared_folder`
- Enter the docker container as root user
```
docker exec -u root -it hpc-courses /bin/bash
```
- Extract the archive `lis_2.0.34.tar` by typing 
```
tar xvf lis_2.0.34.tar -C /
```
- Load the new LIS module: `module load lis/2.0.34`
- Compile the LIS implementations using the Fortran compiler 
```
gfortran test1.c -I${mkLisInc} -L${mkLisLib} -llis -o test1
```
- Run the examples by typing
```
./test1 testmat2.mtx 2 sol.mtx hist.txt -p saamg -i gmres
```

## Exercise 2

Assess the performances of the SA-AMG preconditioner on the matrices considered in the previous examples.

```
./test1 bcsstm12.mtx 2 sol.mtx hist.txt -i gmres
./test1 bcsstm12.mtx 2 sol.mtx hist.txt -i gmres -p saamg
./test1 bcsstm12.mtx 2 sol.mtx hist.txt -i gmres -p saamg -saamg_unsym true
./test1 bcsstm12.mtx 2 sol.mtx hist.txt -i gmres -p saamg -adds true -adds_iter 3

./test1 nos1.mtx 1 sol.txt hist.txt -i cg
./test1 nos1.mtx 1 sol.txt hist.txt -i cg -p saamg
./test1 nos1.mtx 1 sol.txt hist.txt -i cg -p saamg -saamg_theta 1.0
./test1 nos1.mtx 1 sol.txt hist.txt -i cg -p saamg -adds true

gfortran test2.c -I${mkLisInc} -L${mkLisLib} -llis -o test2
./test2 100 100 1 sol.mtx hist.txt -i cg
./test2 100 100 1 sol.mtx hist.txt -i cg -p saamg
./test2 100 100 1 sol.mtx hist.txt -i bicgstab -p saamg
./test2 100 100 1 sol.mtx hist.txt -i bicgstab -p saamg -adds true
```

## Hand-made Gradient method with Eigen

In the folder `iter_sol++` we create a new file called `grad.hpp` where we write the implementation of the preconditioned Gradient iterative method for solving linear system as a `c++` class template. The preamble and input parameter is the same as the one for the `cg.hpp` function. 

```
namespace LinearAlgebra
{
template <class Matrix, class Vector, class Preconditioner>
int GRAD(const Matrix &A, Vector &x, const Vector &b, const Preconditioner &M,
   int &max_iter, typename Vector::Scalar &tol)
{
  using Real = typename Matrix::Scalar;
  Real   resid;
  Vector q(b.size());
  Vector z(b.size());
  Real   alpha, rho;

  Real   normb = b.norm();
  Vector r = b - A * x;
  if(normb == 0.0) normb = 1;

  if((resid = r.norm() / normb) <= tol)
    {
      tol = resid;
      max_iter = 0;
      return 0;
    }

  for(int i = 1; i <= max_iter; i++)
    {
      z = M.solve(r);
      rho = r.dot(z);
      q = A * z;
      alpha = rho / z.dot(q);

      x += alpha * z;
      r -= alpha * q;

      if((resid = r.norm() / normb) <= tol)
        {
          tol = resid;
          max_iter = i;
          return 0;
        }
    }

  tol = resid;
  return 1;
}
} // namespace LinearAlgebra      
```

### Exercise 3: Test the preconditioned Gradient method

Test the gradient method on a tri-diagonal linear system assuming that the exact solution as all the coefficients equal to 1. Compare the results obtained with the Gradient method and the Gmres method.

```
#include <cstdlib>                      
#include <iostream>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

using std::endl;
using std::cout;

#include "grad.hpp"
#include "gmres.hpp"

int main(int argc, char** argv)
{
  using namespace LinearAlgebra;
  // Some useful alias
  using SpMat=Eigen::SparseMatrix<double>;
  using SpVec=Eigen::VectorXd;

  int n = 400;
  SpMat A(n,n);                       // define matrix
  for (int i=0; i<n; i++) {
      A.coeffRef(i, i) = 2.0*(i+1);
      if(i>0) A.coeffRef(i, i-1) -= i;
      if(i<n-1) A.coeffRef(i, i+1) -= (i+1);
  }

  double tol = 1.e-6;                  // Convergence tolerance
  int result, maxit = 10000;           // Maximum iterations

  std::cout<<"Matrix size:"<<A.rows()<<"X"<<A.cols()<<endl;
  std::cout<<"Non zero entries:"<<A.nonZeros()<<endl;

  // Create Rhs b
  SpVec e = SpVec::Ones(A.rows());
  SpVec b = A*e;
  SpVec x(A.rows());
  Eigen::DiagonalPreconditioner<double> D(A);

  // with Gradient Method
  x=0*x;
  result = GRAD(A, x, b, D, maxit, tol);   
  cout << "iterations performed: " << maxit << endl;
  cout << "tolerance achieved  : " << tol << endl;
  std::cout << "Error norm: "<<(x-e).norm()<< endl;

  return result;
}
```

### Exercise 4: Compare the Gradient method 

Compare the results obtained with the Gradient method and the Gmres method assuming that the exact solution as all the coefficients equal to 1. Consider a non-symmetric matricx of size $1000\times 1000$ constrcted as in Exercise 1 of the exam simulation.

```
#include <cstdlib>                     
#include <iostream>                      
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include "gmres.hpp"          
#include "grad.hpp"

int main(int argc, char** argv)
{
  using namespace LinearAlgebra;
  using SpMat=Eigen::SparseMatrix<double>;
  using SpVec=Eigen::VectorXd;

  int n = 1000;
  SpMat A(n,n);                       
  for (int i=0; i<n; i++) {
      A.coeffRef(i, i) = 2.0*(i+1);
      A.coeffRef(i,n-i-1) = -i;
      if(i>0) A.coeffRef(i, i-1) += i;
      if(i<n-1) A.coeffRef(i, i+1) += i+1;
  }

  double tol = 1.e-12;                 // Convergence tolerance
  int result, maxit = 1000;            // Maximum iterations
  int restart = 50;                    // Restart for gmres
  std::cout << "Matrix size: " << A.rows() << "X" << A.cols() << std::endl;
  std::cout << "Non zero entries: " << A.nonZeros() << std::endl;

  // Create Rhs b
  SpVec e = SpVec::Ones(A.rows());
  SpVec b = A*e;
  SpVec x(A.rows());
  Eigen::DiagonalPreconditioner<double> D(A); // Create diagonal preconditioner

  // Solve with GMRES method with restart
  result = GMRES(A, x, b, D, restart, maxit, tol);
  std::cout << "GMRES with restart "    << std::endl;
  std::cout << "iterations performed: " << maxit << std::endl;
  std::cout << "tolerance achieved  : " << tol << std::endl;
  std::cout << "Error:                " << (x-e).norm()<< std::endl;

  // Solve with GMRES method without restart
  x=0*x; restart = 1000; maxit = 1000; tol = 1.e-12;
  result = GMRES(A, x, b, D, restart, maxit, tol);    
  std::cout << "GMRES without restart " << std::endl;
  std::cout << "iterations performed: " << maxit << std::endl;
  std::cout << "tolerance achieved  : " << tol << std::endl;
  std::cout << "Error norm: "<<(x-e).norm()<< std::endl;

  // Solve with GRADIENT Method
  x=0*x; maxit = 10000; tol = 1.e-6;
  result = GRAD(A, x, b, D, maxit, tol);
  std::cout << "Gradient method "    << std::endl;
  std::cout << "iterations performed: " << maxit << std::endl;
  std::cout << "tolerance achieved  : " << tol << std::endl;
  std::cout << "Error:                " << (x-e).norm()<< std::endl;

  return result;
}
```