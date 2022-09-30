# Numerical linear algebra with Eigen and LIS libraries

Read the documentation available at

- [Eigen user guide](https://eigen.tuxfamily.org/dox/GettingStarted.html),
- [LIS documentation](https://www.ssisc.org/lis/index.en.html).

## 1. Basic linear algebra with Eigen

Eigen is a high-level C++ library of template headers for linear algebra, matrix and vector operations, geometrical transformations, numerical solvers and related algorithms.

### 1.1 A simple example

In the following example we declare a 3-by-3 matrix m which is initialized using the 'Random()' method with random values between -1 and 1. The next line applies a linear mapping such that the values are between 0 and 20. 

The next line of the main function introduces a new type: 'VectorXd'. This represents a (column) vector of arbitrary size. Here, the vector is created to contain 3 coefficients which are left uninitialized. The one but last line uses the so-called comma-initializer and the final line of the program multiplies the matrix m with the vector v and outputs the result.

```
#include <iostream>
#include <Eigen/Dense>
 
using Eigen::MatrixXd;
using Eigen::VectorXd;
 
int main()
{
  MatrixXd m = MatrixXd::Random(3,3);
  m = (m + MatrixXd::Constant(3,3,1.0)) * 10;
  std::cout << "m =" << std::endl << m << std::endl;
  VectorXd v(3);
  v << 1, 0, 0;
  std::cout << "m * v =" << std::endl << m * v << std::endl;
}
```

1. Using VS Code, open the shared folder and create a file `eigen-test1.cpp` with the content of the previous example

2. Change the current directory to the shared folder `cd /home/jellyfish/shared-folder`.

3. In the container, make sure the Eigen module is loaded by typing: `module list`.

4. Compile and run the test.

### 1.2 Linear systems 

We aim to solve a linear system of equations `Ax = b`, where `A` and `b` are matrices (b could be a vector, as a special case). In Eigen we can choose between various decompositions, depending on the properties of your matrix `A`, and depending on the desired speed and accuracy. 

In this example, the `colPivHouseholderQr()` method returns an object of class ColPivHouseholderQR, which is a QR decomposition with column pivoting.

```
#include <iostream>
#include <Eigen/Dense>
 
int main()
{
   Eigen::Matrix3f A;
   Eigen::Vector3f b;
   A << 1,2,3,  4,5,6,  7,8,10;
   b << 3, 3, 4;
   std::cout << "Here is the matrix A:\n" << A << std::endl;
   std::cout << "Here is the vector b:\n" << b << std::endl;
   Eigen::Vector3f x = A.colPivHouseholderQr().solve(b);
   std::cout << "The solution is:\n" << x << std::endl;
}
```

A table of some other decomposition methods available in Eigen is reported [here](https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html). 
If you know more about the properties of your matrix, you can use the above table to select the best method. For example, a good choice for solving linear systems with a non-symmetric matrix of full rank is `PartialPivLU`. If you know that your matrix is also symmetric and positive definite, the above table says that a very good choice is the LLT or LDLT decomposition. 

### 1.3 Compute the error

Only you know what error margin you want to allow for a solution to be considered valid. So Eigen lets you do this computation for yourself, if you want to, as in this example:

```
#include <iostream>
#include <Eigen/Dense>
 
using Eigen::MatrixXd;
 
int main()
{
   MatrixXd A = MatrixXd::Random(100,100);
   MatrixXd b = MatrixXd::Random(100,50);
   MatrixXd x = A.fullPivLu().solve(b);
   double relative_error = (A*x - b).norm() / b.norm(); // norm() is L2 norm
   std::cout << "The relative error is:\n" << relative_error << std::endl;
}
```

### 1.4 Sparse linear systems

The class `SparseMatrix` is the main sparse matrix representation of Eigen's sparse module; it offers high performance and low memory usage. It implements a more versatile variant of the widely-used Compressed Column (or Row) Storage scheme.

The `SparseMatrix` and `SparseVector` classes take three template arguments: the scalar type (e.g., double) the storage order (ColMajor or RowMajor, the default is ColMajor) the inner index type (default is int). As for dense Matrix objects, constructors takes the size of the object. Here are some examples:

```
SparseMatrix<std::complex<float> > mat(1000,2000);   
// declares a 1000x2000 column-major compressed sparse matrix of complex<float>

SparseMatrix<double,RowMajor> mat(1000,2000);              
// declares a 1000x2000 row-major compressed sparse matrix of double

SparseVector<std::complex<float> > vec(1000);              
// declares a column sparse vector of complex<float> of size 1000

SparseVector<double,RowMajor> vec(1000);                   
// declares a row sparse vector of double of size 1000
```

In Eigen, there are several methods available to solve linear systems when the coefficient matrix is sparse. Because of the special representation of this class of matrices, special care should be taken in order to get a good performance. [This page](https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html) lists the sparse solvers available in Eigen. All the available solvers follow the same general concept.

```
#include <Eigen/RequiredModuleName>
// ...
SparseMatrix<double> A;
// fill A
VectorXd b, x;
// fill b
// solve Ax = b
SolverClassName<SparseMatrix<double> > solver;
solver.compute(A);
if(solver.info()!=Success) {
  // decomposition failed
  return;
}
x = solver.solve(b);
if(solver.info()!=Success) {
  // solving failed
  return;
}
// solve for another right hand side:
x1 = solver.solve(b1);
```

A simple example:

```
#include <iostream>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
    SparseMatrix<double> mat(10,10);                              // define matrix
    for (int i=0; i<10; i++) {
        mat.coeffRef(i, i) = 1.0;
    }

    VectorXd b = VectorXd::Constant(mat.rows(), 1);                // define right-hand side

    // Solving 
    SimplicialLDLT<Eigen::SparseMatrix<double> > solver(mat);       // performs LDLT factorization 
    solver.compute(mat);
    if(solver.info()!=Success) {                                    // first sanity check 
        cout << "cannot factorize the matrix" << endl;              // decomposition failed
        return 0;
    }
    
    VectorXd x = solver.solve(b);                                   // solving
    cout << x << endl;                                              // display solution
    return 0;    
}
```

### 1.5 Load and export sparse matrices

To export your matrices and right-hand-side vectors in the matrix-market format, we can use the `unsupported SparseExtra` module.

```
#include <unsupported/Eigen/SparseExtra>
...
Eigen::saveMarket(A, "filename.mtx");
Eigen::saveMarket(A, "filename_SPD.mtx", Eigen::Symmetric);   // if A is symmetric-positive-definite
Eigen::saveMarketVector(B, "filename_b.mtx");
```

To load a matrix in the matrix market format, follow the instructions below:

- in the terminal, use `wget` to download a matrix from the matrix market, e.g. `wget https://math.nist.gov/pub/MatrixMarket2/NEP/mhd/mhd416a.mtx.gz`.

- unzip the file by typing `gzip -dk mhd416a.mtx.gz`

- in Eigen, include the `unsupported SparseExtra` module and use 

```
SparseMatrix<double> mat;
loadMarket(mat, "gr_30_30.mtx");
```

### Exercises

- 1. Compile and run the following example 

```
#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
    // Load matrix
    SparseMatrix<double> mat;
    loadMarket(mat, "mhd416a.mtx");

    VectorXd xe = VectorXd::Constant(mat.rows(), 1);                // define exact solution
    VectorXd b = mat*xe;                                            // compute right-hand side
    cout << b << endl;

    return 0;    
}
```

- 2. Using `wget` download and unzip the matrix [here](https://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/laplace/gr_30_30.mtx.gz). Take as exact solution a vector `xe` defined as in the previous example and compute the right-hand side `b`. Solve the resulting linear system using the sparse solver in Eigen (SparseLU, SparseLDLT, etc..). Then, compute and display the relative error between the exact solution `xe` and the approximated solution.

### Solution:

```
#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
    // Load matrix
    SparseMatrix<double> mat;
    loadMarket(mat, "gr_30_30.mtx");

    VectorXd xe = VectorXd::Constant(mat.rows(), 1);                // define exact solution
    VectorXd b = mat*xe;                                            // compute right-hand side

    // Solving 
    // SimplicialLDLT<Eigen::SparseMatrix<double> > solver(mat);       // performs LDLT factorization of mat
    SparseLU<SparseMatrix<double> > solver;                         // performs LU factorization of mat
    solver.compute(mat);
    if(solver.info()!=Success) {                                    // first sanity check 
        cout << "cannot factorize the matrix" << endl;              // decomposition failed
        return 0;
    }
    
    VectorXd x = solver.solve(b);                                   // use the factorization to solve for the given right-hand side

    if(solver.info()!=Success) {                                    // second sanity check
        cout << "cannot solve the linear system" << endl;           // solving failed
        return 0;
    }

    double relative_error = (x-xe).norm()/(xe).norm();              // compute and print relative error 
    cout << relative_error << endl;

    return 0;    
}
```


## 2. Basic linear algebra with LIS

Lis (Library of Iterative Solvers for linear systems) is a parallel software library for solving discretized linear equations and eigenvalue problems that arise in the numerical solution of partial differential equations using iterative methods. 

As a first example of the usage of LIS we aim to run the first test reported in the user guide. To do so, follow the next steps:

- Download the [latest version](https://www.ssisc.org/lis/dl/lis-2.0.34.zip) of LIS and copy the unzipped folder into you shared folder;

- Move into the test folder: `cd lis-2.0.34/test/`

- Compile the file `test1.c` by typing 

```
mpicc -DUSE_MPI -I${mkLisInc} -L${mkLisLib} -llis test1.c -o test1
```

- Run the test using the command 

```
mpirun -n 4 ./test1 testmat0.mtx 1 sol.txt hist.txt
```

### Exercises: 

- Run in the proper way `test1` in order to solve the second exercise of the previous section

- Following the instructions available on the Lis user-guide, compile and run test2.c