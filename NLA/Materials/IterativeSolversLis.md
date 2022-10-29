# Iterative solvers with LIS

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
./test1 testmat0.mtx 1 sol.txt hist.txt
```

- To run the code with multiprocessors using mpi, type 

```
mpirun -n 4 ./test1 testmat0.mtx 1 sol.txt hist.txt
```


## Exercise 1

- 1. Using `wget` download and unzip the matrix [here](https://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/laplace/gr_30_30.mtx.gz). Take as exact solution a vector with all the coefficients equal to 1. 

```
wget https://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/laplace/gr_30_30.mtx.gz
gzip gr_30_30.mtx.gz

mpirun -n 4 ./test1 gr_30_30.mtx 2 sol.txt hist.txt -tol 1.0e-14
```

- 2. Solve the resulting linear system using the LIS library. Explore different iterative solvers (Jacobi, Gauss-Seidel, Conjugate Gradient...).

```
mpirun -n 4 ./test1 gr_30_30.mtx 2 sol.txt hist.txt -i cg 

mpirun -n 4 ./test1 gr_30_30.mtx 2 sol.txt hist.txt -i jacobi -maxiter 2000

mpirun -n 4 ./test1 gr_30_30.mtx 2 sol.txt hist.txt -i gs -tol 1.0e-10
```

- 3. Set different `options` (tolerance, maximum number of iterations, restart...)

```
mpirun -n 4 ./test1 gr_30_30.mtx 2 sol.txt hist.txt -i bicgstab -maxiter 100

mpirun -n 4 ./test1 gr_30_30.mtx 2 sol.txt hist.txt -i gmres -restart 20

mpirun -n 4 ./test1 gr_30_30.mtx 2 sol.txt hist.txt -i bicg -p jacobi
```


## Exercise 2

- 1. Following the instructions available on the LIS user-guide, compile and run test4.c

- 2. Modify the implementation by changing the size of the linear system to 120 and by setting the conjugate gradient method as iterative solver

- 3. Print the relative error

```
#include <stdio.h>
#include "lis.h"
LIS_INT main(int argc, char* argv[])
{
    LIS_Comm comm;  
    LIS_MATRIX A;
    LIS_VECTOR b,x,u;
    LIS_SOLVER solver;
    int nprocs,my_rank;
    LIS_INT err,i,n,gn,is,ie,iter;
    LIS_REAL resid;

    n = 120;
    lis_initialize(&argc, &argv);
    comm = LIS_COMM_WORLD;

#ifdef USE_MPI
    MPI_Comm_size(comm,&nprocs);
    MPI_Comm_rank(comm,&my_rank);
#else
    nprocs  = 1;
    my_rank = 0;
#endif

    lis_printf(comm,"\n");
    lis_printf(comm,"number of processes = %d\n",nprocs);

#ifdef _OPENMP
    lis_printf(comm,"max number of threads = %d\n",omp_get_num_procs());
    lis_printf(comm,"number of threads = %d\n",omp_get_max_threads());
#endif

    lis_matrix_create(comm,&A); 
    err = lis_matrix_set_size(A,0,n);
    CHKERR(err);
    lis_matrix_get_size(A,&n,&gn);
    lis_matrix_get_range(A,&is,&ie);
    for(i=is;i<ie;i++)
    {
        if( i>0   )  lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0,A);
        if( i<gn-1 ) lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.0,A);
        lis_matrix_set_value(LIS_INS_VALUE,i,i,2.0,A);
    }
    lis_matrix_set_type(A,LIS_MATRIX_CSR);
    lis_matrix_assemble(A);

    lis_vector_duplicate(A,&u);
    lis_vector_duplicate(A,&b);
    lis_vector_duplicate(A,&x);
    lis_vector_set_all(1.0,u);
    lis_matvec(A,u,b);
    lis_solver_create(&solver);
    lis_solver_set_option("-i cg -p jacobi -tol 1e-14",solver);
    err = lis_solver_set_optionC(solver);
    CHKERR(err);    
    lis_solve(A,b,x,solver);
    lis_solver_get_iter(solver,&iter);
    lis_solver_get_residualnorm(solver,&resid);
    lis_printf(comm,"number of iterations = %D\n",iter);
    lis_printf(comm,"relative residual = %e\n",(double)resid);
    lis_printf(comm,"\n");
    lis_vector_print(x);

    lis_matrix_destroy(A);
    lis_vector_destroy(b);
    lis_vector_destroy(x);
    lis_vector_destroy(u);
    lis_solver_destroy(solver);
    lis_finalize();
    return 0;
}
```

- 4. Repeat the exercise using the Eigen sparse module. Try both the `SparseLU` and the `ConjugateGradient` solvers

```
#include <iostream>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
    int n = 120;	
    SparseMatrix<double> mat(n,n);                           // define matrix
    for (int i=0; i<n; i++) {
        mat.coeffRef(i, i) = 2.0;
	if(i>0) mat.coeffRef(i, i-1) = -1.0;
        if(i<n-1) mat.coeffRef(i, i+1) = -1.0;	
    }

    VectorXd xe = VectorXd::Constant(mat.rows(), 1);         // define sol
    VectorXd b = mat*xe;                                     // compute rhs

    // Solving 
    //SparseLU<SparseMatrix<double> > solver;      
    ConjugateGradient<SparseMatrix<double>, Lower|Upper> solver;
    solver.compute(mat);
    VectorXd x = solver.solve(b);
    std::cout << "#iterations:     " << solver.iterations() << std::endl;
    cout << x << endl;                                       // display sol

    double relative_error = (x-xe).norm()/(xe).norm();       // compute err 
    cout << relative_error << endl;
    return 0;    
}
```

## Exercise 3

- 1. Following the instructions available on the LIS user-guide, compile and run test5.c

- 2. Set n = 100 and test different values of `gamma` and different iterative solvers

```
./test5 100 0.1 -i jacobi

mpirun -n 2 ./test5 100 -1.0 -i gs 

mpirun -n 2 ./test5 100 12.0 -i gmres -restart 30  

mpirun -n 4 ./test5 100 9.9 -i bicgstab -p sainv
```

- 3. Repeat the exercise using the Eigen `BiCGSTAB` solver

```
#include <iostream>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
    int n = 100;
    double gam = 0.1;    
    SparseMatrix<double> mat(n,n);                       // define matrix
    for (int i=0; i<n; i++) {
        mat.coeffRef(i, i) = 2.0;
	if(i>1) mat.coeffRef(i, i-2) = gam;
        if(i<n-1) mat.coeffRef(i, i+1) = 1.0;	
    }

    VectorXd xe = VectorXd::Constant(mat.rows(), 1);     // define sol
    VectorXd b = mat*xe;                                 // compute rhs

    // Solving 
    BiCGSTAB<SparseMatrix<double> > solver;
    solver.compute(mat);
    VectorXd x = solver.solve(b);
    std::cout << "#iterations:   " << solver.iterations() << std::endl;
    cout << x << endl;                                   // display sol

    double relative_error = (x-xe).norm()/(xe).norm();   // compute err 
    cout << relative_error << endl;
    return 0;    
}
```


## Exercise 4

- Repeat Exercise 1 by considering some of the matrices available [here](https://sparse.tamu.edu/?per_page=All)

```
wget https://suitesparse-collection-website.herokuapp.com/MM/HB/bcsstm12.tar.gz
tar -xf bcsstm12.tar.gz
mv bcsstm12/bcsstm12.mtx .
rm -rf bcsstm12 bcsstm12.tar.gz 

mpirun -n 4 ./test1 bcsstm12.mtx 2 sol.txt hist.txt -i bicgstab -maxiter 5000 -tol 1e-11

mpirun -n 4 ./test1 bcsstm12.mtx 2 sol.txt hist.txt -i gmres -p sainv

mpirun -n 4 ./test1 bcsstm12.mtx 1 sol.txt hist.txt -i bicg -p ilu
```

The importance of choosing a proper preconditioning technique can be observed by testing the solvers with and without preconditioners (cf. `-p option`).