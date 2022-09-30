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

- 2. Solve the resulting linear system using the LIS library. Explore different iterative solvers (Jacobi, Gauss-Seidel, Conjugate Gradient...).

- 3. Set different `options` (tolerance, maximum number of iterations, restart...)


## Exercise 2

- 1. Following the instructions available on the LIS user-guide, compile and run test4.c

- 2. Modify the implementation by changing the size of the linear system to 120 and by setting the conjugate gradient method as iterative solver

- 3. Print the relative error

- 4. Repeat the exercise using the Eigen sparse module. Try both the `SparseLU` and the `ConjugateGradient` solvers


## Exercise 3

- 1. Following the instructions available on the LIS user-guide, compile and run test5.c

- 2. Set n = 100 and test different values of `gamma` and different iterative solvers

- 3. Repeat the exercise using the Eigen `SparseLU` and `BiCGSTAB` solvers


## Exercise 4

- Repeat Exercise 1 on some of the matrices available [here](https://sparse.tamu.edu/?per_page=All)