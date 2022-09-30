#include <cmath>
#include <vector>

namespace LinearAlgebra
{
template <class Matrix, class Vector, class Preconditioner>
int jacobi(const Matrix &A, Vector &x, const Vector &b, const Preconditioner &D,
   int &max_iter, typename Vector::Scalar &tol)

   {
  using Real = typename Matrix::Scalar;
  Real   resid;
  Vector p(b.size());
  Vector z(b.size());
  Vector q(b.size());
  //Real   alpha, beta, rho;
 // Real   rho_1(0.0);

  Real   normb = b.norm();
  Vector r = b - A * x;

  if(normb == 0.0)
    normb = 1;

  if((resid = r.norm() / normb) <= tol)
    {
      tol = resid;
      max_iter = 0;
      return 0;
    }

  for(int i = 1; i <= max_iter; i++)
    {
      
      if((resid = r.norm() / normb) <= tol)
        {
          x=D.solve(b-A*x)+x;
          resid=(b-A*x).norm();
          tol = resid;
          max_iter = i;

          return 0;
        }

    }

  tol = resid;
  return 1;
}
} // namespace LinearAlgebra