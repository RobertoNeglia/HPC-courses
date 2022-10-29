//*****************************************************************
// Iterative template routine -- CGS
//
// CGS solves the unsymmetric linear system Ax = b
// using the Conjugate Gradient Squared method
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//
//*****************************************************************

namespace LinearAlgebra
{
template <class Matrix, class Vector, class Preconditioner>
int CGS(const Matrix &A, Vector &x, const Vector &b, const Preconditioner &M,
    int &max_iter, typename Vector::Scalar &tol)
{
  using Real = typename Vector::Scalar;
  Real   resid(0.0);
  Real   rho_1(0.0), rho_2(0.0), alpha(0.0), beta(0.0);
  Vector p(b.size());
  Vector phat(b.size()), q(b.size()), qhat(b.size());
  Vector vhat(b.size()), u(b.size()), uhat(b.size());

  Real   normb = b.norm();
  Vector r = b - A * x;
  Vector rtilde = r;

  if(normb == 0.0)
    normb = 1.0;

  if((resid = r.norm() / normb) <= tol)
    {
      tol = resid;
      max_iter = 0;
      return 0;
    }

  for(int i = 1; i <= max_iter; i++)
    {
      rho_1 = r.dot(rtilde);
      if(rho_1 == 0.0)
        {
          tol = r.norm() / normb;
          return 2;
        }
      if(i == 1)
        {
          u = r;
          p = u;
        }
      else
        {
          beta = rho_1 / rho_2;
          u = r + beta * q;
          p = u + beta * (q + beta * p);
        }
      phat = M.solve(p);
      vhat = A * phat;
      alpha = rho_1 / vhat.dot(rtilde);
      q = u - alpha * vhat;
      uhat = M.solve(u + q);
      x += alpha * uhat;
      qhat = A * uhat;
      r -= alpha * qhat;
      rho_2 = rho_1;
      if((resid = r.norm() / normb) < tol)
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
