//*****************************************************************
// Iterative template routine -- BiCGSTAB
//
// BiCGSTAB solves the unsymmetric linear system Ax = b
// using the Preconditioned BiConjugate Gradient Stabilized method
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
int BiCGSTAB(const Matrix &A, Vector &x, const Vector &b, const Preconditioner &M,
         int &max_iter, typename Vector::Scalar &tol)
{
  using Real = typename Vector::Scalar;

  Real   resid;
  Real   rho_1(0.), rho_2(0.), alpha(0.), beta(0.), omega(0.);
  Vector p, phat, s, shat, t, v;

  Real   normb = b.norm();
  Vector r = b - A * x;
  Vector rtilde = r;

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
      rho_1 = r.dot(rtilde);
      if(rho_1 == 0)
        {
          tol = r.norm() / normb;
          max_iter = i;
	  return 2;
        }
      if(i == 1)
        p = r;
      else
        {
          beta = (rho_1 / rho_2) * (alpha / omega);
          p = r + beta * (p - omega * v);
        }
      phat = M.solve(p);
      v = A * phat;
      alpha = rho_1 / v.dot(rtilde);
      s = r - alpha * v;
      if((resid = s.norm() / normb) < tol)
        {
          x += alpha * phat;
	  max_iter = i;
          tol = resid;
          return 0;
        }
      shat = M.solve(s);
      t = A * shat;
      omega = s.dot(t) / t.dot(t);
      x += alpha * phat + omega * shat;
      r = s - omega * t;

      rho_2 = rho_1;
      if((resid = r.norm() / normb) < tol)
        {
          tol = resid;
          max_iter = i;
          return 0;
        }
      if(omega == 0)
        {
          tol = r.norm() / normb;
	  max_iter = i;
          return 3;
        }
    }

  tol = resid;
  return 1;
}
} // namespace LinearAlgebra
