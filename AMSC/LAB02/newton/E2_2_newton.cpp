#include "my_newton.hpp"
#include <iostream>

using std::cin;
using std::cout;
using std::endl;

void
print_vector(const std::vector<double> &v)
{
  for(double d : v)
    cout << d << endl;
}

double
f(double x)
{
  return x * x - 2;
}

double
df(double x)
{
  return 2 * x;
}

int
main(int argc, char **argv)
{
  double              x0;
  std::vector<double> zeroes(2);
  int                 i;
  cout << "Ready to find zeroes of the function f = x^2 - 2..." << endl;
  NewtonSolver n(f, df, 1000, 1e-6, 1e-6);

  cout << "f is a second degree function, so there are two zeroes." << endl;

  i = 0;
  std::vector<double> xs;
  while(i < 2)
    {
      cout << "Please insert an initial guess: ";
      cin >> x0;

      cout << endl << "Ok, starting iterations..." << endl;

      n.solve(x0);

      cout << "Number of iterations: " << n.getIter() << endl;
      cout << "Final residual: " << n.getRes() << endl;
      xs = n.getHistory();
      cout << "The zero is at: " << xs[xs.size() - 1] << endl;
      //   print_vector(xs);

      zeroes.emplace_back(xs[xs.size() - 1]);
      if(i == 1 && zeroes.at(0) == zeroes.at(1))
        i--;
      i++;
    }
}