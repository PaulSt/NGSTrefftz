#include "helpers.hpp"

namespace ngfem
{

  double Horner (Vector<double> a, int x)
  {
    int deg = a.Size () - 1;
    double result = a[deg];
    for (int i = deg - 1; i >= 0; --i)
      result = result * x + a[i];
    return result;
  }

  double MultiHorner (Vector<double> a, Vector<double> x, int deg) {}

}
