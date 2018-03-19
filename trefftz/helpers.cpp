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

  double MultiHorner (Matrix<int> multiind, Vector<double> coeff,
                      Vector<double> point)
  {
    double result = 0;
    Vector<int> k (maxj (multiind));
    cout << k << endl;
    return result;
  }

  Vector<int> maxj (Matrix<int> multiind)
  {
    Vector<int> result (multiind.Height ());
    for (int i = multiind.Height () - 1; i >= 0; i--)
      {
        for (int d = multiind.Width () - 1; d >= 0; d--)
          {
            if (multiind (i, d) - multiind (i - 1, d) != 0)
              {
                result (i) = d + 1;
                break;
              }
          }
      }
    return result;
  }

}
