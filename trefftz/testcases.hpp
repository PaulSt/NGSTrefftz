#include <comp.hpp> // provides FESpace, ...
#include <h1lofe.hpp>
#include <regex>
#include <fem.hpp>
#include <multigrid.hpp>

namespace ngcomp
{

  template <int D> Vec<D + 2> simplesin (Vec<D + 1> p, double wavespeed)
  {
    double x = p[0];
    double t = p[D];
    Vec<D + 2> sol;
    int k = 1;
    if (D == 1)
      {
        sol[0] = sin (k * (wavespeed * t + x));
        sol[1] = k * cos (k * (wavespeed * t + x));
        sol[2] = wavespeed * k * cos (k * (wavespeed * t + x));
      }
    else if (D == 2)
      {
        double y = p[1];
        double sq = sqrt (0.5);
        sol[0] = sin (wavespeed * t + sq * (x + y));
        sol[1] = sq * cos (wavespeed * t + sq * (x + y));
        sol[2] = sq * cos (wavespeed * t + sq * (x + y));
        sol[3] = wavespeed * cos (wavespeed * t + sq * (x + y));
      }
    else if (D == 3)
      {
        double y = p[1];
        double z = p[2];
        double sq = sqrt (1.0 / 3.0);
        sol[0] = sin (wavespeed * t + sq * (x + y + z));
        sol[1] = sq * cos (wavespeed * t + sq * (x + y + z));
        sol[2] = sq * cos (wavespeed * t + sq * (x + y + z));
        sol[3] = sq * cos (wavespeed * t + sq * (x + y + z));
        sol[4] = wavespeed * cos (wavespeed * t + sq * (x + y + z));
      }
    return sol;
  }

  template <int D> Vec<D + 2> gausspw (Vec<D + 1> p, double wavespeed)
  {
    double x = p[0];
    double t = p[D];
    Vec<D + 2> sol;
    int k = 1;
    if (D == 1)
      {
        sol[0] = exp (-100 * ((x - 0.5) * (x - 0.5)));
        sol[1] = -200 * (x - 0.5) * sol[0];
        sol[2] = 0;
      }
    else if (D == 2)
      {
        double y = p[1];
        sol[0] = exp (-100 * ((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)));
        sol[1] = -200 * (x - 0.5) * sol[0];
        sol[2] = -200 * (y - 0.5) * sol[0];
        sol[3] = 0;
      }
    else if (D == 3)
      {
        double y = p[1];
        double z = p[2];
        sol[0] = exp (-100
                      * ((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)
                         + (z - 0.5) * (z - 0.5)));
        sol[1] = -200 * (x - 0.5) * sol[0];
        sol[2] = -200 * (y - 0.5) * sol[0];
        sol[3] = -200 * (z - 0.5) * sol[0];
        sol[4] = 0;
      }
    return sol;
  }

  template <int D> Vec<D + 2> vertgausspw (Vec<D + 1> p, double wavespeed)
  {
    double x = p[0];
    double t = p[D];
    Vec<D + 2> sol;
    int k = 1;
    if (D == 1)
      {
        sol[0] = exp (-100 * ((x - 0.5) * (x - 0.5)));
        sol[1] = -200 * (x - 0.5) * sol[0];
        sol[2] = 0;
      }
    else if (D == 2)
      {
        double y = p[1];
        sol[0] = exp (-100 * ((x - 0.25) * (x - 0.25)));
        sol[1] = -200 * (x - 0.25) * sol[0];
        sol[2] = 0;
        sol[3] = 0;
      }
    else if (D == 3)
      {
        double y = p[1];
        double z = p[2];
        sol[0] = exp (-100 * ((x - 0.25) * (x - 0.25)));
        sol[1] = -200 * (x - 0.25) * sol[0];
        sol[2] = 0;
        sol[3] = 0;
        sol[4] = 0;
      }
    return sol;
  }

  template <int D> Vec<D + 2> standingwave (Vec<D + 1> p, double wavespeed)
  {
    double x = p[0];
    double t = p[D];
    Vec<D + 2> sol;
    int k = 1;
    double sq = sqrt (D);
    if (D == 1)
      {
        sol[0]
            = cos (M_PI * x) * sin (M_PI * t * wavespeed * sq) / (sq * M_PI);
        sol[1] = -sin (M_PI * x) * sin (M_PI * t * wavespeed * sq) / sq;
        sol[2] = cos (M_PI * x) * cos (M_PI * t * wavespeed * sq) * wavespeed;
      }
    else if (D == 2)
      {
        double y = p[1];
        sol[0] = cos (M_PI * x) * cos (M_PI * y)
                 * sin (M_PI * t * wavespeed * sq) / (sq * M_PI);
        sol[1] = -sin (M_PI * x) * cos (M_PI * y)
                 * sin (M_PI * t * wavespeed * sq) / sq;
        sol[2] = -cos (M_PI * x) * sin (M_PI * y)
                 * sin (M_PI * t * wavespeed * sq) / sq;
        sol[3] = cos (M_PI * x) * cos (M_PI * y)
                 * cos (M_PI * t * wavespeed * sq) * wavespeed;
      }
    else if (D == 3)
      {
        double y = p[1];
        double z = p[2];
        sol[0] = cos (M_PI * x) * cos (M_PI * y) * cos (M_PI * z)
                 * sin (M_PI * t * wavespeed * sq) / (sq * M_PI);
        sol[1] = -sin (M_PI * x) * cos (M_PI * y) * cos (M_PI * z)
                 * sin (M_PI * t * wavespeed * sq) / sq;
        sol[2] = -cos (M_PI * x) * sin (M_PI * y) * cos (M_PI * z)
                 * sin (M_PI * t * wavespeed * sq) / sq;
        sol[3] = -cos (M_PI * x) * cos (M_PI * y) * sin (M_PI * z)
                 * sin (M_PI * t * wavespeed * sq) / sq;
        sol[4] = cos (M_PI * x) * cos (M_PI * y) * cos (M_PI * z)
                 * cos (M_PI * t * wavespeed * sq) * wavespeed;
      }
    return sol;
  }

  template <int D>
  Vec<D + 2> TestSolution (Vec<D + 1> p, double wavespeed, char const *solname)
  {
    Vec<D + 2> sol;
    if (strcmp (solname, "gausspw") == 0)
      sol = gausspw<D> (p, wavespeed);
    else if (strcmp (solname, "standingwave") == 0)
      sol = standingwave<D> (p, wavespeed);
    else if (strcmp (solname, "vertgausspw") == 0)
      sol = vertgausspw<D> (p, wavespeed);
    else
      sol = simplesin<D> (p, wavespeed);
    return sol;
  }
}
