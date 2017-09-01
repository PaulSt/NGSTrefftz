namespace ngfem
{
  constexpr int ipow (int base, int exp, int result = 1)
  {
    return exp < 1 ? result
                   : ipow (base * base, exp / 2,
                           (exp % 2) ? result * base : result);
  }

  constexpr int BinCoeff (int n, int k)
  {
    return round (tgamma (n + 1) / (tgamma (k + 1) * tgamma (n - k + 1)));
  }

}
