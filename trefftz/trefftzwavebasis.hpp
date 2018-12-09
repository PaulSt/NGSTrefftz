#ifndef FILE_TREFFTZBASIS_HPP
#define FILE_TREFFTZBASIS_HPP

#include <fem.hpp>
#include "helpers.hpp"

namespace ngfem
{
    template<int D>
    void TB_inner(int ord, Matrix<> &trefftzbasis, Vec<D, int> coeffnum, int basis, int dim, int &tracker);

    template<int D>
    const Matrix<>*  TB(int ord);

    template<int D>
    int IndexMap2(Vec<D, int> index, int ord);
}

#endif // FILE_TrefftzBasis_hpp
