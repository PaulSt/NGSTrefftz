#include "helpers.hpp"

namespace ngfem
{
    void MatToCSR(Matrix<> mat, CSR &sparsemat)
    {
        int spsize = 0;
        sparsemat[0].SetSize(0);
        sparsemat[1].SetSize(0);
        sparsemat[2].SetSize(0);
        for(int i=0;i<mat.Height();i++)
        {
            // guarantee sparsemat[0].Size() to be nbasis,
            // could be shorter but easier for mat*vec
            sparsemat[0].Append(spsize);
            for(int j=0;j<mat.Width();j++)
            {
                if(mat(i,j))
                {
                    spsize++;
                    sparsemat[1].Append(j);
                    sparsemat[2].Append(mat(i,j));
                }
            }
        }
        sparsemat[0].Append(spsize);
    };
}
