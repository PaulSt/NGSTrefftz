#include "svdtrefftz.hpp"
//#include "../fem/integratorcf.hpp"
//#include "../fem/h1lofe.hpp"

namespace ngcomp
{
  void LapackSVD (SliceMatrix<> A, SliceMatrix<double, ColMajor> U,
                  SliceMatrix<double, ColMajor> V)
  {
    static Timer t ("LapackSVD");
    RegionTimer reg (t);
    ngbla::integer m = A.Width (), n = A.Height ();
    // Matrix<> U(m), V(n);
    Vector<> S (min (n, m));
    Array<double> work (n * m + 100);
    ngbla::integer info;
    char jobu = 'A', jobv = 'A';
    ngbla::integer lda = A.Dist (), ldu = U.Dist (), ldv = V.Dist ();
    ngbla::integer lwork = work.Size ();

    dgesvd_ (&jobu, &jobv, &m, &n, A.Data (), &lda, S.Data (), U.Data (), &ldu,
             V.Data (), &ldv, work.Data (), &lwork, &info);
    // cout << "info = " << info << endl;
    // if (n <= 100)
    // cout << "S = " << S << endl;
    A.Diag (0) = S;
  }

  shared_ptr<BaseMatrix> SVDTrefftz (shared_ptr<SumOfIntegrals> bf,
                                     shared_ptr<FESpace> fes, double eps)
  {
    LocalHeap lh (1000 * 1000 * 1000);

    auto ma = fes->GetMeshAccess ();
    // cout << *bf;

    // const BitArray & freedofs = *fes->GetFreeDofs();

    Array<shared_ptr<BilinearFormIntegrator>> bfis[4]; // VOL, BND, ...
    // auto gfvec = gf->GetVector().FV<double>();
    // gfvec = 0.0;

    for (auto icf : bf->icfs)
      {
        auto &dx = icf->dx;
        bfis[dx.vb] += make_shared<SymbolicBilinearFormIntegrator> (
            icf->cf, dx.vb, dx.element_vb);
      }

    // Array<FlatMatrix<double> > elmats(ma->GetNE());
    Array<Matrix<double>> elmats;
    Array<Array<int>> newdofs;
    int newdofscounter = 0;
    // Array<int> widths;
    // size_t Tndofs = 0;
    int max_elsperrow = 0;
    int Pwidth = 0;

    ma->IterateElements (
        VOL, lh,
        [&] (auto ei, LocalHeap &mlh)
        // for(auto ei : ma->Elements())
        {
          HeapReset hr (lh);
          // ElementId ei(VOL, el);
          // cout << "ei = " << ei << endl;
          Array<DofId> dofs;
          fes->GetDofNrs (ei, dofs);
          // cout << "dofs = " << dofs << endl;
          auto &trafo = ma->GetTrafo (ei, lh);
          auto &fel = fes->GetFE (ei, lh);

          FlatMatrix<> elmat (dofs.Size (), lh);
          FlatMatrix<> elmati (dofs.Size (), lh);
          elmat = 0.0;
          for (auto &bfi : bfis[VOL])
            {
              bfi->CalcElementMatrix (fel, trafo, elmati, lh);
              elmat += elmati;
              // bfi -> CalcElementMatrixAdd(fel, trafo, elmat, lh);
            }
          // cout << ei << " elmat = " << elmat << endl;
          // cout << elmat << endl;
          Matrix<double, ColMajor> U (dofs.Size (), dofs.Size ()),
              V (dofs.Size (), dofs.Size ());
          // CalcSVD(elmat,U,V);
          LapackSVD (elmat, U, V);
          // int nnz=0;
          Array<int> newdof;
          for (auto sv : elmat.Diag ())
            if (sv < eps)
              newdof.Append (newdofscounter++);
          Pwidth += newdof.Size ();
          newdofs.Append (newdof);
          // max_elsperrow =
          // max_elsperrow<newdof.Size()?newdof.Size():max_elsperrow; cout <<
          // newdof.Size(); cout << nnz <<  "SVD" << elmat.Diag() << endl << U
          // << endl << V << endl; FlatMatrix<> P = U.Cols(nnz,dofs.Size());
          // cout <<"P"<<endl<< P << endl;
          // cout << "U"<< endl<<U <<endl;
          // cout << "V**T"<< endl<<V <<endl;
          // cout << elmat << endl;
          // cout << endl << "U";
          // for(int i=0;i<dofs.Size();i++)
          //{
          // cout << endl;
          // for(int j=0;j<dofs.Size();j++)
          //{
          // cout << U(i,j) << " ";
          // }
          //}
          // cout << endl;
          elmats.Append (U.Cols (dofs.Size () - newdof.Size (), dofs.Size ()));
          // widths.Append(dofs.Size()-nnz);
          // Tndofs += dofs.Size()-nnz;
          //(elmats[ei.Nr()]
        });

    int nnewndof = newdofs.Last ().Last ();

    TableCreator<int> creator (ma->GetNE (VOL));
    for (; !creator.Done (); creator++)
      ma->IterateElements (VOL, lh, [&] (auto ei, LocalHeap &mlh) {
        Array<DofId> dnums;
        fes->GetDofNrs (ei, dnums);
        for (DofId d : dnums)
          if (IsRegularDof (d))
            creator.Add (ei.Nr (), d);
      });
    auto table = creator.MoveTable ();
    cout << table << endl;
    TableCreator<int> creator2 (ma->GetNE (VOL));
    for (; !creator2.Done (); creator2++)
      ma->IterateElements (VOL, lh, [&] (auto ei, LocalHeap &mlh) {
        for (int d : newdofs[ei.Nr ()])
          if (IsRegularDof (d))
            creator2.Add (ei.Nr (), d);
      });
    auto table2 = creator2.MoveTable ();
    cout << table2 << endl;

    // SparseMatrix<double> P(fes->GetNDof(), max_elsperrow);
    SparseMatrix<double> P (fes->GetNDof (), Pwidth, table, table2, false);

    for (auto ei : ma->Elements ())
      {
        Array<DofId> dofs;
        fes->GetDofNrs (ei, dofs);
        P.AddElementMatrix (dofs, newdofs[ei.Nr ()], elmats[ei.Nr ()]);
      }
    return make_shared<SparseMatrix<double>> (P);
  }
}

#ifdef NGS_PYTHON
//#include <python_ngstd.hpp>
//#include <comp.hpp>
//#include <fem.hpp>
// using namespace ngfem;
void ExportSVDTrefftz (py::module m)
{
  m.def (
      "SVDTrefftz",
      [] (shared_ptr<ngfem::SumOfIntegrals> bf,
          shared_ptr<ngcomp::FESpace> fes) -> shared_ptr<ngcomp::BaseMatrix> {
        return SVDTrefftz (bf, fes);
      },
      py::arg ("bf"), py::arg ("fes"));
}
#endif // NGS_PYTHON
