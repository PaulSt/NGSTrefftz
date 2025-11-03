#include "condensedg.hpp"

namespace ngcomp
{
  void GetSubMatrix (shared_ptr<BaseMatrix> mat, FlatArray<int> drow,
                     FlatArray<int> dcol, FlatMatrix<> out)
  {
    // auto colnr = mat->GetColIndices();
    // for (int i = 0; i < drow.Size (); i++)
    // for (size_t si = mat->First[drow[i]], j = 0; si < mat->First[drow[i] +
    // 1]; si++) if (dcol.Contains (colnr[si])) out (i, j++) = data[si];
    auto smat = dynamic_pointer_cast<SparseMatrix<double>> (mat);
    for (size_t i = 0; i < drow.Size (); i++)
      for (size_t j = 0; j < dcol.Size (); j++)
        out (i, j) = smat->operator() (drow[i], dcol[j]);
  }

  // template <typename T>
  // inline void AInvBt (ngbla::FlatMatrix<T> a, ngbla::FlatMatrix<T> b)
  //{
  // LapackAInvBt (a, b, 'N');
  //}

  // inline void AInvBt (ngbla::FlatMatrix<double> a, ngbla::FlatMatrix<double>
  // b)
  //{
  // ArrayMem<int, 100> p (a.Height ());
  // CalcLU (a, p);
  // SolveTransFromLU (a, p, Trans (b));
  //}

  inline void AInvBt (FlatMatrix<double> a, FlatMatrix<double> b)
  {
    // TODO: improve this
    Matrix<> ainv = a;
    CalcInverse (ainv, INVERSE_LIB::INV_LAPACK);
    Matrix<> c = b * ainv;
    b = c;
  }

  Array<int>
  GetElNeighbours (shared_ptr<ngcomp::MeshAccess> ma, ElementId elnr)
  {
    auto fnums = ma->GetElFacets (elnr);
    Array<int> els;
    for (int i : fnums)
      {
        Array<int> elnums;
        ma->GetFacetElements (i, elnums);
        for (int elnr2 : elnums)
          {
            if (elnr2 != (int)elnr.Nr ())
              els.Append (elnr2);
          }
      }
    return els;
  }

  shared_ptr<BaseMatrix>
  CondenseDG (shared_ptr<BaseMatrix> mat, shared_ptr<BaseVector> vec,
              shared_ptr<FESpace> fes)
  {
    auto smat = dynamic_pointer_cast<SparseMatrix<double>> (mat);
    auto svec = vec; // dynamic_pointer_cast<VVector<double>> (vec);
    if (!smat || !svec)
      throw Exception ("CondenseDG: matrix must be of type "
                       "SparseMatrix and vector must be of type VVector");

    static Timer sc ("CondenseDG");
    RegionTimer reg (sc);
    LocalHeap lh (1000 * 1000 * 1000);

    auto ma = fes->GetMeshAccess ();
    size_t ne = ma->GetNE (VOL);

    // prepare output matrix
    Table<int> table;
    TableCreator<int> creator (ne);
    size_t vndof = 0;
    for (; !creator.Done (); creator++)
      {
        for (auto ei : ma->Elements (VOL))
          {
            Array<DofId> dnums;
            fes->GetDofNrs (ei, dnums, VISIBLE_DOF);
            for (DofId d : dnums)
              if (IsRegularDof (d))
                {
                  creator.Add (ei.Nr (), d);
                  if (creator.GetMode () == 2)
                    vndof++;
                }

            Array<int> els = GetElNeighbours (ma, ei);
            int nels = els.Size ();
            for (int i = 0; i < nels; i++)
              {
                int elnr = els[i];
                dnums.SetSize0 ();
                fes->GetDofNrs (ElementId (elnr), dnums, VISIBLE_DOF);
                for (DofId d : dnums)
                  if (IsRegularDof (d))
                    creator.Add (ei.Nr (), d);

                Array<int> els2 = GetElNeighbours (ma, ElementId (elnr));
                for (auto elnr2 : els2)
                  {
                    if (elnr2 == (int)ei.Nr () || els.Contains (elnr2))
                      continue;
                    els.Append (elnr2);
                    dnums.SetSize0 ();
                    fes->GetDofNrs (ElementId (elnr2), dnums, VISIBLE_DOF);
                    for (DofId d : dnums)
                      if (IsRegularDof (d))
                        creator.Add (ei.Nr (), d);
                  }
              }
          }
      }
    table = creator.MoveTable ();

    auto PP = make_shared<SparseMatrix<double>> (
        vndof, vndof, table, table, false); // TODO: support symmetric
    PP->SetZero ();
    // finished output matrix

    bool elim_only_hidden = true;

    ma->IterateElements (VOL, lh, [&] (auto ei, LocalHeap &mlh) {
      Array<DofId> dofs1;
      fes->GetDofNrs (ei, dofs1);
      Array<int> idofs1 (dofs1.Size (), mlh), odofs1 (dofs1.Size (), mlh);
      idofs1.SetSize0 ();
      odofs1.SetSize0 ();
      auto ctype = elim_only_hidden ? HIDDEN_DOF : CONDENSABLE_DOF;
      for (auto i : dofs1)
        {
          auto ct = fes->GetDofCouplingType (i);
          if (ct & ctype)
            idofs1.AppendHaveMem (i);
          else if (ct != UNUSED_DOF)
            odofs1.AppendHaveMem (i);
        }

      FlatMatrix<> DD (idofs1.Size (), idofs1.Size (), mlh);
      GetSubMatrix (smat, idofs1, idofs1, DD);

      FlatVector<> dd (idofs1.Size (), mlh);
      svec->GetIndirect (idofs1, dd);

      Array<int> els = GetElNeighbours (ma, ei);
      els.Append (ei.Nr ());

      for (auto elnr : els)
        {
          Array<DofId> dofs2;
          fes->GetDofNrs (ngfem::ElementId (elnr), dofs2);
          Array<int> idofs2 (dofs2.Size (), mlh), odofs2 (dofs2.Size (), mlh);
          idofs2.SetSize0 ();
          odofs2.SetSize0 ();
          auto ctype = elim_only_hidden ? HIDDEN_DOF : CONDENSABLE_DOF;
          for (auto i : dofs2)
            {
              auto ct = fes->GetDofCouplingType (i);
              if (ct & ctype)
                idofs2.AppendHaveMem (i);
              else if (ct != UNUSED_DOF)
                odofs2.AppendHaveMem (i);
            }

          FlatMatrix<> AA (odofs1.Size (), odofs2.Size (), mlh);
          GetSubMatrix (smat, odofs1, odofs2, AA);
          PP->AddElementMatrix (odofs1, odofs2, AA, true);

          FlatMatrix<> BB (odofs2.Size (), idofs1.Size (), mlh);
          GetSubMatrix (smat, odofs2, idofs1, BB);
          AInvBt (DD, BB); // b <--- b d^-1 TODO: compute only once Dinv

          FlatVector<> vv (odofs2.Size (), mlh);
          vv = -1.0 * BB * dd;
          svec->AddIndirect (odofs2, vv, true);

          // Array<int> *odofs[2] = { &odofs1, &odofs2 };
          for (auto elnr2 : els)
            {
              Array<DofId> dofs3;
              fes->GetDofNrs (ngfem::ElementId (elnr2), dofs3);
              Array<int> idofs3 (dofs3.Size (), mlh),
                  odofs3 (dofs3.Size (), mlh);
              idofs3.SetSize0 ();
              odofs3.SetSize0 ();
              auto ctype = elim_only_hidden ? HIDDEN_DOF : CONDENSABLE_DOF;
              for (auto i : dofs3)
                {
                  auto ct = fes->GetDofCouplingType (i);
                  if (ct & ctype)
                    idofs3.AppendHaveMem (i);
                  else if (ct != UNUSED_DOF)
                    odofs3.AppendHaveMem (i);
                }

              FlatMatrix<> CC (idofs1.Size (), odofs3.Size (), mlh);
              GetSubMatrix (smat, idofs1, odofs3, CC);

              FlatMatrix<> BDC (odofs2.Size (), odofs3.Size (), mlh);
              BDC = 0.0;
              // GetSubMatrix<double> (smat, odofs2, odofs3, AA);
              SubAB (BB, CC, BDC); // BDC is -b c = -b d^-1 c why no SubABt?
              PP->AddElementMatrix (odofs2, odofs3, BDC, true);
              // smat->AddElementMatrix (odofs2, odofs3, BDC, true);
            }
        }
    });

    //*smat.get() = move(PP);
    // SparseMatrixTM<double> SPP (move(PP));
    // swap (SPP, PP);
    return PP;
  }
}

////////////////////////// python interface ///////////////////////////

#ifdef NGS_PYTHON

void ExportCondenseDG (py::module m)
{
  m.def (
      "CondenseDG",
      [] (shared_ptr<ngcomp::BaseMatrix> mat,
          shared_ptr<ngcomp::BaseVector> vec,
          shared_ptr<ngcomp::FESpace> fes) -> shared_ptr<ngcomp::BaseMatrix> {
        return ngcomp::CondenseDG (mat, vec, fes);
      },
      R"mydelimiter(
      hello
            )mydelimiter",
      py::arg ("mat"), py::arg ("vec"), py::arg ("ncondense"));
}
#endif // NGS_PYTHON
