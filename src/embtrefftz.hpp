#ifndef FILE_SVDTREFFTZ_HPP
#define FILE_SVDTREFFTZ_HPP
#include <comp.hpp>
#include <python_comp.hpp>
#include <fem.hpp>
#include <integratorcf.hpp>
#include <variant>
#include <bla.hpp>

namespace ngcomp
{
  template <class SCAL>
  std::tuple<shared_ptr<BaseMatrix>, shared_ptr<BaseVector>>
  EmbTrefftz (shared_ptr<SumOfIntegrals> bf, shared_ptr<FESpace> fes,
              shared_ptr<SumOfIntegrals> lf, double eps,
              shared_ptr<FESpace> fes_test, int tndof,
              std::map<std::string, Vector<SCAL>> *stats = nullptr);

  template <class SCAL>
  int LocalTrefftzEmb (Ngs_Element ei, FlatMatrix<SCAL> elmat,
                       FlatMatrix<SCAL, ColMajor> U,
                       FlatMatrix<SCAL, ColMajor> Vt,
                       Array<shared_ptr<BilinearFormIntegrator>> *bfis,
                       shared_ptr<FESpace> fes, double eps,
                       shared_ptr<FESpace> test_fes, int tndof,
                       LocalHeap &mlh);

  class EmbTrefftzFESpace
      : public L2HighOrderFESpace //, public
                                  //std::enable_shared_from_this<EmbTrefftzFESpace>
  {

    // shared_ptr<SumOfIntegrals> bf;
    // shared_ptr<FESpace> fesforet;
    Vector<shared_ptr<Matrix<double>>> ETmats;
    Vector<shared_ptr<Vector<double>>> ETvecs;
    // shared_ptr<BaseVector> uf;
    shared_ptr<VVector<double>> uf;

  public:
    EmbTrefftzFESpace (shared_ptr<MeshAccess> ama, const Flags &flags,
                       bool parseflags = false)
        : L2HighOrderFESpace (ama, flags, parseflags)
    {
      name = "EmbTrefftzFESpace";
      type = "embt";
      needs_transform_vec = true;

      // cout << "DGJUMPS " << this->dgjumps << endl;

      // fesforet = make_shared<L2HighOrderFESpace>(ma,flags,true);
    }

    shared_ptr<BaseVector>
    SetOp (shared_ptr<SumOfIntegrals> bf, shared_ptr<SumOfIntegrals> lf,
           double eps, shared_ptr<FESpace> test_fes, int tndof)
    {
      VVector<double> lfvec (this->GetNDof ());
      static Timer timer ("EmbTrefftz: SetOp");
      RegionTimer reg (timer);
      LocalHeap lh (1000 * 1000 * 1000);

      if (eps == 0 && tndof == 0 && test_fes == nullptr)
        throw Exception ("Need to specify eps, tndof, or test_fes");

      shared_ptr<FESpace> use_as_l2 = dynamic_pointer_cast<FESpace> (
          const_cast<EmbTrefftzFESpace *> (this)->shared_from_this ());
      bool mixed_mode = true;
      if (test_fes == nullptr)
        {
          mixed_mode = false;
          test_fes = use_as_l2;
        }

      Array<shared_ptr<BilinearFormIntegrator>> bfis[4]; // VOL, BND, ...
      for (auto icf : bf->icfs)
        {
          auto &dx = icf->dx;
          bfis[dx.vb] += icf->MakeBilinearFormIntegrator ();
        }

      Array<shared_ptr<LinearFormIntegrator>> lfis[4];
      if (lf)
        {
          ETvecs.SetSize (ma->GetNE (VOL));
          for (auto icf : lf->icfs)
            {
              auto &dx = icf->dx;
              lfis[dx.vb] += icf->MakeLinearFormIntegrator ();
            }
        }

      auto ma = this->GetMeshAccess ();
      ETmats.SetSize (ma->GetNE (VOL));

      ma->IterateElements (VOL, lh, [&] (auto ei, LocalHeap &mlh) {
        bool definedhere = false;
        for (auto icf : bf->icfs)
          {
            if (icf->dx.vb == VOL)
              if ((!icf->dx.definedonelements)
                  || (icf->dx.definedonelements->Test (ei.Nr ())))
                definedhere = true;
          }
        if (!definedhere)
          return; // escape lambda

        Array<DofId> test_dofs;
        test_fes->GetDofNrs (ei, test_dofs);
        Array<DofId> dofs;
        this->GetDofNrs (ei, dofs);

        FlatMatrix<double> elmat (test_dofs.Size (), dofs.Size (), mlh);
        FlatMatrix<double, ColMajor> U (test_dofs.Size (), mlh),
            Vt (dofs.Size (), mlh);

        int nz = LocalTrefftzEmb (ei, elmat, U, Vt, bfis, use_as_l2, eps,
                                  test_fes, tndof, mlh);

        Array<int> dnums;
        this->GetDofNrs (ei, dnums);

        if (lf)
          {
            int nnz = dofs.Size () - nz;
            auto &test_fel = test_fes->GetFE (ei, mlh);
            auto &trafo = ma->GetTrafo (ei, mlh);
            FlatVector<double> elvec (test_dofs.Size (), mlh),
                elveci (test_dofs.Size (), mlh);
            elvec = 0.0;
            for (auto &lfi : lfis[VOL])
              {
                if (lfi->DefinedOnElement (ei.Nr ()))
                  {
                    auto &mapped_trafo = trafo.AddDeformation (
                        lfi->GetDeformation ().get (), mlh);
                    lfi->CalcElementVector (test_fel, mapped_trafo, elveci,
                                            mlh);
                    elvec += elveci;
                  }
              }
            Matrix<double> Ut = Trans (U).Rows (0, nnz);
            Matrix<double> V = Trans (Vt).Cols (0, nnz);
            Matrix<double> SigI (nnz, nnz);

            SigI = static_cast<double> (0.0);
            for (int i = 0; i < nnz; i++)
              SigI (i, i) = 1.0 / elmat (i, i);
            Matrix<double> elinverse = V * SigI * Ut;
            Vector<double> elinv = elinverse * elvec;
            ETvecs[ei.Nr ()] = make_shared<Vector<double>> (elinv);

            for (int k = 0; k < dnums.Size (); k++)
              if (IsRegularDof (dnums[k]))
                lfvec (dnums[k]) = elinv (k);
          }

        // for(int i=0;i<nz;i++) // permute rows such that zeros are on
        // bottom.. is there a better way? Vt.Row(i)=Vt.Row(i+1);
        // Vt.Rows(nz,dofs.Size()) = 0;
        Vt.Rows (0, dofs.Size () - nz) = 0;
        ETmats[ei.Nr ()] = make_shared<Matrix<double>> (Vt);
        // cout << "Vt" << endl;
        // cout <<  Vt << endl;

        // cout << dnums << endl;
        for (int i = 0; i < dofs.Size () - nz; i++)
          free_dofs->Clear (dnums[i]);
      });

      // return make_shared<VVector<double>>(lfvec);
      uf = make_shared<VVector<double>> (lfvec);
      return uf;
    }

    void VTransformMR (ElementId ei, const SliceMatrix<double> mat,
                       TRANSFORM_TYPE type) const
    {
      static Timer timer ("EmbTrefftz: MTransform");
      RegionTimer reg (timer);

      cout << "mtrans " << type << endl;
      // Array<int> dnums, dnums1, dnums2, elnums, fnums, vnums1, vnums2;
      // this->GetDofNrs (ei, dnums1);

      // cout << "=========" << endl;
      // cout << "elnr " << ei.Nr() << " VorB " << ei.VB() << " type " << type
      // << endl; cout << mat.Height() << " x " << mat.Width() << endl; cout <<
      // dnums1.Size() << endl; cout << "=========" << endl; cout << mat <<
      // endl; cout << (*ETmats[ei.Nr()]).Height() << " x " <<
      // (*ETmats[ei.Nr()]).Width() << endl; mat=0;
      if (type == TRANSFORM_MAT_LEFT)
        mat = (*ETmats[ei.Nr ()]) * mat;
      if (type == TRANSFORM_MAT_RIGHT)
        mat = mat * Trans (*ETmats[ei.Nr ()]);
      if (type == TRANSFORM_MAT_LEFT_RIGHT)
        {
          mat = mat * Trans (*ETmats[ei.Nr ()]);
          mat = (*ETmats[ei.Nr ()]) * mat;
        }
      // cout << mat << endl;

      // TRANSFORM_MAT_LEFT = 1,
      // TRANSFORM_MAT_RIGHT = 2,
      // TRANSFORM_MAT_LEFT_RIGHT = 3,
      // TRANSFORM_RHS = 4,
      // TRANSFORM_SOL = 8,
      // TRANSFORM_SOL_INVERSE = 16
    }

    virtual void VTransformVR (ElementId ei, const SliceVector<double> vec,
                               TRANSFORM_TYPE type) const
    {
      static Timer timer ("EmbTrefftz: VTransform");
      RegionTimer reg (timer);

      cout << "vtrans " << type << endl;
      // cout << (*ETmats[ei.Nr()]).Height() << " x " <<
      // (*ETmats[ei.Nr()]).Width() << endl; cout << vec.Size() << endl;
      if (type == TRANSFORM_RHS)
        {
          Vector<double> new_vec (vec.Size ());
          // cout << vec << endl;
          // new_vec = (*ETmats[ei.Nr()]) * vec;
          // vec = new_vec;
        }
      else if (type == TRANSFORM_SOL)
        {
          // cout << vec << endl;
          Vector<double> new_vec (vec.Size ());
          // Vector<double> uf_loc(vec.Size());
          // Array<int> dnums;
          // this->GetDofNrs (ei, dnums);
          // for (int k = 0; k < dnums.Size(); k++)
          // uf_loc(k) = (*uf)(dnums[k]);
          // uf_loc(k) = (*dynamic_pointer_cast<VVector<double>(uf))(dnums[k]);

          // new_vec = Trans(*ETmats[ei.Nr()]) * vec + lfvec.Range();
          // cout << "#######################" << endl;
          // cout << (*ETmats[ei.Nr()]) << endl;
          // cout << vec << endl;
          // new_vec = Trans(*ETmats[ei.Nr()]) * vec;
          // cout << new_vec << endl;
          // new_vec += uf_loc;
          // vec = new_vec;
          // cout << vec << endl;
        }
      else
        cout << "transform " << type << " nothing here" << endl;
    }
  };

}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportEmbTrefftz (py::module m);
#endif // NGS_PYTHON

#endif
