#ifndef FUNCTOREXAMPLE_H
#define FUNCTOREXAMPLE_H

#include "headers.hh"

#include <vector>

class FunctorExample
{
   public:
      FunctorExample();
      FunctorExample(const Char_t *Name, Double_t XMin, Double_t XMax, TF1 *fyp, TF1 *fyHe);
      FunctorExample(const FunctorExample &);

      virtual ~FunctorExample();

      Double_t operator()(Double_t *x, Double_t *par);
      Double_t Eval(Double_t x);

      // display R_sum
      Double_t RSum(Double_t *x, Double_t *par);
      TF1 *GetRSumTF1Pointer() { return _f_R_sum; }

      TF1 *GetTF1Pointer()  { return _func; }

      void Print();

      void SetProtonFlux(TF1 *flux_p) { _flux_p = flux_p; }
      void SetHeliumFluxReplicate(TF1 *flux_He) { _flux_He = flux_He; }
      void AddElementFlux(TF1 *flux_el, Double_t A) { _flux_elem.push_back(flux_el); _A.push_back(A); }
      void AddElementFlux_Ave(TF1 *flux_el, Double_t A) { _flux_elem_ave.push_back(flux_el); _A_ave.push_back(A); }

   private:
      UInt_t _npars;
      UInt_t _nparsp;
      UInt_t _nparsHe;

      TF1 *_func; //!
      TF1 *_fyp;
      TF1 *_fyHe;
      // TGraph *_R_sum; // R(P), ratio of heavier/helium
      TF1 *_f_R_sum; //!, display R_sum

      TF1 *_flux_p, *_flux_He;
      std::vector<TF1 *> _flux_elem;
      std::vector<TF1 *> _flux_elem_ave; // average flux used to compute R(P)
      std::vector<Double_t> _A;
      std::vector<Double_t> _A_ave; // average A[i] used to compute R(P)

      ClassDef(FunctorExample,1)
};

class NormSpline
{
   public:
      NormSpline();
      NormSpline(const Char_t *Name, Double_t XMin, Double_t XMax, TF1 *fSpline);
      NormSpline(const NormSpline &);

      virtual ~NormSpline();

      Double_t operator()(Double_t *x, Double_t *par);
      Double_t Eval(Double_t x);

      TF1 *GetTF1Pointer()  { return _func; }

      void Print();

   private:
      UInt_t _npars;
      UInt_t _y0par_index;

      TF1 *_func; //!
      TF1 *_f_spline;

      ClassDef(NormSpline, 1)

};

#endif
