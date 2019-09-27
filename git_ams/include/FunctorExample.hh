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

      TF1 *GetTF1Pointer()  { return _func; }

      void Print();

      void SetProtonFlux(TF1 *flux_p) { _flux_p = flux_p; }
      void AddElementFlux(TF1 *flux_el, Double_t A) { _flux_elem.push_back(flux_el); _A.push_back(A); }

   private:
      UInt_t _npars;
      UInt_t _nparsp;
      UInt_t _nparsHe;

      TF1 *_func; //!
      TF1 *_fyp;
      TF1 *_fyHe;

      TF1 *_flux_p;
      std::vector<TF1 *> _flux_elem;
      std::vector<Double_t> _A;

      ClassDef(FunctorExample,1)
};

#endif

// Notes
/* 
TF1 *fyfp = NeutronMonitors::YieldFunctions::CreateFunction("fyfp", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mishev13H1, Particle::PROTON, Energy::RIGIDITY); 
TF1 *fyfHe = NeutronMonitors::YieldFunctions::CreateFunction("fyfHe", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mishev13He4, Particle::HELIUM4, Energy::RIGIDITY); 
FunctorExample *fe = new FunctorExample("fyf", 0.1, 3e3, fyfp, fyfHe)
TF1 *fyf = fe->GetTF1Pointer()
HistTools::PrintFunction(fyfp)
HistTools::PrintFunction(fyfHe)
fe->Print()
fyf->Eval(1)
.qqqqqqqqqqq
TF1 *fyfp = NeutronMonitors::YieldFunctions::CreateFunction("fyfp", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mishev13H1, Particle::PROTON, Energy::RIGIDITY); 
TF1 *fyfHe = NeutronMonitors::YieldFunctions::CreateFunction("fyfHe", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mishev13He4, Particle::HELIUM4, Energy::RIGIDITY); 
FunctorExample *fe = new FunctorExample("fyf", 0.1, 3e3, fyfp, fyfHe)
fyfp->Eval(1)
fyfHe->Eval(1)
TF1 *fyf = fe->GetTF1Pointer()
fyf->Print()
HistTools::PrintFunction(fyf)
fyf->Eval(1)
double yp = fyfp->Eval(100)
double yHe = fyfHe->Eval(100)
fyf->Eval(100)
fyf->Eval(100) - (yp + yHe)
fyfp->SetParameter(12, 3)
yp = fyfp->Eval(100)
fyfp->SetParameter(16, 2)
yp = fyfp->Eval(100)
fyf->Eval(100)
fe->SetProtonFlux(flux_p)
.q
TF1 *fyfp = NeutronMonitors::YieldFunctions::CreateFunction("fyfp", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mishev13H1, Particle::PROTON, Energy::RIGIDITY); 
TF1 *fyfHe = NeutronMonitors::YieldFunctions::CreateFunction("fyfHe", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mishev13He4, Particle::HELIUM4, Energy::RIGIDITY); 
FunctorExample *fe = new FunctorExample("fyf", 0.1, 3e3, fyfp, fyfHe)
TF1 *fyf = fe->GetTF1Pointer()
fe->Eval(10)
.q
TF1 *fyfp = NeutronMonitors::YieldFunctions::CreateFunction("fyfp", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mishev13H1, Particle::PROTON, Energy::RIGIDITY); 
TF1 *fyfHe = NeutronMonitors::YieldFunctions::CreateFunction("fyfHe", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mishev13He4, Particle::HELIUM4, Energy::RIGIDITY); 
FunctorExample *fe = new FunctorExample("fyf", 0.1, 3e3, fyfp, fyfHe)
TF1 *fyf = fe->GetTF1Pointer()
TF1 *fluxp = new TF1("fluxp", "[0]*pow(x,-[1])", 0.1, 3e3)
fluxp->SetParameters(1e4, 2.8)
fluxp->Draw(); gPad->SetLog(); gPad->SetLogy()
fluxp->Draw(); gPad->SetLogx(); gPad->SetLogy()
fyf->Eval(1)
fe->SetProtonFlux(fluxp)
fyf->Eval(1)
*/
