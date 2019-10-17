#include "FunctorExample.hh"

#include "HistTools.hh"
#include "debug.hh"

ClassImp(FunctorExample)

using namespace std;

// default constructor
FunctorExample::FunctorExample() :
   _npars(0), _func(NULL)
{
}

// explicit constructor
FunctorExample::FunctorExample(const Char_t *Name, Double_t XMin, Double_t XMax, TF1 *fyp, TF1 *fyHe, TH1 *J_sum) :
   _func(NULL), _fyp(fyp), _fyHe(fyHe), _J_sum(J_sum),  
   _flux_p(NULL)
{
   _nparsp = fyp->GetNpar();
   _nparsHe = fyHe->GetNpar();
   _npars = _nparsp + _nparsHe; 

   _func = new TF1(Name, this, XMin, XMax, _npars);
   _func->SetNpx(1000);

   HistTools::CopyParameters(_fyp, _func); // { _fyp_par0, _fyp_par1, ..., _fyp_npars-1, 0, ... }
   HistTools::CopyParameters(_fyHe, _func, _nparsp);// { _fyp_par0, _fyp_par1, ..., _fyp_npars-1, _fyHe_par0, ... }
}


// copy constructor
FunctorExample::FunctorExample(const FunctorExample &fe)
{
   Double_t xmin, xmax;
   fe._func->GetRange(xmin, xmax);
   _func = new TF1("", this, xmin, xmax, fe._func->GetNpar());
   _func->SetNpx(fe._func->GetNpx());
   HistTools::CopyParameters(fe._func, _func);
}


// destructor
FunctorExample::~FunctorExample()
{
   if (_func != NULL) delete _func;
}

/*

// needed by TF, for us 
Double_t FunctorExample::operator()(Double_t *x, Double_t *par)
{
   Double_t xx = x[0]; 

   _fyp->SetParameters(par);
   _fyHe->SetParameters(par+_nparsp);

   Double_t f = _flux_p->Eval(xx)*_fyp->Eval(xx);
   for (int i = 0; i < _A.size(); ++i)
   {
      f += _A[i]/4.*_flux_elem[i]->Eval(xx)*_fyHe->Eval(xx); 
   }

   return f;
} 

*/

// needed by TF, for reproduction 
Double_t FunctorExample::operator()(Double_t *x, Double_t *par) 
{
   Double_t xx = x[0]; 

   _fyp->SetParameters(par);
   _fyHe->SetParameters(par+_nparsp); 

   Double_t f = _flux_p->Eval(xx)*_fyp->Eval(xx); 
   for (int i = 0; i < _A.size(); ++i)
   {
      f += _A[i]/4.*_flux_elem[i]->Eval(xx)*_fyHe->Eval(xx)*(1.+_J_sum->GetBinContent(_J_sum->FindBin(x[0]))); 
   } 

   return f;
}

Double_t FunctorExample::Eval(Double_t x)
{
   return (*this)(&x, _func->GetParameters());
}

void FunctorExample::Print()
{
   printf("[%p] FunctorExample", this);
   printf(" _npars = %u\n", _npars);

   if (_flux_p == NULL) // error
   printf(" Elements = %d\n", (int) _A.size());
   

   HistTools::PrintFunction(_fyp);
   HistTools::PrintFunction(_func);
}
