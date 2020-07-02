#include "FunctorExample.hh"

#include "HistTools.hh"
#include "debug.hh"

ClassImp(FunctorExample)
ClassImp(NormSpline)

using namespace std;

// default constructor
FunctorExample::FunctorExample() :
   _npars(0), _func(NULL)
{
}

// explicit constructor
FunctorExample::FunctorExample(const Char_t *Name, Double_t XMin, Double_t XMax, TF1 *fyp, TF1 *fyHe) :
   _func(NULL), _fyp(fyp), _fyHe(fyHe),
   _flux_p(NULL), _flux_p_low(NULL), _flux_He(NULL) 
{
   _nparsp = fyp->GetNpar();
   _nparsHe = fyHe->GetNpar();
   _npars = _nparsp + _nparsHe;

   _func = new TF1(Name, this, XMin, XMax, _npars); 
   _func->SetNpx(500);

   // display R_sum
   _f_R_sum = new TF1(Form("%s_R_sum", Name), this, &FunctorExample::RSum, XMin, XMax, 0);
   _f_R_sum->SetNpx(500);

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

// needed by TF, for reproduction
Double_t FunctorExample::operator()(Double_t *x, Double_t *par)
{
   Double_t xx = x[0];

   Double_t xmin, xmax;
   _flux_p->GetRange(xmin, xmax); 

   _fyp->SetParameters(par);
   _fyHe->SetParameters(par+_nparsp);

   Double_t R_sum = RSum2(x, par); // Rsum w/ the 1.746 term, Rsum2 w/o the 1.746 term. 

   Double_t f;

   if (xx<xmin) f = _flux_p_low->Eval(xx)*_fyp->Eval(xx); 
   if (xx>=xmin) f = _flux_p->Eval(xx)*_fyp->Eval(xx); 

   for (int i = 0; i < (int)_A.size(); ++i)
   {

	if (xx<xmin) f += _flux_elem_low[i]->Eval(xx)*_fyHe->Eval(xx)*(1.+R_sum); 
	if (xx>=xmin) f += _flux_elem[i]->Eval(xx)*_fyHe->Eval(xx)*(1.+R_sum);

        //if (f!=0 && _flux_elem[i]->Eval(xx)!=0) f += _flux_elem[i]->Eval(xx)*_fyHe->Eval(xx)*(1.+R_sum);
   	//if (f==0 && _flux_elem[i]->Eval(xx)!=0) f += _flux_elem[i]->Eval(xx)*_fyHe->Eval(xx); 
   	//if (f==0 && _flux_elem[i]->Eval(xx)==0) f += _flux_He->Eval(xx)*_fyHe->Eval(xx)*R_sum;
   } 

   for (int i = 0; i < (int)_A_above.size(); ++i)
   {
	f += _flux_above[i]->Eval(xx)*_fyHe->Eval(xx); 
   }

   return f;
}

/*
// operator OG
Double_t FunctorExample::operator()(Double_t *x, Double_t *par)
{
   Double_t xx = x[0];

   _fyp->SetParameters(par);
   _fyHe->SetParameters(par+_nparsp);

   Double_t R_sum = RSum(x, par);

   Double_t f = _flux_p->Eval(xx)*_fyp->Eval(xx);
   for (int i = 0; i < (int)_A.size(); ++i)
   {
        if (f!=0 && _flux_elem[i]->Eval(xx)!=0) f += _flux_elem[i]->Eval(xx)*_fyHe->Eval(xx)*(1.+R_sum);
   if (f==0 && _flux_elem[i]->Eval(xx)!=0) f += _flux_elem[i]->Eval(xx)*_fyHe->Eval(xx);
   if (f==0 && _flux_elem[i]->Eval(xx)==0) f += _flux_He->Eval(xx)*_fyHe->Eval(xx)*R_sum;
   }

   return f;
}
*/ 

Double_t FunctorExample::Eval(Double_t x)
{
   return (*this)(&x, _func->GetParameters());
}

// display R_sum
Double_t FunctorExample::RSum(Double_t *x, Double_t *par)
{
   Double_t xx = x[0];

   Double_t R_sum = 0.;

   // starts from He+1
   for (int i = 2; i < (int)_A_ave.size(); ++i)
   {
      R_sum += (_A_ave[i]/4.)*(_flux_elem_ave[i]->Eval(xx)/_flux_elem_ave[1]->Eval(xx));
   }
   R_sum = R_sum + 1.746*(_A_ave[7]/4.)*_flux_elem_ave[7]->Eval(xx)/_flux_elem_ave[1]->Eval(xx); 

   return R_sum;
} 

// estimated F(R,t)/<F(R,t)> ratio, replacing _ave w/ the actual ratio values 
Double_t FunctorExample::RSum2(Double_t *x, Double_t *par)
{
   Double_t xx = x[0];

   Double_t R_sum = 0.;

   // starts from He+1
   for (int i = 2; i < (int)_A_ave.size(); ++i)
   {
      R_sum += (_A_ave[i]/4.)*(_flux_elem_ave[i]->Eval(xx)/_flux_elem_ave[1]->Eval(xx)); 
   }
   // R_sum = R_sum + 1.746*(_A_ave[7]/4.)*_flux_elem_ave[7]->Eval(xx)/_flux_elem_ave[1]->Eval(xx); 

   return R_sum;
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

// Modify Spline

// default constructor
NormSpline::NormSpline() :
   _npars(3), _func(NULL)
{
}

// explicit constructor
NormSpline::NormSpline(const Char_t *Name, Double_t XMin, Double_t XMax, TF1 *fSpline) :
   _npars(3), _func(NULL), _f_spline(fSpline)
{
   _func = new TF1(Name, this, XMin, XMax, _npars);
   _func->SetNpx(1000);
   _func->SetParName(0, "norm");
   _func->SetParName(1, "gb");
   _func->SetParName(2, "y0");

    int nnodes = (_f_spline->GetNpar() - 2) / 2;
    _y0par_index = 2 + nnodes; // gb, ge, { x nodes }
}

// copy constructor
NormSpline::NormSpline(const NormSpline &fsp)
{
   Double_t xmin, xmax;
   fsp._func->GetRange(xmin, xmax);
   _func = new TF1("", this, xmin, xmax, fsp._func->GetNpar());
   _func->SetNpx(fsp._func->GetNpx());
   HistTools::CopyParameters(fsp._func, _func);
}


// destructor
NormSpline::~NormSpline()
{
   if (_func != NULL) delete _func;
}

// (x-nodes, y-parameters)
Double_t NormSpline::operator()(Double_t *x, Double_t *par)
{
   Double_t xx = x[0];

   Double_t norm = par[0];
   Double_t gb = par[1];
   Double_t y0 = par[2];

    _f_spline->SetParameter(0, gb);
    _f_spline->SetParameter(_y0par_index, y0);

   return norm * _f_spline->Eval(xx);
}

Double_t NormSpline::Eval(Double_t x)
{
   return (*this)(&x, _func->GetParameters());
}

void NormSpline::Print()
{
   printf("[%p] NormSpline", this);
   printf(" _npars = %u\n", _npars);

   HistTools::PrintFunction(_f_spline);
   HistTools::PrintFunction(_func);
}
