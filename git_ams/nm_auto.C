// Neutron Monitor Yield Function Analysis
// By Lingqiang He 09/18/2019

#include "TFile.h"
#include "TTree.h"
#include <set>
#include <array>
#include <string>
#include <FitTools.hh> 
#include <iostream>
#include <fstream>

#include "commonlib/include/Experiments.hh"
#include "commonlib/include/HistTools.hh"
#include "commonlib/include/Spline.hh"
#include "commonlib/include/DateTimeTools.hh"
#include "commonlib/include/debug.hh"
#include "include/FunctorExample.hh" 

using namespace std;

const int n_ams = 8; 
const int n_ele = 24; // number of ACE elements 
const int n_total = 28; // number of total elements 
const int nBRs = 79; 

const char *AMS_Element[n_ele] = { "p", "he", "li", "be", "b", "c", "n", "o", "f" };
const char *ACE_Element[n_ele] = { "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "Va", "Cr", "Mn", "Fe", "Co", "Ni" };
const char *Element[n_total] = { "p", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "Va", "Cr", "Mn", "Fe", "Co", "Ni" };
Particle::Type ACE_Isotope[n_ele] = { Particle::BORON11, Particle::CARBON12, Particle::NITROGEN15, Particle::OXYGEN16, Particle::FLUORINE19, Particle::NEON20, Particle::SODIUM23, Particle::MAGNESIUM24, Particle::ALUMINUM27, Particle::SILICON28, Particle::PHOSPHORUS31, Particle::SULFUR32, Particle::CHLORINE35, Particle::ARGON36, Particle::POTASSIUM41, Particle::CALCIUM40, Particle::SCANDIUM45, Particle::TITANIUM46, Particle::VANADIUM51, Particle::CHROMIUM52, Particle::MANGANESE55, Particle::IRON56, Particle::COBALT59, Particle::NICKEL60 };

Particle::Type compare_isotope[24][10] = {
      { Particle::BORON10, Particle::BORON11 }, // sec, sec
      { Particle::CARBON12,  Particle::CARBON13, Particle::CARBON14 }, // pri, sec, abs
      { Particle::NITROGEN14, Particle::NITROGEN15 }, // mix, sec
      { Particle::OXYGEN16, Particle::OXYGEN17, Particle::OXYGEN18 }, // pri, sec, mix
      { Particle::FLUORINE18, Particle::FLUORINE19 }, // abs, sec 
      { Particle::NEON20, Particle::NEON21, Particle::NEON22 }, // pri, sec, mix
      { Particle::SODIUM22, Particle::SODIUM23, Particle::SODIUM24 }, // abs, mix, abs
      { Particle::MAGNESIUM24, Particle::MAGNESIUM25, Particle::MAGNESIUM26 }, // pri, mix, mix
      { Particle::ALUMINUM26, Particle::ALUMINUM27 }, // sec, mix
      { Particle::SILICON28, Particle::SILICON29, Particle::SILICON30, Particle::SILICON31, Particle::SILICON32 }, // pri, mix, mix, abs, abs
      { Particle::PHOSPHORUS29, Particle::PHOSPHORUS31, Particle::PHOSPHORUS32, Particle::PHOSPHORUS33 }, // abs, mix, abs, abs
      { Particle::SULFUR30, Particle::SULFUR31, Particle::SULFUR32, Particle::SULFUR33, Particle::SULFUR34, Particle::SULFUR35, Particle::SULFUR36 }, // abs, abs, pri, sec, mix, abs, mix
      { Particle::CHLORINE33, Particle::CHLORINE35, Particle::CHLORINE36, Particle::CHLORINE37 }, // abs, sec, sec, sec
      { Particle::ARGON34, Particle::ARGON36, Particle::ARGON37, Particle::ARGON38, Particle::ARGON39, Particle::ARGON40, Particle::ARGON41, Particle::ARGON42 }, // abs, mix, sec, mix, abs, sec, abs, abs
      { Particle::POTASSIUM37, Particle::POTASSIUM39, Particle::POTASSIUM40, Particle::POTASSIUM41 }, // abs, sec, sec, sec
      { Particle::CALCIUM40, Particle::CALCIUM41, Particle::CALCIUM42, Particle::CALCIUM43, Particle::CALCIUM44, Particle::CALCIUM45, Particle::CALCIUM46, Particle::CALCIUM47, Particle::CALCIUM48 }, // pri, sec, sec, sec, sec, abs, sec, abs, sec
      { Particle::SCANDIUM43, Particle::SCANDIUM44, Particle::SCANDIUM45, Particle::SCANDIUM46, Particle::SCANDIUM47, Particle::SCANDIUM48 }, // abs, abs, sec, abs, abs, abs
      { Particle::TITANIUM44, Particle::TITANIUM45, Particle::TITANIUM46, Particle::TITANIUM47, Particle::TITANIUM48, Particle::TITANIUM49, Particle::TITANIUM50 }, // sec, abs, sec, sec, sec, sec, sec
      { Particle::VANADIUM48, Particle::VANADIUM49, Particle::VANADIUM50, Particle::VANADIUM51 }, // abs, sec, sec, sec
      { Particle::CHROMIUM50, Particle::CHROMIUM51, Particle::CHROMIUM52, Particle::CHROMIUM53, Particle::CHROMIUM54 }, // sec, sec, mix, sec, sec
      { Particle::MANGANESE52, Particle::MANGANESE53, Particle::MANGANESE54, Particle::MANGANESE55 }, // abs, sec, sec, mix
      { Particle::IRON54, Particle::IRON55, Particle::IRON56, Particle::IRON57, Particle::IRON58, Particle::IRON59, Particle::IRON60 }, // mix, sec, pri, pri, pri, abs, abs
      { Particle::COBALT56, Particle::COBALT57, Particle::COBALT58, Particle::COBALT59, Particle::COBALT60 }, // abs, sec, abs, pri, abs
      { Particle::NICKEL56, Particle::NICKEL58, Particle::NICKEL59, Particle::NICKEL60, Particle::NICKEL61, Particle::NICKEL62, Particle::NICKEL63, Particle::NICKEL64 } // abs, pri, sec, pri, pri, pri, abs, pri
};

const char *Element_Type[n_total] = {

	"pri", // proton
	"pri", // He
	"sec", // Li
	"sec", // Be
	"sec", // B
	"pri", // C
	"mix", // N
	"pri", // O
	"sec", // F
	"mix", // Ne
	"mix", // Na
	"mix", // Mg
	"mix", // Al
	"pri", // Si
	"mix", // P
	"pri", // S
	"sec", // Cl
	"mix", // Ar
	"sec", // K
	"mix", // Ca
	"sec", // Sc
	"sec", // Ti
	"sec", // Va
	"mix", // Cr
	"mix", // Mn
	"pri", // Fe
	"mix", // Co
	"pri"  // Ni

}; 

int i_temp[n_total] = {

	1, // proton
	1, // He
	0, // Li
	0, // Be
	0, // B
	1, // C
	2, // N
	1, // O
	0, // F
	2, // Ne
	2, // Na
	2, // Mg
	2, // Al
	1, // Si
	2, // P
	2, // S
	0, // Cl
	2, // Ar
	0, // K
	2, // Ca
	0, // Sc
	0, // Ti
	0, // Va
	2, // Cr
	2, // Mn
	1, // Fe
	2, // Co
	1  // Ni

}; 

const Double_t A[n_total] = {

      1.01961, // H
      3.85, // He
      6.52381, // Li
      7.93, // Be
      10.6897, // B
      12.0654, // C
      14.55, // N
      16.36, // O
      19, // F
      20.7, // Ne
      23, //Na
      24.58, // Mg
      26.95, // Al
      28.25, // Si
      31., // P
      32.45, // S
      35.75, // Cl
      37.3, // Ar
      40.3, // K
      41.46, // Ca
      45., // Sc
      47.07, // Ti
      50.2, // V
      52.22, // Cr
      54.22, // Mn
      55.89, // Fe
      58.34, // Co
      59.02, // Ni

   };  

string name_isotope[24][10] = {

      { "B10", "B11" },
      { "C12", "C13", "C14" },
      { "N14", "N15" },
      { "O16", "O17", "O18" },
      { "F18", "F19" },
      { "Ne20", "Ne21", "Ne22" },
      { "Na22", "Na23", "Na24" },
      { "Mg24", "Mg25", "Mg26" },
      { "Al26", "Al27" },
      { "Si28", "Si29", "Si30", "Si31", "Si32" },
      { "P29", "P31", "P32", "P33" },
      { "S30", "S31", "S32", "S33", "S34", "S35", "S36" },
      { "Cl33", "Cl35", "Cl36", "Cl37" },
      { "Ar34", "Ar36", "Ar37", "Ar38", "Ar39", "Ar40", "Ar41", "Ar42" },
      { "K37", "K39", "K40", "K41" },
      { "Ca40", "Ca41", "Ca42", "Ca43", "Ca44", "Ca45", "Ca46", "Ca47", "Ca48" },
      { "Sc43", "Sc44", "Sc45", "Sc46", "Sc47", "Sc48" },
      { "Ti44", "Ti45", "Ti46", "Ti47", "Ti48", "Ti49", "Ti50" },
      { "V48", "V49", "V50", "V51" },
      { "Cr50", "Cr51", "Cr52", "Cr53", "Cr54" },
      { "Mn52", "Mn53", "Mn54", "Mn55" },
      { "Fe54", "Fe55", "Fe56", "Fe57", "Fe58", "Fe59", "Fe60" },
      { "Co56", "Co57", "Co58", "Co59", "Co60" },
      { "Ni56", "Ni58", "Ni59", "Ni60", "Ni61", "Ni62", "Ni63", "Ni64" },
};

int isotope_size[24] = { 2, 3, 2, 3, 2, 3, 3, 3, 2, 5, 4, 7, 4, 8, 4, 9, 6, 7, 4, 5, 4, 7, 5, 8 }; // change when compare_isotope[][] is changed 

// list of functions 
UInt_t UTimeToBR(Long64_t utime);
UInt_t UBRToTime(int BR);
double *get_EMed(const char *element);
double *get_kin_bins(const char *element);
double *get_spall_corr(const char *element);
double *get_spall_corr_unc(const char *element);
double compare_sig(TF1 *fit1, TF1 *fit2); // obtain sigma value between two fits

TGraphAsymmErrors *get_ace_graph(const char *element, UInt_t iBin, UInt_t nBRs); // flux in Kinetic Energy over time 
TGraphAsymmErrors *get_ace_average_graph(const char *element, UInt_t *BRs, UInt_t nBRs); // flux in Kinetic over energy bins 

void nm_auto(); 
void nm_reproduce(); // reproduce the results from Koldobisky
TH1 *heavier_he(); // return ratio heavier/helium

// #### main function ####
void nm_auto(){

	// fmul = f_fit[i] X fyf_pHE 

	gSystem->mkdir("data/nm/NMYF", true); 

	TF1 *f_fit[n_total]; 

	TFile *fin = new TFile("data/ACE/contribute2/h_contribute.root"); 

	TCanvas *c1 = new TCanvas("c1", "NMYF", 900, 2700); 
	c1->Divide(1, 3); 

	TF1 *fyfp = NeutronMonitors::YieldFunctions::CreateFunction("fyfp", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mishev13H1, Particle::PROTON, Energy::RIGIDITY);
	TF1 *fyfHe = NeutronMonitors::YieldFunctions::CreateFunction("fyfHe", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mishev13He4, Particle::HELIUM4, Energy::RIGIDITY); 

	for (int i=0; i<n_total; i++){

		int nnodes = i<4?9:7; 

		TF1 *fyf_pHe = i==0 ? (TF1 *) fyfp->Clone("fyf_p") : (TF1 *) fyfHe->Clone("fyf_He"); // p or He 

		TFile *file = new TFile(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[i_temp[i]], nnodes)); // load ACE Template combined fit 		

		if (i<8) { 

			Spline *sp_comb = new Spline(Form("fit_%s", Element[i]), nnodes, Spline::LogLog | Spline::PowerLaw); 
			f_fit[i] = sp_comb->GetTF1Pointer();  // real function 
			TF1 *fit_comb = (TF1*) fin->Get(Form("fit_%s", Element[i]))->Clone(Form("fit_%s", Element[i]));
	
			HistTools::CopyParameters(fit_comb, f_fit[i]); 
			double x1, x2;
			fit_comb->GetRange(x1,x2); 
			f_fit[i]->SetRange(x1,x2); 

		} else if (i>=8) {

			int nnodes = 7;
	
			TF1 *f_extract = (TF1 *) fin->Get(Form("fit_%s", Element[i]))->Clone(Form("f_%s_extract", Element[i])); // extract the multiply const 	
			Spline *sp_comb = new Spline("sp_comb", nnodes, Spline::LogLog | Spline::PowerLaw); 
			f_fit[i] = sp_comb->GetTF1Pointer();  // real function 
			TF1 *fit_comb = (TF1*) file->Get("fit_both")->Clone("fit_both"); 
	
			HistTools::CopyParameters(fit_comb, f_fit[i]); 
			double x1, x2;
			fit_comb->GetRange(x1,x2); 
			f_fit[i]->SetRange(x1,x2); 

			f_fit[i] = HistTools::CombineTF1Const(f_fit[i], f_extract->GetParameter(1), HistTools::MultiplyConst, Form("fit_%s", Element[i])); 

		}

		HistTools::PrintFunction(f_fit[i]); 

		TF1 *fmul = HistTools::CombineTF1(f_fit[i], fyf_pHe, HistTools::Multiply, "fmul");

		c1->cd(1);
		fyf_pHe->Draw(); 
		gPad->SetLogx(); 
		gPad->SetLogy();
		fyf_pHe->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV"));

		c1->cd(2);
		f_fit[i]->Draw(); 
		gPad->SetLogx(); 
		gPad->SetLogy(); 
		f_fit[i]->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV"));
 
		c1->cd(3); 
		fmul->Draw(); 
		gPad->SetLogx(); 
		//gPad->SetLogy();
		fmul->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV"));

		// fyf_pHe->Eval(1e3);
		// fmul->Integral(0.1, 1e4); 
		// fmul->Integral(0.1, 1e4)/fmul->Integral(0.1, 3e3); 
		HistTools::CombineTF1Const(fmul, A[i]/1, HistTools::MultiplyConst, "fmulconst"); 
		// HistTools::PrintFunction(fmul); 

		c1->Print(Form("data/nm/NMYF/fmul_%s.png", Element[i])); 

	}

} 

// Reconstruct the NM count with a simple cosmic ray contribution model from Koldobisky's assumption 
void nm_reproduce(){

	int nnodes_ams = 6; 	

	TGraph *N_t = new TGraph(); 

	Experiments::DataPath = "data"; 

	TH1 *J_sum = heavier_he(); 

	J_sum->Print("range");  

	const int nBRs = Experiments::Info[Experiments::AMS02].Dataset[1].nMeasurements; 

	TF1 *fyfp = NeutronMonitors::YieldFunctions::CreateFunction("fyfp", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mishev13H1, Particle::PROTON, Energy::RIGIDITY);
	TF1 *fyfhe = NeutronMonitors::YieldFunctions::CreateFunction("fyfhe", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mishev13He4, Particle::HELIUM4, Energy::RIGIDITY); 

	//TF1 *fyfp = NeutronMonitors::YieldFunctions::CreateFunction("fyfp", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mangeard16H1, Particle::PROTON, Energy::RIGIDITY);
	//TF1 *fyfhe = NeutronMonitors::YieldFunctions::CreateFunction("fyfhe", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mangeard16He4, Particle::HELIUM4, Energy::RIGIDITY); 

	HistTools::PrintFunction(fyfp);
	HistTools::PrintFunction(fyfhe); 

	TCanvas *c1 = new TCanvas("c1", "Estimated NM Count Rate", 2700, 900); 

	TF1 *f_BR_p[nBRs];   
	TF1 *f_BR_he[nBRs];  

	TFile *fit_result = new TFile(Form("data/amsfit/fit_result_node%d.root", nnodes_ams));

	for (int i=0; i<nBRs; i++){

		Spline *sp_p = new Spline(Form("f_BR_p_%d", i), nnodes_ams, Spline::LogLog | Spline::PowerLaw); 
		f_BR_p[i] = sp_p->GetTF1Pointer();  // real function 
		TF1 *fit_p = (TF1*) fit_result->Get(Form("fsp_BR_p_%02d", i))->Clone(Form("f_BR_p_%d", i)); 
	
		HistTools::CopyParameters(fit_p, f_BR_p[i]); 
		double x1, x2;
		fit_p->GetRange(x1,x2); 
		f_BR_p[i]->SetRange(x1,x2); 

		Spline *sp_he = new Spline(Form("f_BR_he_%d", i), nnodes_ams, Spline::LogLog | Spline::PowerLaw); 
		f_BR_he[i] = sp_he->GetTF1Pointer();  // real function 
		TF1 *fit_he = (TF1*) fit_result->Get(Form("fsp_BR_he_%02d", i))->Clone(Form("f_BR_he_%d", i)); 
	
		HistTools::CopyParameters(fit_he, f_BR_he[i]); 
		double xx1, xx2; 
		fit_he->GetRange(xx1,xx2); 
		f_BR_he[i]->SetRange(xx1, xx2); 	

		FunctorExample *fe = new FunctorExample("fe", 0.1, 3e3, fyfp, fyfhe, J_sum); 

		fe->SetProtonFlux(f_BR_p[i]); 
		fe->AddElementFlux(f_BR_he[i], A[1]); 

		TF1 *f = fe->GetTF1Pointer(); 

		// fe->Print(); 

		double f_check = 0.; 
		double *bl = HistTools::BuildLogBins(0.1, 3e3, 10); // check the integral by separation into few integrals 
			
		for (int k=0; k<10; k++){ 
		
			f_check += f->Integral(bl[k], bl[k+1]);  
			if (i==0) printf("k = %d, lower limit = %4.1f, upper limit = %4.1f, f = %10.4f \n", k, bl[k], bl[k+1], f->Integral(0.1, bl[k]));  
				 
		}		

		printf("Time = %d, N_t = %10.4f,  f-f_check = %10.4f (good if =0) \n", UBRToTime(i+2426), f->Integral(0.1, 3e3), f->Integral(0.1, 3e3)-f_check); 

		N_t->SetPoint(i, UBRToTime(i+2426), f->Integral(0.1, 3e3)); // this internally makes a loop in the range Rmin to Rmax, and calls FunctorExample::operator() at every step, computing the integral as the sum of all the steps  

		// break; 

	} 

	c1->cd(1); 
	HistTools::SetStyle(N_t, kBlue, kFullCircle, 0.9, 1, 1); 

	N_t->GetXaxis()->SetTimeDisplay(1);
  	N_t->GetXaxis()->SetTimeFormat("%m-%y");
	N_t->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	N_t->SetTitle("; ; Estimated NM Count Rate"); 
	N_t->Draw("APL"); 

	double sum_N = 0., average_N = 0.;  
	for (int i=0; i<nBRs; i++){
		double x_N, y_N; 
		N_t->GetPoint(i, x_N, y_N); 
		sum_N += y_N; 
	}
	average_N = sum_N/nBRs; 
	TGraph *N_norm_t = new TGraph(); 
	for (int i=0; i<nBRs; i++){
		double x_N, y_N; 
		N_t->GetPoint(i, x_N, y_N); 
		N_norm_t->SetPoint(i, UBRToTime(i+2426), y_N/average_N);  
	}
	// PRINT_GRAPH(N_t); 

	TCanvas *c2 = new TCanvas("c2", "Estimated NM Count Rate (Normalized)", 2700, 900); 
	c2->cd(1); 

	TLegend *legend2 = new TLegend(0.1,0.8,0.28,0.9); // left, down, right, top 
	legend2->AddEntry(N_norm_t, "Mi13", "l"); 

	HistTools::SetStyle(N_norm_t, kBlue, kFullCircle, 0.9, 1, 1); 

	N_norm_t->GetXaxis()->SetTimeDisplay(1);
  	N_norm_t->GetXaxis()->SetTimeFormat("%m-%y");
	N_norm_t->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	N_norm_t->SetTitle("; ; Normalized Estimated NM Count Rate"); 

	N_norm_t->Draw("APL"); 
	legend2->Draw("SAME"); 

	c2->Print("data/nm/reproduce/estimated_nm_count.png"); 

}

// Return ratio heavier/helium
TH1 *heavier_he(){

	int nnodes_ams = 6; 

	gSystem->mkdir("data/nm/reproduce", true); 

	TF1 *f_fit[n_total]; 

	TFile *fin = new TFile("data/ACE/contribute2/h_contribute.root");  

	TF1 *fyfp = NeutronMonitors::YieldFunctions::CreateFunction("fyfp", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mishev13H1, Particle::PROTON, Energy::RIGIDITY);
	TF1 *fyfHe = NeutronMonitors::YieldFunctions::CreateFunction("fyfHe", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mishev13He4, Particle::HELIUM4, Energy::RIGIDITY);

	TFile *file0 = new TFile(Form("data/amsfit/fit_result_node%d.root", nnodes_ams)); 

	TH1 *J_int[n_ams]; // integrated ams flux 
	TH1 *J_ratio[n_ams]; // A[i] * J_i / J_He / 4.  
	TH1 *J_sum; // R(P), summed ratio as function of rigidity   

	// 0,  1,  2,  3, 4, 5, 6, 7 
	// p, he, li, be, b, c, n, o
	for (int i=0; i<n_ams; i++){
		
		J_int[i] = (TH1*) file0->Get(Form("h_%s", AMS_Element[1]))->Clone(Form("h_%s", AMS_Element[1]));
		
		Spline *sp_ams = new Spline("sp_ams", nnodes_ams, Spline::LogLog | Spline::PowerLaw); 
		f_fit[i] = sp_ams->GetTF1Pointer();  // real function 
		TF1 *fit_ams = (TF1*) file0->Get(Form("fsp_%s", AMS_Element[i]))->Clone(Form("f_%s", AMS_Element[i])); 
	
		HistTools::CopyParameters(fit_ams, f_fit[i]); 
		double x1, x2;
		fit_ams->GetRange(x1,x2); 
		f_fit[i]->SetRange(x1,x2); 
 
		for (int j=1; j<=J_int[i]->GetNbinsX(); j++){ 

      			Double_t R = J_int[i]->GetBinLowEdge(j);
     			Double_t w = J_int[i]->GetBinWidth(j);
      			Double_t flux = f_fit[i]->Integral(R, R+w)/w; 
			J_int[i]->SetBinContent(j, flux);
			J_int[i]->SetBinError(j, 0.); 
 
		} 
		if (i==1) {
			J_sum = (TH1*) file0->Get(Form("h_%s", AMS_Element[i]))->Clone(Form("h_%s", AMS_Element[i])); 
			// J_sum->Print("range"); 
			for (int j=1; j<=J_int[1]->GetNbinsX(); j++){
				J_sum->SetBinContent(j, 0.);
				J_sum->SetBinError(j, 0.);  
			}
		} 
		//J_int[i]->Print("range");  
	}

	// printf("Helium Bins = %d \n", J_int[1]->GetNbinsX()); 

	for (int i=2; i<n_ams; i++){
		
		J_ratio[i] = (TH1*) J_int[1]->Clone(Form("h_%s", AMS_Element[i])); 

		// J_ratio[i]->Print("range"); 

		for (int j=1; j<=J_int[1]->GetNbinsX(); j++){
			Double_t dM = J_int[i]->GetBinError(j), M = J_int[i]->GetBinContent(j); 
			Double_t dN = J_int[1]->GetBinError(j), N = J_int[1]->GetBinContent(j); 
			J_ratio[i]->SetBinContent(j, (A[i]/4.)*(M/N));			 
			J_ratio[i]->SetBinError(j, (A[i]/4.)*(M/N)*sqrt((dM/M)*(dM/M)+(dN/N)*(dN/N))); 
		}	 

		for (int j=1; j<=J_int[1]->GetNbinsX(); j++){ 
			J_sum->AddBinContent(j, J_ratio[i]->GetBinContent(j)); 
			// printf("J_ratio Element %s = %10.4f \n", Element[i], J_ratio[i]->GetBinContent(j)); 
		}
	} 

	// Z>8 Elements  
	for (int j=1; j<=J_int[1]->GetNbinsX(); j++){ 

			Double_t R = J_int[1]->GetBinLowEdge(j);
     			Double_t w = J_int[1]->GetBinWidth(j);
      			Double_t flux = f_fit[7]->Integral(R, R+w)/w; 

			J_sum->AddBinContent(j, (1.746*A[7]/4.)*(flux/J_int[1]->GetBinContent(j))); 
		}

	J_sum->Print("range"); 

	TCanvas *c1 = new TCanvas("c1", "NM", 2700, 900); 
	// c1->Divide(1, 3);
	c1->cd(1);

	gPad->SetLogx(); 
	//gPad->SetLogy();
	gPad->SetGrid();

	HistTools::SetStyle(J_sum, kBlue, kFullCircle, 0.9, 1, 1);

	J_sum->GetYaxis()->SetRangeUser(0.3, 0.6);
	J_sum->GetXaxis()->SetRangeUser(0.01, 1210);
	J_sum->SetTitle("summed ratio as function of rigidity"); 
	J_sum->SetYTitle("ratio heavier/helium"); 
	J_sum->Draw("P"); 

	c1->Print("data/nm/reproduce/ratio_heavier_helium.png"); 

	return J_sum; 
}

UInt_t UTimeToBR(Long64_t utime){
	Long64_t first_BR = -4351622400; // Unix time corresponding to the first Bartels rotation: Feb. 8th, 1832
	return (utime - first_BR) / (27*86400) + 1;
}

UInt_t UBRToTime(int BR){
	Long64_t first_BR = -4351622400; // Unix time corresponding to the first Bartels rotation: Feb. 8th, 1832 
	return (BR - 1) * (27*86400) + first_BR;
}

// flux over time (might be incorrect) 
TGraphAsymmErrors *get_ace_graph(const char *element, UInt_t iBin, UInt_t nBRs){

   double *kin_bins = get_kin_bins(element);
   double *SpallCorr = get_spall_corr(element);
   double *SpallCorrUnc = get_spall_corr_unc(element);
   //double *EMed = get_EMed(element);

   TFile *_file2 = new TFile(Form("data/ACE/%s_97226_18359.root", element));
   TTree *ace=(TTree*)_file2->Get("ace");

   float F[7], C[7]; // must be initialized 
   Long64_t utime;
   Float_t livetime;
	
   ace->SetBranchAddress("F", F);
   ace->SetBranchAddress("C", C);
   //ace->SetBranchAddress("start_utime", &utime);	
   //ace->SetBranchAddress("livetime", &livetime);

   TGraphAsymmErrors *graph = new TGraphAsymmErrors(nBRs); 

   for (int iBR=2426; iBR<=2493; iBR++){

	if (iBR != 2472 && iBR != 2473){

		ace->GetEntry(iBR-2426); 

      		double sys_err = F[iBin] * sqrt(8e-4 + SpallCorrUnc[iBin]*SpallCorrUnc[iBin]/SpallCorr[iBin]/SpallCorr[iBin]);
		double stat_err = F[iBin]/sqrt(C[iBin]); 
		double err = sqrt(stat_err*stat_err + sys_err*sys_err);

		graph->SetPoint(iBR-2426, UBRToTime(iBR), F[iBin]); 
		graph->SetPointError(iBR-2426, 0., 0., err, err); 
	}

   }

   graph->Set(65);

   _file2->Close();
   
   return graph; 
}

// flux over energy bins
TGraphAsymmErrors *get_ace_average_graph(const char *element, UInt_t *BRs, UInt_t nBRs){

   double *kin_bins = get_kin_bins(element);
   double *SpallCorr = get_spall_corr(element);
   double *SpallCorrUnc = get_spall_corr_unc(element);
   double *EMed = get_EMed(element);

   TFile *_file2 = new TFile(Form("data/ACE/%s_97226_18359.root", element));
   TTree *ace=(TTree*)_file2->Get("ace");

   float F[7], C[7]; // must be initialized 
   Long64_t utime;
   Float_t livetime;
	
   ace->SetBranchAddress("F", F);
   ace->SetBranchAddress("C", C);
   ace->SetBranchAddress("start_utime", &utime);	
   ace->SetBranchAddress("livetime", &livetime);

   const UShort_t nbins = 7;
   TGraphAsymmErrors *graph = new TGraphAsymmErrors(nbins); 

   for (UShort_t ibin = 0; ibin < nbins; ++ibin)
   {
      		graph->SetPoint(ibin, EMed[ibin], 0.);
      		graph->SetPointError(ibin, EMed[ibin] - kin_bins[2*ibin], kin_bins[2*ibin+1] - EMed[ibin], 0., 0.);
   }

   Double_t dT = 0.;
   for (UInt_t ibr = 0; ibr < nBRs; ++ibr)
   {
      if (BRs[ibr] >= ace->GetEntries()) continue; 

      ace->GetEntry(BRs[ibr]);
      dT += livetime;

      for (UShort_t ibin = 0; ibin < nbins; ++ibin)
      {
         	graph->GetY()[ibin] += livetime*F[ibin];
         	graph->GetEYlow()[ibin] += C[ibin];
      }
   }
   for (UShort_t ibin = 0; ibin < nbins; ++ibin)
   {
      		Double_t flux   = graph->GetY()[ibin]/dT;
      		Double_t counts = graph->GetEYlow()[ibin];
      		Double_t syst   = flux * sqrt(8e-4 + SpallCorrUnc[ibin]*SpallCorrUnc[ibin]/SpallCorr[ibin]/SpallCorr[ibin]);
      		Double_t stat   = flux/sqrt(counts);
      		Double_t err    = sqrt(stat*stat + syst*syst);

      		graph->GetY()[ibin]      = flux;
      		graph->GetEYlow()[ibin]  = err;
      		graph->GetEYhigh()[ibin] = err;
   }

   _file2->Close();
   return graph;
}

double compare_sig(TF1 *fit1, TF1 *fit2){ 

			// w/ nnodes = fit1 
			double chi2_1 = fit1->GetChisquare(); // chi-squared = sum of the square of the residuals, where residual = (data - fit)/data_error
			double dof_1 = fit1->GetNDF(); // number of degrees of freedom = number of data points - number of free parameters in the fit
			double norm_chi2_1 = chi2_1/dof_1; // normalized chi-squared

			// w/ nnodes = fit2 
			double chi2_2 = fit2->GetChisquare(); // chi-squared = sum of the square of the residuals, where residual = (data - fit)/data_error
			double dof_2 = fit2->GetNDF(); // number of degrees of freedom = number of data points - number of free parameters in the fit
			double norm_chi2_2 = chi2_2/dof_2; // normalized chi-squared

			double pvalue = TMath::Prob(chi2_1 - chi2_2, dof_1 - dof_2);
			double sigma = TMath::Sqrt2()*TMath::ErfcInverse(pvalue);	
				
			return sigma; 
}


double *get_kin_bins(const char *element)
{
   const int nBins = 14;
   static double kin_bins_B[nBins] = { 51.4, 65.8, 67.3, 90.4, 91.5, 110.6, 111.6, 128.4, 129.2, 144.6, 145.3, 159.4, 160.2, 173.7 };
   static double kin_bins_C[nBins] = { 59.0, 75.6, 77.2, 103.8, 105.1, 127.3, 128.3, 147.9, 148.8, 166.6, 167.4, 183.9, 184.8, 200.4 };
   static double kin_bins_N[nBins] = { 63.2, 81.0, 82.8, 111.4, 112.7, 136.6, 137.7, 158.8, 159.8, 179.0, 179.9, 197.6, 198.6, 215.5 };
   static double kin_bins_O[nBins] = { 69.4, 89.0, 91.0, 122.5, 124.0, 150.3, 151.6, 174.9, 176.0, 197.3, 198.3, 218.0, 219.1, 237.9 };
   static double kin_bins_F[nBins] = { 72.0, 92.4, 94.4, 127.2, 128.8, 156.2, 157.5, 181.8, 182.9, 205.1, 206.2, 226.8, 227.9, 247.5 };
   static double kin_bins_Ne[nBins] = { 77.1, 99.0, 101.2, 136.4, 138.1, 167.6, 169.1, 195.2, 196.5, 220.4, 221.6, 243.8, 245.0, 266.3 };
   static double kin_bins_Na[nBins] = { 81.0, 104.1, 106.4, 143.5, 145.3, 176.5, 178.0, 205.6, 206.9, 232.3, 233.5, 257.0, 258.4, 280.9 };
   static double kin_bins_Mg[nBins] = { 86.3, 110.9, 113.4, 153.1, 155.0, 188.4, 190.1, 219.7, 221.1, 248.4, 249.7, 275.0, 276.4, 300.7 };
   static double kin_bins_Al[nBins] = { 89.4, 115.0, 117.5, 158.8, 160.8, 195.5, 197.2, 228.1, 229.5, 257.9, 259.2, 285.6, 287.1, 312.4 };
   static double kin_bins_Si[nBins] = { 94.8, 121.9, 124.7, 168.6, 170.7, 207.8, 209.6, 242.6, 244.1, 274.5, 275.9, 304.2, 305.7, 332.9 };
   static double kin_bins_P[nBins] = { 97.1, 124.9, 127.7, 172.8, 175.0, 213.0, 214.9, 248.8, 250.4, 281.6, 283.1, 312.2, 313.8, 341.7 };
   static double kin_bins_S[nBins] = { 101.6, 130.8, 133.7, 181.2, 183.4, 223.5, 225.4, 261.2, 262.8, 295.8, 297.3, 328.1, 329.8, 359.3 };
   static double kin_bins_Cl[nBins] = { 103.2, 133.0, 136.0, 184.3, 186.6, 227.4, 229.4, 265.8, 267.5, 301.1, 302.7, 334.0, 335.8, 365.9 };
   static double kin_bins_Ar[nBins] = { 107.7, 138.8, 141.9, 192.4, 194.9, 237.7, 239.7, 278.0, 279.8, 315.1, 316.7, 349.7, 351.6, 383.3 };
   static double kin_bins_K[nBins] = { 110.0, 141.8, 145.0, 196.8, 199.3, 243.1, 245.3, 284.5, 286.3, 322.5, 324.2, 358.1, 360.0, 392.6 };
   static double kin_bins_Ca[nBins] = { 113.3, 146.2, 149.5, 203.1, 205.6, 251.0, 253.2, 293.9, 295.8, 333.3, 335.1, 370.2, 372.2, 406.0 };
   static double kin_bins_Sc[nBins] = { 114.7, 148.1, 151.4, 205.8, 208.4, 254.4, 256.7, 297.9, 299.9, 338.0, 339.8, 375.5, 377.5, 411.9 };
   static double kin_bins_Ti[nBins] = { 117.8, 152.1, 155.5, 211.5, 214.1, 261.6, 263.9, 306.5, 308.5, 347.9, 349.7, 386.6, 388.7, 424.2 };
   static double kin_bins_V[nBins] = { 120.2, 155.3, 158.8, 216.0, 218.8, 267.4, 269.8, 313.4, 315.4, 355.8, 357.7, 395.5, 397.6, 434.1 };
   static double kin_bins_Cr[nBins] = { 123.4, 159.6, 163.2, 222.2, 225.0, 275.2, 277.6, 322.6, 324.8, 366.5, 368.4, 407.6, 409.8, 447.5 };
   static double kin_bins_Mn[nBins] = { 126.0, 163.0, 166.7, 227.0, 230.0, 281.3, 283.9, 330.0, 332.2, 375.0, 377.0, 417.2, 419.4, 458.2 };
   static double kin_bins_Fe[nBins] = { 129.1, 167.0, 170.8, 232.9, 235.9, 288.7, 291.3, 338.9, 341.1, 385.2, 387.3, 428.7, 431.0, 471.0 };
   static double kin_bins_Co[nBins] = { 131.7, 170.6, 174.5, 238.0, 241.0, 295.2, 297.9, 346.6, 348.9, 394.1, 396.3, 438.8, 441.1, 482.2 };
   static double kin_bins_Ni[nBins] = { 136.2, 176.5, 180.5, 246.4, 249.6, 305.9, 308.7, 359.4, 361.8, 408.9, 411.2, 455.5, 458.0, 500.8 };
   // repeat for all elements
   
   	if (!strcmp(element, "B")) return kin_bins_B;
   else if (!strcmp(element, "C")) return kin_bins_C;
   else if (!strcmp(element, "N")) return kin_bins_N;
   else if (!strcmp(element, "O")) return kin_bins_O;
   else if (!strcmp(element, "F")) return kin_bins_F;
   else if (!strcmp(element, "Ne")) return kin_bins_Ne;
   else if (!strcmp(element, "Na")) return kin_bins_Na;
   else if (!strcmp(element, "Mg")) return kin_bins_Mg;
   else if (!strcmp(element, "Al")) return kin_bins_Al;
   else if (!strcmp(element, "Si")) return kin_bins_Si;
   else if (!strcmp(element, "P")) return kin_bins_P;
   else if (!strcmp(element, "S")) return kin_bins_S;
   else if (!strcmp(element, "Cl")) return kin_bins_Cl;
   else if (!strcmp(element, "Ar")) return kin_bins_Ar;
   else if (!strcmp(element, "K")) return kin_bins_K;
   else if (!strcmp(element, "Ca")) return kin_bins_Ca;
   else if (!strcmp(element, "Sc")) return kin_bins_Sc;
   else if (!strcmp(element, "Ti")) return kin_bins_Ti;
   else if (!strcmp(element, "Va")) return kin_bins_V;
   else if (!strcmp(element, "Cr")) return kin_bins_Cr;
   else if (!strcmp(element, "Mn")) return kin_bins_Mn;
   else if (!strcmp(element, "Fe")) return kin_bins_Fe;
   else if (!strcmp(element, "Co")) return kin_bins_Co;
   else if (!strcmp(element, "Ni")) return kin_bins_Ni;
}

double *get_spall_corr(const char *element)
{
	const int nBins = 7; 
	static double SpallCorr_B[nBins] = { 0.955, 0.929, 0.895, 0.863, 0.832, 0.802, 0.773 };
	static double SpallCorr_C[nBins] = { 0.953, 0.926, 0.891, 0.857, 0.825, 0.794, 0.764 };
	static double SpallCorr_N[nBins] = { 0.949, 0.920, 0.883, 0.847, 0.813, 0.781, 0.750 };
	static double SpallCorr_O[nBins] = { 0.947, 0.917, 0.878, 0.842, 0.807, 0.773, 0.741 };
	static double SpallCorr_F[nBins] = { 0.943, 0.911, 0.870, 0.832, 0.795, 0.760, 0.726 };
	static double SpallCorr_Ne[nBins] = { 0.940, 0.908, 0.866, 0.826, 0.788, 0.752, 0.718 };
	static double SpallCorr_Na[nBins] = { 0.938, 0.904, 0.860, 0.819, 0.780, 0.743, 0.708 };
	static double SpallCorr_Mg[nBins] = { 0.936, 0.901, 0.857, 0.815, 0.775, 0.737, 0.701 };
	static double SpallCorr_Al[nBins] = { 0.933, 0.897, 0.851, 0.808, 0.767, 0.728, 0.691 };
	static double SpallCorr_Si[nBins] = { 0.932, 0.895, 0.848, 0.804, 0.763, 0.724, 0.686 };
	static double SpallCorr_P[nBins] = { 0.929, 0.891, 0.843, 0.797, 0.754, 0.714, 0.676 };
	static double SpallCorr_S[nBins] = { 0.928, 0.889, 0.840, 0.794, 0.751, 0.710, 0.671 };
	static double SpallCorr_Cl[nBins] = { 0.925, 0.885, 0.834, 0.787, 0.742, 0.700, 0.660 };
	static double SpallCorr_Ar[nBins] = { 0.922, 0.881, 0.830, 0.782, 0.737, 0.694, 0.654 };
	static double SpallCorr_K[nBins] = { 0.920, 0.878, 0.826, 0.776, 0.730, 0.687, 0.646 };
	static double SpallCorr_Ca[nBins] = { 0.917, 0.874, 0.821, 0.771, 0.724, 0.681, 0.639 };
	static double SpallCorr_Sc[nBins] = { 0.915, 0.871, 0.817, 0.765, 0.718, 0.673, 0.631 };
	static double SpallCorr_Ti[nBins] = { 0.913, 0.869, 0.813, 0.761, 0.713, 0.667, 0.625 };
	static double SpallCorr_V[nBins] = { 0.911, 0.866, 0.809, 0.756, 0.707, 0.661, 0.618 };
	static double SpallCorr_Cr[nBins] = { 0.910, 0.864, 0.807, 0.753, 0.703, 0.657, 0.614 };
	static double SpallCorr_Mn[nBins] = { 0.907, 0.861, 0.802, 0.748, 0.698, 0.651, 0.607 };
	static double SpallCorr_Fe[nBins] = { 0.906, 0.859, 0.799, 0.744, 0.693, 0.646, 0.602 };
	static double SpallCorr_Co[nBins] = { 0.904, 0.856, 0.796, 0.740, 0.689, 0.641, 0.596 };
	static double SpallCorr_Ni[nBins] = { 0.904, 0.855, 0.795, 0.739, 0.687, 0.639, 0.595 };

	if (!strcmp(element, "B")) return SpallCorr_B;
   else if (!strcmp(element, "C")) return SpallCorr_C;
   else if (!strcmp(element, "N")) return SpallCorr_N;
   else if (!strcmp(element, "O")) return SpallCorr_O;
   else if (!strcmp(element, "F")) return SpallCorr_F;
   else if (!strcmp(element, "Ne")) return SpallCorr_Ne;
   else if (!strcmp(element, "Na")) return SpallCorr_Na;
   else if (!strcmp(element, "Mg")) return SpallCorr_Mg;
   else if (!strcmp(element, "Al")) return SpallCorr_Al;
   else if (!strcmp(element, "Si")) return SpallCorr_Si;
   else if (!strcmp(element, "P")) return SpallCorr_P;
   else if (!strcmp(element, "S")) return SpallCorr_S;
   else if (!strcmp(element, "Cl")) return SpallCorr_Cl;
   else if (!strcmp(element, "Ar")) return SpallCorr_Ar;
   else if (!strcmp(element, "K")) return SpallCorr_K;
   else if (!strcmp(element, "Ca")) return SpallCorr_Ca;
   else if (!strcmp(element, "Sc")) return SpallCorr_Sc;
   else if (!strcmp(element, "Ti")) return SpallCorr_Ti;
   else if (!strcmp(element, "Va")) return SpallCorr_V;
   else if (!strcmp(element, "Cr")) return SpallCorr_Cr;
   else if (!strcmp(element, "Mn")) return SpallCorr_Mn;
   else if (!strcmp(element, "Fe")) return SpallCorr_Fe;
   else if (!strcmp(element, "Co")) return SpallCorr_Co;
   else if (!strcmp(element, "Ni")) return SpallCorr_Ni;
}

double *get_spall_corr_unc(const char *element)
{
	const int nBins = 7; 
	static double SpallCorrUnc_B[nBins] = { 0.005, 0.007, 0.011, 0.015, 0.019, 0.022, 0.026 };
	static double SpallCorrUnc_C[nBins] = { 0.005, 0.008, 0.012, 0.016, 0.019, 0.023, 0.027 };
	static double SpallCorrUnc_N[nBins] = { 0.005, 0.008, 0.013, 0.017, 0.021, 0.025, 0.029 };
	static double SpallCorrUnc_O[nBins] = { 0.005, 0.009, 0.013, 0.017, 0.022, 0.026, 0.030 };
	static double SpallCorrUnc_F[nBins] = { 0.006, 0.009, 0.014, 0.019, 0.023, 0.028, 0.033 };
	static double SpallCorrUnc_Ne[nBins] = { 0.006, 0.010, 0.015, 0.019, 0.024, 0.029, 0.034 };
	static double SpallCorrUnc_Na[nBins] = { 0.006, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035 };
	static double SpallCorrUnc_Mg[nBins] = { 0.007, 0.010, 0.016, 0.021, 0.026, 0.031, 0.036 };
	static double SpallCorrUnc_Al[nBins] = { 0.007, 0.011, 0.016, 0.022, 0.027, 0.032, 0.038 };
	static double SpallCorrUnc_Si[nBins] = { 0.007, 0.011, 0.017, 0.022, 0.027, 0.033, 0.038 };
	static double SpallCorrUnc_P[nBins] = { 0.007, 0.012, 0.017, 0.023, 0.029, 0.034, 0.040 };
	static double SpallCorrUnc_S[nBins] = { 0.008, 0.012, 0.018, 0.023, 0.029, 0.035, 0.041 };
	static double SpallCorrUnc_Cl[nBins] = { 0.008, 0.012, 0.018, 0.024, 0.030, 0.036, 0.042 };
	static double SpallCorrUnc_Ar[nBins] = { 0.008, 0.013, 0.019, 0.025, 0.031, 0.037, 0.043 };
	static double SpallCorrUnc_K[nBins] = { 0.008, 0.013, 0.019, 0.026, 0.032, 0.038, 0.045 };
	static double SpallCorrUnc_Ca[nBins] = { 0.009, 0.013, 0.020, 0.026, 0.033, 0.039, 0.046 };
	static double SpallCorrUnc_Sc[nBins] = { 0.009, 0.014, 0.020, 0.027, 0.034, 0.040, 0.047 };
	static double SpallCorrUnc_Ti[nBins] = { 0.009, 0.014, 0.021, 0.028, 0.034, 0.041, 0.048 };
	static double SpallCorrUnc_V[nBins] = { 0.009, 0.015, 0.021, 0.028, 0.035, 0.042, 0.049 };
	static double SpallCorrUnc_Cr[nBins] = { 0.010, 0.015, 0.022, 0.029, 0.036, 0.043, 0.050 };
	static double SpallCorrUnc_Mn[nBins] = { 0.010, 0.015, 0.022, 0.029, 0.037, 0.044, 0.051 };
	static double SpallCorrUnc_Fe[nBins] = { 0.010, 0.015, 0.023, 0.030, 0.037, 0.045, 0.052 };
	static double SpallCorrUnc_Co[nBins] = { 0.010, 0.016, 0.023, 0.031, 0.038, 0.045, 0.053 };
	static double SpallCorrUnc_Ni[nBins] = { 0.010, 0.016, 0.023, 0.031, 0.038, 0.046, 0.053 };

	if (!strcmp(element, "B")) return SpallCorrUnc_B;
   else if (!strcmp(element, "C")) return SpallCorrUnc_C;
   else if (!strcmp(element, "N")) return SpallCorrUnc_N;
   else if (!strcmp(element, "O")) return SpallCorrUnc_O;
   else if (!strcmp(element, "F")) return SpallCorrUnc_F;
   else if (!strcmp(element, "Ne")) return SpallCorrUnc_Ne;
   else if (!strcmp(element, "Na")) return SpallCorrUnc_Na;
   else if (!strcmp(element, "Mg")) return SpallCorrUnc_Mg;
   else if (!strcmp(element, "Al")) return SpallCorrUnc_Al;
   else if (!strcmp(element, "Si")) return SpallCorrUnc_Si;
   else if (!strcmp(element, "P")) return SpallCorrUnc_P;
   else if (!strcmp(element, "S")) return SpallCorrUnc_S;
   else if (!strcmp(element, "Cl")) return SpallCorrUnc_Cl;
   else if (!strcmp(element, "Ar")) return SpallCorrUnc_Ar;
   else if (!strcmp(element, "K")) return SpallCorrUnc_K;
   else if (!strcmp(element, "Ca")) return SpallCorrUnc_Ca;
   else if (!strcmp(element, "Sc")) return SpallCorrUnc_Sc;
   else if (!strcmp(element, "Ti")) return SpallCorrUnc_Ti;
   else if (!strcmp(element, "Va")) return SpallCorrUnc_V;
   else if (!strcmp(element, "Cr")) return SpallCorrUnc_Cr;
   else if (!strcmp(element, "Mn")) return SpallCorrUnc_Mn;
   else if (!strcmp(element, "Fe")) return SpallCorrUnc_Fe;
   else if (!strcmp(element, "Co")) return SpallCorrUnc_Co;
   else if (!strcmp(element, "Ni")) return SpallCorrUnc_Ni;
}

double *get_EMed(const char *element)
{
	const int nBins = 7; 
	static double EMed_B[nBins] = {59.6, 79.7, 102.0, 121.1, 138.2, 154.0, 168.6}; 
	static double EMed_C[nBins] = {68.3, 91.5, 117.3, 139.3, 159.1, 177.4, 194.5};
	static double EMed_N[nBins] = {73.3, 98.1, 125.9, 149.6, 171.0, 190.7, 209.2};
	static double EMed_O[nBins] = {80.4, 107.8, 138.4, 164.7, 188.4, 210.3, 230.8};
	static double EMed_F[nBins] = {83.5, 112.0, 143.8, 171.1, 195.9, 218.7, 240.0};
	static double EMed_Ne[nBins] = {89.5, 120.1, 154.4, 183.9, 210.6, 235.3, 258.4};
	static double EMed_Na[nBins] = {94.0, 126.2, 162.4, 193.5, 221.7, 247.8, 272.3};
	static double EMed_Mg[nBins] = {100.2, 134.7, 173.4, 206.8, 237.1, 265.2, 291.5};
	static double EMed_Al[nBins] = {103.8, 139.6, 179.8, 214.5, 246.1, 275.3, 302.8};
	static double EMed_Si[nBins] = {110.1, 148.2, 191.1, 228.1, 261.8, 293.1, 322.6};
	static double EMed_P[nBins] = {112.7, 151.8, 195.9, 233.9, 268.6, 300.8, 331.1};
	static double EMed_S[nBins] = {118.2, 159.4, 205.8, 245.9, 282.5, 316.6, 348.7};
	static double EMed_Cl[nBins] = {120.2, 162.1, 209.4, 250.3, 287.7, 322.4, 355.1};
	static double EMed_Ar[nBins] = {125.0, 168.8, 218.1, 260.9, 300.0, 336.4, 370.8};
	static double EMed_K[nBins] = {127.9, 172.8, 223.4, 267.4, 307.5, 344.9, 380.3};
	static double EMed_Ca[nBins] = {131.6, 177.9, 230.1, 275.6, 317.1, 355.9, 392.4};
	static double EMed_Sc[nBins] = {133.5, 180.5, 233.7, 279.9, 322.2, 361.6, 398.8};
	static double EMed_Ti[nBins] = {137.1, 185.5, 240.3, 287.9, 331.6, 372.3, 410.8};
	static double EMed_V[nBins] = {139.9, 189.5, 245.5, 294.3, 339.1, 380.8, 420.3};
	static double EMed_Cr[nBins] = {144.0, 195.1, 253.0, 303.5, 349.8, 393.0, 434.0};
	static double EMed_Mn[nBins] = {146.8, 199.1, 258.3, 309.9, 357.3, 401.6, 443.5};
	static double EMed_Fe[nBins] = {150.4, 204.1, 265.0, 318.1, 366.9, 412.6, 455.9};
	static double EMed_Co[nBins] = {153.6, 208.5, 270.9, 325.3, 375.4, 422.3, 466.7};
	static double EMed_Ni[nBins] = {158.9, 215.9, 280.7, 337.3, 389.5, 438.4, 484.7}; 

	if (!strcmp(element, "B")) return EMed_B;
   else if (!strcmp(element, "C")) return EMed_C;
   else if (!strcmp(element, "N")) return EMed_N;
   else if (!strcmp(element, "O")) return EMed_O;
   else if (!strcmp(element, "F")) return EMed_F;
   else if (!strcmp(element, "Ne")) return EMed_Ne;
   else if (!strcmp(element, "Na")) return EMed_Na;
   else if (!strcmp(element, "Mg")) return EMed_Mg;
   else if (!strcmp(element, "Al")) return EMed_Al;
   else if (!strcmp(element, "Si")) return EMed_Si;
   else if (!strcmp(element, "P")) return EMed_P;
   else if (!strcmp(element, "S")) return EMed_S;
   else if (!strcmp(element, "Cl")) return EMed_Cl;
   else if (!strcmp(element, "Ar")) return EMed_Ar;
   else if (!strcmp(element, "K")) return EMed_K;
   else if (!strcmp(element, "Ca")) return EMed_Ca;
   else if (!strcmp(element, "Sc")) return EMed_Sc;
   else if (!strcmp(element, "Ti")) return EMed_Ti;
   else if (!strcmp(element, "Va")) return EMed_V;
   else if (!strcmp(element, "Cr")) return EMed_Cr;
   else if (!strcmp(element, "Mn")) return EMed_Mn;
   else if (!strcmp(element, "Fe")) return EMed_Fe;
   else if (!strcmp(element, "Co")) return EMed_Co;
   else if (!strcmp(element, "Ni")) return EMed_Ni;
}



