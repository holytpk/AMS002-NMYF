// Neutron Monitor Yield Function Analysis, replacing Koldobskiy's fluxes w/ our estimated fluxes 
// By Lingqiang He 06/21/2020

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
const int nBRs = Experiments::Info[Experiments::AMS02].Dataset[1].nMeasurements;  

const int F34 = 3; 

const int nNMs_useful = 15;
const char *NM_useful[nNMs_useful+1] = { "OULU", "PSNM", "MXCO", "HRMS", "JUNG", "JUNG1", "NEWK", "KIEL2", "APTY", "FSMT", "NAIN", "PWNK", "THUL", "SOPB", "SOPO" };

const int nNMs_Koldob = 7;
const char *NM_Koldob[nNMs_Koldob] = { "OULU", "HRMS", "MOSC", "APTY", "NEWK", "INVK", "ATHN" }; // add ATHN for analysis 
const double R_cutoff[nNMs_Koldob] = { 0.8, 4.58, 2.43, 0.65, 2.4, 0.3, 8.53 }; 
const int NM_tubes[nNMs_Koldob] = { 9, 12, 24, 18, 9, 18, 6 }; 

const double Rcut[nNMs_useful+1] = { 0.81, 16.80, 8.28, 4.58, 4.49, 4.49, 2.40, 2.36, 0.65, 0.30, 0.30, 0.30, 0.30, 0.10, 0.10 }; // Rigidity Cutoff in GV
const double alt[nNMs_useful+1] = { 15.0, 2565.0, 2274.0, 26.0, 3570.0, 3475.0, 50.0, 54.0, 181.0, 180.0, 46.0, 53.0, 26.0, 2820.0, 2820.0 }; // NM altitude

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

int choose_BO(const char *element){

	int temp = 0; 
	int pri = 3; //O, pri 
	int sec = 0; //B, sec or mix 

	for (int i=0; i<n_total; i++){
		if (!strcmp(element, Form("%s", Element[i]))) temp = i_temp[i]; 
	}

	if (temp==1) return pri;
	else if (temp==0 || temp==2) return sec; 
	  
} 

int isotope_size[24] = { 2, 3, 2, 3, 2, 3, 3, 3, 2, 5, 4, 7, 4, 8, 4, 9, 6, 7, 4, 5, 4, 7, 5, 8 }; // change when compare_isotope[][] is changed 

// list of functions 
UInt_t UTimeToBR(Long64_t utime);
UInt_t UBRToTime(int BR);
double *get_EMed(const char *element); // get energy bin mid point 
double *get_kin_bins(const char *element); // get kinetic energy bins 
double *get_spall_corr(const char *element); // get spallation correlation 
double *get_spall_corr_unc(const char *element); // get spallation correlation uncertainty
double compare_sig(TF1 *fit1, TF1 *fit2); // get sigma value between two fits
double get_Rcutoff(const char *NM); // get R_cutoff
int get_NM_tubes(const char *NM); // get # of NM tubes 

double get_k_BR_mean(TGraph *k_sf, const char *option); 
TGraphErrors *get_k_sf_ratio( TGraphErrors *k_sf1, TGraphErrors *k_sf2 ); 

TGraphAsymmErrors *get_ace_graph(const char *element, UInt_t iBin, UInt_t nBRs); // flux in Kinetic Energy over time 
TGraphAsymmErrors *get_ace_average_graph(const char *element, UInt_t *BRs, UInt_t nBRs); // flux in Kinetic over energy bins 

void nm_auto(); 
void nm_all_N_t(); 
void nm_FY(); 
void nm_plot_F34(const char *NM); 
void nm_plot_F34_auto(); 
void nm_plot_dis(const char *NM, const char *YF); 
TGraph *nm_reproduce1(const char *NM, const char *option1, const char *option2); // reproduce the results from Koldobisky, Mi13/Ma16. option1 = "k_norm", "k_sf" or "N_t", option2 = "Mi13" "Ma16"  
TGraph *nm_reproduce2(const char *NM, const char *option1, const char *option2); // reproduce the results from Koldobisky, but separate each contribution for p, He & elements above He 
										 // option 2 = "p" "He" "above" ""(all)
TH1 *heavier_he(); // return ratio heavier/helium

void nm_auto(){

	Debug::Enable(Debug::ALL);

	gSystem->mkdir("data/nm/reproduce2", true); 
 
	TLegend *legend1 = new TLegend(0.62,0.8,0.9,0.9); // left, down, right, top 
	legend1->SetNColumns(3); 
	
	TLegend *legend2 = new TLegend(0.32,0.8,0.68,0.9); // left, down, right, top 
	legend2->SetNColumns(3); 

	TCanvas *c = new TCanvas("c", "Estimated Scaling Factor", 2400, 1800); 
	c->Divide(1, 2); 

	TGraphErrors *k_sf1[nNMs_Koldob];
	TGraphErrors *k_sf2[nNMs_Koldob];

	gStyle->SetPalette(109); 

	// int F34 = 3;   

	gROOT->ProcessLine(Form(".> data/nm/reproduce2/F%d_k_BR_mean.txt", F34)); 

	for (int i=0; i<nNMs_Koldob; ++i){ 

		// Mi13
		k_sf1[i] = (TGraphErrors *) nm_reproduce1(NM_Koldob[i], "k_norm", "Mi13"); // Mi13

		legend1->AddEntry(k_sf1[i], Form("%s", NM_Koldob[i]), "l"); 

		c->cd(1);

		HistTools::SetStyle(k_sf1[i], HistTools::GetColorPalette(i, nNMs_Koldob), kFullCircle, 0.7, 1, 1);
		if (i==0) HistTools::SetStyle(k_sf1[i], kPink, kFullCircle, 0.7, 1, 1);

		k_sf1[i]->GetXaxis()->SetTimeDisplay(1);
  		k_sf1[i]->GetXaxis()->SetTimeFormat("%m-%y");
		k_sf1[i]->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
		k_sf1[i]->GetYaxis()->SetRangeUser(0.965, 1.035); 
		k_sf1[i]->SetTitle("Yield Function Mi13; ; Estimated NM Scaling Factor (Normalized)"); 
		 
		// PRINT_GRAPH(k_sf[i]); 

		if (i==0) k_sf1[i]->Draw("APL"); 

		// break; 
		if (i>0) k_sf1[i]->Draw("PLSAME"); 

		// Ma16
		k_sf2[i] = (TGraphErrors *) nm_reproduce1(NM_Koldob[i], "k_norm", "Ma16"); // Ma16

		legend2->AddEntry(k_sf2[i], Form("%s", NM_Koldob[i]), "l"); 

		//nm_reproduce1(NM_Koldob[i], "k_norm", "CM12"); 
		//nm_reproduce1(NM_Koldob[i], "k_norm", "CD00");  

		c->cd(2);

		HistTools::SetStyle(k_sf2[i], HistTools::GetColorPalette(i, nNMs_Koldob), kFullCircle, 0.7, 1, 1);
		if (i==0) HistTools::SetStyle(k_sf2[i], kPink, kFullCircle, 0.7, 1, 1);

		k_sf2[i]->GetXaxis()->SetTimeDisplay(1);
  		k_sf2[i]->GetXaxis()->SetTimeFormat("%m-%y");
		k_sf2[i]->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
		k_sf2[i]->GetYaxis()->SetRangeUser(0.965, 1.035); 
		k_sf2[i]->SetTitle("Yield Function Ma16; ; Estimated NM Scaling Factor (Normalized)"); 
		 
		// PRINT_GRAPH(k_sf2[i]); 

		if (i==0) k_sf2[i]->Draw("APL"); 

		// break; 
		if (i>0) k_sf2[i]->Draw("PLSAME"); 

	} 

	//k_sf1[nNMs_Koldob] = (TGraphErrors*) nm_reproduce1("ATHN","k_norm","Mi13");
	//k_sf2[nNMs_Koldob] = (TGraphErrors*) nm_reproduce1("ATHN","k_norm","Ma16"); 
	//nm_reproduce1("ATHN", "N_t", "CM12"); 
	//nm_reproduce1("ATHN", "N_t", "CD00");

	TFile *file_k = new TFile(Form("data/nm/reproduce2/F%d_k_sf.root", F34), "recreate");

	for (int i=0; i<nNMs_Koldob; i++){ 
		
		k_sf1[i]->Write(Form("k_sf1_%s", NM_Koldob[i])); 
		k_sf2[i]->Write(Form("k_sf2_%s", NM_Koldob[i])); 

	}
		
	file_k->Close(); 

	gROOT->ProcessLine(".> ");  

	// compute average NM count rates  
	// 1, my Mi13
	// 2, my Ma16

	TGraphErrors *k_ave1 = new TGraphErrors();
	TGraphErrors *k_ave2 = new TGraphErrors(); 

	for (int iBR=0; iBR<nBRs; iBR++){

		Double_t x1=0, y1=0, x2=0, y2=0; 

		Double_t ave1=0, sum1=0, ave2=0, sum2=0; 

		for (int i=0; i<nNMs_Koldob; i++){
		
			k_sf1[i]->GetPoint(iBR, x1, y1); 
			k_sf2[i]->GetPoint(iBR, x2, y2); 

			sum1 += y1; 
			sum2 += y2; 

			//printf("i = %d, x1 = %f, y1 = %f, sum1 = %f \n", iBR, x1, y1, sum1);  

		}
		
		ave1 = sum1/nNMs_Koldob;
		ave2 = sum2/nNMs_Koldob;

		Double_t std1=0, std2=0; 

		for (int i=0; i<nNMs_Koldob; i++){

			Double_t x1_i=0, y1_i=0, x2_i=0, y2_i=0; // for std  

			k_sf1[i]->GetPoint(iBR, x1_i, y1_i); 
			k_sf2[i]->GetPoint(iBR, x2_i, y2_i);

			std1 += pow(y1_i-ave1, 2);
			std2 += pow(y2_i-ave2, 2); 

		} 

		std1 = std1/(nNMs_Koldob-1); 
		std1 = sqrt(std1)/sqrt(nNMs_Koldob); 
	
		std2 = std2/(nNMs_Koldob-1); 
		std2 = sqrt(std2)/sqrt(nNMs_Koldob); 

		k_ave1->SetPoint(iBR, x1, ave1); 
		k_ave1->SetPointError(iBR, 0, std1); 
		k_ave2->SetPoint(iBR, x2, ave2); 
		k_ave2->SetPointError(iBR, 0, std2); 

	} 

	c->cd(1);

	HistTools::SetStyle(k_ave1, kBlack, kFullCircle, 0.9, 1, 3);  
	//PRINT_GRAPH(k_ave1);
	k_ave1->Draw("PLSAME"); 
	legend1->AddEntry(k_ave1, "Mean", "l");  
	legend1->Draw("SAME"); 
	c->cd(2);
	HistTools::SetStyle(k_ave2, kBlack, kFullCircle, 0.9, 1, 3);
	//PRINT_GRAPH(k_ave2);
	k_ave2->Draw("PLSAME"); 
	legend2->AddEntry(k_ave2, "Mean", "l");
	legend2->Draw("SAME"); 

	c->Print(Form("data/nm/reproduce2/F%d_estimated_nm_k_all.png", F34)); 

	// compare every single NM k_sf w/ Koldobskiy's result

	TGraphErrors *k_sf1_kol[nNMs_Koldob]; // Mi13
	TGraphErrors *k_sf2_kol[nNMs_Koldob]; // Ma16 

	for (int i=0; i<nNMs_Koldob-1; i++){

	   if (i!=1){

			k_sf1_kol[i] = new TGraphErrors(Form("./data/nm/reproduce/SF1-KOL-%s.dat", NM_Koldob[i]), "%lg %lg", ""); // Mi13
			k_sf2_kol[i] = new TGraphErrors(Form("./data/nm/reproduce/SF2-KOL-%s.dat", NM_Koldob[i]), "%lg %lg", ""); // Ma16

			// convert x points to unix time
			for (int iBR=0; iBR<nBRs; iBR++){

				Double_t x1=0, y1=0, // Ling's
					 x2=0, y2=0, // Koldobskiy's
					 x0=0, y0=0; // fix 
	
				k_sf1_kol[i]->GetPoint(iBR, x1, y1);
				k_sf2_kol[i]->GetPoint(iBR, x2, y2);

				k_sf1[i]->GetPoint(iBR, x0, y0); 

				x1 = x0;
				x2 = x0;  
				
				k_sf1_kol[i]->SetPoint(iBR, x1, y1);
				k_sf2_kol[i]->SetPoint(iBR, x2, y2);
				
			} 

			TGraphErrors *k_sf1_ratio = get_k_sf_ratio( k_sf1[i], k_sf1_kol[i] ); // Mi13	
			TGraphErrors *k_sf2_ratio = get_k_sf_ratio( k_sf2[i], k_sf2_kol[i] ); // Ma16
	
			TLegend *legend5 = new TLegend(0.62,0.8,0.9,0.9); // left, down, right, top 
			TLegend *legend6 = new TLegend(0.62,0.8,0.9,0.9); // left, down, right, top 	

			legend5->AddEntry(k_sf1[i], "Ling's", "l");
			legend5->AddEntry(k_sf1_kol[i], "Koldobskiy's", "l"); 
			legend6->AddEntry(k_sf2[i], "Ling's", "l");
			legend6->AddEntry(k_sf2_kol[i], "Koldobskiy's", "l");  
	
			TCanvas *c01 = new TCanvas("c01", "Comparison of Two Models", 2400, 1800); // Mi13
			c01->Divide(1, 2); 

			TCanvas *c02 = new TCanvas("c02", "Comparison of Two Models", 2400, 2800); // Ma16 
			c02->Divide(1, 2); 
	
			c01->cd(1);
			gPad->SetGrid(); 

			k_sf1[i]->Draw("A PL"); 
			HistTools::SetStyle(k_sf1_kol[i], kBlack, kFullCircle, 0.7, 1, 1); 
			k_sf1_kol[i]->Draw("PL SAME"); 
			//k_sf1_kol[i]->Print("range"); 
			legend5->Draw("SAME"); 

			c01->cd(2); 
			gPad->SetGrid();

			HistTools::SetStyle(k_sf1_ratio, kBlue, kFullCircle, 0.7, 1, 1);

			k_sf1_ratio->SetTitle(Form(" ; ; %s Scaling Factor Ratio (Ling's/Koldobskiy's) w/ Mi13", NM_Koldob[i]));
			k_sf1_ratio->GetXaxis()->SetTimeDisplay(1);
  			k_sf1_ratio->GetXaxis()->SetTimeFormat("%m-%y");
			k_sf1_ratio->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			k_sf1_ratio->GetYaxis()->SetRangeUser(0.8, 1.2);		

			k_sf1_ratio->Draw("A PL"); 

			c02->cd(1);
			gPad->SetGrid();
			
			k_sf2[i]->Draw("A PL"); 
			HistTools::SetStyle(k_sf2_kol[i], kBlack, kFullCircle, 0.7, 1, 1); 
			k_sf2_kol[i]->Draw("PL SAME"); 
			legend6->Draw("SAME");

			c02->cd(2);  
			gPad->SetGrid();

			HistTools::SetStyle(k_sf2_ratio, kBlue, kFullCircle, 0.7, 1, 1);

			k_sf2_ratio->SetTitle(Form(" ; ; %s Scaling Factor Ratio (Ling's/Koldobskiy's) w/ Ma16", NM_Koldob[i]));
			k_sf2_ratio->GetXaxis()->SetTimeDisplay(1);
  			k_sf2_ratio->GetXaxis()->SetTimeFormat("%m-%y");
			k_sf2_ratio->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			k_sf2_ratio->GetYaxis()->SetRangeUser(0.8, 1.2);			
			k_sf2_ratio->Draw("A PL");  

			c01->Print(Form("./data/nm/reproduce2/F%d_k_sf_ratio_%s_Mi13.png", F34, NM_Koldob[i]));
			c02->Print(Form("./data/nm/reproduce2/F%d_k_sf_ratio_%s_Ma16.png", F34, NM_Koldob[i]));
		} 

	   }

	// compare w/ Koldobskiy's NM k_sf Mean

	TGraphErrors *k_ave1_kol = new TGraphErrors("./data/nm/reproduce/k_ave1_kol.dat", "%lg %lg %lg", ""); // Mi13
	TGraphErrors *k_ave2_kol = new TGraphErrors("./data/nm/reproduce/k_ave2_kol.dat", "%lg %lg %lg", ""); // Ma16

	TGraphErrors *k_ave1_ratio = new TGraphErrors(); 
	TGraphErrors *k_ave2_ratio = new TGraphErrors();

	TGraphErrors *k_ave1_kol_err = new TGraphErrors(); 
	TGraphErrors *k_ave2_kol_err = new TGraphErrors(); 

	TLegend *legend3 = new TLegend(0.62,0.8,0.9,0.9); // left, down, right, top 
	TLegend *legend4 = new TLegend(0.62,0.8,0.9,0.9); // left, down, right, top 	

	for (int iBR=0; iBR<nBRs; iBR++){

		Double_t x1=0, y1=0, x2=0, y2=0; 
		Double_t xa=0, ya=0, xb=0, yb=0; // Koldob
		
		k_ave1->GetPoint(iBR, x1, y1); 
		k_ave2->GetPoint(iBR, x2, y2); 
		k_ave1_kol->GetPoint(iBR, xa, ya);
		k_ave2_kol->GetPoint(iBR, xb, yb); 

		k_ave1_ratio->SetPoint(iBR, x1, y1/ya); // Mi13
		k_ave2_ratio->SetPoint(iBR, x2, y2/yb); // Ma16 

		// error propogation example 
		Double_t dy1 = k_ave1->GetErrorY(iBR);  
		Double_t dya = k_ave1_kol->GetErrorY(iBR); 

		k_ave1_ratio->SetPoint(iBR, x1, y1/ya);			 
		k_ave1_ratio->SetPointError(iBR, 0, (y1/ya)*sqrt((dy1/y1)*(dy1/y1)+(dya/ya)*(dya/ya))); 

		Double_t dy2 = k_ave2->GetErrorY(iBR);  
		Double_t dyb = k_ave2_kol->GetErrorY(iBR); 

		k_ave2_ratio->SetPoint(iBR, x2, y2/yb);			 
		k_ave2_ratio->SetPointError(iBR, 0, (y2/yb)*sqrt((dy2/y2)*(dy2/y2)+(dyb/yb)*(dyb/yb))); 

		// error bars from the paper
		k_ave1_kol_err->SetPoint(iBR, x1, 1);
		k_ave1_kol_err->SetPointError(iBR, 0, dya);
			
		k_ave2_kol_err->SetPoint(iBR, x2, 1);
		k_ave2_kol_err->SetPointError(iBR, 0, dyb);

		// printf("iBR = %d, dy1 = %f, dya = %f \n", iBR, dy1, dya); 

	} 

	legend3->AddEntry(k_ave1_ratio, "Mi13", "l");
	// legend3->AddEntry(k_ave1_ratio, "Koldobskiy's STD", "");  
	legend4->AddEntry(k_ave2_ratio, "Ma16", "l"); 

	TCanvas *c1 = new TCanvas("c1", "Comparison of Two Models", 2400, 1800); 
	c1->Divide(1, 2); 
 
	HistTools::SetStyle(k_ave1_ratio, kRed, kFullCircle, 0.9, 1, 2);  
	HistTools::SetStyle(k_ave2_ratio, kBlue, kFullCircle, 0.9, 1, 2);

	HistTools::SetFillStyle(k_ave1_kol_err, kYellow-7, 1001);
	HistTools::SetFillStyle(k_ave2_kol_err, kYellow-7, 1001);
	
	k_ave1_kol_err->SetTitle(" ; ; NM Scaling Factor Ratio (Ling's/Koldobskiy's)"); 
	k_ave2_kol_err->SetTitle(" ; ; NM Scaling Factor Ratio (Ling's/Koldobskiy's)"); 
	
	c1->cd(1);
	gPad->SetGrid(); 

	k_ave1_kol_err->GetXaxis()->SetTimeDisplay(1);
  	k_ave1_kol_err->GetXaxis()->SetTimeFormat("%m-%y");
	k_ave1_kol_err->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	k_ave1_kol_err->GetYaxis()->SetRangeUser(0.985, 1.015); 

	k_ave1_kol_err->Draw("A E3"); // std from Kol 
 
	k_ave1_ratio->Draw("PL SAME");
	legend3->Draw("SAME"); 

	c1->cd(2); 
	gPad->SetGrid(); 

	k_ave2_kol_err->GetXaxis()->SetTimeDisplay(1);
  	k_ave2_kol_err->GetXaxis()->SetTimeFormat("%m-%y");
	k_ave2_kol_err->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	k_ave2_kol_err->GetYaxis()->SetRangeUser(0.985, 1.015);

	k_ave2_kol_err->Draw("A E3"); // std from Kol 

	k_ave2_ratio->Draw("PL SAME");
	legend4->Draw("SAME"); 
	
	c1->Print(Form("./data/nm/reproduce2/F%d_k_sf_ave_ratio.png", F34));

}

void compare_Rsum(){
	
	const char *NM = "OULU"; 
	const char *option2 = "Mi13";  

	TFile *file1 = new TFile(Form("data/nm/reproduce/KOL_%s_%s.root", NM, option2)); 
	TFile *file2 = new TFile(Form("data/nm/reproduce2/F%d_Analysis/F%d_%s_%s.root", F34, F34, NM, option2)); 

	TF1 *f_RP1 = (TF1*) file1->Get("f_RP"); 
	TF1 *f_RP2 = (TF1*) file2->Get("f_RP"); 

	HistTools::SetStyle(f_RP2, kBlue, kFullCircle, 1, 1, 2); 

	TCanvas *c0 = new TCanvas("c0","",1800,900);

	c0->cd(1);
	gPad->SetGrid(); 
	gPad->SetLogx();

	f_RP1->SetTitle("summed ratio as function of rigidity; rigidity (GV); ratio heavier/helium");
	//f_RP1->SetRange(0.1, 3e3); 
	f_RP1->GetXaxis()->SetRangeUser(0.1, 3e3);   
	f_RP1->GetYaxis()->SetRangeUser(0.18, 0.625);
	f_RP1->Draw();  	
	f_RP2->Draw("SAME");  

	TLegend *l0 = new TLegend(0.8,0.8,0.9,0.9); 
	l0->AddEntry(f_RP1, "Koldobskiy's", "l");
	l0->AddEntry(f_RP2, "F3/F4", "l"); 
	l0->Draw("SAME"); 

	c0->Print("data/nm/reproduce2/Rsum_F34_vs_Kol.png"); 

	return 0; 

}

// plot F3, F4 & Koldob N_t
void nm_plot_F34(const char *NM){

	Debug::Enable(Debug::ALL); 

	int f3 = 3; 
	int f4 = 4; 

	TCanvas *c0 = new TCanvas("c0", "", 1800, 900); 
	c0->Divide(2, 2); 

	TFile *file_f3 = new TFile(Form("data/nm/reproduce2/F%d_k_sf.root", f3)); 
	TFile *file_f4 = new TFile(Form("data/nm/reproduce2/F%d_k_sf.root", f4)); 
	
	TFile *file_reproduce1[2]; // i for F3 or F4
	for (int i=0; i<2; ++i) file_reproduce1[i] = new TFile(Form("data/nm/reproduce2/F%d_Analysis/F%d_%s_%s.root", f3, f3, NM, "Mi13")); 
	TFile *file_reproduce2[2];
	for (int i=0; i<2; ++i) file_reproduce2[i] = new TFile(Form("data/nm/reproduce2/F%d_Analysis/F%d_%s_%s.root", f4, f4, NM, "Ma16")); 

	TFile *file_kol1 = new TFile(Form("data/nm/reproduce/KOL_%s_%s.root", NM, "Mi13")); 
	TFile *file_kol2 = new TFile(Form("data/nm/reproduce/KOL_%s_%s.root", NM, "Ma16")); 

	TFile *file_nm = new TFile(Form("./data/nm/NM-%s.root", NM));
	TGraphErrors *N_nm = (TGraphErrors *) file_nm->Get("g_ave"); 

	TGraphErrors *N_t1[2], *N_t2[2]; 
	
	N_t1[0] = (TGraphErrors *) file_reproduce1[0]->Get("N_t"); 
	N_t2[0] = (TGraphErrors *) file_reproduce2[0]->Get("N_t"); 

	N_t1[1] = (TGraphErrors *) file_reproduce1[1]->Get("N_t"); 
	N_t2[1] = (TGraphErrors *) file_reproduce2[1]->Get("N_t");

	TGraphErrors *N_kol1 = (TGraphErrors *) file_kol1->Get("N_t");
	TGraphErrors *N_kol2 = (TGraphErrors *) file_kol2->Get("N_t");	

	TGraphErrors *k_sf1[2], *k_sf2[2]; 

	k_sf1[0] = (TGraphErrors*) file_f3->Get(Form("k_sf1_%s", NM)); 
	k_sf2[0] = (TGraphErrors*) file_f3->Get(Form("k_sf2_%s", NM)); 

	k_sf1[1] = (TGraphErrors*) file_f4->Get(Form("k_sf1_%s", NM)); 
	k_sf2[1] = (TGraphErrors*) file_f4->Get(Form("k_sf2_%s", NM)); 

	TGraphErrors *N_ratio1[3], *N_ratio2[3]; // F3, F4, Koldob vs. Observed 

	N_ratio1[0] = new TGraphErrors(); 
	N_ratio2[0] = new TGraphErrors(); 
	N_ratio1[1] = new TGraphErrors(); 
	N_ratio2[1] = new TGraphErrors();
	N_ratio1[2] = new TGraphErrors(); 
	N_ratio2[2] = new TGraphErrors();  

	N_ratio1[0]->GetXaxis()->SetTimeDisplay(1);
  	N_ratio1[0]->GetXaxis()->SetTimeFormat("%m-%y");
	N_ratio1[0]->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	N_ratio2[0]->GetXaxis()->SetTimeDisplay(1);
  	N_ratio2[0]->GetXaxis()->SetTimeFormat("%m-%y");
	N_ratio2[0]->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 

	for (int iBR=0; iBR<N_t1[0]->GetN(); ++iBR){

		N_ratio1[0]->SetPoint(iBR, N_t1[0]->GetX()[iBR], N_t1[0]->GetY()[iBR]/N_nm->GetY()[iBR]); 
		N_ratio2[0]->SetPoint(iBR, N_t2[0]->GetX()[iBR], N_t2[0]->GetY()[iBR]/N_nm->GetY()[iBR]); 
		N_ratio1[1]->SetPoint(iBR, N_t1[1]->GetX()[iBR], N_t1[1]->GetY()[iBR]/N_nm->GetY()[iBR]); 
		N_ratio2[1]->SetPoint(iBR, N_t2[1]->GetX()[iBR], N_t2[1]->GetY()[iBR]/N_nm->GetY()[iBR]);
		N_ratio1[2]->SetPoint(iBR, N_kol1->GetX()[iBR], N_kol1->GetY()[iBR]/N_nm->GetY()[iBR]); 
		N_ratio2[2]->SetPoint(iBR, N_kol2->GetX()[iBR], N_kol2->GetY()[iBR]/N_nm->GetY()[iBR]);
	}	

	c0->cd(1);
	gPad->SetGrid();  

	TLegend *l1 = new TLegend(0.8,0.8,0.9,0.9); 
	l1->AddEntry(N_t1[0], "F3", "pl"); 
	l1->AddEntry(N_t1[1], "F4", "pl"); 
	l1->AddEntry(N_kol1, "KOL", "pl"); 
	
	HistTools::SetStyle(N_t1[0], kBlue, kFullCircle, 0.8, 1, 1); 
	HistTools::SetStyle(N_t1[1], kPink, kFullCircle, 0.8, 1, 1); 
	HistTools::SetStyle(N_kol1, kBlack, kFullCircle, 0.8, 1, 1); 

	N_t1[0]->SetTitle(Form("Mi13; ;%s Estimated Count Rate", NM));  
	N_t1[0]->Draw("APL"); 
	N_t1[0]->GetYaxis()->SetRangeUser(0.9*N_t1[0]->GetMean(2), 1.15*N_t1[0]->GetMean(2)); 
	N_t1[1]->Draw("PL SAME"); 
	N_kol1->Draw("PL SAME");
	l1->Draw("SAME");   

	c0->cd(2);
	gPad->SetGrid();

	TLegend *l2 = new TLegend(0.8,0.8,0.9,0.9); 
	l2->AddEntry(N_t2[0], "F3", "pl"); 
	l2->AddEntry(N_t2[1], "F4", "pl"); 
	l2->AddEntry(N_kol2, "KOL", "pl"); 

	HistTools::SetStyle(N_t2[0], kBlue, kFullCircle, 0.8, 1, 1); 
	HistTools::SetStyle(N_t2[1], kPink, kFullCircle, 0.8, 1, 1); 
	HistTools::SetStyle(N_kol2, kBlack, kFullCircle, 0.8, 1, 1);
 
	N_t2[0]->SetTitle(Form("Ma16; ;%s Estimated Count Rate", NM)); 
	N_t2[0]->Draw("APL"); 
	N_t2[0]->GetYaxis()->SetRangeUser(0.9*N_t2[0]->GetMean(2), 1.15*N_t2[0]->GetMean(2)); 
	N_t2[1]->Draw("PL SAME"); 
	N_kol2->Draw("PL SAME"); 
	l2->Draw("SAME");
	
	c0->cd(3);
	gPad->SetGrid(); 

	TLegend *l3 = new TLegend(0.8,0.8,0.9,0.9); 
	l3->AddEntry(N_ratio1[0], "F3", "pl"); 
	l3->AddEntry(N_ratio1[1], "F4", "pl"); 
	l3->AddEntry(N_ratio1[2], "KOL", "pl"); 

	HistTools::SetStyle(N_ratio1[0], kBlue, kFullCircle, 0.8, 1, 1); 
	HistTools::SetStyle(N_ratio1[1], kPink, kFullCircle, 0.8, 1, 1); 
	HistTools::SetStyle(N_ratio1[2], kBlack, kFullCircle, 0.8, 1, 1);

	N_ratio1[0]->SetTitle("Mi13; ; Estimated Count Rate/Observed Count Rate"); 
	N_ratio1[0]->Draw("APL"); 
	N_ratio1[0]->GetYaxis()->SetRangeUser(0.9*N_ratio1[0]->GetMean(2), 1.15*N_ratio1[0]->GetMean(2)); 
	N_ratio1[1]->Draw("PL SAME");
	N_ratio1[2]->Draw("PL SAME"); 
	l3->Draw("SAME"); 

	c0->cd(4);
	gPad->SetGrid();

	TLegend *l4 = new TLegend(0.8,0.8,0.9,0.9); 
	l4->AddEntry(N_ratio2[0], "F3", "pl"); 
	l4->AddEntry(N_ratio2[1], "F4", "pl"); 
	l4->AddEntry(N_ratio2[2], "KOL", "pl"); 

	HistTools::SetStyle(N_ratio2[0], kBlue, kFullCircle, 0.8, 1, 1); 
	HistTools::SetStyle(N_ratio2[1], kPink, kFullCircle, 0.8, 1, 1); 
	HistTools::SetStyle(N_ratio2[2], kBlack, kFullCircle, 0.8, 1, 1);

	N_ratio2[0]->SetTitle("Ma16; ; Estimated Count Rate/Observed Count Rate");
	N_ratio2[0]->Draw("APL");
	N_ratio2[0]->GetYaxis()->SetRangeUser(0.9*N_ratio2[0]->GetMean(2), 1.15*N_ratio2[0]->GetMean(2)); 
	N_ratio2[1]->Draw("PL SAME"); 
	N_ratio2[2]->Draw("PL SAME"); 
	l4->Draw("SAME"); 

	c0->Print(Form("data/nm/reproduce2/plot_N_t_%s.png", NM)); 

}

void nm_plot_F34_auto(){ 

	for (int k=0; k<nNMs_Koldob; ++k){
	
		nm_plot_F34(Form("%s", NM_Koldob[k]));
		//nm_plot_dis(Form("%s", NM_Koldob[k]), "Mi13"); 

	} 

} 

// check discrepancy
// plot F3/F4 ratio at 1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000 GV. 
void nm_plot_dis(const char *NM, const char *YF){
 
	int nnodes = 7; 
	int nnodes_ams = 7; 

	TF1 *fyfp;
	TF1 *fyfhe; 

	double R_cha[11] = { 1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000 }; 

	if (!strcmp(YF, "Mi13")){
		fyfp = NeutronMonitors::YieldFunctions::CreateFunction("fyfp", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mishev13H1, Particle::PROTON, Energy::RIGIDITY);
		fyfhe = NeutronMonitors::YieldFunctions::CreateFunction("fyfhe", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mishev13He4, Particle::HELIUM4, Energy::RIGIDITY); 
	}else if (!strcmp(YF, "Ma16")){
		fyfp = NeutronMonitors::YieldFunctions::CreateFunction("fyfp", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mangeard16H1, Particle::PROTON, Energy::RIGIDITY);
		fyfhe = NeutronMonitors::YieldFunctions::CreateFunction("fyfhe", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mangeard16He4, Particle::HELIUM4, Energy::RIGIDITY); 
	}else if (!strcmp(YF, "CM12")){
		fyfp = NeutronMonitors::YieldFunctions::CreateFunction("fyfp", 0.1, 3e3, NeutronMonitors::YieldFunctions::CaballeroLopezMoraal12H1, Particle::PROTON, Energy::RIGIDITY);
		fyfhe = NeutronMonitors::YieldFunctions::CreateFunction("fyfhe", 0.1, 3e3, NeutronMonitors::YieldFunctions::CaballeroLopezMoraal12He4, Particle::HELIUM4, Energy::RIGIDITY); 
	}else if (!strcmp(YF, "CD00")){
		fyfp = NeutronMonitors::YieldFunctions::CreateFunction("fyfp", 0.1, 3e3, NeutronMonitors::YieldFunctions::ClemDorman00H1, Particle::PROTON, Energy::RIGIDITY);
		fyfhe = NeutronMonitors::YieldFunctions::CreateFunction("fyfhe", 0.1, 3e3, NeutronMonitors::YieldFunctions::ClemDorman00He4, Particle::HELIUM4, Energy::RIGIDITY); 
	}

	TCanvas *c3 = new TCanvas("c3", "", 1800, 900);
	c3->Divide(1, 1);  

	TCanvas *c4 = new TCanvas("c4", "", 1800, 900); 
	c4->Divide(1, 1); 

	for (int k=0; k<24; ++k){ 

		int iBR_true=0; 
		for (int i=0; i<nBRs; i++){

			int f34 = 3;

			TFile *file0 = new TFile(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[choose_BO(ACE_Element[k])], nnodes)); // BO temp 

			Spline *sp_comb = new Spline("sp_comb", nnodes, Spline::LogLog | Spline::PowerLaw);
			TF1 *fsp_comb = sp_comb->GetTF1Pointer();  
			TF1 *fit_comb = (TF1*) file0->Get("fit_both"); 
			HistTools::CopyParameters(fit_comb, fsp_comb); // error 
			double x1, x2;
			fit_comb->GetRange(x1,x2);
			fsp_comb->SetRange(x1,x2);

			// load ACE BR 
			TFile *file2 = new TFile(Form("data/ACE/fill/%s_fill.root", ACE_Element[k]));
			TH1 *h_ace_BR = (TH1*) file2->Get(Form("h_rig_%s_BR%d", ACE_Element[k], 2426+iBR_true)); 

			// rescaled combined fit 
			TH1 *h_ratio;
			h_ratio = (TH1*) h_ace_BR->Clone("h_ratio"); // ACE Ave/Temp 	
			h_ratio->Divide(fsp_comb); 

			// HistTools::SetStyle(h_ratio, kPink, kFullCircle, 1.4, 1, 1);

			double ratio_sum=0; // compute average of h_ratio manually  
			for(int nbin=0;nbin<14;++nbin){
				ratio_sum += h_ratio->GetBinContent(nbin);
				//printf("ratio_sum = %0.6f \n", ratio_sum);
			}
			double ratio_ave = ratio_sum/7;	

			// rescaled combined fit 
			TF1 *f_fit1 = HistTools::CombineTF1Const(fsp_comb, ratio_ave, HistTools::MultiplyConst, "rescaled_fit", 0.1, 3e3); 
			TF1 *f_fit2 = HistTools::CombineTF1Const(fsp_comb, ratio_ave, HistTools::MultiplyConst, "rescaled_fit", 0.1, 3e3); 
			//if (k!=22){ 

			// load F3 He(R,t)/<F(R,t)> 
			TFile *file_ratio1 = new TFile(Form("data/ACE/fill/F%d_%s.root", 3, ACE_Element[k])); // load ratio of F(R,t)/<F(R,t) fit
			Spline *sp_he = new Spline("sp_he", nnodes_ams, Spline::LogLog | Spline::PowerLaw);
			TF1 *fsp_he = sp_he->GetTF1Pointer();  
			TF1 *fit_he = (TF1*) file_ratio1->Get(Form("fit_he_BR%d", 2426+iBR_true)); 

			HistTools::CopyParameters(fit_he, fsp_he); 
			x1=0, x2=0; 
			fit_he->GetRange(x1,x2);
			fsp_he->SetRange(x1,x2);

			Spline *sp_he_ave = new Spline("sp_he_ave", nnodes_ams, Spline::LogLog | Spline::PowerLaw);
			TF1 *fsp_he_ave = sp_he_ave->GetTF1Pointer();  
			TF1 *fit_he_ave = (TF1*) file_ratio1->Get(Form("fit_he_ave_BR%d", 2426+iBR_true)); 

			HistTools::CopyParameters(fit_he_ave, fsp_he_ave);
			// double x1, x2;
			fit_he_ave->GetRange(x1,x2);
			fsp_he_ave->SetRange(x1,x2);

			TF1 *fit_ratio1 = HistTools::CombineTF1(fsp_he, fsp_he_ave, HistTools::Divide, "fsp_he2"); // He(R,t) Fit vs. <He(R)> Fit
			f_fit1 = HistTools::CombineTF1(f_fit1, fit_ratio1, HistTools::Multiply, "f_fit", 0.1, 3e3);

			// load F4 F(R,t)/<F(R,t)> ratio 
			TFile *file_ratio2 = new TFile(Form("data/ACE/fill/F%d_%s.root", 4, ACE_Element[k])); // load ratio of F(R,t)/<F(R,t) fit
			TF1 *fit_ratio2 = (TF1*) file_ratio2->Get(Form("fit_ratio_BR%d", 2426+iBR_true));	
			f_fit2 = HistTools::CombineTF1(f_fit2, fit_ratio2, HistTools::Multiply, "f_fit", 0.1, 3e3);

			TGraphErrors *g_f34_ratio = new TGraphErrors(); 

			for (int j=0; j<11; ++j){
				g_f34_ratio->SetPoint(j, R_cha[j], fit_ratio1->Eval(R_cha[j])/fit_ratio2->Eval(R_cha[j])); 
			} 
	
			c3->cd(1);
			gPad->SetGrid();
			//gPad->SetLogy();	
			gPad->SetLogx(); 

			g_f34_ratio->GetXaxis()->SetRangeUser(0.1, 3e3); 
			g_f34_ratio->GetYaxis()->SetNdivisions(20); 
			g_f34_ratio->SetTitle(Form("%s Flux_F3/Flux_F4 at Characteristic Rigidities; ; ratio of F3/F4", ACE_Element[k])); 
			g_f34_ratio->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV")); 

			HistTools::SetStyle(g_f34_ratio, kBlue, kFullCircle, 2, 1, 2); 
			g_f34_ratio->Draw("AP"); 

			if (i==0) c3->Print(Form("data/nm/reproduce2/F%d_discrepency/ratioF3F4_%s_%s.pdf(", f34, ACE_Element[k], YF), "pdf"); 	
			if (i>0 && i<nBRs-1) c3->Print(Form("data/nm/reproduce2/F%d_discrepency/ratioF3F4_%s_%s.pdf", f34, ACE_Element[k], YF), "pdf"); 	
			if (i==nBRs-1) c3->Print(Form("data/nm/reproduce2/F%d_discrepency/ratioF3F4_%s_%s.pdf)", f34, ACE_Element[k], YF), "pdf");  		

			file_ratio1->Close(); 
			file_ratio2->Close(); 
			// } 
/*
			// discrepancy vs. sensitivity check

			TF1 *f_yf_divide1 = HistTools::CombineTF1(f_fit1, fyfhe, HistTools::Multiply, "f_fit", 0.1, 3e3); 
			f_yf_divide1 = HistTools::CombineTF1Const(f_yf_divide1, 1/f_yf_divide1->Integral(0.1, 3e3), HistTools::MultiplyConst, "f_yf_divide", 0.1, 3e3);
			
			TF1 *f_yf_divide2 = HistTools::CombineTF1(f_fit2, fyfhe, HistTools::Multiply, "f_fit", 0.1, 3e3); 
			f_yf_divide2 = HistTools::CombineTF1Const(f_yf_divide2, 1/f_yf_divide2->Integral(0.1, 3e3), HistTools::MultiplyConst, "f_yf_divide", 0.1, 3e3);

			c4->cd(1); 

			gPad->SetGrid(); 
			gPad->SetLogy();
			gPad->SetLogx(); 
			gPad->SetMargin(0.12, 0.04, 0.1, 0.08);	
	
			TLegend *l1 = new TLegend(0.62,0.8,0.9,0.9); // left, down, right, top 
			l1->AddEntry(f_yf_divide1, "F3", "l");
			l1->AddEntry(f_yf_divide2, "F4", "l");

			f_yf_divide1->Draw(); 
			f_yf_divide2->Draw("SAME"); 
			l1->Draw("SAME"); 
			HistTools::SetStyle(f_yf_divide1, kRed, kFullCircle, 0.9, 1, 1); 
			HistTools::SetStyle(f_yf_divide2, kBlue, kFullCircle, 0.9, 1, 1); 
			f_yf_divide1->SetRange(0.8, 3e3);
			f_yf_divide1->GetYaxis()->SetRangeUser(1e-4, 1e-1); 
			f_yf_divide1->SetTitle(Form("F_%s*Y_He Helium Discrepancy vs. Sensitivity (F3 & F4);%s; Y*F/Integral BR%d", ACE_Element[k], Unit::GetEnergyLabel("GV"), iBR_true+2426));

			if (i==0) c4->Print(Form("data/nm/reproduce2/F%d_discrepency/YF_divide_%s_%s.pdf(", f34, ACE_Element[k], YF), "pdf"); 	
			if (i>0 && i<nBRs-1) c4->Print(Form("data/nm/reproduce2/F%d_discrepency/YF_divide_%s_%s.pdf", f34, ACE_Element[k], YF), "pdf"); 	
			if (i==nBRs-1) c4->Print(Form("data/nm/reproduce2/F%d_discrepency/YF_divide_%s_%s.pdf)", f34, ACE_Element[k], YF), "pdf");  

*/ 

	   		if (i+2426==2472-1) iBR_true += 3; 
			else iBR_true ++;

			file0->Close();
			file2->Close(); 
		} 
	}
}

void nm_plot_k(){

	Debug::Enable(Debug::ALL); 

	// int F34 = 3; 

	gSystem->mkdir(Form("data/nm/reproduce2/F%d_Analysis", F34), true);

	TFile *file_k = new TFile(Form("data/nm/reproduce2/F%d_k_sf.root", F34)); 

	TGraphErrors *N_t1[nNMs_Koldob]; 
	TGraphErrors *N_t2[nNMs_Koldob]; 

	TGraphErrors *k_sf1[nNMs_Koldob]; // Mi13
	TGraphErrors *k_sf2[nNMs_Koldob]; // Ma16

	TGraphErrors *k_sf1_kol[nNMs_Koldob]; // Mi13
	TGraphErrors *k_sf2_kol[nNMs_Koldob]; // Ma16 

	TGraphErrors *k_ave1 = new TGraphErrors();
	TGraphErrors *k_ave2 = new TGraphErrors(); 

	TFile *file_nm[nNMs_Koldob];
	TGraph *N_nm[nNMs_Koldob];  

	TCanvas *c0 = new TCanvas("c0", "Comparison of Three Models", 2400, 1800); 
	c0->Divide(1, 2); 

	for (int iBR=0; iBR<nBRs; iBR++){

		Double_t x1=0, y1=0, x2=0, y2=0; 

		Double_t ave1=0, sum1=0, ave2=0, sum2=0; 

		for (int i=0; i<nNMs_Koldob; ++i){
		
			k_sf1[i]->GetPoint(iBR, x1, y1); 
			k_sf2[i]->GetPoint(iBR, x2, y2);  

			sum1 += y1; 
			sum2 += y2; 

			//printf("i = %d, x1 = %f, y1 = %f, sum1 = %f \n", iBR, x1, y1, sum1);  

		}  
		
		ave1 = sum1/nNMs_Koldob;
		ave2 = sum2/nNMs_Koldob;

		Double_t std1=0, std2=0; 

		for (int i=0; i<nNMs_Koldob; i++){

			Double_t x1_i=0, y1_i=0, x2_i=0, y2_i=0; // for std  

			k_sf1[i]->GetPoint(iBR, x1_i, y1_i); 
			k_sf2[i]->GetPoint(iBR, x2_i, y2_i);

			std1 += pow(y1_i-ave1, 2);
			std2 += pow(y2_i-ave2, 2); 

		} 

		std1 = std1/(nNMs_Koldob-1); 
		std1 = sqrt(std1)/sqrt(nNMs_Koldob); 
	
		std2 = std2/(nNMs_Koldob-1); 
		std2 = sqrt(std2)/sqrt(nNMs_Koldob); 

		k_ave1->SetPoint(iBR, x1, ave1); 
		k_ave1->SetPointError(iBR, 0, std1); 
		k_ave2->SetPoint(iBR, x2, ave2); 
		k_ave2->SetPointError(iBR, 0, std2); 

	} 
		
	for (int i=0; i<nNMs_Koldob-1; i++){

	   // if (i!=1){

			k_sf1[i] = (TGraphErrors*) file_k->Get(Form("k_sf1_%s", NM_Koldob[i])); 
			k_sf2[i] = (TGraphErrors*) file_k->Get(Form("k_sf2_%s", NM_Koldob[i])); 

			k_sf1_kol[i] = new TGraphErrors(Form("./data/nm/reproduce/SF1-KOL-%s.dat", NM_Koldob[i]), "%lg %lg", ""); // Mi13
			k_sf2_kol[i] = new TGraphErrors(Form("./data/nm/reproduce/SF2-KOL-%s.dat", NM_Koldob[i]), "%lg %lg", ""); // Ma16 

			// convert x points to unix time
			for (int iBR=0; iBR<nBRs; iBR++){

				Double_t x1=0, y1=0, // Ling's
					 x2=0, y2=0, // Koldobskiy's
					 x0=0, y0=0; // fix 
	
				k_sf1_kol[i]->GetPoint(iBR, x1, y1);
				k_sf2_kol[i]->GetPoint(iBR, x2, y2);

				k_sf1[i]->GetPoint(iBR, x0, y0); 

				x1 = x0;
				x2 = x0;  
				
				k_sf1_kol[i]->SetPoint(iBR, x1, y1);
				k_sf2_kol[i]->SetPoint(iBR, x2, y2);
				
			} 

			TGraphErrors *k_sf1_ratio = get_k_sf_ratio( k_sf1[i], k_sf1_kol[i] ); // Mi13	
			TGraphErrors *k_sf2_ratio = get_k_sf_ratio( k_sf2[i], k_sf2_kol[i] ); // Ma16
	
			TLegend *legend5 = new TLegend(0.62,0.8,0.9,0.9); // left, down, right, top 
			TLegend *legend6 = new TLegend(0.62,0.8,0.9,0.9); // left, down, right, top 	

			legend5->AddEntry(k_sf1[i], "Ling's", "l");
			legend5->AddEntry(k_sf1_kol[i], "Koldobskiy's", "l"); 
			legend6->AddEntry(k_sf2[i], "Ling's", "l");
			legend6->AddEntry(k_sf2_kol[i], "Koldobskiy's", "l");  
	
			TCanvas *c01 = new TCanvas("c01", "Comparison of Two Models", 2400, 1800); // Mi13
			c01->Divide(1, 2); 

			TCanvas *c02 = new TCanvas("c02", "Comparison of Two Models", 2400, 2800); // Ma16 
			c02->Divide(1, 2); 
	
			c01->cd(1);
			gPad->SetGrid(); 

			k_sf1[i]->Draw("A PL"); 
			//k_sf1[i]->Print("range"); 
			HistTools::SetStyle(k_sf1_kol[i], kBlack, kFullCircle, 0.7, 1, 1); 
			k_sf1_kol[i]->Draw("PL SAME"); 
			//k_sf1_kol[i]->Print("range"); 
			legend5->Draw("SAME"); 

			c01->cd(2); 
			gPad->SetGrid();

			HistTools::SetStyle(k_sf1_ratio, kBlue, kFullCircle, 0.7, 1, 1);

			k_sf1_ratio->SetTitle(Form(" ; ; %s Scaling Factor Ratio (Ling's/Koldobskiy's) w/ Mi13", NM_Koldob[i]));
			k_sf1_ratio->GetXaxis()->SetTimeDisplay(1);
  			k_sf1_ratio->GetXaxis()->SetTimeFormat("%m-%y");
			k_sf1_ratio->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			k_sf1_ratio->GetYaxis()->SetRangeUser(0.5, 1.5);		

			k_sf1_ratio->Draw("A PL"); 

			c02->cd(1);
			gPad->SetGrid();
			
			k_sf2[i]->Draw("A PL"); 
			HistTools::SetStyle(k_sf2_kol[i], kBlack, kFullCircle, 0.7, 1, 1); 
			k_sf2_kol[i]->Draw("PL SAME"); 
			legend6->Draw("SAME");

			c02->cd(2);  
			gPad->SetGrid();

			HistTools::SetStyle(k_sf2_ratio, kBlue, kFullCircle, 0.7, 1, 1);

			k_sf2_ratio->SetTitle(Form(" ; ; %s Scaling Factor Ratio (Ling's/Koldobskiy's) w/ Ma16", NM_Koldob[i]));
			k_sf2_ratio->GetXaxis()->SetTimeDisplay(1);
  			k_sf2_ratio->GetXaxis()->SetTimeFormat("%m-%y");
			k_sf2_ratio->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			k_sf2_ratio->GetYaxis()->SetRangeUser(0.5, 1.5);			
			k_sf2_ratio->Draw("A PL");  

			//c01->Print(Form("data/nm/reproduce2/F%d_Analysis/F%d_k_sf_ratio_%s_Mi13.png", F34, F34, NM_Koldob[i]));
			//c02->Print(Form("data/nm/reproduce2/F%d_Analysis/F%d_k_sf_ratio_%s_Ma16.png", F34, F34, NM_Koldob[i]));
		// } 

	   }

	// compare w/ Koldobskiy's NM k_sf Mean

	TGraphErrors *k_ave1_kol = new TGraphErrors("./data/nm/reproduce2/k_ave1_kol.dat", "%lg %lg %lg", ""); // Mi13
	TGraphErrors *k_ave2_kol = new TGraphErrors("./data/nm/reproduce2/k_ave2_kol.dat", "%lg %lg %lg", ""); // Ma16

	TGraphErrors *k_ave1_ratio = new TGraphErrors(); 
	TGraphErrors *k_ave2_ratio = new TGraphErrors();

	TGraphErrors *k_ave1_kol_err = new TGraphErrors(); 
	TGraphErrors *k_ave2_kol_err = new TGraphErrors(); 

	TLegend *legend3 = new TLegend(0.62,0.8,0.9,0.9); // left, down, right, top 
	TLegend *legend4 = new TLegend(0.62,0.8,0.9,0.9); // left, down, right, top 	

	for (int iBR=0; iBR<nBRs; iBR++){

		Double_t x1=0, y1=0, x2=0, y2=0; 
		Double_t xa=0, ya=0, xb=0, yb=0; // Koldob
		
		k_ave1->GetPoint(iBR, x1, y1); 
		k_ave2->GetPoint(iBR, x2, y2); 
		k_ave1_kol->GetPoint(iBR, xa, ya);
		k_ave2_kol->GetPoint(iBR, xb, yb); 

		k_ave1_ratio->SetPoint(iBR, x1, y1/ya); // Mi13
		k_ave2_ratio->SetPoint(iBR, x2, y2/yb); // Ma16 

		// error propogation example 
		Double_t dy1 = k_ave1->GetErrorY(iBR);  
		Double_t dya = k_ave1_kol->GetErrorY(iBR); 

		k_ave1_ratio->SetPoint(iBR, x1, y1/ya);			 
		k_ave1_ratio->SetPointError(iBR, 0, (y1/ya)*sqrt((dy1/y1)*(dy1/y1)+(dya/ya)*(dya/ya))); 

		Double_t dy2 = k_ave2->GetErrorY(iBR);  
		Double_t dyb = k_ave2_kol->GetErrorY(iBR); 

		k_ave2_ratio->SetPoint(iBR, x2, y2/yb);			 
		k_ave2_ratio->SetPointError(iBR, 0, (y2/yb)*sqrt((dy2/y2)*(dy2/y2)+(dyb/yb)*(dyb/yb))); 

		// error bars from the paper
		k_ave1_kol_err->SetPoint(iBR, x1, 1);
		k_ave1_kol_err->SetPointError(iBR, 0, dya);
			
		k_ave2_kol_err->SetPoint(iBR, x2, 1);
		k_ave2_kol_err->SetPointError(iBR, 0, dyb);

		// printf("iBR = %d, dy1 = %f, dya = %f \n", iBR, dy1, dya); 

	} 

	legend3->AddEntry(k_ave1_ratio, "Mi13", "l");
	// legend3->AddEntry(k_ave1_ratio, "Koldobskiy's STD", "");  
	legend4->AddEntry(k_ave2_ratio, "Ma16", "l"); 

	TCanvas *c1 = new TCanvas("c1", "Comparison of Two Models", 2400, 1800); 
	c1->Divide(1, 2); 
 
	HistTools::SetStyle(k_ave1_ratio, kRed, kFullCircle, 0.9, 1, 2);  
	HistTools::SetStyle(k_ave2_ratio, kBlue, kFullCircle, 0.9, 1, 2);

	HistTools::SetFillStyle(k_ave1_kol_err, kYellow-7, 1001);
	HistTools::SetFillStyle(k_ave2_kol_err, kYellow-7, 1001);
	
	k_ave1_kol_err->SetTitle(" ; ; NM Scaling Factor Ratio (Ling's/Koldobskiy's)"); 
	k_ave2_kol_err->SetTitle(" ; ; NM Scaling Factor Ratio (Ling's/Koldobskiy's)"); 
	
	c1->cd(1);
	gPad->SetGrid(); 

	k_ave1_kol_err->GetXaxis()->SetTimeDisplay(1);
  	k_ave1_kol_err->GetXaxis()->SetTimeFormat("%m-%y");
	k_ave1_kol_err->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	k_ave1_kol_err->GetYaxis()->SetRangeUser(0.985, 1.015); 

	k_ave1_kol_err->Draw("A E3"); // std from Kol 
 
	k_ave1_ratio->Draw("PL SAME");
	legend3->Draw("SAME"); 

	c1->cd(2); 
	gPad->SetGrid(); 

	k_ave2_kol_err->GetXaxis()->SetTimeDisplay(1);
  	k_ave2_kol_err->GetXaxis()->SetTimeFormat("%m-%y");
	k_ave2_kol_err->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	k_ave2_kol_err->GetYaxis()->SetRangeUser(0.985, 1.015);

	k_ave2_kol_err->Draw("A E3"); // std from Kol 

	k_ave2_ratio->Draw("PL SAME");
	legend4->Draw("SAME"); 
	
	//c1->Print(Form("data/nm/reproduce2/F%d_Analysis/F%d_k_sf_ave_ratio.png", F34, F34));
}

// plot all N(t) in one canvas
// plot N(t) by p, He, element above He contribution 
void nm_all_N_t(){

	Debug::Enable(Debug::ALL);

	gSystem->mkdir("data/nm/reproduce2", true); 

	TGraph *N_nm[nNMs_Koldob];
	TGraph *N_nm_ave[nNMs_Koldob]; 
	TGraphErrors *N_simu1[nNMs_Koldob]; //Mi13
	TGraphErrors *N_simu2[nNMs_Koldob]; //Ma16 

	TGraph *N_simu_p[nNMs_Koldob]; 
	TGraph *N_simu_He[nNMs_Koldob]; 
	TGraph *N_simu_above[nNMs_Koldob]; 
	
	for (int i=0; i<nNMs_Koldob; i++){

		TFile *nm_data = new TFile(Form("./data/nm/NM-%s.root", NM_Koldob[i]));
		N_nm[i] = (TGraph*) nm_data->Get("g");
		N_nm_ave[i] = (TGraph*) nm_data->Get("g_ave");

		N_simu1[i] = (TGraphErrors*) nm_reproduce1(NM_Koldob[i], "N_t", "Mi13"); 
		N_simu2[i] = (TGraphErrors*) nm_reproduce1(NM_Koldob[i], "N_t", "Ma16"); 

		N_simu_p[i] = (TGraph*) nm_reproduce2(NM_Koldob[i], "N_t", "p"); 
		N_simu_p[i]->SetTitle(Form(" ; ; %s NM Simulated Partial Count Rate (p)", NM_Koldob[i])); 
		N_simu_He[i] = (TGraph*) nm_reproduce2(NM_Koldob[i], "N_t", "He");  
		N_simu_He[i]->SetTitle(Form(" ; ; %s NM Simulated Partial Count Rate (He)", NM_Koldob[i])); 
		N_simu_above[i] = (TGraph*) nm_reproduce2(NM_Koldob[i], "N_t", "above"); 
		N_simu_above[i]->SetTitle(Form(" ; ; %s NM Simulated Partial Count Rate (above He)", NM_Koldob[i]));  

		TLegend *legend = new TLegend(0.8,0.8,0.9,0.9); // left, down, right, top 

		TCanvas *c = new TCanvas("c","", 2400, 1600); 

		c->cd(1);
		N_nm[i]->SetTitle(Form("%s Count Rate Data & Simuations;Time; Count Rate", NM_Koldob[i]));
		N_nm[i]->GetXaxis()->SetTimeDisplay(1);
 		N_nm[i]->GetXaxis()->SetTimeFormat("%m-%y");
 		N_nm[i]->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00");
		N_nm[i]->GetYaxis()->SetRangeUser(0., N_nm[i]->Eval(UBRToTime(2426))*1.8); 
		N_nm[i]->Draw("ALP"); 
		N_nm_ave[i]->SetMarkerColor(kPink);
		N_nm_ave[i]->SetMarkerStyle(kFullCircle);
		N_nm_ave[i]->Draw("PSAME"); 
		N_simu1[i]->SetMarkerColor(kBlue-4);
		N_simu1[i]->SetMarkerStyle(kFullCircle);
		N_simu1[i]->Draw("PSAME");
		N_simu2[i]->SetMarkerColor(kBlue+4);
		N_simu2[i]->SetMarkerStyle(kFullCircle);
		N_simu2[i]->Draw("PSAME"); 

		legend->AddEntry(N_nm[i], "Raw", "lp"); 
		legend->AddEntry(N_nm_ave[i], "Averaged", "p"); 
		legend->AddEntry(N_simu1[i], "Mi13", "p"); 
		legend->AddEntry(N_simu2[i], "Ma16", "p"); 

		legend->Draw("SAME");  

		c->Print(Form("./data/nm/reproduce2/N_all_%s.png", NM_Koldob[i])); 

		TFile *NM_simu = new TFile(Form("./data/nm/NM_simu_%s.root", NM_Koldob[i]), "recreate");
		
		N_simu_p[i]->Write("N_simu_p");
		N_simu_He[i]->Write("N_simu_He");
		N_simu_above[i]->Write("N_simu_above");
		
		NM_simu->Write();
		NM_simu->Close(); 
	
	} 

}

// compute F(R,t)*Y(R)
void nm_FY(){

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

		// HistTools::PrintFunction(f_fit[i]); 

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

// Reconstruct the NM count with a simple cosmic ray contribution model from Koldobisky's assumption, using Mi13
// option1 = "k_norm" or "k_sf" or "N_t", option2 = "Mi13" or "Ma16" 
// using F3 fits 
TGraph *nm_reproduce1(const char *NM, const char *option1, const char *option2){
	
	Debug::Enable(Debug::ALL); 

	gSystem->mkdir("data/nm/reproduce2", true);

	gSystem->mkdir(Form("data/nm/reproduce2/F%d_Analysis", F34), true);

	int nnodes_ams = 6; 	
	int nnodes_ams2 = 7; // fsp_he in F3 
	double n_tubes = (double) get_NM_tubes(NM)/6; 
	double R_cutoff = get_Rcutoff(NM); 

	// printf("n_tubes (%s, %s) = %2.2f \n", NM, option2, n_tubes); 

	TGraph *N_t = new TGraph();

	Experiments::DataPath = "data"; 

	TF1 *f_fit[n_total]; 

	// TH1 *J_sum = heavier_he(); 
	// J_sum->Print("range"); 

	TF1 *fyfp;
	TF1 *fyfhe; 

	if (!strcmp(option2, "Mi13")){
		fyfp = NeutronMonitors::YieldFunctions::CreateFunction("fyfp", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mishev13H1, Particle::PROTON, Energy::RIGIDITY);
		fyfhe = NeutronMonitors::YieldFunctions::CreateFunction("fyfhe", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mishev13He4, Particle::HELIUM4, Energy::RIGIDITY); 
	}else if (!strcmp(option2, "Ma16")){
		fyfp = NeutronMonitors::YieldFunctions::CreateFunction("fyfp", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mangeard16H1, Particle::PROTON, Energy::RIGIDITY);
		fyfhe = NeutronMonitors::YieldFunctions::CreateFunction("fyfhe", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mangeard16He4, Particle::HELIUM4, Energy::RIGIDITY); 
	}else if (!strcmp(option2, "CM12")){
		fyfp = NeutronMonitors::YieldFunctions::CreateFunction("fyfp", 0.1, 3e3, NeutronMonitors::YieldFunctions::CaballeroLopezMoraal12H1, Particle::PROTON, Energy::RIGIDITY);
		fyfhe = NeutronMonitors::YieldFunctions::CreateFunction("fyfhe", 0.1, 3e3, NeutronMonitors::YieldFunctions::CaballeroLopezMoraal12He4, Particle::HELIUM4, Energy::RIGIDITY); 
	}else if (!strcmp(option2, "CD00")){
		fyfp = NeutronMonitors::YieldFunctions::CreateFunction("fyfp", 0.1, 3e3, NeutronMonitors::YieldFunctions::ClemDorman00H1, Particle::PROTON, Energy::RIGIDITY);
		fyfhe = NeutronMonitors::YieldFunctions::CreateFunction("fyfhe", 0.1, 3e3, NeutronMonitors::YieldFunctions::ClemDorman00He4, Particle::HELIUM4, Energy::RIGIDITY); 
	}

	//HistTools::PrintFunction(fyfp);
	//HistTools::PrintFunction(fyfhe); 

	TF1 *f_BR_p[nBRs];   
	TF1 *f_BR_he[nBRs];  

	TCanvas *c3 = new TCanvas("c3", "R_sum", 1800, 900); 

	int nnodes = 7; 

	TFile *file0 = new TFile(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[3], nnodes)); // O temp 

	Spline *sp_O = new Spline("sp_O", nnodes, Spline::LogLog | Spline::PowerLaw);
	TF1 *fsp_O = sp_O->GetTF1Pointer();  
	TF1 *fit_O = (TF1*) file0->Get("fit_both"); 

	HistTools::CopyParameters(fit_O, fsp_O); // error 
	double x1, x2;
	fit_O->GetRange(x1,x2);
	fsp_O->SetRange(x1,x2);

	TFile *fit_result = new TFile(Form("data/amsfit/fit_result_node%d.root", nnodes_ams));
	TFile *nm_data = new TFile(Form("./data/nm/NM-%s.root", NM));
	TGraph *N_nm = (TGraph*) nm_data->Get("g_ave"); 

	TGraph *k_sf = new TGraph(); // k = N_t/N_nm, scaling factor of NM stations 

	// int F34 = 3;  

	TF1 *f_RP; 

	gSystem->mkdir(Form("data/nm/reproduce2/F%d_discrepency", F34), true);

	int iBR_true=0; 
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

		FunctorExample *fe = new FunctorExample("fe", 0.1, 3e3, fyfp, fyfhe); 

		TF1 *f0 = new TF1("f0", "10e-25", 0.1, 3e3); // empty 
	
		for(int k=0; k<n_ams; k++){

			Spline *sp_ams = new Spline("sp_ams", nnodes_ams, Spline::LogLog | Spline::PowerLaw); 
			f_fit[k] = sp_ams->GetTF1Pointer();  // real function 
			TF1 *fit_ams = (TF1*) fit_result->Get(Form("fsp_%s", AMS_Element[k]))->Clone(Form("f_%s", AMS_Element[k])); 
	
			HistTools::CopyParameters(fit_ams, f_fit[k]); 
			double x1, x2;
			fit_ams->GetRange(x1,x2); 
			f_fit[k]->SetRange(x1,x2); 

			//fe->AddElementFlux(f0, A[k]); 

			fe->AddElementFlux_Ave(f_fit[k], A[k]);  
			// if (i==nBRs-1) HistTools::PrintFunction(fit_ams);   
		}  

		TF1 *f_BR_p_low = HistTools::CombineTF1Const(fsp_O, f_BR_p[i]->Eval(2.)/fsp_O->Eval(2.), HistTools::MultiplyConst, "rescaled_fit", 0.1, 3e3); 
		TF1 *f_BR_he_low = HistTools::CombineTF1Const(fsp_O, f_BR_he[i]->Eval(2.)/fsp_O->Eval(2.), HistTools::MultiplyConst, "rescaled_fit", 0.1, 3e3);

/* TEST THE EXTRAPOLATION 
		f_BR_he[i]->SetRange(0.1, 3e3); 
		f_BR_he[i]->Draw(); 
		f_BR_he_low->SetLineColor(kGreen); 
		f_BR_he_low->Draw("SAME"); 

		TLegend *l_he = new TLegend(0.9,0.7,0.7,0.9); 
		l_he->AddEntry(f_BR_he[i], "He Flux", "l"); 
		l_he->AddEntry(f_BR_he_low, "Rescaled O Temp", "l"); 
		l_he->Draw("SAME"); 

		f_BR_p[i]->SetRange(0.1, 3e3); 
		f_BR_p[i]->Draw(); 
		f_BR_p_low->SetLineColor(kGreen); 
		f_BR_p_low->Draw("SAME"); 

		TLegend *l_p = new TLegend(0.9,0.7,0.7,0.9); 
		l_p->AddEntry(f_BR_he[i], "Proton Flux", "l"); 
		l_p->AddEntry(f_BR_he_low, "Rescaled O Temp", "l"); 
		l_p->Draw("SAME"); 

		gPad->SetLogy();
		gPad->SetLogx(); 

		//double xmin, xmax; 
		//f_BR_he[i]->GetRange(xmin, xmax); 
		//printf(" (%0.4f, %0.4f) \n", xmin, xmax); 
*/

		fe->SetProtonFlux(f_BR_p[i]); 
		if (F34==3) fe->SetProtonFluxLow(f_BR_p[i]); 
		else if (F34==4) fe->SetProtonFluxLow(f_BR_p_low); 
		fe->AddElementFlux(f_BR_he[i], A[1]); 
		if (F34==3) fe->AddElementFluxLow(f_BR_he[i]); 
		else if (F34==4) fe->AddElementFluxLow(f_BR_he_low); 

		for(int k=4; k<n_ele; ++k){

		 //  if (k==22){

			TFile *file0 = new TFile(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[choose_BO(ACE_Element[k])], nnodes)); // BO temp 

			Spline *sp_comb = new Spline("sp_comb", nnodes, Spline::LogLog | Spline::PowerLaw);
			TF1 *fsp_comb = sp_comb->GetTF1Pointer();  
			TF1 *fit_comb = (TF1*) file0->Get("fit_both"); 

			HistTools::CopyParameters(fit_comb, fsp_comb); // error 
			double x1, x2;
			fit_comb->GetRange(x1,x2);
			fsp_comb->SetRange(x1,x2);
	
			// load ACE BR 
			TFile *file2 = new TFile(Form("data/ACE/fill/%s_fill.root", ACE_Element[k]));
			TH1 *h_ace_BR = (TH1*) file2->Get(Form("h_rig_%s_BR%d", ACE_Element[k], 2426+iBR_true)); 

			// rescaled combined fit 
			TH1 *h_ratio;
			h_ratio = (TH1*) h_ace_BR->Clone("h_ratio"); // ACE Ave/Temp 	
			h_ratio->Divide(fsp_comb); 

			// HistTools::SetStyle(h_ratio, kPink, kFullCircle, 1.4, 1, 1);

			double ratio_sum=0; // compute average of h_ratio manually  
			for(int nbin=0;nbin<14;++nbin){
				ratio_sum += h_ratio->GetBinContent(nbin); 
				//printf("ratio_sum = %0.6f \n", ratio_sum);
			}
			double ratio_ave = ratio_sum/7;		

			// rescaled combined fit 
			f_fit[k] = HistTools::CombineTF1Const(fsp_comb, ratio_ave, HistTools::MultiplyConst, "rescaled_fit", 0.1, 3e3); 

			fe->AddElementFlux_Ave(f_fit[k], A[k]); 

			if (F34==3){

				// load F3 He(R,t)/<F(R,t)> 
				TFile *file_ratio1 = new TFile(Form("data/ACE/fill/F%d_%s.root", F34, ACE_Element[k])); // load ratio of F(R,t)/<F(R,t) fit
				Spline *sp_he1 = new Spline("sp_he1", nnodes_ams, Spline::LogLog | Spline::PowerLaw);
				TF1 *fsp_he1 = sp_he1->GetTF1Pointer();  
				TF1 *fit_he1 = (TF1*) file_ratio1->Get(Form("fit_he_BR%d", 2426+iBR_true)); 

				HistTools::CopyParameters(fit_he1, fsp_he1); 
				x1=0, x2=0; 
				fit_he1->GetRange(x1,x2);
				fsp_he1->SetRange(x1,x2);

				Spline *sp_he_ave = new Spline("sp_he_ave", nnodes_ams, Spline::LogLog | Spline::PowerLaw);
				TF1 *fsp_he_ave = sp_he_ave->GetTF1Pointer();  
				TF1 *fit_he_ave = (TF1*) file_ratio1->Get(Form("fit_he_ave_BR%d", 2426+iBR_true)); 

				HistTools::CopyParameters(fit_he_ave, fsp_he_ave);
				// double x1, x2;
				fit_he_ave->GetRange(x1,x2);
				fsp_he_ave->SetRange(x1,x2);

				TF1 *fit_ratio = HistTools::CombineTF1(fsp_he1, fsp_he_ave, HistTools::Divide, "fsp_he2"); // He(R,t) Fit vs. <He(R)> Fit
				f_fit[k] = HistTools::CombineTF1(f_fit[k], fit_ratio, HistTools::Multiply, "f_fit", 0.1, 3e3);
				file_ratio1->Close(); 

			} else if (F34==4){
				
				// load F4 F(R,t)/<F(R,t)> ratio 
				TFile *file_ratio = new TFile(Form("data/ACE/fill/F%d_%s.root", F34, ACE_Element[k])); // load ratio of F(R,t)/<F(R,t) fit
				TF1 *fit_ratio = (TF1*) file_ratio->Get(Form("fit_ratio_BR%d", 2426+iBR_true));	
				f_fit[k] = HistTools::CombineTF1(f_fit[k], fit_ratio, HistTools::Multiply, "f_fit", 0.1, 3e3);

			} 

			fe->AddElementFluxAbove(f_fit[k], A[k]); 

			file0->Close(); 
			file2->Close(); 

		//   }

		} 
 
		TF1 *f = fe->GetTF1Pointer(); 

		if (i==nBRs-1){
			c3->cd(1);
			gPad->SetGrid(); 
			gPad->SetLogx();
			f_RP = fe->GetRSumTF1Pointer(); 

			f_RP->SetTitle("summed ratio as function of rigidity; rigidity (GV); ratio heavier/helium");
			//f_RP->GetXaxis()->SetRangeUser(1.0, 1.15e3);   
			//f_RP->GetYaxis()->SetRangeUser(0.3, 0.625);
			f_RP->Draw();  
		}
		// fe->Print(); 

		double l_max = 1e4; 
		double f_check = 0.; 
		double *bl = HistTools::BuildLogBins(R_cutoff, l_max, 10); // check the integral by separation into few integrals 
			
		for (int k=0; k<10; k++){ 
		
			//f_check += f->Integral(bl[k], bl[k+1]);  
			// if (i==0) printf("k = %d, lower limit = %4.1f, upper limit = %4.1f, f = %10.4f \n", k, bl[k], bl[k+1], f->Integral(R_cutoff, bl[k]));  
				 
		}		

		// printf("Time = %d, N_t = %10.4f,  f/f_check = %10.4f (good if =1) \n", UBRToTime(i+2426), f->Integral(0.1, l_max), f->Integral(0.1, l_max)/f_check);  

		N_t->SetPoint(i, UBRToTime(iBR_true+2426), n_tubes*f->Integral(R_cutoff, l_max)); // this internally makes a loop in the range Rmin to Rmax, and calls FunctorExample::operator() at every step, computing the integral as the sum of all the steps  

		// break; 

	   	if (i+2426==2472-1) iBR_true += 3; 
		else iBR_true ++; 

	}  

	// N_t->Print("range"); 

	c3->Print("data/nm/reproduce2/R_sum_all_template.png"); 

	// PRINT_GRAPH(N_t); 

	TCanvas *c1 = new TCanvas("c1", "Estimated NM Count Rate", 2700, 900); 
	c1->cd(1); 

	TLegend *legend2 = new TLegend(0.1,0.7,0.28,0.9); // left, down, right, top 
	legend2->AddEntry(N_t, "Mi13", "l"); 

	HistTools::SetStyle(N_t, kBlue, kFullCircle, 0.9, 1, 1); 

	N_t->GetXaxis()->SetTimeDisplay(1);
  	N_t->GetXaxis()->SetTimeFormat("%m-%y");
	N_t->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	N_t->SetTitle("; ; Estimated NM Count Rate"); 

	N_t->Draw("APL"); 
	legend2->Draw("SAME"); 

	c1->Print(Form("data/nm/reproduce2/F%d_estimated_nm_count_%s.png", F34, NM)); 

	iBR_true = 0; 
	for (int i=0; i<nBRs; i++){

		double x1, y1, x2, y2; 

		N_t->GetPoint(i, x1, y1);
		N_nm->GetPoint(i, x2, y2); 

		k_sf->SetPoint(i, UBRToTime(iBR_true+2426), y1/y2); 

		// printf("%s, iBR = %d, N_nm(t)=%f, N(t)=%f, k_sf=%10.4f \n", NM, i, y2, y1, y1/y2);

	   	if (i+2426==2472-1) iBR_true += 3; 
		else iBR_true ++;  
	}

	double sum_k = 0., average_k = 0.;  
	for (int i=0; i<nBRs; i++){
		double x_k, y_k; 
		k_sf->GetPoint(i, x_k, y_k); 
		sum_k += y_k; 
	}
	average_k = sum_k/nBRs; 

	TGraph *k_norm = new TGraph(); // normalized k 

	iBR_true = 0; 
	for (int i=0; i<nBRs; i++){
		double x_k, y_k; 
		k_sf->GetPoint(i, x_k, y_k); 
		k_norm->SetPoint(i, UBRToTime(iBR_true+2426), y_k/average_k); 

	   	if (i+2426==2472-1) iBR_true += 3; 
		else iBR_true ++;  
	}

	TCanvas *c2 = new TCanvas("c2", "Estimated NM Scaling Factor (Normalized)", 2700, 900); 

	c2->cd(1); 
	HistTools::SetStyle(k_norm, kBlue, kFullCircle, 0.9, 1, 1); 
	k_norm->GetXaxis()->SetTimeDisplay(1);
  	k_norm->GetXaxis()->SetTimeFormat("%m-%y");
	k_norm->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	k_norm->SetTitle("; ; Estimated NM Scaling Factor (Normalized)"); 
	k_norm->Draw("APL"); 

	c2->Print(Form("data/nm/reproduce2/F%d_estimated_nm_k_%s.png", F34, NM)); 

	// compute BR mean + std

	double k_BR_mean = get_k_BR_mean( k_sf, "mean" );
	double k_BR_std = get_k_BR_mean( k_sf, "std" );  

	printf("(%s, %s) k_BR_mean = %f, k_BR_std = %f \n", NM, option2, k_BR_mean, k_BR_std);  

	TFile *file_reproduce = new TFile(Form("data/nm/reproduce2/F%d_Analysis/F%d_%s_%s.root", F34, F34, NM, option2), "recreate"); 

	k_norm->Write("k_norm");
	k_sf->Write("k_sf"); 
	N_t->Write("N_t"); 
	f_RP->Write("f_RP");

	file_reproduce->Close(); 

	if (!strcmp(option1, "k_norm")) return k_norm; 
	if (!strcmp(option1, "k_sf")) return k_sf; 
	if (!strcmp(option1, "N_t")) return N_t; 
	//if (!strcmp(option1, "YF_divide")) return YF_divide; 

}	

// Reconstruct the NM count with a simple cosmic ray contribution model from Koldobisky's assumption w/ YF Mi13, but separate for p, He, elements above He contributions
TGraph *nm_reproduce2(const char *NM, const char* option1, const char* option2){
	
	Debug::Enable(Debug::ALL); 

	int nnodes_ams = 6; 
	double n_tubes = (double) get_NM_tubes(NM)/6; 
	double R_cutoff = get_Rcutoff(NM); 	

	printf("n_tubes (%s, %s) = %2.2f \n", NM, option2, n_tubes); 

	TGraph *N_t = new TGraph(); 

	Experiments::DataPath = "data"; 

	TF1 *f_fit[n_total]; 

	// TH1 *J_sum = heavier_he(); 
	// J_sum->Print("range"); 

	TF1 *fyfp = NeutronMonitors::YieldFunctions::CreateFunction("fyfp", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mangeard16H1, Particle::PROTON, Energy::RIGIDITY);
	TF1 *fyfhe = NeutronMonitors::YieldFunctions::CreateFunction("fyfhe", 0.1, 3e3, NeutronMonitors::YieldFunctions::Mangeard16He4, Particle::HELIUM4, Energy::RIGIDITY); 

	//HistTools::PrintFunction(fyfp);
	//HistTools::PrintFunction(fyfhe); 

	TF1 *f_BR_p[nBRs];   
	TF1 *f_BR_he[nBRs];  

	TCanvas *c3 = new TCanvas("c3", "R_sum", 1800, 900); 

	TFile *fit_result = new TFile(Form("data/amsfit/fit_result_node%d.root", nnodes_ams)); // ams fits 
	TFile *nm_data = new TFile(Form("./data/nm/NM-%s.root", NM));
	TGraph *N_nm = (TGraph*) nm_data->Get("g_ave"); 

	TGraph *k_sf = new TGraph(); // k = N_t/N_nm, scaling factor of NM stations 

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

		//FunctorExample *fe = new FunctorExample("fe", 0.1, 3e3, fyfp, fyfhe); 
		FunctorExample *fe_p = new FunctorExample("fe_p", 0.1, 3e3, fyfp, fyfhe); // p only 
		FunctorExample *fe_he = new FunctorExample("fe_he", 0.1, 3e3, fyfp, fyfhe); // he only 
		FunctorExample *fe_above = new FunctorExample("fe_above", 0.1, 3e3, fyfp, fyfhe); // all other elements above only 
	
		for(int k=0; k<n_ams; k++){

			Spline *sp_ams = new Spline("sp_ams", nnodes_ams, Spline::LogLog | Spline::PowerLaw); 
			f_fit[k] = sp_ams->GetTF1Pointer();  // real function 
			TF1 *fit_ams = (TF1*) fit_result->Get(Form("fsp_%s", AMS_Element[k]))->Clone(Form("f_%s", AMS_Element[k])); 
	
			HistTools::CopyParameters(fit_ams, f_fit[k]); 
			double x1, x2;
			fit_ams->GetRange(x1,x2); 
			f_fit[k]->SetRange(x1,x2); 
 
			//fe->AddElementFlux_Ave(f_fit[k], A[k]); 
			fe_p->AddElementFlux_Ave(f_fit[k], A[k]); 
			fe_he->AddElementFlux_Ave(f_fit[k], A[k]);
			fe_above->AddElementFlux_Ave(f_fit[k], A[k]); 
			// if (i==nBRs-1) HistTools::PrintFunction(fit_ams);   
		} 

		TF1 *f0 = new TF1("f0", "0", 0.1, 3e3);  

		fe_p->SetProtonFlux(f_BR_p[i]); 
		fe_p->AddElementFlux(f0, A[1]);

		fe_he->SetProtonFlux(f0); 
		fe_he->AddElementFlux(f_BR_he[i], A[1]);

		fe_above->SetProtonFlux(f0); 
		fe_above->SetHeliumFluxReplicate(f_BR_he[i]); 
		fe_above->AddElementFlux(f0, A[1]);

		TF1 *f;
		if (!strcmp(option2, "p")) f = fe_p->GetTF1Pointer();
		if (!strcmp(option2, "He")) f = fe_he->GetTF1Pointer();
		if (!strcmp(option2, "above")) f = fe_above->GetTF1Pointer();		
 
/*
		if (i==nBRs-1){
			c3->cd(1);
			gPad->SetLogx();
			TF1 *f_RP = fe->GetRSumTF1Pointer(); 

			f_RP->SetTitle("summed ratio as function of rigidity; rigidity (GV); ratio heavier/helium");  
			f_RP->GetXaxis()->SetRangeUser(1.0, 1.15e3);  
			f_RP->GetYaxis()->SetRangeUser(0.3, 0.625); 
			f_RP->Draw(); 
		}
		// fe->Print(); 
*/

		double l_max = 1e4; 
		double f_check = 0.; 
		double *bl = HistTools::BuildLogBins(R_cutoff, l_max, 10); // check the integral by separation into few integrals 
			
		for (int k=0; k<10; k++){ 
		
			f_check += f->Integral(bl[k], bl[k+1]);  
			// if (i==0) printf("k = %d, lower limit = %4.1f, upper limit = %4.1f, f = %10.4f \n", k, bl[k], bl[k+1], f->Integral(R_cutoff, bl[k]));  
				 
		}		

		// printf("Time = %d, N_t = %10.4f,  f/f_check = %10.4f (good if =1) \n", UBRToTime(i+2426), f->Integral(0.1, l_max), f->Integral(0.1, l_max)/f_check);   

		N_t->SetPoint(i, UBRToTime(i+2426), n_tubes*f->Integral(R_cutoff, l_max)); // this internally makes a loop in the range Rmin to Rmax, and calls FunctorExample::operator() at every step, computing the integral as the sum of all the steps  

		// break; 

	} 

	// c3->Print("data/nm/reproduce2/ratio_heavier_helium.png"); 

	// PRINT_GRAPH(N_t); 

	TCanvas *c1 = new TCanvas("c1", "Estimated NM Count Rate", 2700, 900); 
	c1->cd(1); 

	TLegend *legend2 = new TLegend(0.1,0.7,0.28,0.9); // left, down, right, top 
	legend2->AddEntry(N_t, "Mi13", "l"); 

	HistTools::SetStyle(N_t, kBlue, kFullCircle, 0.9, 1, 1); 

	N_t->GetXaxis()->SetTimeDisplay(1);
  	N_t->GetXaxis()->SetTimeFormat("%m-%y");
	N_t->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	N_t->SetTitle("; ; Estimated NM Count Rate"); 

	N_t->Draw("APL"); 
	legend2->Draw("SAME"); 

	// c1->Print(Form("data/nm/reproduce2/estimated_nm_count_%s.png", NM)); 

	for (int i=0; i<nBRs; i++){

		double x1, y1, x2, y2; 

		N_t->GetPoint(i, x1, y1);
		N_nm->GetPoint(i, x2, y2); 

		k_sf->SetPoint(i, UBRToTime(i+2426), y1/y2); 

		// printf("iBR = %d, N_nm(t)=%10.4f, N(t)=%10.4f, k_sf=%10.4f \n", i, y2, y1, y1/y2); 
	}

	double sum_k = 0., average_k = 0.;  
	for (int i=0; i<nBRs; i++){
		double x_k, y_k; 
		k_sf->GetPoint(i, x_k, y_k); 
		sum_k += y_k; 
	}
	average_k = sum_k/nBRs; 

	TGraph *k_norm = new TGraph(); // normalized k 

	for (int i=0; i<nBRs; i++){
		double x_k, y_k; 
		k_sf->GetPoint(i, x_k, y_k); 
		k_norm->SetPoint(i, UBRToTime(i+2426), y_k/average_k);  
	}

	TCanvas *c2 = new TCanvas("c2", "Estimated NM Scaling Factor (Normalized)", 2700, 900); 

	c2->cd(1); 
	HistTools::SetStyle(k_norm, kBlue, kFullCircle, 0.9, 1, 1); 
	k_norm->GetXaxis()->SetTimeDisplay(1);
  	k_norm->GetXaxis()->SetTimeFormat("%m-%y");
	k_norm->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	k_norm->SetTitle("; ; Estimated NM Scaling Factor (Normalized)"); 
	k_norm->Draw("APL"); 

	// c2->Print(Form("data/nm/reproduce2/estimated_nm_k_%s.png", NM)); 

	if (!strcmp(option1, "k_norm")) return k_norm; 
	if (!strcmp(option1, "k_sf")) return k_sf; 
	if (!strcmp(option1, "N_t")) return N_t;  
}

// Return ratio heavier/helium
TH1 *heavier_he(){

	int nnodes_ams = 6; 

	gSystem->mkdir("data/nm/reproduce2", true); 

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

	// J_sum->Print("range"); 

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

	// c1->Print("data/nm/reproduce2/ratio_heavier_helium.png"); 

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

   TFile *_file2 = new TFile(Form("data/ACE/%s_97226_18359.root",  element)); 
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

double get_Rcutoff(const char *NM){	
	for (int i=0; i<nNMs_Koldob+1; i++){
		if (!strcmp(NM, Form("%s", NM_Koldob[i]))) return R_cutoff[i];
	}
}

TGraphErrors *get_k_sf_ratio( TGraphErrors *k_sf1, TGraphErrors *k_sf2 ){

	TGraphErrors *k_sf_ratio = new TGraphErrors(); 

	for (int iBR=0; iBR<nBRs; iBR++){

		Double_t x1=0, y1=0, // Ling's 
			 x2=0, y2=0; // Koldob 		
		
		k_sf1->GetPoint(iBR, x1, y1); 
		k_sf2->GetPoint(iBR, x2, y2);

		k_sf_ratio->SetPoint(iBR, x1, y1/y2); 

		// error propogation 
		Double_t dy1 = k_sf1->GetErrorY(iBR)==-1.0 ? 0 : k_sf1->GetErrorY(iBR);  
		Double_t dy2 = k_sf2->GetErrorY(iBR)==-1.0 ? 0 : k_sf2->GetErrorY(iBR); 

		k_sf_ratio->SetPoint(iBR, x1, y1/y2);			 
		k_sf_ratio->SetPointError(iBR, 0, (y1/y2)*sqrt((dy1/y1)*(dy1/y1)+(dy2/y2)*(dy2/y2))); 

		// printf("iBR = %d, dy1 = %f, dy2 = %f \n", iBR, dy1, dy2); 

	} 

	return k_sf_ratio; 	
} 

double get_k_BR_mean(TGraph *k_sf, const char *option){

	double mean = 0, std = 0;

	for (int iBR=0; iBR<nBRs; iBR++){
	
		Double_t x=0, y=0; 

		k_sf->GetPoint(iBR, x, y);
		mean += y; 

	} 
	
	mean = mean/nBRs; 

	for (int iBR=0; iBR<nBRs; iBR++){

		Double_t x=0, y=0; 

		k_sf->GetPoint(iBR, x, y);
		std += pow((y-mean), 2);	

	}

	std = std/(nBRs-1); 	
	std = sqrt(std/nBRs);

	if (!strcmp(option, "mean")) return mean; 
	if (!strcmp(option, "std")) return std; 
}

int get_NM_tubes(const char *NM){

	for (int i=0; i<nNMs_Koldob+1; i++){
		if (!strcmp(NM, Form("%s", NM_Koldob[i]))) return NM_tubes[i];
	}

}
