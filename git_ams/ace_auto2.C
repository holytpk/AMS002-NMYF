// Manipulation of ACE/CRIS data by Lingqiang He, updated on 09/05/2019
// for i in {5..28}; do sed -r -e'1,5d' -e's/^ +//' -e's/ +/ /g' cris_energy_bands.txt | egrep "^$((i))" | cut -d' ' -f3,4 | paste -sd' ' | sed -r 's/ /, /g' >> ace_energy_band.dat; done

//CRIS: BR 2240-2529

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

using namespace std;

const int n_ele = 24; // number of ACE elements 
const int n_total = 28; // number of total elements 

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

int get_index(const char *element){

	int index = 0; 
	for (int i=0; i<n_total; i++){
		if (!strcmp(element, Form("%s", Element[i]))) index = i;
	}

	return index; 
}

int get_temp_int(const char *element){

	int temp = 0; 
	for (int i=0; i<n_total; i++){
		if (!strcmp(element, Form("%s", Element[i]))) temp = i_temp[i];
	}

	return temp; 
}

double pow(double factor){

	double product = factor*factor;

	return product; 
}

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
const char *get_template(const char *element); 
int get_ams_data_value(const char *element); 
TGraph *get_norm_graph(TGraph *g); 
TGraphErrors *get_norm_graph(TGraphErrors *g); 

void ace_fill(const char *element, Particle::Type isotope); // fill ACE/CRIS data into root histograms, then convert to GV/m^2 and save 
TGraphErrors *ace_rescale_BR(const char *element, Particle::Type isotope); // rescale ACE BR points to the level of ACE BR Average, returns g_ratio_time[iset]
void ace_rescale_BR_averaged(const char *element, Particle::Type isotope); // improved ver. for BCNO. Plots spectral indices vs. rigidity and time. y
void ace_fake_td_ams(const char *element, Particle::Type isotope); // make a fake time-dependent TOA Flux Model with (AMS Helium BR)*(<AMS>_C/<AMS>_He)   
void ace_fake_td_ams_v2(const char *element, Particle::Type isotope); // fake_td_ams but plot He(R,t) Fit vs. <He(R)> Fit
void ace_fake_td_ams_v3(const char *element, Particle::Type isotope); // use O pars for C, B pars for N. 
void ace_fake_td_ams_v4(const char *element, Particle::Type isotope); // par0==linear fit, par1==average over time; use O pars for C, B pars for N. 
void plot_fit_pars(); // plot fit par_BCN vs. par_O  
void compare_fake_flux( const char *element, Particle::Type isotope ); // plot F1, F2, F3 

void ace_all_average(); // plot averaged flux over energy bins for all elements
void ace_convert(const char *element, Particle::Type isotope); // convert h_ene into h_rig and also plot flux and normalized flux over time
void ace_fitboth(int nnodes); // fit ACE&AMS combined data 
void compare_nodes(int k); // compare among spline fit results for 5-9 nodes 
double compare_sig(TF1 *fit1, TF1 *fit2); // obtain sigma value between two fits
void ace_extend(); // divide remaining ACE data that is not measured by AMS by combined fit to find a flat residual for ACE assumption at high energy region
		   // divide elements with Z>8 ACE data by combinced C cit 
void ace_extend2(); // fit_comb of BNO / C
void ace_extend3(); // group ratios from extend() part II and check isotope assumptions 
void ace_extend4(); // plot average relative absolute difference, chi-2, max variation of residuals vs. element Z, element Z/A, difference (Z/A)_element - (Z/A)_C
void ace_extend5(int temp); // extend low energy profile of p, He, Li, Be AMS data through BCNO templates and rescaled few times 
// please input template, B = 0, C = 1, N = 2, o = 3 
void ace_contribution(); // plot the contributions by all estimated element fluxes to the total cosmic ray flux   
void ace_contribution2(); // plot the optimized contribution (proper template for each element) table  
TH1* ace_contribute(int temp, const char *option); // reconstructe ACE fit for all elements through C combined fit template, check the contribution of each element by plotting the ratio of CR Flux/total CR Flux 
// option "n" - normal total; option "mw" - mass-weighted total.
void ace_contribute2(); // recompute the total cosmic ray flux with selected template for each element according to its category   

TH1* ave_hist( TH1D **h_set, int nBRs ); // average a stack of histograms 
TGraphErrors *ave_grapherrors( TGraphErrors *g_set[], int nset ); // average a stack of TGraphErrors

// double *getmean_ace( TGraphErrors *g, const char *option )
TGraphAsymmErrors *get_ace_graph(const char *element, UInt_t iBin, UInt_t nBRs); // flux in Kinetic Energy over time 
TGraphAsymmErrors *get_ace_average_graph(const char *element, UInt_t *BRs, UInt_t nBRs); // flux in Kinetic over energy bins 

void ace_auto(const char *operation){
 
	Experiments::DataPath = "data";

	//gROOT->SetBatch();
	//gROOT->ProcessLine(".L ace_auto.C");

	gSystem->mkdir("data/ACE/fill", true);
	gSystem->mkdir("data/ACE/extend", true);
	gSystem->mkdir("data/ACE/extend2", true);
	gSystem->mkdir("data/ACE/compare", true);
	gSystem->mkdir("data/ACE/convert", true);
	gSystem->mkdir("data/ACE/contribute", true);

	gSystem->Setenv("TZ", "UCT"); 
	gStyle->SetTimeOffset(0);
	gStyle->SetOptStat(0); 
	gStyle->SetNumberContours(99); 
	gStyle->SetPalette(55); 

	const int nAce = 28-5+1; 
	// const char *elements[nAce] = { "B11", "C12", "N15", "O16", "F19", "Ne20", "Na23", "Mg24", "Al27", "Si28", "P31", "S32", "Cl35", "Ar36", "K41", "Ca40", "Sc45", "Ti46", "Va51", "Cr52", "Mn55", "Fe56", "Co59", "Ni60" }; list of most abundant isotopes
	
	if (strcmp(operation, "fill") == 0){

		// gROOT->ProcessLine(".> data/ACE/fill/fill_all.txt"); 

		for (int i=0; i<24;i++){
			//ace_fill( ACE_Element[i], ACE_Isotope[i] );
		} 

		for (int i=0; i<4;i++){
			//ace_rescale_BR( ACE_Element[i], ACE_Isotope[i] ); 
			//ace_rescale_BR_averaged( ACE_Element[i], ACE_Isotope[i] ); 
			//ace_fake_td_ams( ACE_Element[i], ACE_Isotope[i] );
			//ace_fake_td_ams_v2( ACE_Element[i], ACE_Isotope[i] );
			//ace_fake_td_ams_v3( ACE_Element[i], ACE_Isotope[i] ); 
			//ace_fake_td_ams_v4( ACE_Element[i], ACE_Isotope[i] ); 

			//compare_fake_flux( ACE_Element[i], ACE_Isotope[i] ); 
		}

		for (int i=4; i<24; i++){

			if (i==22) continue; 
			else {	
				//ace_rescale_BR( ACE_Element[i], ACE_Isotope[i] ); 
				//ace_rescale_BR_averaged( ACE_Element[i], ACE_Isotope[i] ); 
				//ace_fake_td_ams( ACE_Element[i], ACE_Isotope[i] );
				ace_fake_td_ams_v2( ACE_Element[i], ACE_Isotope[i] );
				//ace_fake_td_ams_v3( ACE_Element[i], ACE_Isotope[i] ); 
				//ace_fake_td_ams_v4( ACE_Element[i], ACE_Isotope[i] ); 

				//compare_fake_flux( ACE_Element[i], ACE_Isotope[i] ); 
			} 
		}
		//compare_fake_flux( ACE_Element[7], ACE_Isotope[7] ); 

		// TLegend *legend = new TLegend(0.1,0.6,0.28,0.9); 

		/*
		for (int i=0; i<n_ele; ++i){	

			TGraph *g_ratio_ave = new TGraph(Form("./data/ACE/fill/rescaling_factor_%s.dat", ACE_Element[i])); 
			HistTools::SetStyle(g_ratio_ave, HistTools::GetColorPalette(i, n_ele), kFullCircle, 0.65, 1, 1); 

			TGraph *g_ratio_ave_norm = get_norm_graph(g_ratio_ave); 

			g_ratio_ave->Print("range"); 

			g_ratio_ave_norm->GetYaxis()->SetRangeUser(0.007, 2.5); 
			g_ratio_ave_norm->SetTitle("Normalized Time-dependent Rescaling Factor; ; "); 

			g_ratio_ave_norm->GetXaxis()->SetTimeDisplay(1);
			g_ratio_ave_norm->GetXaxis()->SetTimeFormat("%Y-%m");
			g_ratio_ave_norm->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g_ratio_ave_norm->GetXaxis()->SetTitleSize(0.7);

			legend->AddEntry(g_ratio_ave_norm, Form("%s", ACE_Element[i]), "pl");
			legend->SetNColumns(3); 

			if (i==4){ 
				gPad->SetGrid();
				g_ratio_ave_norm->Draw("APL"); 
				legend->Draw("SAME");  
			}
			if (i>4){
				gPad->SetGrid();
				g_ratio_ave_norm->Draw("PLSAME");
				legend->Draw("SAME"); 
			} 
		}
		
		c0->Print("./data/ACE/fill/time_dependent_rescaling_factor.png");
		*/
  
		// gROOT->ProcessLine(".> ");

	} else if (strcmp(operation, "convert") == 0){

		gROOT->ProcessLine(".> data/ACE/convert/rigidity_all.txt");
	
		ace_convert( "B", Particle::BORON11 );
		ace_convert( "C", Particle::CARBON12 );
		ace_convert( "N", Particle::NITROGEN15 );
		ace_convert( "O", Particle::OXYGEN16 );
		ace_convert( "F", Particle::FLUORINE19 );
		ace_convert( "Ne", Particle::NEON20 );
		ace_convert( "Na", Particle::SODIUM23 );
		ace_convert( "Mg", Particle::MAGNESIUM24 );
		ace_convert( "Al", Particle::ALUMINUM27 );
		ace_convert( "Si", Particle::SILICON28 );
		ace_convert( "P", Particle::PHOSPHORUS31 );
		ace_convert( "S", Particle::SULFUR32 );
		ace_convert( "Cl", Particle::CHLORINE35 );
		ace_convert( "Ar", Particle::ARGON36 );
		ace_convert( "K", Particle::POTASSIUM41 );
		ace_convert( "Ca", Particle::CALCIUM40 );
		ace_convert( "Sc", Particle::SCANDIUM45 );
		ace_convert( "Ti", Particle::TITANIUM46 );
		ace_convert( "Va", Particle::VANADIUM51 );
		ace_convert( "Cr", Particle::CHROMIUM52 );
		ace_convert( "Mn", Particle::MANGANESE55 );
		ace_convert( "Fe", Particle::IRON56 );
		ace_convert( "Co", Particle::COBALT59 );
		ace_convert( "Ni", Particle::NICKEL60 ); 
	
		gROOT->ProcessLine(".> ");
		//gROOT->Reset(); 
	} else if (strcmp(operation, "average") == 0){
		ace_all_average();
	} else if (strcmp(operation, "compare") == 0){ 

	   gROOT->ProcessLine(".> data/ACE/compare/sigmatest_allnodes.txt");
	   for(int k=0;k<4;k++){
		compare_nodes(k);
	   }
	   gROOT->ProcessLine(".> "); 

	} else if (strcmp(operation, "fitboth") == 0){

	   gROOT->ProcessLine(".> data/ACE/compare/allnodes.txt");
	   for(int i=5;i<=9;i++){
		//gROOT->ProcessLine(Form(".> data/ACE/fit_both_%dnodes.txt", nnodes));
		int nnodes= i; // 8 nodes for best stats  
		
		ace_fitboth(nnodes);
	   }
	   gROOT->ProcessLine(".> "); 
	} else if (strcmp(operation, "extend") == 0){
		ace_extend();
	} else if (strcmp(operation, "extend2") == 0){
		ace_extend2();
	} else if (strcmp(operation, "extend3") == 0){
		ace_extend3();
	} else if (strcmp(operation, "extend4") == 0){
		ace_extend4();
	} else if (strcmp(operation, "compare_fake") == 0){

		for (int i=0; i<4;i++){ 
			compare_fake_flux( ACE_Element[i], ACE_Isotope[i] ); 
		} 

	} 
} 

// Fill CRIS Data into Bins
void ace_fill(const char *element, Particle::Type isotope){

	double *kin_bins = get_kin_bins(element);
        double *SpallCorr = get_spall_corr(element);
	double *SpallCorrUnc = get_spall_corr_unc(element);
	double *EMed = get_EMed(element);

	gSystem->mkdir("data/ACE/fill", true);
	
	TFile *_file0 = new TFile(Form("data/ACE/%s_97226_18359.root", element));
	TTree *ace=(TTree*)_file0->Get("ace");
	
	float F[7], C[7]; // must be initialized 
	Long64_t utime;
	Float_t livetime;
	
	ace->SetBranchAddress("F", F);
	ace->SetBranchAddress("C", C);
	ace->SetBranchAddress("start_utime", &utime);	
	ace->SetBranchAddress("livetime", &livetime);		

	TFile file1(Form("data/ACE/fill/%s_fill.root", element), "RECREATE");

	TCanvas *c0 = new TCanvas("c0", "", 800, 600);
	c0->cd(1);
	// fill the bins for each BR 
	for (int k=0; k<ace->GetEntries(); ++k){

		TH1F *h_ene = new TH1F("h_ene","", 13, kin_bins);		

		ace->GetEntry(k); //Get the nth entry of TTree!! 						

		for (int i=0; i<h_ene->GetNbinsX(); ++i) { 
			if (i%2==0) { 
				double sys_err = F[i/2] * sqrt(8e-4 + SpallCorrUnc[i/2]*SpallCorrUnc[i/2]/SpallCorr[i/2]/SpallCorr[i/2]);
				double stat_err = F[i/2]/sqrt(C[i/2]); 
				double tot_err = sqrt(stat_err*stat_err + sys_err*sys_err);
				h_ene->SetBinContent(i+1, F[i/2]); 

				if (C[i/2]!=0){
					h_ene->SetBinError(i+1, tot_err); 
				} else if (C[i/2]==0){
					continue;
				}
			} 
		}

		h_ene->SetMarkerStyle(kFullCircle);
		//HistTools::SetColors(h_ene, 290, kFullCircle, 1.4);
		gPad->SetGrid(); 
		h_ene->SetTitle(Form("%s BR-%d Energy Spectrum; Energy (MeV/nuc); Flux (/(cm^2 sr s)(MeV/nuc)", element, UTimeToBR(utime)));
		//h_ene->Draw("PSAME"); 

		TH1 *h_rig = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene, isotope, "MeV/n cm", "GV m", "_rig");

		h_rig->SetTitle(Form("%s BR-%d Rigidity Spectrum; ; ", element, UTimeToBR(utime)));
		h_rig->SetXTitle(Unit::GetEnergyLabel("GV"));
  		h_rig->SetYTitle(Unit::GetDifferentialFluxLabel("GV m"));

		h_ene->Write(Form("h_kin_%s_BR%d", element, UTimeToBR(utime)));
		h_rig->Write(Form("h_rig_%s_BR%d", element, UTimeToBR(utime))); 

	}

	file1.Write();
	file1.Close();

	_file0->Close();

}

// & Rescale ACE BR points to the level of ACE BR Average (for convinience, I just merge these functions) 
TGraphErrors *ace_rescale_BR(const char *element, Particle::Type isotope){

	gStyle->SetOptStat(0);

 	Experiments::DataPath = "data";	
   	int data_value[n_ele] = {0, 18, 23, 27, 31, 20, 43, 22}; // p, He, Li, Be, B, C, N, O

	double *kin_bins = get_kin_bins(element);
        double *SpallCorr = get_spall_corr(element);
	double *SpallCorrUnc = get_spall_corr_unc(element);
	double *EMed = get_EMed(element);

	const UInt_t FirstACEBR = 2240;
	vector<UInt_t> BRs;
	// we stop at BR 2493, which ends on 2016/05/23, just 3 days before the end of the data taking period for AMS nuclei
	for (UInt_t br=2426; br<=2493; ++br) { 
		if (br != 2472 && br != 2473) BRs.push_back(br-FirstACEBR); 
	}

	TH1 *h_ams = Experiments::GetMeasurementHistogram(Experiments::AMS02, data_value[1], 0); // load AMS data for a given element
	TH1 *h_ace_ene = HistTools::GraphToHist(get_ace_average_graph( element , &BRs[0], BRs.size() ), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
	TH1 *h_ace = HistTools::TransformEnergyAndDifferentialFluxNew(h_ace_ene, isotope, "MeV/n cm", "GV m", "_rig"); // load averaged ACE data for the same element, converted in rigidity

	UShort_t namsbins = h_ams->GetNbinsX();
	UShort_t nacebins = h_ace->GetNbinsX();
	double R1 = h_ace->GetBinLowEdge(1);
	double R2 = h_ams->GetBinLowEdge(namsbins+1);

	int nnodes = 7;

	TFile *file0 = new TFile(Form("data/ACE/compare/fit_%s_%dnodes.root", get_template(element), nnodes));  

	gSystem->mkdir("data/ACE/fill", true);
	
	TFile *_file0 = new TFile(Form("data/ACE/%s_97226_18359.root", element));
	TTree *ace=(TTree*)_file0->Get("ace");
	
	float F[7], C[7]; // must be initialized 
	Long64_t utime;
	Float_t livetime;
	
	ace->SetBranchAddress("F", F);
	ace->SetBranchAddress("C", C);
	ace->SetBranchAddress("start_utime", &utime);	
	ace->SetBranchAddress("livetime", &livetime);		

	TCanvas *c0 = new TCanvas("c0", "", 800, 600);
	c0->cd(1);
	// fill the bins for each BR 

	TH1 *h_a1 = new TH1D("", "", 200, 0, 200); 
	TH1 *h_a2 = new TH1D("", "", 200, 0, 200); 

	// gROOT->ProcessLine(Form(".> data/ACE/fill/rescaling_factor_%s.dat", element)); 

	TCanvas *c1 = new TCanvas("c1", "", 1600, 900); 
	c1->Divide(1, 2);  

	TFile *file3 = new TFile(Form("data/ACE/fill/%s_fill.root", element)); 
	TFile *file1 = new TFile(Form("data/ACE/fill/F1_%s.root", element), "RECREATE");

	TGraphErrors *g_ratio_time[7]; 
	TGraphErrors *g_ratio_time_norm[7]; 
	TGraphErrors *g_ratio_ace_fit[7];
	// TGraphErrors *g_ratio_ams_fit[7];

	for (int i=0; i<7; ++i){
		g_ratio_time[i] = new TGraphErrors(); 
	}

	TGraphErrors *g_ratio_ace[7];
	TGraphErrors *g_ratio_ams[7]; 
	TGraphErrors *g_residual_ace[7];
	TGraphErrors *g_resierr_ace[7]; 
	TGraphErrors *g_residual_norm_ace[7];
	TGraphErrors *g_residual_ams[h_ams->GetNbinsX()]; 
	//TGraphErrors *g_resierr_ams[h_ams->GetNbinsX()]; 
	TGraphErrors *g_residual_norm_ams[h_ams->GetNbinsX()];

	for (int bin=0; bin <= h_ace->GetNbinsX(); ++bin){ 
		if (bin%2==0){  
			g_ratio_ace[bin/2] = new TGraphErrors(); 
			g_ratio_ace_fit[bin/2] = new TGraphErrors(); 
			g_residual_ace[bin/2] = new TGraphErrors();  
			g_resierr_ace[bin/2] = new TGraphErrors(); 
		} 
	}  

	for (int k=0; k<ace->GetEntries(); ++k){

		ace->GetEntry(k); //Get the nth entry of TTree!! 						

		TH1 *h_rig = (TH1*) file3->Get(Form("h_rig_%s_BR%d", element, UTimeToBR(utime))); 
			
		// average of h_rig 
		double rig_sum=0; // compute average of h_rig manually  
		for(int nbin=0;nbin<14;++nbin){
			rig_sum += h_rig->GetBinContent(nbin);
			//printf("ratio_sum = %0.6f \n", ratio_sum);
		}
		double rig_ave = rig_sum/7;

		// printf("rig_ave = %10.4f \n", rig_ave);  

		h_rig->SetTitle(Form("%s BR-%d Rigidity Spectrum; ; ", element, UTimeToBR(utime)));
		h_rig->SetXTitle(Unit::GetEnergyLabel("GV"));
  		h_rig->SetYTitle(Unit::GetDifferentialFluxLabel("GV m"));

		// h_ene->Write(Form("h_kin_%s_BR%d", element, UTimeToBR(utime)));
		// h_rig->Write(Form("h_rig_%s_BR%d", element, UTimeToBR(utime))); 

		// rescale ACE BR to ACE Averaged Magnitude 
		Spline *sp_comb = new Spline("sp_comb", nnodes, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_comb = sp_comb->GetTF1Pointer();  
		TF1 *fit_comb = (TF1*) file0->Get("fit_both"); 

		HistTools::CopyParameters(fit_comb, fsp_comb); // error 
		double x1, x2;
		fit_comb->GetRange(x1,x2);
		fsp_comb->SetRange(x1,x2);

		double *nodes = sp_comb->GetXNodes(); 

		TLine *l1 = new TLine( nodes[0], 0, nodes[0], h_rig->GetMaximum()*2); 
		TLine *l2 = new TLine( nodes[1], 0, nodes[1], h_rig->GetMaximum()*2); 

		TLegend *legend = new TLegend(0.62,0.8,0.9,0.9); 

		l1->SetLineColor(kGreen-3);
		l2->SetLineColor(kGreen-3);
		l1->SetLineStyle(2);
		l2->SetLineStyle(2); 		

		TH1 *h_ratio = (TH1 *) h_rig->Clone("h_ratio");

		h_ratio->Divide(fsp_comb);

		HistTools::SetStyle(h_ratio, kRed, kFullCircle, 0.9, 1, 1);

		double ratio_sum=0; // compute average of h_ratio manually  
		for(int nbin=0;nbin<14;++nbin){
			ratio_sum += h_ratio->GetBinContent(nbin);
			//printf("ratio_sum = %0.6f \n", ratio_sum);
		}
		double ratio_ave = ratio_sum/7;

		// printf("%d %0.6f \n", (UInt_t) utime, ratio_ave);
		//HistTools::PrintFunction(fit_comb);
			
		double scale = 1./ratio_ave;
		//if (get_index(element)<4) h_ratio->Scale(scale);	

		TH1 *h_rig_rescale = (TH1 *) h_rig->Clone("h_rig_rescale");
		h_rig_rescale->Divide(fsp_comb); 
		HistTools::SetStyle(h_rig_rescale, kRed, kFullCircle, 1.4, 1, 1);

		TF1 *rescaled_fit = HistTools::CombineTF1Const(fsp_comb, ratio_ave, HistTools::MultiplyConst, "rescaled_fit", R1, R2); 
		TF1 *fit_ratio = HistTools::CombineTF1Const(fsp_comb, ratio_ave, HistTools::MultiplyConst, "rescaled_fit", R1, R2); 
		fit_ratio = HistTools::CombineTF1(fit_ratio, fsp_comb, HistTools::Divide, "fit_ratio", R1, R2);  

		// rescaled_fit->Print(); 

		TObjArray data; 
		data.Add(h_rig_rescale);
		// data.Add(h_ratio0);

		// fit AMS and ACE data at the same time
		vector<double> rigmin, rigmax, chi2norm(2);
		ROOT::Fit::Fitter fitter;
		//HistTools::PrintFunction(fit);
		FitTools::SetCommonFitterOptions(fitter);
		FitTools::FitCombinedData(data, fit_ratio, "I", rigmin, rigmax, chi2norm, fitter, 3); 

		TH1 *h_fitres[2]; 
		TH1 *h_fiterr[2]; 

		for (UShort_t i = 0; i < 1; ++i)
		{
  			TH1 *hist = HistTools::ToHist(data[i]);
			h_fiterr[i] = (TH1*) HistTools::GetFitError( (TH1*) data[i], fit_ratio, "_fiterr", false, false, true, 10, DBL_MAX, 0.68, &fitter);
			//h_fiterr[i]->Print("range");  
   			h_fitres[i] = (TH1 *) HistTools::GetResiduals(hist, fit_ratio, "_fitres", false, false, true, 5, 1, 0.68, &fitter);
			HistTools::CopyStyle(hist, h_fitres[i]);

   			UShort_t ndf  = hist->GetNbinsX();
   			Double_t chi2 = chi2norm[i]*ndf; 
   			// printf(" %d nodes ### %-34s   chi2/ndf=%6.2f/%-2u   chi2norm=%5.2f   prob.=%5.2f%%\n", nnodes, hist->GetTitle(), chi2, ndf, chi2norm[i], TMath::Prob(chi2, ndf)*1e2);	
			
			if (i==0){ 
				//g_chi2_ace->SetPoint(iBR, UBRToTime(iBR_true+2426), chi2/ndf);
				for (int bin=0; bin <= h_fitres[i]->GetNbinsX(); ++bin){ 
					if (bin%2==0){ 
					 
						g_ratio_time[bin/2]->SetPoint(k, utime, h_ratio->GetBinContent(bin+1));
						g_ratio_time[bin/2]->SetPointError(k, 0, h_ratio->GetBinError(bin+1));
						g_ratio_ace[bin/2]->SetPoint(k, utime, h_rig_rescale->GetBinContent(bin+1));
						g_ratio_ace[bin/2]->SetPointError(k, 0, h_rig_rescale->GetBinError(bin+1)); 
						g_ratio_ace_fit[bin/2]->SetPoint(k, utime, fit_ratio->Eval(h_rig_rescale->GetBinCenter(bin+1))); 
						g_ratio_ace_fit[bin/2]->SetPointError(k, 0, h_fiterr[i]->GetBinError(bin+1)); 
						g_residual_ace[bin/2]->SetPoint(k, utime, h_fitres[i]->GetBinContent(bin+1));
						g_residual_ace[bin/2]->SetPointError(k, 0, h_fitres[i]->GetBinError(bin+1)); 

						double df = h_fitres[i]->GetBinContent(bin+1)*hist->GetBinContent(bin+1)/hist->GetBinError(bin+1); 
						//printf(" df = %0.4f, a= %0.4f, b=%0.4f \n", df, h_fitres[i]->GetBinContent(bin+1)*h_fiterr[i]->GetBinContent(bin+1), h_fiterr[i]->GetBinError(bin+1)); 
						double ddf = df*sqrt(pow(h_fitres[i]->GetBinError(bin+1)/h_fitres[i]->GetBinContent(bin+1))+pow(h_fiterr[i]->GetBinError(bin+1)/h_fiterr[i]->GetBinContent(bin+1))); // d(data-fit) 
						g_resierr_ace[bin/2]->SetPoint(k, utime, df); 
						g_resierr_ace[bin/2]->SetPointError(k, 0, ddf); 

						// printf("bin/2 = %d, iBR = %d, time = %10.4u, fitres = %10.4f \n", bin/2, iBR, UBRToTime(iBR_true+2426), h_fitres[0]->GetBinContent(bin+1)); 
					} 
				}   
			} 
		
		}
		
		c1->cd(1);
		gPad->SetGrid(); 
		//gPad->SetLogy();	

		h_a1->SetTitle(Form("%s BR-%d Rigidity Spectrum; ; ", element, UTimeToBR(utime)));
		h_a1->SetXTitle(Unit::GetEnergyLabel("GV"));
  		h_a1->SetYTitle(Unit::GetDifferentialFluxLabel("GV m"));

		//if (get_temp_int(element)<4) h_a1->GetYaxis()->SetRangeUser(h_rig->GetMinimum()*0.5, h_rig->GetMaximum()*2); 
		h_a1->GetXaxis()->SetRangeUser(0., 3.); 	

		//h_rig->GetYaxis()->SetRangeUser(0.0001, 2.0); 
		h_rig->GetXaxis()->SetRangeUser(0., 3.);

		legend->AddEntry(rescaled_fit, "rescaled model", "l"); 
		legend->AddEntry(h_rig, "ACE flux", "p"); 
	
		h_a1->Draw("E1X0");  
		h_rig->Draw("E1X0 SAME"); 

		//rescaled_fit->SetRange(0., 3.); 
		rescaled_fit->Draw("SAME"); 
		legend->Draw("SAME"); 
		l1->Draw();
		l2->Draw();   

		c1->cd(2); 
		gPad->SetGrid(); 

		//if (get_temp_int(element)<4) h_a2->GetYaxis()->SetRangeUser(0, 2); 
		h_a2->GetXaxis()->SetRangeUser(0., 3.); 

		h_a2->SetTitle("");
		h_a2->SetXTitle(Unit::GetEnergyLabel("GV"));
		h_a2->SetYTitle(Form("%s BR-%d Flux / %s Model", element, UTimeToBR(utime), get_template(element))); 
		//h_ratio->GetXaxis()->SetRangeUser(0, 3.0); 

		h_a2->Draw("E1X0");
		h_ratio->Draw("HIST P SAME"); 

		if (k==0) c1->Print(Form("./data/ACE/fill/compare_%s_flux_model.pdf(", element), "pdf"); 
		if (k>0 && k<ace->GetEntries()) c1->Print(Form("./data/ACE/fill/compare_%s_flux_model.pdf", element), "pdf"); 
		if (k==ace->GetEntries()-1) c1->Print(Form("./data/ACE/fill/compare_%s_flux_model.pdf)", element), "pdf"); 

		rescaled_fit->Write(Form("rescaled_fit_BR%d", UTimeToBR(utime))); 
		fit_ratio->Write(Form("fit_ratio_BR%d", UTimeToBR(utime))); 
		h_rig->Write(Form("h_rig_%d", UTimeToBR(utime)));
		h_rig_rescale->Write(Form("h_ratio_BR%d", UTimeToBR(utime))); 
		h_fitres[0]->Write(Form("h_fitres_ace_BR%d", UTimeToBR(utime)));		 

	}

	// gROOT->ProcessLine(".> "); 

	// gStyle->SetPalette(109); 

	for (int bin=0; bin <= h_ace->GetNbinsX(); ++bin){
 
		if (bin%2==0) {  

			g_ratio_time_norm[bin/2] = get_norm_graph( g_ratio_time[bin/2] ); 
		
			for (int iBR=0; iBR<g_ratio_time[bin/2]->GetN(); ++iBR){
				Double_t x=0, y=0; 
				g_ratio_time_norm[bin/2]->GetPoint(iBR, x, y); 
				g_ratio_time_norm[bin/2]->SetPoint(iBR, x, y-1); 			 
			}

			//c4->cd(1);
			//gPad->SetGrid(); 

			//TLegend *l_comp1 = new TLegend(0.62, 0.8, 0.9, 1.0); 
			//l_comp1->AddEntry(g_ratio_ave_norm, "Normalized Averaged ACE Data/Template-1", "PL"); 
			//l_comp1->AddEntry(g_residual_ace[bin/2], Form("ACE Fitting Residuals (%0.4f GV)", h_ace_rig_ave->GetBinCenter(bin+1)), "PL"); 

			g_residual_ace[bin/2]->SetTitle(Form(" ; ; ACE %s Data/Template-1 vs. Fitting Residuals", element));  
			g_residual_ace[bin/2]->GetXaxis()->SetTimeDisplay(1);
			g_residual_ace[bin/2]->GetXaxis()->SetTimeFormat("%m-%y");
			g_residual_ace[bin/2]->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g_residual_ace[bin/2]->GetXaxis()->SetTitleSize(0.7);
			//if (get_temp_int(element)<4) g_residual_ace[bin/2]->GetYaxis()->SetRangeUser(-0.25, 0.3);
			//g_residual_ace[bin/2]->Draw("APL");
			// g_ratio_ave_norm->Draw("PL SAME"); 
			//g_ratio_time_norm[bin/2]->Draw("PL SAME"); 
			//l_comp1->Draw("SAME"); 		
					
			g_ratio_time[bin/2]->Write(Form("g_ratio_time_%d", bin/2)); 
			g_ratio_ace[bin/2]->Write(Form("g_ratio_ace_%d", bin/2)); 
			g_ratio_ace_fit[bin/2]->Write(Form("g_ratio_ace_fit_%d", bin/2)); 
			// g_ratio_ace_fit[bin/2]->Print("range"); 
			g_residual_ace[bin/2]->Write(Form("g_residual_ace_%d", bin/2)); 
			g_resierr_ace[bin/2]->Write(Form("g_resierr_ace_%d", bin/2)); 
			g_ratio_time_norm[bin/2]->Write(Form("g_ratio_time_norm_%d", bin/2)); 
		} 
	} 
 
	file1->Close();

	TCanvas *c[7];

	TH1D *h_ene = HistTools::GraphToHist( get_ace_average_graph( element, &BRs[0], BRs.size()), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
	TH1 *h_ave = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene, isotope, "MeV/n cm", "GV m", "_rig"); 

	TFile *file2 = new TFile(Form("data/ACE/fill/%s_fill2.root", element), "RECREATE");

	for (int i=0; i<7; ++i){

		c[i] = new TCanvas("c1", "", 800, 600);
		c[i]->cd(1);

		// TLegend *legend = new TLegend(0.1,0.8,0.28,0.9);

		HistTools::SetStyle(g_ratio_time[i], HistTools::GetColorPalette(i, 7), kFullCircle, 0.8, 1, 1); 
		g_ratio_time[i]->GetXaxis()->SetTimeDisplay(1);
		g_ratio_time[i]->GetXaxis()->SetTimeFormat("%Y-%m");
		g_ratio_time[i]->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
		g_ratio_time[i]->GetXaxis()->SetTitleSize(0.7); 

		g_ratio_time[i]->SetTitle(Form(" ; ; Ratio of ACE %s Flux and %s Model (%0.4f GV)", element, get_template(element), h_ave->GetBinLowEdge(i*2))); 

		// legend->AddEntry(g_ratio_time[i], Form("%0.4f GV", h_ave->GetBinLowEdge(i*2)), "p"); 

		gPad->SetGrid();
		g_ratio_time[i]->Draw("APL");

		// SetDirectory(file1);
		// g_ratio_time[i]->Print("range"); 
		g_ratio_time[i]->Write(Form("g_ratio_time_%d", i));

		// legend->Draw("SAME"); 
		
		if (i==0) c[i]->Print(Form("./data/ACE/fill/ratio_%s_flux_model.pdf(", element), "pdf"); 
		if (i>0 && i<6) c[i]->Print(Form("./data/ACE/fill/ratio_%s_flux_model.pdf", element), "pdf"); 
		if (i==6) c[i]->Print(Form("./data/ACE/fill/ratio_%s_flux_model.pdf)", element), "pdf"); 

	}

	file2->Close();

	_file0->Close();

	return 0; 

}

void ace_rescale_BR_averaged(const char *element, Particle::Type isotope){

	gStyle->SetOptStat(0); 

 	Experiments::DataPath = "data";	
   	int data_value[n_ele] = {0, 18, 23, 27, 31, 20, 43, 22}; // p, He, Li, Be, B, C, N, O

	double *kin_bins = get_kin_bins(element);
        double *SpallCorr = get_spall_corr(element);
	double *SpallCorrUnc = get_spall_corr_unc(element);
	double *EMed = get_EMed(element);

	const UInt_t FirstACEBR = 2240;
	vector<UInt_t> BRs;
	// we stop at BR 2493, which ends on 2016/05/23, just 3 days before the end of the data taking period for AMS nuclei
	for (UInt_t br=2426; br<=2493; ++br) { 
		if (br != 2472 && br != 2473) BRs.push_back(br-FirstACEBR); 
	}

	TH1 *h_ams = Experiments::GetMeasurementHistogram(Experiments::AMS02, get_ams_data_value(element), 0); // load AMS data for a given element
	TH1 *h_ene_ave = HistTools::GraphToHist(get_ace_average_graph( element , &BRs[0], BRs.size() ), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
	TH1 *h_rig_ave = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene_ave, isotope, "MeV/n cm", "GV m", "_rig"); // load averaged ACE data for the same element, converted in rigidity
	HistTools::SetStyle(h_ams, kBlue, kFullCircle, 1.4, 1, 1);
	HistTools::SetStyle(h_rig_ave, kBlack, kFullCircle, 1.4, 1, 1);

	UShort_t namsbins = h_ams->GetNbinsX();
	UShort_t nacebins = h_rig_ave->GetNbinsX();
	double R1 = h_rig_ave->GetBinLowEdge(1);
	double R2 = h_ams->GetBinLowEdge(namsbins+1);

	int nnodes = 7; 

	TFile *file = new TFile(Form("data/ACE/compare/fit_%s_%dnodes.root", get_template(element), nnodes)); 

	gSystem->mkdir("data/ACE/fill", true);	

	TFile *_file0 = new TFile(Form("data/ACE/%s_97226_18359.root", element));
	TTree *ace=(TTree*)_file0->Get("ace");
	
	float F[7], C[7]; // must be initialized 
	Long64_t utime;
	Float_t livetime;
	
	ace->SetBranchAddress("F", F);
	ace->SetBranchAddress("C", C);
	ace->SetBranchAddress("start_utime", &utime);	
	ace->SetBranchAddress("livetime", &livetime);	

	TCanvas *c0 = new TCanvas("c0", "", 800, 600);
	c0->cd(1);
	// fill the bins for each BR 

	TH1 *h_a1 = new TH1D("", "", 200, 0, 200); 
	TH1 *h_a2 = new TH1D("", "", 200, 0, 200); 
	TH1 *h_a3 = new TH1D("", "", 300, 0, 300);
	TH1 *h_a4 = new TH1D("", "", 300, 0, 300); 

	// gROOT->ProcessLine(Form(".> data/ACE/fill/rescaling_factor_%s.dat", element)); 

	TCanvas *c1 = new TCanvas("c1", "", 1600, 900); 
	c1->Divide(1, 2);   

	TCanvas *c2 = new TCanvas("c2", "", 1600, 900);
	c2->Divide(1, 1);

	// create ACE spectral indcies vs. time 
	TGraphAsymmErrors *gspind_ace_t[5]; 
	double x_spec=0, x_sp[5], y_spec=0; // get spectral indicies 

	for (int i=0; i<5; ++i){
		gspind_ace_t[i] = new TGraphAsymmErrors(ace->GetEntries()); 
		HistTools::SetStyle(gspind_ace_t[i], kPink, kFullCircle, 0.9, 1, 1); 
	} 

	// plot BR-averaged data vs. rescaled template for each BR 
	for (int k=0; k<ace->GetEntries(); ++k){

		TH1F *h_ene = new TH1F("h_ene","", 13, kin_bins);		

		ace->GetEntry(k); //Get the nth entry of TTree!! 						

		for (int i=0; i<h_ene->GetNbinsX(); ++i) { 
			if (i%2==0) { 
				double sys_err = F[i/2] * sqrt(8e-4 + SpallCorrUnc[i/2]*SpallCorrUnc[i/2]/SpallCorr[i/2]/SpallCorr[i/2]);
				double stat_err = F[i/2]/sqrt(C[i/2]); 
				double tot_err = sqrt(stat_err*stat_err + sys_err*sys_err);
				h_ene->SetBinContent(i+1, F[i/2]); 

				if (C[i/2]!=0){
					h_ene->SetBinError(i+1, tot_err); 
				} else if (C[i/2]==0){
					continue;
				}
			} 
		}

		TH1 *h_rig = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene, isotope, "MeV/n cm", "GV m", "_rig"); // per BR
		HistTools::SetStyle(h_rig, kPink, kFullCircle, 1.4, 1, 1); 

		TGraphAsymmErrors *gspind_ace = HistTools::GetSpectralIndex(h_rig, 5, 2); // get spectral indices 
		TGraphAsymmErrors *gspind_ams = HistTools::GetSpectralIndex(h_ams, 4, 1); 

		// average of h_rig
		double rig_sum=0; // compute average of h_rig manually  
		for(int nbin=0;nbin<14;++nbin){
			rig_sum += h_rig->GetBinContent(nbin);
			//printf("ratio_sum = %0.6f \n", ratio_sum);
		}
		double rig_ave = rig_sum/7;

		// printf("rig_ave = %10.4f \n", rig_ave);  

		h_rig_ave->SetTitle(Form("%s BR-%d Rigidity Spectrum; ; ", element,  UTimeToBR(utime)));
		h_rig_ave->SetXTitle(Unit::GetEnergyLabel("GV"));
  		h_rig_ave->SetYTitle(Unit::GetDifferentialFluxLabel("GV m"));

		// rescale ACE BR to ACE Averaged Magnitude 
		Spline *sp_comb = new Spline("sp_comb", nnodes, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_comb = sp_comb->GetTF1Pointer();  
		TF1 *fit_comb = (TF1*) file->Get("fit_both"); 

		HistTools::CopyParameters(fit_comb, fsp_comb); // error 
		double x1, x2;
		fit_comb->GetRange(x1,x2);
		fsp_comb->SetRange(x1,x2);

		// Modify Nodes of The Spline 
		//NormSpline *ns = new NormSpline("f_norm_spline", R1, R2, fsp_comb); // normalized spline with first node and derivative free
		//TF1 *f_norm_spline = ns->GetTF1Pointer();
		//ns->Print(); 
		// set initial parameters: average value of current BR ACE flux, same first derivative and first Y node as template
		//f_norm_spline->SetParameters(rig_ave, fsp_comb->GetParameter(0), fsp_comb->GetParameter(8)); // replace 8 with appropriate index of first Y node
		//h_rig->Fit(f_norm_spline, "NQ"); // fit current BR ACE flux
		//f_norm_spline->Print(); 

		double *nodes = sp_comb->GetXNodes(); 

		TLine *l1 = new TLine( nodes[0], 0, nodes[0], h_rig->GetMaximum()*3); 
		TLine *l2 = new TLine( nodes[1], 0, nodes[1], h_rig->GetMaximum()*3); 

		TLegend *legend = new TLegend(0.62,0.8,0.9,1.0);
		TLegend *legend2 = new TLegend(0.62,0.8,0.9,0.9);  

		l1->SetLineColor(kGreen-3);
		l2->SetLineColor(kGreen-3);
		l1->SetLineStyle(2);
		l2->SetLineStyle(2); 		

		TH1 *h_ratio = (TH1 *) h_rig->Clone("h_ratio");
		TH1 *h_ratio2 = (TH1 *) h_ams->Clone("h_ratio2"); 		

		h_ratio->Divide(fsp_comb); 
		h_ratio2->Divide(fsp_comb); 

		HistTools::SetStyle(h_ratio, kPink, kFullCircle, 1.4, 1, 1);
		HistTools::SetStyle(h_ratio2, kBlue, kFullCircle, 1.4, 1, 1);

		double ratio_sum=0; // compute average of h_ratio manually  
		for(int nbin=0;nbin<14;++nbin){
			ratio_sum += h_ratio->GetBinContent(nbin);
			//printf("ratio_sum = %0.6f \n", ratio_sum);
		}
		double ratio_ave = ratio_sum/7;

		// printf("%0.6f \n", ratio_ave);
		//HistTools::PrintFunction(fit_comb);
			
		double scale = 1./ratio_ave;
		//if (get_index(element)<4) h_ratio->Scale(scale);	

		TF1 *rescaled_fit = HistTools::CombineTF1Const(fsp_comb, ratio_ave, HistTools::MultiplyConst, "rescaled_fit", R1, R2); 

		// rescaled_fit->Print(); 
		
		c1->cd(1);
		gPad->SetGrid(); 
		gPad->SetLogy();	

		h_a1->SetTitle(Form("%s BR-%d Model Rigidity Spectrum; ; ", element, UTimeToBR(utime)));
		h_a1->SetXTitle(Unit::GetEnergyLabel("GV"));
  		h_a1->SetYTitle(Unit::GetDifferentialFluxLabel("GV m"));

		//if (get_temp_int(element)<4) h_a1->GetYaxis()->SetRangeUser(0.1, h_rig_ave->GetMaximum()*3); 
		h_a1->GetXaxis()->SetRangeUser(0., 10.); 	

		//h_rig_ave->GetYaxis()->SetRangeUser(0.0001, 2.0); 
		h_rig_ave->GetXaxis()->SetRangeUser(0., 10.);

		legend->AddEntry(rescaled_fit, "rescaled template", "l"); 
		legend->AddEntry(h_rig_ave, "ACE Mean flux", "p"); 
		legend->AddEntry(h_rig, "ACE BR Flux", "p"); 
		legend->AddEntry(h_ams, "AMS Integrated flux", "p"); 
		legend->SetTextSize(0.05); 

		legend2->AddEntry(h_ratio, "ratio of ACE BR Flux wrt the template ", "p");
		// legend2->AddEntry(h_ratio2, "ratio wrt Integrated AMS Flux", "p"); 
	
		h_a1->Draw("E1X0");  
		rescaled_fit->SetRange(0., 10.); 
		rescaled_fit->Draw("SAME"); 
		h_rig_ave->Draw("E1X0 SAME"); 
		h_rig->Draw("E1X0 SAME"); 
		h_ams->Draw("E1X0 SAME"); 
		legend->Draw("SAME"); 
		l1->Draw();
		l2->Draw();   

		c1->cd(2); 
		gPad->SetGrid(); 

		//if (get_temp_int(element)<4) h_a2->GetYaxis()->SetRangeUser(0.5, 1.5); 
		h_a2->GetXaxis()->SetRangeUser(0., 10.); 

		h_a2->SetTitle("");
		h_a2->SetXTitle(Unit::GetEnergyLabel("GV"));
		h_a2->SetYTitle(Form("%s BR-Averaged Flux / %s Model", element, get_template(element))); 
		//h_ratio->GetXaxis()->SetRangeUser(0, 3.0); 

		h_a2->Draw("E1X0"); 
		h_ratio->Draw("E1X0 SAME"); 
		// h_ratio2->Draw("E1X0 SAME"); 
		legend2->Draw("SAME"); 

		if (k==0) c1->Print(Form("./data/ACE/fill/compare2_%s_flux_model.pdf(", element), "pdf"); 
		if (k>0 && k<ace->GetEntries()) c1->Print(Form("./data/ACE/fill/compare2_%s_flux_model.pdf", element), "pdf"); 
		if (k==ace->GetEntries()-1) c1->Print(Form("./data/ACE/fill/compare2_%s_flux_model.pdf)", element), "pdf"); 

		c2->cd(1);
		gPad->SetGrid();
	
		TLegend *l_both = new TLegend(0.62,0.8,0.9,0.9); 
		l_both->AddEntry(gspind_ace, "ACE Spectral Indices", "PL");
		l_both->AddEntry(gspind_ams, "AMS Spectral Indices", "PL"); 		
	
		//if (get_temp_int(element)<4) gspind_ace->GetYaxis()->SetRangeUser(-4, 7); 
		TAxis *axis1 = gspind_ace->GetXaxis(); 
		axis1->SetLimits(0., 10.); 

		//h_a3->Draw("E1X0"); 
		HistTools::SetStyle(gspind_ace, kPink, kFullCircle, 1.4, 1, 1); 
		gspind_ace->SetTitle(Form("; ; ACE %s Spectral Indices BR-%d", element, UTimeToBR(utime)));
		gspind_ace->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV"));
		gspind_ace->Draw("APL"); 
		
		// h_a3->Draw("E1X0"); 

		//h_a4->GetYaxis()->SetRangeUser(-5, 0); 
		TAxis *axis2 = gspind_ams->GetXaxis();

		axis2->SetLimits(0., 10.); 

		//h_a4->Draw("E1X0"); 
		HistTools::SetStyle(gspind_ams, kBlue, kFullCircle, 1.4, 1, 1); 
		gspind_ams->SetTitle(Form("; ; AMS %s Spectral Indices", element)); 
		gspind_ams->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV"));
		gspind_ams->Draw("PL SAME"); 
		l_both->Draw("SAME"); 

		if (k==0) c2->Print(Form("./data/ACE/fill/spind_%s_flux_model.pdf(", element), "pdf"); 
		if (k>0 && k<ace->GetEntries()) c2->Print(Form("./data/ACE/fill/spind_%s_flux_model.pdf", element), "pdf"); 
		if (k==ace->GetEntries()-1){
			// gspind_ace->Print();
			// cout << " " << endl; 
			// gspind_ams->Print();  
			c2->Print(Form("./data/ACE/fill/spind_%s_flux_model.pdf)", element), "pdf"); 
		}

		for (int i=0; i<5; ++i){ 
			gspind_ace->GetPoint(i, x_spec, y_spec);
			x_sp[i] = x_spec; 
			// printf("x_spec = %10.4f, y_spec=%10.4f \n", x_spec, y_spec);  
			gspind_ace_t[i]->SetPoint(k, utime, y_spec); 
			gspind_ace_t[i]->SetPointError(k, 0, 0, gspind_ace->GetErrorYlow(i), gspind_ace->GetErrorYhigh(i)); 
		}
	} // end of BR loop 

	TCanvas *c3 = new TCanvas("", "", 1800, 600); 

	for (int i=0; i<5; ++i){

		c3->cd(1); 
		gPad->SetGrid(); 

		gspind_ace_t[i]->SetTitle(Form(" ; ; ACE %s Spectral Indices at R = %7.3f GV", element, x_sp[i])); 
		gspind_ace_t[i]->GetXaxis()->SetTimeDisplay(1);
		gspind_ace_t[i]->GetXaxis()->SetTimeFormat("%Y-%m");
		gspind_ace_t[i]->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
		gspind_ace_t[i]->GetXaxis()->SetTitleSize(0.7); 	

	 	//if (get_temp_int(element)<4) gspind_ace_t[i]->GetYaxis()->SetRangeUser(-4, 7); 
		gspind_ace_t[i]->Draw("APL"); 

		//gspind_ace_t[i]->Print(); 
		//cout << " " << endl; 
		if (i==0) c3->Print(Form("./data/ACE/fill/spind(t)_%s_flux_model.pdf(", element), "pdf");
		if (i>0 && i<4) c3->Print(Form("./data/ACE/fill/spind(t)_%s_flux_model.pdf", element), "pdf");
		if (i==4) c3->Print(Form("./data/ACE/fill/spind(t)_%s_flux_model.pdf)", element), "pdf"); 
	} 

	// gROOT->ProcessLine(".> "); 

}

void ace_fake_td_ams(const char *element, Particle::Type isotope){

	gStyle->SetOptStat(0); 
	// gStyle->SetTitleSize(26,"t");

	gSystem->mkdir("data/ACE/fill/fake_td_ams", true);	

 	Experiments::DataPath = "data";	
   	int data_value[n_ele] = {0, 18, 23, 27, 31, 20, 43, 22}; // p, He, Li, Be, B, C, N, O

	double *kin_bins = get_kin_bins(element);
        double *SpallCorr = get_spall_corr(element);
	double *SpallCorrUnc = get_spall_corr_unc(element);
	double *EMed = get_EMed(element);

	const UInt_t FirstACEBR = 2240;
	vector<UInt_t> BRs;
	// we stop at BR 2493, which ends on 2016/05/23, just 3 days before the end of the data taking period for AMS nuclei
	for (UInt_t br=2426; br<=2493; ++br) { 
		if (br != 2472 && br != 2473) BRs.push_back(br-FirstACEBR); 
	}

	// load AMS He BR fluxes 
	const int nBRs = Experiments::Info[Experiments::AMS02].Dataset[1].nMeasurements; // read the number of BRs
	TH1D **h_BR_he = Experiments::GetDatasetHistograms(Experiments::AMS02, 4);
	// Create a new set of histograms which match the AMS monthly bins with the AMS integrated bins 

	TH1 *h_ams = Experiments::GetMeasurementHistogram(Experiments::AMS02, get_ams_data_value(element), 0); // load AMS data for a given template 
	TH1 *h_ams_new = (TH1D*) h_BR_he[0]->Clone("h_ams_new"); 
	// TH1 *h_he_int = (TH1*) ave_hist( h_BR_he, nBRs);  
	TH1 *h_he = Experiments::GetMeasurementHistogram(Experiments::AMS02, 18, 0);
	TH1 *h_he_int = (TH1D*) h_BR_he[0]->Clone("h_ams_new"); 

	// create h_ams_new and h_he_int to match the bins  

	for (int bin=1; bin <= h_ams_new->GetNbinsX(); ++bin){ 
	   if (!strcmp(element, "B") || !strcmp(element, "C")){ 
		h_ams_new->SetBinContent(bin, h_ams->GetBinContent(bin)); 
		h_ams_new->SetBinError(bin, h_ams->GetBinError(bin));
	   }
	   if (!strcmp(element, "N") || !strcmp(element, "O")){ 
		h_ams_new->SetBinContent(bin, h_ams->GetBinContent(bin-1)); 
		h_ams_new->SetBinError(bin, h_ams->GetBinError(bin-1)); 
	   }else{ 
		h_ams_new->SetBinContent(bin, h_ams->GetBinContent(bin)); 
		h_ams_new->SetBinError(bin, h_ams->GetBinError(bin)); 
	   }
	}

	if (!strcmp(element, "N") || !strcmp(element, "O")){
		h_ams_new->SetBinContent(1, 0); 
		h_ams_new->SetBinError(1, 0);
	}
 
	for (int bin=1; bin<=h_BR_he[0]->GetNbinsX(); ++bin) { 
		h_he_int->SetBinContent(bin, h_he->GetBinContent(bin)); 
		h_he_int->SetBinError(bin, h_he->GetBinError(bin)); 
	}

	//return 0; 

	// h_ams->Print("range"); 
	// h_ams_new->Print("range"); 

	// load ACE average 
	TH1 *h_ace_ene_ave = HistTools::GraphToHist(get_ace_average_graph( element , &BRs[0], BRs.size() ), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
	TH1 *h_ace_rig_ave = HistTools::TransformEnergyAndDifferentialFluxNew(h_ace_ene_ave, isotope, "MeV/n cm", "GV m", "_rig"); // load averaged ACE data for the same element, converted in rigidity
	HistTools::SetStyle(h_ams_new, kBlue, kFullCircle, 1.4, 1, 1);
	HistTools::SetStyle(h_ace_rig_ave, kBlack, kFullCircle, 1.4, 1, 1);

	UShort_t namsbins = h_ams_new->GetNbinsX(); 
	UShort_t nacebins = h_ace_rig_ave->GetNbinsX(); 
	double R1 = h_ace_rig_ave->GetBinLowEdge(1);
	double R2 = h_ams->GetBinLowEdge(namsbins+1);

	int nnodes = 9; // combined template 
	int nnodes_ams = 7; 

	printf(" template = %s, %d \n", element, get_ams_data_value(element));  

	// load combined fit template 
	TFile *file0 = new TFile(Form("data/ACE/compare/fit_%s_%dnodes.root", get_template(element), nnodes));  

	Spline *sp_comb = new Spline("sp_comb", nnodes, Spline::LogLog | Spline::PowerLaw);
	TF1 *fsp_comb = sp_comb->GetTF1Pointer();  
	TF1 *fit_comb = (TF1*) file0->Get("fit_both"); 

	HistTools::CopyParameters(fit_comb, fsp_comb); // error 
	double x1, x2;
	fit_comb->GetRange(x1,x2);
	fsp_comb->SetRange(x1,x2);

	// load AMS He combined, modified spline  
	TFile *file1 = new TFile(Form("data/amsfit/fit_result_node%d.root", nnodes_ams));

	// load ACE BR 
	TFile *file2 = new TFile(Form("data/ACE/fill/%s_fill.root", element));
	TFile *file2a = new TFile(Form("data/ACE/fill/%s_fill2.root", element));  

	TGraphErrors *g_ratio_time[7]; 
	TGraphErrors *g_ratio_time_norm[7]; 

	for (int i=0; i<7; ++i){
		g_ratio_time[i] = (TGraphErrors*) file2a->Get(Form("g_ratio_time_%d", i)); 	
		HistTools::SetStyle(g_ratio_time[i], kBlue, kFullCircle, 1.4, 1, 1);
	}

	// load AMS He Model by Rescaled C Template
	TFile *file3 = new TFile(Form("data/ACE/extend2/fit_C_temp_he_%dnodes.root", nnodes));

	TCanvas *c0 = new TCanvas("c0", "", 1600, 900); // fluxes & flux/template
	c0->Divide(1, 2); 

	TCanvas *c1 = new TCanvas("c1", "", 1600, 900); // He ratio, F(R,t)/<F(R)> & spectral indices 
	c1->Divide(1, 2); 
	
	TCanvas *c2 = new TCanvas("c2", "", 1600, 900); // Modeling F(R,t)/<F(R)> 
							// fit/data of F(R,t)/<F(R)> 
	c2->Divide(1, 3); 				// spectral indices of F(R,t)/<F(R)>, here is defined by h_ratio (ACE) and h_ratio0 (AMS)

	TCanvas *c3 = new TCanvas("c3", "", 1600, 900); // plot chi2/ndf vs. BR for both ACE & AMS F(R,t)/<F(R)> fits  
	c3->Divide(2, 1);  

	TCanvas *c4 = new TCanvas("c4", "", 1600, 900); // plot for F(R,t)/<F(R)> fit residuals for each ACE bin 
	c4->Divide(1, 1); 

	TCanvas *c4_2 = new TCanvas("c4_2", "", 1600, 900); // c4 but for each AMS bin 
	c4_2->Divide(1, 1); 
	
	TCanvas *c5 = new TCanvas("c5", "", 1600, 900); // plot parameter [0] [1] vs. BR for the fit 
	c5->Divide(1, 3); 

	TH1 *h_a1 = new TH1D("", "", 3000, 0, 3000); 
	TH1 *h_a2 = new TH1D("", "", 3000, 0, 3000);
	TH1 *h_a3 = new TH1D("", "", 3000, 0, 3000);
	TH1 *h_a4 = new TH1D("", "", 3000, 0, 3000);	
	TH1 *h_a5 = new TH1D("", "", 3000, 0, 3000); 

	vector<double> last_pars;
	vector<double> last_pars2 = {0.5, 0.5};  

	TGraph *g_chi2_ace = new TGraph(); 
	TGraph *g_chi2_ams = new TGraph(); 
	TGraph *g_chi2_ace2 = new TGraph(); 
	TGraph *g_chi2_ams2 = new TGraph(); 
	TGraphErrors *g_pars[2];

	for (int ipar=0; ipar<2; ++ipar){
		g_pars[ipar] = new TGraphErrors(); 
	} 

	TGraphErrors *g_ratio_ace[7];	
	TGraphErrors *g_ratio_ace_fit[7]; 
	TGraphErrors *g_ratio_ams[h_ams_new->GetNbinsX()]; 
	TGraphErrors *g_ratio_ams_fit[h_ams_new->GetNbinsX()]; 

	TGraphErrors *g_ratio_ave = new TGraphErrors();
	TGraphErrors *g_ratio_ave_norm; // normalized averaged ratio of data/template-1 vs. time  

	TGraphErrors *g_residual_ace[7];
	TGraphErrors *g_resierr_ace[7]; 
	TGraphErrors *g_residual_norm_ace[7];
	TGraphErrors *g_residual_ams[h_ams_new->GetNbinsX()]; 
	TGraphErrors *g_residual_norm_ams[h_ams_new->GetNbinsX()];

	for (int bin=0; bin <= h_ace_rig_ave->GetNbinsX(); ++bin){ 
		if (bin%2==0){  
			g_ratio_ace[bin/2] = new TGraphErrors(nBRs); 
			g_ratio_ace_fit[bin/2] = new TGraphErrors(nBRs); 
			g_residual_ace[bin/2] = new TGraphErrors(nBRs);  
			g_resierr_ace[bin/2] = new TGraphErrors(nBRs); 
			HistTools::SetStyle(g_residual_ace[bin/2], kPink, kFullCircle, 1.4, 1, 1);
			HistTools::SetStyle(g_ratio_ave, kBlue+1, kFullCircle, 1.4, 1, 1);
		} 
	}  
 
	for (int bin=0; bin <= h_ams_new->GetNbinsX(); ++bin){ 
		g_ratio_ams[bin] = new TGraphErrors(nBRs); 
		g_ratio_ams_fit[bin] = new TGraphErrors(nBRs); 
		g_residual_ams[bin] = new TGraphErrors(nBRs);  
		HistTools::SetStyle(g_residual_ams[bin], kPink, kFullCircle, 1.4, 1, 1); 
	}   

	TFile *file = new TFile(Form("data/ACE/fill/F2_%s.root", element), "RECREATE");
 
	//return 0; 

	int iBR_true = 0;
	for (int iBR=0; iBR<nBRs; ++iBR){

		// AMS_He Monthly Spline 
		Spline *sp_he = new Spline("sp_he", nnodes_ams, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_he = sp_he->GetTF1Pointer();  
		TF1 *fit_he = (TF1*) file1->Get(Form("fsp_BR_he_%02d", iBR)); 

		HistTools::CopyParameters(fit_he, fsp_he); // error 
		fit_he->GetRange(x1,x2);
		fsp_he->SetRange(x1,x2);

		// AMS_He Integrated Spline  
		Spline *sp_he_ave = new Spline("sp_he_ave", nnodes_ams, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_he_ave = sp_he_ave->GetTF1Pointer();  
		TF1 *fit_he_ave = (TF1*) file1->Get("fsp_he"); 

		HistTools::CopyParameters(fit_he_ave, fsp_he_ave); // error 
		x1=0, x2=0;
		fit_he_ave->GetRange(x1,x2);
		fsp_he_ave->SetRange(x1,x2); 

		// AMS_He Model by Rescaled C ACE+AMS Combined Template  
		Spline *sp_he_temp = new Spline("sp_he_temp", nnodes, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_he_temp = sp_he_temp->GetTF1Pointer();  
		TF1 *fit_he_temp = (TF1*) file3->Get("fsp_he"); 

		HistTools::CopyParameters(fit_he_temp, fsp_he_temp); // error 
		x1=0, x2=0; 
		fit_he_temp->GetRange(x1,x2);
		fsp_he_temp->SetRange(x1,x2);  

		TF1 *fsp_he2 = HistTools::CombineTF1(fsp_he, fsp_he_ave, HistTools::Divide, "fsp_he2"); // He(R,t) Fit vs. <He(R)> Fit
		TF1 *fspind_he = HistTools::GetSpectralIndex(fsp_he2, "fspind_he", R1, R2); 
		fspind_he->SetRange(h_he_int->GetBinCenter(2), x2);

		// AMS_He(R,t)/<He(R)>
		TH1 *h_ratio0 = (TH1D*) h_BR_he[iBR]->Clone("h_ratio0"); 
		h_ratio0->Divide(h_he_int); 

		// refill AMS BR He to match the bins  
		HistTools::SetStyle(h_BR_he[iBR], kBlue, kFullCircle, 1.1, 1, 1); 

		TH1 *h_ace_BR = (TH1*) file2->Get(Form("h_rig_%s_BR%d", element, 2426+iBR_true)); 

		// rescaled <AMS_C> by AMS He BR / <AMS_He>   
		TH1 *h_ams_BR_fake = (TH1D*) h_ams_new->Clone("h_ams_BR_fake");

		h_ams_BR_fake->Multiply(h_BR_he[iBR]); 
		h_ams_BR_fake->Divide(h_he_int); 
 
		// h_ams_BR_fake->Print("range"); 

/*
		fsp_comb->Draw();
		h_ace_rig_ave->Draw("E1X0 SAME");
		h_ams_new->Draw("E1X0 SAME");  

		fsp_he_ave->Draw("SAME");
		fsp_he_ave->SetLineColor(kGreen);		
		HistTools::SetStyle(h_ace_rig_ave, kBlack, kFullCircle, 0.9, 1, 1);
		HistTools::SetStyle(h_ams_new, kBlack, kFullCircle, 0.9, 1, 1);

		break; 
*/  

		// rescaled combined fit 
		TH1 *h_ratio;
		if (get_index(element) <4) h_ratio = (TH1 *) h_ace_BR->Clone("h_ratio"); // ACE BR/Temp 	
		if (get_index(element)>=4) h_ratio = (TH1 *) h_ace_rig_ave->Clone("h_ratio"); // ACE Ave/Temp 	
		h_ratio->Divide(fsp_comb); 

		HistTools::SetStyle(h_ratio, kPink, kFullCircle, 1.4, 1, 1);

		double ratio_sum=0; // compute average of h_ratio manually  
		for(int nbin=0;nbin<14;++nbin){
			ratio_sum += h_ratio->GetBinContent(nbin);
			//printf("ratio_sum = %0.6f \n", ratio_sum);
		}
		double ratio_ave = ratio_sum/7;

		// printf("%0.6f \n", ratio_ave);
		//HistTools::PrintFunction(fsp_comb);	
			
		double scale = 1./ratio_ave;
		//if (get_index(element) <4) h_ratio->Scale(scale);	

		// rescaled combined fit 
		TF1 *rescaled_fit = HistTools::CombineTF1Const(fsp_comb, ratio_ave, HistTools::MultiplyConst, "rescaled_fit", R1, R2); 
		// rescaled_fit = HistTools::CombineTF1(rescaled_fit, fsp_comb, HistTools::Divide, "fit_ratio", R1, R2);  

		TH1 *h_ratio1 = (TH1D*) h_ace_BR->Clone("h_ratio1"); // spind model  
		if (get_index(element) <4) h_ratio1->Divide(fsp_comb); 
		if (get_index(element)>=4) h_ratio1->Divide(rescaled_fit); 

		TH1 *h_ratio2 = (TH1D*) h_ams_BR_fake->Clone("h_ratio2"); // spind model 
		h_ratio2->Divide(h_ams_new);
   
		// h_ratio2->Print("range");

		// estimated AMS / resclaed combined fit
		TH1 *h_ratio3 = (TH1D *) h_ams_BR_fake->Clone("h_ratio3");		

		h_ratio3->Divide(rescaled_fit); 
		HistTools::SetStyle(h_ratio3, kPink, kFullCircle, 1.4, 1, 1);

		c0->cd(1);
		gPad->SetGrid(); 
		gPad->SetLogy();
		gPad->SetLogx(); 	

		TLegend *legend = new TLegend(0.62,0.75,0.9,0.9);
		legend->AddEntry(h_ams, Form("Integrated AMS %s Flux", element), "p"); 
		legend->AddEntry(h_ace_BR, Form("ACE %s BR Flux", element), "p"); 
		legend->AddEntry(h_ams_BR_fake, Form("Estimated AMS %s BR Flux", element), "p"); 
		legend->AddEntry(rescaled_fit, "rescaled combined template", "l"); 

		h_a1->SetTitle(Form("%s BR-%d Model Rigidity Spectrum; ; ", element, 2426+iBR_true));
		h_a1->SetXTitle(Unit::GetEnergyLabel("GV"));
  		h_a1->SetYTitle(Unit::GetDifferentialFluxLabel("GV m"));

		//if (get_temp_int(element)<4) h_a1->GetYaxis()->SetRangeUser(1e-3, 1e1); 
		// h_a1->GetXaxis()->SetRangeUser(0.01, 3000.); 

		TAxis *axis1 = h_a1->GetXaxis(); 
		axis1->SetLimits(0.7, 60.); 
	
		h_a1->Draw("E1X0"); 
	
		HistTools::SetStyle(h_ams, kPink-3, kFullCircle, 1.4, 1, 1); 
		HistTools::SetStyle(h_ace_BR, kRed, kFullCircle, 1.4, 1, 1);
		// HistTools::SetStyle(h_ams_BR_he, kBlue, kFullCircle, 1.4, 1, 1);
		HistTools::SetStyle(h_ams_BR_fake, kBlue, kFullCircle, 1.4, 1, 1);

		// fsp_he->Draw("SAME"); 
		rescaled_fit->Draw("SAME"); 
		h_ams->Draw("E1X0 SAME"); 
		h_ace_BR->Draw("E1X0 SAME"); 
		// h_ams_BR_he->Draw("E1X0 SAME"); 
		h_ams_BR_fake->Draw("E1X0 SAME"); 

		//fit->SetLineColor(kBlue); 
		//fit->Draw("SAME"); 
		legend->Draw("SAME");   

		c0->cd(2); 
		gPad->SetGrid(); 
		gPad->SetLogx(); 

		TAxis *axis2 = h_a2->GetXaxis(); 
		axis2->SetLimits(0.7, 60.); 

		//if (get_temp_int(element)<4) h_a2->GetYaxis()->SetRangeUser(0.5, 1.5); 
		h_a2->GetXaxis()->SetRangeUser(0.7, 60.); 

		h_a2->SetTitle("");
		h_a2->SetXTitle(Unit::GetEnergyLabel("GV"));
		h_a2->SetYTitle(Form("Estimated BR Fluxes / Rescaled %s Model", element)); 
		h_a2->SetTitleSize(0.3,"t"); 

		h_a2->Draw("E1X0"); 
		h_ratio->Draw("E1X0 SAME"); 
		h_ratio3->Draw("E1X0 SAME"); 

		if (iBR==0) { c0->Print(Form("./data/ACE/fill/fake_td_ams/%s_flux_model.pdf(", element), "pdf"); } 
		if (iBR>0 && iBR<nBRs-1) { c0->Print(Form("./data/ACE/fill/fake_td_ams/%s_flux_model.pdf", element), "pdf"); } 
		if (iBR==nBRs-1){
			c0->Print(Form("./data/ACE/fill/fake_td_ams/%s_flux_model.pdf)", element), "pdf"); 
		} 

/*
		// plot AMS Integrated / Rescaled Combined Fit Template
		c1->cd(1); 
		gPad->SetGrid(); 
		gPad->SetLogx(); 

		TH1 *h_ratio3 = (TH1*) h_ams_new->Clone("h_ams_new"); 
		h_ratio3->Divide(fsp_comb); 

		TAxis *axis3 = h_ratio3->GetXaxis(); 
		axis3->SetLimits(0.7, 60.);

		HistTools::SetStyle(h_ratio3, kBlack, kFullCircle, 1.4, 1, 1);
		h_ratio3->SetXTitle(Unit::GetEnergyLabel("GV"));
		h_ratio3->SetTitle(Form("Template Test; ; AMS Integrated %s Flux / %s Template", element, element)); 
		h_ratio3->Draw("E1X0"); 
*/

		// plot the ratio of AMS_He_BR / <AMS_He>
		c1->cd(1); 
		gPad->SetGrid(); 
		gPad->SetLogx(); 

		//if (get_temp_int(element)<4) h_a3->GetYaxis()->SetRangeUser(0.5, 1.5); 
		h_a3->GetXaxis()->SetLimits(0.7, 60);		

		TLegend *l1 = new TLegend(0.62,0.8,0.9,0.9); 
		l1->AddEntry(h_ratio0, "AMS He(R,t)/<He(R)>", "PL"); 	

		h_a3->Draw("E1X0"); 
		h_a3->SetXTitle(Unit::GetEnergyLabel("GV"));
		h_a3->SetTitle(Form("; ; BR-%d He(R,t)/<He(R)>", iBR_true+2426)); 
		HistTools::SetStyle(h_ratio0, kBlue, kFullCircle, 1.4, 1, 1);
		h_ratio0->Draw("E1X0 SAME"); 

		// plot spectral indices
		TGraphAsymmErrors *gspind_ace = HistTools::GetSpectralIndex(h_ace_BR, 5, 2); // get spectral indices 
		TGraphAsymmErrors *gspind_ams = HistTools::GetSpectralIndex(h_ams_BR_fake, 4, 1); 

		c1->cd(2); 
		gPad->SetGrid(); 
		gPad->SetLogx();  

		TLegend *l_both = new TLegend(0.62,0.8,0.9,1.); 
		l_both->AddEntry(gspind_ace, "ACE Spectral Indices", "PL");
		l_both->AddEntry(gspind_ams, "AMS Spectral Indices", "PL"); 		
	
		//if (get_temp_int(element)<4) gspind_ace->GetYaxis()->SetRangeUser(-3, 3); 
		gspind_ace->GetXaxis()->SetLimits(0.7, 60);

		HistTools::SetStyle(gspind_ace, kPink, kFullCircle, 1.4, 1, 1); 
		gspind_ace->SetTitle(Form("; ; ACE %s Spectral Indices BR-%d", element, 2426+iBR_true));
		gspind_ace->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV"));
		gspind_ace->Draw("APL"); 

		//h_a4->Draw("E1X0"); 
		HistTools::SetStyle(gspind_ams, kBlue, kFullCircle, 1.4, 1, 1); 
		gspind_ams->SetTitle(Form("; ; AMS %s Spectral Indices", element)); 
		gspind_ams->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV"));
		gspind_ams->GetXaxis()->SetTitleSize(1.5); 
		gspind_ams->Draw("PL SAME"); 
		// fspind->Draw("SAME"); 
		l_both->Draw("SAME"); 

		if (iBR==0) { c1->Print(Form("./data/ACE/fill/fake_td_ams/fake_ams_spind_%s.pdf(", element), "pdf"); }
		if (iBR>0 && iBR<nBRs-1) { c1->Print(Form("./data/ACE/fill/fake_td_ams/fake_ams_spind_%s.pdf", element), "pdf"); }
		if (iBR==nBRs-1){
			// gspind_ace->Print("range");
			// cout << " " << endl; 
			// gspind_ams->Print("range"); 
			c1->Print(Form("./data/ACE/fill/fake_td_ams/fake_ams_spind_%s.pdf)", element), "pdf"); 
		} 

		// find spectral indices of h_ratio + h_ratio0 
		TGraphAsymmErrors *gspind_ace2 = HistTools::GetSpectralIndex(h_ratio1, 5, 2); // get spectral indices 
		TGraphAsymmErrors *gspind_ams2 = HistTools::GetSpectralIndex(h_ratio2, 4, 1); 

		HistTools::SetStyle(gspind_ace2, kPink, kFullCircle, 1.4, 1, 1); 
		HistTools::SetStyle(gspind_ams2, kBlue, kFullCircle, 1.4, 1, 1); 
		HistTools::SetStyle(fspind_he, kBlack, kFullCircle, 0.8, 1, 1);

		// fit to F(R,t)/<F(R)>
		TF1 *fit_ratio = new TF1("f_ratio", "1+[0]*exp(-[1]*x)", 0.7, 60);
		//TF1 *fit_ratio2 = new TF1("f_ratio2", "1+[0]*exp(-[1]*log(x))", 0.7, 60);

		double F1 = h_ratio1->GetBinContent(h_ratio1->FindBin(R1)), F2 = h_ratio0->GetBinContent(h_ratio0->GetNbinsX()); 
		double A1 = log(abs(F1-1)), A2 = log(abs(F2-1)); 
		TF1 *f1 = new TF1("f1","abs(sin(x)/x)*sqrt(x)",0,1);
		double rho = f1->GetRandom();
		double N = f1->GetRandom(); 
		
		// double rho = (R2-R1)/(A1-A2)
		// double N = 0.5*(abs(F1-1)/exp(-R1/rho) + abs(F2-1)/exp(-R2/rho)); 

		last_pars = {N, rho}; 
		last_pars2 = {N, rho};
 
		printf(" rho = %0.4f, N = %0.20f \n ", rho, N); 
		// printf(" rho = %0.4f, N = %0.20f, F1 = %0.4f, F2 = %0.4f, (F1-1)/exp(-R1/rho) = %0.6f, (F2-1)/exp(-R2/rho) = %0.6f \n", rho, N, F1, F2, (F1-1)/exp(-R1/rho), (F2-1)/exp(-R2/rho));

		for (int ipar=0; ipar<fit_ratio->GetNpar(); ++ipar){ 
			fit_ratio->SetParameter(ipar, last_pars[ipar]); 
			//fit_ratio2->SetParameter(ipar, last_pars2[ipar]); 
		}   

		printf("par0 = %0.4f, par1 = %0.4f \n", fit_ratio->GetParameter(0), fit_ratio->GetParameter(1)); 

		//h_ratio1->Print("range");
		//h_ratio0->Print("range");  

		TObjArray data; 
		data.Add(h_ratio1);
		data.Add(h_ratio0); 

		// fit AMS and ACE data at the same time
		vector<double> rigmin, rigmax, chi2norm(2);
		ROOT::Fit::Fitter fitter;
		//HistTools::PrintFunction(fit);
		FitTools::SetCommonFitterOptions(fitter);
		FitTools::FitCombinedData(data, fit_ratio, "I", rigmin, rigmax, chi2norm, fitter, 3); 

		//vector<double> rigmin2, rigmax2, chi2norm2(2);
		//ROOT::Fit::Fitter fitter2;
		//FitTools::SetCommonFitterOptions(fitter2);
		//FitTools::FitCombinedData(data, fit_ratio2, "I", rigmin2, rigmax2, chi2norm2, fitter2, 3); 
 
		for (int ipar=0; ipar<fit_ratio->GetNpar(); ++ipar){
			// last_pars[ipar] = fit_ratio->GetParameter(ipar); 
			// last_pars2[ipar] = fit_ratio2->GetParameter(ipar); 
			if (abs(fit_ratio->GetParError(ipar)/fit_ratio->GetParameter(ipar))<1){ 
				g_pars[ipar]->SetPoint(iBR, UBRToTime(iBR_true+2426), fit_ratio->GetParameter(ipar)); 
				g_pars[ipar]->SetPointError(iBR, 0, fit_ratio->GetParError(ipar)); 
			} else {
				g_pars[ipar]->SetPoint(iBR, UBRToTime(iBR_true+2426), 0); 
				g_pars[ipar]->SetPointError(iBR, 0, 0); 
			}
		}  

		TF1 *fspind_ratio = HistTools::GetSpectralIndex(fit_ratio, "fspind_ratio", R1, R2); 
		//TF1 *fspind_ratio2 = HistTools::GetSpectralIndex(fit_ratio2, "fspind_ratio2", R1, R2); 

		//h_ams_BR_fake->Print("range"); 
		c2->cd(1); 
		gPad->SetGrid(); 
		gPad->SetLogx();
		gPad->SetBottomMargin(0.01); 

		//if (get_temp_int(element)<4) h_a4->GetYaxis()->SetRangeUser(0.5, 2.2); 
		h_a4->GetXaxis()->SetLimits(0.7, 60); 
	
		TLegend *l2 = new TLegend(0.62,0.8,0.9,1.); 
		l2->AddEntry(h_ratio1, Form("ACE %s(R,t)/<%s(R)>", element, element), "PL");
		l2->AddEntry(h_ratio0, "AMS He(R,t)/<He(R)>", "PL");  
		l2->AddEntry(fit_ratio, Form("1+%6.3fexp(-%6.3fR)", fit_ratio->GetParameter(0), fit_ratio->GetParameter(1)), "L"); 
		//l2->AddEntry(fit_ratio2, Form("1+%6.3fexp(-%6.3fln(R))", fit_ratio2->GetParameter(0), fit_ratio2->GetParameter(1)), "L"); 

		h_a4->Draw("E1X0"); 
		h_a4->SetXTitle(Unit::GetEnergyLabel("GV"));
		h_a4->SetTitle(Form("; ; BR-%d %s(R,t)/<%s(R)>", iBR_true+2426, element, element));
		// h_a4->SetTitleSize(0.5,"y"); 
		HistTools::SetStyle(h_ratio1, kPink, kFullCircle, 1.4, 1, 1);
		HistTools::SetStyle(h_ratio0, kBlue, kFullCircle, 1.4, 1, 1);
		//fit_ratio2->SetLineColor(kGreen); 
		fit_ratio->Draw("SAME");
		//fit_ratio2->Draw("SAME"); 
		h_ratio1->Draw("E1X0 SAME");
		h_ratio0->Draw("E1X0 SAME");
		l2->Draw("SAME"); 

		TH1 *h_fiterr[2]; 
		TH1D *h_fitres[2];
		//TH1D *h_fitres2[2]; 
		TLegend *l_chi2 = new TLegend(0.62,0.8,0.9,1.0);

		for (UShort_t i = 0; i < 2; ++i)
		{
  			TH1 *hist = HistTools::ToHist(data[i]);
			h_fiterr[i] = (TH1*) HistTools::GetFitError( (TH1*) data[i], fit_ratio, "_fiterr", false, false, true, 10, DBL_MAX, 0.68, &fitter);
			//h_fiterr[i]->Print("range"); 
   			h_fitres[i] = (TH1D *)HistTools::GetResiduals(hist, fit_ratio, "_fitres", false, true, true, 5, 1, 0.68, &fitter);
			HistTools::CopyStyle(hist, h_fitres[i]);

   			UShort_t ndf  = hist->GetNbinsX(); 
   			Double_t chi2 = chi2norm[i]*ndf; 
   			// printf(" %d nodes ### %-34s   chi2/ndf=%6.2f/%-2u   chi2norm=%5.2f   prob.=%5.2f%%\n", nnodes, hist->GetTitle(), chi2, ndf, chi2norm[i], TMath::Prob(chi2, ndf)*1e2);

			l_chi2->AddEntry(h_fitres[i], Form("1+[0]*exp(-[1]*R) chi2/ndf=%6.2f/%-2u", chi2, ndf), "P");  

			if (i==0){ 
				g_chi2_ace->SetPoint(iBR, UBRToTime(iBR_true+2426), chi2/ndf);
				for (int bin=0; bin <= h_fitres[i]->GetNbinsX(); ++bin){ 
					if (bin%2==0){ 
						g_ratio_ace[bin/2]->SetPoint(iBR, UBRToTime(iBR_true+2426), h_ratio1->GetBinContent(bin+1)); 
						g_ratio_ace[bin/2]->SetPointError(iBR, 0, h_ratio1->GetBinError(bin+1)); 
						g_ratio_ace_fit[bin/2]->SetPoint(iBR, UBRToTime(iBR_true+2426), h_fiterr[i]->GetBinContent(bin+1)); 
						g_ratio_ace_fit[bin/2]->SetPointError(iBR, 0, h_fiterr[i]->GetBinError(bin+1)); 
						g_residual_ace[bin/2]->SetPoint(iBR, UBRToTime(iBR_true+2426), h_fitres[i]->GetBinContent(bin+1));
						g_residual_ace[bin/2]->SetPointError(iBR, 0, h_fitres[i]->GetBinError(bin+1)); 

						double df = h_fitres[i]->GetBinContent(bin+1)*hist->GetBinContent(bin+1)/hist->GetBinError(bin+1); // (data-fit)/error_data
						double ddf = df*sqrt(pow(h_fitres[i]->GetBinError(bin+1)/h_fitres[i]->GetBinContent(bin+1))+pow(h_fiterr[i]->GetBinError(bin+1)/h_fiterr[i]->GetBinContent(bin+1))); // d[(data-fit)/erro_data]  
						//printf(" df = %0.4f, a= %0.4f, b=%0.4f \n", df, h_fitres[i]->GetBinContent(bin+1)*h_fiterr[i]->GetBinContent(bin+1), h_fiterr[i]->GetBinError(bin+1)); 
						g_resierr_ace[bin/2]->SetPoint(iBR, UBRToTime(iBR_true+2426), df); 
						g_resierr_ace[bin/2]->SetPointError(iBR, 0, ddf); 
 
						// printf("bin/2 = %d, iBR = %d, time = %10.4u, fitres = %10.4f \n", bin/2, iBR, UBRToTime(iBR_true+2426), h_fitres[0]->GetBinContent(bin+1)); 
					} 
				}   
			} 

			if (i==1){ 
				g_chi2_ams->SetPoint(iBR, UBRToTime(iBR_true+2426), chi2/ndf);
				for (int bin=0; bin <= h_ams_new->GetNbinsX()-1; ++bin){ 
					g_ratio_ams[bin]->SetPoint(iBR, UBRToTime(iBR_true+2426), h_ratio0->GetBinContent(bin+1));
					g_ratio_ams[bin]->SetPointError(iBR, 0, h_ratio0->GetBinError(bin+1));  
					g_ratio_ams_fit[bin]->SetPoint(iBR, UBRToTime(iBR_true+2426), h_fiterr[i]->GetBinContent(bin+1)); 
					g_ratio_ams_fit[bin]->SetPointError(iBR, 0, h_fiterr[i]->GetBinError(bin+1)); 
					g_residual_ams[bin]->SetPoint(iBR, UBRToTime(iBR_true+2426), h_fitres[i]->GetBinContent(bin+1));
					g_residual_ams[bin]->SetPointError(iBR, 0, h_fitres[i]->GetBinError(bin+1)); 
					// printf("bin = %d, iBR = %d, time = %10.4u, fitres = %10.4f \n", bin, iBR, UBRToTime(iBR_true+2426), h_fitres[1]->GetBinContent(bin+1)); 
				}  
			} 	
		}	

/*
		for (UShort_t i = 0; i < 2; ++i)
		{
  			TH1 *hist = HistTools::ToHist(data[i]);
   			h_fitres2[i] = (TH1D *)HistTools::GetResiduals(hist, fit_ratio2, "_fitres2", false, true, true, 5, 1, 0.68, &fitter2);
			HistTools::CopyStyle(hist, h_fitres2[i]);

   			UShort_t ndf  = hist->GetNbinsX();
   			Double_t chi2 = chi2norm2[i]*ndf; 
   			// printf(" %d nodes ### %-34s   chi2/ndf=%6.2f/%-2u   chi2norm=%5.2f   prob.=%5.2f%%\n", nnodes, hist->GetTitle(), chi2, ndf, chi2norm[i], TMath::Prob(chi2, ndf)*1e2); 

			l_chi2->AddEntry(h_fitres2[i], Form("1+[0]*exp(-[1]*ln(R)) chi2/ndf=%6.2f/%-2u", chi2, ndf), "P"); 

			if (i==0){ 
				g_chi2_ace2->SetPoint(iBR, UBRToTime(iBR_true+2426), chi2/ndf);
			} 

			if (i==1){ 
				g_chi2_ams2->SetPoint(iBR, UBRToTime(iBR_true+2426), chi2/ndf); 
			} 	
		}
*/

		TH1 *h_comp1 = (TH1*) h_ace_BR->Clone("h_comp1");
		h_comp1->Divide(rescaled_fit); 

		double comp1_sum=0, dcomp1_sum=0; // compute average ratio of data/template-1 of ACE F(R,t)/<F(R)> manually 		

		for(int nbin=0;nbin<h_ace_rig_ave->GetNbinsX();++nbin){
			comp1_sum += h_comp1->GetBinContent(nbin); 
			dcomp1_sum += h_comp1->GetBinError(nbin)*h_comp1->GetBinError(nbin); 
			//printf("comp1_sum = %0.6f \n", comp1_sum);
		}
		double comp1_ave = comp1_sum/7; 
		double dcomp1_ave = sqrt(dcomp1_sum)/7; 

		g_ratio_ave->SetPoint(iBR, UBRToTime(iBR_true+2426), comp1_ave-1); 
		g_ratio_ave->SetPointError(iBR, 0, dcomp1_ave); 

		c2->cd(2);
		gPad->SetGrid();
		gPad->SetLogx(); 
		gPad->SetTopMargin(0.01); 

		h_a5->GetYaxis()->SetRangeUser(-0.15, 0.15)  ; 
		h_a5->GetXaxis()->SetLimits(0.7, 60); 

		h_a5->Draw("E1X0"); 
		h_a5->SetXTitle(Unit::GetEnergyLabel("GV"));
		h_a5->SetTitle("; ; Data/Fit-1"); 

		//h_fitres2[0]->Draw("E1X0 SAME"); 
		//h_fitres2[1]->Draw("E1X0 SAME");

		h_fitres[0]->Draw("E1X0 SAME");
		h_fitres[1]->Draw("E1X0 SAME");

		//HistTools::SetStyle(h_fitres2[0], kGreen-1, kFullCircle, 1.4, 1, 1);
		//HistTools::SetStyle(h_fitres2[1], kGreen-1, kFullCircle, 1.4, 1, 1); 

		l_chi2->Draw("SAME"); 		

		c2->cd(3);
		gPad->SetGrid();
		gPad->SetLogx(); 

		TLegend *l_both3 = new TLegend(0.62,0.7,0.8,1.0); 
		l_both3->AddEntry(gspind_ace2, Form("ACE %s(R,t)/<%s(R)> Spectral Indices", element, element), "PL");
		l_both3->AddEntry(gspind_ams2, Form("AMS %s(R,t)/<%s(R)> Spectral Indices", element, element), "PL"); 
		// l_both2->AddEntry(fspind_he, "AMS He(R,t)_Fit/temp<He(R)>_Fit Spectral Indices", "L");  	
		l_both3->AddEntry(fspind_ratio, "Spectral Index of 1+[0]*exp(-[1]*R) Fit to He(R,t)/<He(R)>", "L");	
		l_both3->AddEntry(fspind_he, "AMS He(R,t)_Fit/<He(R)>_Fit Spectral Indices", "L");  
		//l_both3->AddEntry(fspind_ratio2, "Spectral Index of 1+[0]*exp(-[1]*ln(R)) Fit to He(R,t)/<He(R)>", "L"); 

		//if (get_temp_int(element)<4) gspind_ace2->GetYaxis()->SetRangeUser(-2, 1.3); 
		TAxis *axis5 = gspind_ace2->GetXaxis(); 
		axis5->SetLimits(0.7, 60.); 

		gspind_ace2->SetTitle(Form("; ; Spectral Indices Model of %s F(R,t)/<F(R)> at BR-%d", element, 2426+iBR_true));
		gspind_ace2->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV"));
		gspind_ace2->Draw("APL");  
		fspind_ratio->Draw("SAME"); 
		//fspind_ratio2->Draw("SAME"); 
		//fspind_ratio2->SetLineColor(kGreen); 
		fspind_he->Draw("SAME");
		gspind_ams2->Draw("PL SAME"); 
		gspind_ace2->Draw("PL SAME");
		l_both3->Draw("SAME");   

		rescaled_fit->Write(Form("rescaled_fit_BR%d", 2426+iBR_true)); 
		fit_ratio->Write(Form("fit_ratio_BR%d", 2426+iBR_true)); 
		//fit_ratio2->Write(Form("fit_ratio2_BR%d", 2426+iBR_true)); 
		h_ratio1->Write(Form("h_ratio1_BR%d", 2426+iBR_true));
		h_ratio0->Write(Form("h_ratio0_BR%d", 2426+iBR_true)); 
		h_fitres[0]->Write(Form("h_fitres_ace_BR%d", 2426+iBR_true));
		h_fitres[1]->Write(Form("h_fitres_ams_BR%d", 2426+iBR_true)); 

		if (iBR==0) c2->Print(Form("./data/ACE/fill/fake_td_ams/spind_model_%s_temp.pdf(", element), "pdf"); 
		if (iBR>0 && iBR<nBRs-1) c2->Print(Form("./data/ACE/fill/fake_td_ams/spind_model_%s_temp.pdf", element), "pdf"); 
		if (iBR==nBRs-1){
			//gspind_ace2->Print("range");
			cout << " " << endl; 
			//gspind_ams2->Print("range");
			cout << " " << endl; 
			// h_ratio2->Print("range");
			c2->Print(Form("./data/ACE/fill/fake_td_ams/spind_model_%s_temp.pdf)", element), "pdf"); 
		}

		printf(" BR = %d \n", iBR_true);  

	   	if (iBR+2426==2472-1) iBR_true += 3; 
		else iBR_true ++; 
	}  

	g_chi2_ace->Write("g_chi2_ace"); 
	g_chi2_ams->Write("g_chi2_ams"); 
			
	c3->cd(1);
	gPad->SetGrid(); 
	HistTools::SetStyle(g_chi2_ace, kPink, kFullCircle, 1.4, 1, 1); 
	HistTools::SetStyle(g_chi2_ace2, kBlue, kFullCircle, 1.4, 1, 1); 
	TLegend *l_c3_1 = new TLegend(0.62,0.8,0.9,0.9); 
	l_c3_1->AddEntry(g_chi2_ace, "1+[0]*exp(-[1]*R)", "PL"); 
	l_c3_1->AddEntry(g_chi2_ace2, "1+[0]*exp(-[1]*ln(R))", "PL"); 
	g_chi2_ace->SetTitle(Form("; ; ACE %s(R,t)/<%s(R)> Fitting Chi2/NDF", element, element)); 
	g_chi2_ace->GetXaxis()->SetTimeDisplay(1);
	g_chi2_ace->GetXaxis()->SetTimeFormat("%m-%y"); 
	g_chi2_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	g_chi2_ace->GetXaxis()->SetTitleSize(0.7); 
	// //if (get_temp_int(element)<4) g_chi2_ace->GetYaxis()->SetRangeUser(0., 4.5); 
	g_chi2_ace->Draw("APL"); 
	g_chi2_ace2->Draw("PL SAME"); 
	l_c3_1->Draw("SAME");

	c3->cd(2); 
	gPad->SetGrid(); 
	HistTools::SetStyle(g_chi2_ams, kPink, kFullCircle, 1.4, 1, 1); 
	HistTools::SetStyle(g_chi2_ams2, kBlue, kFullCircle, 1.4, 1, 1); 
	TLegend *l_c3_2 = new TLegend(0.62,0.8,0.9,0.9); 
	l_c3_2->AddEntry(g_chi2_ams, "1+[0]*exp(-[1]*R)", "PL");
	l_c3_2->AddEntry(g_chi2_ams2, "1+[0]*exp(-[1]*ln(R))", "PL");
	g_chi2_ams->SetTitle(Form("; ; AMS %s(R,t)/<%s(R)> Fitting Chi2/NDF", element, element)); 
	g_chi2_ams->GetXaxis()->SetTimeDisplay(1);
	g_chi2_ams->GetXaxis()->SetTimeFormat("%m-%y");
	g_chi2_ams->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	g_chi2_ams->GetXaxis()->SetTitleSize(0.7);
	// //if (get_temp_int(element)<4) g_chi2_ams->GetYaxis()->SetRangeUser(0., 2.2); 
	g_chi2_ams->Draw("APL");
	g_chi2_ams2->Draw("PL SAME");  
	l_c3_2->Draw("SAME"); 

	g_chi2_ace->Print();
	g_chi2_ams->Print(); 

	c3->Print(Form("./data/ACE/fill/fake_td_ams/chi2_vs_BR_%s.png", element)); 

	g_ratio_ave_norm = get_norm_graph( g_ratio_ave ); 

	for (int iBR=0; iBR<nBRs; ++iBR){
				Double_t x=0, y=0; 
				g_ratio_ave_norm->GetPoint(iBR, x, y); 
				g_ratio_ave_norm->SetPoint(iBR, x, y-1); 			 
	}	

	TGraphErrors *g_ave_residual_ace_rig = new TGraphErrors(); // time-averaged residual vs. rigidity
	TGraphErrors *g_ave_residual_ams_rig = new TGraphErrors(); // time-averaged residual vs. rigidity  

	for (int bin=0; bin <= h_ace_rig_ave->GetNbinsX(); ++bin){
 
		if (bin%2==0) {  

			g_ratio_time_norm[bin/2] = get_norm_graph( g_ratio_time[bin/2] ); 
		
			for (int iBR=0; iBR<g_ratio_time_norm[bin/2]->GetN(); ++iBR){
				Double_t x=0, y=0; 
				g_ratio_time_norm[bin/2]->GetPoint(iBR, x, y); 
				g_ratio_time_norm[bin/2]->SetPoint(iBR, x, y-1); 			 
			}

			g_ave_residual_ace_rig->SetPoint(bin/2, h_ace_rig_ave->GetBinCenter(bin+1), g_residual_ace[bin/2]->GetMean(2) ); 
			g_ave_residual_ace_rig->SetPointError(bin/2, 0, g_residual_ace[bin/2]->GetRMS(2) ); 

			c4->cd(1);
			gPad->SetGrid(); 

			TLegend *l_comp1 = new TLegend(0.62, 0.8, 0.9, 1.0); 
			l_comp1->AddEntry(g_ratio_ave_norm, "Normalized Averaged ACE Data/Template-1", "PL"); 
			l_comp1->AddEntry(g_residual_ace[bin/2], Form("ACE Fitting Residuals (%0.4f GV)", h_ace_rig_ave->GetBinCenter(bin+1)), "PL"); 

			g_residual_ace[bin/2]->SetTitle(Form(" ; ; ACE %s Data/Template-1 vs. Fitting Residuals", element));  
			g_residual_ace[bin/2]->GetXaxis()->SetTimeDisplay(1);
			g_residual_ace[bin/2]->GetXaxis()->SetTimeFormat("%m-%y");
			g_residual_ace[bin/2]->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g_residual_ace[bin/2]->GetXaxis()->SetTitleSize(0.7);
			//if (get_temp_int(element)<4) g_residual_ace[bin/2]->GetYaxis()->SetRangeUser(-0.25, 0.3); 
			g_residual_ace[bin/2]->Draw("APL");
			// g_ratio_ave_norm->Draw("PL SAME"); 
			g_ratio_time_norm[bin/2]->Draw("PL SAME"); 
			l_comp1->Draw("SAME"); 		

			if (bin==0) c4->Print(Form("./data/ACE/fill/fake_td_ams/fitres_ace_vs_BR_%s.pdf(", element), "pdf"); 
			if (bin>0 && bin<h_ace_rig_ave->GetNbinsX()-1) c4->Print(Form("./data/ACE/fill/fake_td_ams/fitres_ace_vs_BR_%s.pdf", element), "pdf");  
			if (bin==h_ace_rig_ave->GetNbinsX()-1){
				c4->Print(Form("./data/ACE/fill/fake_td_ams/fitres_ace_vs_BR_%s.pdf)", element), "pdf"); 
			} 
			
			g_ratio_ace[bin/2]->Write(Form("g_ratio_ace_%d", bin/2)); 
			g_ratio_ace_fit[bin/2]->Write(Form("g_ratio_ace_fit_%d", bin/2)); 
			//g_ratio_ace_fit[bin/2]->Print("range"); 
			g_residual_ace[bin/2]->Write(Form("g_residual_ace_%d", bin/2)); 
			g_resierr_ace[bin/2]->Write(Form("g_resierr_ace_%d", bin/2)); 
			//g_resierr_ace[bin/2]->Print("range"); 

		} 
	} 

	for (int bin=0; bin <= h_ams_new->GetNbinsX()-1; ++bin){ 

			g_ave_residual_ams_rig->SetPoint(bin, h_ams_new->GetBinCenter(bin+1), g_residual_ams[bin]->GetMean(2) ); 
			g_ave_residual_ams_rig->SetPointError(bin, 0, g_residual_ams[bin]->GetRMS(2) ); 

			c4_2->cd(1); 
			gPad->SetGrid(); 

			g_residual_ams[bin]->SetTitle(Form("AMS %s(R,t)/<%s(R)>; ; Fitting Residuals (%0.4f GV)", element, element, h_ams_new->GetBinCenter(bin+1)));  
			g_residual_ams[bin]->GetXaxis()->SetTimeDisplay(1);
			g_residual_ams[bin]->GetXaxis()->SetTimeFormat("%m-%y");
			g_residual_ams[bin]->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g_residual_ams[bin]->GetXaxis()->SetTitleSize(0.7);
			////if (get_temp_int(element)<4) g_residual_ams[bin]->GetYaxis()->SetRangeUser(-0.06, 0.1);  
			// g_residual_ams[bin]->Print("range"); 
			g_residual_ams[bin]->Draw("APL"); 

			if (bin==0) c4_2->Print(Form("./data/ACE/fill/fake_td_ams/fitres_ams_vs_BR_%s.pdf(", element), "pdf"); 
			if (bin>0 && bin<h_ams_new->GetNbinsX()-1) c4_2->Print(Form("./data/ACE/fill/fake_td_ams/fitres_ams_vs_BR_%s.pdf", element), "pdf");  
			if (bin==h_ams_new->GetNbinsX()-1){
				c4_2->Print(Form("./data/ACE/fill/fake_td_ams/fitres_ams_vs_BR_%s.pdf)", element), "pdf"); 
			} 
			
			g_ratio_ams[bin]->Write(Form("g_ratio_ams_%d", bin)); 
			g_ratio_ams_fit[bin]->Write(Form("g_ratio_ams_fit_%d", bin)); 
			//g_ratio_ams_fit[bin]->Print("range"); 
			g_residual_ams[bin]->Write(Form("g_residual_ams_%d", bin)); 

			//g_ratio_ams[bin]->Print("range");
			//g_residual_ams[bin]->Print("range"); 

	} 

	// TGraphErrors *g_mean_residual_ace = ave_grapherrors( g_residual_ace, 7 ); 
	// TGraphErrors *g_mean_residual_ams = ave_grapherrors( g_residual_ams, h_ams_new->GetNbinsX()-1 ); 

	g_ave_residual_ace_rig->Write("g_ave_residual_ace_rig"); 
	g_ave_residual_ams_rig->Write("g_ave_residual_ams_rig"); 

	file->Close(); 

	for (int ipar=0; ipar<2; ++ipar){
		c5->cd(ipar+1); 
		gPad->SetGrid(); 
		g_pars[ipar]->Draw("APL"); 
		//g_pars[ipar]->Print("range"); 
		HistTools::SetStyle(g_pars[ipar], kPink, kFullCircle, 1.4, 1, 1); 
		g_pars[ipar]->GetXaxis()->SetTimeDisplay(1);
		g_pars[ipar]->GetXaxis()->SetTimeFormat("%m-%y");
		g_pars[ipar]->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
		g_pars[ipar]->GetXaxis()->SetTitleSize(0.7);
		////if (get_temp_int(element)<4) g_pars[ipar]->GetYaxis()->SetRangeUser(-3., 3.); 
		if (ipar==0) g_pars[ipar]->SetTitle(Form("1+[0]*exp(-[1]*R); ; Parameter [%d]", ipar));
		else g_pars[ipar]->SetTitle(Form(" ; ; Parameter [%d]", ipar));
	}		

	TGraphErrors *g_pars_ratio = new TGraphErrors(); 
	for (int i=0; i<g_pars[0]->GetN(); ++i){
		double ratio_par = g_pars[1]->GetY()[i]/g_pars[0]->GetY()[i]; 
		double ratio_par_err = g_pars[1]->GetY()[i]/g_pars[0]->GetY()[i]*sqrt(pow(g_pars[1]->GetEY()[i]/g_pars[1]->GetY()[i], 2)+pow(g_pars[0]->GetEY()[i]/g_pars[0]->GetY()[i], 2));  
		if (abs(ratio_par_err/ratio_par) < 1){ 
			g_pars_ratio->SetPoint(i, g_pars[0]->GetX()[i], g_pars[1]->GetY()[i]/g_pars[0]->GetY()[i]); 
			g_pars_ratio->SetPointError(i, 0, g_pars[1]->GetY()[i]/g_pars[0]->GetY()[i]*sqrt(pow(g_pars[1]->GetEY()[i]/g_pars[1]->GetY()[i], 2)+pow(g_pars[0]->GetEY()[i]/g_pars[0]->GetY()[i], 2))); 
		} else {
			g_pars_ratio->SetPoint(i, g_pars[0]->GetX()[i], 0); 
			g_pars_ratio->SetPointError(i, 0, 0); 
		}
	}

	c5->cd(3);
	gPad->SetGrid();
	HistTools::SetStyle(g_pars_ratio, kPink, kFullCircle, 1.4, 1, 1);
	g_pars_ratio->GetXaxis()->SetTimeDisplay(1);
	g_pars_ratio->GetXaxis()->SetTimeFormat("%m-%y");
	g_pars_ratio->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	g_pars_ratio->GetXaxis()->SetTitleSize(0.7);
	////if (get_temp_int(element)<4) g_pars_ratio->GetYaxis()->SetRangeUser(-15., 30.); 
	g_pars_ratio->SetTitle(Form(" ; ; [%d]/[%d]", 1, 0)); 
	g_pars_ratio->Draw("AP"); 

	// g_pars_ratio->Print("range"); 

	c5->Print(Form("./data/ACE/fill/fake_td_ams/pars_vs_BR_%s.png", element)); 

	TCanvas *c6 = new TCanvas("c6", "", 1600, 900); // plot time-averaged fit residuals to F(R,t)/<F(R)> vs. rigidity 
	c6->Divide(1, 1);  

	c6->cd(1); 
	gPad->SetGrid(); 
	gPad->SetLogx(); 

	HistTools::CopyStyle( h_ace_rig_ave, g_ave_residual_ace_rig ); 
	HistTools::CopyStyle( h_ams_new, g_ave_residual_ams_rig );
	////if (get_temp_int(element)<4) g_ave_residual_ace_rig->GetYaxis()->SetRangeUser(-0.085, 0.085); 
	g_ave_residual_ace_rig->GetXaxis()->SetLimits(0.7, 60);
	g_ave_residual_ace_rig->SetTitle(Form(" ; ; Time-averaged Fit Residuals of %s(R,t)/<%s(R)>", element, element)); 
	g_ave_residual_ace_rig->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV")); 

	TLegend *l_both4 = new TLegend(0.62,0.8,0.9,0.9); 
	l_both4->AddEntry(g_ave_residual_ace_rig, "ACE", "PL");
	l_both4->AddEntry(g_ave_residual_ams_rig, "AMS", "PL"); 

	g_ave_residual_ace_rig->Draw("APL"); 
	g_ave_residual_ams_rig->Draw("PLSAME"); 
	l_both4->Draw("SAME"); 

	c6->Print(Form("./data/ACE/fill/fake_td_ams/time_ave_residual_vs_rig_%s.png", element)); 

	return 0; 
}

void plot_fit_pars(){

	gStyle->SetOptStat(0);
	const int nBRs = Experiments::Info[Experiments::AMS02].Dataset[1].nMeasurements; // read the number of BRs

	TCanvas *c0 = new TCanvas("c0", "", 1600, 900); // BCN/O 
	c0->Divide(1, 2);

	TCanvas *c1 = new TCanvas("c1", "", 1600, 900); // BCN/ACE_BR
	c1->Divide(1, 2);

	TCanvas *c2 = new TCanvas("c2", "", 1600, 900);
	c2->Divide(1, 1); 

	TLegend *l0 = new TLegend(0.62,0.8,0.9,0.9);
	TLegend *l1 = new TLegend(0.62,0.8,0.9,0.9); 

	TF1 *f_corr2[4][7]; 
	TH1 *h_correrr[4][7]; 
	TH1 *h_corr2[4]; 

	// pars vs time 	
	for (int i=0; i<4; ++i){ 

		const char *element = ACE_Element[i]; 
		Particle::Type isotope = ACE_Isotope[i]; 

		double *kin_bins = get_kin_bins(element);
       		double *SpallCorr = get_spall_corr(element);
		double *SpallCorrUnc = get_spall_corr_unc(element);
		double *EMed = get_EMed(element);

		const UInt_t FirstACEBR = 2240;
		vector<UInt_t> BRs;
		// we stop at BR 2493, which ends on 2016/05/23, just 3 days before the end of the data taking period for AMS nuclei
		for (UInt_t br=2426; br<=2493; ++br) { 
			if (br != 2472 && br != 2473) BRs.push_back(br-FirstACEBR); 
		}

		TH1 *h_ace_ene = HistTools::GraphToHist(get_ace_average_graph( element , &BRs[0], BRs.size() ), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
		TH1 *h_ace_rig = HistTools::TransformEnergyAndDifferentialFluxNew(h_ace_ene, isotope, "MeV/n cm", "GV m", "_rig"); // load averaged ACE data for the same element, converted in rigidity

		TFile *file2 = new TFile(Form("data/ACE/fill/%s_fill.root", element));

		TGraphErrors *par0_X = new TGraphErrors(Form("data/ACE/fill/fake_td_ams_v2/fit_pars_%s.dat", element), "%lg %lg %lg %*s %*s"); 
		TGraphErrors *par0_O = new TGraphErrors(Form("data/ACE/fill/fake_td_ams_v2/fit_pars_%s.dat", "O"), "%lg %lg %lg %*s %*s"); 

		TGraphErrors *par1_X = new TGraphErrors(Form("data/ACE/fill/fake_td_ams_v2/fit_pars_%s.dat", element), "%lg %*s %*s %lg %lg"); 
		TGraphErrors *par1_O = new TGraphErrors(Form("data/ACE/fill/fake_td_ams_v2/fit_pars_%s.dat", "O"), "%lg %*s %*s %lg %lg"); 	

		TGraphErrors *par0_ratio = new TGraphErrors();
		TGraphErrors *par1_ratio = new TGraphErrors(); 

		TGraphErrors *par0_ACE[7]; 
		TGraphErrors *par1_ACE[7];

		h_corr2[i] = new TH1D("", "", 7, 0, 7); 

		int nbins =7; 
		for (int bin=0; bin<nbins; ++bin){

			TLegend *l2 = new TLegend(0.62,0.8,0.9,0.9);
			TLegend *l3 = new TLegend(0.62,0.8,0.9,0.9);

			par0_ACE[bin] = new TGraphErrors();
			par1_ACE[bin] = new TGraphErrors(); 

			TGraphErrors *par1_par0_sc = new TGraphErrors(); // scatter plot, include every BR in one plot, separate for each bin 

			TGraphErrors *par0_ACE_sc = new TGraphErrors();
			TGraphErrors *par1_ACE_sc = new TGraphErrors();

			HistTools::SetStyle(par0_ACE[bin], kBlue-2, kFullCircle, 0.9, 1, 1);
			HistTools::SetStyle(par1_ACE[bin], kBlue-2, kFullCircle, 0.9, 1, 1);

			int iBR_true = 0;
			for (int iBR=0; iBR<nBRs; ++iBR){ 

				TH1 *h_ace_BR = (TH1*) file2->Get(Form("h_rig_%s_BR%d", element, 2426+iBR_true));

				par0_ACE[bin]->SetPoint(iBR, par0_X->GetX()[iBR], par0_X->GetY()[iBR]/h_ace_BR->GetBinContent(2*bin+1)); 
				par0_ACE[bin]->SetPointError(iBR, 0, par0_X->GetY()[iBR]/h_ace_BR->GetBinContent(2*bin+1)*sqrt(pow(par0_X->GetEY()[iBR]/par0_X->GetY()[iBR], 2)+pow(h_ace_BR->GetBinError(2*bin+1)/h_ace_BR->GetBinContent(2*bin+1), 2)));

				par1_ACE[bin]->SetPoint(iBR, par1_X->GetX()[iBR], par1_X->GetY()[iBR]/h_ace_BR->GetBinContent(2*bin+1)); 
				par1_ACE[bin]->SetPointError(iBR, 0, par1_X->GetY()[iBR]/h_ace_BR->GetBinContent(2*bin+1)*sqrt(pow(par1_X->GetEY()[iBR]/par1_X->GetY()[iBR], 2)+pow(h_ace_BR->GetBinError(2*bin+1)/h_ace_BR->GetBinContent(2*bin+1), 2)));

				par1_par0_sc->SetPoint(iBR, par0_X->GetY()[iBR], par1_X->GetY()[iBR]); 
				par1_par0_sc->SetPointError(iBR, par0_X->GetEY()[iBR], par1_X->GetEY()[iBR]); 	

				par0_ACE_sc->SetPoint(iBR, h_ace_BR->GetBinContent(2*bin+1), par0_X->GetY()[iBR]); 
				par0_ACE_sc->SetPointError(iBR, h_ace_BR->GetBinError(2*bin+1), par0_X->GetEY()[iBR]); 

				par1_ACE_sc->SetPoint(iBR, h_ace_BR->GetBinContent(2*bin+1), par1_X->GetY()[iBR]); 
				par1_ACE_sc->SetPointError(iBR, h_ace_BR->GetBinError(2*bin+1), par1_X->GetEY()[iBR]); 				
		
				if (iBR+2426==2472-1) iBR_true += 3; 
				else iBR_true ++; 
			} 

			l2->AddEntry(par0_ACE[bin], Form("%0.2f GV", h_ace_rig->GetBinCenter(2*bin+1)), "PL"); 
			l3->AddEntry(par0_ACE[bin], Form("%0.2f GV", h_ace_rig->GetBinCenter(2*bin+1)), "PL"); 

			c1->cd(1);
			gPad->SetGrid();

			par0_ACE[bin]->SetTitle(" 1+[0]*exp(-[1]*R) ; ; F(R,t)/<F(R)> Fit Par_BCN[0]/ACE(R,t)"); 
			par0_ACE[bin]->GetXaxis()->SetTimeDisplay(1);
			par0_ACE[bin]->GetXaxis()->SetTimeFormat("%m-%y");
			par0_ACE[bin]->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			par0_ACE[bin]->GetXaxis()->SetTitleSize(0.7);
			// par0_ACE[bin]->GetYaxis()->SetRangeUser(-0.1, 2.); 
			par0_ACE[bin]->Draw("APL");
			l2->Draw("SAME");
	
			c1->cd(2);
			gPad->SetGrid();

			par1_ACE[bin]->SetTitle("; ; F(R,t)/<F(R)> Fit Par_BCN[1]/ACE(R,t)"); 
			par1_ACE[bin]->GetXaxis()->SetTimeDisplay(1);
			par1_ACE[bin]->GetXaxis()->SetTimeFormat("%m-%y");
			par1_ACE[bin]->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			par1_ACE[bin]->GetXaxis()->SetTitleSize(0.7);
			// par0_ACE[bin]->GetYaxis()->SetRangeUser(-0.1, 2.); 
			par1_ACE[bin]->Draw("APL");
			l3->Draw("SAME"); 

			if (bin==0) c1->Print(Form("./data/ACE/fill/fake_td_ams_v2/fit_pars_%s_vs_ACE.pdf(", element), "pdf"); 
			if (bin>0 && bin<nbins-1) c1->Print(Form("./data/ACE/fill/fake_td_ams_v2/fit_pars_%s_vs_ACE.pdf", element), "pdf"); 
			if (bin==nbins-1){
				c1->Print(Form("./data/ACE/fill/fake_td_ams_v2/fit_pars_%s_vs_ACE.pdf)", element), "pdf"); 
			}

			c2->cd(1);  
			gPad->SetGrid(); 

			TF1 *f_rand = new TF1("f_rand","abs(sin(x)/x)*sqrt(x)",0,1); // random parameters 
			double mrand = f_rand->GetRandom();
			double nrand = f_rand->GetRandom(); 

			double corr = par1_par0_sc->GetCorrelationFactor();

			double corr_low=0, corr_up=0; 
			// these two variables will contain the lower and upper value of the correlation, accounting for statistical uncertainty at a given confidence level defined below
	
			double cl = TMath::Erf(1/sqrt(2)); 
			// confidence level corresponding to 1 sigma, roughly 0.68, which means ~68% probability

			double corr_pvalue = HistTools::CorrelationZTest(corr, par1_par0_sc->GetN(), corr_low, corr_up, cl); 
			// the p-value is a number that tells you how much is probable that the correlation factor you computed is different from zero simply by chance, and not by some physical reason

			double corr_sign = TMath::ErfcInverse(corr_pvalue); 
			// get the statistical significance of the correlation being different from zero: the significance is given in unit of sigma. Any number above 3 is very significant, i.e. we can claim that most probably the computed correlation factor is different from zero not simply by chance.

			TF1 *f_corr1 = new TF1("f_corr1", "[0]*x+[1]", -1, 2);
			f_corr1->SetParameters(10, 1.5); 
			par1_par0_sc->Fit(f_corr1, "NQ"); 

			TLegend *l_corr1 = new TLegend(0.62,0.8,0.9,0.9); 
			l_corr1->AddEntry(f_corr1, Form("%0.3fpm%0.3fx+(%0.3fpm%0.3f)", f_corr1->GetParameter(0), f_corr1->GetParError(0), f_corr1->GetParameter(1), f_corr1->GetParError(1))); 
					
			par1_par0_sc->SetTitle(Form(" %s Flux 1+[0]*exp(-[1]*R) Fit ( Corr. F=%0.3f, %0.2f GV ); Par_%s[0]; Par_%s[1]", element, corr, h_ace_rig->GetBinCenter(2*bin+1), element, element)); 
			//par1_par0_sc->GetXaxis()->SetTitleSize(1.3);
			HistTools::SetStyle(par1_par0_sc, kBlue, kFullCircle, 1.2, 1, 1);
			par1_par0_sc->SetLineStyle(2); 
			par1_par0_sc->SetLineWidth(1.0); 
			par1_par0_sc->GetYaxis()->SetRangeUser(-0.5, 10); 
			par1_par0_sc->GetXaxis()->SetRangeUser(-1, 2); 
			par1_par0_sc->Draw("AP"); 
			f_corr1->Draw("SAME"); 
			par1_par0_sc->Draw("PSAME"); 
			l_corr1->Draw("SAME"); 

			// par1_par0_sc->Print("range"); 

			if (bin==0) c2->Print(Form("./data/ACE/fill/fake_td_ams_v2/fit_pars_sca_%s_1.pdf(", element), "pdf"); 
			if (bin>0 && bin<nbins-1) c2->Print(Form("./data/ACE/fill/fake_td_ams_v2/fit_pars_sca_%s_1.pdf", element), "pdf"); 
			if (bin==nbins-1){
				c2->Print(Form("./data/ACE/fill/fake_td_ams_v2/fit_pars_sca_%s_1.pdf)", element), "pdf"); 
			} 

			c2->cd(1);
			gPad->SetGrid();

			corr = par0_ACE_sc->GetCorrelationFactor(); 
			corr_low=0, corr_up=0; 
			cl = TMath::Erf(1/sqrt(2)); 
			corr_pvalue = HistTools::CorrelationZTest(corr, par0_ACE_sc->GetN(), corr_low, corr_up, cl); 
			corr_sign = TMath::ErfcInverse(corr_pvalue); 
		
			f_corr2[i][bin] = new TF1("f_corr2[i][bin]", "[0]*x+[1]", -1, 2);
			f_corr2[i][bin]->SetParameters(10, 1); 
			par0_ACE_sc->Fit(f_corr2[i][bin], "NQ"); 

			TLegend *l_corr2 = new TLegend(0.62,0.8,0.9,0.9); 
			l_corr2->AddEntry(f_corr2[i][bin], Form("%0.3fpm%0.3fx+(%0.3fpm%0.3f)", f_corr2[i][bin]->GetParameter(0), f_corr2[i][bin]->GetParError(0), f_corr2[i][bin]->GetParameter(1), f_corr2[i][bin]->GetParameter(1)));

			// printf(" corr = %0.3f \n", corr); 
			h_corr2[i]->SetBinContent(bin+1, corr); 

			// TH1 *h_par0_ACE_sc = HistTools::ToHist(par0_ACE_sc); 
			h_correrr[i][bin] = (TH1*) HistTools::GetFitError( par0_ACE_sc, f_corr2[i][bin], "_f_correrr", false, false, true, 10, DBL_MAX, 0.68);
			//h_correrr[i][bin]->Print("range"); 

			par0_ACE_sc->SetTitle(Form(" %s Flux 1+[0]*exp(-[1]*R) Fit (Corr. F=%0.3f, %0.2f GV); ; Par_%s[0]", element, corr, h_ace_rig->GetBinCenter(2*bin+1), element));
			par0_ACE_sc->GetXaxis()->SetTitle(Form("ACE %s %s", element, Unit::GetDifferentialFluxLabel("GV m"))); 
			//par0_ACE_sc->GetXaxis()->SetTitleSize(1.3);
			HistTools::SetStyle(par0_ACE_sc, kRed, kFullCircle, 1.2, 1, 1);
			par0_ACE_sc->SetLineStyle(2); 
			par0_ACE_sc->SetLineWidth(1.0); 
			par0_ACE_sc->GetYaxis()->SetRangeUser(-0.5, 10); 
			par0_ACE_sc->Draw("AP");
			f_corr2[i][bin]->Draw("SAME"); 
			par0_ACE_sc->Draw("PSAME"); 
			l_corr2->Draw("SAME"); 			

			if (bin==0) c2->Print(Form("./data/ACE/fill/fake_td_ams_v2/fit_pars_sca_%s_2.pdf(", element), "pdf"); 
			if (bin>0 && bin<nbins-1) c2->Print(Form("./data/ACE/fill/fake_td_ams_v2/fit_pars_sca_%s_2.pdf", element), "pdf"); 
			if (bin==nbins-1){
				c2->Print(Form("./data/ACE/fill/fake_td_ams_v2/fit_pars_sca_%s_2.pdf)", element), "pdf"); 
			} 

			c2->cd(1);
			gPad->SetGrid(); 

			corr = par1_ACE_sc->GetCorrelationFactor(); 
			corr_low=0, corr_up=0; 
			cl = TMath::Erf(1/sqrt(2)); 
			corr_pvalue = HistTools::CorrelationZTest(corr, par1_ACE_sc->GetN(), corr_low, corr_up, cl); 
			corr_sign = TMath::ErfcInverse(corr_pvalue); 

			TF1 *f_corr3 = new TF1("f_corr3", "[0]*x+[1]", -1, 2);
			f_corr3->SetParameters(2, 0.3); 
			par1_ACE_sc->Fit(f_corr3, "NQ"); 

			TLegend *l_corr3 = new TLegend(0.62,0.8,0.9,0.9); 
			l_corr3->AddEntry(f_corr3, Form("%0.3fpm%0.3fx+(%0.3fpm%0.3f)", f_corr3->GetParameter(0), f_corr3->GetParError(0), f_corr3->GetParameter(1), f_corr3->GetParameter(1)));

			par1_ACE_sc->SetTitle(Form(" %s Flux 1+[0]*exp(-[1]*R) Fit (Corr. F=%0.3f, %0.2f GV); ; Par_%s[1]", element, corr, h_ace_rig->GetBinCenter(2*bin+1), element));  
			par1_ACE_sc->GetXaxis()->SetTitle(Form("ACE %s %s", element, Unit::GetDifferentialFluxLabel("GV m")));
			//par1_ACE_sc->GetXaxis()->SetTitleSize(1.3);
			HistTools::SetStyle(par1_ACE_sc, kRed, kFullCircle, 1.2, 1, 1);
			par1_ACE_sc->SetLineStyle(2);
			par1_ACE_sc->SetLineWidth(1.0);  
			par1_ACE_sc->GetYaxis()->SetRangeUser(-0.5, 10); 
			par1_ACE_sc->Draw("AP"); 
			f_corr3->Draw("SAME"); 
			par1_ACE_sc->Draw("PSAME");
			l_corr3->Draw("SAME"); 

			if (bin==0) c2->Print(Form("./data/ACE/fill/fake_td_ams_v2/fit_pars_sca_%s_3.pdf(", element), "pdf"); 
			if (bin>0 && bin<nbins-1) c2->Print(Form("./data/ACE/fill/fake_td_ams_v2/fit_pars_sca_%s_3.pdf", element), "pdf"); 
			if (bin==nbins-1){
				c2->Print(Form("./data/ACE/fill/fake_td_ams_v2/fit_pars_sca_%s_3.pdf)", element), "pdf"); 
			} 
		} 

		h_corr2[i]->Print("range"); 
		
		HistTools::SetStyle(par0_ratio, HistTools::GetColorPalette(3*i, 9), kFullCircle, 0.9, 1, 1);
		HistTools::SetStyle(par1_ratio, HistTools::GetColorPalette(3*i, 9), kFullCircle, 0.9, 1, 1);

		if (i<3) l0->AddEntry(par0_ratio, Form("Par_%s[0]/Par_O[0]", element), "PL"); 
		if (i<3) l1->AddEntry(par1_ratio, Form("Par_%s[1]/Par_O[1]", element), "PL"); 

		for (int iBR=0; iBR<nBRs; ++iBR){

			par0_ratio->SetPoint(iBR, par0_X->GetX()[iBR], par0_X->GetY()[iBR]/par0_O->GetY()[iBR]);
			par0_ratio->SetPointError(iBR, 0, par0_X->GetY()[iBR]/par0_O->GetY()[iBR]*sqrt(pow(par0_X->GetEY()[iBR]/par0_X->GetY()[iBR], 2)+pow(par0_O->GetEY()[iBR]/par0_O->GetY()[iBR], 2)));

			par1_ratio->SetPoint(iBR, par1_X->GetX()[iBR], par1_X->GetY()[iBR]/par1_O->GetY()[iBR]);
			par1_ratio->SetPointError(iBR, 0, par1_X->GetY()[iBR]/par1_O->GetY()[iBR]*sqrt(pow(par1_X->GetEY()[iBR]/par1_X->GetY()[iBR], 2)+pow(par1_O->GetEY()[iBR]/par1_O->GetY()[iBR], 2)));		

		}	

		// par0_X->Print("range");
		// par0_O->Print("range");  

		// par0_ratio->Print("range");
		// par1_ratio->Print("range"); 

		printf(" \n");

		c0->cd(1);
		gPad->SetGrid();

		par0_ratio->SetTitle(" 1+[0]*exp(-[1]*R) ; ; F(R,t)/<F(R)> Fit Par_BCN[0]/Par_O[0]"); 
		par0_ratio->GetXaxis()->SetTimeDisplay(1);
		par0_ratio->GetXaxis()->SetTimeFormat("%m-%y");
		par0_ratio->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
		par0_ratio->GetXaxis()->SetTitleSize(0.7);
		par0_ratio->GetYaxis()->SetRangeUser(-0.1, 2.); 
		if (i==0) par0_ratio->Draw("APL");
		if (i<3) par0_ratio->Draw("PL SAME"); 
		if (i==2) l0->Draw("SAME"); 
		
		c0->cd(2); 
		gPad->SetGrid(); 

		par1_ratio->SetTitle("; ; F(R,t)/<F(R)> Fit Par_BCN[1]/Par_O[1]"); 
		par1_ratio->GetXaxis()->SetTimeDisplay(1);
		par1_ratio->GetXaxis()->SetTimeFormat("%m-%y");
		par1_ratio->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
		par1_ratio->GetXaxis()->SetTitleSize(0.7);
		par1_ratio->GetYaxis()->SetRangeUser(-0.1, 2.5); 
		if (i==0) par1_ratio->Draw("APL"); 
		if (i<3) par1_ratio->Draw("PL SAME"); 
		if (i==2) l1->Draw("SAME"); 

	} 

	TFile *file1 = new TFile("./data/ACE/fill/fake_td_ams_v2/par0_fit.root", "RECREATE");

	for (int i=0; i<4; ++i){
		const char *element = ACE_Element[i]; 
		h_corr2[i]->Write(Form("h_corr2_%s", element)); 
		f_corr2[i][h_corr2[i]->GetMaximumBin()-1]->Write(Form("f_corr2_%s", element));
		h_correrr[i][h_corr2[i]->GetMaximumBin()-1]->Write(Form("h_correrr_%s", element)); 
		printf("best at bin # %d, f_corr par0 = %0.3f, par1 = %0.3f  \n ", h_corr2[i]->GetMaximumBin()-1, f_corr2[i][h_corr2[i]->GetMaximumBin()-1]->GetParameter(0) , f_corr2[i][h_corr2[i]->GetMaximumBin()-1]->GetParameter(1)); 
	} 
	file1->Close(); 

	for (int i=0; i<4; ++i){  

		// load ACE BR 
		TFile *file2 = new TFile(Form("data/ACE/fill/%s_fill.root", ACE_Element[i]));

		gROOT->ProcessLine(Form(".> data/ACE/fill/fake_td_ams_v4/compare_par0_%s.dat", ACE_Element[i]));

		int iBR_true = 0;
		for (int iBR=0; iBR<nBRs; ++iBR){

			TH1 *h_ace_BR = (TH1*) file2->Get(Form("h_rig_%s_BR%d", ACE_Element[i], 2426+iBR_true)); 
			TGraph *g_compare_par0 = new TGraph(); // compare par0 from each energy bin f_corr. 

			for (int bin=0; bin<7; ++bin){

				h_ace_BR->GetBinContent(2*bin+1); 
				g_compare_par0->SetPoint(bin, bin, f_corr2[i][bin]->Eval( h_ace_BR->GetBinContent(2*bin+1) )  ); 
				//printf("corr = %0.8f \n", f_corr2[i][bin]->Eval( h_ace_BR->GetBinContent(2*bin+1)) );  
				printf( "ACE_Flux=%0.4f, par0=%0.4f \n", h_ace_BR->GetBinContent(2*bin+1), f_corr2[i][bin]->Eval( h_ace_BR->GetBinContent(2*bin+1) ) ); 
			} 
	
			TCanvas *c_comp_par0 = new TCanvas("","",1600,900); 
			c_comp_par0->Divide(1, 1);
	
			c_comp_par0->cd(1); 	
			HistTools::SetStyle(g_compare_par0, kBlue, kFullCircle, 0.9, 1, 1); 
			//g_compare_par0->Print(); 	
			g_compare_par0->Draw("AP"); 

			printf( "BR=%d, RMS=%0.4f \n \n", 2426+iBR_true, g_compare_par0->GetRMS(2) ); 

			if (iBR==0) c_comp_par0->Print(Form("./data/ACE/fill/fake_td_ams_v4/compare_par0_%s.pdf(", ACE_Element[i]), "pdf"); 
			if (iBR>0 && iBR<nBRs-1) c_comp_par0->Print(Form("./data/ACE/fill/fake_td_ams_v4/compare_par0_%s.pdf", ACE_Element[i]), "pdf"); 
			if (iBR==nBRs-1) c_comp_par0->Print(Form("./data/ACE/fill/fake_td_ams_v4/compare_par0_%s.pdf)", ACE_Element[i]), "pdf"); 
	
	   		if (iBR+2426==2472-1) iBR_true += 3; 
			else iBR_true ++; 
		}
		gROOT->ProcessLine(".>");
	} 

	c0->Print("data/ACE/fill/fake_td_ams_v2/plot_fit_pars.png");

	TGraphErrors *par0_X[4]; // X stands for an element 
	TGraphErrors *par1_X[4];  

	TGraphErrors *par0_ratio_CO = new TGraphErrors();
	TGraphErrors *par1_ratio_CO = new TGraphErrors(); 

	TGraphErrors *par0_ratio_NB = new TGraphErrors();
	TGraphErrors *par1_ratio_NB = new TGraphErrors(); 


	c1 = new TCanvas("c1", "", 1600, 900);
	c1->Divide(2, 1); 

	for (int i=0; i<4; ++i){ 

		par0_X[i] = new TGraphErrors(Form("data/ACE/fill/fake_td_ams_v2/fit_pars_%s.dat", ACE_Element[i]), "%lg %lg %lg %*s %*s"); 
		par1_X[i] = new TGraphErrors(Form("data/ACE/fill/fake_td_ams_v2/fit_pars_%s.dat", ACE_Element[i]), "%lg %*s %*s %lg %lg"); 

		// Double_t par1_edges[nBRs+1];  
	
		//for (int iBR=0; iBR<nBRs; ++iBR){
		//	par1_edges[iBR] = par1_X[i]->GetY()[iBR]/par1_X[i]->GetMean(); 
		//} 

		TH1D *h_par1_norm = new TH1D("h_par1_norm", "", 40., 0., 10.); // par1(t)/<par1(t)>
		TH1D *h_par1_var_err = new TH1D("h_par1_var_err", "", 40., -20., 20.); // (pars[1](t)-<pars[1]>)/err_pars[1](t) 

		for (int iBR=0; iBR<nBRs; ++iBR){
			h_par1_norm->Fill(par1_X[i]->GetY()[iBR]/par1_X[i]->GetMean(2)); 
			h_par1_var_err->Fill( (par1_X[i]->GetY()[iBR]-par1_X[i]->GetMean(2))/par1_X[i]->GetEY()[iBR] ); 
			// printf("%0.4f ## %0.4f \n", par1_X[i]->GetMean(2), par1_X[i]->GetY()[iBR]/par1_X[i]->GetMean(2));  
			// printf("%0.4f ## %0.4f \n", par1_X[i]->GetMean(2), (par1_X[i]->GetY()[iBR]-par1_X[i]->GetMean(2))/par1_X[i]->GetEY()[iBR]);  
		}  

		c1->cd(1); 
		gPad->SetGrid();
		gPad->SetMargin(0.12, 0.04, 0.1, 0.08);

		TF1 *f_par1 = new TF1("f_par1", "gaus", 0., 10.); 
		f_par1->SetParameters(h_par1_norm->GetMaximum(), h_par1_norm->GetMean(), h_par1_norm->GetRMS() ); 
		h_par1_norm->Fit("f_par1");  

		TLegend *l_par1 = new TLegend(0.42,0.85,0.95,0.95); 
		l_par1->AddEntry(f_par1, Form("%0.2f*exp(-0.5*((x-%0.3f)/%0.3f)^2)", f_par1->GetParameter(0), f_par1->GetParameter(1), f_par1->GetParameter(2)), "l"); 
	
		h_par1_norm->SetFillColor(kBlue-2); 
		h_par1_norm->SetTitle(" ; par1(t)/<par1(t)>; counts"); 
		h_par1_norm->Draw("HIST"); 
		f_par1->Draw("SAME"); 
		l_par1->Draw("SAME"); 
		// h_par1_norm->Print("range"); 

		c1->cd(2); 
		gPad->SetGrid(); 
		gPad->SetMargin(0.12, 0.04, 0.1, 0.08); 
	
		TF1 *f_par1_var_err = new TF1("f_par1_var_err", "gaus", -20., 20.); 
		f_par1_var_err->SetParameters(h_par1_var_err->GetMaximum(), h_par1_var_err->GetMean(), h_par1_var_err->GetRMS() ); 
		h_par1_var_err->Fit("f_par1_var_err");  

		TLegend *l_par1_var_err = new TLegend(0.42,0.85,0.95,0.95); 
		l_par1_var_err->AddEntry(f_par1_var_err, Form("%0.2f*exp(-0.5*((x-%0.3f)/%0.3f)^2)", f_par1_var_err->GetParameter(0), f_par1_var_err->GetParameter(1), f_par1_var_err->GetParameter(2)), "l"); 
	
		h_par1_var_err->SetFillColor(kBlue-2); 
		h_par1_var_err->SetTitle(" ; (pars[1](t)-<pars[1]>)/err_pars[1](t); counts"); 
		h_par1_var_err->Draw("HIST"); 
		f_par1_var_err->Draw("SAME"); 
		l_par1_var_err->Draw("SAME"); 

		c1->Print(Form("data/ACE/fill/fake_td_ams_v2/plot_fit_pars_h_par1_var_%s.png", ACE_Element[i])); 

	}

	for (int iBR=0; iBR<nBRs; ++iBR){

		par0_ratio_CO->SetPoint(iBR, par0_X[1]->GetX()[iBR], par0_X[1]->GetY()[iBR]/par0_X[3]->GetY()[iBR]);
		par0_ratio_CO->SetPointError(iBR, 0, par0_X[1]->GetY()[iBR]/par0_X[3]->GetY()[iBR]*sqrt(pow(par0_X[1]->GetEY()[iBR]/par0_X[1]->GetY()[iBR], 2)+pow(par0_X[3]->GetEY()[iBR]/par0_X[3]->GetY()[iBR], 2)));

		par1_ratio_CO->SetPoint(iBR, par1_X[1]->GetX()[iBR], par1_X[1]->GetY()[iBR]/par1_X[3]->GetY()[iBR]);
		par1_ratio_CO->SetPointError(iBR, 0, par1_X[1]->GetY()[iBR]/par1_X[3]->GetY()[iBR]*sqrt(pow(par1_X[1]->GetEY()[iBR]/par1_X[1]->GetY()[iBR], 2)+pow(par1_X[3]->GetEY()[iBR]/par1_X[3]->GetY()[iBR], 2)));		

		par0_ratio_NB->SetPoint(iBR, par0_X[2]->GetX()[iBR], par0_X[2]->GetY()[iBR]/par0_X[0]->GetY()[iBR]);
		par0_ratio_NB->SetPointError(iBR, 0, par0_X[2]->GetY()[iBR]/par0_X[0]->GetY()[iBR]*sqrt(pow(par0_X[2]->GetEY()[iBR]/par0_X[2]->GetY()[iBR], 2)+pow(par0_X[0]->GetEY()[iBR]/par0_X[0]->GetY()[iBR], 2)));

		par1_ratio_NB->SetPoint(iBR, par1_X[2]->GetX()[iBR], par1_X[2]->GetY()[iBR]/par1_X[0]->GetY()[iBR]);
		par1_ratio_NB->SetPointError(iBR, 0, par1_X[2]->GetY()[iBR]/par1_X[0]->GetY()[iBR]*sqrt(pow(par1_X[2]->GetEY()[iBR]/par1_X[2]->GetY()[iBR], 2)+pow(par1_X[0]->GetEY()[iBR]/par1_X[0]->GetY()[iBR], 2)));

	}

	HistTools::SetStyle(par0_ratio_CO, kBlue-2, kFullCircle, 0.9, 1, 1);
	HistTools::SetStyle(par1_ratio_CO, kBlue-2, kFullCircle, 0.9, 1, 1);
	HistTools::SetStyle(par0_ratio_NB, kBlue-2, kFullCircle, 0.9, 1, 1);
	HistTools::SetStyle(par1_ratio_NB, kBlue-2, kFullCircle, 0.9, 1, 1);

	c0 = new TCanvas("c0", "", 1600, 900);
	c0->Divide(1, 2); 

	// reuse of the canvas
	c0->cd(1);
	
	par0_ratio_CO->SetTitle(" 1+[0]*exp(-[1]*R) ; ; F(R,t)/<F(R)> Fit Par_C[0]/Par_O[0]"); 
	par0_ratio_CO->GetXaxis()->SetTimeDisplay(1);
	par0_ratio_CO->GetXaxis()->SetTimeFormat("%m-%y");
	par0_ratio_CO->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	par0_ratio_CO->GetXaxis()->SetTitleSize(0.7);
	par0_ratio_CO->GetYaxis()->SetRangeUser(-0.5, 10); 
	par0_ratio_CO->Draw("APL");

	c0->cd(2);
	
	par1_ratio_CO->SetTitle(" ; ; F(R,t)/<F(R)> Fit Par_C[1]/Par_O[1]"); 
	par1_ratio_CO->GetXaxis()->SetTimeDisplay(1);
	par1_ratio_CO->GetXaxis()->SetTimeFormat("%m-%y");
	par1_ratio_CO->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	par1_ratio_CO->GetXaxis()->SetTitleSize(0.7);
	par1_ratio_CO->GetYaxis()->SetRangeUser(0, 2); 
	par1_ratio_CO->Draw("APL"); 

	c0->Print("data/ACE/fill/fake_td_ams_v2/plot_fit_pars_CO.png");

	// reuse of the canvas
	c0->cd(1);
	
	par0_ratio_NB->SetTitle(" 1+[0]*exp(-[1]*R) ; ; F(R,t)/<F(R)> Fit Par_N[0]/Par_B[0]"); 
	par0_ratio_NB->GetXaxis()->SetTimeDisplay(1);
	par0_ratio_NB->GetXaxis()->SetTimeFormat("%m-%y");
	par0_ratio_NB->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	par0_ratio_NB->GetXaxis()->SetTitleSize(0.7);
	par0_ratio_NB->GetYaxis()->SetRangeUser(-0.5, 10);
	par0_ratio_NB->Draw("APL");

	c0->cd(2);
	
	par1_ratio_NB->SetTitle(" ; ; F(R,t)/<F(R)> Fit Par_N[1]/Par_B[1]"); 
	par1_ratio_NB->GetXaxis()->SetTimeDisplay(1);
	par1_ratio_NB->GetXaxis()->SetTimeFormat("%m-%y");
	par1_ratio_NB->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	par1_ratio_NB->GetXaxis()->SetTitleSize(0.7);
	par1_ratio_NB->GetYaxis()->SetRangeUser(0, 2); 
	par1_ratio_NB->Draw("APL"); 

	c0->Print("data/ACE/fill/fake_td_ams_v2/plot_fit_pars_NB.png");

} 

// F3 
void ace_fake_td_ams_v2(const char *element, Particle::Type isotope){

	gStyle->SetOptStat(0); 
	// gStyle->SetTitleSize(26,"t");

	gSystem->mkdir("data/ACE/fill/fake_td_ams_v2", true);	

 	Experiments::DataPath = "data";	
   	int data_value[n_ele] = {0, 18, 23, 27, 31, 20, 43, 22}; // p, He, Li, Be, B, C, N, O

	double *kin_bins = get_kin_bins(element);
        double *SpallCorr = get_spall_corr(element);
	double *SpallCorrUnc = get_spall_corr_unc(element);
	double *EMed = get_EMed(element);

	const UInt_t FirstACEBR = 2240;
	vector<UInt_t> BRs;
	// we stop at BR 2493, which ends on 2016/05/23, just 3 days before the end of the data taking period for AMS nuclei
	for (UInt_t br=2426; br<=2493; ++br) { 
		if (br != 2472 && br != 2473) BRs.push_back(br-FirstACEBR); 
	}

	// load AMS He BR fluxes 
	const int nBRs = Experiments::Info[Experiments::AMS02].Dataset[1].nMeasurements; // read the number of BRs
	TH1D **h_BR_he = Experiments::GetDatasetHistograms(Experiments::AMS02, 4);
	// Create a new set of histograms which match the AMS monthly bins with the AMS integrated bins 

	TH1 *h_ams = Experiments::GetMeasurementHistogram(Experiments::AMS02, get_ams_data_value(element), 0); // load AMS data for a given template 
	TH1 *h_ams_new = (TH1D*) h_BR_he[0]->Clone("h_ams_new"); 
	// TH1 *h_he_int = (TH1*) ave_hist( h_BR_he, nBRs);  
	TH1 *h_he = Experiments::GetMeasurementHistogram(Experiments::AMS02, 18, 0);
	TH1 *h_he_int = (TH1D*) h_BR_he[0]->Clone("h_ams_new");  

	// create h_ams_new and h_he_int to match the bins  

	for (int bin=1; bin <= h_ams_new->GetNbinsX(); ++bin){ 
	   if (!strcmp(element, "B") || !strcmp(element, "C")){ 
		h_ams_new->SetBinContent(bin, h_ams->GetBinContent(bin)); 
		h_ams_new->SetBinError(bin, h_ams->GetBinError(bin));
	   }
	   if (!strcmp(element, "N") || !strcmp(element, "O")){ 
		h_ams_new->SetBinContent(bin, h_ams->GetBinContent(bin-1)); 
		h_ams_new->SetBinError(bin, h_ams->GetBinError(bin-1)); 
	   }else{ 
		h_ams_new->SetBinContent(bin, h_ams->GetBinContent(bin)); 
		h_ams_new->SetBinError(bin, h_ams->GetBinError(bin)); 
	   }
	}
	if (!strcmp(element, "N") || !strcmp(element, "O")){
		h_ams_new->SetBinContent(1, 0); 
		h_ams_new->SetBinError(1, 0);
	}
 
	for (int bin=1; bin<=h_BR_he[0]->GetNbinsX(); ++bin) { 
		h_he_int->SetBinContent(bin, h_he->GetBinContent(bin)); 
		h_he_int->SetBinError(bin, h_he->GetBinError(bin)); 
	}

	// h_ams->Print("range"); 
	// h_ams_new->Print("range"); 

	// load ACE average 
	TH1 *h_ace_ene_ave = HistTools::GraphToHist(get_ace_average_graph( element , &BRs[0], BRs.size() ), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
	TH1 *h_ace_rig_ave = HistTools::TransformEnergyAndDifferentialFluxNew(h_ace_ene_ave, isotope, "MeV/n cm", "GV m", "_rig"); // load averaged ACE data for the same element, converted in rigidity
	HistTools::SetStyle(h_ams_new, kBlue, kFullCircle, 1.4, 1, 1);
	HistTools::SetStyle(h_ace_rig_ave, kBlack, kFullCircle, 1.4, 1, 1);

	UShort_t namsbins = h_ams_new->GetNbinsX(); 
	UShort_t nacebins = h_ace_rig_ave->GetNbinsX(); 
	double R1 = h_ace_rig_ave->GetBinLowEdge(1);
	double R2 = h_ams->GetBinLowEdge(namsbins+1);

	int nnodes = 9; // combined template 
	int nnodes_ams = 7; 

	printf(" template = %s, %d \n", element, get_ams_data_value(element)); 

	// load combined fit template 
	TFile *file0 = new TFile(Form("data/ACE/compare/fit_%s_%dnodes.root", get_template(element), nnodes));  

	Spline *sp_comb = new Spline("sp_comb", nnodes, Spline::LogLog | Spline::PowerLaw);
	TF1 *fsp_comb = sp_comb->GetTF1Pointer();  
	TF1 *fit_comb = (TF1*) file0->Get("fit_both"); 

	HistTools::CopyParameters(fit_comb, fsp_comb); // error 
	double x1, x2;
	fit_comb->GetRange(x1,x2);

	x1=R1, x2=R2; 
	fsp_comb->SetRange(x1,x2);

	printf("R1=%0.4f, R2=%0.4f", R1, R2); 

	// load AMS He combined, modified spline  
	TFile *file1 = new TFile(Form("data/amsfit/fit_result_node%d.root", nnodes_ams));

	// load ACE BR 
	TFile *file2 = new TFile(Form("data/ACE/fill/%s_fill.root", element));
	TFile *file2a = new TFile(Form("data/ACE/fill/%s_fill2.root", element));  

	TGraphErrors *g_ratio_time[7]; 
	TGraphErrors *g_ratio_time_norm[7]; 

	for (int i=0; i<7; ++i){
		g_ratio_time[i] = (TGraphErrors*) file2a->Get(Form("g_ratio_time_%d", i)); 	
		HistTools::SetStyle(g_ratio_time[i], kBlue, kFullCircle, 1.4, 1, 1);
	}

	// load AMS He Model by Rescaled C Template
	TFile *file3 = new TFile(Form("data/ACE/extend2/fit_C_temp_he_%dnodes.root", nnodes));

	TCanvas *c0 = new TCanvas("c0", "", 1600, 900); // fluxes & flux/template
	c0->Divide(1, 2); 

	TCanvas *c1 = new TCanvas("c1", "", 1600, 900); // He ratio, F(R,t)/<F(R)> & spectral indices 
	c1->Divide(1, 2); 
	
	TCanvas *c2 = new TCanvas("c2", "", 1600, 900); // Modeling F(R,t)/<F(R)> 
							// fit/data of F(R,t)/<F(R)> 
	c2->Divide(1, 3); 				// spectral indices of F(R,t)/<F(R)>, here is defined by h_ratio (ACE) and h_ratio0 (AMS)

	TCanvas *c3 = new TCanvas("c3", "", 1600, 900); // plot chi2/ndf vs. BR for both ACE & AMS F(R,t)/<F(R)> fits  
	c3->Divide(2, 1);  

	TCanvas *c4 = new TCanvas("c4", "", 1600, 900); // plot for F(R,t)/<F(R)> fit residuals for each ACE bin 
	c4->Divide(1, 1); 

	TCanvas *c4_2 = new TCanvas("c4_2", "", 1600, 900); // c4 but for each AMS bin 
	c4_2->Divide(1, 1); 
	
	TCanvas *c5 = new TCanvas("c5", "", 1600, 900); // plot parameter [0] [1] vs. BR for the fit 
	c5->Divide(1, 3); 

	TH1 *h_a1 = new TH1D("", "", 3000, 0, 3000); 
	TH1 *h_a2 = new TH1D("", "", 3000, 0, 3000);
	TH1 *h_a3 = new TH1D("", "", 3000, 0, 3000);
	TH1 *h_a4 = new TH1D("", "", 3000, 0, 3000);	
	TH1 *h_a5 = new TH1D("", "", 3000, 0, 3000); 

	vector<double> last_pars;
	vector<double> last_pars2 = {0.5, 0.5};  

	TGraph *g_chi2_ace = new TGraph(); 
	TGraph *g_chi2_ams = new TGraph(); 
	TGraph *g_chi2_ace2 = new TGraph(); 
	TGraph *g_chi2_ams2 = new TGraph(); 
	TGraphErrors *g_pars[2];

	for (int ipar=0; ipar<2; ++ipar){
		g_pars[ipar] = new TGraphErrors(); 
	} 

	TGraphErrors *g_ratio_ace[7];
	TGraphErrors *g_ratio_ace_fit[7]; 
	TGraphErrors *g_ratio_ams[h_ams_new->GetNbinsX()]; 
	TGraphErrors *g_ratio_ams_fit[h_ams_new->GetNbinsX()]; 

	TGraphErrors *g_ratio_ave = new TGraphErrors();
	TGraphErrors *g_ratio_ave_norm; // normalized averaged ratio of data/template-1 vs. time  

	TGraphErrors *g_residual_ace[7];
	TGraphErrors *g_resierr_ace[7]; 
	TGraphErrors *g_residual_norm_ace[7];
	TGraphErrors *g_residual_ams[h_ams_new->GetNbinsX()]; 
	TGraphErrors *g_residual_norm_ams[h_ams_new->GetNbinsX()];

	for (int bin=0; bin <= h_ace_rig_ave->GetNbinsX(); ++bin){ 
		if (bin%2==0){  
			g_ratio_ace[bin/2] = new TGraphErrors(nBRs); 
			g_ratio_ace_fit[bin/2] = new TGraphErrors(nBRs); 
			g_residual_ace[bin/2] = new TGraphErrors(nBRs); 
			g_resierr_ace[bin/2] = new TGraphErrors(nBRs);  
			HistTools::SetStyle(g_residual_ace[bin/2], kPink, kFullCircle, 1.4, 1, 1);
			HistTools::SetStyle(g_ratio_ave, kBlue+1, kFullCircle, 1.4, 1, 1);
		} 
	}  
 
	for (int bin=0; bin <= h_ams_new->GetNbinsX(); ++bin){ 
		g_ratio_ams[bin] = new TGraphErrors(nBRs); 
		g_ratio_ams_fit[bin] = new TGraphErrors(nBRs); 
		g_residual_ams[bin] = new TGraphErrors(nBRs);  
		HistTools::SetStyle(g_residual_ams[bin], kPink, kFullCircle, 1.4, 1, 1); 
	}   

	gROOT->ProcessLine(Form(".> data/ACE/fill/fake_td_ams_v2/fit_pars_%s.dat", element)); 

	TFile *file = new TFile(Form("data/ACE/fill/F3_%s.root", element), "RECREATE");
 
	int iBR_true = 0;
	for (int iBR=0; iBR<nBRs; ++iBR){

		// AMS_He Monthly Spline 
		Spline *sp_he = new Spline("sp_he", nnodes_ams, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_he = sp_he->GetTF1Pointer();  
		TF1 *fit_he = (TF1*) file1->Get(Form("fsp_BR_he_%02d", iBR)); 

		HistTools::CopyParameters(fit_he, fsp_he); // error 
		fit_he->GetRange(x1,x2);
		fsp_he->SetRange(x1,x2);

		// AMS_He Integrated Spline  
		Spline *sp_he_ave = new Spline("sp_he_ave", nnodes_ams, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_he_ave = sp_he_ave->GetTF1Pointer();  
		TF1 *fit_he_ave = (TF1*) file1->Get("fsp_he"); 

		HistTools::CopyParameters(fit_he_ave, fsp_he_ave); // error 
		x1=0, x2=0;
		fit_he_ave->GetRange(x1,x2);
		fsp_he_ave->SetRange(x1,x2); 

		// AMS_He Model by Rescaled C ACE+AMS Combined Template  
		Spline *sp_he_temp = new Spline("sp_he_temp", nnodes, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_he_temp = sp_he_temp->GetTF1Pointer();  
		TF1 *fit_he_temp = (TF1*) file3->Get("fsp_he"); 

		HistTools::CopyParameters(fit_he_temp, fsp_he_temp); // error 
		x1=0, x2=0; 
		fit_he_temp->GetRange(x1,x2);
		fsp_he_temp->SetRange(x1,x2);  

		TF1 *fsp_he2 = HistTools::CombineTF1(fsp_he, fsp_he_ave, HistTools::Divide, "fsp_he2"); // He(R,t) Fit vs. <He(R)> Fit
		fsp_he2->SetRange(x1, x2);
		TF1 *fspind_he = HistTools::GetSpectralIndex(fsp_he2, "fspind_he", x1, x2); 
		fspind_he->SetRange(x1, x2); 

		// AMS_He(R,t)/<He(R)>
		TH1 *h_ratio0 = (TH1D*) h_BR_he[iBR]->Clone("h_ratio0"); 
		h_ratio0->Divide(h_he_int); 

		// refill AMS BR He to match the bins  
		HistTools::SetStyle(h_BR_he[iBR], kBlue, kFullCircle, 1.1, 1, 1); 

		TH1 *h_ace_BR = (TH1*) file2->Get(Form("h_rig_%s_BR%d", element, 2426+iBR_true)); 

		// rescaled <AMS_C> by AMS He BR / <AMS_He>   
		TH1 *h_ams_BR_fake = (TH1D*) h_ams_new->Clone("h_ams_BR_fake");

		h_ams_BR_fake->Multiply(h_BR_he[iBR]); 
		h_ams_BR_fake->Divide(h_he_int); 
 
		// h_ams_BR_fake->Print("range"); 

		// rescaled combined fit 
		TH1 *h_ratio;
		if (get_index(element) <4) h_ratio = (TH1 *) h_ace_BR->Clone("h_ratio"); // ACE BR/Temp 	
		if (get_index(element)>=4) h_ratio = (TH1 *) h_ace_rig_ave->Clone("h_ratio"); // ACE Ave/Temp 	
		h_ratio->Divide(fsp_comb); 

		HistTools::SetStyle(h_ratio, kPink, kFullCircle, 1.4, 1, 1);

		double ratio_sum=0; // compute average of h_ratio manually  
		for(int nbin=0;nbin<14;++nbin){
			ratio_sum += h_ratio->GetBinContent(nbin);
			//printf("ratio_sum = %0.6f \n", ratio_sum);
		}
		double ratio_ave = ratio_sum/7;

		// printf("%0.6f \n", ratio_ave);
		//HistTools::PrintFunction(fsp_comb);
			
		double scale = 1./ratio_ave;
		if (get_index(element)<4) h_ratio->Scale(scale);	

		// rescaled combined fit 
		TF1 *rescaled_fit = HistTools::CombineTF1Const(fsp_comb, ratio_ave, HistTools::MultiplyConst, "rescaled_fit", R1, R2); 
		//rescaled_fit = HistTools::CombineTF1(rescaled_fit, fsp_comb, HistTools::Divide, "fit_ratio", R1, R2);  

		TH1 *h_ratio1 = (TH1D*) h_ace_BR->Clone("h_ratio1"); // spind model  
		if (get_index(element) <4) h_ratio1->Divide(fsp_comb); 
		if (get_index(element)>=4) h_ratio1->Divide(rescaled_fit); 

		TH1 *h_ratio2 = (TH1D*) h_ams_BR_fake->Clone("h_ratio2"); // spind model 
		h_ratio2->Divide(h_ams_new);
   
		// h_ratio2->Print("range");

		// estimated AMS / resclaed combined fit
		TH1 *h_ratio3 = (TH1D *) h_ams_BR_fake->Clone("h_ratio3");		

		h_ratio3->Divide(rescaled_fit); 
		HistTools::SetStyle(h_ratio3, kPink, kFullCircle, 1.4, 1, 1);

		c0->cd(1);
		gPad->SetGrid(); 
		gPad->SetLogy();
		gPad->SetLogx(); 	

		TLegend *legend = new TLegend(0.62,0.75,0.9,0.9);
		legend->AddEntry(h_ams, Form("Integrated AMS %s Flux", element), "p"); 
		legend->AddEntry(h_ace_BR, Form("ACE %s BR Flux", element), "p"); 
		legend->AddEntry(h_ams_BR_fake, Form("Estimated AMS %s BR Flux", element), "p"); 
		legend->AddEntry(rescaled_fit, "rescaled combined template", "l"); 

		h_a1->SetTitle(Form("%s BR-%d Model Rigidity Spectrum; ; ", element, 2426+iBR_true));
		h_a1->SetXTitle(Unit::GetEnergyLabel("GV"));
  		h_a1->SetYTitle(Unit::GetDifferentialFluxLabel("GV m"));

		////if (get_temp_int(element)<4) h_a1->GetYaxis()->SetRangeUser(1e-3, 1e1); 
		// h_a1->GetXaxis()->SetRangeUser(0.01, 3000.); 

		TAxis *axis1 = h_a1->GetXaxis(); 
		axis1->SetLimits(0.7, 60.); 
	
		h_a1->Draw("E1X0"); 
	
		HistTools::SetStyle(h_ams, kPink-3, kFullCircle, 1.4, 1, 1); 
		HistTools::SetStyle(h_ace_BR, kRed, kFullCircle, 1.4, 1, 1);
		// HistTools::SetStyle(h_ams_BR_he, kBlue, kFullCircle, 1.4, 1, 1);
		HistTools::SetStyle(h_ams_BR_fake, kBlue, kFullCircle, 1.4, 1, 1);

		// fsp_he->Draw("SAME"); 
		rescaled_fit->Draw("SAME"); 
		h_ams->Draw("E1X0 SAME"); 
		h_ace_BR->Draw("E1X0 SAME"); 
		// h_ams_BR_he->Draw("E1X0 SAME"); 
		h_ams_BR_fake->Draw("E1X0 SAME"); 

		//fit->SetLineColor(kBlue); 
		//fit->Draw("SAME"); 
		legend->Draw("SAME");   

		c0->cd(2); 
		gPad->SetGrid(); 
		gPad->SetLogx(); 

		TAxis *axis2 = h_a2->GetXaxis(); 
		axis2->SetLimits(0.7, 60.); 

		////if (get_temp_int(element)<4) h_a2->GetYaxis()->SetRangeUser(0.5, 1.5); 
		h_a2->GetXaxis()->SetRangeUser(0.7, 60.); 

		h_a2->SetTitle("");
		h_a2->SetXTitle(Unit::GetEnergyLabel("GV"));
		h_a2->SetYTitle(Form("Estimated BR Fluxes / Rescaled %s Model", element)); 
		h_a2->SetTitleSize(0.3,"t"); 

		h_a2->Draw("E1X0"); 
		h_ratio->Draw("E1X0 SAME"); 
		h_ratio3->Draw("E1X0 SAME"); 

		if (iBR==0) { c0->Print(Form("./data/ACE/fill/fake_td_ams_v2/%s_flux_model.pdf(", element), "pdf"); } 
		if (iBR>0 && iBR<nBRs-1) { c0->Print(Form("./data/ACE/fill/fake_td_ams_v2/%s_flux_model.pdf", element), "pdf"); } 
		if (iBR==nBRs-1){
			c0->Print(Form("./data/ACE/fill/fake_td_ams_v2/%s_flux_model.pdf)", element), "pdf"); 
		} 

/*
		// plot AMS Integrated / Rescaled Combined Fit Template
		c1->cd(1); 
		gPad->SetGrid(); 
		gPad->SetLogx(); 

		TH1 *h_ratio3 = (TH1*) h_ams_new->Clone("h_ams_new"); 
		h_ratio3->Divide(fsp_comb); 

		TAxis *axis3 = h_ratio3->GetXaxis(); 
		axis3->SetLimits(0.7, 60.);

		HistTools::SetStyle(h_ratio3, kBlack, kFullCircle, 1.4, 1, 1);
		h_ratio3->SetXTitle(Unit::GetEnergyLabel("GV"));
		h_ratio3->SetTitle(Form("Template Test; ; AMS Integrated %s Flux / %s Template", element, element)); 
		h_ratio3->Draw("E1X0"); 
*/

		// plot the ratio of AMS_He_BR / <AMS_He>
		c1->cd(1); 
		gPad->SetGrid(); 
		gPad->SetLogx(); 

		////if (get_temp_int(element)<4) h_a3->GetYaxis()->SetRangeUser(0.5, 1.5); 
		h_a3->GetXaxis()->SetLimits(0.7, 60);		

		TLegend *l1 = new TLegend(0.62,0.8,0.9,0.9); 
		l1->AddEntry(h_ratio0, "AMS He(R,t)/<He(R)>", "PL"); 	

		h_a3->Draw("E1X0"); 
		h_a3->SetXTitle(Unit::GetEnergyLabel("GV"));
		h_a3->SetTitle(Form("; ; BR-%d He(R,t)/<He(R)>", iBR_true+2426)); 
		HistTools::SetStyle(h_ratio0, kBlue, kFullCircle, 1.4, 1, 1);
		h_ratio0->Draw("E1X0 SAME"); 

		// plot spectral indices
		TGraphAsymmErrors *gspind_ace = HistTools::GetSpectralIndex(h_ace_BR, 5, 2); // get spectral indices 
		TGraphAsymmErrors *gspind_ams = HistTools::GetSpectralIndex(h_ams_BR_fake, 4, 1); 

		c1->cd(2); 
		gPad->SetGrid(); 
		gPad->SetLogx();  

		TLegend *l_both = new TLegend(0.62,0.8,0.9,1.); 
		l_both->AddEntry(gspind_ace, "ACE Spectral Indices", "PL"); 
		l_both->AddEntry(gspind_ams, "AMS Spectral Indices", "PL"); 		
	
		////if (get_temp_int(element)<4) gspind_ace->GetYaxis()->SetRangeUser(-3, 3); 
		gspind_ace->GetXaxis()->SetLimits(0.7, 60);

		HistTools::SetStyle(gspind_ace, kPink, kFullCircle, 1.4, 1, 1); 
		gspind_ace->SetTitle(Form("; ; ACE %s Spectral Indices BR-%d", element, 2426+iBR_true));
		gspind_ace->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV"));
		gspind_ace->Draw("APL"); 

		//h_a4->Draw("E1X0"); 
		HistTools::SetStyle(gspind_ams, kBlue, kFullCircle, 1.4, 1, 1); 
		gspind_ams->SetTitle(Form("; ; AMS %s Spectral Indices", element)); 
		gspind_ams->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV"));
		gspind_ams->GetXaxis()->SetTitleSize(1.5); 
		gspind_ams->Draw("PL SAME"); 
		// fspind->Draw("SAME"); 
		l_both->Draw("SAME"); 

		if (iBR==0) { c1->Print(Form("./data/ACE/fill/fake_td_ams_v2/fake_ams_spind_%s.pdf(", element), "pdf"); }
		if (iBR>0 && iBR<nBRs-1) { c1->Print(Form("./data/ACE/fill/fake_td_ams_v2/fake_ams_spind_%s.pdf", element), "pdf"); }
		if (iBR==nBRs-1){
			c1->Print(Form("./data/ACE/fill/fake_td_ams_v2/fake_ams_spind_%s.pdf)", element), "pdf"); 
		} 

		// find spectral indices of h_ratio + h_ratio0 
		TGraphAsymmErrors *gspind_ace2 = HistTools::GetSpectralIndex(h_ratio1, 5, 2); // get spectral indices 
		TGraphAsymmErrors *gspind_ams2 = HistTools::GetSpectralIndex(h_ratio2, 4, 1); 

		HistTools::SetStyle(gspind_ace2, kPink, kFullCircle, 1.4, 1, 1); 
		HistTools::SetStyle(gspind_ams2, kBlue, kFullCircle, 1.4, 1, 1); 
		HistTools::SetStyle(fspind_he, kRed, kFullCircle, 0.8, 1, 1);

		// fit to F(R,t)/<F(R)>
		TF1 *fit_ratio = new TF1("f_ratio", "1+[0]*exp(-[1]*x)", 0.7, 60);
		TF1 *fit_ratio2 = new TF1("f_ratio2", "1+[0]*exp(-[1]*log(x))", 0.7, 60);

		double F1 = h_ratio1->GetBinContent(h_ratio1->FindBin(R1)), F2 = h_ratio0->GetBinContent(h_ratio0->GetNbinsX()); 
		double A1 = log(abs(F1-1)), A2 = log(abs(F2-1)); 
		TF1 *f1 = new TF1("f1","abs(sin(x)/x)*sqrt(x)",0,1);
		double rho = f1->GetRandom();
		double N = f1->GetRandom(); 

		last_pars = {N, rho}; 
		last_pars2 = {N, rho};

		for (int ipar=0; ipar<fit_ratio->GetNpar(); ++ipar){ 
			fit_ratio->SetParameter(ipar, last_pars[ipar]); 
			fit_ratio2->SetParameter(ipar, last_pars2[ipar]); 
		}   

		TObjArray data; 
		data.Add(h_ratio1);
		data.Add(h_ratio0);

		// fsp_he2->Print("V"); 

		// compute residuals w/ He Fit
		vector<double> rigmin, rigmax, chi2norm(2);
		ROOT::Fit::Fitter fitter;
		FitTools::SetCommonFitterOptions(fitter);
		FitTools::FitCombinedData(data, fsp_he2, "I", rigmin, rigmax, chi2norm, fitter, 3); 

		// fit AMS and ACE data at the same time
		vector<double> rigmin1, rigmax1, chi2norm1(2);
		ROOT::Fit::Fitter fitter1;
		FitTools::SetCommonFitterOptions(fitter1);
		FitTools::FitCombinedData(data, fit_ratio, "I", rigmin1, rigmax1, chi2norm1, fitter1, 3); 

		vector<double> rigmin2, rigmax2, chi2norm2(2);
		ROOT::Fit::Fitter fitter2;
		FitTools::SetCommonFitterOptions(fitter2);
		FitTools::FitCombinedData(data, fit_ratio2, "I", rigmin2, rigmax2, chi2norm2, fitter2, 3); 
 
		for (int ipar=0; ipar<fit_ratio->GetNpar(); ++ipar){
			if (abs(fit_ratio->GetParError(ipar)/fit_ratio->GetParameter(ipar))<1){ 
				g_pars[ipar]->SetPoint(iBR, UBRToTime(iBR_true+2426), fit_ratio->GetParameter(ipar)); 
				g_pars[ipar]->SetPointError(iBR, 0, fit_ratio->GetParError(ipar)); 
			} else {
				TLatex *ltx_pars = new TLatex(g_pars[ipar]->GetX()[iBR], g_pars[ipar]->GetY()[iBR],"annotation"); 
				g_pars[ipar]->GetListOfFunctions()->Add(ltx_pars); 

				g_pars[ipar]->SetPoint(iBR, UBRToTime(iBR_true+2426), 0); 
				g_pars[ipar]->SetPointError(iBR, 0, 0); 
			}
		}  

		TF1 *fspind_ratio = HistTools::GetSpectralIndex(fit_ratio, "fspind_ratio", R1, R2); 
		TF1 *fspind_ratio2 = HistTools::GetSpectralIndex(fit_ratio2, "fspind_ratio2", R1, R2); 

		c2->cd(1); 
		gPad->SetGrid(); 
		gPad->SetLogx();
		gPad->SetBottomMargin(0.01); 

		////if (get_temp_int(element)<4) h_a4->GetYaxis()->SetRangeUser(0.5, 2.2); 
		h_a4->GetXaxis()->SetLimits(0.7, 60); 
	
		TLegend *l2 = new TLegend(0.62,0.8,0.9,1.); 
		l2->AddEntry(h_ratio1, Form("ACE %s(R,t)/<%s(R)>", element, element), "PL");
		l2->AddEntry(h_ratio0, "AMS He(R,t)/<He(R)>", "PL");  
		l2->AddEntry(fsp_he2, "He(R,t)Fit/<He(R)>Fit", "L"); 

		h_a4->Draw("E1X0"); 
		h_a4->SetXTitle(Unit::GetEnergyLabel("GV"));
		h_a4->SetTitle(Form("; ; BR-%d %s(R,t)/<%s(R)>", iBR_true+2426, element, element));
		// h_a4->SetTitleSize(0.5,"y"); 
		HistTools::SetStyle(h_ratio1, kPink, kFullCircle, 1.4, 1, 1);
		HistTools::SetStyle(h_ratio0, kBlue, kFullCircle, 1.4, 1, 1);
		//fit_ratio2->SetLineColor(kGreen); 
		//fit_ratio->Draw("SAME");
		//fit_ratio2->Draw("SAME");
		fsp_he2->Draw("SAME"); 
		h_ratio1->Draw("E1X0 SAME");
		h_ratio0->Draw("E1X0 SAME"); 
		l2->Draw("SAME"); 

		TH1D *h_fitres[2]; 
		TH1D *h_fitres1[2];
		TH1D *h_fitres2[2]; 
		TLegend *l_chi2 = new TLegend(0.62,0.8,0.9,1.0); 

		for (UShort_t i = 0; i < 2; ++i)
		{
  			TH1 *hist = HistTools::ToHist(data[i]);
			TH1 *h_fiterr = (TH1*) HistTools::GetFitError(hist, fsp_he2, "_fiterr", false, false, true, 10, DBL_MAX, 0.68, &fitter);
			//h_fiterr->Print("range"); 
   			h_fitres[i] = (TH1D *)HistTools::GetResiduals(hist, fsp_he2, "_fitres", false, true, true, 5, 1, 0.68, &fitter);
			//h_fitres[i]->Print("range"); 
			HistTools::CopyStyle(hist, h_fitres[i]);

   			UShort_t ndf  = hist->GetNbinsX();
   			Double_t chi2 = chi2norm[i]*ndf; 
 
			l_chi2->AddEntry(h_fitres[i], Form("He(R,t)Fit/<He(R)>Fit chi2/ndf=%6.2f/%-2u", chi2, ndf), "P");  

			if (i==0){ 
				g_chi2_ace->SetPoint(iBR, UBRToTime(iBR_true+2426), chi2/ndf);
				for (int bin=0; bin <= h_fitres[i]->GetNbinsX(); ++bin){ 
					if (bin%2==0){ 
						g_ratio_ace[bin/2]->SetPoint(iBR, UBRToTime(iBR_true+2426), h_ratio1->GetBinContent(bin+1)); 
						g_ratio_ace[bin/2]->SetPointError(iBR, 0, h_ratio1->GetBinError(bin+1)); 
						g_ratio_ace_fit[bin/2]->SetPoint(iBR, UBRToTime(iBR_true+2426), fsp_he2->Eval(h_ratio1->GetBinCenter(bin+1))); 
						g_ratio_ace_fit[bin/2]->SetPointError(iBR, 0, h_fiterr->GetBinError(bin+1)); 
						g_residual_ace[bin/2]->SetPoint(iBR, UBRToTime(iBR_true+2426), h_fitres[i]->GetBinContent(bin+1)); 
						g_residual_ace[bin/2]->SetPointError(iBR, 0, h_fitres[i]->GetBinError(bin+1)); 

						double df = h_fitres[i]->GetBinContent(bin+1)*hist->GetBinContent(bin+1)/hist->GetBinError(bin+1); 
						double ddf = df*sqrt(pow(h_fitres[i]->GetBinError(bin+1)/h_fitres[i]->GetBinContent(bin+1))+pow(hist->GetBinError(bin+1)/hist->GetBinContent(bin+1))); // d(data-fit) 
						//printf(" df = %0.4f, a= %0.4f, b=%0.4f \n", df, h_fitres[i]->GetBinContent(bin+1)*h_fiterr->GetBinContent(bin+1), h_fiterr->GetBinError(bin+1)); 

						g_resierr_ace[bin/2]->SetPoint(iBR, UBRToTime(iBR_true+2426), df); 
						g_resierr_ace[bin/2]->SetPointError(iBR, 0, ddf); 

						// printf("bin/2 = %d, iBR = %d, time = %10.4u, fitres = %10.4f \n", bin/2, iBR, UBRToTime(iBR_true+2426), h_fitres[0]->GetBinContent(bin+1)); 
					} 
				}   
			} 

			if (i==1){ 
				g_chi2_ams->SetPoint(iBR, UBRToTime(iBR_true+2426), chi2/ndf);
				for (int bin=0; bin <= h_ams_new->GetNbinsX()-1; ++bin){ 
					g_ratio_ams[bin]->SetPoint(iBR, UBRToTime(iBR_true+2426), h_ratio0->GetBinContent(bin+1));
					g_ratio_ams[bin]->SetPointError(iBR, 0, h_ratio0->GetBinError(bin+1)); 
					g_ratio_ams_fit[bin]->SetPoint(iBR, UBRToTime(iBR_true+2426), fsp_he2->Eval(h_ratio0->GetBinCenter(bin+1))); 
					g_ratio_ams_fit[bin]->SetPointError(iBR, 0, h_fiterr->GetBinError(bin+1)); 
					g_residual_ams[bin]->SetPoint(iBR, UBRToTime(iBR_true+2426), h_fitres[i]->GetBinContent(bin+1));
					g_residual_ams[bin]->SetPointError(iBR, 0, h_fitres[i]->GetBinError(bin+1)); 
					// printf("bin = %d, iBR = %d, time = %10.4u, fitres = %10.4f \n", bin, iBR, UBRToTime(iBR_true+2426), h_fitres[1]->GetBinContent(bin+1)); 
				}  
			} 	

		} 

		for (UShort_t i = 0; i < 2; ++i)
		{
  			TH1 *hist = HistTools::ToHist(data[i]);
   			h_fitres1[i] = (TH1D *)HistTools::GetResiduals(hist, fit_ratio, "_fitres1", false, true, true, 5, 1, 0.68, &fitter1);
			HistTools::CopyStyle(hist, h_fitres1[i]);

   			UShort_t ndf  = hist->GetNbinsX();
   			Double_t chi2 = chi2norm1[i]*ndf; 
 
			// l_chi2->AddEntry(h_fitres1[i], Form("1+[0]*exp(-[1]*R) chi2/ndf=%6.2f/%-2u", chi2, ndf), "P");

			if (i==0){ 
				g_chi2_ace2->SetPoint(iBR, UBRToTime(iBR_true+2426), chi2/ndf);
			} 

			if (i==1){ 
				g_chi2_ams2->SetPoint(iBR, UBRToTime(iBR_true+2426), chi2/ndf);
			}   

		}	

		for (UShort_t i = 0; i < 2; ++i)
		{
  			TH1 *hist = HistTools::ToHist(data[i]);
   			h_fitres2[i] = (TH1D *)HistTools::GetResiduals(hist, fit_ratio2, "_fitres2", false, true, true, 5, 1, 0.68, &fitter2);
			HistTools::CopyStyle(hist, h_fitres2[i]);

   			UShort_t ndf  = hist->GetNbinsX();
   			Double_t chi2 = chi2norm2[i]*ndf; 
   	
			// l_chi2->AddEntry(h_fitres2[i], Form("1+[0]*exp(-[1]*ln(R)) chi2/ndf=%6.2f/%-2u", chi2, ndf), "P"); 	
		}

		// %time %par0 %par0_er %par1 %par1_er  
		printf("%d %0.6f %0.6f %0.6f %0.6f \n", UBRToTime(iBR_true+2426), fit_ratio->GetParameter(0), fit_ratio->GetParError(0), fit_ratio->GetParameter(1), fit_ratio->GetParError(1)); 

		TH1 *h_comp1 = (TH1*) h_ace_BR->Clone("h_comp1");
		h_comp1->Divide(rescaled_fit); 

		double comp1_sum=0, dcomp1_sum=0; // compute average ratio of data/template-1 of ACE F(R,t)/<F(R)> manually 		

		for(int nbin=0;nbin<h_ace_rig_ave->GetNbinsX();++nbin){
			comp1_sum += h_comp1->GetBinContent(nbin); 
			dcomp1_sum += h_comp1->GetBinError(nbin)*h_comp1->GetBinError(nbin); 
			//printf("comp1_sum = %0.6f \n", comp1_sum);
		}
		double comp1_ave = comp1_sum/7; 
		double dcomp1_ave = sqrt(dcomp1_sum)/7; 

		g_ratio_ave->SetPoint(iBR, UBRToTime(iBR_true+2426), comp1_ave-1); 
		g_ratio_ave->SetPointError(iBR, 0, dcomp1_ave); 

		c2->cd(2);
		gPad->SetGrid();
		gPad->SetLogx(); 
		gPad->SetTopMargin(0.01); 

		////if (get_temp_int(element)<4) h_a5->GetYaxis()->SetRangeUser(-0.15, 0.15)  ; 
		h_a5->GetXaxis()->SetLimits(0.7, 60); 

		h_a5->Draw("E1X0"); 
		h_a5->SetXTitle(Unit::GetEnergyLabel("GV"));
		h_a5->SetTitle("; ; Data/He_Fit_Ratio-1"); 

		h_fitres[0]->Draw("E1X0 SAME");
		h_fitres[1]->Draw("E1X0 SAME");

		l_chi2->Draw("SAME"); 		

		c2->cd(3);
		gPad->SetGrid();
		gPad->SetLogx(); 

		TLegend *l_both3 = new TLegend(0.62,0.7,0.8,1.0); 
		l_both3->AddEntry(gspind_ace2, Form("ACE %s(R,t)/<%s(R)> Spectral Indices", element, element), "PL");
		l_both3->AddEntry(gspind_ams2, Form("AMS %s(R,t)/<%s(R)> Spectral Indices", element, element), "PL"); 
		// l_both3->AddEntry(fspind_he, "AMS He(R,t)_Fit/temp<He(R)>_Fit Spectral Indices", "L");  	
		// l_both3->AddEntry(fspind_ratio, "Spectral Index of 1+[0]*exp(-[1]*R) Fit to He(R,t)/<He(R)>", "L");	
		l_both3->AddEntry(fspind_he, "AMS He(R,t)_Fit/<He(R)>_Fit Spectral Indices", "L");  
		// l_both3->AddEntry(fspind_ratio2, "Spectral Index of 1+[0]*exp(-[1]*ln(R)) Fit to He(R,t)/<He(R)>", "L"); 

		//if (get_temp_int(element)<4) gspind_ace2->GetYaxis()->SetRangeUser(-2, 1.3); 
		TAxis *axis5 = gspind_ace2->GetXaxis(); 
		axis5->SetLimits(0.7, 60.); 

		gspind_ace2->SetTitle(Form("; ; Spectral Indices Model of %s F(R,t)/<F(R)> at BR-%d", element, 2426+iBR_true));
		gspind_ace2->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV"));
		gspind_ace2->Draw("APL");  
		//fspind_ratio->Draw("SAME"); 
		//fspind_ratio2->Draw("SAME"); 
		//fspind_ratio2->SetLineColor(kGreen); 
		fspind_he->Draw("SAME");
		gspind_ams2->Draw("PL SAME"); 
		gspind_ace2->Draw("PL SAME");
		l_both3->Draw("SAME"); 

		// break; 
		rescaled_fit->Write(Form("rescaled_fit_BR%d", 2426+iBR_true)); 
		fit_he->Write(Form("fit_he_BR%d", 2426+iBR_true)); 
		//fsp_he2->Write(Form("fsp_he2_BR%d", 2426+iBR_true)); 
		fit_he_ave->Write(Form("fit_he_ave_BR%d", 2426+iBR_true)); 
		fit_ratio->Write(Form("fit_ratio_BR%d", 2426+iBR_true)); 
		h_ratio1->Write(Form("h_ratio1_BR%d", 2426+iBR_true));
		h_ratio0->Write(Form("h_ratio0_BR%d", 2426+iBR_true)); 
		h_fitres[0]->Write(Form("h_fitres_ace_BR%d", 2426+iBR_true));
		h_fitres[1]->Write(Form("h_fitres_ams_BR%d", 2426+iBR_true)); 

		if (iBR==0) c2->Print(Form("./data/ACE/fill/fake_td_ams_v2/spind_model_%s_temp.pdf(", element), "pdf"); 
		if (iBR>0 && iBR<nBRs-1) c2->Print(Form("./data/ACE/fill/fake_td_ams_v2/spind_model_%s_temp.pdf", element), "pdf"); 
		if (iBR==nBRs-1){
			//gspind_ace2->Print("range");
			cout << " " << endl; 
			//gspind_ams2->Print("range");
			cout << " " << endl; 
			// h_ratio2->Print("range");
			c2->Print(Form("./data/ACE/fill/fake_td_ams_v2/spind_model_%s_temp.pdf)", element), "pdf"); 
		}

	   	if (iBR+2426==2472-1) iBR_true += 3; 
		else iBR_true ++; 
	}
			
	c3->cd(1);
	gPad->SetGrid(); 
	HistTools::SetStyle(g_chi2_ace, kPink, kFullCircle, 1.4, 1, 1); 
	HistTools::SetStyle(g_chi2_ace2, kBlue, kFullCircle, 1.4, 1, 1); 
	TLegend *l_c3_1 = new TLegend(0.62,0.8,0.9,0.9); 
	l_c3_1->AddEntry(g_chi2_ace, "He(R,t)_Fit/<He(R)>_Fit", "PL"); 
	l_c3_1->AddEntry(g_chi2_ace2, "1+[0]*exp(-[1]*R)", "PL"); 
	g_chi2_ace->SetTitle(Form("; ; ACE %s(R,t)/<%s(R)> Fitting Chi2/NDF", element, element)); 
	g_chi2_ace->GetXaxis()->SetTimeDisplay(1);
	g_chi2_ace->GetXaxis()->SetTimeFormat("%m-%y"); 
	g_chi2_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	g_chi2_ace->GetXaxis()->SetTitleSize(0.7); 
	// //if (get_temp_int(element)<4) g_chi2_ace->GetYaxis()->SetRangeUser(0., 4.5); 
	g_chi2_ace->Draw("APL"); 
	g_chi2_ace2->Draw("PL SAME"); 
	l_c3_1->Draw("SAME");

	c3->cd(2); 
	gPad->SetGrid(); 
	HistTools::SetStyle(g_chi2_ams, kPink, kFullCircle, 1.4, 1, 1); 
	HistTools::SetStyle(g_chi2_ams2, kBlue, kFullCircle, 1.4, 1, 1); 
	TLegend *l_c3_2 = new TLegend(0.62,0.8,0.9,0.9); 
	l_c3_2->AddEntry(g_chi2_ams, "He(R,t)_Fit/<He(R)>_Fit", "PL");
	l_c3_2->AddEntry(g_chi2_ams2, "1+[0]*exp(-[1]*R)", "PL");
	g_chi2_ams->SetTitle(Form("; ; AMS %s(R,t)/<%s(R)> Fitting Chi2/NDF", element, element)); 
	g_chi2_ams->GetXaxis()->SetTimeDisplay(1);
	g_chi2_ams->GetXaxis()->SetTimeFormat("%m-%y");
	g_chi2_ams->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	g_chi2_ams->GetXaxis()->SetTitleSize(0.7);
	// //if (get_temp_int(element)<4) g_chi2_ams->GetYaxis()->SetRangeUser(0., 2.2); 
	g_chi2_ams->Draw("APL");
	g_chi2_ams2->Draw("PL SAME");  
	l_c3_2->Draw("SAME"); 

	c3->Print(Form("./data/ACE/fill/fake_td_ams_v2/chi2_vs_BR_%s.png", element)); 

	g_ratio_ave_norm = get_norm_graph( g_ratio_ave ); 

	for (int iBR=0; iBR<nBRs; ++iBR){
				Double_t x=0, y=0; 
				g_ratio_ave_norm->GetPoint(iBR, x, y); 
				g_ratio_ave_norm->SetPoint(iBR, x, y-1); 			 
	}	

	TGraphErrors *g_ave_residual_ace_rig = new TGraphErrors(); // time-averaged residual vs. rigidity
	TGraphErrors *g_ave_residual_ams_rig = new TGraphErrors(); // time-averaged residual vs. rigidity 

	for (int bin=0; bin <= h_ace_rig_ave->GetNbinsX(); ++bin){
 
		if (bin%2==0) {  

			g_ratio_time_norm[bin/2] = get_norm_graph( g_ratio_time[bin/2] ); 
		
			for (int iBR=0; iBR<g_ratio_time_norm[bin/2]->GetN(); ++iBR){
				Double_t x=0, y=0; 
				g_ratio_time_norm[bin/2]->GetPoint(iBR, x, y); 
				g_ratio_time_norm[bin/2]->SetPoint(iBR, x, y-1); 			 
			}

			g_ave_residual_ace_rig->SetPoint(bin/2, h_ace_rig_ave->GetBinCenter(bin+1), g_residual_ace[bin/2]->GetMean(2) ); 
			g_ave_residual_ace_rig->SetPointError(bin/2, 0, g_residual_ace[bin/2]->GetRMS(2) ); 

			c4->cd(1);
			gPad->SetGrid(); 

			TLegend *l_comp1 = new TLegend(0.62, 0.8, 0.9, 1.0); 
			l_comp1->AddEntry(g_ratio_ave_norm, "Normalized Averaged ACE Data/Template-1", "PL"); 
			l_comp1->AddEntry(g_residual_ace[bin/2], Form("ACE He_Fit Residuals (%0.4f GV)", h_ace_rig_ave->GetBinCenter(bin+1)), "PL"); 

			g_residual_ace[bin/2]->SetTitle(Form(" ; ; ACE %s Data/Template-1 vs. Fitting Residuals", element));  
			g_residual_ace[bin/2]->GetXaxis()->SetTimeDisplay(1);
			g_residual_ace[bin/2]->GetXaxis()->SetTimeFormat("%m-%y");
			g_residual_ace[bin/2]->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g_residual_ace[bin/2]->GetXaxis()->SetTitleSize(0.7);
			//if (get_temp_int(element)<4) g_residual_ace[bin/2]->GetYaxis()->SetRangeUser(-0.5, 0.5);
			g_residual_ace[bin/2]->Draw("APL");
			// g_ratio_ave_norm->Draw("PL SAME"); 
			g_ratio_time_norm[bin/2]->Draw("PL SAME"); 
			l_comp1->Draw("SAME"); 		

			if (bin==0) c4->Print(Form("./data/ACE/fill/fake_td_ams_v2/fitres_ace_vs_BR_%s.pdf(", element), "pdf"); 
			if (bin>0 && bin<h_ace_rig_ave->GetNbinsX()-1) c4->Print(Form("./data/ACE/fill/fake_td_ams_v2/fitres_ace_vs_BR_%s.pdf", element), "pdf");  
			if (bin==h_ace_rig_ave->GetNbinsX()-1){
				c4->Print(Form("./data/ACE/fill/fake_td_ams_v2/fitres_ace_vs_BR_%s.pdf)", element), "pdf"); 
			} 

			g_ratio_ace[bin/2]->Write(Form("g_ratio_ace_%d", bin/2)); 
			g_ratio_ace_fit[bin/2]->Write(Form("g_ratio_ace_fit_%d", bin/2)); 
			g_residual_ace[bin/2]->Write(Form("g_residual_ace_%d", bin/2)); 
			g_resierr_ace[bin/2]->Write(Form("g_resierr_ace_%d", bin/2)); 
		} 
	} 

	for (int bin=0; bin <= h_ams_new->GetNbinsX()-1; ++bin){ 

			g_ave_residual_ams_rig->SetPoint(bin, h_ams_new->GetBinCenter(bin+1), g_residual_ams[bin]->GetMean(2) ); 
			g_ave_residual_ams_rig->SetPointError(bin, 0, g_residual_ams[bin]->GetRMS(2) ); 

			c4_2->cd(1); 
			gPad->SetGrid(); 

			g_residual_ams[bin]->SetTitle(Form("AMS %s(R,t)/<%s(R)>; ; Fitting Residuals (%0.4f GV)", "He_Fit", "He_Fit", h_ams_new->GetBinCenter(bin+1)));  
			g_residual_ams[bin]->GetXaxis()->SetTimeDisplay(1);
			g_residual_ams[bin]->GetXaxis()->SetTimeFormat("%m-%y");
			g_residual_ams[bin]->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g_residual_ams[bin]->GetXaxis()->SetTitleSize(0.7);
			//if (get_temp_int(element)<4) g_residual_ams[bin]->GetYaxis()->SetRangeUser(-0.1, 0.05);  
			// g_residual_ams[bin]->Print("range"); 
			g_residual_ams[bin]->Draw("APL"); 

			if (bin==0) c4_2->Print(Form("./data/ACE/fill/fake_td_ams_v2/fitres_ams_vs_BR_%s.pdf(", "He_fit"), "pdf"); 
			if (bin>0 && bin<h_ams_new->GetNbinsX()-1) c4_2->Print(Form("./data/ACE/fill/fake_td_ams_v2/fitres_ams_vs_BR_%s.pdf", "He_fit"), "pdf");  
			if (bin==h_ams_new->GetNbinsX()-1){
				c4_2->Print(Form("./data/ACE/fill/fake_td_ams_v2/fitres_ams_vs_BR_%s.pdf)", "He_fit"), "pdf"); 
			} 

			g_ratio_ams[bin]->Write(Form("g_ratio_ams_%d", bin)); 
			g_ratio_ams_fit[bin]->Write(Form("g_ratio_ams_fit_%d", bin)); 
			// g_ratio_ams_fit[bin]->Print("range"); 
			g_residual_ams[bin]->Write(Form("g_residual_ams_%d", bin)); 
	} 

	g_chi2_ace->Write("g_chi2_ace"); 
	g_chi2_ams->Write("g_chi2_ams"); 

	g_ave_residual_ace_rig->Write("g_ave_residual_ace_rig"); 
	g_ave_residual_ams_rig->Write("g_ave_residual_ams_rig"); 

	file->Close(); 

	gROOT->ProcessLine(".>"); 

	g_chi2_ace->Print("range");
	g_chi2_ams->Print("range"); 

	// TGraphErrors *g_mean_residual_ace = ave_grapherrors( g_residual_ace, 7 ); 
	// TGraphErrors *g_mean_residual_ams = ave_grapherrors( g_residual_ams, h_ams_new->GetNbinsX()-1 ); 

	// g_ave_residual_ace_rig->Print("range");
	// g_ave_residual_ams_rig->Print("range"); 

	for (int ipar=0; ipar<2; ++ipar){
		c5->cd(ipar+1); 
		gPad->SetGrid(); 
		g_pars[ipar]->Draw("APL"); 
		// g_pars[ipar]->Print("range"); 
		HistTools::SetStyle(g_pars[ipar], kPink, kFullCircle, 1.4, 1, 1); 
		g_pars[ipar]->GetXaxis()->SetTimeDisplay(1);
		g_pars[ipar]->GetXaxis()->SetTimeFormat("%m-%y");
		g_pars[ipar]->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
		g_pars[ipar]->GetXaxis()->SetTitleSize(0.7);
		//if (get_temp_int(element)<4) g_pars[ipar]->GetYaxis()->SetRangeUser(-3., 3.); 
		if (ipar==0) g_pars[ipar]->SetTitle(Form("1+[0]*exp(-[1]*R); ; Parameter [%d]", ipar));
		else g_pars[ipar]->SetTitle(Form(" ; ; Parameter [%d]", ipar));
	}		

	TGraphErrors *g_pars_ratio = new TGraphErrors(); 
	for (int i=0; i<g_pars[0]->GetN(); ++i){
		double ratio_par = g_pars[1]->GetY()[i]/g_pars[0]->GetY()[i]; 
		double ratio_par_err = g_pars[1]->GetY()[i]/g_pars[0]->GetY()[i]*sqrt(pow(g_pars[1]->GetEY()[i]/g_pars[1]->GetY()[i], 2)+pow(g_pars[0]->GetEY()[i]/g_pars[0]->GetY()[i], 2));  
		//if (abs(ratio_par_err/ratio_par) < 1){ 
		g_pars_ratio->SetPoint(i, g_pars[0]->GetX()[i], g_pars[1]->GetY()[i]/g_pars[0]->GetY()[i]); 
		g_pars_ratio->SetPointError(i, 0, g_pars[1]->GetY()[i]/g_pars[0]->GetY()[i]*sqrt(pow(g_pars[1]->GetEY()[i]/g_pars[1]->GetY()[i], 2)+pow(g_pars[0]->GetEY()[i]/g_pars[0]->GetY()[i], 2))); 
		//} else {
		//	TLatex *ltx_pars = new TLatex(g_pars_ratio->GetX()[i], g_pars_ratio->GetY()[i],"outlier"); 
		//	g_pars_ratio->GetListOfFunctions()->Add(ltx_pars);
		//	g_pars_ratio->SetPoint(i, g_pars[0]->GetX()[i], 0); 
		//	g_pars_ratio->SetPointError(i, 0, 0); 
		//}
		// printf ("par1 = %0.4f, par1err = %0.4f, par0 = %0.4f, par0err = %0.4f \n", g_pars[1]->GetY()[i], g_pars[1]->GetEY()[i], g_pars[0]->GetY()[i], g_pars[0]->GetEY()[i]); 
	}

	c5->cd(3);
	gPad->SetGrid();
	HistTools::SetStyle(g_pars_ratio, kPink, kFullCircle, 1.4, 1, 1);
	g_pars_ratio->GetXaxis()->SetTimeDisplay(1);
	g_pars_ratio->GetXaxis()->SetTimeFormat("%m-%y");
	g_pars_ratio->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	g_pars_ratio->GetXaxis()->SetTitleSize(0.7);
	//if (get_temp_int(element)<4) g_pars_ratio->GetYaxis()->SetRangeUser(-15., 30.); 
	g_pars_ratio->SetTitle(Form(" ; ; [%d]/[%d]", 1, 0)); 
	g_pars_ratio->Draw("APL"); 

	g_pars_ratio->Print("range"); 

	c5->Print(Form("./data/ACE/fill/fake_td_ams_v2/pars_vs_BR_%s.png", element)); 

	TCanvas *c6 = new TCanvas("c6", "", 1600, 900); // plot time-averaged fit residuals to F(R,t)/<F(R)> vs. rigidity 
	c6->Divide(1, 1);  

	c6->cd(1); 
	gPad->SetGrid(); 
	gPad->SetLogx(); 

	HistTools::CopyStyle( h_ace_rig_ave, g_ave_residual_ace_rig ); 
	HistTools::CopyStyle( h_ams_new, g_ave_residual_ams_rig );
	//if (get_temp_int(element)<4) g_ave_residual_ace_rig->GetYaxis()->SetRangeUser(-0.085, 0.085); 
	g_ave_residual_ace_rig->GetXaxis()->SetLimits(0.7, 60);
	g_ave_residual_ace_rig->SetTitle(Form(" ; ; Time-averaged Fit Residuals of He(R,t)_Fit/<He(R)>_Fit vs. ACE %s Flux", element)); 
	g_ave_residual_ace_rig->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV")); 

	TLegend *l_both4 = new TLegend(0.62,0.8,0.9,0.9); 
	l_both4->AddEntry(g_ave_residual_ace_rig, "ACE", "PL");
	l_both4->AddEntry(g_ave_residual_ams_rig, "AMS", "PL"); 

	g_ave_residual_ace_rig->Draw("APL");
	g_ave_residual_ams_rig->Draw("PLSAME"); 
	l_both4->Draw("SAME"); 

	c6->Print(Form("./data/ACE/fill/fake_td_ams_v2/time_ave_residual_vs_rig_%s.png", element)); 

	return 0; 
}

// F4  
void ace_fake_td_ams_v3(const char *element, Particle::Type isotope){

	gStyle->SetOptStat(0); 
	// gStyle->SetTitleSize(26,"t");

	gSystem->mkdir("data/ACE/fill/fake_td_ams_v3", true);	

 	Experiments::DataPath = "data";	
   	int data_value[n_ele] = {0, 18, 23, 27, 31, 20, 43, 22}; // p, He, Li, Be, B, C, N, O

	TGraphErrors *par0_X[4];
	TGraphErrors *par1_X[4]; 

	for (int i=0; i<4; ++i){ 

		par0_X[i] = new TGraphErrors(Form("data/ACE/fill/fake_td_ams_v2/fit_pars_%s.dat", ACE_Element[i]), "%lg %lg %lg %*s %*s"); 
		par1_X[i] = new TGraphErrors(Form("data/ACE/fill/fake_td_ams_v2/fit_pars_%s.dat", ACE_Element[i]), "%lg %*s %*s %lg %lg"); 

	}

	double *kin_bins = get_kin_bins(element);
        double *SpallCorr = get_spall_corr(element);
	double *SpallCorrUnc = get_spall_corr_unc(element);
	double *EMed = get_EMed(element);

	const UInt_t FirstACEBR = 2240;
	vector<UInt_t> BRs;
	// we stop at BR 2493, which ends on 2016/05/23, just 3 days before the end of the data taking period for AMS nuclei
	for (UInt_t br=2426; br<=2493; ++br) { 
		if (br != 2472 && br != 2473) BRs.push_back(br-FirstACEBR); 
	}

	// load AMS He BR fluxes 
	const int nBRs = Experiments::Info[Experiments::AMS02].Dataset[1].nMeasurements; // read the number of BRs
	TH1D **h_BR_he = Experiments::GetDatasetHistograms(Experiments::AMS02, 4);
	// Create a new set of histograms which match the AMS monthly bins with the AMS integrated bins 

	TH1 *h_ams = Experiments::GetMeasurementHistogram(Experiments::AMS02, get_ams_data_value(element), 0); // load AMS data for a given template 
	TH1 *h_ams_new = (TH1D*) h_BR_he[0]->Clone("h_ams_new"); 
	// TH1 *h_he_int = (TH1*) ave_hist( h_BR_he, nBRs);  
	TH1 *h_he = Experiments::GetMeasurementHistogram(Experiments::AMS02, 18, 0);
	TH1 *h_he_int = (TH1D*) h_BR_he[0]->Clone("h_ams_new");  

	// create h_ams_new and h_he_int to match the bins  

	for (int bin=1; bin <= h_ams_new->GetNbinsX(); ++bin){ 
	   if (!strcmp(element, "B") || !strcmp(element, "C")){ 
		h_ams_new->SetBinContent(bin, h_ams->GetBinContent(bin)); 
		h_ams_new->SetBinError(bin, h_ams->GetBinError(bin));
	   }
	   if (!strcmp(element, "N") || !strcmp(element, "O")){ 
		h_ams_new->SetBinContent(bin, h_ams->GetBinContent(bin-1)); 
		h_ams_new->SetBinError(bin, h_ams->GetBinError(bin-1)); 
	   }
	}
	if (!strcmp(element, "N") || !strcmp(element, "O")){
		h_ams_new->SetBinContent(1, 0); 
		h_ams_new->SetBinError(1, 0);
	}
 
	for (int bin=1; bin<=h_BR_he[0]->GetNbinsX(); ++bin) { 
		h_he_int->SetBinContent(bin, h_he->GetBinContent(bin)); 
		h_he_int->SetBinError(bin, h_he->GetBinError(bin)); 
	}

	// h_ams->Print("range"); 
	// h_ams_new->Print("range"); 

	// load ACE average 
	TH1 *h_ace_ene_ave = HistTools::GraphToHist(get_ace_average_graph( element , &BRs[0], BRs.size() ), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
	TH1 *h_ace_rig_ave = HistTools::TransformEnergyAndDifferentialFluxNew(h_ace_ene_ave, isotope, "MeV/n cm", "GV m", "_rig"); // load averaged ACE data for the same element, converted in rigidity
	HistTools::SetStyle(h_ams_new, kBlue, kFullCircle, 1.4, 1, 1);
	HistTools::SetStyle(h_ace_rig_ave, kBlack, kFullCircle, 1.4, 1, 1);

	UShort_t namsbins = h_ams_new->GetNbinsX(); 
	UShort_t nacebins = h_ace_rig_ave->GetNbinsX(); 
	double R1 = h_ace_rig_ave->GetBinLowEdge(1);
	double R2 = h_ams->GetBinLowEdge(namsbins+1);

	int nnodes = 9; // combined template 
	int nnodes_ams = 7; 

	printf(" template = %s, %d \n", element, get_ams_data_value(element)); 

	// load combined fit template 
	TFile *file0 = new TFile(Form("data/ACE/compare/fit_%s_%dnodes.root", get_template(element), nnodes));  

	Spline *sp_comb = new Spline("sp_comb", nnodes, Spline::LogLog | Spline::PowerLaw);
	TF1 *fsp_comb = sp_comb->GetTF1Pointer();  
	TF1 *fit_comb = (TF1*) file0->Get("fit_both"); 

	HistTools::CopyParameters(fit_comb, fsp_comb); // error 
	double x1, x2;
	fit_comb->GetRange(x1,x2);
	fsp_comb->SetRange(x1,x2);

	// load AMS He combined, modified spline  
	TFile *file1 = new TFile(Form("data/amsfit/fit_result_node%d.root", nnodes_ams));

	// load ACE BR 
	TFile *file2 = new TFile(Form("data/ACE/fill/%s_fill.root", element));
	TFile *file2a = new TFile(Form("data/ACE/fill/%s_fill2.root", element));  

	TGraphErrors *g_ratio_time[7]; 
	TGraphErrors *g_ratio_time_norm[7]; 

	for (int i=0; i<7; ++i){
		g_ratio_time[i] = (TGraphErrors*) file2a->Get(Form("g_ratio_time_%d", i)); 	
		HistTools::SetStyle(g_ratio_time[i], kBlue, kFullCircle, 1.4, 1, 1);
	}

	// load AMS He Model by Rescaled C Template
	TFile *file3 = new TFile(Form("data/ACE/extend2/fit_C_temp_he_%dnodes.root", nnodes));

	TCanvas *c0 = new TCanvas("c0", "", 1600, 900); // fluxes & flux/template
	c0->Divide(1, 2); 

	TCanvas *c1 = new TCanvas("c1", "", 1600, 900); // He ratio, F(R,t)/<F(R)> & spectral indices 
	c1->Divide(1, 2); 
	
	TCanvas *c2 = new TCanvas("c2", "", 1600, 900); // Modeling F(R,t)/<F(R)> 
							// fit/data of F(R,t)/<F(R)> 
	c2->Divide(1, 3); 				// spectral indices of F(R,t)/<F(R)>, here is defined by h_ratio (ACE) and h_ratio0 (AMS)

	TCanvas *c3 = new TCanvas("c3", "", 1600, 900); // plot chi2/ndf vs. BR for both ACE & AMS F(R,t)/<F(R)> fits  
	c3->Divide(2, 1);  

	TCanvas *c4 = new TCanvas("c4", "", 1600, 900); // plot for F(R,t)/<F(R)> fit residuals for each ACE bin 
	c4->Divide(1, 1); 

	TCanvas *c4_2 = new TCanvas("c4_2", "", 1600, 900); // c4 but for each AMS bin 
	c4_2->Divide(1, 1); 
	
	TCanvas *c5 = new TCanvas("c5", "", 1600, 900); // plot parameter [0] [1] vs. BR for the fit 
	c5->Divide(1, 3); 

	TH1 *h_a1 = new TH1D("", "", 3000, 0, 3000); 
	TH1 *h_a2 = new TH1D("", "", 3000, 0, 3000);
	TH1 *h_a3 = new TH1D("", "", 3000, 0, 3000);
	TH1 *h_a4 = new TH1D("", "", 3000, 0, 3000);	
	TH1 *h_a5 = new TH1D("", "", 3000, 0, 3000); 

	vector<double> last_pars;
	vector<double> last_pars2 = {0.5, 0.5};  

	TGraph *g_chi2_ace = new TGraph(); 
	TGraph *g_chi2_ams = new TGraph(); 
	TGraph *g_chi2_ace2 = new TGraph(); 
	TGraph *g_chi2_ams2 = new TGraph(); 
	TGraphErrors *g_pars[2];

	for (int ipar=0; ipar<2; ++ipar){
		g_pars[ipar] = new TGraphErrors(); 

	}

	TGraphErrors *g_ratio_ace[7]; 
	TGraphErrors *g_ratio_ace_fit[7];
	TGraphErrors *g_ratio_ams[h_ams_new->GetNbinsX()]; 
	TGraphErrors *g_ratio_ams_fit[h_ams_new->GetNbinsX()]; 

	TGraphErrors *g_ratio_ave = new TGraphErrors();
	TGraphErrors *g_ratio_ave_norm; // normalized averaged ratio of data/template-1 vs. time  

	TGraphErrors *g_residual_ace[7];
	TGraphErrors *g_resierr_ace[7]; 
	TGraphErrors *g_residual_norm_ace[7];
	TGraphErrors *g_residual_ams[h_ams_new->GetNbinsX()]; 
	TGraphErrors *g_residual_norm_ams[h_ams_new->GetNbinsX()];

	for (int bin=0; bin <= h_ace_rig_ave->GetNbinsX(); ++bin){ 
		if (bin%2==0){  
			g_ratio_ace[bin/2] = new TGraphErrors(nBRs); 
			g_ratio_ace_fit[bin/2] = new TGraphErrors(nBRs); 
			g_residual_ace[bin/2] = new TGraphErrors(nBRs);  
			g_resierr_ace[bin/2] = new TGraphErrors(nBRs); 
			HistTools::SetStyle(g_residual_ace[bin/2], kPink, kFullCircle, 1.4, 1, 1);
			HistTools::SetStyle(g_ratio_ave, kBlue+1, kFullCircle, 1.4, 1, 1);
		} 
	}  
 
	for (int bin=0; bin <= h_ams_new->GetNbinsX(); ++bin){
		g_ratio_ams[bin] = new TGraphErrors(nBRs);  
		g_ratio_ams_fit[bin] = new TGraphErrors(nBRs); 
		g_residual_ams[bin] = new TGraphErrors(nBRs);  
		HistTools::SetStyle(g_residual_ams[bin], kPink, kFullCircle, 1.4, 1, 1); 
	}   

	TFile *file_f4 = new TFile(Form("data/ACE/fill/F4_%s.root", element), "recreate"); 
 
	int iBR_true = 0;
	for (int iBR=0; iBR<nBRs; ++iBR){

		// AMS_He Monthly Spline 
		Spline *sp_he = new Spline("sp_he", nnodes_ams, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_he = sp_he->GetTF1Pointer();  
		TF1 *fit_he = (TF1*) file1->Get(Form("fsp_BR_he_%02d", iBR)); 

		HistTools::CopyParameters(fit_he, fsp_he); // error 
		fit_he->GetRange(x1,x2);
		fsp_he->SetRange(x1,x2);

		// AMS_He Integrated Spline  
		Spline *sp_he_ave = new Spline("sp_he_ave", nnodes_ams, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_he_ave = sp_he_ave->GetTF1Pointer();  
		TF1 *fit_he_ave = (TF1*) file1->Get("fsp_he"); 

		HistTools::CopyParameters(fit_he_ave, fsp_he_ave); // error 
		x1=0, x2=0;
		fit_he_ave->GetRange(x1,x2);
		fsp_he_ave->SetRange(x1,x2); 

		// AMS_He Model by Rescaled C ACE+AMS Combined Template  
		Spline *sp_he_temp = new Spline("sp_he_temp", nnodes, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_he_temp = sp_he_temp->GetTF1Pointer();  
		TF1 *fit_he_temp = (TF1*) file3->Get("fsp_he"); 

		HistTools::CopyParameters(fit_he_temp, fsp_he_temp); // error 
		x1=0, x2=0; 
		fit_he_temp->GetRange(x1,x2);
		fsp_he_temp->SetRange(x1,x2);  

		TF1 *fsp_he2 = HistTools::CombineTF1(fsp_he, fsp_he_ave, HistTools::Divide, "fsp_he2"); // He(R,t) Fit vs. <He(R)> Fit
		TF1 *fspind_he = HistTools::GetSpectralIndex(fsp_he2, "fspind_he", R1, R2); 
		fspind_he->SetRange(h_he_int->GetBinCenter(2), x2);

		// AMS_He(R,t)/<He(R)>
		TH1 *h_ratio0 = (TH1D*) h_BR_he[iBR]->Clone("h_ratio0"); 
		h_ratio0->Divide(h_he_int); 

		// refill AMS BR He to match the bins  
		HistTools::SetStyle(h_BR_he[iBR], kBlue, kFullCircle, 1.1, 1, 1); 

		TH1 *h_ace_BR = (TH1*) file2->Get(Form("h_rig_%s_BR%d", element, 2426+iBR_true)); 

		// rescaled <AMS_C> by AMS He BR / <AMS_He>   
		TH1 *h_ams_BR_fake = (TH1D*) h_ams_new->Clone("h_ams_BR_fake");

		h_ams_BR_fake->Multiply(h_BR_he[iBR]); 
		h_ams_BR_fake->Divide(h_he_int); 
 
		// h_ams_BR_fake->Print("range");  

		// rescaled combined fit 
		TH1 *h_ratio;
		if (get_index(element) <4) h_ratio = (TH1 *) h_ace_BR->Clone("h_ratio"); // ACE BR/Temp 	
		if (get_index(element)>=4) h_ratio = (TH1 *) h_ace_rig_ave->Clone("h_ratio"); // ACE Ave/Temp 	
		h_ratio->Divide(fsp_comb); 

		HistTools::SetStyle(h_ratio, kPink, kFullCircle, 1.4, 1, 1);

		double ratio_sum=0; // compute average of h_ratio manually  
		for(int nbin=0;nbin<14;++nbin){
			ratio_sum += h_ratio->GetBinContent(nbin);
			//printf("ratio_sum = %0.6f \n", ratio_sum);
		}
		double ratio_ave = ratio_sum/7;

		// printf("%0.6f \n", ratio_ave);
		//HistTools::PrintFunction(fsp_comb);
			
		double scale = 1./ratio_ave;
		if (get_index(element)<4) h_ratio->Scale(scale);	

		// rescaled combined fit 
		TF1 *rescaled_fit = HistTools::CombineTF1Const(fsp_comb, ratio_ave, HistTools::MultiplyConst, "rescaled_fit", R1, R2); 
		//if (get_index(element)<4) rescaled_fit = HistTools::CombineTF1(rescaled_fit, fsp_comb, HistTools::Divide, "fit_ratio", R1, R2);  

		TH1 *h_ratio1 = (TH1D*) h_ace_BR->Clone("h_ratio1"); // spind model  
		if (get_index(element) <4) h_ratio1->Divide(fsp_comb); 
		if (get_index(element)>=4) h_ratio1->Divide(rescaled_fit); 

		TH1 *h_ratio2 = (TH1D*) h_ams_BR_fake->Clone("h_ratio2"); // spind model 
		h_ratio2->Divide(h_ams_new);
   
		// h_ratio2->Print("range");

		// estimated AMS / resclaed combined fit
		TH1 *h_ratio3 = (TH1D *) h_ams_BR_fake->Clone("h_ratio3");		

		h_ratio3->Divide(rescaled_fit); 
		HistTools::SetStyle(h_ratio3, kPink, kFullCircle, 1.4, 1, 1);

		c0->cd(1);
		gPad->SetGrid(); 
		gPad->SetLogy();
		gPad->SetLogx(); 	

		TLegend *legend = new TLegend(0.62,0.75,0.9,0.9);
		legend->AddEntry(h_ams, Form("Integrated AMS %s Flux", element), "p"); 
		legend->AddEntry(h_ace_BR, Form("ACE %s BR Flux", element), "p"); 
		legend->AddEntry(h_ams_BR_fake, Form("Estimated AMS %s BR Flux", element), "p"); 
		legend->AddEntry(rescaled_fit, "rescaled combined template", "l"); 

		h_a1->SetTitle(Form("%s BR-%d Model Rigidity Spectrum; ; ", element, 2426+iBR_true));
		h_a1->SetXTitle(Unit::GetEnergyLabel("GV"));
  		h_a1->SetYTitle(Unit::GetDifferentialFluxLabel("GV m"));

		//if (get_temp_int(element)<4) h_a1->GetYaxis()->SetRangeUser(1e-3, 1e1); 
		// h_a1->GetXaxis()->SetRangeUser(0.01, 3000.); 

		TAxis *axis1 = h_a1->GetXaxis(); 
		axis1->SetLimits(0.7, 60.); 
	
		h_a1->Draw("E1X0"); 
	
		HistTools::SetStyle(h_ams, kPink-3, kFullCircle, 1.4, 1, 1); 
		HistTools::SetStyle(h_ace_BR, kRed, kFullCircle, 1.4, 1, 1);
		// HistTools::SetStyle(h_ams_BR_he, kBlue, kFullCircle, 1.4, 1, 1);
		HistTools::SetStyle(h_ams_BR_fake, kBlue, kFullCircle, 1.4, 1, 1);

		// fsp_he->Draw("SAME"); 
		rescaled_fit->Draw("SAME"); 
		h_ams->Draw("E1X0 SAME"); 
		h_ace_BR->Draw("E1X0 SAME"); 
		// h_ams_BR_he->Draw("E1X0 SAME"); 
		h_ams_BR_fake->Draw("E1X0 SAME"); 

		//fit->SetLineColor(kBlue); 
		//fit->Draw("SAME"); 
		legend->Draw("SAME");   

		c0->cd(2); 
		gPad->SetGrid(); 
		gPad->SetLogx(); 

		TAxis *axis2 = h_a2->GetXaxis(); 
		axis2->SetLimits(0.7, 60.); 

		//if (get_temp_int(element)<4) h_a2->GetYaxis()->SetRangeUser(0.5, 1.5); 
		h_a2->GetXaxis()->SetRangeUser(0.7, 60.); 

		h_a2->SetTitle("");
		h_a2->SetXTitle(Unit::GetEnergyLabel("GV"));
		h_a2->SetYTitle(Form("Estimated BR Fluxes / Rescaled %s Model", element)); 
		h_a2->SetTitleSize(0.3,"t"); 

		h_a2->Draw("E1X0"); 
		h_ratio->Draw("E1X0 SAME"); 
		h_ratio3->Draw("E1X0 SAME"); 

		if (iBR==0) { c0->Print(Form("./data/ACE/fill/fake_td_ams_v3/%s_flux_model.pdf(", element), "pdf"); } 
		if (iBR>0 && iBR<nBRs-1) { c0->Print(Form("./data/ACE/fill/fake_td_ams_v3/%s_flux_model.pdf", element), "pdf"); } 
		if (iBR==nBRs-1){
			c0->Print(Form("./data/ACE/fill/fake_td_ams_v3/%s_flux_model.pdf)", element), "pdf"); 
		} 

		// plot the ratio of AMS_He_BR / <AMS_He>
		c1->cd(1); 
		gPad->SetGrid(); 
		gPad->SetLogx(); 

		//if (get_temp_int(element)<4) h_a3->GetYaxis()->SetRangeUser(0.5, 1.5); 
		h_a3->GetXaxis()->SetLimits(0.7, 60);		

		TLegend *l1 = new TLegend(0.62,0.8,0.9,0.9); 
		l1->AddEntry(h_ratio0, "AMS He(R,t)/<He(R)>", "PL"); 	

		h_a3->Draw("E1X0"); 
		h_a3->SetXTitle(Unit::GetEnergyLabel("GV"));
		h_a3->SetTitle(Form("; ; BR-%d He(R,t)/<He(R)>", iBR_true+2426)); 
		HistTools::SetStyle(h_ratio0, kBlue, kFullCircle, 1.4, 1, 1);
		h_ratio0->Draw("E1X0 SAME"); 

		// plot spectral indices
		TGraphAsymmErrors *gspind_ace = HistTools::GetSpectralIndex(h_ace_BR, 5, 2); // get spectral indices 
		TGraphAsymmErrors *gspind_ams = HistTools::GetSpectralIndex(h_ams_BR_fake, 4, 1); 

		c1->cd(2); 
		gPad->SetGrid(); 
		gPad->SetLogx();  

		TLegend *l_both = new TLegend(0.62,0.8,0.9,1.); 
		l_both->AddEntry(gspind_ace, "ACE Spectral Indices", "PL");
		l_both->AddEntry(gspind_ams, "AMS Spectral Indices", "PL"); 		
	
		//if (get_temp_int(element)<4) gspind_ace->GetYaxis()->SetRangeUser(-3, 3); 
		gspind_ace->GetXaxis()->SetLimits(0.7, 60);

		HistTools::SetStyle(gspind_ace, kPink, kFullCircle, 1.4, 1, 1); 
		gspind_ace->SetTitle(Form("; ; ACE %s Spectral Indices BR-%d", element, 2426+iBR_true));
		gspind_ace->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV"));
		gspind_ace->Draw("APL"); 

		//h_a4->Draw("E1X0"); 
		HistTools::SetStyle(gspind_ams, kBlue, kFullCircle, 1.4, 1, 1); 
		gspind_ams->SetTitle(Form("; ; AMS %s Spectral Indices", element)); 
		gspind_ams->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV"));
		gspind_ams->GetXaxis()->SetTitleSize(1.5); 
		gspind_ams->Draw("PL SAME"); 
		// fspind->Draw("SAME"); 
		l_both->Draw("SAME"); 

		if (iBR==0) { c1->Print(Form("./data/ACE/fill/fake_td_ams_v3/fake_ams_spind_%s.pdf(", element), "pdf"); }
		if (iBR>0 && iBR<nBRs-1) { c1->Print(Form("./data/ACE/fill/fake_td_ams_v3/fake_ams_spind_%s.pdf", element), "pdf"); }
		if (iBR==nBRs-1){
			// gspind_ace->Print("range");
			// cout << " " << endl; 
			// gspind_ams->Print("range"); 
			c1->Print(Form("./data/ACE/fill/fake_td_ams_v3/fake_ams_spind_%s.pdf)", element), "pdf"); 
		} 

		// find spectral indices of h_ratio + h_ratio0 
		TGraphAsymmErrors *gspind_ace2 = HistTools::GetSpectralIndex(h_ratio1, 5, 2); // get spectral indices 
		TGraphAsymmErrors *gspind_ams2 = HistTools::GetSpectralIndex(h_ratio2, 4, 1); 

		HistTools::SetStyle(gspind_ace2, kPink, kFullCircle, 1.4, 1, 1); 
		HistTools::SetStyle(gspind_ams2, kBlue, kFullCircle, 1.4, 1, 1); 
		HistTools::SetStyle(fspind_he, kBlack, kFullCircle, 0.8, 1, 1);

		// fit to F(R,t)/<F(R)>
		TF1 *fit_ratio = new TF1("f_ratio", "1+[0]*exp(-[1]*x)", 0.7, 60);
		TF1 *fit_ratio2 = new TF1("f_ratio2", "1+[0]*exp(-[1]*log(x))", 0.7, 60);

		double F1 = h_ratio1->GetBinContent(h_ratio1->FindBin(R1)), F2 = h_ratio0->GetBinContent(h_ratio0->GetNbinsX()); 
		double A1 = log(abs(F1-1)), A2 = log(abs(F2-1)); 
		TF1 *f1 = new TF1("f1","abs(sin(x)/x)*sqrt(x)",0,1);
		double rho = f1->GetRandom();
		double N = f1->GetRandom(); 
		
		// double rho = (R2-R1)/(A1-A2)
		// double N = 0.5*(abs(F1-1)/exp(-R1/rho) + abs(F2-1)/exp(-R2/rho)); 
		
		last_pars = {N, rho}; 
		last_pars2 = {N, rho};

		// printf(" rho = %0.4f, N = %0.20f \n ", rho, N); 
		// printf(" rho = %0.4f, N = %0.20f, F1 = %0.4f, F2 = %0.4f, (F1-1)/exp(-R1/rho) = %0.6f, (F2-1)/exp(-R2/rho) = %0.6f \n", rho, N, F1, F2, (F1-1)/exp(-R1/rho), (F2-1)/exp(-R2/rho));

		for (int ipar=0; ipar<fit_ratio->GetNpar(); ++ipar){ 
			fit_ratio->SetParameter(ipar, last_pars[ipar]); 
			fit_ratio2->SetParameter(ipar, last_pars2[ipar]); 
		}   

		TObjArray data; 
		data.Add(h_ratio1);
		data.Add(h_ratio0);

		// fit AMS and ACE data at the same time
		vector<double> rigmin, rigmax, chi2norm(2);
		ROOT::Fit::Fitter fitter;
		//HistTools::PrintFunction(fit);
		FitTools::SetCommonFitterOptions(fitter);
		FitTools::FitCombinedData(data, fit_ratio, "I", rigmin, rigmax, chi2norm, fitter, 3); 

		vector<double> rigmin2, rigmax2, chi2norm2(2);
		ROOT::Fit::Fitter fitter2;
		FitTools::SetCommonFitterOptions(fitter2);
		FitTools::FitCombinedData(data, fit_ratio2, "I", rigmin2, rigmax2, chi2norm2, fitter2, 3); 
 
		for (int ipar=0; ipar<fit_ratio->GetNpar(); ++ipar){
			g_pars[ipar]->SetPoint(iBR, UBRToTime(iBR_true+2426), fit_ratio->GetParameter(ipar)); 
			g_pars[ipar]->SetPointError(iBR, 0, fit_ratio->GetParError(ipar));   
		}  

		TF1 *fspind_ratio = HistTools::GetSpectralIndex(fit_ratio, "fspind_ratio", R1, R2); 
		TF1 *fspind_ratio2 = HistTools::GetSpectralIndex(fit_ratio2, "fspind_ratio2", R1, R2); 

		//h_ams_BR_fake->Print("range"); 
		c2->cd(1); 
		gPad->SetGrid(); 
		gPad->SetLogx();
		gPad->SetBottomMargin(0.01); 

		//if (get_temp_int(element)<4) h_a4->GetYaxis()->SetRangeUser(0.5, 2.2); 
		h_a4->GetXaxis()->SetLimits(0.7, 60); 
	
		TLegend *l2 = new TLegend(0.62,0.8,0.9,1.); 
		l2->AddEntry(h_ratio1, Form("ACE %s(R,t)/<%s(R)>", element, element), "PL");
		l2->AddEntry(h_ratio0, "AMS He(R,t)/<He(R)>", "PL");  
		l2->AddEntry(fit_ratio, Form("1+%6.3fexp(-%6.3fR)", fit_ratio->GetParameter(0), fit_ratio->GetParameter(1)), "L"); 
		l2->AddEntry(fit_ratio2, Form("1+%6.3fexp(-%6.3fln(R))", fit_ratio2->GetParameter(0), fit_ratio2->GetParameter(1)), "L"); 

		h_a4->Draw("E1X0"); 
		h_a4->SetXTitle(Unit::GetEnergyLabel("GV"));
		h_a4->SetTitle(Form("; ; BR-%d %s(R,t)/<%s(R)>", iBR_true+2426, element, element));
		// h_a4->SetTitleSize(0.5,"y"); 
		HistTools::SetStyle(h_ratio1, kPink, kFullCircle, 1.4, 1, 1);
		HistTools::SetStyle(h_ratio0, kBlue, kFullCircle, 1.4, 1, 1);
		fit_ratio2->SetLineColor(kGreen); 
		fit_ratio->Draw("SAME");
		fit_ratio2->Draw("SAME"); 
		h_ratio1->Draw("E1X0 SAME");
		h_ratio0->Draw("E1X0 SAME");
		l2->Draw("SAME"); 

		TH1D *h_fitres[2];
		TH1D *h_fitres2[2]; 
		TLegend *l_chi2 = new TLegend(0.62,0.8,0.9,1.0); 

		if (!strcmp(element, "N")){
			last_pars = {par0_X[0]->GetY()[iBR], par1_X[0]->GetY()[iBR]};
			for (int ipar=0; ipar<fit_ratio->GetNpar(); ++ipar){ 
				fit_ratio->SetParameter(ipar, last_pars[ipar]); 
			}
		} 

		if (!strcmp(element, "C")){
			last_pars = {par0_X[3]->GetY()[iBR], par1_X[3]->GetY()[iBR]};   
			for (int ipar=0; ipar<fit_ratio->GetNpar(); ++ipar){ 
				fit_ratio->SetParameter(ipar, last_pars[ipar]); 
			}			
		} 

		else if (get_index(element)>3){
			last_pars = {par0_X[choose_BO(element)]->GetY()[iBR], par1_X[choose_BO(element)]->GetY()[iBR]};
			for (int ipar=0; ipar<fit_ratio->GetNpar(); ++ipar){ 
				fit_ratio->SetParameter(ipar, last_pars[ipar]); 
			}
		}

		for (UShort_t i = 0; i < 2; ++i)
		{
  			TH1 *hist = HistTools::ToHist(data[i]);
			TH1 *h_fiterr = (TH1*) HistTools::GetFitError( (TH1*) data[i], fit_ratio, "_fiterr", false, false, true, 10, DBL_MAX, 0.68, &fitter);
   			h_fitres[i] = (TH1D *)HistTools::GetResiduals(hist, fit_ratio, "_fitres", false, true, true, 5, 1, 0.68, &fitter);
			HistTools::CopyStyle(hist, h_fitres[i]);

   			UShort_t ndf  = hist->GetNbinsX();
   			Double_t chi2 = hist->Chisquare(fit_ratio, "R"); 
   			// printf(" %d nodes ### %-34s   chi2/ndf=%6.2f/%-2u   chi2norm=%5.2f   prob.=%5.2f%%\n", nnodes, hist->GetTitle(), chi2, ndf, chi2norm[i], TMath::Prob(chi2, ndf)*1e2);

			l_chi2->AddEntry(h_fitres[i], Form("1+[0]*exp(-[1]*R) chi2/ndf=%6.2f/%-2u", chi2, ndf), "P");  

			if (i==0){ 
				g_chi2_ace->SetPoint(iBR, UBRToTime(iBR_true+2426), chi2/ndf);
				for (int bin=0; bin <= h_fitres[i]->GetNbinsX(); ++bin){ 
					if (bin%2==0){ 
						g_ratio_ace[bin/2]->SetPoint(iBR, UBRToTime(iBR_true+2426), h_ratio1->GetBinContent(bin+1)); 
						g_ratio_ace[bin/2]->SetPointError(iBR, 0, h_ratio1->GetBinError(bin+1)); 
						g_ratio_ace_fit[bin/2]->SetPoint(iBR, UBRToTime(iBR_true+2426), fit_ratio->Eval(h_ratio1->GetBinCenter(bin+1))); 
						g_ratio_ace_fit[bin/2]->SetPointError(iBR, 0, h_fiterr->GetBinError(bin+1)); 
						g_residual_ace[bin/2]->SetPoint(iBR, UBRToTime(iBR_true+2426), h_fitres[i]->GetBinContent(bin+1)); 
						g_residual_ace[bin/2]->SetPointError(iBR, 0, h_fitres[i]->GetBinError(bin+1)); 

						double df = h_fitres[i]->GetBinContent(bin+1)*hist->GetBinContent(bin+1)/hist->GetBinError(bin+1); 
						double ddf = df*sqrt(pow(h_fitres[i]->GetBinError(bin+1)/h_fitres[i]->GetBinContent(bin+1))+pow(h_fiterr->GetBinError(bin+1)/h_fiterr->GetBinContent(bin+1))); // d(data-fit) 
						g_resierr_ace[bin/2]->SetPoint(iBR, UBRToTime(iBR_true+2426), df); 
						g_resierr_ace[bin/2]->SetPointError(iBR, 0, ddf); 

						// printf("bin/2 = %d, iBR = %d, time = %10.4u, fitres = %10.4f \n", bin/2, iBR, UBRToTime(iBR_true+2426), h_fitres[0]->GetBinContent(bin+1)); 
					} 
				}   
			} 

			if (i==1){  
				g_chi2_ams->SetPoint(iBR, UBRToTime(iBR_true+2426), chi2/ndf);
				for (int bin=0; bin <= h_ams_new->GetNbinsX()-1; ++bin){ 
					g_ratio_ams[bin]->SetPoint(iBR, UBRToTime(iBR_true+2426), h_ratio0->GetBinContent(bin+1)); 
					g_ratio_ams[bin]->SetPointError(iBR, 0, h_ratio0->GetBinError(bin+1));
					g_ratio_ams_fit[bin]->SetPoint(iBR, UBRToTime(iBR_true+2426), fit_ratio->Eval(h_ratio0->GetBinCenter(bin+1))); 
					g_ratio_ams_fit[bin]->SetPointError(iBR, 0, h_fiterr->GetBinError(bin+1)); 
					g_residual_ams[bin]->SetPoint(iBR, UBRToTime(iBR_true+2426), h_fitres[i]->GetBinContent(bin+1));
					g_residual_ams[bin]->SetPointError(iBR, 0, h_fitres[i]->GetBinError(bin+1)); 
					// printf("bin = %d, iBR = %d, time = %10.4u, fitres = %10.4f \n", bin, iBR, UBRToTime(iBR_true+2426), h_fitres[1]->GetBinContent(bin+1)); 
				}  
			} 	
		}	

		for (UShort_t i = 0; i < 2; ++i)
		{
  			TH1 *hist = HistTools::ToHist(data[i]);
   			h_fitres2[i] = (TH1D *)HistTools::GetResiduals(hist, fit_ratio2, "_fitres2", false, true, true, 5, 1, 0.68, &fitter2);
			HistTools::CopyStyle(hist, h_fitres2[i]);

   			UShort_t ndf  = hist->GetNbinsX();
   			Double_t chi2 = chi2norm2[i]*ndf; 
   			// printf(" %d nodes ### %-34s   chi2/ndf=%6.2f/%-2u   chi2norm=%5.2f   prob.=%5.2f%%\n", nnodes, hist->GetTitle(), chi2, ndf, chi2norm[i], TMath::Prob(chi2, ndf)*1e2); 

			l_chi2->AddEntry(h_fitres2[i], Form("1+[0]*exp(-[1]*ln(R)) chi2/ndf=%6.2f/%-2u", chi2, ndf), "P"); 

			if (i==0){ 
				g_chi2_ace2->SetPoint(iBR, UBRToTime(iBR_true+2426), chi2);
			} 

			if (i==1){ 
				g_chi2_ams2->SetPoint(iBR, UBRToTime(iBR_true+2426), chi2);
			} 	
		}

		TH1 *h_comp1 = (TH1*) h_ace_BR->Clone("h_comp1");
		h_comp1->Divide(rescaled_fit); 

		double comp1_sum=0, dcomp1_sum=0; // compute average ratio of data/template-1 of ACE F(R,t)/<F(R)> manually 		

		for(int nbin=0;nbin<h_ace_rig_ave->GetNbinsX();++nbin){
			comp1_sum += h_comp1->GetBinContent(nbin); 
			dcomp1_sum += h_comp1->GetBinError(nbin)*h_comp1->GetBinError(nbin); 
			//printf("comp1_sum = %0.6f \n", comp1_sum);
		}
		double comp1_ave = comp1_sum/7; 
		double dcomp1_ave = sqrt(dcomp1_sum)/7; 

		g_ratio_ave->SetPoint(iBR, UBRToTime(iBR_true+2426), comp1_ave-1); 
		g_ratio_ave->SetPointError(iBR, 0, dcomp1_ave); 

		c2->cd(2);
		gPad->SetGrid();
		gPad->SetLogx(); 
		gPad->SetTopMargin(0.01); 

		//if (get_temp_int(element)<4) h_a5->GetYaxis()->SetRangeUser(-0.15, 0.15)  ; 
		h_a5->GetXaxis()->SetLimits(0.7, 60); 

		h_a5->Draw("E1X0"); 
		h_a5->SetXTitle(Unit::GetEnergyLabel("GV"));
		h_a5->SetTitle("; ; Data/Fit-1"); 

		h_fitres2[0]->Draw("E1X0 SAME"); 
		h_fitres2[1]->Draw("E1X0 SAME");

		h_fitres[0]->Draw("E1X0 SAME");
		h_fitres[1]->Draw("E1X0 SAME");

		HistTools::SetStyle(h_fitres2[0], kGreen-1, kFullCircle, 1.4, 1, 1);
		HistTools::SetStyle(h_fitres2[1], kGreen-1, kFullCircle, 1.4, 1, 1); 

		l_chi2->Draw("SAME"); 		

		c2->cd(3);
		gPad->SetGrid();
		gPad->SetLogx(); 

		TLegend *l_both3 = new TLegend(0.62,0.7,0.8,1.0); 
		l_both3->AddEntry(gspind_ace2, Form("ACE %s(R,t)/<%s(R)> Spectral Indices", element, element), "PL");
		l_both3->AddEntry(gspind_ams2, Form("AMS %s(R,t)/<%s(R)> Spectral Indices", element, element), "PL"); 
		// l_both2->AddEntry(fspind_he, "AMS He(R,t)_Fit/temp<He(R)>_Fit Spectral Indices", "L");  	
		l_both3->AddEntry(fspind_ratio, "Spectral Index of 1+[0]*exp(-[1]*R) Fit to He(R,t)/<He(R)>", "L");	
		l_both3->AddEntry(fspind_he, "AMS He(R,t)_Fit/<He(R)>_Fit Spectral Indices", "L");  
		l_both3->AddEntry(fspind_ratio2, "Spectral Index of 1+[0]*exp(-[1]*ln(R)) Fit to He(R,t)/<He(R)>", "L"); 

		//if (get_temp_int(element)<4) gspind_ace2->GetYaxis()->SetRangeUser(-2, 1.3); 
		TAxis *axis5 = gspind_ace2->GetXaxis(); 
		axis5->SetLimits(0.7, 60.); 

		gspind_ace2->SetTitle(Form("; ; Spectral Indices Model of %s F(R,t)/<F(R)> at BR-%d", element, 2426+iBR_true));
		gspind_ace2->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV"));
		gspind_ace2->Draw("APL");  
		fspind_ratio->Draw("SAME"); 
		fspind_ratio2->Draw("SAME"); 
		fspind_ratio2->SetLineColor(kGreen); 
		fspind_he->Draw("SAME");
		gspind_ams2->Draw("PL SAME"); 
		gspind_ace2->Draw("PL SAME");
		l_both3->Draw("SAME"); 

		// break; 

		rescaled_fit->Write(Form("rescaled_fit_BR%d", 2426+iBR_true)); 
		fit_he->Write(Form("fit_he_BR%d", 2426+iBR_true)); 
		//fsp_he2->Write(Form("fsp_he2_BR%d", 2426+iBR_true)); 
		fit_ratio->Write(Form("fit_ratio_BR%d", 2426+iBR_true)); 
		fit_he_ave->Write(Form("fit_he_ave_BR%d", 2426+iBR_true)); 
		h_ratio1->Write(Form("h_ratio1_BR%d", 2426+iBR_true));
		h_ratio0->Write(Form("h_ratio0_BR%d", 2426+iBR_true)); 
		h_fitres[0]->Write(Form("h_fitres_ace_BR%d", 2426+iBR_true));
		h_fitres[1]->Write(Form("h_fitres_ams_BR%d", 2426+iBR_true));

		if (iBR==0) c2->Print(Form("./data/ACE/fill/fake_td_ams_v3/spind_model_%s_temp.pdf(", element), "pdf"); 
		if (iBR>0 && iBR<nBRs-1) c2->Print(Form("./data/ACE/fill/fake_td_ams_v3/spind_model_%s_temp.pdf", element), "pdf"); 
		if (iBR==nBRs-1){ 
			//gspind_ace2->Print("range");
			cout << " " << endl; 
			//gspind_ams2->Print("range");
			cout << " " << endl; 
			// h_ratio2->Print("range");
			c2->Print(Form("./data/ACE/fill/fake_td_ams_v3/spind_model_%s_temp.pdf)", element), "pdf"); 
		}

	   	if (iBR+2426==2472-1) iBR_true += 3; 
		else iBR_true ++; 
	}
			
	c3->cd(1);
	gPad->SetGrid(); 
	HistTools::SetStyle(g_chi2_ace, kPink, kFullCircle, 1.4, 1, 1); 
	HistTools::SetStyle(g_chi2_ace2, kBlue, kFullCircle, 1.4, 1, 1); 
	TLegend *l_c3_1 = new TLegend(0.62,0.8,0.9,0.9); 
	l_c3_1->AddEntry(g_chi2_ace, "1+[0]*exp(-[1]*R)", "PL"); 
	l_c3_1->AddEntry(g_chi2_ace2, "1+[0]*exp(-[1]*ln(R))", "PL"); 
	g_chi2_ace->SetTitle(Form("; ; ACE %s(R,t)/<%s(R)> Fitting Chi2/NDF", element, element)); 
	g_chi2_ace->GetXaxis()->SetTimeDisplay(1);
	g_chi2_ace->GetXaxis()->SetTimeFormat("%m-%y"); 
	g_chi2_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	g_chi2_ace->GetXaxis()->SetTitleSize(0.7); 
	////if (get_temp_int(element)<4) g_chi2_ace->GetYaxis()->SetRangeUser(0., 4.5); 
	g_chi2_ace->Draw("APL"); 
	g_chi2_ace2->Draw("PL SAME"); 
	l_c3_1->Draw("SAME");

	c3->cd(2); 
	gPad->SetGrid(); 
	HistTools::SetStyle(g_chi2_ams, kPink, kFullCircle, 1.4, 1, 1); 
	HistTools::SetStyle(g_chi2_ams2, kBlue, kFullCircle, 1.4, 1, 1); 
	TLegend *l_c3_2 = new TLegend(0.62,0.8,0.9,0.9); 
	l_c3_2->AddEntry(g_chi2_ams, "1+[0]*exp(-[1]*R)", "PL");
	l_c3_2->AddEntry(g_chi2_ams2, "1+[0]*exp(-[1]*ln(R))", "PL");
	g_chi2_ams->SetTitle(Form("; ; AMS %s(R,t)/<%s(R)> Fitting Chi2/NDF", element, element)); 
	g_chi2_ams->GetXaxis()->SetTimeDisplay(1);
	g_chi2_ams->GetXaxis()->SetTimeFormat("%m-%y");
	g_chi2_ams->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	g_chi2_ams->GetXaxis()->SetTitleSize(0.7);
	// //if (get_temp_int(element)<4) g_chi2_ams->GetYaxis()->SetRangeUser(0., 2.2); 
	g_chi2_ams->Draw("APL");
	g_chi2_ams2->Draw("PL SAME");  
	l_c3_2->Draw("SAME"); 

	c3->Print(Form("./data/ACE/fill/fake_td_ams_v3/chi2_vs_BR_%s.png", element)); 

	g_ratio_ave_norm = get_norm_graph( g_ratio_ave ); 

	for (int iBR=0; iBR<nBRs; ++iBR){
				Double_t x=0, y=0; 
				g_ratio_ave_norm->GetPoint(iBR, x, y); 
				g_ratio_ave_norm->SetPoint(iBR, x, y-1); 			 
	}	

	TGraphErrors *g_ave_residual_ace_rig = new TGraphErrors(); // time-averaged residual vs. rigidity
	TGraphErrors *g_ave_residual_ams_rig = new TGraphErrors(); // time-averaged residual vs. rigidity 

	for (int bin=0; bin <= h_ace_rig_ave->GetNbinsX(); ++bin){
 
		if (bin%2==0) {  

			g_ratio_time_norm[bin/2] = get_norm_graph( g_ratio_time[bin/2] ); 
		
			for (int iBR=0; iBR<g_ratio_time_norm[bin/2]->GetN(); ++iBR){
				Double_t x=0, y=0; 
				g_ratio_time_norm[bin/2]->GetPoint(iBR, x, y); 
				g_ratio_time_norm[bin/2]->SetPoint(iBR, x, y-1); 			 
			}

			g_ave_residual_ace_rig->SetPoint(bin/2, h_ace_rig_ave->GetBinCenter(bin+1), g_residual_ace[bin/2]->GetMean(2) ); 
			g_ave_residual_ace_rig->SetPointError(bin/2, 0, g_residual_ace[bin/2]->GetRMS(2) ); 

			c4->cd(1);
			gPad->SetGrid(); 

			TLegend *l_comp1 = new TLegend(0.62, 0.8, 0.9, 1.0); 
			l_comp1->AddEntry(g_ratio_ave_norm, "Normalized Averaged ACE Data/Template-1", "PL"); 
			l_comp1->AddEntry(g_residual_ace[bin/2], Form("ACE Fitting Residuals (%0.4f GV)", h_ace_rig_ave->GetBinCenter(bin+1)), "PL"); 

			g_residual_ace[bin/2]->SetTitle(Form(" ; ; ACE %s Data/Template-1 vs. Fitting Residuals", element));  
			g_residual_ace[bin/2]->GetXaxis()->SetTimeDisplay(1);
			g_residual_ace[bin/2]->GetXaxis()->SetTimeFormat("%m-%y");
			g_residual_ace[bin/2]->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g_residual_ace[bin/2]->GetXaxis()->SetTitleSize(0.7);
			//if (get_temp_int(element)<4) g_residual_ace[bin/2]->GetYaxis()->SetRangeUser(-0.25, 0.3);
			g_residual_ace[bin/2]->Draw("APL");
			// g_ratio_ave_norm->Draw("PL SAME"); 
			g_ratio_time_norm[bin/2]->Draw("PL SAME"); 
			l_comp1->Draw("SAME"); 		

			if (bin==0) c4->Print(Form("./data/ACE/fill/fake_td_ams_v3/fitres_ace_vs_BR_%s.pdf(", element), "pdf"); 
			if (bin>0 && bin<h_ace_rig_ave->GetNbinsX()-1) c4->Print(Form("./data/ACE/fill/fake_td_ams_v3/fitres_ace_vs_BR_%s.pdf", element), "pdf");  
			if (bin==h_ace_rig_ave->GetNbinsX()-1){
				c4->Print(Form("./data/ACE/fill/fake_td_ams_v3/fitres_ace_vs_BR_%s.pdf)", element), "pdf"); 
			} 
		
			g_ratio_ace[bin/2]->Write(Form("g_ratio_ace_%d", bin/2)); 
			g_ratio_ace_fit[bin/2]->Write(Form("g_ratio_ace_fit_%d", bin/2)); 
			g_residual_ace[bin/2]->Write(Form("g_residual_ace_%d", bin/2)); 
			g_resierr_ace[bin/2]->Write(Form("g_resierr_ace_%d", bin/2)); 
		} 
	} 

	for (int bin=0; bin <= h_ams_new->GetNbinsX()-1; ++bin){ 

			g_ave_residual_ams_rig->SetPoint(bin, h_ams_new->GetBinCenter(bin+1), g_residual_ams[bin]->GetMean(2) ); 
			g_ave_residual_ams_rig->SetPointError(bin, 0, g_residual_ams[bin]->GetRMS(2) ); 

			c4_2->cd(1); 
			gPad->SetGrid(); 

			g_residual_ams[bin]->SetTitle(Form("AMS %s(R,t)/<%s(R)>; ; Fitting Residuals (%0.4f GV)", element, element, h_ams_new->GetBinCenter(bin+1)));  
			g_residual_ams[bin]->GetXaxis()->SetTimeDisplay(1);
			g_residual_ams[bin]->GetXaxis()->SetTimeFormat("%m-%y");
			g_residual_ams[bin]->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g_residual_ams[bin]->GetXaxis()->SetTitleSize(0.7);
			//if (get_temp_int(element)<4) g_residual_ams[bin]->GetYaxis()->SetRangeUser(-0.06, 0.1);  
			// g_residual_ams[bin]->Print("range"); 
			g_residual_ams[bin]->Draw("APL"); 

			if (bin==0) c4_2->Print(Form("./data/ACE/fill/fake_td_ams_v3/fitres_ams_vs_BR_%s.pdf(", element), "pdf"); 
			if (bin>0 && bin<h_ams_new->GetNbinsX()-1) c4_2->Print(Form("./data/ACE/fill/fake_td_ams_v3/fitres_ams_vs_BR_%s.pdf", element), "pdf");  
			if (bin==h_ams_new->GetNbinsX()-1){
				c4_2->Print(Form("./data/ACE/fill/fake_td_ams_v3/fitres_ams_vs_BR_%s.pdf)", element), "pdf"); 
			} 

			g_ratio_ams[bin]->Write(Form("g_ratio_ams_%d", bin)); 
			g_ratio_ams_fit[bin]->Write(Form("g_ratio_ams_fit_%d", bin)); 
			g_residual_ams[bin]->Write(Form("g_residual_ams_%d", bin)); 
	} 

	for (int ipar=0; ipar<2; ++ipar){
		c5->cd(ipar+1); 
		gPad->SetGrid(); 
		g_pars[ipar]->Draw("APL"); 
		HistTools::SetStyle(g_pars[ipar], kPink, kFullCircle, 1.4, 1, 1); 
		g_pars[ipar]->GetXaxis()->SetTimeDisplay(1);
		g_pars[ipar]->GetXaxis()->SetTimeFormat("%m-%y");
		g_pars[ipar]->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
		g_pars[ipar]->GetXaxis()->SetTitleSize(0.7);
		//if (get_temp_int(element)<4) g_pars[ipar]->GetYaxis()->SetRangeUser(-3., 3.); 
		if (ipar==0) g_pars[ipar]->SetTitle(Form("1+[0]*exp(-[1]*R); ; Parameter [%d]", ipar));
		else g_pars[ipar]->SetTitle(Form(" ; ; Parameter [%d]", ipar));
	}		

	TGraphErrors *g_pars_ratio = new TGraphErrors(); 
	for (int i=0; i<g_pars[0]->GetN(); ++i){
		g_pars_ratio->SetPoint(i, g_pars[0]->GetX()[i], g_pars[1]->GetY()[i]/g_pars[0]->GetY()[i]); 
		g_pars_ratio->SetPointError(i, 0, g_pars[1]->GetY()[i]/g_pars[0]->GetY()[i]*sqrt(pow(g_pars[1]->GetEY()[i]/g_pars[1]->GetY()[i], 2)+pow(g_pars[0]->GetEY()[i]/g_pars[0]->GetY()[i], 2))); 
	}

	c5->cd(3);
	gPad->SetGrid();
	HistTools::SetStyle(g_pars_ratio, kPink, kFullCircle, 1.4, 1, 1);
	g_pars_ratio->GetXaxis()->SetTimeDisplay(1);
	g_pars_ratio->GetXaxis()->SetTimeFormat("%m-%y");
	g_pars_ratio->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	g_pars_ratio->GetXaxis()->SetTitleSize(0.7);
	//if (get_temp_int(element)<4) g_pars_ratio->GetYaxis()->SetRangeUser(-15., 30.); 
	g_pars_ratio->SetTitle(Form(" ; ; [%d]/[%d]", 1, 0)); 
	g_pars_ratio->Draw("AP"); 

	// g_pars_ratio->Print("range"); 

	c5->Print(Form("./data/ACE/fill/fake_td_ams_v3/pars_vs_BR_%s.png", element)); 

	TCanvas *c6 = new TCanvas("c6", "", 1600, 900); // plot time-averaged fit residuals to F(R,t)/<F(R)> vs. rigidity 
	c6->Divide(1, 1);  

	c6->cd(1); 
	gPad->SetGrid(); 
	gPad->SetLogx(); 

	HistTools::CopyStyle( h_ace_rig_ave, g_ave_residual_ace_rig ); 
	HistTools::CopyStyle( h_ams_new, g_ave_residual_ams_rig );
	////if (get_temp_int(element)<4) g_ave_residual_ace_rig->GetYaxis()->SetRangeUser(-0.15, 0.15); 
	g_ave_residual_ace_rig->GetXaxis()->SetLimits(0.7, 60);
	g_ave_residual_ace_rig->SetTitle(Form(" ; ; Time-averaged Fit Residuals of %s(R,t)/<%s(R)>", element, element)); 
	g_ave_residual_ace_rig->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV")); 

	TLegend *l_both4 = new TLegend(0.62,0.8,0.9,0.9); 
	l_both4->AddEntry(g_ave_residual_ace_rig, "ACE", "PL");
	l_both4->AddEntry(g_ave_residual_ams_rig, "AMS", "PL"); 

	g_ave_residual_ace_rig->Draw("APL");
	g_ave_residual_ams_rig->Draw("PLSAME"); 
	l_both4->Draw("SAME"); 

	g_chi2_ace->Write("g_chi2_ace"); 
	g_chi2_ams->Write("g_chi2_ams"); 

	g_ave_residual_ace_rig->Write("g_ave_residual_ace_rig"); 
	g_ave_residual_ams_rig->Write("g_ave_residual_ams_rig"); 

	file_f4->Close(); 

	c6->Print(Form("./data/ACE/fill/fake_td_ams_v3/time_ave_residual_vs_rig_%s.png", element)); 

	return 0; 
} 

// F5  
void ace_fake_td_ams_v4(const char *element, Particle::Type isotope){

	gStyle->SetOptStat(0); 
	// gStyle->SetTitleSize(26,"t");

	gSystem->mkdir("data/ACE/fill/fake_td_ams_v4", true);	

 	Experiments::DataPath = "data";	
   	int data_value[n_ele] = {0, 18, 23, 27, 31, 20, 43, 22}; // p, He, Li, Be, B, C, N, O

	TF1 *f_corr[4]; 
	TH1 *h_corr[4]; 
	TH1 *h_correrr[4]; 
	TGraphErrors *g_correrr[4]; 
	int best_bin[4]; 

	TGraphErrors *par0_X[4];
	TGraphErrors *par1_X[4]; 

	TFile *file5 = new TFile("./data/ACE/fill/fake_td_ams_v2/par0_fit.root");

	for (int i=0; i<4; ++i){ 

		f_corr[i] = (TF1 *) file5->Get(Form("f_corr2_%s", ACE_Element[i]));  
		h_corr[i] = (TH1 *) file5->Get(Form("h_corr2_%s", ACE_Element[i])); 
		h_correrr[i] = (TH1 *) file5->Get(Form("h_correrr_%s", ACE_Element[i])); 
		g_correrr[i] = (TGraphErrors *) file5->Get(Form("h_correrr_%s", ACE_Element[i])); 
		best_bin[i] = h_corr[i]->GetMaximumBin()-1; 

		par0_X[i] = new TGraphErrors(Form("data/ACE/fill/fake_td_ams_v2/fit_pars_%s.dat", ACE_Element[i]), "%lg %lg %lg %*s %*s"); 
		par1_X[i] = new TGraphErrors(Form("data/ACE/fill/fake_td_ams_v2/fit_pars_%s.dat", ACE_Element[i]), "%lg %*s %*s %lg %lg"); 

	}

	double *kin_bins = get_kin_bins(element);
        double *SpallCorr = get_spall_corr(element);
	double *SpallCorrUnc = get_spall_corr_unc(element);
	double *EMed = get_EMed(element);

	const UInt_t FirstACEBR = 2240;
	vector<UInt_t> BRs;
	// we stop at BR 2493, which ends on 2016/05/23, just 3 days before the end of the data taking period for AMS nuclei
	for (UInt_t br=2426; br<=2493; ++br) { 
		if (br != 2472 && br != 2473) BRs.push_back(br-FirstACEBR); 
	}

	// load AMS He BR fluxes 
	const int nBRs = Experiments::Info[Experiments::AMS02].Dataset[1].nMeasurements; // read the number of BRs
	TH1D **h_BR_he = Experiments::GetDatasetHistograms(Experiments::AMS02, 4);
	// Create a new set of histograms which match the AMS monthly bins with the AMS integrated bins 

	TH1 *h_ams = Experiments::GetMeasurementHistogram(Experiments::AMS02, get_ams_data_value(element), 0); // load AMS data for a given template 
	TH1 *h_ams_new = (TH1D*) h_BR_he[0]->Clone("h_ams_new"); 
	// TH1 *h_he_int = (TH1*) ave_hist( h_BR_he, nBRs);  
	TH1 *h_he = Experiments::GetMeasurementHistogram(Experiments::AMS02, 18, 0);
	TH1 *h_he_int = (TH1D*) h_BR_he[0]->Clone("h_ams_new");  

	// create h_ams_new and h_he_int to match the bins  

	for (int bin=1; bin <= h_ams_new->GetNbinsX(); ++bin){ 
	   if (!strcmp(element, "B") || !strcmp(element, "C")){ 
		h_ams_new->SetBinContent(bin, h_ams->GetBinContent(bin)); 
		h_ams_new->SetBinError(bin, h_ams->GetBinError(bin));
	   }
	   if (!strcmp(element, "N") || !strcmp(element, "O")){ 
		h_ams_new->SetBinContent(bin, h_ams->GetBinContent(bin-1)); 
		h_ams_new->SetBinError(bin, h_ams->GetBinError(bin-1)); 
	   }
	}
	if (!strcmp(element, "N") || !strcmp(element, "O")){
		h_ams_new->SetBinContent(1, 0); 
		h_ams_new->SetBinError(1, 0);
	}
 
	for (int bin=1; bin<=h_BR_he[0]->GetNbinsX(); ++bin) { 
		h_he_int->SetBinContent(bin, h_he->GetBinContent(bin)); 
		h_he_int->SetBinError(bin, h_he->GetBinError(bin)); 
	}

	// h_ams->Print("range"); 
	// h_ams_new->Print("range"); 

	// load ACE average 
	TH1 *h_ace_ene_ave = HistTools::GraphToHist(get_ace_average_graph( element , &BRs[0], BRs.size() ), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
	TH1 *h_ace_rig_ave = HistTools::TransformEnergyAndDifferentialFluxNew(h_ace_ene_ave, isotope, "MeV/n cm", "GV m", "_rig"); // load averaged ACE data for the same element, converted in rigidity
	HistTools::SetStyle(h_ams_new, kBlue, kFullCircle, 1.4, 1, 1);
	HistTools::SetStyle(h_ace_rig_ave, kBlack, kFullCircle, 1.4, 1, 1);

	UShort_t namsbins = h_ams_new->GetNbinsX(); 
	UShort_t nacebins = h_ace_rig_ave->GetNbinsX(); 
	double R1 = h_ace_rig_ave->GetBinLowEdge(1);
	double R2 = h_ams->GetBinLowEdge(namsbins+1);

	int nnodes = 9; // combined template 
	int nnodes_ams = 7; 

	printf(" template = %s, %d \n", element, get_ams_data_value(element)); 

	// load combined fit template 
	TFile *file0 = new TFile(Form("data/ACE/compare/fit_%s_%dnodes.root", get_template(element), nnodes));  

	Spline *sp_comb = new Spline("sp_comb", nnodes, Spline::LogLog | Spline::PowerLaw);
	TF1 *fsp_comb = sp_comb->GetTF1Pointer();  
	TF1 *fit_comb = (TF1*) file0->Get(Form("fit_both")); 

	HistTools::CopyParameters(fit_comb, fsp_comb); // error 
	double x1, x2;
	fit_comb->GetRange(x1,x2);
	fsp_comb->SetRange(x1,x2);

	// load AMS He combined, modified spline  
	TFile *file1 = new TFile(Form("data/amsfit/fit_result_node%d.root", nnodes_ams));

	// load ACE BR 
	TFile *file2 = new TFile(Form("data/ACE/fill/%s_fill.root", element));
	TFile *file2a = new TFile(Form("data/ACE/fill/%s_fill2.root", element));  

	TGraphErrors *g_ratio_time[7]; 
	TGraphErrors *g_ratio_time_norm[7]; 

	for (int i=0; i<7; ++i){
		g_ratio_time[i] = (TGraphErrors*) file2a->Get(Form("g_ratio_time_%d", i)); 	
		HistTools::SetStyle(g_ratio_time[i], kBlue, kFullCircle, 1.4, 1, 1);
	}

	// load AMS He Model by Rescaled C Template
	TFile *file3 = new TFile(Form("data/ACE/extend2/fit_C_temp_he_%dnodes.root", nnodes));

	TCanvas *c0 = new TCanvas("c0", "", 1600, 900); // fluxes & flux/template
	c0->Divide(1, 2); 

	TCanvas *c1 = new TCanvas("c1", "", 1600, 900); // He ratio, F(R,t)/<F(R)> & spectral indices 
	c1->Divide(1, 2); 
	
	TCanvas *c2 = new TCanvas("c2", "", 1600, 900); // Modeling F(R,t)/<F(R)> 
							// fit/data of F(R,t)/<F(R)> 
	c2->Divide(1, 3); 				// spectral indices of F(R,t)/<F(R)>, here is defined by h_ratio (ACE) and h_ratio0 (AMS)

	TCanvas *c3 = new TCanvas("c3", "", 1600, 900); // plot chi2/ndf vs. BR for both ACE & AMS F(R,t)/<F(R)> fits  
	c3->Divide(2, 1);  

	TCanvas *c4 = new TCanvas("c4", "", 1600, 900); // plot for F(R,t)/<F(R)> fit residuals for each ACE bin 
	c4->Divide(1, 1); 

	TCanvas *c4_2 = new TCanvas("c4_2", "", 1600, 900); // c4 but for each AMS bin 
	c4_2->Divide(1, 1); 
	
	TCanvas *c5 = new TCanvas("c5", "", 1600, 900); // plot parameter [0] [1] vs. BR for the fit 
	c5->Divide(1, 3); 

	TH1 *h_a1 = new TH1D("", "", 3000, 0, 3000); 
	TH1 *h_a2 = new TH1D("", "", 3000, 0, 3000);
	TH1 *h_a3 = new TH1D("", "", 3000, 0, 3000);
	TH1 *h_a4 = new TH1D("", "", 3000, 0, 3000);	
	TH1 *h_a5 = new TH1D("", "", 3000, 0, 3000); 

	vector<double> last_pars;
	vector<double> last_pars2 = {0.5, 0.5};  

	TGraph *g_chi2_ace = new TGraph(); 
	TGraph *g_chi2_ams = new TGraph(); 
	TGraph *g_chi2_ace2 = new TGraph(); 
	TGraph *g_chi2_ams2 = new TGraph(); 
	TGraphErrors *g_pars[2];

	for (int ipar=0; ipar<2; ++ipar){
		g_pars[ipar] = new TGraphErrors(); 

	}

	TGraphErrors *g_ratio_ace[7]; 
	TGraphErrors *g_ratio_ace_fit[7];
	TGraphErrors *g_ratio_ams[h_ams_new->GetNbinsX()]; 
	TGraphErrors *g_ratio_ams_fit[h_ams_new->GetNbinsX()]; 

	TGraphErrors *g_ratio_ave = new TGraphErrors();
	TGraphErrors *g_ratio_ave_norm; // normalized averaged ratio of data/template-1 vs. time  

	TGraphErrors *g_residual_ace[7];
	TGraphErrors *g_resierr_ace[7]; 
	TGraphErrors *g_residual_norm_ace[7];
	TGraphErrors *g_residual_ams[h_ams_new->GetNbinsX()]; 
	TGraphErrors *g_residual_norm_ams[h_ams_new->GetNbinsX()];

	for (int bin=0; bin <= h_ace_rig_ave->GetNbinsX(); ++bin){ 
		if (bin%2==0){  
			g_ratio_ace[bin/2] = new TGraphErrors(nBRs); 
			g_ratio_ace_fit[bin/2] = new TGraphErrors(nBRs); 
			g_residual_ace[bin/2] = new TGraphErrors(nBRs); 
			g_resierr_ace[bin/2] = new TGraphErrors(nBRs);  
			HistTools::SetStyle(g_residual_ace[bin/2], kPink, kFullCircle, 1.4, 1, 1);
			HistTools::SetStyle(g_ratio_ave, kBlue+1, kFullCircle, 1.4, 1, 1);
		} 
	}  
 
	for (int bin=0; bin <= h_ams_new->GetNbinsX(); ++bin){
		g_ratio_ams[bin] = new TGraphErrors(nBRs);  
		g_ratio_ams_fit[bin] = new TGraphErrors(nBRs); 
		g_residual_ams[bin] = new TGraphErrors(nBRs);  
		HistTools::SetStyle(g_residual_ams[bin], kPink, kFullCircle, 1.4, 1, 1); 
	}    

	TFile *file_f5 = new TFile(Form("data/ACE/fill/F5_%s.root", element), "recreate"); 
 
	int iBR_true = 0;
	for (int iBR=0; iBR<nBRs; ++iBR){

		// AMS_He Monthly Spline 
		Spline *sp_he = new Spline("sp_he", nnodes_ams, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_he = sp_he->GetTF1Pointer();  
		TF1 *fit_he = (TF1*) file1->Get(Form("fsp_BR_he_%02d", iBR)); 

		HistTools::CopyParameters(fit_he, fsp_he); // error 
		fit_he->GetRange(x1,x2);
		fsp_he->SetRange(x1,x2);

		// AMS_He Integrated Spline  
		Spline *sp_he_ave = new Spline("sp_he_ave", nnodes_ams, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_he_ave = sp_he_ave->GetTF1Pointer();  
		TF1 *fit_he_ave = (TF1*) file1->Get("fsp_he"); 

		HistTools::CopyParameters(fit_he_ave, fsp_he_ave); // error 
		x1=0, x2=0;
		fit_he_ave->GetRange(x1,x2);
		fsp_he_ave->SetRange(x1,x2); 

		// AMS_He Model by Rescaled C ACE+AMS Combined Template  
		Spline *sp_he_temp = new Spline("sp_he_temp", nnodes, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_he_temp = sp_he_temp->GetTF1Pointer();  
		TF1 *fit_he_temp = (TF1*) file3->Get("fsp_he"); 

		HistTools::CopyParameters(fit_he_temp, fsp_he_temp); // error 
		x1=0, x2=0; 
		fit_he_temp->GetRange(x1,x2);
		fsp_he_temp->SetRange(x1,x2);  

		TF1 *fsp_he2 = HistTools::CombineTF1(fsp_he, fsp_he_ave, HistTools::Divide, "fsp_he2"); // He(R,t) Fit vs. <He(R)> Fit
		TF1 *fspind_he = HistTools::GetSpectralIndex(fsp_he2, "fspind_he", R1, R2); 
		fspind_he->SetRange(h_he_int->GetBinCenter(2), x2);

		// AMS_He(R,t)/<He(R)>
		TH1 *h_ratio0 = (TH1D*) h_BR_he[iBR]->Clone("h_ratio0"); 
		h_ratio0->Divide(h_he_int); 

		// refill AMS BR He to match the bins  
		HistTools::SetStyle(h_BR_he[iBR], kBlue, kFullCircle, 1.1, 1, 1); 

		TH1 *h_ace_BR = (TH1*) file2->Get(Form("h_rig_%s_BR%d", element, 2426+iBR_true)); 

		// rescaled <AMS_C> by AMS He BR / <AMS_He>   
		TH1 *h_ams_BR_fake = (TH1D*) h_ams_new->Clone("h_ams_BR_fake");

		h_ams_BR_fake->Multiply(h_BR_he[iBR]); 
		h_ams_BR_fake->Divide(h_he_int); 
 
		// h_ams_BR_fake->Print("range");  

		// rescaled combined fit 
		TH1 *h_ratio;
		if (get_index(element) <4) h_ratio = (TH1 *) h_ace_BR->Clone("h_ratio"); // ACE BR/Temp 	
		if (get_index(element)>=4) h_ratio = (TH1 *) h_ace_rig_ave->Clone("h_ratio"); // ACE Ave/Temp 	
		h_ratio->Divide(fsp_comb);  

		HistTools::SetStyle(h_ratio, kPink, kFullCircle, 1.4, 1, 1);

		double ratio_sum=0; // compute average of h_ratio manually  
		for(int nbin=0;nbin<14;++nbin){
			ratio_sum += h_ratio->GetBinContent(nbin);
			//printf("ratio_sum = %0.6f \n", ratio_sum);
		}
		double ratio_ave = ratio_sum/7;

		// printf("%0.6f \n", ratio_ave);
		//HistTools::PrintFunction(fsp_comb);
			
		double scale = 1./ratio_ave;
		if (get_index(element)<4) h_ratio->Scale(scale);	

		// rescaled combined fit 
		TF1 *rescaled_fit = HistTools::CombineTF1Const(fsp_comb, ratio_ave, HistTools::MultiplyConst, "rescaled_fit", R1, R2); 
		//if (get_index(element)<4) rescaled_fit = HistTools::CombineTF1(rescaled_fit, fsp_comb, HistTools::Divide, "fit_ratio", R1, R2);  

		TH1 *h_ratio1 = (TH1D*) h_ace_BR->Clone("h_ratio1"); // spind model 
		if (get_index(element) <4) h_ratio1->Divide(fsp_comb); 
		if (get_index(element)>=4) h_ratio1->Divide(rescaled_fit); 

		TH1 *h_ratio2 = (TH1D*) h_ams_BR_fake->Clone("h_ratio2"); // spind model 
		h_ratio2->Divide(h_ams_new);
   
		// h_ratio2->Print("range");

		// estimated AMS / resclaed combined fit
		TH1 *h_ratio3 = (TH1D *) h_ams_BR_fake->Clone("h_ratio3");		

		h_ratio3->Divide(rescaled_fit); 
		HistTools::SetStyle(h_ratio3, kPink, kFullCircle, 1.4, 1, 1);

		c0->cd(1);
		gPad->SetGrid(); 
		gPad->SetLogy();
		gPad->SetLogx(); 	

		TLegend *legend = new TLegend(0.62,0.75,0.9,0.9);
		legend->AddEntry(h_ams, Form("Integrated AMS %s Flux", element), "p"); 
		legend->AddEntry(h_ace_BR, Form("ACE %s BR Flux", element), "p"); 
		legend->AddEntry(h_ams_BR_fake, Form("Estimated AMS %s BR Flux", element), "p"); 
		legend->AddEntry(rescaled_fit, "rescaled combined template", "l"); 

		h_a1->SetTitle(Form("%s BR-%d Model Rigidity Spectrum; ; ", element, 2426+iBR_true));
		h_a1->SetXTitle(Unit::GetEnergyLabel("GV"));
  		h_a1->SetYTitle(Unit::GetDifferentialFluxLabel("GV m"));

		//if (get_temp_int(element)<4) h_a1->GetYaxis()->SetRangeUser(1e-3, 1e1); 
		// h_a1->GetXaxis()->SetRangeUser(0.01, 3000.); 

		TAxis *axis1 = h_a1->GetXaxis(); 
		axis1->SetLimits(0.7, 60.); 
	
		h_a1->Draw("E1X0"); 
	
		HistTools::SetStyle(h_ams, kPink-3, kFullCircle, 1.4, 1, 1); 
		HistTools::SetStyle(h_ace_BR, kRed, kFullCircle, 1.4, 1, 1);
		// HistTools::SetStyle(h_ams_BR_he, kBlue, kFullCircle, 1.4, 1, 1);
		HistTools::SetStyle(h_ams_BR_fake, kBlue, kFullCircle, 1.4, 1, 1);

		// fsp_he->Draw("SAME"); 
		rescaled_fit->Draw("SAME"); 
		h_ams->Draw("E1X0 SAME"); 
		h_ace_BR->Draw("E1X0 SAME"); 
		// h_ams_BR_he->Draw("E1X0 SAME"); 
		h_ams_BR_fake->Draw("E1X0 SAME"); 

		//fit->SetLineColor(kBlue); 
		//fit->Draw("SAME"); 
		legend->Draw("SAME");   

		c0->cd(2); 
		gPad->SetGrid(); 
		gPad->SetLogx(); 

		TAxis *axis2 = h_a2->GetXaxis(); 
		axis2->SetLimits(0.7, 60.); 

		//if (get_temp_int(element)<4) h_a2->GetYaxis()->SetRangeUser(0.5, 1.5); 
		h_a2->GetXaxis()->SetRangeUser(0.7, 60.); 

		h_a2->SetTitle("");
		h_a2->SetXTitle(Unit::GetEnergyLabel("GV"));
		h_a2->SetYTitle(Form("Estimated BR Fluxes / Rescaled %s Model", element)); 
		h_a2->SetTitleSize(0.3,"t"); 

		h_a2->Draw("E1X0"); 
		h_ratio->Draw("E1X0 SAME"); 
		h_ratio3->Draw("E1X0 SAME"); 

		if (iBR==0) { c0->Print(Form("./data/ACE/fill/fake_td_ams_v4/%s_flux_model.pdf(", element), "pdf"); } 
		if (iBR>0 && iBR<nBRs-1) { c0->Print(Form("./data/ACE/fill/fake_td_ams_v4/%s_flux_model.pdf", element), "pdf"); } 
		if (iBR==nBRs-1){
			c0->Print(Form("./data/ACE/fill/fake_td_ams_v4/%s_flux_model.pdf)", element), "pdf"); 
		} 

		// plot the ratio of AMS_He_BR / <AMS_He>
		c1->cd(1); 
		gPad->SetGrid(); 
		gPad->SetLogx(); 

		//if (get_temp_int(element)<4) h_a3->GetYaxis()->SetRangeUser(0.5, 1.5); 
		h_a3->GetXaxis()->SetLimits(0.7, 60);		

		TLegend *l1 = new TLegend(0.62,0.8,0.9,0.9); 
		l1->AddEntry(h_ratio0, "AMS He(R,t)/<He(R)>", "PL"); 	

		h_a3->Draw("E1X0"); 
		h_a3->SetXTitle(Unit::GetEnergyLabel("GV"));
		h_a3->SetTitle(Form("; ; BR-%d He(R,t)/<He(R)>", iBR_true+2426)); 
		HistTools::SetStyle(h_ratio0, kBlue, kFullCircle, 1.4, 1, 1);
		h_ratio0->Draw("E1X0 SAME"); 

		// plot spectral indices
		TGraphAsymmErrors *gspind_ace = HistTools::GetSpectralIndex(h_ace_BR, 5, 2); // get spectral indices 
		TGraphAsymmErrors *gspind_ams = HistTools::GetSpectralIndex(h_ams_BR_fake, 4, 1); 

		c1->cd(2); 
		gPad->SetGrid(); 
		gPad->SetLogx();  

		TLegend *l_both = new TLegend(0.62,0.8,0.9,1.); 
		l_both->AddEntry(gspind_ace, "ACE Spectral Indices", "PL");
		l_both->AddEntry(gspind_ams, "AMS Spectral Indices", "PL"); 		
	
		//if (get_temp_int(element)<4) gspind_ace->GetYaxis()->SetRangeUser(-3, 3); 
		gspind_ace->GetXaxis()->SetLimits(0.7, 60);

		HistTools::SetStyle(gspind_ace, kPink, kFullCircle, 1.4, 1, 1); 
		gspind_ace->SetTitle(Form("; ; ACE %s Spectral Indices BR-%d", element, 2426+iBR_true));
		gspind_ace->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV"));
		gspind_ace->Draw("APL"); 

		//h_a4->Draw("E1X0"); 
		HistTools::SetStyle(gspind_ams, kBlue, kFullCircle, 1.4, 1, 1); 
		gspind_ams->SetTitle(Form("; ; AMS %s Spectral Indices", element)); 
		gspind_ams->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV"));
		gspind_ams->GetXaxis()->SetTitleSize(1.5); 
		gspind_ams->Draw("PL SAME"); 
		// fspind->Draw("SAME"); 
		l_both->Draw("SAME"); 

		if (iBR==0) { c1->Print(Form("./data/ACE/fill/fake_td_ams_v4/fake_ams_spind_%s.pdf(", element), "pdf"); }
		if (iBR>0 && iBR<nBRs-1) { c1->Print(Form("./data/ACE/fill/fake_td_ams_v4/fake_ams_spind_%s.pdf", element), "pdf"); }
		if (iBR==nBRs-1){
			// gspind_ace->Print("range");
			// cout << " " << endl; 
			// gspind_ams->Print("range"); 
			c1->Print(Form("./data/ACE/fill/fake_td_ams_v4/fake_ams_spind_%s.pdf)", element), "pdf"); 
		} 

		// find spectral indices of h_ratio + h_ratio0 
		TGraphAsymmErrors *gspind_ace2 = HistTools::GetSpectralIndex(h_ratio1, 5, 2); // get spectral indices 
		TGraphAsymmErrors *gspind_ams2 = HistTools::GetSpectralIndex(h_ratio2, 4, 1); 

		HistTools::SetStyle(gspind_ace2, kPink, kFullCircle, 1.4, 1, 1); 
		HistTools::SetStyle(gspind_ams2, kBlue, kFullCircle, 1.4, 1, 1); 
		HistTools::SetStyle(fspind_he, kBlack, kFullCircle, 0.8, 1, 1);

		// fit to F(R,t)/<F(R)>
		TF1 *fit_ratio = new TF1("f_ratio", "1+[0]*exp(-[1]*x)", 0.7, 60);
		TF1 *fit_ratio2 = new TF1("f_ratio2", "1+[0]*exp(-[1]*log(x))", 0.7, 60);

		double F1 = h_ratio1->GetBinContent(h_ratio1->FindBin(R1)), F2 = h_ratio0->GetBinContent(h_ratio0->GetNbinsX()); 
		double A1 = log(abs(F1-1)), A2 = log(abs(F2-1)); 
		TF1 *f1 = new TF1("f1","abs(sin(x)/x)*sqrt(x)",0,1);
		double rho = f1->GetRandom();
		double N = f1->GetRandom(); 
		
		// double rho = (R2-R1)/(A1-A2)
		// double N = 0.5*(abs(F1-1)/exp(-R1/rho) + abs(F2-1)/exp(-R2/rho)); 
		
		last_pars = {N, rho}; 
		last_pars2 = {N, rho};

		// printf(" rho = %0.4f, N = %0.20f \n ", rho, N); 
		// printf(" rho = %0.4f, N = %0.20f, F1 = %0.4f, F2 = %0.4f, (F1-1)/exp(-R1/rho) = %0.6f, (F2-1)/exp(-R2/rho) = %0.6f \n", rho, N, F1, F2, (F1-1)/exp(-R1/rho), (F2-1)/exp(-R2/rho));

		for (int ipar=0; ipar<fit_ratio->GetNpar(); ++ipar){ 
			fit_ratio->SetParameter(ipar, last_pars[ipar]); 
			fit_ratio2->SetParameter(ipar, last_pars2[ipar]); 
		}   

		TObjArray data; 
		data.Add(h_ratio1);
		data.Add(h_ratio0);

		// fit AMS and ACE data at the same time
		vector<double> rigmin, rigmax, chi2norm(2);
		ROOT::Fit::Fitter fitter;
		//HistTools::PrintFunction(fit);
		FitTools::SetCommonFitterOptions(fitter);
		FitTools::FitCombinedData(data, fit_ratio, "I", rigmin, rigmax, chi2norm, fitter, 3); 

		vector<double> rigmin2, rigmax2, chi2norm2(2);
		ROOT::Fit::Fitter fitter2;
		FitTools::SetCommonFitterOptions(fitter2);
		FitTools::FitCombinedData(data, fit_ratio2, "I", rigmin2, rigmax2, chi2norm2, fitter2, 3); 
 
		for (int ipar=0; ipar<fit_ratio->GetNpar(); ++ipar){
			g_pars[ipar]->SetPoint(iBR, UBRToTime(iBR_true+2426), fit_ratio->GetParameter(ipar)); 
			g_pars[ipar]->SetPointError(iBR, 0, fit_ratio->GetParError(ipar));   
		}  

		// get parameters from par0 f_corr 
		int index = choose_BO(element); // bin index
		int shift = 1; 
		fit_ratio->SetParameters(f_corr[index]->Eval(h_ace_BR->GetBinContent(2*best_bin[index]+shift)), par1_X[index]->GetMean(2)); // F5 set par0 to the f_corr fit 
		fit_ratio->SetParError(0, g_correrr[index]->GetEY()[h_correrr[index]->FindBin(h_ace_BR->GetBinContent(2*best_bin[index]+shift))-1] );
		fit_ratio->SetParError(1, par1_X[index]->GetRMS(2));

		//printf("h_correrr bin # = %d \n", h_correrr[0]->GetBinContent(h_correrr[0]->FindBin(h_ace_BR->GetBinContent(best_bin[0]+2)))); 
		//h_ace_BR->Print("range");
		//h_correrr[0]->Print("range");    

		TF1 *fspind_ratio = HistTools::GetSpectralIndex(fit_ratio, "fspind_ratio", R1, R2); 
		TF1 *fspind_ratio2 = HistTools::GetSpectralIndex(fit_ratio2, "fspind_ratio2", R1, R2); 

		//h_ams_BR_fake->Print("range"); 
		c2->cd(1); 
		gPad->SetGrid(); 
		gPad->SetLogx();
		gPad->SetBottomMargin(0.01); 

		//if (get_temp_int(element)<4) h_a4->GetYaxis()->SetRangeUser(0.5, 2.2); 
		h_a4->GetXaxis()->SetLimits(0.7, 60); 
	
		TLegend *l2 = new TLegend(0.62,0.8,0.9,1.); 
		l2->AddEntry(h_ratio1, Form("ACE %s(R,t)/<%s(R)>", element, element), "PL");
		l2->AddEntry(h_ratio0, "AMS He(R,t)/<He(R)>", "PL");  
		l2->AddEntry(fit_ratio, Form("1+%6.3fexp(-%6.3fR)", fit_ratio->GetParameter(0), fit_ratio->GetParameter(1)), "L"); 
		l2->AddEntry(fit_ratio2, Form("1+%6.3fexp(-%6.3fln(R))", fit_ratio2->GetParameter(0), fit_ratio2->GetParameter(1)), "L"); 

		h_a4->Draw("E1X0"); 
		h_a4->SetXTitle(Unit::GetEnergyLabel("GV"));
		h_a4->SetTitle(Form("; ; BR-%d %s(R,t)/<%s(R)>", iBR_true+2426, element, element));
		// h_a4->SetTitleSize(0.5,"y"); 
		HistTools::SetStyle(h_ratio1, kPink, kFullCircle, 1.4, 1, 1);
		HistTools::SetStyle(h_ratio0, kBlue, kFullCircle, 1.4, 1, 1);
		fit_ratio2->SetLineColor(kGreen); 
		fit_ratio->Draw("SAME");
		fit_ratio2->Draw("SAME"); 
		h_ratio1->Draw("E1X0 SAME");
		h_ratio0->Draw("E1X0 SAME");
		l2->Draw("SAME"); 

		TH1D *h_fitres[2];
		TH1D *h_fitres2[2]; 
		TLegend *l_chi2 = new TLegend(0.62,0.8,0.9,1.0); 

		for (UShort_t i = 0; i < 2; ++i)
		{
  			TH1 *hist = HistTools::ToHist(data[i]);
			TH1 *h_fiterr = (TH1*) HistTools::GetFitError( (TH1*) data[i], fit_ratio, "_fiterr", false, false, true, 10, DBL_MAX, 0.68, &fitter);
   			h_fitres[i] = (TH1D *) HistTools::GetResiduals(hist, fit_ratio, "_fitres", false, true, true, 5, 1, 0.68, &fitter);
			HistTools::CopyStyle(hist, h_fitres[i]);

   			UShort_t ndf  = hist->GetNbinsX();
   			Double_t chi2 = hist->Chisquare(fit_ratio, "R"); 
   			// printf(" %d nodes ### %-34s   chi2/ndf=%6.2f/%-2u   chi2norm=%5.2f   prob.=%5.2f%%\n", nnodes, hist->GetTitle(), chi2, ndf, chi2norm[i], TMath::Prob(chi2, ndf)*1e2);

			l_chi2->AddEntry(h_fitres[i], Form("1+[0]*exp(-[1]*R) chi2/ndf=%6.2f/%-2u", chi2, ndf), "P");  

			if (i==0){ 
				g_chi2_ace->SetPoint(iBR, UBRToTime(iBR_true+2426), chi2/ndf);
				for (int bin=0; bin <= h_fitres[i]->GetNbinsX(); ++bin){ 
					if (bin%2==0){ 
						g_ratio_ace[bin/2]->SetPoint(iBR, UBRToTime(iBR_true+2426), h_ratio1->GetBinContent(bin+1)); 
						g_ratio_ace[bin/2]->SetPointError(iBR, 0, h_ratio1->GetBinError(bin+1)); 
						g_ratio_ace_fit[bin/2]->SetPoint(iBR, UBRToTime(iBR_true+2426), fit_ratio->Eval(h_ratio1->GetBinCenter(bin+1))); 
						g_ratio_ace_fit[bin/2]->SetPointError(iBR, 0, h_fiterr->GetBinError(bin+1)); 
						g_residual_ace[bin/2]->SetPoint(iBR, UBRToTime(iBR_true+2426), h_fitres[i]->GetBinContent(bin+1)); 
						g_residual_ace[bin/2]->SetPointError(iBR, 0, h_fitres[i]->GetBinError(bin+1)); 

						double df = h_fitres[i]->GetBinContent(bin+1)*hist->GetBinContent(bin+1)/hist->GetBinError(bin+1); 
						double ddf = df*sqrt(pow(h_fitres[i]->GetBinError(bin+1)/h_fitres[i]->GetBinContent(bin+1))+pow(h_fiterr->GetBinError(bin+1)/h_fiterr->GetBinContent(bin+1))); // d(data-fit) 
						g_resierr_ace[bin/2]->SetPoint(iBR, UBRToTime(iBR_true+2426), df); 
						g_resierr_ace[bin/2]->SetPointError(iBR, 0, ddf); 

						// printf("bin/2 = %d, iBR = %d, time = %10.4u, fitres = %10.4f \n", bin/2, iBR, UBRToTime(iBR_true+2426), h_fitres[0]->GetBinContent(bin+1)); 
					} 
				}   
			} 

			if (i==1){  
				g_chi2_ams->SetPoint(iBR, UBRToTime(iBR_true+2426), chi2/ndf);
				for (int bin=0; bin <= h_ams_new->GetNbinsX()-1; ++bin){ 
					g_ratio_ams[bin]->SetPoint(iBR, UBRToTime(iBR_true+2426), h_ratio0->GetBinContent(bin+1)); 
					g_ratio_ams[bin]->SetPointError(iBR, 0, h_ratio0->GetBinError(bin+1));
					g_ratio_ams_fit[bin]->SetPoint(iBR, UBRToTime(iBR_true+2426), fit_ratio->Eval(h_ratio0->GetBinCenter(bin+1))); 
					g_ratio_ams_fit[bin]->SetPointError(iBR, 0, h_fiterr->GetBinError(bin+1)); 
					g_residual_ams[bin]->SetPoint(iBR, UBRToTime(iBR_true+2426), h_fitres[i]->GetBinContent(bin+1));
					g_residual_ams[bin]->SetPointError(iBR, 0, h_fitres[i]->GetBinError(bin+1)); 
					// printf("bin = %d, iBR = %d, time = %10.4u, fitres = %10.4f \n", bin, iBR, UBRToTime(iBR_true+2426), h_fitres[1]->GetBinContent(bin+1)); 
				}  
			} 	
		}	

		for (UShort_t i = 0; i < 2; ++i)
		{
  			TH1 *hist = HistTools::ToHist(data[i]);
   			h_fitres2[i] = (TH1D *)HistTools::GetResiduals(hist, fit_ratio2, "_fitres2", false, true, true, 5, 1, 0.68, &fitter2);
			HistTools::CopyStyle(hist, h_fitres2[i]);

   			UShort_t ndf  = hist->GetNbinsX();
   			Double_t chi2 = chi2norm2[i]*ndf; 
   			// printf(" %d nodes ### %-34s   chi2/ndf=%6.2f/%-2u   chi2norm=%5.2f   prob.=%5.2f%%\n", nnodes, hist->GetTitle(), chi2, ndf, chi2norm[i], TMath::Prob(chi2, ndf)*1e2); 

			l_chi2->AddEntry(h_fitres2[i], Form("1+[0]*exp(-[1]*ln(R)) chi2/ndf=%6.2f/%-2u", chi2, ndf), "P"); 

			if (i==0){ 
				g_chi2_ace2->SetPoint(iBR, UBRToTime(iBR_true+2426), chi2);
			} 

			if (i==1){ 
				g_chi2_ams2->SetPoint(iBR, UBRToTime(iBR_true+2426), chi2);
			} 	
		}

		TH1 *h_comp1 = (TH1*) h_ace_BR->Clone("h_comp1");
		h_comp1->Divide(rescaled_fit); 

		double comp1_sum=0, dcomp1_sum=0; // compute average ratio of data/template-1 of ACE F(R,t)/<F(R)> manually 		

		for(int nbin=0;nbin<h_ace_rig_ave->GetNbinsX();++nbin){
			comp1_sum += h_comp1->GetBinContent(nbin); 
			dcomp1_sum += h_comp1->GetBinError(nbin)*h_comp1->GetBinError(nbin); 
			//printf("comp1_sum = %0.6f \n", comp1_sum);
		}
		double comp1_ave = comp1_sum/7; 
		double dcomp1_ave = sqrt(dcomp1_sum)/7; 

		g_ratio_ave->SetPoint(iBR, UBRToTime(iBR_true+2426), comp1_ave-1); 
		g_ratio_ave->SetPointError(iBR, 0, dcomp1_ave); 

		c2->cd(2);
		gPad->SetGrid();
		gPad->SetLogx(); 
		gPad->SetTopMargin(0.01); 

		//if (get_temp_int(element)<4) h_a5->GetYaxis()->SetRangeUser(-0.15, 0.15)  ; 
		h_a5->GetXaxis()->SetLimits(0.7, 60); 

		h_a5->Draw("E1X0"); 
		h_a5->SetXTitle(Unit::GetEnergyLabel("GV"));
		h_a5->SetTitle("; ; Data/Fit-1"); 

		h_fitres2[0]->Draw("E1X0 SAME"); 
		h_fitres2[1]->Draw("E1X0 SAME");

		h_fitres[0]->Draw("E1X0 SAME");
		h_fitres[1]->Draw("E1X0 SAME");

		HistTools::SetStyle(h_fitres2[0], kGreen-1, kFullCircle, 1.4, 1, 1);
		HistTools::SetStyle(h_fitres2[1], kGreen-1, kFullCircle, 1.4, 1, 1); 

		l_chi2->Draw("SAME"); 		

		c2->cd(3);
		gPad->SetGrid();
		gPad->SetLogx(); 

		TLegend *l_both3 = new TLegend(0.62,0.7,0.8,1.0); 
		l_both3->AddEntry(gspind_ace2, Form("ACE %s(R,t)/<%s(R)> Spectral Indices", element, element), "PL");
		l_both3->AddEntry(gspind_ams2, Form("AMS %s(R,t)/<%s(R)> Spectral Indices", element, element), "PL"); 
		// l_both2->AddEntry(fspind_he, "AMS He(R,t)_Fit/temp<He(R)>_Fit Spectral Indices", "L");  	
		l_both3->AddEntry(fspind_ratio, "Spectral Index of 1+[0]*exp(-[1]*R) Fit to He(R,t)/<He(R)>", "L");	
		l_both3->AddEntry(fspind_he, "AMS He(R,t)_Fit/<He(R)>_Fit Spectral Indices", "L");  
		l_both3->AddEntry(fspind_ratio2, "Spectral Index of 1+[0]*exp(-[1]*ln(R)) Fit to He(R,t)/<He(R)>", "L"); 

		//if (get_temp_int(element)<4) gspind_ace2->GetYaxis()->SetRangeUser(-2, 1.3); 
		TAxis *axis5 = gspind_ace2->GetXaxis(); 
		axis5->SetLimits(0.7, 60.); 

		gspind_ace2->SetTitle(Form("; ; Spectral Indices Model of %s F(R,t)/<F(R)> at BR-%d", element, 2426+iBR_true));
		gspind_ace2->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV"));
		gspind_ace2->Draw("APL");  
		fspind_ratio->Draw("SAME"); 
		fspind_ratio2->Draw("SAME"); 
		fspind_ratio2->SetLineColor(kGreen); 
		fspind_he->Draw("SAME");
		gspind_ams2->Draw("PL SAME"); 
		gspind_ace2->Draw("PL SAME");
		l_both3->Draw("SAME"); 

		// break; 

		rescaled_fit->Write(Form("rescaled_fit_BR%d", 2426+iBR_true)); 
		fit_he->Write(Form("fit_he_BR%d", 2426+iBR_true)); 
		//fsp_he2->Write(Form("fsp_he2_BR%d", 2426+iBR_true)); 
		fit_ratio->Write(Form("fit_ratio_BR%d", 2426+iBR_true)); 
		fit_he_ave->Write(Form("fit_he_ave_BR%d", 2426+iBR_true)); 
		h_ratio1->Write(Form("h_ratio1_BR%d", 2426+iBR_true));
		h_ratio0->Write(Form("h_ratio0_BR%d", 2426+iBR_true)); 
		h_fitres[0]->Write(Form("h_fitres_ace_BR%d", 2426+iBR_true));
		h_fitres[1]->Write(Form("h_fitres_ams_BR%d", 2426+iBR_true));

		if (iBR==0) c2->Print(Form("./data/ACE/fill/fake_td_ams_v4/spind_model_%s_temp.pdf(", element), "pdf"); 
		if (iBR>0 && iBR<nBRs-1) c2->Print(Form("./data/ACE/fill/fake_td_ams_v4/spind_model_%s_temp.pdf", element), "pdf"); 
		if (iBR==nBRs-1){ 
			//gspind_ace2->Print("range");
			cout << " " << endl; 
			//gspind_ams2->Print("range");
			cout << " " << endl; 
			// h_ratio2->Print("range");
			c2->Print(Form("./data/ACE/fill/fake_td_ams_v4/spind_model_%s_temp.pdf)", element), "pdf"); 
		}

	   	if (iBR+2426==2472-1) iBR_true += 3; 
		else iBR_true ++; 
	}
			
	c3->cd(1);
	gPad->SetGrid(); 
	HistTools::SetStyle(g_chi2_ace, kPink, kFullCircle, 1.4, 1, 1); 
	HistTools::SetStyle(g_chi2_ace2, kBlue, kFullCircle, 1.4, 1, 1); 
	TLegend *l_c3_1 = new TLegend(0.62,0.8,0.9,0.9); 
	l_c3_1->AddEntry(g_chi2_ace, "1+[0]*exp(-[1]*R)", "PL"); 
	l_c3_1->AddEntry(g_chi2_ace2, "1+[0]*exp(-[1]*ln(R))", "PL"); 
	g_chi2_ace->SetTitle(Form("; ; ACE %s(R,t)/<%s(R)> Fitting Chi2/NDF", element, element)); 
	g_chi2_ace->GetXaxis()->SetTimeDisplay(1);
	g_chi2_ace->GetXaxis()->SetTimeFormat("%m-%y"); 
	g_chi2_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	g_chi2_ace->GetXaxis()->SetTitleSize(0.7); 
	// //if (get_temp_int(element)<4) g_chi2_ace->GetYaxis()->SetRangeUser(0., 4.5); 
	g_chi2_ace->Draw("APL"); 
	g_chi2_ace2->Draw("PL SAME"); 
	l_c3_1->Draw("SAME");

	c3->cd(2); 
	gPad->SetGrid(); 
	HistTools::SetStyle(g_chi2_ams, kPink, kFullCircle, 1.4, 1, 1); 
	HistTools::SetStyle(g_chi2_ams2, kBlue, kFullCircle, 1.4, 1, 1); 
	TLegend *l_c3_2 = new TLegend(0.62,0.8,0.9,0.9); 
	l_c3_2->AddEntry(g_chi2_ams, "1+[0]*exp(-[1]*R)", "PL");
	l_c3_2->AddEntry(g_chi2_ams2, "1+[0]*exp(-[1]*ln(R))", "PL");
	g_chi2_ams->SetTitle(Form("; ; AMS %s(R,t)/<%s(R)> Fitting Chi2/NDF", element, element)); 
	g_chi2_ams->GetXaxis()->SetTimeDisplay(1);
	g_chi2_ams->GetXaxis()->SetTimeFormat("%m-%y");
	g_chi2_ams->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	g_chi2_ams->GetXaxis()->SetTitleSize(0.7);
	// //if (get_temp_int(element)<4) g_chi2_ams->GetYaxis()->SetRangeUser(0., 2.2); 
	g_chi2_ams->Draw("APL");
	g_chi2_ams2->Draw("PL SAME");  
	l_c3_2->Draw("SAME"); 

	c3->Print(Form("./data/ACE/fill/fake_td_ams_v4/chi2_vs_BR_%s.png", element)); 

	g_ratio_ave_norm = get_norm_graph( g_ratio_ave ); 

	for (int iBR=0; iBR<nBRs; ++iBR){
				Double_t x=0, y=0; 
				g_ratio_ave_norm->GetPoint(iBR, x, y); 
				g_ratio_ave_norm->SetPoint(iBR, x, y-1); 			 
	}	

	TGraphErrors *g_ave_residual_ace_rig = new TGraphErrors(); // time-averaged residual vs. rigidity
	TGraphErrors *g_ave_residual_ams_rig = new TGraphErrors(); // time-averaged residual vs. rigidity 

	for (int bin=0; bin <= h_ace_rig_ave->GetNbinsX(); ++bin){
 
		if (bin%2==0) {  

			g_ratio_time_norm[bin/2] = get_norm_graph( g_ratio_time[bin/2] ); 
		
			for (int iBR=0; iBR<g_ratio_time_norm[bin/2]->GetN(); ++iBR){
				Double_t x=0, y=0; 
				g_ratio_time_norm[bin/2]->GetPoint(iBR, x, y); 
				g_ratio_time_norm[bin/2]->SetPoint(iBR, x, y-1); 			 
			}

			g_ave_residual_ace_rig->SetPoint(bin/2, h_ace_rig_ave->GetBinCenter(bin+1), g_residual_ace[bin/2]->GetMean(2) ); 
			g_ave_residual_ace_rig->SetPointError(bin/2, 0, g_residual_ace[bin/2]->GetRMS(2) ); 

			c4->cd(1);
			gPad->SetGrid(); 

			TLegend *l_comp1 = new TLegend(0.62, 0.8, 0.9, 1.0); 
			l_comp1->AddEntry(g_ratio_ave_norm, "Normalized Averaged ACE Data/Template-1", "PL"); 
			l_comp1->AddEntry(g_residual_ace[bin/2], Form("ACE Fitting Residuals (%0.4f GV)", h_ace_rig_ave->GetBinCenter(bin+1)), "PL"); 

			g_residual_ace[bin/2]->SetTitle(Form(" ; ; ACE %s Data/Template-1 vs. Fitting Residuals", element));  
			g_residual_ace[bin/2]->GetXaxis()->SetTimeDisplay(1);
			g_residual_ace[bin/2]->GetXaxis()->SetTimeFormat("%m-%y");
			g_residual_ace[bin/2]->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g_residual_ace[bin/2]->GetXaxis()->SetTitleSize(0.7);
			//if (get_temp_int(element)<4) g_residual_ace[bin/2]->GetYaxis()->SetRangeUser(-0.25, 0.3);
			g_residual_ace[bin/2]->Draw("APL");
			// g_ratio_ave_norm->Draw("PL SAME"); 
			g_ratio_time_norm[bin/2]->Draw("PL SAME"); 
			l_comp1->Draw("SAME"); 		

			if (bin==0) c4->Print(Form("./data/ACE/fill/fake_td_ams_v4/fitres_ace_vs_BR_%s.pdf(", element), "pdf"); 
			if (bin>0 && bin<h_ace_rig_ave->GetNbinsX()-1) c4->Print(Form("./data/ACE/fill/fake_td_ams_v4/fitres_ace_vs_BR_%s.pdf", element), "pdf");  
			if (bin==h_ace_rig_ave->GetNbinsX()-1){
				c4->Print(Form("./data/ACE/fill/fake_td_ams_v4/fitres_ace_vs_BR_%s.pdf)", element), "pdf"); 
			} 
		
			g_ratio_ace[bin/2]->Write(Form("g_ratio_ace_%d", bin/2)); 
			g_ratio_ace_fit[bin/2]->Write(Form("g_ratio_ace_fit_%d", bin/2)); 
			g_residual_ace[bin/2]->Write(Form("g_residual_ace_%d", bin/2)); 
			g_resierr_ace[bin/2]->Write(Form("g_resierr_ace_%d", bin/2)); 
		} 
	} 

	for (int bin=0; bin <= h_ams_new->GetNbinsX()-1; ++bin){ 

			g_ave_residual_ams_rig->SetPoint(bin, h_ams_new->GetBinCenter(bin+1), g_residual_ams[bin]->GetMean(2) ); 
			g_ave_residual_ams_rig->SetPointError(bin, 0, g_residual_ams[bin]->GetRMS(2) ); 

			c4_2->cd(1); 
			gPad->SetGrid(); 

			g_residual_ams[bin]->SetTitle(Form("AMS %s(R,t)/<%s(R)>; ; Fitting Residuals (%0.4f GV)", element, element, h_ams_new->GetBinCenter(bin+1)));  
			g_residual_ams[bin]->GetXaxis()->SetTimeDisplay(1);
			g_residual_ams[bin]->GetXaxis()->SetTimeFormat("%m-%y");
			g_residual_ams[bin]->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g_residual_ams[bin]->GetXaxis()->SetTitleSize(0.7);
			//if (get_temp_int(element)<4) g_residual_ams[bin]->GetYaxis()->SetRangeUser(-0.06, 0.1);  
			// g_residual_ams[bin]->Print("range"); 
			g_residual_ams[bin]->Draw("APL"); 

			if (bin==0) c4_2->Print(Form("./data/ACE/fill/fake_td_ams_v4/fitres_ams_vs_BR_%s.pdf(", element), "pdf"); 
			if (bin>0 && bin<h_ams_new->GetNbinsX()-1) c4_2->Print(Form("./data/ACE/fill/fake_td_ams_v4/fitres_ams_vs_BR_%s.pdf", element), "pdf");  
			if (bin==h_ams_new->GetNbinsX()-1){
				c4_2->Print(Form("./data/ACE/fill/fake_td_ams_v4/fitres_ams_vs_BR_%s.pdf)", element), "pdf"); 
			} 

			g_ratio_ams[bin]->Write(Form("g_ratio_ams_%d", bin)); 
			g_ratio_ams_fit[bin]->Write(Form("g_ratio_ams_fit_%d", bin)); 
			g_residual_ams[bin]->Write(Form("g_residual_ams_%d", bin)); 
	} 

	for (int ipar=0; ipar<2; ++ipar){
		c5->cd(ipar+1); 
		gPad->SetGrid(); 
		g_pars[ipar]->Draw("APL"); 
		HistTools::SetStyle(g_pars[ipar], kPink, kFullCircle, 1.4, 1, 1); 
		g_pars[ipar]->GetXaxis()->SetTimeDisplay(1);
		g_pars[ipar]->GetXaxis()->SetTimeFormat("%m-%y");
		g_pars[ipar]->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
		g_pars[ipar]->GetXaxis()->SetTitleSize(0.7);
		//if (get_temp_int(element)<4) g_pars[ipar]->GetYaxis()->SetRangeUser(-3., 3.); 
		if (ipar==0) g_pars[ipar]->SetTitle(Form("1+[0]*exp(-[1]*R); ; Parameter [%d]", ipar));
		else g_pars[ipar]->SetTitle(Form(" ; ; Parameter [%d]", ipar));
	}		

	TGraphErrors *g_pars_ratio = new TGraphErrors(); 
	for (int i=0; i<g_pars[0]->GetN(); ++i){
		g_pars_ratio->SetPoint(i, g_pars[0]->GetX()[i], g_pars[1]->GetY()[i]/g_pars[0]->GetY()[i]); 
		g_pars_ratio->SetPointError(i, 0, g_pars[1]->GetY()[i]/g_pars[0]->GetY()[i]*sqrt(pow(g_pars[1]->GetEY()[i]/g_pars[1]->GetY()[i], 2)+pow(g_pars[0]->GetEY()[i]/g_pars[0]->GetY()[i], 2))); 
	}

	c5->cd(3);
	gPad->SetGrid();
	HistTools::SetStyle(g_pars_ratio, kPink, kFullCircle, 1.4, 1, 1);
	g_pars_ratio->GetXaxis()->SetTimeDisplay(1);
	g_pars_ratio->GetXaxis()->SetTimeFormat("%m-%y");
	g_pars_ratio->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
	g_pars_ratio->GetXaxis()->SetTitleSize(0.7);
	//if (get_temp_int(element)<4) g_pars_ratio->GetYaxis()->SetRangeUser(-15., 30.); 
	g_pars_ratio->SetTitle(Form(" ; ; [%d]/[%d]", 1, 0)); 
	g_pars_ratio->Draw("AP"); 

	// g_pars_ratio->Print("range"); 

	c5->Print(Form("./data/ACE/fill/fake_td_ams_v4/pars_vs_BR_%s.png", element)); 

	TCanvas *c6 = new TCanvas("c6", "", 1600, 900); // plot time-averaged fit residuals to F(R,t)/<F(R)> vs. rigidity 
	c6->Divide(1, 1);  

	c6->cd(1); 
	gPad->SetGrid(); 
	gPad->SetLogx(); 

	HistTools::CopyStyle( h_ace_rig_ave, g_ave_residual_ace_rig ); 
	HistTools::CopyStyle( h_ams_new, g_ave_residual_ams_rig );
	////if (get_temp_int(element)<4) g_ave_residual_ace_rig->GetYaxis()->SetRangeUser(-0.15, 0.15); 
	g_ave_residual_ace_rig->GetXaxis()->SetLimits(0.7, 60);
	g_ave_residual_ace_rig->SetTitle(Form(" ; ; Time-averaged Fit Residuals of %s(R,t)/<%s(R)>", element, element)); 
	g_ave_residual_ace_rig->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV")); 

	TLegend *l_both4 = new TLegend(0.62,0.8,0.9,0.9); 
	l_both4->AddEntry(g_ave_residual_ace_rig, "ACE", "PL");
	l_both4->AddEntry(g_ave_residual_ams_rig, "AMS", "PL"); 

	g_ave_residual_ace_rig->Draw("APL");
	g_ave_residual_ams_rig->Draw("PLSAME"); 
	l_both4->Draw("SAME"); 

	g_chi2_ace->Write("g_chi2_ace"); 
	g_chi2_ams->Write("g_chi2_ams"); 

	g_ave_residual_ace_rig->Write("g_ave_residual_ace_rig"); 
	g_ave_residual_ams_rig->Write("g_ave_residual_ams_rig"); 

	file_f5->Close(); 

	c6->Print(Form("./data/ACE/fill/fake_td_ams_v4/time_ave_residual_vs_rig_%s.png", element)); 

	return 0; 
} 

void compare_fake_flux( const char *element, Particle::Type isotope ){

	int nnodes = 7; 
 	Experiments::DataPath = "data";	
	TFile *file = new TFile(Form("data/ACE/compare/fit_%s_%dnodes.root", get_template(element), nnodes)); 
	
	gStyle->SetOptStat(0);
	const int nBRs = Experiments::Info[Experiments::AMS02].Dataset[1].nMeasurements; // read the number of BRs
	TH1D **h_BR_he = Experiments::GetDatasetHistograms(Experiments::AMS02, 4);

   	int data_value[n_ele] = {0, 18, 23, 27, 31, 20, 43, 22}; // p, He, Li, Be, B, C, N, O

	double *kin_bins = get_kin_bins(element);
        double *SpallCorr = get_spall_corr(element);
	double *SpallCorrUnc = get_spall_corr_unc(element);
	double *EMed = get_EMed(element);

	const UInt_t FirstACEBR = 2240;
	vector<UInt_t> BRs;
	// we stop at BR 2493, which ends on 2016/05/23, just 3 days before the end of the data taking period for AMS nuclei
	for (UInt_t br=2426; br<=2493; ++br) { 
		if (br != 2472 && br != 2473) BRs.push_back(br-FirstACEBR); 
	}

	TH1 *h_ams = Experiments::GetMeasurementHistogram(Experiments::AMS02, data_value[1], 0); // load AMS data for a given element
	TH1 *h_ace_ene = HistTools::GraphToHist(get_ace_average_graph( element , &BRs[0], BRs.size() ), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
	TH1 *h_ace = HistTools::TransformEnergyAndDifferentialFluxNew(h_ace_ene, isotope, "MeV/n cm", "GV m", "_rig"); // load averaged ACE data for the same element, converted in rigidity

	UShort_t namsbins = h_ams->GetNbinsX();
	UShort_t nacebins = h_ace->GetNbinsX();
	double R1 = h_ace->GetBinLowEdge(1);
	double R2 = h_ams->GetBinLowEdge(namsbins+1);

	TFile *file_f1 = new TFile(Form("data/ACE/fill/F1_%s.root", element)); 
	TFile *file_f2 = new TFile(Form("data/ACE/fill/F2_%s.root", element)); 
	TFile *file_f3 = new TFile(Form("data/ACE/fill/F3_%s.root", element)); 
	TFile *file_f4 = new TFile(Form("data/ACE/fill/F4_%s.root", element));
	TFile *file_f5 = new TFile(Form("data/ACE/fill/F5_%s.root", element)); 

	//TCanvas *c0 = new TCanvas("c0", "", 1600, 900); 
	//c0->Divide(1, 2);

	int i_skip = 5; 

	TCanvas *c1 = new TCanvas("c1", "", 1600, 900); // P1
	c1->Divide(5, 2); 

	TCanvas *c2 = new TCanvas("c2", "", 1600, 900); // P2-ACE
	c2->Divide(5, 2); 

	TCanvas *c3 = new TCanvas("c3", "", 1600, 900); // P2-AMS
	c3->Divide(5, 2); 

	TCanvas *c4 = new TCanvas("c4", "", 1600, 900); // P3
	c4->Divide(5, 1); 

	TCanvas *c5 = new TCanvas("c5", "", 1600, 900); // P4
	c5->Divide(5, 2); 


	// P1F123
	int iBR_true=0; 
	for (int iBR=0; iBR<nBRs; ++iBR){

	   //F1
	   for (int k=0; k<1; ++k){

		TH1* h_rig = (TH1*) file_f1->Get(Form("h_rig_%d", 2426+iBR_true));

		TH1 *h_a = new TH1D("", "", 3000, 0, 3000); 
		TH1 *h_b = new TH1D("", "", 3000, 0, 3000); 

		// average of h_rig 
		double rig_sum=0; // compute average of h_rig manually  
		for(int nbin=0;nbin<14;++nbin){
			rig_sum += h_rig->GetBinContent(nbin);
			//printf("ratio_sum = %0.6f \n", ratio_sum);
		}
		double rig_ave = rig_sum/7; 

		// rescale ACE BR to ACE Averaged Magnitude 
		Spline *sp_comb = new Spline("sp_comb", nnodes, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_comb = sp_comb->GetTF1Pointer();  
		TF1 *fit_comb = (TF1*) file->Get("fit_both"); 

		HistTools::CopyParameters(fit_comb, fsp_comb); // error 
		double x1, x2;
		fit_comb->GetRange(x1,x2);
		fsp_comb->SetRange(x1,x2);	

		TH1 *h_ratio1 = (TH1 *) h_rig->Clone("h_ratio1");

		h_ratio1->Divide(fsp_comb);

		HistTools::SetStyle(h_rig, kBlack, kFullCircle, 0.9, 1, 1);
		HistTools::SetStyle(h_ratio1, kBlack, kFullCircle, 0.9, 1, 1);

		double ratio_sum=0; // compute average of h_ratio1 manually  
		for(int nbin=0;nbin<14;++nbin){
			ratio_sum += h_ratio1->GetBinContent(nbin);
			//printf("ratio_sum = %0.6f \n", ratio_sum);
		}
		double ratio_ave = ratio_sum/7;

		// printf("%d %0.6f \n", (UInt_t) utime, ratio_ave);
		//HistTools::PrintFunction(fit_comb);
			
		double scale = 1./ratio_ave;	

		TH1 *h_rig_rescale = (TH1 *) h_rig->Clone("h_rig_rescale");
		h_rig_rescale->Divide(fsp_comb); 
		HistTools::SetStyle(h_rig_rescale, kBlack, kFullCircle, 0.8, 1, 1);

		// TF1 *rescaled_fit = HistTools::CombineTF1Const(fsp_comb, ratio_ave, HistTools::MultiplyConst, "rescaled_fit", R1, R2); 

		TF1 *fit_ratio = HistTools::CombineTF1Const(fsp_comb, ratio_ave, HistTools::MultiplyConst, "rescaled_fit", R1, R2); 
		fit_ratio = HistTools::CombineTF1(fit_ratio, fsp_comb, HistTools::Divide, "fit_ratio", R1, R2);  

		// rescaled_fit->Print(); 

		TObjArray data; 
		data.Add(h_rig_rescale);
		// data.Add(h_ratio0);

		// fit AMS and ACE data at the same time
		vector<double> rigmin, rigmax, chi2norm(2);
		ROOT::Fit::Fitter fitter;
		//HistTools::PrintFunction(fit);
		FitTools::SetCommonFitterOptions(fitter);
		FitTools::FitCombinedData(data, fit_ratio, "I", rigmin, rigmax, chi2norm, fitter, 3); 

		TH1D *h_fitres[2]; 
		h_fitres[0] = (TH1D *) file_f1->Get(Form("h_fitres_ace_BR%d", 2426+iBR_true)); 

		c1->cd(1);
		gPad->SetGrid(); 
		gPad->SetLogx();
		gPad->SetMargin(0.12, 0.04, 0.1, 0.08);	

		h_a->SetTitle(Form("F1; ; BR-%d %s(R,t)/<%s(R)>", 2426+iBR_true, element, element));
		h_a->SetXTitle(Unit::GetEnergyLabel("GV"));

		h_a->GetYaxis()->SetRangeUser(-0.3, 3.); 
		h_a->GetXaxis()->SetLimits(0.7, 60); 	

		TLegend *l2 = new TLegend(0.3,0.7,0.8,0.9);
		l2->AddEntry(fit_ratio, Form("averaged ratio=%0.4f", ratio_ave), "l"); 
		l2->AddEntry(h_rig, Form("ACE %s(R,t)/<%s(R)>", element, element), "p"); 
	
		h_a->Draw("E1X0"); 
		//rescaled_fit->SetRange(0., 3.); 
		fit_ratio->Draw("SAME"); 
		h_rig_rescale->Draw("E1X0 SAME"); 
		l2->Draw("SAME"); 

		c1->cd(6); 
		gPad->SetGrid(); 
		gPad->SetLogx();
		gPad->SetMargin(0.12, 0.04, 0.2, 0.01);

		h_b->GetYaxis()->SetRangeUser(-0.5, 0.5)  ; 
		h_b->GetXaxis()->SetLimits(0.7, 60); 

		h_b->SetTitle("");
		h_b->SetXTitle(Unit::GetEnergyLabel("GV"));
		h_b->SetYTitle("Data/Model-1"); 
		//h_ratio1->GetXaxis()->SetRangeUser(0, 3.); 

		h_b->Draw("E1X0");
		HistTools::SetStyle(h_fitres[0], kBlack, kFullCircle, 0.8, 1, 1);
		h_fitres[0]->Draw("E1X0 SAME"); 
		// h_fitres[0]->Print("range"); 
	   }
	 
	   //F2
	   for (int k=0; k<1; ++k){

		TGraphErrors *g_chi2_ace = (TGraphErrors *) file_f2->Get("g_chi2_ace"); 
		TGraphErrors *g_chi2_ams = (TGraphErrors *) file_f2->Get("g_chi2_ams");

		TH1 *h_a = new TH1D("", "", 3000, 0, 3000); 
		TH1 *h_b = new TH1D("", "", 3000, 0, 3000); 

		TF1 *fit_ratio = (TF1*) file_f2->Get(Form("fit_ratio_BR%d", 2426+iBR_true)); 
		TH1 *h_ratio0 = (TH1*) file_f2->Get(Form("h_ratio0_BR%d", 2426+iBR_true)); 
		TH1 *h_ratio1 = (TH1*) file_f2->Get(Form("h_ratio1_BR%d", 2426+iBR_true)); 

		TH1 *h_fitres[2]; 
		h_fitres[0] = (TH1*) file_f2->Get(Form("h_fitres_ace_BR%d", 2426+iBR_true)); 
		h_fitres[1] = (TH1*) file_f2->Get(Form("h_fitres_ams_BR%d", 2426+iBR_true));

		c1->cd(2); 
		gPad->SetGrid(); 
		gPad->SetLogx();
		gPad->SetMargin(0.12, 0.04, 0.1, 0.08);	

		h_a->GetYaxis()->SetRangeUser(0, 3.5); 
		h_a->GetXaxis()->SetLimits(0.7, 60); 
	
		TLegend *l2 = new TLegend(0.3,0.7,0.8,0.9); 
		l2->AddEntry(h_ratio1, Form("ACE %s(R,t)/<%s(R)>", element, element), "PL");
		l2->AddEntry(h_ratio0, "AMS He(R,t)/<He(R)>", "PL");  
		l2->AddEntry(fit_ratio, Form("1+%0.2fpm%0.2fexp(-%0.2fpm%0.2fR)", fit_ratio->GetParameter(0), fit_ratio->GetParError(0), fit_ratio->GetParameter(1), fit_ratio->GetParError(1)), "L"); 
		// l2->AddEntry(fit_ratio2, Form("1+%6.3fexp(-%6.3fln(R))", fit_ratio2->GetParameter(0), fit_ratio2->GetParameter(1)), "L"); 

		h_a->Draw("E1X0"); 
		h_a->SetTitle(Form("F2; ; BR-%d %s(R,t)/<%s(R)>", iBR_true+2426, element, element));
		h_a->SetXTitle(Unit::GetEnergyLabel("GV"));
		// h_a->SetTitleSize(0.5,"y"); 
		HistTools::SetStyle(h_ratio1, kBlack, kFullCircle, 0.8, 1, 1);
		HistTools::SetStyle(h_ratio0, kBlack, kFullSquare, 0.8, 1, 1);
		fit_ratio->Draw("SAME");
		h_ratio1->Draw("E1X0 SAME"); 
		h_ratio0->Draw("E1X0 SAME");
		l2->Draw("SAME"); 

		c1->cd(7);
		gPad->SetGrid();
		gPad->SetLogx();  
		gPad->SetMargin(0.12, 0.04, 0.2, 0.01);

		h_b->GetYaxis()->SetRangeUser(-0.5, 0.5)  ; 
		h_b->GetXaxis()->SetLimits(0.7, 60); 

		TLegend *l_chi2 = new TLegend(0.62,0.8,0.9,0.9); 
		l_chi2->AddEntry(h_ratio1, Form("chi2=%0.4f", g_chi2_ace->GetY()[iBR]), "PL");
		l_chi2->AddEntry(h_ratio0, Form("chi2=%0.4f", g_chi2_ams->GetY()[iBR]), "PL");  

		h_b->Draw("E1X0"); 
		h_b->SetTitle("; ; Data/Fit-1"); 
		h_b->SetXTitle(Unit::GetEnergyLabel("GV"));

		HistTools::SetStyle(h_fitres[0], kBlack, kFullCircle, 0.8, 1, 1);
		HistTools::SetStyle(h_fitres[1], kBlack, kFullSquare, 0.8, 1, 1);

		h_fitres[0]->Draw("E1X0 SAME");
		h_fitres[1]->Draw("E1X0 SAME");
		l_chi2->Draw("SAME"); 
	   }

	   //F3
	   for (int k=0; k<1; ++k){

		int nnodes_ams = 7; 

		TGraphErrors *g_chi2_ace = (TGraphErrors *) file_f3->Get("g_chi2_ace"); 
		TGraphErrors *g_chi2_ams = (TGraphErrors *) file_f3->Get("g_chi2_ams"); 

		TH1 *h_a = new TH1D("", "", 3000, 0, 3000); 
		TH1 *h_b = new TH1D("", "", 3000, 0, 3000); 

		Spline *sp_he = new Spline("sp_he", nnodes_ams, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_he = sp_he->GetTF1Pointer();  
		TF1 *fit_he = (TF1*) file_f3->Get(Form("fit_he_BR%d", 2426+iBR_true)); 

		HistTools::CopyParameters(fit_he, fsp_he);
		double x1, x2;
		fit_he->GetRange(x1,x2);
		fsp_he->SetRange(x1,x2);

		Spline *sp_he_ave = new Spline("sp_he_ave", nnodes_ams, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_he_ave = sp_he_ave->GetTF1Pointer();  
		TF1 *fit_he_ave = (TF1*) file_f3->Get(Form("fit_he_ave_BR%d", 2426+iBR_true)); 

		HistTools::CopyParameters(fit_he_ave, fsp_he_ave);
		// double x1, x2;
		fit_he_ave->GetRange(x1,x2);
		fsp_he_ave->SetRange(x1,x2);

		TF1 *fsp_he2 = HistTools::CombineTF1(fsp_he, fsp_he_ave, HistTools::Divide, "fsp_he2"); // He(R,t) Fit vs. <He(R)> Fit
		fsp_he2->SetRange(0.7, 60); 

		TH1 *h_ratio0 = (TH1*) file_f3->Get(Form("h_ratio0_BR%d", 2426+iBR_true)); 
		TH1 *h_ratio1 = (TH1*) file_f3->Get(Form("h_ratio1_BR%d", 2426+iBR_true)); 

		TH1 *h_fitres[2]; 
		h_fitres[0] = (TH1*) file_f3->Get(Form("h_fitres_ace_BR%d", 2426+iBR_true)); 
		h_fitres[1] = (TH1*) file_f3->Get(Form("h_fitres_ams_BR%d", 2426+iBR_true)); 

		c1->cd(3); 
		gPad->SetGrid(); 
		gPad->SetLogx();
		gPad->SetMargin(0.12, 0.04, 0.1, 0.08);	

		////if (get_temp_int(element)<4) h_a->GetYaxis()->SetRangeUser(0.5, 2.2); 
		h_a->GetYaxis()->SetRangeUser(0, 3.5); 
		h_a->GetXaxis()->SetLimits(0.7, 60); 
	
		TLegend *l2 = new TLegend(0.3,0.7,0.8,0.9); 
		l2->AddEntry(h_ratio1, Form("ACE %s(R,t)/<%s(R)>", element, element), "PL");
		l2->AddEntry(h_ratio0, "AMS He(R,t)/<He(R)>", "PL");  
		l2->AddEntry(fsp_he2, "He(R,t)_Fit/<He(R)>_Fit", "L");
		// l2->AddEntry(fit_ratio2, Form("1+%6.3fexp(-%6.3fln(R))", fit_ratio2->GetParameter(0), fit_ratio2->GetParameter(1)), "L"); 

		h_a->Draw("E1X0"); 
		h_a->SetTitle(Form("F3; ; BR-%d %s(R,t)/<%s(R)>", iBR_true+2426, element, element));
		h_a->SetXTitle(Unit::GetEnergyLabel("GV"));
		// h_a->SetTitleSize(0.5,"y"); 
		HistTools::SetStyle(h_ratio1, kBlack, kFullCircle, 0.8, 1, 1);
		HistTools::SetStyle(h_ratio0, kBlack, kFullSquare, 0.8, 1, 1);
		fsp_he2->Draw("SAME");
		h_ratio1->Draw("E1X0 SAME"); 
		h_ratio0->Draw("E1X0 SAME");
		l2->Draw("SAME"); 

		c1->cd(8);
		gPad->SetGrid();
		gPad->SetLogx();  
		gPad->SetMargin(0.12, 0.04, 0.2, 0.01);

		////if (get_temp_int(element)<4) h_b->GetYaxis()->SetRangeUser(-0.5, 0.5)  ; 
		h_b->GetYaxis()->SetRangeUser(-0.5, 0.5)  ; 
		h_b->GetXaxis()->SetLimits(0.7, 60); 

		TLegend *l_chi2 = new TLegend(0.62,0.8,0.9,0.9); 
		l_chi2->AddEntry(h_ratio1, Form("He_Fit(R,t)/<He_Fit(R)> chi2=%0.4f", g_chi2_ace->GetY()[iBR]), "PL");
		l_chi2->AddEntry(h_ratio0, Form("He_Fit(R,t)/<He_Fit(R)> chi2=%0.4f", g_chi2_ams->GetY()[iBR]), "PL");  

		h_b->Draw("E1X0"); 
		h_b->SetTitle("; ; Data/Fit-1"); 
		h_b->SetXTitle(Unit::GetEnergyLabel("GV"));

		HistTools::SetStyle(h_fitres[0], kBlack, kFullCircle, 0.8, 1, 1);
		HistTools::SetStyle(h_fitres[1], kBlack, kFullSquare, 0.8, 1, 1);

		h_fitres[0]->Draw("E1X0 SAME");
		h_fitres[1]->Draw("E1X0 SAME");
		l_chi2->Draw("SAME"); 
	   }
	
	   //F4
	   if (!strcmp(element, "C") || !strcmp(element, "N") || get_index(element)>3 ){

		TGraphErrors *g_chi2_ace = (TGraphErrors *) file_f4->Get("g_chi2_ace"); 
		TGraphErrors *g_chi2_ams = (TGraphErrors *) file_f4->Get("g_chi2_ams");

		TH1 *h_a = new TH1D("", "", 3000, 0, 3000); 
		TH1 *h_b = new TH1D("", "", 3000, 0, 3000); 

		TF1 *fit_ratio = (TF1*) file_f4->Get(Form("fit_ratio_BR%d", 2426+iBR_true)); 
		TH1 *h_ratio0 = (TH1*) file_f4->Get(Form("h_ratio0_BR%d", 2426+iBR_true)); 
		TH1 *h_ratio1 = (TH1*) file_f4->Get(Form("h_ratio1_BR%d", 2426+iBR_true)); 

		TH1 *h_fitres[2]; 
		h_fitres[0] = (TH1*) file_f4->Get(Form("h_fitres_ace_BR%d", 2426+iBR_true)); 
		h_fitres[1] = (TH1*) file_f4->Get(Form("h_fitres_ams_BR%d", 2426+iBR_true));

		c1->cd(4); 
		gPad->SetGrid(); 
		gPad->SetLogx();
		gPad->SetMargin(0.12, 0.04, 0.1, 0.08);	

		////if (get_temp_int(element)<4) h_a->GetYaxis()->SetRangeUser(0.5, 2.2); 
		h_a->GetYaxis()->SetRangeUser(0, 3.5); 
		h_a->GetXaxis()->SetLimits(0.7, 60); 
	
		TLegend *l2 = new TLegend(0.3,0.7,0.8,0.9); 
		l2->AddEntry(h_ratio1, Form("ACE %s(R,t)/<%s(R)>", element, element), "PL");
		l2->AddEntry(h_ratio0, "AMS He(R,t)/<He(R)>", "PL");  
		l2->AddEntry(fit_ratio, Form("1+%0.2fpm%0.2fexp(-%0.2fpm%0.2fR)", fit_ratio->GetParameter(0), fit_ratio->GetParError(0), fit_ratio->GetParameter(1), fit_ratio->GetParError(1)), "L"); 
		// l2->AddEntry(fit_ratio2, Form("1+%6.3fexp(-%6.3fln(R))", fit_ratio2->GetParameter(0), fit_ratio2->GetParameter(1)), "L"); 

		h_a->Draw("E1X0"); 
		h_a->SetTitle(Form("F4; ; BR-%d %s(R,t)/<%s(R)>", iBR_true+2426, element, element));
		h_a->SetXTitle(Unit::GetEnergyLabel("GV"));
		// h_a->SetTitleSize(0.5,"y"); 
		HistTools::SetStyle(h_ratio1, kBlack, kFullCircle, 0.8, 1, 1);
		HistTools::SetStyle(h_ratio0, kBlack, kFullSquare, 0.8, 1, 1);
		fit_ratio->Draw("SAME");
		h_ratio1->Draw("E1X0 SAME"); 
		h_ratio0->Draw("E1X0 SAME");
		l2->Draw("SAME"); 

		c1->cd(9);
		gPad->SetGrid();
		gPad->SetLogx();  
		gPad->SetMargin(0.12, 0.04, 0.2, 0.01);

		////if (get_temp_int(element)<4) h_b->GetYaxis()->SetRangeUser(-0.5, 0.5)  ; 
		h_b->GetYaxis()->SetRangeUser(-0.5, 0.5)  ; 
		h_b->GetXaxis()->SetLimits(0.7, 60); 

		TLegend *l_chi2 = new TLegend(0.62,0.8,0.9,0.9); 
		l_chi2->AddEntry(h_ratio1, Form("chi2=%0.4f", g_chi2_ace->GetY()[iBR]), "PL");
		l_chi2->AddEntry(h_ratio0, Form("chi2=%0.4f", g_chi2_ams->GetY()[iBR]), "PL");  

		h_b->Draw("E1X0"); 
		h_b->SetTitle("; ; Data/Fit-1"); 
		h_b->SetXTitle(Unit::GetEnergyLabel("GV"));

		HistTools::SetStyle(h_fitres[0], kBlack, kFullCircle, 0.8, 1, 1);
		HistTools::SetStyle(h_fitres[1], kBlack, kFullSquare, 0.8, 1, 1);

		h_fitres[0]->Draw("E1X0 SAME");
		h_fitres[1]->Draw("E1X0 SAME");
		l_chi2->Draw("SAME"); 

	   }

	   //F5
	   for( int k=0; k<1; ++k ){

	   	TGraphErrors *g_chi2_ace = (TGraphErrors *) file_f5->Get("g_chi2_ace"); 
		TGraphErrors *g_chi2_ams = (TGraphErrors *) file_f5->Get("g_chi2_ams");  

		TH1 *h_a = new TH1D("", "", 3000, 0, 3000); 
		TH1 *h_b = new TH1D("", "", 3000, 0, 3000); 

		TF1 *fit_ratio = (TF1*) file_f5->Get(Form("fit_ratio_BR%d", 2426+iBR_true)); 
		TH1 *h_ratio0 = (TH1*) file_f5->Get(Form("h_ratio0_BR%d", 2426+iBR_true)); 
		TH1 *h_ratio1 = (TH1*) file_f5->Get(Form("h_ratio1_BR%d", 2426+iBR_true)); 

		TH1 *h_fitres[2]; 
		h_fitres[0] = (TH1*) file_f5->Get(Form("h_fitres_ace_BR%d", 2426+iBR_true)); 
		h_fitres[1] = (TH1*) file_f5->Get(Form("h_fitres_ams_BR%d", 2426+iBR_true));

		c1->cd(5); 
		gPad->SetGrid(); 
		gPad->SetLogx();
		gPad->SetMargin(0.12, 0.04, 0.1, 0.08);	

		h_a->GetYaxis()->SetRangeUser(0, 3.5);  
		h_a->GetXaxis()->SetLimits(0.7, 60); 
	
		TLegend *l2 = new TLegend(0.3,0.7,0.8,0.9); 
		l2->AddEntry(h_ratio1, Form("ACE %s(R,t)/<%s(R)>", element, element), "PL");
		l2->AddEntry(h_ratio0, "AMS He(R,t)/<He(R)>", "PL");  
		l2->AddEntry(fit_ratio, Form("1+%0.2fpm%0.3fexp(-%0.2fpm%0.2fR)", fit_ratio->GetParameter(0), fit_ratio->GetParError(0), fit_ratio->GetParameter(1), fit_ratio->GetParError(1)), "L"); 
		// l2->AddEntry(fit_ratio2, Form("1+%6.3fexp(-%6.3fln(R))", fit_ratio2->GetParameter(0), fit_ratio2->GetParameter(1)), "L"); 

		h_a->Draw("E1X0"); 
		h_a->SetTitle(Form("F5; ; BR-%d %s(R,t)/<%s(R)>", iBR_true+2426, element, element));
		h_a->SetXTitle(Unit::GetEnergyLabel("GV"));
		// h_a->SetTitleSize(0.5,"y"); 
		HistTools::SetStyle(h_ratio1, kBlack, kFullCircle, 0.8, 1, 1);
		HistTools::SetStyle(h_ratio0, kBlack, kFullSquare, 0.8, 1, 1);
		fit_ratio->Draw("SAME");
		h_ratio1->Draw("E1X0 SAME"); 
		h_ratio0->Draw("E1X0 SAME");
		l2->Draw("SAME"); 

		c1->cd(10);
		gPad->SetGrid();
		gPad->SetLogx();  
		gPad->SetMargin(0.12, 0.04, 0.2, 0.01);

		h_b->GetYaxis()->SetRangeUser(-0.5, 0.5)  ; 
		h_b->GetXaxis()->SetLimits(0.7, 60); 

		TLegend *l_chi2 = new TLegend(0.62,0.8,0.9,0.9); 
		l_chi2->AddEntry(h_ratio1, Form("chi2=%0.4f", g_chi2_ace->GetY()[iBR]), "PL");
		l_chi2->AddEntry(h_ratio0, Form("chi2=%0.4f", g_chi2_ams->GetY()[iBR]), "PL");  

		h_b->Draw("E1X0"); 
		h_b->SetTitle("; ; Data/Fit-1"); 
		h_b->SetXTitle(Unit::GetEnergyLabel("GV"));

		HistTools::SetStyle(h_fitres[0], kBlack, kFullCircle, 0.8, 1, 1);
		HistTools::SetStyle(h_fitres[1], kBlack, kFullSquare, 0.8, 1, 1);

		h_fitres[0]->Draw("E1X0 SAME");
		h_fitres[1]->Draw("E1X0 SAME");
		l_chi2->Draw("SAME"); 
	
	   }

	   // break; 

	   if (iBR==0) c1->Print(Form("./data/ACE/fill/compare_fake_%s_p1.pdf(", element), "pdf"); 
	   if (iBR>0 && nBRs-1) c1->Print(Form("./data/ACE/fill/compare_fake_%s_p1.pdf", element), "pdf");  
	   if (iBR==nBRs-1) c1->Print(Form("./data/ACE/fill/compare_fake_%s_p1.pdf)", element), "pdf"); 

	   if (iBR+2426==2472-1) iBR_true += 3; 
		else iBR_true ++; 
	} 

	//P2F123-ACE
	for (int bin=0; bin<7; ++bin){

		//F1
		for (int k=0; k<1; ++k){ 

			TGraphErrors *g_residual_ace = (TGraphErrors *) file_f1->Get(Form("g_residual_ace_%d", bin)); 
			TGraphErrors *g_resierr_ace = (TGraphErrors *) file_f1->Get(Form("g_resierr_ace_%d", bin)); 
			TGraphErrors *g_ratio_ace = (TGraphErrors *) file_f1->Get(Form("g_ratio_ace_%d", bin)); 
			TGraphErrors *g_ratio_ace_fit = (TGraphErrors *) file_f1->Get(Form("g_ratio_ace_fit_%d", bin)); 

			c2->cd(1);
			gPad->SetGrid();
			gPad->SetMargin(0.12, 0.04, 0.1, 0.08);

			TLegend *l_a = new TLegend(0.62, 0.8, 0.9, 0.9); 
			l_a->AddEntry(g_ratio_ace, "Averaged ACE Data/Template", "PL"); 
			l_a->AddEntry(g_ratio_ace_fit, "Fit", "PL"); 

			g_ratio_ace->GetXaxis()->SetTimeDisplay(1);
			g_ratio_ace->GetXaxis()->SetTimeFormat("%m-%y");
			g_ratio_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g_ratio_ace->GetXaxis()->SetLimits(UBRToTime(2426-6), UBRToTime(2426+nBRs+6)); 
			g_ratio_ace->GetYaxis()->SetRangeUser(-0.3, 3); 
		
			HistTools::SetStyle(g_ratio_ace, kPink, kFullCircle, 0.8, 1, 1);
			g_ratio_ace->SetTitle(Form("F1; ; ACE %s Data/Template", element)); 
			g_ratio_ace->Draw("AP"); 
			HistTools::SetStyle(g_ratio_ace_fit, kGreen-3, kFullSquare, 0.6, 1, 1);
			g_ratio_ace_fit->SetLineColor(kGreen); 
			g_ratio_ace_fit->SetFillStyle(1001); 
			g_ratio_ace_fit->Draw("3 SAME"); 

			TGraph *g_ratio_ace_fiterr = (TGraph *) g_ratio_ace_fit->Clone();
			g_ratio_ace_fiterr->Draw("LX SAME"); 
			
			l_a->Draw("SAME"); 

			c2->cd(6);
			gPad->SetGrid(); 
			gPad->SetMargin(0.12, 0.04, 0.2, 0.01);

			TLegend *l_b = new TLegend(0.42, 0.8, 0.9, 0.9); 
			l_b->AddEntry(g_residual_ace, Form("ACE He_Fit Residuals (%0.4f GV)", h_ace->GetBinCenter(bin+1)), "PL"); 

			g_residual_ace->SetTitle(Form(" ; ; ACE %s Fitting Residuals", element));  
			g_residual_ace->GetXaxis()->SetTimeDisplay(1);
			g_residual_ace->GetXaxis()->SetTimeFormat("%m-%y");
			g_residual_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			//g_residual_ace->GetXaxis()->SetTitleSize(0.7);
			g_residual_ace->GetYaxis()->SetRangeUser(-1, 2); 
			g_residual_ace->GetXaxis()->SetLimits(UBRToTime(2426-6), UBRToTime(2426+nBRs+6)); 
			HistTools::SetStyle(g_residual_ace, kPink, kFullCircle, 0.8, 1, 1); 
			g_residual_ace->Draw("APL");
			l_b->Draw("SAME"); 		

		} 

		//F2
		for (int k=0; k<1; ++k){ 

			TGraphErrors *g_residual_ace = (TGraphErrors *) file_f2->Get(Form("g_residual_ace_%d", bin)); 
			TGraphErrors *g_resierr_ace = (TGraphErrors *) file_f2->Get(Form("g_resierr_ace_%d", bin)); 
			TGraphErrors *g_ratio_ace = (TGraphErrors *) file_f2->Get(Form("g_ratio_ace_%d", bin)); 
			TGraphErrors *g_ratio_ace_fit = (TGraphErrors *) file_f2->Get(Form("g_ratio_ace_fit_%d", bin)); 

			c2->cd(2);
			gPad->SetGrid();
			gPad->SetMargin(0.12, 0.04, 0.1, 0.08);

			TLegend *l_a = new TLegend(0.62, 0.8, 0.9, 0.9); 
			l_a->AddEntry(g_ratio_ace, "ACE Data/Template", "PL"); 
			l_a->AddEntry(g_ratio_ace_fit, "Fit", "PL"); 

			g_ratio_ace->SetTitle(Form("F2; ; ACE %s Data/Template", element));  
			g_ratio_ace->GetXaxis()->SetTimeDisplay(1);
			g_ratio_ace->GetXaxis()->SetTimeFormat("%m-%y");
			g_ratio_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00");  
			g_ratio_ace->GetYaxis()->SetRangeUser(-0.2, 2.2); 
		
			HistTools::SetStyle(g_ratio_ace, kPink, kFullCircle, 0.8, 1, 1);
			g_ratio_ace->Draw("AP"); 
			//g_ratio_ace_fit->Draw("PSAME"); 
			HistTools::SetStyle(g_ratio_ace_fit, kGreen-3, kFullSquare, 0.6, 1, 1);
			g_ratio_ace_fit->SetLineColor(kGreen); 
			g_ratio_ace_fit->SetFillStyle(1001); 
			g_ratio_ace_fit->Draw("3 SAME"); 
			l_a->Draw("SAME"); 

			//g_ratio_ace->Print("range"); 
			//g_ratio_ace_fit->Print("range"); 

			TGraph *g_ratio_ace_fiterr = (TGraph *) g_ratio_ace_fit->Clone();
			g_ratio_ace_fiterr->Draw("LX SAME");

			c2->cd(7);
			gPad->SetGrid(); 
			gPad->SetMargin(0.12, 0.04, 0.2, 0.01);

			TLegend *l_b = new TLegend(0.42, 0.8, 0.9, 0.9); 
			l_b->AddEntry(g_residual_ace, Form("ACE He_Fit Residuals (%0.4f GV)", h_ace->GetBinCenter(bin+1)), "PL"); 

			g_residual_ace->SetTitle(Form(" ; ; ACE %s Fitting Residuals", element));  
			g_residual_ace->GetXaxis()->SetTimeDisplay(1);
			g_residual_ace->GetXaxis()->SetTimeFormat("%m-%y");
			g_residual_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			//g_residual_ace->GetXaxis()->SetTitleSize(0.7);
			g_residual_ace->GetYaxis()->SetRangeUser(-1, 2);
			HistTools::SetStyle(g_residual_ace, kPink, kFullCircle, 0.8, 1, 1); 
			g_residual_ace->Draw("APL");
			l_b->Draw("SAME"); 

			c2->cd(12); 
			gPad->SetGrid(); 
			gPad->SetMargin(0.12, 0.04, 0.2, 0.01);

			TLegend *l_c = new TLegend(0.42, 0.8, 0.9, 0.9); 
			l_c->AddEntry(g_resierr_ace, Form("ACE He_Fit Residuals (%0.4f GV)", h_ace->GetBinCenter(bin+1)), "PL"); 

			g_resierr_ace->SetTitle(Form(" ; ; ACE %s (Data-Fit)/FitError", element));  
			g_resierr_ace->GetXaxis()->SetTimeDisplay(1);
			g_resierr_ace->GetXaxis()->SetTimeFormat("%m-%y");
			g_resierr_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			//g_resierr_ace->GetXaxis()->SetTitleSize(0.7);
			//g_resierr_ace->GetYaxis()->SetRangeUser(-1, 2); 
			g_resierr_ace->GetXaxis()->SetLimits(UBRToTime(2426-6), UBRToTime(2426+nBRs+6)); 
			HistTools::SetStyle(g_resierr_ace, kPink, kFullCircle, 0.8, 1, 1); 
			g_resierr_ace->Draw("APL");

			//g_resierr_ace->Print("range"); 

			l_c->Draw("SAME"); 		
	
		} 

		//F3
		for (int k=0; k<1; ++k){ 

			TGraphErrors *g_residual_ace = (TGraphErrors *) file_f3->Get(Form("g_residual_ace_%d", bin)); 
			TGraphErrors *g_resierr_ace = (TGraphErrors *) file_f3->Get(Form("g_resierr_ace_%d", bin)); 
			TGraphErrors *g_ratio_ace = (TGraphErrors *) file_f3->Get(Form("g_ratio_ace_%d", bin));	
			TGraphErrors *g_ratio_ace_fit = (TGraphErrors *) file_f3->Get(Form("g_ratio_ace_fit_%d", bin));  

			c2->cd(3);
			gPad->SetGrid();
			gPad->SetMargin(0.12, 0.04, 0.1, 0.08);

			TLegend *l_a = new TLegend(0.62, 0.8, 0.9, 0.9); 
			l_a->AddEntry(g_ratio_ace, "ACE Data/Template", "PL"); 
			l_a->AddEntry(g_ratio_ace_fit, "Fit", "PL"); 

			g_ratio_ace->SetTitle(Form("F3; ; ACE %s Data/Template", element));  
			g_ratio_ace->GetXaxis()->SetTimeDisplay(1);
			g_ratio_ace->GetXaxis()->SetTimeFormat("%m-%y");
			g_ratio_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			////if (get_temp_int(element)<4) g_ratio_ace->GetYaxis()->SetLimits(0.2, 2.2); 
			g_ratio_ace->GetYaxis()->SetRangeUser(-0.2, 2.2); 
		
			HistTools::SetStyle(g_ratio_ace, kPink, kFullCircle, 0.8, 1, 1);
			g_ratio_ace->Draw("AP"); 
			HistTools::SetStyle(g_ratio_ace_fit, kGreen-3, kFullSquare, 0.6, 1, 1);
			g_ratio_ace_fit->SetLineColor(kGreen); 
			g_ratio_ace_fit->SetFillStyle(1001); 
			g_ratio_ace_fit->Draw("3 SAME"); 
			l_a->Draw("SAME"); 

			TGraph *g_ratio_ace_fiterr = (TGraph *) g_ratio_ace_fit->Clone();
			g_ratio_ace_fiterr->Draw("LX SAME");

			c2->cd(8);
			gPad->SetGrid(); 
			gPad->SetMargin(0.12, 0.04, 0.2, 0.01);

			TLegend *l_b = new TLegend(0.42, 0.8, 0.9, 0.9); 
			l_b->AddEntry(g_residual_ace, Form("ACE He_Fit Residuals (%0.4f GV)", h_ace->GetBinCenter(bin+1)), "PL"); 

			g_residual_ace->SetTitle(Form(" ; ; ACE %s Fitting Residuals", element));  
			g_residual_ace->GetXaxis()->SetTimeDisplay(1);
			g_residual_ace->GetXaxis()->SetTimeFormat("%m-%y");
			g_residual_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			//g_residual_ace->GetXaxis()->SetTitleSize(0.7);
			////if (get_temp_int(element)<4) g_residual_ace->GetYaxis()->SetRangeUser(-0.5, 0.5) ;
			g_residual_ace->GetYaxis()->SetRangeUser(-1, 2); 
			HistTools::SetStyle(g_residual_ace, kPink, kFullCircle, 0.8, 1, 1); 
			g_residual_ace->Draw("APL");
			l_b->Draw("SAME"); 

			c2->cd(13); 
			gPad->SetGrid(); 
			gPad->SetMargin(0.12, 0.04, 0.2, 0.01);

			TLegend *l_c = new TLegend(0.42, 0.8, 0.9, 0.9); 
			l_c->AddEntry(g_resierr_ace, Form("ACE He_Fit Residuals (%0.4f GV)", h_ace->GetBinCenter(bin+1)), "PL"); 

			g_resierr_ace->SetTitle(Form(" ; ; ACE %s (Data-Fit)/FitError", element));  
			g_resierr_ace->GetXaxis()->SetTimeDisplay(1);
			g_resierr_ace->GetXaxis()->SetTimeFormat("%m-%y");
			g_resierr_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			//g_resierr_ace->GetXaxis()->SetTitleSize(0.7);
			//g_resierr_ace->GetYaxis()->SetRangeUser(-1, 2); 
			g_resierr_ace->GetXaxis()->SetLimits(UBRToTime(2426-6), UBRToTime(2426+nBRs+6)); 
			HistTools::SetStyle(g_resierr_ace, kPink, kFullCircle, 0.8, 1, 1); 
			g_resierr_ace->Draw("APL");

			g_resierr_ace->Print("range"); 

			l_c->Draw("SAME"); 	
	
		} 

		//F4
	 	if (!strcmp(element, "C") || !strcmp(element, "N") || get_index(element)>3){

			TGraphErrors *g_residual_ace = (TGraphErrors *) file_f4->Get(Form("g_residual_ace_%d", bin));
			TGraphErrors *g_resierr_ace = (TGraphErrors *) file_f4->Get(Form("g_resierr_ace_%d", bin));  
			TGraphErrors *g_ratio_ace = (TGraphErrors *) file_f4->Get(Form("g_ratio_ace_%d", bin)); 
			TGraphErrors *g_ratio_ace_fit = (TGraphErrors *) file_f4->Get(Form("g_ratio_ace_fit_%d", bin));

			c2->cd(4);
			gPad->SetGrid();
			gPad->SetMargin(0.12, 0.04, 0.1, 0.08);

			TLegend *l_a = new TLegend(0.62, 0.8, 0.9, 0.9); 
			l_a->AddEntry(g_ratio_ace, "ACE Data/Template", "PL"); 
			l_a->AddEntry(g_ratio_ace_fit, "Fit", "PL"); 

			g_ratio_ace->SetTitle(Form("F4; ; ACE %s Data/Template", element));  
			g_ratio_ace->GetXaxis()->SetTimeDisplay(1);
			g_ratio_ace->GetXaxis()->SetTimeFormat("%m-%y");
			g_ratio_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			////if (get_temp_int(element)<4) g_ratio_ace->GetYaxis()->SetLimits(0.2, 2.2); 
			g_ratio_ace->GetYaxis()->SetRangeUser(-0.2, 2.2); 
		
			HistTools::SetStyle(g_ratio_ace, kPink, kFullCircle, 0.8, 1, 1); 
			g_ratio_ace->Draw("AP"); 
			HistTools::SetStyle(g_ratio_ace_fit, kGreen-3, kFullSquare, 0.6, 1, 1);
			g_ratio_ace_fit->SetLineColor(kGreen); 
			g_ratio_ace_fit->SetFillStyle(1001); 
			g_ratio_ace_fit->Draw("3 SAME"); 
			l_a->Draw("SAME"); 

			TGraph *g_ratio_ace_fiterr = (TGraph *) g_ratio_ace_fit->Clone();
			g_ratio_ace_fiterr->Draw("LX SAME");

			c2->cd(9);
			gPad->SetGrid(); 
			gPad->SetMargin(0.12, 0.04, 0.2, 0.01);

			TLegend *l_b = new TLegend(0.42, 0.8, 0.9, 0.9); 
			l_b->AddEntry(g_residual_ace, Form("ACE He_Fit Residuals (%0.4f GV)", h_ace->GetBinCenter(bin+1)), "PL"); 

			g_residual_ace->SetTitle(Form(" ; ; ACE %s Fitting Residuals", element));  
			g_residual_ace->GetXaxis()->SetTimeDisplay(1);
			g_residual_ace->GetXaxis()->SetTimeFormat("%m-%y");
			g_residual_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			//g_residual_ace->GetXaxis()->SetTitleSize(0.7);
			////if (get_temp_int(element)<4) g_residual_ace->GetYaxis()->SetRangeUser(-0.5, 0.5) ;
			g_residual_ace->GetYaxis()->SetRangeUser(-1, 2);
			HistTools::SetStyle(g_residual_ace, kPink, kFullCircle, 0.8, 1, 1); 
			g_residual_ace->Draw("APL");
			l_b->Draw("SAME"); 

			c2->cd(14); 
			gPad->SetGrid(); 
			gPad->SetMargin(0.12, 0.04, 0.2, 0.01);

			TLegend *l_c = new TLegend(0.42, 0.8, 0.9, 0.9); 
			l_c->AddEntry(g_resierr_ace, Form("ACE He_Fit Residuals (%0.4f GV)", h_ace->GetBinCenter(bin+1)), "PL"); 

			g_resierr_ace->SetTitle(Form(" ; ; ACE %s (Data-Fit)/FitError", element));  
			g_resierr_ace->GetXaxis()->SetTimeDisplay(1);
			g_resierr_ace->GetXaxis()->SetTimeFormat("%m-%y");
			g_resierr_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			//g_resierr_ace->GetXaxis()->SetTitleSize(0.7);
			//g_resierr_ace->GetYaxis()->SetRangeUser(-1, 2); 
			g_resierr_ace->GetXaxis()->SetLimits(UBRToTime(2426-6), UBRToTime(2426+nBRs+6)); 
			HistTools::SetStyle(g_resierr_ace, kPink, kFullCircle, 0.8, 1, 1); 
			g_resierr_ace->Draw("APL");

			//g_resierr_ace->Print("range"); 

			l_c->Draw("SAME"); 
		}

		//F5
		for (int k=0; k<1; ++k){ 

			TGraphErrors *g_residual_ace = (TGraphErrors *) file_f5->Get(Form("g_residual_ace_%d", bin)); 
			TGraphErrors *g_resierr_ace = (TGraphErrors *) file_f5->Get(Form("g_resierr_ace_%d", bin)); 
			TGraphErrors *g_ratio_ace = (TGraphErrors *) file_f5->Get(Form("g_ratio_ace_%d", bin)); 
			TGraphErrors *g_ratio_ace_fit = (TGraphErrors *) file_f5->Get(Form("g_ratio_ace_fit_%d", bin));

			c2->cd(5);
			gPad->SetGrid();
			gPad->SetMargin(0.12, 0.04, 0.1, 0.08);

			TLegend *l_a = new TLegend(0.62, 0.8, 0.9, 0.9); 
			l_a->AddEntry(g_ratio_ace, "ACE Data/Template", "PL"); 
			l_a->AddEntry(g_ratio_ace_fit, "Fit", "PL"); 

			g_ratio_ace->SetTitle(Form("F5; ; ACE %s Data/Template", element));  
			g_ratio_ace->GetXaxis()->SetTimeDisplay(1);
			g_ratio_ace->GetXaxis()->SetTimeFormat("%m-%y");
			g_ratio_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			////if (get_temp_int(element)<4) g_ratio_ace->GetYaxis()->SetLimits(0.2, 2.2); 
			g_ratio_ace->GetYaxis()->SetRangeUser(-0.2, 2.2); 
		
			HistTools::SetStyle(g_ratio_ace, kPink, kFullCircle, 0.8, 1, 1); 
			g_ratio_ace->Draw("AP"); 
			HistTools::SetStyle(g_ratio_ace_fit, kGreen-3, kFullSquare, 0.6, 1, 1);
			g_ratio_ace_fit->SetLineColor(kGreen); 
			g_ratio_ace_fit->SetFillStyle(1001); 
			g_ratio_ace_fit->Draw("3 SAME"); 
			l_a->Draw("SAME"); 

			TGraph *g_ratio_ace_fiterr = (TGraph *) g_ratio_ace_fit->Clone();
			g_ratio_ace_fiterr->Draw("LX SAME");

			c2->cd(10);
			gPad->SetGrid(); 
			gPad->SetMargin(0.12, 0.04, 0.2, 0.01);

			TLegend *l_b = new TLegend(0.42, 0.8, 0.9, 0.9); 
			l_b->AddEntry(g_residual_ace, Form("ACE He_Fit Residuals (%0.4f GV)", h_ace->GetBinCenter(bin+1)), "PL"); 

			g_residual_ace->SetTitle(Form(" ; ; ACE %s Fitting Residuals", element));  
			g_residual_ace->GetXaxis()->SetTimeDisplay(1);
			g_residual_ace->GetXaxis()->SetTimeFormat("%m-%y");
			g_residual_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			//g_residual_ace->GetXaxis()->SetTitleSize(0.7);
			////if (get_temp_int(element)<4) g_residual_ace->GetYaxis()->SetRangeUser(-0.5, 0.5) ;
			g_residual_ace->GetYaxis()->SetRangeUser(-1, 2) ; 
			HistTools::SetStyle(g_residual_ace, kPink, kFullCircle, 0.8, 1, 1); 
			g_residual_ace->Draw("APL");
			l_b->Draw("SAME"); 

			c2->cd(15); 
			gPad->SetGrid(); 
			gPad->SetMargin(0.12, 0.04, 0.2, 0.01);

			TLegend *l_c = new TLegend(0.42, 0.8, 0.9, 0.9); 
			l_c->AddEntry(g_resierr_ace, Form("ACE He_Fit Residuals (%0.4f GV)", h_ace->GetBinCenter(bin+1)), "PL"); 

			g_resierr_ace->SetTitle(Form(" ; ; ACE %s (Data-Fit)/FitError", element));  
			g_resierr_ace->GetXaxis()->SetTimeDisplay(1);
			g_resierr_ace->GetXaxis()->SetTimeFormat("%m-%y");
			g_resierr_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			//g_resierr_ace->GetXaxis()->SetTitleSize(0.7);
			//g_resierr_ace->GetYaxis()->SetRangeUser(-1, 2); 
			g_resierr_ace->GetXaxis()->SetLimits(UBRToTime(2426-6), UBRToTime(2426+nBRs+6)); 
			HistTools::SetStyle(g_resierr_ace, kPink, kFullCircle, 0.8, 1, 1); 
			g_resierr_ace->Draw("APL");
	
			//g_resierr_ace->Print("range"); 
	
			l_c->Draw("SAME"); 
		}

		if (bin==0) c2->Print(Form("./data/ACE/fill/compare_fake_%s_p2_ace.pdf(", element), "pdf"); 
		if (bin>0 && bin<7-1) c2->Print(Form("./data/ACE/fill/compare_fake_%s_p2_ace.pdf", element), "pdf");  
		if (bin==7-1){
			c2->Print(Form("./data/ACE/fill/compare_fake_%s_p2_ace.pdf)", element), "pdf"); 
		} 

	}

	//P2F123-ACE (Data-Fit)/Error_Data
	for (int bin=0; bin<7; ++bin){

		//F1
		for (int k=0; k<1; ++k){ 

			TGraphErrors *g_residual_ace = (TGraphErrors *) file_f1->Get(Form("g_residual_ace_%d", bin)); 
			TGraphErrors *g_resierr_ace = (TGraphErrors *) file_f1->Get(Form("g_resierr_ace_%d", bin)); 
			TGraphErrors *g_ratio_ace = (TGraphErrors *) file_f1->Get(Form("g_ratio_ace_%d", bin)); 
			TGraphErrors *g_ratio_ace_fit = (TGraphErrors *) file_f1->Get(Form("g_ratio_ace_fit_%d", bin));
	
			TH1 *h_resierr = new TH1D("","", 40, -5, 7);  
			for (int i=0; i<g_resierr_ace->GetN(); ++i){
				h_resierr->Fill(g_resierr_ace->GetY()[i]); 
			} 
			for (int i=0; i<g_resierr_ace->GetN(); ++i){
				g_resierr_ace->SetPointError(i, 0, h_resierr->GetRMS(2)); 
			} 

			c2->cd(1);
			gPad->SetGrid(); 
			gPad->SetMargin(0.12, 0.04, 0.1, 0.08);

			TLegend *l_c = new TLegend(0.42, 0.8, 0.9, 0.9); 
			l_c->AddEntry(g_resierr_ace, Form("%0.4f GV", h_ace->GetBinCenter(bin+1)), "PL"); 

			g_resierr_ace->SetTitle(Form(" ; ; ACE %s (Data-Fit)/FitError", element));  
			g_resierr_ace->GetXaxis()->SetTimeDisplay(1);
			g_resierr_ace->GetXaxis()->SetTimeFormat("%m-%y");
			g_resierr_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			//g_resierr_ace->GetXaxis()->SetTitleSize(0.7);
			//g_resierr_ace->GetYaxis()->SetRangeUser(-1, 2); 
			g_resierr_ace->GetXaxis()->SetLimits(UBRToTime(2426-6), UBRToTime(2426+nBRs+6)); 
			HistTools::SetStyle(g_resierr_ace, kPink, kFullCircle, 0.8, 1, 1); 
			g_resierr_ace->Draw("APL");

			//g_resierr_ace->Print("range"); 

			l_c->Draw("SAME");  
	
			c2->cd(6);
			gPad->SetGrid(); 
			gPad->SetMargin(0.12, 0.04, 0.2, 0.01);

			TLegend *l_d = new TLegend(0.42, 0.8, 0.9, 0.9); 
			l_d->AddEntry(h_resierr, Form("%0.4f GV", h_ace->GetBinCenter(bin+1)), "PL"); 

			h_resierr->SetTitle(Form(" ; ACE %s (Data-Fit)/FitError; counts", element)); 
			h_resierr->SetFillColor(kBlue-2);
			h_resierr->Draw("HIST"); 
			l_d->Draw("SAME"); 

						
		} 

		//F2
		for (int k=0; k<1; ++k){ 

			TGraphErrors *g_residual_ace = (TGraphErrors *) file_f2->Get(Form("g_residual_ace_%d", bin)); 
			TGraphErrors *g_resierr_ace = (TGraphErrors *) file_f2->Get(Form("g_resierr_ace_%d", bin)); 
			TGraphErrors *g_ratio_ace = (TGraphErrors *) file_f2->Get(Form("g_ratio_ace_%d", bin)); 
			TGraphErrors *g_ratio_ace_fit = (TGraphErrors *) file_f2->Get(Form("g_ratio_ace_fit_%d", bin)); 

			TH1 *h_resierr = new TH1D("","", 40, -5, 7);  
			for (int i=0; i<g_resierr_ace->GetN(); ++i){
				h_resierr->Fill(g_resierr_ace->GetY()[i]); 
			} 
			for (int i=0; i<g_resierr_ace->GetN(); ++i){
				g_resierr_ace->SetPointError(i, 0, h_resierr->GetRMS(2)); 
			} 

			c2->cd(2); 
			gPad->SetGrid(); 
			gPad->SetMargin(0.12, 0.04, 0.1, 0.08);

			TLegend *l_c = new TLegend(0.42, 0.8, 0.9, 0.9); 
			l_c->AddEntry(g_resierr_ace, Form("%0.4f GV", h_ace->GetBinCenter(bin+1)), "PL"); 

			g_resierr_ace->SetTitle(Form(" ; ; ACE %s (Data-Fit)/FitError", element));  
			g_resierr_ace->GetXaxis()->SetTimeDisplay(1);
			g_resierr_ace->GetXaxis()->SetTimeFormat("%m-%y");
			g_resierr_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			//g_resierr_ace->GetXaxis()->SetTitleSize(0.7);
			//g_resierr_ace->GetYaxis()->SetRangeUser(-1, 2); 
			g_resierr_ace->GetXaxis()->SetLimits(UBRToTime(2426-6), UBRToTime(2426+nBRs+6)); 
			HistTools::SetStyle(g_resierr_ace, kPink, kFullCircle, 0.8, 1, 1); 
			g_resierr_ace->Draw("APL");

			//g_resierr_ace->Print("range"); 

			l_c->Draw("SAME"); 

			c2->cd(7);
			gPad->SetGrid(); 
			gPad->SetMargin(0.12, 0.04, 0.2, 0.01);

			TLegend *l_d = new TLegend(0.42, 0.8, 0.9, 0.9); 
			l_d->AddEntry(h_resierr, Form("%0.4f", h_ace->GetBinCenter(bin+1)), "PL"); 

			h_resierr->SetTitle(Form(" ; ACE %s (Data-Fit)/FitError; counts", element)); 
			h_resierr->SetFillColor(kBlue-2);
			h_resierr->Draw("HIST"); 
			l_d->Draw("SAME"); 		
	
		} 

		//F3
		for (int k=0; k<1; ++k){ 

			TGraphErrors *g_residual_ace = (TGraphErrors *) file_f3->Get(Form("g_residual_ace_%d", bin)); 
			TGraphErrors *g_resierr_ace = (TGraphErrors *) file_f3->Get(Form("g_resierr_ace_%d", bin)); 
			TGraphErrors *g_ratio_ace = (TGraphErrors *) file_f3->Get(Form("g_ratio_ace_%d", bin));	
			TGraphErrors *g_ratio_ace_fit = (TGraphErrors *) file_f3->Get(Form("g_ratio_ace_fit_%d", bin));  

			TH1 *h_resierr = new TH1D("","", 40, -5, 7);  
			for (int i=0; i<g_resierr_ace->GetN(); ++i){
				h_resierr->Fill(g_resierr_ace->GetY()[i]); 
			} 
			for (int i=0; i<g_resierr_ace->GetN(); ++i){
				g_resierr_ace->SetPointError(i, 0, h_resierr->GetRMS(2)); 
			} 

			c2->cd(3); 
			gPad->SetGrid(); 
			gPad->SetMargin(0.12, 0.04, 0.1, 0.08);

			TLegend *l_c = new TLegend(0.42, 0.8, 0.9, 0.9); 
			l_c->AddEntry(g_resierr_ace, Form("%0.4f GV", h_ace->GetBinCenter(bin+1)), "PL"); 

			g_resierr_ace->SetTitle(Form(" ; ; ACE %s (Data-Fit)/FitError", element));  
			g_resierr_ace->GetXaxis()->SetTimeDisplay(1);
			g_resierr_ace->GetXaxis()->SetTimeFormat("%m-%y");
			g_resierr_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			//g_resierr_ace->GetXaxis()->SetTitleSize(0.7);
			//g_resierr_ace->GetYaxis()->SetRangeUser(-1, 2); 
			g_resierr_ace->GetXaxis()->SetLimits(UBRToTime(2426-6), UBRToTime(2426+nBRs+6)); 
			HistTools::SetStyle(g_resierr_ace, kPink, kFullCircle, 0.8, 1, 1); 
			g_resierr_ace->Draw("APL");

			//g_resierr_ace->Print("range"); 

			l_c->Draw("SAME"); 
	
			c2->cd(8);
			gPad->SetGrid(); 
			gPad->SetMargin(0.12, 0.04, 0.2, 0.01);

			TLegend *l_d = new TLegend(0.42, 0.8, 0.9, 0.9); 
			l_d->AddEntry(h_resierr, Form("%0.4f GV", h_ace->GetBinCenter(bin+1)), "PL"); 

			h_resierr->SetTitle(Form(" ; ACE %s (Data-Fit)/FitError; counts", element)); 
			h_resierr->SetFillColor(kBlue-2);
			h_resierr->Draw("HIST"); 
			l_d->Draw("SAME"); 
	
		} 

		//F4
	 	if (!strcmp(element, "C") || !strcmp(element, "N") || get_index(element)>3){

			TGraphErrors *g_residual_ace = (TGraphErrors *) file_f4->Get(Form("g_residual_ace_%d", bin));
			TGraphErrors *g_resierr_ace = (TGraphErrors *) file_f4->Get(Form("g_resierr_ace_%d", bin));  
			TGraphErrors *g_ratio_ace = (TGraphErrors *) file_f4->Get(Form("g_ratio_ace_%d", bin)); 
			TGraphErrors *g_ratio_ace_fit = (TGraphErrors *) file_f4->Get(Form("g_ratio_ace_fit_%d", bin));

			TH1 *h_resierr = new TH1D("","", 40, -5, 7);  
			for (int i=0; i<g_resierr_ace->GetN(); ++i){
				h_resierr->Fill(g_resierr_ace->GetY()[i]); 
			} 
			for (int i=0; i<g_resierr_ace->GetN(); ++i){
				g_resierr_ace->SetPointError(i, 0, h_resierr->GetRMS(2)); 
			} 

			c2->cd(4); 
			gPad->SetGrid(); 
			gPad->SetMargin(0.12, 0.04, 0.1, 0.08);

			TLegend *l_c = new TLegend(0.42, 0.8, 0.9, 0.9); 
			l_c->AddEntry(g_resierr_ace, Form("%0.4f GV", h_ace->GetBinCenter(bin+1)), "PL"); 

			g_resierr_ace->SetTitle(Form(" ; ; ACE %s (Data-Fit)/FitError", element));  
			g_resierr_ace->GetXaxis()->SetTimeDisplay(1);
			g_resierr_ace->GetXaxis()->SetTimeFormat("%m-%y");
			g_resierr_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			//g_resierr_ace->GetXaxis()->SetTitleSize(0.7);
			//g_resierr_ace->GetYaxis()->SetRangeUser(-1, 2); 
			g_resierr_ace->GetXaxis()->SetLimits(UBRToTime(2426-6), UBRToTime(2426+nBRs+6)); 
			HistTools::SetStyle(g_resierr_ace, kPink, kFullCircle, 0.8, 1, 1); 
			g_resierr_ace->Draw("APL");

			//g_resierr_ace->Print("range"); 

			l_c->Draw("SAME"); 

			c2->cd(9);
			gPad->SetGrid(); 
			gPad->SetMargin(0.12, 0.04, 0.2, 0.01);

			TLegend *l_d = new TLegend(0.42, 0.8, 0.9, 0.9); 
			l_d->AddEntry(h_resierr, Form("%0.4f GV", h_ace->GetBinCenter(bin+1)), "PL"); 

			h_resierr->SetTitle(Form(" ; ACE %s (Data-Fit)/FitError; counts", element)); 
			h_resierr->SetFillColor(kBlue-2);
			h_resierr->Draw("HIST"); 
			l_d->Draw("SAME"); 
		}

		//F5
		for (int k=0; k<1; ++k){ 

			TGraphErrors *g_residual_ace = (TGraphErrors *) file_f5->Get(Form("g_residual_ace_%d", bin)); 
			TGraphErrors *g_resierr_ace = (TGraphErrors *) file_f5->Get(Form("g_resierr_ace_%d", bin)); 
			TGraphErrors *g_ratio_ace = (TGraphErrors *) file_f5->Get(Form("g_ratio_ace_%d", bin)); 
			TGraphErrors *g_ratio_ace_fit = (TGraphErrors *) file_f5->Get(Form("g_ratio_ace_fit_%d", bin));

			TH1 *h_resierr = new TH1D("","", 40, -5, 7);  
			for (int i=0; i<g_resierr_ace->GetN(); ++i){
				h_resierr->Fill(g_resierr_ace->GetY()[i]); 
			} 
			for (int i=0; i<g_resierr_ace->GetN(); ++i){
				g_resierr_ace->SetPointError(i, 0, h_resierr->GetRMS(2)); 
			} 

			c2->cd(5); 
			gPad->SetGrid(); 
			gPad->SetMargin(0.12, 0.04, 0.1, 0.08); 

			TLegend *l_c = new TLegend(0.42, 0.8, 0.9, 0.9); 
			l_c->AddEntry(g_resierr_ace, Form("%0.4f GV", h_ace->GetBinCenter(bin+1)), "PL"); 

			g_resierr_ace->SetTitle(Form(" ; ; ACE %s (Data-Fit)/FitError", element));  
			g_resierr_ace->GetXaxis()->SetTimeDisplay(1);
			g_resierr_ace->GetXaxis()->SetTimeFormat("%m-%y");
			g_resierr_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			//g_resierr_ace->GetXaxis()->SetTitleSize(0.7);
			//g_resierr_ace->GetYaxis()->SetRangeUser(-1, 2); 
			g_resierr_ace->GetXaxis()->SetLimits(UBRToTime(2426-6), UBRToTime(2426+nBRs+6)); 
			HistTools::SetStyle(g_resierr_ace, kPink, kFullCircle, 0.8, 1, 1); 
			g_resierr_ace->Draw("APL");
	
			//g_resierr_ace->Print("range"); 
	
			l_c->Draw("SAME"); 

			c2->cd(10);
			gPad->SetGrid(); 
			gPad->SetMargin(0.12, 0.04, 0.2, 0.01);

			TLegend *l_d = new TLegend(0.42, 0.8, 0.9, 0.9); 
			l_d->AddEntry(h_resierr, Form("%0.4f GV", h_ace->GetBinCenter(bin+1)), "PL"); 

			h_resierr->SetTitle(Form(" ; ACE %s (Data-Fit)/FitError; counts", element)); 
			h_resierr->SetFillColor(kBlue-2);
			h_resierr->Draw("HIST"); 
			l_d->Draw("SAME"); 
		}

		if (bin==0) c2->Print(Form("./data/ACE/fill/compare_fake_%s_p2_df.pdf(", element), "pdf"); 
		if (bin>0 && bin<7-1) c2->Print(Form("./data/ACE/fill/compare_fake_%s_p2_df.pdf", element), "pdf");  
		if (bin==7-1){
			c2->Print(Form("./data/ACE/fill/compare_fake_%s_p2_df.pdf)", element), "pdf"); 
		} 

	}

	//P2F23-AMS
	for (int bin=0; bin<h_BR_he[0]->GetNbinsX(); ++bin){

/*
		//F1
		for (int k=0; k<1; ++k){ 

			TGraphErrors *g_residual_ams = (TGraphErrors *) file_f1->Get(Form("g_residual_ams_%d", bin)); 
			//TGraphErrors *g_ratio_ams = (TGraphErrors *) file_f1->Get(Form("g_ratio_ams_%d", bin)); 

			c3->cd(1);
			gPad->SetGrid();

			TLegend *l_a = new TLegend(0.62, 0.8, 0.9, 0.9); 
			l_a->AddEntry(g_ratio_ams, "Averaged AMS Data/Template", "PL"); 

			g_ratio_ams->GetXaxis()->SetTimeDisplay(1);
			g_ratio_ams->GetXaxis()->SetTimeFormat("%m-%y");
			g_ratio_ams->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g_ratio_ams->GetXaxis()->SetLimits(UBRToTime(2426-6), UBRToTime(2426+nBRs+6)); 
			//g_ratio_ams->GetYaxis()->SetRangeUser(0.2, 2.2); 
		
			HistTools::SetStyle(g_ratio_ams, kBlue, kFullCircle, 1.1, 1, 1);
			g_ratio_ams->SetTitle(Form("F1; ; AMS %s Data/Template-1", element)); 
			g_ratio_ams->Draw("APL"); 
			l_a->Draw("SAME"); 

			c3->cd(4);
			gPad->SetGrid(); 

			TLegend *l_b = new TLegend(0.62, 0.8, 0.9, 0.9); 
			l_b->AddEntry(g_residual_ams, Form("AMS He_Fit Residuals (%0.4f GV)", h_ams->GetBinCenter(bin+1)), "PL"); 

			g_residual_ams->SetTitle(Form(" ; ; AMS %s Fitting Residuals", element));  
			g_residual_ams->GetXaxis()->SetTimeDisplay(1);
			g_residual_ams->GetXaxis()->SetTimeFormat("%m-%y");
			g_residual_ams->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			//g_residual_ams->GetXaxis()->SetTitleSize(0.7);
			g_residual_ams->GetYaxis()->SetRangeUser(-0.5, 0.5) ;
			g_residual_ams->GetXaxis()->SetLimits(UBRToTime(2426-6), UBRToTime(2426+nBRs+6)); 
			HistTools::SetStyle(g_residual_ams, kPink, kFullCircle, 1.1, 1, 1); 
			g_residual_ams->Draw("APL");
			l_b->Draw("SAME"); 		

		} 
*/

		//F2
		for (int k=0; k<1; ++k){ 

			TGraphErrors *g_residual_ams = (TGraphErrors *) file_f2->Get(Form("g_residual_ams_%d", bin)); 
			TGraphErrors *g_ratio_ams = (TGraphErrors *) file_f2->Get(Form("g_ratio_ams_%d", bin)); 
			TGraphErrors *g_ratio_ams_fit = (TGraphErrors *) file_f2->Get(Form("g_ratio_ams_fit_%d", bin)); 

			c3->cd(2);
			gPad->SetGrid();
			gPad->SetMargin(0.13, 0.04, 0.1, 0.08);

			TLegend *l_a = new TLegend(0.62, 0.8, 0.9, 0.9); 
			l_a->AddEntry(g_ratio_ams, "AMS Data/Template", "PL"); 
			l_a->AddEntry(g_ratio_ams_fit, "Fit", "PL"); 

			g_ratio_ams->SetTitle(Form("F2; ; AMS %s Data/Template", element));  
			g_ratio_ams->GetXaxis()->SetTimeDisplay(1);
			g_ratio_ams->GetXaxis()->SetTimeFormat("%m-%y");
			g_ratio_ams->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00");
			g_ratio_ams->GetYaxis()->SetRangeUser(0.7, 1.3);  
		
			HistTools::SetStyle(g_ratio_ams, kBlue, kFullCircle, 0.8, 1, 1); 
			g_ratio_ams->Draw("AP"); 
			HistTools::SetStyle(g_ratio_ams_fit, kGreen-3, kFullSquare, 0.6, 1, 1);
			g_ratio_ams_fit->SetLineColor(kGreen); 
			g_ratio_ams_fit->SetFillStyle(1001); 
			g_ratio_ams_fit->Draw("3 SAME"); 
			// g_ratio_ams_fit->Print("range"); 
			l_a->Draw("SAME"); 

			TGraph *g_ratio_ams_fiterr = (TGraph *) g_ratio_ams_fit->Clone();
			g_ratio_ams_fiterr->Draw("LX SAME");

			c3->cd(7);
			gPad->SetGrid(); 
			gPad->SetMargin(0.13, 0.04, 0.2, 0.01);

			TLegend *l_b = new TLegend(0.42, 0.8, 0.9, 0.9); 
			l_b->AddEntry(g_residual_ams, Form("AMS He_Fit Residuals (%0.3f GV)", h_ams->GetBinCenter(bin+1)), "PL"); 

			g_residual_ams->SetTitle(Form(" ; ; AMS %s Fitting Residuals", element));  
			g_residual_ams->GetXaxis()->SetTimeDisplay(1);
			g_residual_ams->GetXaxis()->SetTimeFormat("%m-%y");
			g_residual_ams->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			//g_residual_ams->GetXaxis()->SetTitleSize(0.7);
			g_residual_ams->GetYaxis()->SetRangeUser(-0.5, 0.5)  ;
			HistTools::SetStyle(g_residual_ams, kBlue, kFullCircle, 0.8, 1, 1); 
			g_residual_ams->Draw("APL"); 
			l_b->Draw("SAME"); 		
	
		} 

		//F3
		for (int k=0; k<1; ++k){ 

			TGraphErrors *g_residual_ams = (TGraphErrors *) file_f3->Get(Form("g_residual_ams_%d", bin)); 
			TGraphErrors *g_ratio_ams = (TGraphErrors *) file_f3->Get(Form("g_ratio_ams_%d", bin)); 
			TGraphErrors *g_ratio_ams_fit = (TGraphErrors *) file_f3->Get(Form("g_ratio_ams_fit_%d", bin)); 

			c3->cd(3);
			gPad->SetGrid();
			gPad->SetMargin(0.13, 0.04, 0.1, 0.08);

			TLegend *l_a = new TLegend(0.62, 0.8, 0.9, 0.9); 
			l_a->AddEntry(g_ratio_ams, "AMS Data/Template", "PL"); 
			l_a->AddEntry(g_ratio_ams_fit, "Fit", "PL"); 

			g_ratio_ams->SetTitle(Form("F3; ; AMS %s Data/Template", element));
			g_ratio_ams->GetXaxis()->SetTimeDisplay(1);
			g_ratio_ams->GetXaxis()->SetTimeFormat("%m-%y");
			g_ratio_ams->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g_ratio_ams->GetYaxis()->SetRangeUser(0.7, 1.3);  
		
			HistTools::SetStyle(g_ratio_ams, kBlue, kFullCircle, 0.8, 1, 1);
			g_ratio_ams->Draw("AP"); 
			HistTools::SetStyle(g_ratio_ams_fit, kGreen-3, kFullSquare, 0.6, 1, 1);
			g_ratio_ams_fit->SetLineColor(kGreen); 
			g_ratio_ams_fit->SetFillStyle(1001); 
			g_ratio_ams_fit->Draw("3 SAME");
			// g_ratio_ams_fit->Print("range"); 
			l_a->Draw("SAME"); 

			TGraph *g_ratio_ams_fiterr = (TGraph *) g_ratio_ams_fit->Clone();
			g_ratio_ams_fiterr->Draw("LX SAME");

			c3->cd(8);
			gPad->SetGrid(); 
			gPad->SetMargin(0.13, 0.04, 0.2, 0.01);

			TLegend *l_b = new TLegend(0.42, 0.8, 0.9, 0.9); 
			l_b->AddEntry(g_residual_ams, Form("AMS He_Fit Residuals (%0.3f GV)", h_ams->GetBinCenter(bin+1)), "PL"); 

			g_residual_ams->SetTitle(Form(" ; ; AMS %s Fitting Residuals", element));  
			g_residual_ams->GetXaxis()->SetTimeDisplay(1);
			g_residual_ams->GetXaxis()->SetTimeFormat("%m-%y");
			g_residual_ams->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			//g_residual_ams->GetXaxis()->SetTitleSize(0.7);
			g_residual_ams->GetYaxis()->SetRangeUser(-0.5, 0.5)  ;
			HistTools::SetStyle(g_residual_ams, kBlue, kFullCircle, 0.8, 1, 1); 
			g_residual_ams->Draw("APL");
			l_b->Draw("SAME"); 	
	
		} 

		//F4
	 	if (!strcmp(element, "C") || !strcmp(element, "N") || get_index(element)>3){

			TGraphErrors *g_residual_ams = (TGraphErrors *) file_f4->Get(Form("g_residual_ams_%d", bin)); 
			TGraphErrors *g_ratio_ams = (TGraphErrors *) file_f4->Get(Form("g_ratio_ams_%d", bin)); 
			TGraphErrors *g_ratio_ams_fit = (TGraphErrors *) file_f4->Get(Form("g_ratio_ams_fit_%d", bin)); 

			c3->cd(4);
			gPad->SetGrid();
			gPad->SetMargin(0.13, 0.04, 0.1, 0.08);

			TLegend *l_a = new TLegend(0.62, 0.8, 0.9, 0.9); 
			l_a->AddEntry(g_ratio_ams, "AMS Data/Template", "PL"); 
			l_a->AddEntry(g_ratio_ams_fit, "Fit", "PL"); 

			g_ratio_ams->SetTitle(Form("F4; ; AMS %s Data/Template", element));
			g_ratio_ams->GetXaxis()->SetTimeDisplay(1);
			g_ratio_ams->GetXaxis()->SetTimeFormat("%m-%y");
			g_ratio_ams->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g_ratio_ams->GetYaxis()->SetRangeUser(0.7, 1.3);  
		
			HistTools::SetStyle(g_ratio_ams, kBlue, kFullCircle, 0.8, 1, 1);
			g_ratio_ams->Draw("AP"); 	
			HistTools::SetStyle(g_ratio_ams_fit, kGreen-3, kFullSquare, 0.6, 1, 1);
			g_ratio_ams_fit->SetLineColor(kGreen); 
			g_ratio_ams_fit->SetFillStyle(1001); 
			g_ratio_ams_fit->Draw("3 SAME");
			l_a->Draw("SAME"); 

			TGraph *g_ratio_ams_fiterr = (TGraph *) g_ratio_ams_fit->Clone();
			g_ratio_ams_fiterr->Draw("LX SAME");

			c3->cd(9);
			gPad->SetGrid(); 
			gPad->SetMargin(0.13, 0.04, 0.2, 0.01);

			TLegend *l_b = new TLegend(0.42, 0.8, 0.9, 0.9); 
			l_b->AddEntry(g_residual_ams, Form("AMS He_Fit Residuals (%0.3f GV)", h_ams->GetBinCenter(bin+1)), "PL"); 

			g_residual_ams->SetTitle(Form(" ; ; AMS %s Fitting Residuals", element));  
			g_residual_ams->GetXaxis()->SetTimeDisplay(1);
			g_residual_ams->GetXaxis()->SetTimeFormat("%m-%y");
			g_residual_ams->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			//g_residual_ams->GetXaxis()->SetTitleSize(0.7);
			g_residual_ams->GetYaxis()->SetRangeUser(-0.5, 0.5)  ;
			HistTools::SetStyle(g_residual_ams, kBlue, kFullCircle, 0.8, 1, 1); 
			g_residual_ams->Draw("APL");
			l_b->Draw("SAME"); 
		}

		//F5
		for (int k=0; k<1; ++k){ 

			TGraphErrors *g_residual_ams = (TGraphErrors *) file_f5->Get(Form("g_residual_ams_%d", bin)); 
			TGraphErrors *g_ratio_ams = (TGraphErrors *) file_f5->Get(Form("g_ratio_ams_%d", bin)); 
			TGraphErrors *g_ratio_ams_fit = (TGraphErrors *) file_f5->Get(Form("g_ratio_ams_fit_%d", bin)); 

			c3->cd(5);
			gPad->SetGrid();
			gPad->SetMargin(0.13, 0.04, 0.1, 0.08);

			TLegend *l_a = new TLegend(0.62, 0.8, 0.9, 0.9); 
			l_a->AddEntry(g_ratio_ams, "AMS Data/Template", "PL"); 
			l_a->AddEntry(g_ratio_ams_fit, "Fit", "PL"); 

			g_ratio_ams->SetTitle(Form("F5; ; AMS %s Data/Template", element));
			g_ratio_ams->GetXaxis()->SetTimeDisplay(1);
			g_ratio_ams->GetXaxis()->SetTimeFormat("%m-%y");
			g_ratio_ams->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g_ratio_ams->GetYaxis()->SetRangeUser(0.7, 1.3);  
		
			HistTools::SetStyle(g_ratio_ams, kBlue, kFullCircle, 0.8, 1, 1);
			g_ratio_ams->Draw("AP"); 	
			HistTools::SetStyle(g_ratio_ams_fit, kGreen-3, kFullSquare, 0.6, 1, 1);
			g_ratio_ams_fit->SetLineColor(kGreen); 
			g_ratio_ams_fit->SetFillStyle(1001); 
			g_ratio_ams_fit->Draw("3 SAME");
			l_a->Draw("SAME"); 

			TGraph *g_ratio_ams_fiterr = (TGraph *) g_ratio_ams_fit->Clone();
			g_ratio_ams_fiterr->Draw("LX SAME");

			c3->cd(10);
			gPad->SetGrid(); 
			gPad->SetMargin(0.13, 0.04, 0.2, 0.01);

			TLegend *l_b = new TLegend(0.42, 0.8, 0.9, 0.9); 
			l_b->AddEntry(g_residual_ams, Form("AMS He_Fit Residuals (%0.3f GV)", h_ams->GetBinCenter(bin+1)), "PL"); 

			g_residual_ams->SetTitle(Form(" ; ; AMS %s Fitting Residuals", element));  
			g_residual_ams->GetXaxis()->SetTimeDisplay(1);
			g_residual_ams->GetXaxis()->SetTimeFormat("%m-%y");
			g_residual_ams->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			//g_residual_ams->GetXaxis()->SetTitleSize(0.7);
			g_residual_ams->GetYaxis()->SetRangeUser(-0.5, 0.5)  ;
			HistTools::SetStyle(g_residual_ams, kBlue, kFullCircle, 0.8, 1, 1); 
			g_residual_ams->Draw("APL");
			l_b->Draw("SAME"); 
		}

		if (bin==0) c3->Print(Form("./data/ACE/fill/compare_fake_%s_p2_ams.pdf(", element), "pdf"); 
		if (bin>0 && bin<h_BR_he[0]->GetNbinsX()-1) c3->Print(Form("./data/ACE/fill/compare_fake_%s_p2_ams.pdf", element), "pdf");  
		if (bin==h_BR_he[0]->GetNbinsX()-1){
			c3->Print(Form("./data/ACE/fill/compare_fake_%s_p2_ams.pdf)", element), "pdf"); 
		} 
	}

	//P3F23
	for (int i=0; i<1; ++i){

		//F1
		//for (int k=0; k<1; ++k){ 

		//} 


		//F2
		for (int k=0; k<1; ++k){ 

			TGraphErrors *g_ave_residual_ace_rig = (TGraphErrors *) file_f2->Get("g_ave_residual_ace_rig"); 
			TGraphErrors *g_ave_residual_ams_rig = (TGraphErrors *) file_f2->Get("g_ave_residual_ams_rig"); 

			c4->cd(2); 
			gPad->SetGrid(); 
			gPad->SetLogx();
			gPad->SetMargin(0.14, 0.1, 0.1, 0.1);	 

			HistTools::SetStyle(g_ave_residual_ace_rig, kBlack, kFullCircle, 1.1, 1, 1); 
			HistTools::SetStyle(g_ave_residual_ams_rig, kBlack, kFullSquare, 1.1, 1, 1); 

			//if (get_index(element)<4) g_ave_residual_ace_rig->GetYaxis()->SetRangeUser(-0.5, 0.5); 
			g_ave_residual_ace_rig->GetXaxis()->SetLimits(0.7, 60);
			g_ave_residual_ace_rig->SetTitle(Form("F2; ; Time-averaged Fit Residuals of %s(R,t)/<%s(R)>", element, element)); 
			g_ave_residual_ace_rig->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV")); 

			TLegend *l_both4 = new TLegend(0.62,0.8,0.9,0.9); 
			l_both4->AddEntry(g_ave_residual_ace_rig, "ACE", "PL");
			l_both4->AddEntry(g_ave_residual_ams_rig, "AMS", "PL"); 

			g_ave_residual_ace_rig->Draw("APL");
			g_ave_residual_ams_rig->Draw("PLSAME"); 
			l_both4->Draw("SAME"); 
		
		} 

		//F3
		for (int k=0; k<1; ++k){ 

			TGraphErrors *g_ave_residual_ace_rig = (TGraphErrors *) file_f3->Get("g_ave_residual_ace_rig"); 
			TGraphErrors *g_ave_residual_ams_rig = (TGraphErrors *) file_f3->Get("g_ave_residual_ams_rig"); 

			c4->cd(3); 
			gPad->SetGrid(); 
			gPad->SetLogx(); 
			gPad->SetMargin(0.14, 0.1, 0.1, 0.1);	

			HistTools::SetStyle(g_ave_residual_ace_rig, kBlack, kFullCircle, 1.1, 1, 1); 
			HistTools::SetStyle(g_ave_residual_ams_rig, kBlack, kFullSquare, 1.1, 1, 1); 

			//if (get_index(element)<4) g_ave_residual_ace_rig->GetYaxis()->SetRangeUser(-0.5, 0.5); 
			g_ave_residual_ace_rig->GetXaxis()->SetLimits(0.7, 60);
			g_ave_residual_ace_rig->SetTitle(Form("F3; ; Time-averaged Fit Residuals of %s(R,t)/<%s(R)>", element, element)); 
			g_ave_residual_ace_rig->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV")); 

			TLegend *l_both4 = new TLegend(0.62,0.8,0.9,0.9); 
			l_both4->AddEntry(g_ave_residual_ace_rig, "ACE", "PL");
			l_both4->AddEntry(g_ave_residual_ams_rig, "AMS", "PL"); 

			g_ave_residual_ace_rig->Draw("APL");
			g_ave_residual_ams_rig->Draw("PLSAME"); 
			l_both4->Draw("SAME"); 	

		} 
	
		//F4 
	 	if (!strcmp(element, "C") || !strcmp(element, "N") || get_index(element)>3){

			TGraphErrors *g_ave_residual_ace_rig = (TGraphErrors *) file_f4->Get("g_ave_residual_ace_rig"); 
			TGraphErrors *g_ave_residual_ams_rig = (TGraphErrors *) file_f4->Get("g_ave_residual_ams_rig"); 

			c4->cd(4); 
			gPad->SetGrid(); 
			gPad->SetLogx();
			gPad->SetMargin(0.14, 0.1, 0.1, 0.1);	 

			HistTools::SetStyle(g_ave_residual_ace_rig, kBlack, kFullCircle, 1.1, 1, 1); 
			HistTools::SetStyle(g_ave_residual_ams_rig, kBlack, kFullSquare, 1.1, 1, 1); 

			//if (get_index(element)<4) g_ave_residual_ace_rig->GetYaxis()->SetRangeUser(-0.5, 0.5); 
			g_ave_residual_ace_rig->GetXaxis()->SetLimits(0.7, 60);
			g_ave_residual_ace_rig->SetTitle(Form("F4; ; Time-averaged Fit Residuals of %s(R,t)/<%s(R)>", element, element)); 
			g_ave_residual_ace_rig->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV")); 

			TLegend *l_both4 = new TLegend(0.62,0.8,0.9,0.9); 
			l_both4->AddEntry(g_ave_residual_ace_rig, "ACE", "PL");
			l_both4->AddEntry(g_ave_residual_ams_rig, "AMS", "PL"); 

			g_ave_residual_ace_rig->Draw("APL");
			g_ave_residual_ams_rig->Draw("PLSAME"); 
			l_both4->Draw("SAME"); 
		}

		//F5
		for (int k=0; k<1; ++k){ 

			TGraphErrors *g_ave_residual_ace_rig = (TGraphErrors *) file_f5->Get("g_ave_residual_ace_rig"); 
			TGraphErrors *g_ave_residual_ams_rig = (TGraphErrors *) file_f5->Get("g_ave_residual_ams_rig"); 

			c4->cd(5); 
			gPad->SetGrid(); 
			gPad->SetLogx();
			gPad->SetMargin(0.14, 0.1, 0.1, 0.1);	 

			HistTools::SetStyle(g_ave_residual_ace_rig, kBlack, kFullCircle, 1.1, 1, 1); 
			HistTools::SetStyle(g_ave_residual_ams_rig, kBlack, kFullSquare, 1.1, 1, 1); 

			//if (get_index(element)<4) g_ave_residual_ace_rig->GetYaxis()->SetRangeUser(-0.5, 0.5); 
			g_ave_residual_ace_rig->GetXaxis()->SetLimits(0.7, 60);
			g_ave_residual_ace_rig->SetTitle(Form("F5; ; Time-averaged Fit Residuals of %s(R,t)/<%s(R)>", element, element)); 
			g_ave_residual_ace_rig->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV")); 

			TLegend *l_both4 = new TLegend(0.62,0.8,0.9,0.9); 
			l_both4->AddEntry(g_ave_residual_ace_rig, "ACE", "PL");
			l_both4->AddEntry(g_ave_residual_ams_rig, "AMS", "PL"); 

			g_ave_residual_ace_rig->Draw("APL");
			g_ave_residual_ams_rig->Draw("PLSAME"); 
			l_both4->Draw("SAME"); 
		}

		c4->Print(Form("./data/ACE/fill/compare_fake_%s_p3.png", element)); 
		
	}

	//P4F23

	for (int i=0; i<1; ++i){

		// F2
		for (int k=0; k<1; ++k){
			
			TGraphErrors *g_chi2_ace = (TGraphErrors*) file_f2->Get("g_chi2_ace"); 
			TGraphErrors *g_chi2_ams = (TGraphErrors*) file_f2->Get("g_chi2_ams"); 

			c5->cd(2);
			gPad->SetGrid(); 
			HistTools::SetStyle(g_chi2_ace, kPink, kFullCircle, 1., 1, 1); 
			TLegend *l_c5_1 = new TLegend(0.62,0.8,0.9,0.9); 
			l_c5_1->AddEntry(g_chi2_ace, "1+[0]*exp(-[1]*R)", "PL"); 
			g_chi2_ace->SetTitle(Form("F2; ; ACE %s(R,t)/<%s(R)> Fitting Chi2/NDF", element, element)); 
			g_chi2_ace->GetXaxis()->SetTimeDisplay(1);
			g_chi2_ace->GetXaxis()->SetTimeFormat("%m-%y"); 
			g_chi2_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g_chi2_ace->GetXaxis()->SetTitleSize(0.7); 
			//if (get_index(element)<4) g_chi2_ace->GetYaxis()->SetRangeUser(0., 4.5);
			//else g_chi2_ace->GetYaxis()->SetRangeUser(-2, 6);
			g_chi2_ace->Draw("APL"); 
			l_c5_1->Draw("SAME");

			c5->cd(7); 
			gPad->SetGrid(); 
			HistTools::SetStyle(g_chi2_ams, kPink, kFullCircle, 1., 1, 1); 
			TLegend *l_c5_2 = new TLegend(0.62,0.8,0.9,0.9); 
			l_c5_2->AddEntry(g_chi2_ams, "1+[0]*exp(-[1]*R)", "PL");
			g_chi2_ams->SetTitle(Form("; ; AMS %s(R,t)/<%s(R)> Fitting Chi2/NDF", element, element)); 
			g_chi2_ams->GetXaxis()->SetTimeDisplay(1);
			g_chi2_ams->GetXaxis()->SetTimeFormat("%m-%y");
			g_chi2_ams->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g_chi2_ams->GetXaxis()->SetTitleSize(0.7);
			//if (get_index(element)<4) g_chi2_ams->GetYaxis()->SetRangeUser(0., 2.2); 
			g_chi2_ams->Draw("APL");
			l_c5_2->Draw("SAME"); 

		}

		// F3
		for (int k=0; k<1; ++k){

			TGraphErrors *g_chi2_ace = (TGraphErrors*) file_f3->Get("g_chi2_ace"); 
			TGraphErrors *g_chi2_ams = (TGraphErrors*) file_f3->Get("g_chi2_ams"); 

			c5->cd(3);
			gPad->SetGrid(); 
			HistTools::SetStyle(g_chi2_ace, kPink, kFullCircle, 1., 1, 1); 
			TLegend *l_c5_1 = new TLegend(0.62,0.8,0.9,0.9); 
			l_c5_1->AddEntry(g_chi2_ace, "He(R,t)_Fit/<He(R)>_Fit", "PL"); 
			g_chi2_ace->SetTitle(Form("F3; ; ACE %s(R,t)/<%s(R)> Fitting Chi2/NDF", element, element)); 
			g_chi2_ace->GetXaxis()->SetTimeDisplay(1);
			g_chi2_ace->GetXaxis()->SetTimeFormat("%m-%y"); 
			g_chi2_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g_chi2_ace->GetXaxis()->SetTitleSize(0.7); 
			//if (get_index(element)<4) g_chi2_ace->GetYaxis()->SetRangeUser(0., 50.); 
			g_chi2_ace->Draw("APL"); 
			l_c5_1->Draw("SAME");

			c5->cd(8); 
			gPad->SetGrid(); 
			HistTools::SetStyle(g_chi2_ams, kPink, kFullCircle, 1., 1, 1); 
			TLegend *l_c5_2 = new TLegend(0.62,0.8,0.9,0.9); 
			l_c5_2->AddEntry(g_chi2_ams, "He(R,t)_Fit/<He(R)>_Fit", "PL");
			g_chi2_ams->SetTitle(Form("; ; AMS %s(R,t)/<%s(R)> Fitting Chi2/NDF", element, element)); 
			g_chi2_ams->GetXaxis()->SetTimeDisplay(1);
			g_chi2_ams->GetXaxis()->SetTimeFormat("%m-%y");
			g_chi2_ams->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g_chi2_ams->GetXaxis()->SetTitleSize(0.7);
			//if (get_index(element)<4) g_chi2_ams->GetYaxis()->SetRangeUser(0., 2.2); 
			g_chi2_ams->Draw("APL"); 
			l_c5_2->Draw("SAME"); 
		}

		//F4
	 	if (!strcmp(element, "C") || !strcmp(element, "N") || get_index(element)>3){

			TGraphErrors *g_chi2_ace = (TGraphErrors*) file_f4->Get("g_chi2_ace"); 
			TGraphErrors *g_chi2_ams = (TGraphErrors*) file_f4->Get("g_chi2_ams"); 

			c5->cd(4);
			gPad->SetGrid(); 
			HistTools::SetStyle(g_chi2_ace, kPink, kFullCircle, 1., 1, 1); 
			TLegend *l_c5_1 = new TLegend(0.62,0.8,0.9,0.9); 
			l_c5_1->AddEntry(g_chi2_ace, "He(R,t)_Fit/<He(R)>_Fit", "PL"); 
			g_chi2_ace->SetTitle(Form("F4; ; ACE %s(R,t)/<%s(R)> Fitting Chi2/NDF", element, element)); 
			g_chi2_ace->GetXaxis()->SetTimeDisplay(1);
			g_chi2_ace->GetXaxis()->SetTimeFormat("%m-%y"); 
			g_chi2_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g_chi2_ace->GetXaxis()->SetTitleSize(0.7); 
			//if (get_index(element)<4) g_chi2_ace->GetYaxis()->SetRangeUser(0., 4.5); 
			g_chi2_ace->Draw("APL"); 
			l_c5_1->Draw("SAME");

			c5->cd(9); 
			gPad->SetGrid(); 
			HistTools::SetStyle(g_chi2_ams, kPink, kFullCircle, 1., 1, 1); 
			TLegend *l_c5_2 = new TLegend(0.62,0.8,0.9,0.9); 
			l_c5_2->AddEntry(g_chi2_ams, "He(R,t)_Fit/<He(R)>_Fit", "PL");
			g_chi2_ams->SetTitle(Form("; ; AMS %s(R,t)/<%s(R)> Fitting Chi2/NDF", element, element)); 
			g_chi2_ams->GetXaxis()->SetTimeDisplay(1);
			g_chi2_ams->GetXaxis()->SetTimeFormat("%m-%y");
			g_chi2_ams->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g_chi2_ams->GetXaxis()->SetTitleSize(0.7);
			//if (get_index(element)<4) g_chi2_ams->GetYaxis()->SetRangeUser(0., 2.2); 
			g_chi2_ams->Draw("APL"); 
			l_c5_2->Draw("SAME"); 			

		}

		//F5
		for (int k=0; k<1; ++k){

			TGraphErrors *g_chi2_ace = (TGraphErrors*) file_f5->Get("g_chi2_ace"); 
			TGraphErrors *g_chi2_ams = (TGraphErrors*) file_f5->Get("g_chi2_ams"); 

			c5->cd(5); 
			gPad->SetGrid(); 
			HistTools::SetStyle(g_chi2_ace, kPink, kFullCircle, 1., 1, 1); 
			TLegend *l_c5_1 = new TLegend(0.62,0.8,0.9,0.9); 
			l_c5_1->AddEntry(g_chi2_ace, "He(R,t)_Fit/<He(R)>_Fit", "PL"); 
			g_chi2_ace->SetTitle(Form("F4; ; ACE %s(R,t)/<%s(R)> Fitting Chi2/NDF", element, element)); 
			g_chi2_ace->GetXaxis()->SetTimeDisplay(1);
			g_chi2_ace->GetXaxis()->SetTimeFormat("%m-%y"); 
			g_chi2_ace->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g_chi2_ace->GetXaxis()->SetTitleSize(0.7); 
			//if (get_index(element)<4) g_chi2_ace->GetYaxis()->SetRangeUser(0., 4.5); 
			g_chi2_ace->Draw("APL"); 
			l_c5_1->Draw("SAME");

			c5->cd(10); 
			gPad->SetGrid(); 
			HistTools::SetStyle(g_chi2_ams, kPink, kFullCircle, 1., 1, 1); 
			TLegend *l_c5_2 = new TLegend(0.62,0.8,0.9,0.9); 
			l_c5_2->AddEntry(g_chi2_ams, "He(R,t)_Fit/<He(R)>_Fit", "PL");
			g_chi2_ams->SetTitle(Form("; ; AMS %s(R,t)/<%s(R)> Fitting Chi2/NDF", element, element)); 
			g_chi2_ams->GetXaxis()->SetTimeDisplay(1);
			g_chi2_ams->GetXaxis()->SetTimeFormat("%m-%y");
			g_chi2_ams->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g_chi2_ams->GetXaxis()->SetTitleSize(0.7);
			//if (get_index(element)<4) g_chi2_ams->GetYaxis()->SetRangeUser(0., 10.); 
			g_chi2_ams->Draw("APL"); 
			l_c5_2->Draw("SAME"); 			

		}

		c5->Print(Form("./data/ACE/fill/compare_fake_%s_p4.png", element)); 

	}

	return 0; 
}; 

// Plot All Element Averaged Flux over Energy Bins
void ace_all_average(){

   gSystem->mkdir("data/ACE/average", true);

   TCanvas *c1 = new TCanvas("c1","", 2400, 1200);
   c1->Divide(2, 1);
   gStyle->SetPalette(1);

   TLegend *legend1 = new TLegend(0.1,0.6,0.28,0.9); // left, down, right, top
   TLegend *legend2 = new TLegend(0.1,0.6,0.28,0.9); // left, down, right, top

   // create axis histogram
   TH1 *h1 = HistTools::CreateAxis("h1", "haxis1", 0., 500., 7, 1e-10, 0.2e-6, false); 
   TH1 *h2 = HistTools::CreateAxis("h2", "haxis2", 0., 2.5, 7, 0.001, 1., false);

   c1->cd(1);
   gPad->SetLogy();
   h1->Draw();
   h1->SetTitle("All Element Kinetic Energy Spectrum (Averaged by Months)");
   h1->SetXTitle(Unit::GetEnergyLabel("MeV/n"));
   h1->SetYTitle(Unit::GetDifferentialFluxLabel("MeV/n cm"));
   c1->cd(2);
   gPad->SetLogy();
   h2->Draw();
   h2->SetTitle("All Element Rigidity Spectrum (Averaged by Months)");
   h2->SetXTitle(Unit::GetEnergyLabel("GV"));
   h2->SetYTitle(Unit::GetDifferentialFluxLabel("GV m"));

   for (int i=0; i<n_ele; ++i){

	const UInt_t FirstACEBR = 2240;
	vector<UInt_t> BRs;
	// we stop at BR 2493, which ends on 2016/05/23, just 3 days before the end of the data taking period for AMS nuclei
	for (UInt_t br=2426; br<=2493; ++br) { 
		if (br != 2472 && br != 2473) BRs.push_back(br-FirstACEBR); 
	}
	TH1D *h_ene = HistTools::GraphToHist(get_ace_average_graph( ACE_Element[i] , &BRs[0], BRs.size() ), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
	TH1 *h_rig = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene, ACE_Isotope[i], "MeV/n cm", "GV m", "_rig");
				
	HistTools::SetMarkerStyle(h_ene, HistTools::GetColorPalette(i, n_ele), kFullCircle, 0.9);	
	c1->cd(1);
	gPad->SetGrid(); 
	gPad->SetMargin(0.12, 0.08, 0.08, 0.08);
	legend1->AddEntry(h_ene, Form("%s", ACE_Element[i]), "p");
	legend1->SetNColumns(2);
	h_ene->Draw("E1X0 SAME"); 
	legend1->Draw("SAME");	

	HistTools::SetMarkerStyle(h_rig, HistTools::GetColorPalette(i, n_ele), kFullCircle, 0.9);
	c1->cd(2);
	gPad->SetGrid();
	gPad->SetMargin(0.12, 0.08, 0.08, 0.08);
	legend2->AddEntry(h_rig, Form("%s", ACE_Element[i]), "p");
	legend2->SetNColumns(2);
	h_rig->Draw("E1X0 SAME"); 
	legend2->Draw("SAME");

	TFile file2(Form("data/ACE/average/%s_average.root", ACE_Element[i]), "RECREATE");

	h_ene->Write(Form("h_kin_%s_BR_average", ACE_Element[i]));	
	h_rig->Write(Form("h_rig_%s_BR_average", ACE_Element[i]));
	
	file2.Write();
	file2.Close();

   }

   c1->Print("./data/ACE/average/h_kin_all_average.png");

}

// Convert CRIS Data into AMS Structure
void ace_convert(const char *element, Particle::Type isotope){

	double *kin_bins = get_kin_bins(element);
        double *SpallCorr = get_spall_corr(element);
	double *SpallCorrUnc = get_spall_corr_unc(element);
	double *EMed = get_EMed(element);

	gSystem->mkdir("data/ACE/convert/fluxtime", true);
	//gSystem->mkdir("data/ACE/convert/fluxenergy", true);
	gSystem->mkdir("data/ACE/convert/fluxrigidity", true);	
	
	const int nBins = 14;
	
	TCanvas *c2 = new TCanvas("c2","", 800, 500);
	c2->Divide(2, 1);		

	TFile fin(Form("data/ACE/fill/%s_fill.root", element));
	TFile fout(Form("data/ACE/convert/%s_convert.root", element), "RECREATE");

	int utime_0, utime_1; // ams time range  

	time_t *tran_ams = Experiments::GetMeasurementTimeRange(Experiments::AMS02, 1, 0); 
		utime_0 = tran_ams[0];
		tran_ams = Experiments::GetMeasurementTimeRange(Experiments::AMS02, 1, 78);
		utime_1 = tran_ams[1]; 

	// plot flux (MeV/n cm) vs. kinetic energy within AMS range
	for (int k=2240 ; k<=2529; k++){

		// utime (ace), utime_0 (ams), utime_1 (ams)
		// ace range = 2240 ~ 2529
		// ams range = 2426 ~ 2506 
		if ( UBRToTime(k) <= utime_1 && UBRToTime(k) >= utime_0 ){

			TH1F *h = (TH1F*) fin.Get(Form("h_kin_%s_BR%d", element, k));			
		
			HistTools::SetMarkerStyle(h, HistTools::GetColorPalette(k-2426, 81), kFullCircle, 1.1);
			gPad->SetGrid(); 
			//gPad->SetLogx(); 
			//gPad->SetLogy();
			h->GetYaxis()->SetRangeUser(0, h->GetBinContent(1)*4.0);
			h->SetTitle(Form("%s All BR Energy Spectrum; Kinetic Energy (MeV/nuc); Flux (/(cm^2 sr s)(MeV)", element));
			c2->cd(1);
			gPad->SetMargin(0.12, 0.08, 0.08, 0.08);
			h->Draw("E1X0 SAME"); 
		}
			
	}

	// convert to rigidity 
	for (int k=2240 ; k<=2529; k++){

		// utime (ace), utime_0 (ams), utime_1 (ams)
		// ace range = 2240 ~ 2529
		// ams range = 2426 ~ 2506 
		if ( UBRToTime(k) <= utime_1 && UBRToTime(k) >= utime_0 ){

			TH1F *h = (TH1F*) fin.Get(Form("h_kin_%s_BR%d", element, k));			
			TH1 *h_rig = HistTools::TransformEnergyAndDifferentialFluxNew(h, isotope, "MeV/n cm", "GV m", Form("_rig_%s_BR%d", element, k)); // (TH1 *hist, Particle::Type particle, const Char_t *from_flux_unit, const Char_t *to_flux_unit, const Char_t *suffix) 

			for (int i=0; i<h_rig->GetNbinsX(); ++i) { 
				//printf("%s BR=%d [%02u] %.2f-%.2f y=%10.4e dy=%10.4e  \n", element, k, i, h_rig->GetBinLowEdge(i+1), h_rig->GetBinLowEdge(i+2), h_rig->GetBinContent(i+1), h_rig->GetBinError(i+1)); 
			}

			HistTools::SetMarkerStyle(h_rig, HistTools::GetColorPalette(k-2426, 81), kFullCircle, 1.1);
			gPad->SetGrid(); 
			//gPad->SetLogx(); 
			//gPad->SetLogy();
			h_rig->GetYaxis()->SetRangeUser(0, h_rig->GetBinContent(1)*4.0);
			h_rig->SetTitle(Form("%s All BR Energy Spectrum", element));
			h_rig->SetXTitle(Unit::GetEnergyLabel("GV"));
			h_rig->SetYTitle(Unit::GetDifferentialFluxLabel("GV m"));

			c2->cd(2);
			gPad->SetMargin(0.12, 0.08, 0.08, 0.08);
			h_rig->Draw("E1X0 SAME"); 

			h_rig->Write(Form("h_rig_%s_BR%d", element, k)); 
		}
			
	}

	c2->Print(Form("./data/ACE/convert/fluxrigidity/h_rig_%s_all.png", element));

	// plot flux_time 
	TCanvas *c3 = new TCanvas("c3","",800,600);
	TLegend *legend3 = new TLegend(0.1,0.7,0.28,0.9); // left, down, right, top 

	const UInt_t FirstACEBR = 2240;
	vector<UInt_t> BRs;
	// we stop at BR 2493, which ends on 2016/05/23, just 3 days before the end of the data taking period for AMS nuclei
	for (UInt_t br=2426; br<=2493; ++br) { 
		if (br != 2472 && br != 2473) BRs.push_back(br-FirstACEBR); 
	}

	TH1D *h_ene = HistTools::GraphToHist( get_ace_average_graph( element, &BRs[0], BRs.size()), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
	TH1 *h_ave = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene, isotope, "MeV/n cm", "GV m", "_rig");

	for (int i=0; i<nBins; ++i){
		
		if (i%2==0){	
					
			TGraphErrors *g = new TGraphErrors(79);			
			for (int k=2240; k<=2529; k++){	
				if ( k <= UTimeToBR(utime_1) && k >= UTimeToBR(utime_0) && k!=2472 && k!=2473 ){
					TH1 *h_rig = (TH1*) fout.Get(Form("h_rig_%s_BR%d", element, k));			
					g->SetPoint(k-2426, UBRToTime(k), h_rig->GetBinContent(i+1)); 
					g->SetPointError(k-2426, 0, h_rig->GetBinError(i+1)); 
					double x, y;	
					g->GetPoint(k-2426, x, y);
					//printf("%s BR=%d [%02u] %.2f-%.2f x=%.0f y=%10.4e \n", element, k, i, h_rig->GetBinLowEdge(i+1), h_rig->GetBinLowEdge(i+2), x, y); 			 
				} 
			}
				
			//g->Set(79); 
			g->RemovePoint(46); 
			g->RemovePoint(46);

			g->GetXaxis()->SetTimeDisplay(1);
			g->GetXaxis()->SetTimeFormat("%m-%y");
			g->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00"); 
			g->GetXaxis()->SetTitleSize(0.7);			
			g->SetMarkerStyle(kFullCircle);
			g->SetMarkerColor(HistTools::GetColorPalette(i, nBins)); 
			g->SetLineColor(HistTools::GetColorPalette(i, nBins)); 
			g->SetLineWidth(1);
			g->SetMarkerSize(0.7);

			double x, y;	
			g->GetPoint(0, x, y);

			gPad->SetGrid();  
			//gPad->SetLogy();
			g->GetYaxis()->SetRangeUser(y-y*0.99, y+y*2.5);
			g->GetXaxis()->SetRangeUser(UBRToTime(2425), UBRToTime(2507));
			g->SetTitle(Form("%s All Energy Bin Flux Time Series; ;%s", element, Unit::GetDifferentialFluxLabel("GV m"))); 
			//g->Print();
			c3->cd(1);
			
			if (i%2==0) legend3->AddEntry(g, Form("%0.4f GV", h_ave->GetBinLowEdge(i)), "l");

			if (i/2==0){
				g->Draw("ALP");
				legend3->Draw("SAME");
			} else {
				g->Draw("LPSAME");
				legend3->Draw("SAME");
			}
		}
	} 

	c3->Print(Form("./data/ACE/convert/fluxtime/h_fluxtime_%s_all.png", element));	

	// plot normalized flux_time 
	TCanvas *c4 = new TCanvas("c4","",800,600);
	TLegend *legend4 = new TLegend(0.1,0.7,0.28,0.9); // left, down, right, top 

	//TH1D *h_ene = HistTools::GraphToHist( get_ace_average_graph( element, &BRs[0], BRs.size()), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
	//TH1 *h_ave = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene, isotope, "MeV/n cm", "GV m", "_rig");

	for (int i=0; i<nBins/2; ++i){
		printf("i=%d flux_ave=%10.4e dflux_ave=%10.4e \n", i, h_ave->GetBinContent(i+1), h_ave->GetBinError(i+1));
	}

	for (int i=0; i<nBins; ++i){

		if (i%2==0){
					
			//TGraphAsymmErrors *g = get_ace_graph( element, i, BRs.size() );
			TGraphAsymmErrors *g_norm = new TGraphAsymmErrors(BRs.size());	

			for (int iBR=2426; iBR<2493; iBR++){

				if (iBR != 2472 && iBR != 2473){

					TH1 *h_rig = (TH1*) fout.Get(Form("h_rig_%s_BR%d", element, iBR));
					
					g_norm->SetPoint(iBR-2426, UBRToTime(iBR), h_rig->GetBinContent(i+1)/h_ave->GetBinContent(i+1));
					g_norm->SetPointError(iBR-2426, 0., 0., h_rig->GetBinError(i+1)/h_ave->GetBinContent(i+1), h_rig->GetBinError(i+1)/h_ave->GetBinContent(i+1));

				}
			}

			g_norm->RemovePoint(46); 
   			g_norm->RemovePoint(46);
			g_norm->Set(63);

			g_norm->GetXaxis()->SetTimeDisplay(1);
			g_norm->GetXaxis()->SetTimeFormat("%m-%y");
			g_norm->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00");
			g_norm->GetXaxis()->SetTitleSize(0.7);			
			g_norm->SetMarkerStyle(kFullCircle);
			g_norm->SetMarkerColor(HistTools::GetColorPalette(i, nBins));
			g_norm->SetLineColor(HistTools::GetColorPalette(i, nBins));
			g_norm->SetLineWidth(1);
			g_norm->SetMarkerSize(0.7);

			double x, y;	
			g_norm->GetPoint(0, x, y);

			gPad->SetGrid();  
			//gPad->SetLogy();
			g_norm->GetYaxis()->SetRangeUser(y-y*0.99, y+y*0.9);
			g_norm->GetXaxis()->SetRangeUser(UBRToTime(2425), UBRToTime(2494));
			g_norm->SetTitle(Form("%s All Energy Bin Flux Time Series (Normalized by Months); ;Normalized Flux", element)); 

			//g_norm->Print();
			c4->cd(1);
			
			if (i%2==0) legend4->AddEntry(g_norm, Form("%0.4f GV", h_ave->GetBinLowEdge(i)), "l");

			if (i==0){
				g_norm->Draw("ALP");
				legend4->Draw("SAME");
			} else {
				g_norm->Draw("LPSAME");
				legend4->Draw("SAME");
			}
		}
	} 

	c4->Print(Form("./data/ACE/convert/fluxtime/h_norm_fluxtime_%s_all.png", element));	
	
	fout.Write();
	fout.Close();

	fin.Close();
}

// fit both ACE & AMS data
void ace_fitboth(int nnodes){
   ofstream bestnode(Form("data/ACE/compare/%dnodes.txt",nnodes));
	
   Experiments::DataPath = "data";	
   int data_value[n_ele] = {0, 18, 23, 27, 31, 20, 43, 22}; // p, He, Li, Be, B, C, N, O

   const UInt_t FirstACEBR = 2240;
   vector<UInt_t> BRs;
   // we stop at BR 2493, which ends on 2016/05/23, just 3 days before the end of the data taking period for AMS nuclei
   for (UInt_t br=2426; br<=2493; ++br) { 
	   if (br != 2472 && br != 2473) BRs.push_back(br-FirstACEBR); 
   }

   TCanvas *c5 = new TCanvas("c5","compare nodes", 2400, 900);
   c5->Divide(4,3);

   for (int i=0; i<4; i++){

	int nnodes_og=6;
	const char *AMS_Element[n_ele] = { "b", "c", "n", "o", "f" }; // to solve the upper/lower case conflict
	TFile file1(Form("data/amsfit/fit_result_node%d.root", nnodes_og));
	
	Spline *sp_ams = new Spline("sp_ams", nnodes_og, Spline::LogLog | Spline::PowerLaw);
	TF1 *fsp_ams = sp_ams->GetTF1Pointer(); 
	TF1 *fit_ams = (TF1*) file1.Get(Form("fsp_%s", AMS_Element[i]));

	HistTools::CopyParameters(fit_ams, fsp_ams);
	double x1, x2;
	fit_ams->GetRange(x1,x2);
	fsp_ams->SetRange(x1,x2);

	TH1 *h_ams = Experiments::GetMeasurementHistogram(Experiments::AMS02, data_value[i+4], 0); // load AMS data for a given element
	TH1 *h_ene = HistTools::GraphToHist(get_ace_average_graph( ACE_Element[i] , &BRs[0], BRs.size() ), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
	TH1 *h_ace = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene, ACE_Isotope[i], "MeV/n cm", "GV m", "_rig"); // load averaged ACE data for the same element, converted in rigidity
	h_ace->SetTitle(Form("ACE %s Flux in same time span", ACE_Element[i])); 

	UShort_t namsbins = h_ams->GetNbinsX();
	UShort_t nacebins = h_ace->GetNbinsX();
	double R1 = h_ace->GetBinLowEdge(1);
	double R2 = h_ams->GetBinLowEdge(namsbins+1);

	// initializes the X and Y position of the spline nodes
	double *xnodes = HistTools::BuildLogBins(R1, R2, nnodes); // xnodes will be an array of nnodes+1 items
	double *ynodes = new double[nnodes+1];
	UShort_t inode;
	for (inode = 0; inode < nnodes+1; ++inode)
	{
   		if (xnodes[inode] > h_ace->GetBinLowEdge(nacebins+1)) break;
   		ynodes[inode] = h_ace->GetBinContent(h_ace->FindBin(xnodes[inode]));
	}
	for (; inode < nnodes+1; ++inode)
	{
	   double x = xnodes[inode] < h_ams->GetBinLowEdge(1) ? h_ams->GetBinLowEdge(1) : xnodes[inode];
           //xnodes[inode] = x;
	   ynodes[inode] = h_ams->GetBinContent(h_ams->FindBin(x));
	}

	// create spline
	Spline *sp = new Spline("sp", nnodes, Spline::LogLog | Spline::PowerLaw, xnodes, ynodes);
	sp->SetSpectralIndices(-0.2, -2.8, 0, 0, -1, -1); // set initial values and limits for the spectral index at first and last node
	TF1 *fit = sp->GetTF1Pointer();	
	
	sp->Print(); 
	
	// create array of data to be fitted together
	TObjArray data; 
	data.Add(h_ace);
	data.Add(h_ams);

	// fit AMS and ACE data at the same time
	vector<double> rigmin, rigmax, chi2norm(2);
	ROOT::Fit::Fitter fitter;
	//HistTools::PrintFunction(fit);
	FitTools::SetCommonFitterOptions(fitter);
	FitTools::FitCombinedData(data, fit, "I", rigmin, rigmax, chi2norm, fitter, 3); 

	//ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1.E-4); // adjust tolerance
	//ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-2); 

	// compare original AMS fit with ACE/AMS Combined fit 
	TF1 *f_ratio = HistTools::CombineTF1(fsp_ams, fit, HistTools::Divide, "f_ratio", R1, R2); // fit ratio of AMS vs. Combined

	cout << "" << endl;
	//fitter.Result().Print(cout);
	cout << "" << endl;
	TH1D *h_fitres[2];
	//TH1 *h_fiterr[2];
	
	//bestnode.open();

	for (UShort_t i = 0; i < 2; ++i)
	{
  		TH1 *hist = HistTools::ToHist(data[i]);
   		h_fitres[i] = (TH1D *)HistTools::GetResiduals(hist, fit, "_fitres", false, true, true, 5, 1, 0.68, &fitter);
		//h_fiterr[i] = HistTools::GetFitError(hist, fit, "_fiterr", true, false, true, 11, 1.); // this will compute the relative error of the fit, centered in 1 (so 1.01 = 1% upward error; 0.98 = 2% downward error, etc)
		HistTools::CopyStyle(hist, h_fitres[i]);

   		UShort_t ndf  = hist->GetNbinsX();
   		Double_t chi2 = chi2norm[i]*ndf;
   		printf(" %d nodes ### %-34s   chi2/ndf=%6.2f/%-2u   chi2norm=%5.2f   prob.=%5.2f%%\n", nnodes, hist->GetTitle(), chi2, ndf, chi2norm[i], TMath::Prob(chi2, ndf)*1e2);
		
		//bestnode << " " << nnodes << " nodes" << hist->GetTitle() << ", chi2=" << chi2 << ", ndf=" << ndf << ", chi2norm=" << chi2norm[i] << ", prob.=" << TMath::Prob(chi2, ndf)*1e2 << endl;
	}
	
	//bestnode.close();

	cout << " " << endl; 
	
	c5->cd(i+1);

	TH1 *ha = HistTools::CreateAxis("ha", "haxis1", 0.1, 2500., 7, 1e-10, 1e3, false);
	ha->SetXTitle(Unit::GetEnergyLabel("GV"));
  	ha->SetYTitle(Unit::GetDifferentialFluxLabel("GV m")); 
	ha->SetTitle(Form("Averaged ACE and Integrated %s AMS Flux vs. Rigidity", ACE_Element[i]));

	TH1 *ha_res = HistTools::CreateAxis("ha_res", "haxis2", 0.1, 2500., 7, -1, 1, false);
	HistTools::CopyStyle(ha, ha_res);

	TLegend *legend5 = new TLegend(0.1,0.8,0.28,0.9); // left, down, right, top
	legend5->AddEntry(h_ace, Form("ACE %s Flux", ACE_Element[i]), "p");
	legend5->AddEntry(h_ams, Form("AMS Integrated %s Flux", ACE_Element[i]), "p");
	legend5->AddEntry(fit, Form("ACE & AMS Combined %s Flux Fit", ACE_Element[i]), "l"); 
		
	gPad->SetLogy();
	gPad->SetLogx();
	HistTools::SetMarkerStyle(h_ace, HistTools::GetColorPalette(i, 4), kFullCircle, 0.9); 	
	HistTools::SetMarkerStyle(h_ams, kBlue, kFullCircle, 0.9);	 
	HistTools::SetMarkerStyle(h_fitres[0], kBlue, kFullCircle, 0.9);
	HistTools::SetMarkerStyle(h_fitres[1], kBlue, kFullCircle, 0.9);
	ha->Draw("E1X0 SAME");
	h_ace->Draw("E1X0 SAME");
	h_ams->Draw("E1X0 SAME");
	legend5->Draw("SAME");

	TFile file0(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[i], nnodes), "RECREATE"); 
	fit->Write("fit_both");
	h_fitres[0]->Write("h_fitres_ace");
	h_fitres[1]->Write("h_fitres_ams");  
	
	//fit->SetLineColor(kRed);
	//fit->SetLineWidth(1);
	fit->SetRange(0.1,2500.);	
	fit->Draw("SAME");

	c5->cd(i+5);

	gPad->SetLogx();
	ha_res->SetTitle(Form("%s Fit Residuals;Rigidity [GV];", ACE_Element[i]));
	ha_res->Draw("E1X0 SAME");
	h_fitres[0]->Draw("E1X0 SAME");
	h_fitres[1]->Draw("E1X0 SAME");
	
	c5->cd(i+9);
	
	TH1 *ha_ratio = HistTools::CreateAxis("ha_ratio", "haxis_ratio", 0.1, 2500., 7, 0., 2.5, false);
	ha_ratio->SetTitle(Form("AMS %s Fit vs. Combined Fit Ratio;Rigidity [GV];", ACE_Element[i]));
	HistTools::SetMarkerStyle(f_ratio, kRed, kFullCircle, 0.6);
	gPad->SetLogx();
	
	//f_ratio->GetYaxis()->SetRangeUser(0., 1.);
	f_ratio->SetRange(0.1,2500.);
	ha_ratio->Draw("E1X0");
	f_ratio->Draw("SAME");

	file0.Write();
   	file0.Close();

   } // end of BCNO loop

   c5->Print(Form("./data/ACE/compare/compare_BCNO_%dnodes.png", nnodes));

}

void compare_nodes(int k){ 

	TH2D *fit = new TH2D("fit", "ACE/AMS Combined Fit Sigma Test for All Nodes for BCNO",5,0,5,5,0,5); 
	TCanvas *c1 = new TCanvas();
	//c1->Divide(2,2);
	c1->cd(1);

	for (int i=5;i<=9;i++){

	fit->GetXaxis()->SetBinLabel(i-4, Form("%d", i));

	   for (int j=i+1;j<=9;j++){

		TFile file1(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[k], i));
		TFile file2(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[k], j));
	
		TF1 *fit1 = (TF1*) file1.Get("fit_both");
		TF1 *fit2 = (TF1*) file2.Get("fit_both");		
		
		double sigma = compare_sig(fit1, fit2);
		//HistTools::PrintFunction(fit1);	
		printf("%s Fit %d vs. %d, sigma = %0.4f \n", ACE_Element[k], i, j, sigma); 

		fit->SetBinContent(i-4,j-i+1,sigma);
		fit->GetYaxis()->SetBinLabel(j-4, Form("%d", j));
		fit->SetTitle(Form("ACE/AMS Combined Fit Sigma Test for All Nodes for %s", ACE_Element[k]));
		fit->Draw("Colz");

		file1.Close();
		file2.Close();
	   }
	}
}

// make extension assumption for remaining ACE element data that is not measured by AMS
void ace_extend(){

	Debug::Enable(Debug::ALL); 

	const UInt_t FirstACEBR = 2240;
   	vector<UInt_t> BRs;
  	// we stop at BR 2493, which ends on 2016/05/23, just 3 days before the end of the data taking period for AMS nuclei
   	for (UInt_t br=2426; br<=2493; ++br) { 
	   	if (br != 2472 && br != 2473) BRs.push_back(br-FirstACEBR); 
   	}

	for (int i=4; i<24; i++){

		TCanvas *c1 = new TCanvas("c1","f_ratio residuals for remaining elements", 2400, 800);
		c1->Divide(2, 2);

		for (int j=0; j<4; j++){
	
			c1->cd(j+1);

			TH1 *ha_ratio = HistTools::CreateAxis("ha_ratio", Form("Norm %s/%s;Rigidity [GV]; Ratio", ACE_Element[i], ACE_Element[j]), 0.8, 2.5, 7, 0.6, 1.2, false);
			ha_ratio->Draw("E1X0");

			int nnodes = 7;

			TFile file1(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[j], nnodes)); // load combined fit
			TH1 *h_ene = HistTools::GraphToHist(get_ace_average_graph( ACE_Element[i] , &BRs[0], BRs.size() ), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
			TH1 *h_ace = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene, ACE_Isotope[i], "MeV/n cm", "GV m", "_rig"); // load averaged ACE data for the same element in rigidity
	
			Spline *sp_comb = new Spline("sp_comb", nnodes, Spline::LogLog | Spline::PowerLaw);
			TF1 *fsp_comb = sp_comb->GetTF1Pointer();  
			TF1 *fit_comb = (TF1*) file1.Get("fit_both");

			HistTools::CopyParameters(fit_comb, fsp_comb);
			double x1, x2;
			fit_comb->GetRange(x1,x2);
			fsp_comb->SetRange(x1,x2);

			TH1 *h_ratio = (TH1D *)HistTools::GetResiduals(h_ace, fit_comb, "_ratio", false, true, true, 4, 1);
			HistTools::SetStyle(h_ratio, HistTools::GetColorPalette(j, 4), kFullCircle, 0.9, 1, 1);

			double ratio_sum=0; // compute average of h_ratio manually  
			for(int k=0;k<14;k++){
				ratio_sum += h_ratio->GetBinContent(k);
				//printf("ratio_sum = %0.6f \n", ratio_sum);
			}
			double ratio_ave = ratio_sum/h_ratio->GetEntries();
			
			//printf("ratio_ave = %0.6f w/ %0.1f Entries \n", ratio_ave, h_ratio->GetEntries());
			//HistTools::PrintFunction(fit_comb);
			
			double scale = 1./ratio_ave;
			//printf("i=%d, j=%d, scale = %0.6f \n", i, j, scale); 
			h_ratio->Scale(scale);

			//h_ratio->Print("range");

			gPad->SetGrid();
			//gPad->SetLogx();

			h_ratio->SetTitle(Form("Scaled ACE %s Data vs. Combined Fit Ratio;Rigidity [GV];", ACE_Element[i]));
			h_ratio->Draw("E1X0 SAME"); 
						
			file1.Close();
		}

		// break; 
		c1->Print(Form("./data/ACE/extend/ACE_extend_residuals_byBCNO_%s.png", ACE_Element[i]));
	} 

	TCanvas *c2 = new TCanvas("c2","f_ratio panel for remaining elements", 3600, 4800);
	c2->Divide(2, 3);

	int nnodes = 7;
	TFile file2(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[1], nnodes)); // load C combined fit

	// ACE Z>8 element data over C Combined Fit Normalized Ratio for varies isotopes 
	for (int i=4; i<24; i++){
	
		c2->cd(1);

		gPad->SetGrid();
		gPad->SetLogy();
		gPad->SetLogx();

		TH1 *ha = HistTools::CreateAxis("ha", Form("ACE %s Flux & C Combined Fit", ACE_Element[i]), 0.1, 2500., 7, 1e-10, 1e2, false);
		ha->SetXTitle(Unit::GetEnergyLabel("GV"));
  		ha->SetYTitle(Unit::GetDifferentialFluxLabel("GV m"));
		ha->Draw("E1X0");

		c2->cd(2);

		TH1 *ha_ratio = HistTools::CreateAxis("ha_ratio", Form("Normalized %s/%s;Rigidity [GV]; Ratio", ACE_Element[i], ACE_Element[1]), 0.8, 2.4, 7, 0.7, 1.3, false);
		ha_ratio->Draw("E1X0");

		TLegend *legend1 = new TLegend(0.1,0.7,0.2,0.9); // left, down, right, top
		TLegend *legend2 = new TLegend(0.1,0.7,0.48,0.9); // left, down, right, top
		//legend1->SetNColumns(2);
		legend2->SetNColumns(2);

		TH1D *h_chi2 = new TH1D();  
		TH1D *h_var1 = new TH1D();	
		TH1D *h_var1a = new TH1D();
	

		for (int j=0; j<isotope_size[i]; j++){

			TGraph *g_ene = get_ace_average_graph( ACE_Element[i] , &BRs[0], BRs.size() ); 
			TH1 *h_ene = HistTools::GraphToHist(get_ace_average_graph( ACE_Element[i] , &BRs[0], BRs.size() ), DBL_MIN, -DBL_MAX, true, 0.5, 0.); 

			double gx, gy; 

			//cout << " g_ene: " << endl;
			//for(int k=0;k<7;k++){ 
			//	g_ene->GetPoint(k, gx, gy);
			//	printf("%0.18f \n", gy); 
			//}
			//PRINT_GRAPH(g_ene);
			//cout << " " << endl;
			//cout << " h_ene: " << endl;  
			//PRINT_HIST(h_ene);

			//for(int k=0;k<13;k++){ 
			//	printf("%0.18f \n", h_ene->GetBinContent(k+1)); 
			//}

			TH1 *h_ace = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene, compare_isotope[i][j], "MeV/n cm", "GV m", "_rig"); // load averaged ACE data for the same element in rigidity
			HistTools::SetStyle(h_ace, kBlue, kFullCircle, 0.9, 1, 1); 

			Spline *sp_comb = new Spline("sp_comb", nnodes, Spline::LogLog | Spline::PowerLaw);
			TF1 *fsp_comb = sp_comb->GetTF1Pointer();  
			TF1 *fit_comb = (TF1*) file2.Get("fit_both");

			HistTools::CopyParameters(fit_comb, fsp_comb);
			double x1, x2;
			fit_comb->GetRange(x1,x2);
			fsp_comb->SetRange(x1,x2);
	
			c2->cd(1);

			HistTools::SetStyle(h_ace, HistTools::GetColorPalette(j, isotope_size[i]), kFullCircle, 0.9, 1, 1);
			legend1->AddEntry(h_ace, Form("%s", name_isotope[i][j].c_str())); 
			
			h_ace->Draw("E1X0 SAME");
			fit_comb->Draw("SAME");
			legend1->Draw("SAME");

			TH1 *h_res = (TH1D *) HistTools::GetResiduals(h_ace, fit_comb, "_ratio", false, true, true, 4, 1);
			TH1 *h_ratio = (TH1 *) h_ace->Clone("h_ratio");

			h_ratio->Divide(fit_comb);

			HistTools::SetStyle(h_ratio, kRed, kFullCircle, 0.9, 1, 1);

			double ratio_sum=0; // compute average of h_ratio manually  
			for(int k=0;k<14;k++){
				ratio_sum += h_ratio->GetBinContent(k);
				//printf("ratio_sum = %0.6f \n", ratio_sum);
			}
			double ratio_ave = ratio_sum/7;
			
			//printf("ratio_ave = %0.6f w/ %0.1f Entries \n", ratio_ave, h_ratio->GetEntries());
			//HistTools::PrintFunction(fit_comb);
			
			double scale = 1./ratio_ave;
			//printf("i=%d, j=%d, scale = %0.6f \n", i, j, scale); 
			h_ratio->Scale(scale);	

			//h_ratio->Print("range");

			c2->cd(2); 
			gPad->SetGrid();

			HistTools::SetStyle(h_ratio, HistTools::GetColorPalette(j, isotope_size[i]), kFullCircle, 0.9, 1, 1);
			legend2->AddEntry(h_ratio, Form("%s/%s", name_isotope[i][j].c_str(), ACE_Element[1])); 
			legend2->Draw("SAME");

			h_ratio->SetTitle(Form("Scaled Ratio between ACE %s Data & Combined Fit;Rigidity [GV];", ACE_Element[i]));
			h_ratio->Draw("E1X0 SAME"); 

			c2->cd(3);
			gPad->SetGrid();

			int Nbins = 7;

			double mu_abs=0; // average relative absolute difference
			for(int k=0;k<14;k++){
				if(k%2==0) mu_abs += abs(h_ratio->GetBinContent(k+1)-1); 
			
			}
			mu_abs = mu_abs/Nbins; 

			h_var1->SetBinContent(j+1, mu_abs); 
			h_var1->GetXaxis()->SetBinLabel(j+1, name_isotope[i][j].c_str());  
			HistTools::SetStyle(h_var1, kRed, kFullCircle, 0.9, 1, 1);
			h_var1->SetTitle(Form("%s/%s Average Absolute Variation;; (sum of abs(data-1))/N", ACE_Element[i], ACE_Element[1])); 
			h_var1->Draw("HIST P"); 

			c2->cd(4);
			gPad->SetGrid();

			double std_abs=0; 

			for(int k=0;k<14;k++){
				if(k%2==0) std_abs += pow(h_ratio->GetBinContent(k+1)-1-mu_abs, 2); 
				printf("data = %0.7f, mu_abs = %0.7f, std_abs = %0.7f \n", h_ratio->GetBinContent(k+1), mu_abs, std_abs);  
			}

			std_abs = std_abs/(Nbins-1); 
			std_abs = sqrt(std_abs); 

			h_var1a->SetBinContent(j+1, std_abs);
			h_var1a->GetXaxis()->SetBinLabel(j+1, name_isotope[i][j].c_str());  
			HistTools::SetStyle(h_var1a, kRed, kFullCircle, 0.9, 1, 1);
			h_var1a->SetTitle(Form("%s/%s Average Absolute Variation STD;; std_abs", ACE_Element[i], ACE_Element[1])); 
			h_var1a->Draw("HIST P"); 

			c2->cd(5);
			gPad->SetGrid();

			TF1 *f_flat = new TF1("f_flat","[0]"); // create flat ratio line   
			f_flat->SetParameter(0, 1);  	

			h_ratio->Fit(f_flat, "NQ"); 

			double chi2_flat=0; 

			for(int k=0;k<14;k++){
				if(k%2==0) chi2_flat += pow((h_ratio->GetBinContent(k+1)-1)/h_ratio->GetBinError(k+1), 2); 
				//printf("chi2 = %0.7f, data = %0.7f, error = %0.7f \n", chi2_flat, h_ratio->GetBinContent(k+1), h_ratio->GetBinError(k+1)); 
			}

			h_chi2->SetBinContent(j+1, chi2_flat); 

			//h_chi2->Print("range");

			h_chi2->GetXaxis()->SetBinLabel(j+1, name_isotope[i][j].c_str());  
			HistTools::SetStyle(h_chi2, kRed, kFullCircle, 0.9, 1, 1);
			h_chi2->SetTitle(Form("%s/%s Chi-2 Test;;sum of ((data-1)/error)^2", ACE_Element[i], ACE_Element[1])); 
			h_chi2->Draw("HIST P"); 

			c2->cd(6);
			gPad->SetGrid();

		}

		c2->Print(Form("./data/ACE/extend/ACE_extend_ratio_byC_%s.png", ACE_Element[i]));

	} 
	
	file2.Close();	

}

void ace_extend2(){

	const UInt_t FirstACEBR = 2240;
   	vector<UInt_t> BRs;
  	// we stop at BR 2493, which ends on 2016/05/23, just 3 days before the end of the data taking period for AMS nuclei
   	for (UInt_t br=2426; br<=2493; ++br) { 
	   	if (br != 2472 && br != 2473) BRs.push_back(br-FirstACEBR); 
   	}

	TCanvas *c8 = new TCanvas("c8","f_ratio residuals for remaining elements", 2400, 900);
	c8->Divide(2, 2);

	TLegend *legend6 = new TLegend(0.1,0.8,0.24,0.9); // left, down, right, top

	c8->cd(1);

	TH1 *ha = HistTools::CreateAxis("ha", "ACE BCNO Flux", 0.1, 1450., 7, 1e-10, 1e1, false);

	for (int i=0; i<4; i++){ 

		int nnodes = 7;

		TFile file1(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[i], nnodes)); // load combined fit
		TH1 *h_ene = HistTools::GraphToHist(get_ace_average_graph( ACE_Element[i] , &BRs[0], BRs.size() ), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
		TH1 *h_ace = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene, ACE_Isotope[i], "MeV/n cm", "GV m", "_rig"); // load averaged ACE data for the same element in rigidity
	
		Spline *sp_comb = new Spline("sp_comb", nnodes, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_comb = sp_comb->GetTF1Pointer();  
		TF1 *fit_comb = (TF1*) file1.Get("fit_both");

		HistTools::CopyParameters(fit_comb, fsp_comb);
		double x1, x2;
		fit_comb->GetRange(x1,x2);
		fsp_comb->SetRange(x1,x2);

		//HistTools::PrintFunction(fit_comb);
		gPad->SetLogy();
		gPad->SetLogx();
		gPad->SetGrid();

		legend6->AddEntry(fit_comb, Form("ACE %s Combined Fit", ACE_Element[i])); 

		fit_comb->SetLineColor(HistTools::GetColorPalette(i, 4));
		ha->Draw("E1X0 SAME"); 
		fit_comb->Draw("SAME");
		legend6->Draw("SAME");

	}


	int k=0; 
	for (int j=0; j<4; j++){

	   if (j!=1){

		c8->cd(k+2); 
		TH1 *ha_ratio = HistTools::CreateAxis("ha_ratio", Form("%s/%s;Rigidity [GV];", ACE_Element[j], ACE_Element[1]), 0.8, 1.45, 7, 0., 1., false);
		ha_ratio->Draw("E1X0");

		int nnodes = 7;

		TFile file1(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[j], nnodes)); // load combined fit
		TH1 *h_ene = HistTools::GraphToHist(get_ace_average_graph( ACE_Element[j] , &BRs[0], BRs.size() ), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
		TH1 *h_ace = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene, ACE_Isotope[j], "MeV/n cm", "GV m", "_rig"); // load averaged ACE data for the same element in rigidity

		//h_ace->Print("range");
	
		Spline *sp_comb = new Spline("sp_comb", nnodes, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_comb = sp_comb->GetTF1Pointer();  
		TF1 *fit_comb = (TF1*) file1.Get("fit_both");

		HistTools::CopyParameters(fit_comb, fsp_comb);
		double x1, x2;
		fit_comb->GetRange(x1,x2); 
		fsp_comb->SetRange(x1,x2);

		TFile file2(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[1], nnodes)); // load C combined fit
		TH1 *h_ene_C = HistTools::GraphToHist(get_ace_average_graph( ACE_Element[1] , &BRs[0], BRs.size() )); 
		TH1 *h_ace_C = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene_C, ACE_Isotope[1], "MeV/n cm", "GV m", "_rig"); // load averaged ACE data for the same element in rigidity
	
		//h_ace_C->Print("range");

		Spline *sp_comb_C = new Spline("sp_comb_C", nnodes, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_comb_C = sp_comb_C->GetTF1Pointer();  
		TF1 *fit_comb_C = (TF1*) file2.Get("fit_both");

		HistTools::CopyParameters(fit_comb_C, fsp_comb_C);
		double x1_C, x2_C;
		fit_comb_C->GetRange(x1_C,x2_C);
		fsp_comb_C->SetRange(x1_C,x2_C);

		double R1 = h_ace_C->GetBinLowEdge(1);
		double R2 = h_ace->GetBinLowEdge(h_ace->GetNbinsX()+1);
 
		TF1 *f_ratio = HistTools::CombineTF1(fit_comb, fit_comb_C, HistTools::Divide, "f_ratio", R1, R2); // fit ratio of Combined vs. Combined C

		HistTools::PrintFunction(fit_comb);
		HistTools::PrintFunction(fit_comb_C);
		HistTools::PrintFunction(f_ratio);

		gPad->SetGrid();
		//gPad->SetLogx();

		f_ratio->Draw("SAME"); 	
					
		k++;	

		//break; 
	   } else if (j==1) continue; 

	   //break;
	} 	

	c8->Print("./data/ACE/extend2/BNOvsC.png");

}

void ace_extend3(){

	Debug::Enable(Debug::ALL); 

	const UInt_t FirstACEBR = 2240;
   	vector<UInt_t> BRs;
  	// we stop at BR 2493, which ends on 2016/05/23, just 3 days before the end of the data taking period for AMS nuclei
   	for (UInt_t br=2426; br<=2493; ++br) { 
	   	if (br != 2472 && br != 2473) BRs.push_back(br-FirstACEBR); 
	}

	string group[5][6] = {  { "Al", "Ar", "Ca", "Cr", "Na" },
				{ "Cl", "K", "P" }, 
				{ "Co", "Mn", "Ni", "S" }, 
				{ "Mg", "Ne", "Si", "F" }, 
				{ "Fe", "Sc", "Ti", "Va" } }; 

	Particle::Type group_iso[5][6] = { { Particle::ALUMINUM27, Particle::ARGON36, Particle::CALCIUM40, Particle::CHROMIUM52, Particle::SODIUM23 }, 
					   { Particle::CHLORINE35, Particle::POTASSIUM41, Particle::PHOSPHORUS31 }, 
					   { Particle::COBALT59, Particle::MANGANESE55, Particle::NICKEL60, Particle::SULFUR32 }, 
					   { Particle::MAGNESIUM24, Particle::NEON20, Particle::SILICON28, Particle::FLUORINE19 }, 
					   { Particle::IRON56, Particle::SCANDIUM45, Particle::TITANIUM46, Particle::VANADIUM51 } }; 
	
	int size[5] = { 5, 3, 4, 4, 4 }; // update every time when the line above is changed 
	
	string group_abund[5][6] = {  { "Al27", "Ar36", "Ca40", "Cr52", "Na23" },
				{ "Cl35", "K41", "P31" }, 
				{ "Co59", "Mn55", "Ni60", "S32" }, 
				{ "Mg24", "Ne20", "Si28", "F19" }, 
				{ "Fe56", "Sc45", "Ti46", "Va51" } };

	TCanvas *c1 = new TCanvas("c1","f_ratio for Z>8 elements with similar shape of ratio", 2400, 900);
	c1->Divide(3, 2);

	// group similar ratios together 
 	for (int i=0; i<5; i++){

		c1->cd(i+1);

		TH1 *ha_ratio = HistTools::CreateAxis("ha_ratio", Form("Normalized Flux/C Ratio Shape #%d;Rigidity [GV];", i+1), 0.8, 2.5, 7, 0.6, 1.2, false);
		ha_ratio->Draw("E1X0"); 

		TLegend *legend1 = new TLegend(0.1,0.8,0.24,0.9); 

		for (int j=0; j<size[i]; j++){

			int nnodes = 7;

			TFile file1(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[1], nnodes)); // load combined C fit
			TH1 *h_ene = HistTools::GraphToHist(get_ace_average_graph( group[i][j].c_str() , &BRs[0], BRs.size() ), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
			TH1 *h_ace = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene, group_iso[i][j], "MeV/n cm", "GV m", "_rig"); // load averaged ACE data for the same element in rigidity
	
			Spline *sp_comb = new Spline("sp_comb", nnodes, Spline::LogLog | Spline::PowerLaw); 
			TF1 *fsp_comb = sp_comb->GetTF1Pointer();  
			TF1 *fit_comb = (TF1*) file1.Get("fit_both");

			HistTools::CopyParameters(fit_comb, fsp_comb);
			double x1, x2;
			fit_comb->GetRange(x1,x2);
			fsp_comb->SetRange(x1,x2);

			TH1 *h_res = (TH1D *)HistTools::GetResiduals(h_ace, fit_comb, "_ratio", false, true, true, 4, 1);
			HistTools::SetStyle(h_res, HistTools::GetColorPalette(j, size[i]), kFullCircle, 0.9, 1, 1);

			double res_sum=0; // compute average of h_res manually  
			for(int k=0;k<14;k++){
				res_sum += h_res->GetBinContent(k);
				//printf("ratio_sum = %0.6f \n", ratio_sum);
			}
			double res_ave = res_sum/h_res->GetEntries();
			
			double scale = 1./res_ave;
			h_res->Scale(scale); 

			legend1->AddEntry(h_res, Form("%s/C", group[i][j].c_str() )); 

			gPad->SetGrid();
			h_res->Draw("E1X0 SAME");

		}
		//break; 
		legend1->Draw("SAME");
	} 

	c1->Print("./data/ACE/extend/ACE_extend_group_ratio_dividebyC.png");

	// group similar ratios together 
 	for (int i=0; i<5; i++){

		TCanvas *c2 = new TCanvas("c2","f_ratio for Z>8 elements with similar shape of ratio", 2400, 900);
		c2->Divide(3, 2);

		c2->cd(1); 

		TH1 *ha_ratio = HistTools::CreateAxis("ha_ratio", Form("Normalized Flux/C Ratio Shape #%d;Rigidity [GV];", i+1), 0.8, 2.5, 7, 0.6, 1.2, false);
		ha_ratio->Draw("E1X0");

		c2->cd(2); 

		TH1 *ha_res = HistTools::CreateAxis("ha_res", Form("Normalized Flux/C Residual Shape #%d;Rigidity [GV];", i+1), 0.8, 2.5, 7, 0.6, 1.2, false);
		ha_res->Draw("E1X0");

		TLegend *legend1 = new TLegend(0.1,0.8,0.24,0.9); 
		TLegend *legend2 = new TLegend(0.1,0.8,0.24,0.9); 

		TH1D *h_chi2 = new TH1D();  
		TH1D *h_var1 = new TH1D();	
		TH1D *h_var1a = new TH1D();

		for (int j=0; j<size[i]; j++){

			int nnodes = 7; 

			TFile file1(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[1], nnodes)); // load combined fit
			TH1 *h_ene = HistTools::GraphToHist(get_ace_average_graph( group[i][j].c_str() , &BRs[0], BRs.size() ), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
			TH1 *h_ace = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene, group_iso[i][j], "MeV/n cm", "GV m", "_rig"); // load averaged ACE data for the same element in rigidity
	
			Spline *sp_comb = new Spline("sp_comb", nnodes, Spline::LogLog | Spline::PowerLaw); 
			TF1 *fsp_comb = sp_comb->GetTF1Pointer();  
			TF1 *fit_comb = (TF1*) file1.Get("fit_both");

			HistTools::CopyParameters(fit_comb, fsp_comb);
			double x1, x2;
			fit_comb->GetRange(x1,x2);
			fsp_comb->SetRange(x1,x2);

			TH1 *h_ratio = (TH1 *) h_ace->Clone("h_ratio");
			HistTools::SetStyle(h_ratio, HistTools::GetColorPalette(j, size[i]), kFullCircle, 0.9, 1, 1);
			h_ratio->Divide(fit_comb);
			legend1->AddEntry(h_ratio, Form("%s/C", group_abund[i][j].c_str() ));

			TH1 *h_res = (TH1D *)HistTools::GetResiduals(h_ace, fit_comb, "_ratio", false, true, true, 4, 1);
			HistTools::SetStyle(h_res, HistTools::GetColorPalette(j, size[i]), kFullCircle, 0.9, 1, 1);
			legend2->AddEntry(h_res, Form("%s/C", group_abund[i][j].c_str() ));

		//	PRINT_HIST(h_res);

			double ratio_sum=0; // compute average of h_res manually  
			for(int k=0;k<14;k++){
				ratio_sum += h_ratio->GetBinContent(k);
				printf("ratio_sum = %0.6f, # of bins = %f \n", ratio_sum, h_res->GetEntries()); // after divide by the fit, the entries of h_ratio changes  
			}
			double ratio_ave = ratio_sum/7;

			PRINT_HIST(h_ratio);
		
			h_ratio->Scale(1./ratio_ave); 

			double res_sum=0; // compute average of h_ratio manually  
			for(int k=0;k<14;k++){
				res_sum += h_res->GetBinContent(k);
				//printf("ratio_sum = %0.6f \n", ratio_sum);
			}
			double res_ave = res_sum/h_res->GetEntries();
			
			//printf("ratio_ave = %0.6f w/ %0.1f Entries \n", ratio_ave, h_ratio->GetEntries());
			//HistTools::PrintFunction(fit_comb);

			//printf("i=%d, j=%d, scale = %0.6f \n", i, j, scale); 
			h_res->Scale(1./res_ave); 
		
			c2->cd(1); 
			gPad->SetGrid();	

			h_ratio->Draw("E1X0 SAME");
			legend1->Draw("SAME");		

			c2->cd(2);
			gPad->SetGrid();

			h_res->Draw("E1X0 SAME");
			legend2->Draw("SAME"); 

			c2->cd(4);
			gPad->SetGrid();
	
			int Nbins = 7;

			double mu_abs=0; // average relative absolute difference
			for(int k=0;k<14;k++){
				if(k%2==0) mu_abs += abs(h_ratio->GetBinContent(k+1)-1); 
			
			}
			mu_abs = mu_abs/Nbins; 

			h_var1->SetBinContent(j+1, mu_abs); 
			h_var1->GetXaxis()->SetBinLabel(j+1, group_abund[i][j].c_str());  
			HistTools::SetStyle(h_var1, kRed, kFullCircle, 0.9, 1, 1);
			h_var1->SetTitle("Average Absolute Variation;; (sum of abs(data-1))/N"); 
			h_var1->Draw("HIST P"); 

			c2->cd(5);
			gPad->SetGrid();

			double std_abs=0; 

			for(int k=0;k<14;k++){
				if(k%2==0) std_abs += pow(h_ratio->GetBinContent(k+1)-1-mu_abs, 2); 
				//printf("data = %0.7f, mu_abs = %0.7f, std_abs = %0.7f \n", h_ratio->GetBinContent(k+1), mu_abs, std_abs);  
			}

			std_abs = std_abs/(Nbins-1); 
			std_abs = sqrt(std_abs); 

			h_var1a->SetBinContent(j+1, std_abs);
			h_var1a->GetXaxis()->SetBinLabel(j+1, group_abund[i][j].c_str());  
			HistTools::SetStyle(h_var1a, kRed, kFullCircle, 0.9, 1, 1);
			h_var1a->SetTitle("Average Absolute Variation STD;; std_abs"); 
			h_var1a->Draw("HIST P"); 

			c2->cd(6);
			gPad->SetGrid();

			double chi2_flat=0; 

			for(int k=0;k<14;k++){
				if(k%2==0) chi2_flat += pow((h_ratio->GetBinContent(k+1)-1)/h_ratio->GetBinError(k+1), 2); 
				//printf("chi2 = %0.7f, data = %0.7f, error = %0.7f \n", chi2_flat, h_ratio->GetBinContent(k+1), h_ratio->GetBinError(k+1)); 
			}

			h_chi2->SetBinContent(j+1, chi2_flat); 

			h_chi2->GetXaxis()->SetBinLabel(j+1, group_abund[i][j].c_str());  
			HistTools::SetStyle(h_chi2, kRed, kFullCircle, 0.9, 1, 1);
			h_chi2->SetTitle(Form("%s/%s Chi-2 Test;;sum of ((data-1)/error)^2", group[i][j].c_str(), ACE_Element[1])); 
			h_chi2->Draw("HIST P");  
						
			//file1.Close();
		}
 	
		c2->Print(Form("./data/ACE/extend/ACE_extend_group_ratio_dividebyC_shape%d.png", i+1));
	} 

}

void ace_extend4(){

	Debug::Enable(Debug::ALL); 

	const UInt_t FirstACEBR = 2240;
   	vector<UInt_t> BRs;
  	// we stop at BR 2493, which ends on 2016/05/23, just 3 days before the end of the data taking period for AMS nuclei
   	for (UInt_t br=2426; br<=2493; ++br) { 
	   	if (br != 2472 && br != 2473) BRs.push_back(br-FirstACEBR); 
	}

	string group[5][6] = {  { "Al", "Ar", "Ca", "Cr", "Na" },
				{ "Cl", "K", "P" }, 
				{ "Co", "Mn", "Ni", "S" }, 
				{ "Mg", "Ne", "Si", "F" }, 
				{ "Fe", "Sc", "Ti", "Va" } }; 
	
	string group_name[5] = { "Al, Ar, Ca, Cr, Na", "Cl, K, P", "Co, Mn, Ni, S", "Mg, Ne, Si, F", "Fe, Sc, Ti, Va"}; 

	Particle::Type group_iso[5][6] = { { Particle::ALUMINUM27, Particle::ARGON36, Particle::CALCIUM40, Particle::CHROMIUM52, Particle::SODIUM23 }, 
					   { Particle::CHLORINE35, Particle::POTASSIUM41, Particle::PHOSPHORUS31 }, 
					   { Particle::COBALT59, Particle::MANGANESE55, Particle::NICKEL60, Particle::SULFUR32 }, 
					   { Particle::MAGNESIUM24, Particle::NEON20, Particle::SILICON28, Particle::FLUORINE19 }, 
					   { Particle::IRON56, Particle::SCANDIUM45, Particle::TITANIUM46, Particle::VANADIUM51 } }; 
	
	int size[5] = { 5, 3, 4, 4, 4 }; // update every time when the line above is changed 
	
	string group_abund[5][6] = {  { "Al27", "Ar36", "Ca40", "Cr52", "Na23" },
				{ "Cl35", "K41", "P31" }, 
				{ "Co59", "Mn55", "Ni60", "S32" }, 
				{ "Mg24", "Ne20", "Si28", "F19" }, 
				{ "Fe56", "Sc45", "Ti46", "Va51" } };

	int ngroups = 5;  
  
	// over Z
	TGraph *g_rad[ngroups];	
	TGraph *g_radstd[ngroups];
	TGraph *g_chisq[ngroups];
	TGraph *g_maxres[ngroups]; 

	// over A
	TGraph *g_rad1[ngroups];	
	TGraph *g_radstd1[ngroups];
	TGraph *g_chisq1[ngroups];
	TGraph *g_maxres1[ngroups];

	// over Z/A
	TGraph *g_rad2[ngroups];	
	TGraph *g_radstd2[ngroups];
	TGraph *g_chisq2[ngroups];
	TGraph *g_maxres2[ngroups]; 

	// over (Z/A)_element - (Z/A)_C
	TGraph *g_rad3[ngroups];	
	TGraph *g_radstd3[ngroups];
	TGraph *g_chisq3[ngroups];
	TGraph *g_maxres3[ngroups]; 

	TCanvas *c1 = new TCanvas("c1","statistic check for all ratios", 8000, 4500); 
	c1->Divide(4, 4);

	//int startpoint = 4; 

	TH1 *ha1 = HistTools::CreateAxis("ha1", "Average Absolute Variation of All Ratio over C Combined Fit;Z; (sum of abs(data-1))/N", 8, 30, size[0], 0.0, 0.11, false);
	TH1 *ha2 = HistTools::CreateAxis("ha2", "Average Absolute Variation STD of All Ratio over C Combined Fit;Z; std_abs", 8, 30, size[1], 0.0, 0.2, false);
	TH1 *ha3 = HistTools::CreateAxis("ha3", "Chi-2 Test;Z;sum of ((data-1)/error)^2", 8, 30, size[0], 0.0, 50.0, false);
	TH1 *ha4 = HistTools::CreateAxis("ha4", "Max Residual;Z; max residual", 8, 30, size[0], 0.0, 0.2, false); 

	TH1 *ha5 = HistTools::CreateAxis("ha5", "Average Absolute Variation of All Ratio over C Combined Fit;Z; (sum of abs(data-1))/N", 5, 70, size[1], 0.0, 0.11, false);
	TH1 *ha6 = HistTools::CreateAxis("ha6", "Average Absolute Variation STD of All Ratio over C Combined Fit;Z; std_abs", 5, 70, size[1], 0.0, 0.2, false);
	TH1 *ha7 = HistTools::CreateAxis("ha7", "Chi-2 Test;Z;sum of ((data-1)/error)^2", 5, 70, size[1], 0.0, 50.0, false);
	TH1 *ha8 = HistTools::CreateAxis("ha8", "Max Residual;Z; max residual", 5, 70, size[1], 0.0, 0.2, false);

	TH1 *ha9 = HistTools::CreateAxis("ha9", "Average Absolute Variation of All Ratio over C Combined Fit;Z; (sum of abs(data-1))/N", 0.42, 0.54, size[2], 0.0, 0.11, false);
	TH1 *ha10 = HistTools::CreateAxis("ha10", "Average Absolute Variation STD of All Ratio over C Combined Fit;Z; std_abs", 0.42, 0.54, size[2], 0.0, 0.2, false);
	TH1 *ha11 = HistTools::CreateAxis("ha11", "Chi-2 Test;Z;sum of ((data-1)/error)^2", 0.42, 0.54, size[2], 0.0, 50.0, false);
	TH1 *ha12 = HistTools::CreateAxis("ha12", "Max Residual;Z; max residual", 0.42, 0.54, size[2], 0.0, 0.2, false);

	TH1 *ha13 = HistTools::CreateAxis("ha13", "Average Absolute Variation of All Ratio over C Combined Fit;Z; (sum of abs(data-1))/N", -0.07, 0.01, size[3], 0.0, 0.11, false);
	TH1 *ha14 = HistTools::CreateAxis("ha14", "Average Absolute Variation STD of All Ratio over C Combined Fit;Z; std_abs", -0.07, 0.01, size[3], 0.0, 0.2, false);
	TH1 *ha15 = HistTools::CreateAxis("ha15", "Chi-2 Test;Z;sum of ((data-1)/error)^2", -0.07, 0.01, size[3], 0.0, 50.0, false);
	TH1 *ha16 = HistTools::CreateAxis("ha16", "Max Residual;Z; max residual", -0.07, 0.01, size[3], 0.0, 0.2, false);

	TLegend *legend1 = new TLegend(0.1,0.7,0.48,0.9); 
	TLegend *legend2 = new TLegend(0.1,0.7,0.48,0.9); 
	TLegend *legend3 = new TLegend(0.1,0.7,0.48,0.9); 
	TLegend *legend4 = new TLegend(0.1,0.7,0.48,0.9); 
	TLegend *legend5 = new TLegend(0.1,0.7,0.48,0.9); 
	TLegend *legend6 = new TLegend(0.1,0.7,0.48,0.9); 
	TLegend *legend7 = new TLegend(0.1,0.7,0.48,0.9); 
	TLegend *legend8 = new TLegend(0.1,0.7,0.48,0.9); 
	TLegend *legend9 = new TLegend(0.1,0.7,0.48,0.9); 
	TLegend *legend10 = new TLegend(0.1,0.7,0.48,0.9); 
	TLegend *legend11 = new TLegend(0.1,0.7,0.48,0.9); 
	TLegend *legend12 = new TLegend(0.1,0.7,0.48,0.9); 
	TLegend *legend13 = new TLegend(0.1,0.7,0.48,0.9); 
	TLegend *legend14 = new TLegend(0.1,0.7,0.48,0.9); 
	TLegend *legend15 = new TLegend(0.1,0.7,0.48,0.9); 
	TLegend *legend16 = new TLegend(0.1,0.7,0.48,0.9); 

	c1->cd(1); 
	gPad->SetGrid();
	gPad->SetMargin(0.12, 0.08, 0.08, 0.08);
	ha1->Draw("E1X0"); 

	c1->cd(2);
	gPad->SetGrid();
	gPad->SetMargin(0.12, 0.08, 0.08, 0.08);
	ha2->Draw("E1X0");

	c1->cd(3);
	gPad->SetGrid();
	gPad->SetMargin(0.12, 0.08, 0.08, 0.08); 
	ha3->Draw("E1X0");

	c1->cd(4);
	gPad->SetGrid();
	gPad->SetMargin(0.12, 0.08, 0.08, 0.08);
	ha4->Draw("E1X0"); 

	c1->cd(5); 
	gPad->SetGrid();
	gPad->SetMargin(0.12, 0.08, 0.08, 0.08);
	ha5->Draw("E1X0"); 

	c1->cd(6); 
	gPad->SetGrid();
	gPad->SetMargin(0.12, 0.08, 0.08, 0.08);
	ha6->Draw("E1X0"); 

	c1->cd(7); 
	gPad->SetGrid();
	gPad->SetMargin(0.12, 0.08, 0.08, 0.08); 
	ha7->Draw("E1X0"); 

	c1->cd(8); 
	gPad->SetGrid();
	gPad->SetMargin(0.12, 0.08, 0.08, 0.08);
	ha8->Draw("E1X0"); 

	c1->cd(9); 
	gPad->SetGrid();
	gPad->SetMargin(0.12, 0.08, 0.08, 0.08);
	ha9->Draw("E1X0");

	c1->cd(10); 
	gPad->SetGrid();
	gPad->SetMargin(0.12, 0.08, 0.08, 0.08);
	ha10->Draw("E1X0"); 

	c1->cd(11); 
	gPad->SetGrid();
	gPad->SetMargin(0.12, 0.08, 0.08, 0.08);
	ha11->Draw("E1X0"); 

	c1->cd(12); 
	gPad->SetGrid();
	gPad->SetMargin(0.12, 0.08, 0.08, 0.08);
	ha12->Draw("E1X0"); 

	c1->cd(13); 
	gPad->SetGrid();
	gPad->SetMargin(0.12, 0.08, 0.08, 0.08);
	ha13->Draw("E1X0"); 

	c1->cd(14); 
	gPad->SetGrid();
	gPad->SetMargin(0.12, 0.08, 0.08, 0.08);
	ha14->Draw("E1X0"); 

	c1->cd(15); 
	gPad->SetGrid();
	gPad->SetMargin(0.12, 0.08, 0.08, 0.08);
	ha15->Draw("E1X0"); 
	   
	c1->cd(16); 
	gPad->SetGrid();
	gPad->SetMargin(0.12, 0.08, 0.08, 0.08); 
	ha16->Draw("E1X0"); 

	for (int i=0; i<ngroups; i++){
	   
	   g_rad[i] = new TGraph();
	   g_radstd[i] = new TGraph();
	   g_chisq[i] = new TGraph();
	   g_maxres[i] = new TGraph(); 	

	   g_rad1[i] = new TGraph();
	   g_radstd1[i] = new TGraph();
	   g_chisq1[i] = new TGraph();
	   g_maxres1[i] = new TGraph(); 

	   g_rad2[i] = new TGraph();
	   g_radstd2[i] = new TGraph();
	   g_chisq2[i] = new TGraph();
	   g_maxres2[i] = new TGraph(); 
	
	   g_rad3[i] = new TGraph();
	   g_radstd3[i] = new TGraph();
	   g_chisq3[i] = new TGraph();
	   g_maxres3[i] = new TGraph(); 	   

	   for (int j=0; j<size[i]; j++){

		int nnodes = 7; 

		TFile file1(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[1], nnodes)); // load combined C fit
		
		TH1 *h_ene = HistTools::GraphToHist(get_ace_average_graph( group[i][j].c_str(), &BRs[0], BRs.size() ), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
		TH1 *h_ace = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene, group_iso[i][j], "MeV/n cm", "GV m", "_rig"); // load averaged ACE data for the same element in rigidity 
		
		HistTools::SetStyle(h_ace, kBlue, kFullCircle, 0.9, 1, 1); 

		Spline *sp_comb = new Spline("sp_comb", nnodes, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_comb = sp_comb->GetTF1Pointer();  
		TF1 *fit_comb = (TF1*) file1.Get("fit_both"); 

		HistTools::CopyParameters(fit_comb, fsp_comb);
		double x1, x2;
		fit_comb->GetRange(x1,x2);
		fsp_comb->SetRange(x1,x2); 

		TH1 *h_res = (TH1D *) HistTools::GetResiduals(h_ace, fit_comb, "_ratio", false, true, true, 4, 1);
		TH1 *h_ratio = (TH1 *) h_ace->Clone("h_ratio");

		h_ratio->Divide(fit_comb);
		//HistTools::SetStyle(h_ratio, HistTools::GetColorPalette(j, size[i]), kFullCircle, 0.9, 1, 1);

		double ratio_sum=0; // compute average of h_ratio manually  
		for(int k=0;k<14;k++){
			ratio_sum += h_ratio->GetBinContent(k);
			//printf("ratio_sum = %0.6f \n", ratio_sum);
		}
		double ratio_ave = ratio_sum/h_res->GetEntries();
			
		//printf("ratio_ave = %0.6f w/ %0.1f Entries \n", ratio_ave, h_ratio->GetEntries());
		//HistTools::PrintFunction(fit_comb);
			
		double scale = 1./ratio_ave;
		//printf("i=%d, j=%d, scale = %0.6f \n", i, j, scale); 
		h_ratio->Scale(scale);	

		int Nbins = 7;

		double mu_abs=0; // average relative absolute difference
		for(int k=0;k<14;k++){
			if(k%2==0) mu_abs += abs(h_ratio->GetBinContent(k+1)-1); 
			
		}
		mu_abs = mu_abs/Nbins; 

		g_rad[i]->SetPoint(j, Particle::Z[group_iso[i][j]], mu_abs);  
		
		//PRINT_GRAPH(g_rad[i]);

		double std_abs=0; 

		for(int k=0;k<14;k++){
			if(k%2==0) std_abs += pow(h_ratio->GetBinContent(k+1)-1-mu_abs, 2); 
			//printf("data = %0.7f, mu_abs = %0.7f, std_abs = %0.7f \n", h_ratio->GetBinContent(k+1), mu_abs, std_abs);  
		}

		std_abs = std_abs/(Nbins-1); 
		std_abs = sqrt(std_abs); 	

		g_radstd[i]->SetPoint(j, Particle::Z[group_iso[i][j]], std_abs); 

		double chi2_flat=0; 

		for(int k=0;k<14;k++){
			if(k%2==0) chi2_flat += pow((h_ratio->GetBinContent(k+1)-1)/h_ratio->GetBinError(k+1), 2); 
			//printf("chi2 = %0.7f, data = %0.7f, error = %0.7f \n", chi2_flat, h_ratio->GetBinContent(k+1), h_ratio->GetBinError(k+1)); 
		}

		g_chisq[i]->SetPoint(j, Particle::Z[group_iso[i][j]], chi2_flat);  
		//h_chisq->Print("range");	

		double maxres = 0;    
		for(int k=0; k<h_ratio->GetNbinsX(); ++k) {
			if(k%2==0){
				if ( abs(h_ratio->GetBinContent(k+1)-1) > maxres ) maxres = abs(h_ratio->GetBinContent(k+1)-1); 
			}
			//printf("bin = %d, bin value = %0.3f, maxres = %0.3f \n",  k+1, h_ratio->GetBinContent(k+1), maxres);
		} 

		g_maxres[i]->SetPoint(j, Particle::Z[group_iso[i][j]], maxres); 
	
		g_rad1[i]->SetPoint(j, Particle::A[group_iso[i][j]], mu_abs); 
		g_radstd1[i]->SetPoint(j, Particle::A[group_iso[i][j]], std_abs);  
		g_chisq1[i]->SetPoint(j, Particle::A[group_iso[i][j]], chi2_flat);  
		g_maxres1[i]->SetPoint(j, Particle::A[group_iso[i][j]], maxres); 

		g_rad2[i]->SetPoint(j, Particle::Z[group_iso[i][j]]/Particle::A[group_iso[i][j]], mu_abs); 
		g_radstd2[i]->SetPoint(j, Particle::Z[group_iso[i][j]]/Particle::A[group_iso[i][j]], std_abs);  
		g_chisq2[i]->SetPoint(j, Particle::Z[group_iso[i][j]]/Particle::A[group_iso[i][j]], chi2_flat);  
		g_maxres2[i]->SetPoint(j, Particle::Z[group_iso[i][j]]/Particle::A[group_iso[i][j]], maxres); 

		g_rad3[i]->SetPoint(j, Particle::Z[group_iso[i][j]]/Particle::A[group_iso[i][j]]-Particle::Z[ACE_Isotope[1]]/Particle::A[ACE_Isotope[1]], mu_abs); 
		g_radstd3[i]->SetPoint(j, Particle::Z[group_iso[i][j]]/Particle::A[group_iso[i][j]]-Particle::Z[ACE_Isotope[1]]/Particle::A[ACE_Isotope[1]], std_abs);  
		g_chisq3[i]->SetPoint(j, Particle::Z[group_iso[i][j]]/Particle::A[group_iso[i][j]]-Particle::Z[ACE_Isotope[1]]/Particle::A[ACE_Isotope[1]], chi2_flat);  
		g_maxres3[i]->SetPoint(j, Particle::Z[group_iso[i][j]]/Particle::A[group_iso[i][j]]-Particle::Z[ACE_Isotope[1]]/Particle::A[ACE_Isotope[1]], maxres); 
		
	   } 	

	   c1->cd(1); 
	   legend1->AddEntry(g_rad[i], Form("%s", group_name[i].c_str()), "p");
	   legend1->Draw("SAME");

	   HistTools::SetStyle(g_rad[i], HistTools::GetColorPalette(i, size[i]), kFullCircle, 0.9, 1, 1);
	   g_rad[i]->Draw("PSAME");  

	   c1->cd(2);
	   legend2->AddEntry(g_radstd[i], Form("%s ", group_name[i].c_str()), "p");
	   legend2->Draw("SAME");
 
	   HistTools::SetStyle(g_radstd[i], HistTools::GetColorPalette(i, size[i]), kFullCircle, 0.9, 1, 1); 
	   g_radstd[i]->Draw("PSAME"); 

	   c1->cd(3);
	   legend3->AddEntry(g_chisq[i], Form("%s", group_name[i].c_str()), "p");
	   legend3->Draw("SAME");

	   HistTools::SetStyle(g_chisq[i], HistTools::GetColorPalette(i, size[i]), kFullCircle, 0.9, 1, 1); 
	   g_chisq[i]->Draw("PSAME"); 

	   c1->cd(4);
	   legend4->AddEntry(g_maxres[i], Form("%s", group_name[i].c_str()), "p");
	   legend4->Draw("SAME");
   
	   HistTools::SetStyle(g_maxres[i], HistTools::GetColorPalette(i, size[i]), kFullCircle, 0.9, 1, 1);	
	   g_maxres[i]->Draw("PSAME"); 

	   c1->cd(5); 
	   legend5->AddEntry(g_rad1[i], Form("%s", group_name[i].c_str()), "p");
	   legend5->Draw("SAME");

	   HistTools::SetStyle(g_rad1[i], HistTools::GetColorPalette(i, size[i]), kFullCircle, 0.9, 1, 1);
	   g_rad1[i]->Draw("PSAME"); 

	   c1->cd(6);  
	   legend6->AddEntry(g_radstd1[i], Form("%s ", group_name[i].c_str()), "p");
	   legend6->Draw("SAME");

	   HistTools::SetStyle(g_radstd1[i], HistTools::GetColorPalette(i, size[i]), kFullCircle, 0.9, 1, 1); 
	   g_radstd1[i]->Draw("PSAME"); 

	   c1->cd(7); 
	   legend7->AddEntry(g_chisq1[i], Form("%s", group_name[i].c_str()), "p");
	   legend7->Draw("SAME");

	   HistTools::SetStyle(g_chisq1[i], HistTools::GetColorPalette(i, size[i]), kFullCircle, 0.9, 1, 1); 
	   g_chisq1[i]->Draw("PSAME");  

	   c1->cd(8); 
	   legend8->AddEntry(g_maxres1[i], Form("%s", group_name[i].c_str()), "p");
	   legend8->Draw("SAME");

	   HistTools::SetStyle(g_maxres1[i], HistTools::GetColorPalette(i, size[i]), kFullCircle, 0.9, 1, 1);	
	   g_maxres1[i]->Draw("PSAME"); 

	   c1->cd(9); 
	   legend9->AddEntry(g_rad2[i], Form("%s", group_name[i].c_str()), "p");
	   legend9->Draw("SAME");

	   HistTools::SetStyle(g_rad2[i], HistTools::GetColorPalette(i, size[i]), kFullCircle, 0.9, 1, 1);
	   g_rad2[i]->Draw("PSAME");  

	   c1->cd(10);  
	   legend10->AddEntry(g_radstd2[i], Form("%s ", group_name[i].c_str()), "p");
	   legend10->Draw("SAME");

	   HistTools::SetStyle(g_radstd2[i], HistTools::GetColorPalette(i, size[i]), kFullCircle, 0.9, 1, 1); 
	   g_radstd2[i]->Draw("PSAME"); 

	   c1->cd(11);
	   legend11->AddEntry(g_chisq2[i], Form("%s", group_name[i].c_str()), "p");
	   legend11->Draw("SAME"); 
	
	   HistTools::SetStyle(g_chisq2[i], HistTools::GetColorPalette(i, size[i]), kFullCircle, 0.9, 1, 1); 
	   g_chisq2[i]->Draw("PSAME"); 

	   c1->cd(12); 
	   legend12->AddEntry(g_maxres2[i], Form("%s", group_name[i].c_str()), "p");
	   legend12->Draw("SAME"); 

	   HistTools::SetStyle(g_maxres2[i], HistTools::GetColorPalette(i, size[i]), kFullCircle, 0.9, 1, 1);	
	   g_maxres2[i]->Draw("PSAME"); 

	   c1->cd(13);
	   legend13->AddEntry(g_rad3[i], Form("%s", group_name[i].c_str()), "p");
	   legend13->Draw("SAME");  

	   HistTools::SetStyle(g_rad3[i], HistTools::GetColorPalette(i, size[i]), kFullCircle, 0.9, 1, 1);
	   g_rad3[i]->Draw("PSAME");

	   c1->cd(14); 
	   legend14->AddEntry(g_radstd3[i], Form("%s", group_name[i].c_str()), "p");
	   legend14->Draw("SAME");

	   HistTools::SetStyle(g_radstd3[i], HistTools::GetColorPalette(i, size[i]), kFullCircle, 0.9, 1, 1); 
	   g_radstd3[i]->Draw("PSAME"); 

	   c1->cd(15); 
	   legend15->AddEntry(g_chisq3[i], Form("%s", group_name[i].c_str()), "p");
	   legend15->Draw("SAME");

	   HistTools::SetStyle(g_chisq3[i], HistTools::GetColorPalette(i, size[i]), kFullCircle, 0.9, 1, 1); 
	   g_chisq3[i]->Draw("PSAME"); 
	   
	   c1->cd(16);  
	   legend16->AddEntry(g_maxres3[i], Form("%s", group_name[i].c_str()), "p");
	   legend16->Draw("SAME");

	   HistTools::SetStyle(g_maxres3[i], HistTools::GetColorPalette(i, size[i]), kFullCircle, 0.9, 1, 1);	
	   g_maxres3[i]->Draw("PSAME"); 

	} 	
	
	c1->Print(Form("./data/ACE/extend/ACE_extend_last_ratio_test.pdf"));

}

// extend AMS p, He, Li, Be at low energy 
void ace_extend5(int temp){

	Debug::Enable(Debug::ALL); 

	Experiments::DataPath = "data";

	gSystem->mkdir(Form("data/ACE/extend2/%s_temp", ACE_Element[temp]), true);

	int data_value[n_ele] = {0, 18, 23, 27, 31, 20, 43, 22}; 

	const char *AMS_Element2[8] = { "p", "he", "li", "be", "b", "c", "n", "o" }; 
	const char *AMS_Element2_Cap[8] = { "Proton", "He", "Li", "Be", "B", "C", "N", "O" }; 

	const UInt_t FirstACEBR = 2240;
   	vector<UInt_t> BRs;
  	// we stop at BR 2493, which ends on 2016/05/23, just 3 days before the end of the data taking period for AMS nuclei
   	for (UInt_t br=2426; br<=2493; ++br) { 
	   	if (br != 2472 && br != 2473) BRs.push_back(br-FirstACEBR); 
	}

	TH1 *h_ams_C = Experiments::GetMeasurementHistogram(Experiments::AMS02, data_value[1], 0); // load AMS Template data in order to determine Rmax

	TH1 *h_ene_C = HistTools::GraphToHist(get_ace_average_graph( ACE_Element[1] , &BRs[0], BRs.size() ), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
	TH1 *h_ace_C = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene_C, ACE_Isotope[1], "MeV/n cm", "GV m", "_rig"); // load averaged ACE data for the same element in rigidity

	UShort_t namsbins = h_ams_C->GetNbinsX(); 
	UShort_t nacebins = h_ace_C->GetNbinsX();
	
	double R1 = h_ace_C->GetBinLowEdge(1);
	double R2 = h_ams_C->GetBinLowEdge(namsbins+1);

	TCanvas *c1 = new TCanvas("c1", "Extend Fit", 2400, 1800); 
	c1->Divide(1, 2); 

	TCanvas *c2 = new TCanvas("c2", "Global Chi2 and AMS Residual vs. Rescaling Factor", 1200, 900); 
	c2->Divide(1, 2); 

	TGraph *g_chi2_rescale = new TGraph(8);
	TGraph *g_residual_rescale = new TGraph(8);  

	for (int i=0; i<4; i++){

	   if (i==0){

		// p, C template 

		int nnodes = 9; 

		int nnodes_ams = 6; 
		int nnodes_ace = 7; 

		TH1 *ha = HistTools::CreateAxis("ha", "haxis1", 0.1, 2500., 7, 1e-10, 1.2e3, false); 
		ha->SetXTitle(Unit::GetEnergyLabel("GV"));
  		ha->SetYTitle(Unit::GetDifferentialFluxLabel("GV m")); 
		ha->SetTitle(Form("Estimated AMS %s Flux vs. Rigidity", AMS_Element2_Cap[i]));

		TH1 *h_ams = Experiments::GetMeasurementHistogram(Experiments::AMS02, data_value[i], 0); // load AMS data
		HistTools::SetStyle(h_ams, kBlue, kFullCircle, 1.5, 1, 1); 
		TFile file1(Form("data/amsfit/fit_result_node%d.root", nnodes_ams)); // load AMS fit  

		TH1 *h_ams_i = (TH1 *) file1.Get(Form("h_%s", AMS_Element2[i])); // AMS p, He, Li, Be
		HistTools::SetStyle(h_ams_i, HistTools::GetColorPalette(i+1, 4), kFullCircle, 1.5, 1, 1);  

		Spline *sp_ams = new Spline("sp_ams", nnodes_ams, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_ams = sp_ams->GetTF1Pointer(); 
		TF1 *fit_ams = (TF1*) file1.Get(Form("fsp_%s", AMS_Element2[i])); // the fsp in the root file is actually fit 

		HistTools::CopyParameters(fit_ams, fsp_ams); 
		double x1, x2;
		fit_ams->GetRange(x1,x2);
		fsp_ams->SetRange(x1,x2); 

		TFile file2(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[temp], nnodes_ace)); // load ACE Combined Template fit
	
		Spline *sp_comb_temp = new Spline("sp_comb_temp", nnodes_ace, Spline::LogLog | Spline::PowerLaw); 
		TF1 *fsp_comb_temp = sp_comb_temp->GetTF1Pointer();  // real function 
		TF1 *fit_comb_temp = (TF1*) file2.Get("fit_both");
	
		HistTools::CopyParameters(fit_comb_temp, fsp_comb_temp); 
		double x1_temp, x2_temp;
		fit_comb_temp->GetRange(x1_temp,x2_temp);
		fsp_comb_temp->SetRange(x1_temp,x2_temp);
		
		double s_ratio = fsp_ams->Eval(1.0)/fsp_comb_temp->Eval(1.0); // ratio of f_ams(1GV)/f_ace(1GV) 
		//printf("%s s_ratio = %0.4f \n", AMS_Element2_Cap[i], s_ratio); 

		TH1 *h_scale[6];  

		h_scale[0] = (TH1 *) h_ace_C->Clone("h_scale"); 			

		HistTools::SetStyle(h_scale[0], HistTools::GetColorPalette(i, 4), kFullCircle, 1.5, 1, 1); 

		// use first 3 bins only for p 	
		for (int k=6;k<14;k++){
			if (k%2==0){
				h_scale[0]->SetBinContent(k+1, 0);
				h_scale[0]->SetBinError(k+1, 0); 
			} 
		} 

		h_scale[0]->Scale(s_ratio);

		// initializes the X and Y position of the spline nodes
		double *xnodes = HistTools::BuildLogBins(R1, R2, nnodes); // xnodes will be an array of nnodes+1 items
		double *ynodes = new double[nnodes+1];
		UShort_t inode;
		for (inode = 0; inode < nnodes+1; ++inode)
		{
   			if (xnodes[inode] > h_scale[0]->GetBinLowEdge(nacebins+1)) break;
   			ynodes[inode] = h_scale[0]->GetBinContent(h_scale[0]->FindBin(xnodes[inode]));
		}
		for (; inode < nnodes+1; ++inode)
		{
	   		double x = xnodes[inode] < h_ams->GetBinLowEdge(1) ? h_ams->GetBinLowEdge(1) : xnodes[inode];
	   		ynodes[inode] = h_ams->GetBinContent(h_ams->FindBin(x));
		}

		// create spline
		Spline *sp = new Spline("sp", nnodes, Spline::LogLog | Spline::PowerLaw, xnodes, ynodes);
		sp->SetSpectralIndices(3, -2.8, 0, 4, -3.5, -1); // set initial values and limits for the spectral index at first and last node
		TF1 *fit = sp->GetTF1Pointer();	
	
		// create array of data to be fitted together
		TObjArray data; 
		data.Add(h_scale[0]);
		data.Add(h_ams);

		// fit AMS and ACE data at the same time
		vector<double> rigmin, rigmax, chi2norm(2);
		ROOT::Fit::Fitter fitter;
		//HistTools::PrintFunction(fit);
		FitTools::SetCommonFitterOptions(fitter);
		FitTools::FitCombinedData(data, fit, "I", rigmin, rigmax, chi2norm, fitter, 3);  

		// data/fit
		TH1D *h_fitres[2];
 		
		for (UShort_t i = 0; i < 2; ++i)
		{
  			TH1 *hist = HistTools::ToHist(data[i]);
   			h_fitres[i] = (TH1D *)HistTools::GetResiduals(hist, fit, "_fitres", false, true, true, 5, 1, 0.68, &fitter);
			HistTools::CopyStyle(hist, h_fitres[i]);

   			UShort_t ndf  = hist->GetNbinsX();
   			Double_t chi2 = chi2norm[i]*ndf;
   			printf(" %d nodes ### %-34s   chi2/ndf=%6.2f/%-2u   chi2norm=%5.2f   prob.=%5.2f%%\n", nnodes, hist->GetTitle(), chi2, ndf, chi2norm[i], TMath::Prob(chi2, ndf)*1e2);
		
		}

		HistTools::SetMarkerStyle(h_fitres[0], kBlue, kFullCircle, 0.9);
		HistTools::SetMarkerStyle(h_fitres[1], kBlue, kFullCircle, 0.9);

		int jmax = 7; // rescaling index 

		if (temp == 0) jmax = 4;
		if (temp == 1) jmax = 2;
		if (temp == 2) jmax = 4;
		if (temp == 3) jmax = 2;

		for (int j=1; j<jmax; j++){ 

			// data/fit
			TH1D *h_fitres2[2];

			h_fitres2[1] = h_fitres[1]; 

			h_scale[j] = (TH1 *) h_scale[j-1]->Clone("h_scale");			
			HistTools::SetStyle(h_scale[j], HistTools::GetColorPalette(i, 4), kFullCircle, 1.5, 1, 1); 

			h_scale[j]->Scale(1.0 + h_fitres2[1]->GetBinContent(1)); 

			// Fit again !! 

			// create spline
			Spline *sp2 = new Spline("sp2", nnodes, Spline::LogLog | Spline::PowerLaw, xnodes, ynodes);
			sp2->SetSpectralIndices(3, -2.8, 0, 4, -3.5, -1); // set initial values and limits for the spectral index at first and last node
			TF1 *fit2 = sp2->GetTF1Pointer();

			// sp2->Print();  

			// create array of data to be fitted together
			TObjArray data2; 
			data2.Add(h_scale[j]);
			data2.Add(h_ams);

			// fit AMS and ACE data at the same time
			FitTools::FitCombinedData(data2, fit2, "I", rigmin, rigmax, chi2norm, fitter, 3); 

			for (UShort_t i = 0; i < 2; ++i)
			{
  				TH1 *hist = HistTools::ToHist(data2[i]);
   				h_fitres2[i] = (TH1D *)HistTools::GetResiduals(hist, fit2, "_fitres", false, true, true, 5, 1, 0.68, &fitter);
				HistTools::CopyStyle(hist, h_fitres2[i]);

   				UShort_t ndf  = hist->GetNbinsX();
   				Double_t chi2 = chi2norm[i]*ndf;
   				printf("new %d nodes ### %-34s   chi2/ndf=%6.2f/%-2u   chi2norm=%5.2f   prob.=%5.2f%%\n", nnodes, hist->GetTitle(), chi2, ndf, chi2norm[i], TMath::Prob(chi2, ndf)*1e2);
		
			}

			HistTools::SetMarkerStyle(h_fitres2[0], kBlue, kFullCircle, 0.9);
			HistTools::SetMarkerStyle(h_fitres2[1], kBlue, kFullCircle, 0.9);

			TLegend *legend = new TLegend(0.62,0.8,0.9,0.9); // left, down, right, top
			legend->AddEntry(h_scale[j], Form("Estimated AMS %s Low Energy Flux", AMS_Element2_Cap[i]), "p");
			legend->AddEntry(h_ams, Form("AMS Integrated %s Flux", AMS_Element2_Cap[i]), "p");
			legend->AddEntry(fit2, Form("Estimated+Actual AMS Combined %s Flux Reconstruction", AMS_Element2_Cap[i]), "l"); 

			// PRINT_HIST(h_scale); 

			c1->cd(1); 
			gPad->SetLogx();
			gPad->SetLogy();
			gPad->SetGrid();	

			ha->SetTitle(Form("Estimated AMS %s Flux vs. Rigidity ( %s Template, chi2/ndf=%6.2f/%d )", AMS_Element2_Cap[i], ACE_Element[temp], fit2->GetChisquare(), fit2->GetNDF()));	

			ha->Draw("E1X0");
			fit2->Draw("SAME"); 
			h_ams->Draw("E1X0 SAME"); 
			h_scale[j]->Draw("E1X0 SAME");
			legend->Draw("SAME"); 

			c1->cd(2); 
			gPad->SetLogx();
			gPad->SetGrid();
	
			TH1 *h_resaxis = HistTools::CreateAxis("h_resaxis", " ; ; Data / Fit - 1", 0.1, 2500., 7, -0.6, 0.6, false);
			HistTools::CopyStyle(ha, h_resaxis); 
			h_resaxis->SetXTitle(Unit::GetEnergyLabel("GV")); 

			h_resaxis->Draw("E1X0"); 
			h_fitres2[0]->Draw("E1X0 SAME");
			h_fitres2[1]->Draw("E1X0 SAME"); 

			// PRINT_HIST(h_fitres[1]) 
			// PRINT_HIST(h_fitres2[1]) 
	
			g_chi2_rescale->SetPoint(j, 1.0+h_fitres2[1]->GetBinContent(1), fit2->GetChisquare()); 
			g_residual_rescale->SetPoint(j, 1.0+h_fitres2[1]->GetBinContent(1), h_fitres2[1]->GetBinContent(1)); 

			printf(" %s rescaled by %d times fit residual = %10.4f \n", AMS_Element2_Cap[i], j, h_fitres2[1]->GetBinContent(1));  
			printf(" %s rescaled by %d times global fit chi2/ndf = %10.4f/%d \n", AMS_Element2_Cap[i], j, fit2->GetChisquare(), fit2->GetNDF());  	

			TFile file0(Form("data/ACE/extend2/fit_%s_temp_%s_%dnodes.root", ACE_Element[temp], AMS_Element2[i], nnodes), "RECREATE"); 
		
			fit2->Write("fsp_p");
			h_fitres2[0]->Write("h_fitres_ace");
			h_fitres2[1]->Write("h_fitres_ams"); 

			file0.Write();
			file0.Close();

			c1->Print(Form("./data/ACE/extend2/%s_temp/extend_low_%s_rescale_%dtimes.png", ACE_Element[temp], AMS_Element2_Cap[i], j));			

		}
  
	   } else if (i==1) { 

		// He 
		int nnodes = 9; 

		int nnodes_ams = 6; 
		int nnodes_ace = 7; 

		TH1 *ha = HistTools::CreateAxis("ha", "haxis1", 0.1, 2500., 7, 1e-10, 1e3, false); 
		ha->SetXTitle(Unit::GetEnergyLabel("GV"));
  		ha->SetYTitle(Unit::GetDifferentialFluxLabel("GV m")); 
		ha->SetTitle(Form("Estimated AMS %s Flux vs. Rigidity", AMS_Element2_Cap[i]));

		TH1 *h_ams = Experiments::GetMeasurementHistogram(Experiments::AMS02, data_value[i], 0); // load AMS data
		HistTools::SetStyle(h_ams, kBlue, kFullCircle, 1.5, 1, 1); 
		TFile file1(Form("data/amsfit/fit_result_node%d.root", nnodes_ams)); // load AMS fit  

		TH1 *h_ams_i = (TH1 *) file1.Get(Form("h_%s", AMS_Element2[i])); // AMS p, He, Li, Be
		HistTools::SetStyle(h_ams_i, HistTools::GetColorPalette(i+1, 4), kFullCircle, 1.5, 1, 1);  

		Spline *sp_ams = new Spline("sp_ams", nnodes_ams, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_ams = sp_ams->GetTF1Pointer(); 
		TF1 *fit_ams = (TF1*) file1.Get(Form("fsp_%s", AMS_Element2[i])); // the fsp in the root file is actually fit 

		HistTools::CopyParameters(fit_ams, fsp_ams); 
		double x1, x2;
		fit_ams->GetRange(x1,x2);
		fsp_ams->SetRange(x1,x2); 

		TFile file2(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[temp], nnodes_ace)); // load ACE Combined Template fit
	
		Spline *sp_comb_temp = new Spline("sp_comb_temp", nnodes_ace, Spline::LogLog | Spline::PowerLaw); 
		TF1 *fsp_comb_temp = sp_comb_temp->GetTF1Pointer();  // real function 
		TF1 *fit_comb_temp = (TF1*) file2.Get("fit_both");
	
		HistTools::CopyParameters(fit_comb_temp, fsp_comb_temp); 
		double x1_temp, x2_temp;
		fit_comb_temp->GetRange(x1_temp,x2_temp);
		fsp_comb_temp->SetRange(x1_temp,x2_temp);
		
		double s_ratio = fsp_ams->Eval(2.0)/fsp_comb_temp->Eval(2.0); // ratio of f_ams(2GV)/f_ace(2GV) 
		//printf("%s s_ratio = %0.4f \n", AMS_Element2_Cap[i], s_ratio);

		TH1 *h_scale[10];  

		h_scale[0] = (TH1 *) h_ace_C->Clone("h_scale");

		h_scale[0]->Scale(s_ratio); 			

		HistTools::SetStyle(h_scale[0], HistTools::GetColorPalette(i, 4), kFullCircle, 1.5, 1, 1); 

		// initializes the X and Y position of the spline nodes
		double *xnodes = HistTools::BuildLogBins(R1, R2, nnodes-1); // xnodes will be an array of nnodes+1 items
		xnodes[1] = 1;
		// xnodes[8] = 2.5e3;  

		double *ynodes = new double[nnodes+1];
		UShort_t inode;
		for (inode = 0; inode < nnodes+1; ++inode)
		{
   			if (xnodes[inode] > h_scale[0]->GetBinLowEdge(nacebins+1)) break;
   			ynodes[inode] = h_scale[0]->GetBinContent(h_scale[0]->FindBin(xnodes[inode]));
		}
		for (; inode < nnodes+1; ++inode)
		{
	   		double x = xnodes[inode] < h_ams->GetBinLowEdge(1) ? h_ams->GetBinLowEdge(1) : xnodes[inode];
	   		ynodes[inode] = h_ams->GetBinContent(h_ams->FindBin(x));
		}

		// create spline
		Spline *sp = new Spline("sp", nnodes, Spline::LogLog | Spline::PowerLaw, xnodes, ynodes);
		sp->SetSpectralIndices(3, -2.8, 0, 4, -3.5, -1); // set initial values and limits for the spectral index at first and last node
		TF1 *fit = sp->GetTF1Pointer();	
	
		// create array of data to be fitted together
		TObjArray data; 
		data.Add(h_scale[0]);
		data.Add(h_ams);

		PRINT_HIST(h_scale[0])
		PRINT_HIST(h_ams); 

		// fit AMS and ACE data at the same time
		vector<double> rigmin, rigmax, chi2norm(2);
		ROOT::Fit::Fitter fitter;
		HistTools::PrintFunction(fit);
		FitTools::SetCommonFitterOptions(fitter);

		FitTools::FitCombinedData(data, fit, "I", rigmin, rigmax, chi2norm, fitter, 3); 

		// break;

		// data/fit
		TH1D *h_fitres[2];
 		
		for (UShort_t i = 0; i < 2; ++i)
		{
  			TH1 *hist = HistTools::ToHist(data[i]);
   			h_fitres[i] = (TH1D *)HistTools::GetResiduals(hist, fit, "_fitres", false, true, true, 5, 1, 0.68, &fitter);
			HistTools::CopyStyle(hist, h_fitres[i]);

   			UShort_t ndf  = hist->GetNbinsX();
   			Double_t chi2 = chi2norm[i]*ndf;
   			printf("old %d nodes ### %-34s   chi2/ndf=%6.2f/%-2u   chi2norm=%5.2f   prob.=%5.2f%%\n", nnodes, hist->GetTitle(), chi2, ndf, chi2norm[i], TMath::Prob(chi2, ndf)*1e2);
		
		} 

		HistTools::SetMarkerStyle(h_fitres[0], kBlue, kFullCircle, 0.9);  
		HistTools::SetMarkerStyle(h_fitres[1], kBlue, kFullCircle, 0.9);  

		g_chi2_rescale->SetPoint(0, 1, fit->GetChisquare()); 
		g_residual_rescale->SetPoint(0, 1, h_fitres[1]->GetBinContent(1)); 

		int jmax = 7; // rescaling index 

		if (temp == 0) jmax = 11; 
		if (temp == 1) jmax = 2;
		if (temp == 2) jmax = 11; 
		if (temp == 3) jmax = 2;

		for (int j=1; j<jmax; j++){ 

			// data/fit
			TH1D *h_fitres2[2];

			h_fitres2[1] = h_fitres[1]; 

			h_scale[j] = (TH1 *) h_scale[j-1]->Clone("h_scale");			
			HistTools::SetStyle(h_scale[j], HistTools::GetColorPalette(i, 4), kFullCircle, 1.5, 1, 1); 

			h_scale[j]->Scale(1.0 + h_fitres2[1]->GetBinContent(1)); 

			// Fit again !! 

			// create spline
			Spline *sp2 = new Spline("sp2", nnodes, Spline::LogLog | Spline::PowerLaw, xnodes, ynodes);
			sp2->SetSpectralIndices(3, -2.8, 0, 4, -3.5, -1); // set initial values and limits for the spectral index at first and last node
			TF1 *fit2 = sp2->GetTF1Pointer();

			// sp2->Print();  

			// create array of data to be fitted together
			TObjArray data2; 
			data2.Add(h_scale[j]);
			data2.Add(h_ams);

			// fit AMS and ACE data at the same time
			FitTools::FitCombinedData(data2, fit2, "I", rigmin, rigmax, chi2norm, fitter, 3); 

			for (UShort_t i = 0; i < 2; ++i)
			{
  				TH1 *hist = HistTools::ToHist(data2[i]);
   				h_fitres2[i] = (TH1D *)HistTools::GetResiduals(hist, fit2, "_fitres", false, true, true, 5, 1, 0.68, &fitter);
				HistTools::CopyStyle(hist, h_fitres2[i]);

   				UShort_t ndf  = hist->GetNbinsX();
   				Double_t chi2 = chi2norm[i]*ndf;
   				printf("new %d nodes ### %-34s   chi2/ndf=%6.2f/%-2u   chi2norm=%5.2f   prob.=%5.2f%%\n", nnodes, hist->GetTitle(), chi2, ndf, chi2norm[i], TMath::Prob(chi2, ndf)*1e2);
		
			}

			HistTools::SetMarkerStyle(h_fitres2[0], kBlue, kFullCircle, 0.9);
			HistTools::SetMarkerStyle(h_fitres2[1], kBlue, kFullCircle, 0.9);

			TLegend *legend = new TLegend(0.62,0.8,0.9,0.9); // left, down, right, top
			legend->AddEntry(h_scale[j], Form("Estimated AMS %s Low Energy Flux", AMS_Element2_Cap[i]), "p");
			legend->AddEntry(h_ams, Form("AMS Integrated %s Flux", AMS_Element2_Cap[i]), "p");
			legend->AddEntry(fit2, Form("Estimated+Actual AMS Combined %s Flux Reconstruction", AMS_Element2_Cap[i]), "l"); 

			// PRINT_HIST(h_scale); 

			c1->cd(1); 
			gPad->SetLogx();
			gPad->SetLogy();
			gPad->SetGrid();	

			ha->SetTitle(Form("Estimated AMS %s Flux vs. Rigidity ( %s Template, chi2/ndf=%6.2f/%d )", AMS_Element2_Cap[i], ACE_Element[temp], fit2->GetChisquare(), fit2->GetNDF()));	

			ha->Draw("E1X0");
			fit2->Draw("SAME"); 
			h_ams->Draw("E1X0 SAME"); 
			h_scale[j]->Draw("E1X0 SAME");
			legend->Draw("SAME"); 

			c1->cd(2); 
			gPad->SetLogx();
			gPad->SetGrid();
	
			TH1 *h_resaxis = HistTools::CreateAxis("h_resaxis", " ; ; Data / Fit - 1", 0.1, 2500., 7, -0.6, 0.6, false);
			HistTools::CopyStyle(ha, h_resaxis); 
			h_resaxis->SetXTitle(Unit::GetEnergyLabel("GV")); 

			h_resaxis->Draw("E1X0"); 
			h_fitres2[0]->Draw("E1X0 SAME");
			h_fitres2[1]->Draw("E1X0 SAME"); 

			// PRINT_HIST(h_fitres[1]) 
			// PRINT_HIST(h_fitres2[1]) 
	
			g_chi2_rescale->SetPoint(j, 1.0+h_fitres2[1]->GetBinContent(1), fit2->GetChisquare()); 
			g_residual_rescale->SetPoint(j, 1.0+h_fitres2[1]->GetBinContent(1), h_fitres2[1]->GetBinContent(1)); 

			printf(" %s rescaled by %d times fit residual = %10.4f \n", AMS_Element2_Cap[i], j, h_fitres2[1]->GetBinContent(1));  
			printf(" %s rescaled by %d times global fit chi2/ndf = %10.4f/%d \n", AMS_Element2_Cap[i], j, fit2->GetChisquare(), fit2->GetNDF());  	

			TFile file0(Form("data/ACE/extend2/fit_%s_temp_%s_%dnodes.root", ACE_Element[temp], AMS_Element2[i], nnodes), "RECREATE"); 
		
			fit2->Write("fsp_he");
			h_fitres2[0]->Write("h_fitres_ace");
			h_fitres2[1]->Write("h_fitres_ams"); 

			file0.Write();
			file0.Close();

			c1->Print(Form("./data/ACE/extend2/%s_temp/extend_low_%s_rescale_%dtimes.png", ACE_Element[temp], AMS_Element2_Cap[i], j));			

		}

		// PRINT_GRAPH(g_chi2_rescale)
		// PRINT_GRAPH(g_residual_rescale)
		
		c2->cd(1);
		gPad->SetGrid(); 

		HistTools::SetStyle(g_chi2_rescale, kRed, kFullCircle, 1.5, 1, 1); 
		g_chi2_rescale->SetTitle(Form("%s Global Chi2 vs. Rescaling Factor; Rescaling Factor; Chi2", AMS_Element2_Cap[i]));  
		g_chi2_rescale->Draw("AP"); 

		c2->cd(2);
		gPad->SetGrid();
		HistTools::SetStyle(g_residual_rescale, kRed, kFullCircle, 1.5, 1, 1); 
		g_residual_rescale->SetTitle(Form("%s Fitting Residual of 1st Bin of AMS Range vs. Rescaling Factor; Rescaling Factor; Residual", AMS_Element2_Cap[i]));
		g_residual_rescale->Draw("AP");  

		c2->Print(Form("./data/ACE/extend2/%s_temp/extend_chi2_residual_vs_rescaling_factor_%s.png", ACE_Element[temp], AMS_Element2_Cap[i])); 

	   } else if (i==2) { 

		// Li  
		int nnodes = 9; 

		int nnodes_ams = 6; 
		int nnodes_ace = 7;  

		TH1 *ha = HistTools::CreateAxis("ha", "haxis1", 0.1, 2500., 7, 1e-10, 1e3, false); 
		ha->SetXTitle(Unit::GetEnergyLabel("GV"));
  		ha->SetYTitle(Unit::GetDifferentialFluxLabel("GV m")); 
		ha->SetTitle(Form("Estimated AMS %s Flux vs. Rigidity", AMS_Element2_Cap[i]));

		TH1 *h_ams = Experiments::GetMeasurementHistogram(Experiments::AMS02, data_value[i], 0); // load AMS data
		HistTools::SetStyle(h_ams, kBlue, kFullCircle, 1.5, 1, 1); 
		TFile file1(Form("data/amsfit/fit_result_node%d.root", nnodes_ams)); // load AMS fit  

		TH1 *h_ams_i = (TH1 *) file1.Get(Form("h_%s", AMS_Element2[i])); // AMS p, He, Li, Be
		HistTools::SetStyle(h_ams_i, HistTools::GetColorPalette(i+1, 4), kFullCircle, 1.5, 1, 1);  

		Spline *sp_ams = new Spline("sp_ams", nnodes_ams, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_ams = sp_ams->GetTF1Pointer(); 
		TF1 *fit_ams = (TF1*) file1.Get(Form("fsp_%s", AMS_Element2[i])); // the fsp in the root file is actually fit 

		HistTools::CopyParameters(fit_ams, fsp_ams); 
		double x1, x2;
		fit_ams->GetRange(x1,x2);
		fsp_ams->SetRange(x1,x2); 

		TFile file2(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[temp], nnodes_ace)); // load ACE Combined Template fit
	
		Spline *sp_comb_temp = new Spline("sp_comb_temp", nnodes_ace, Spline::LogLog | Spline::PowerLaw); 
		TF1 *fsp_comb_temp = sp_comb_temp->GetTF1Pointer();  // real function 
		TF1 *fit_comb_temp = (TF1*) file2.Get("fit_both");
	
		HistTools::CopyParameters(fit_comb_temp, fsp_comb_temp); 
		double x1_temp, x2_temp;
		fit_comb_temp->GetRange(x1_temp,x2_temp);
		fsp_comb_temp->SetRange(x1_temp,x2_temp);
		
		double s_ratio = fsp_ams->Eval(2.0)/fsp_comb_temp->Eval(2.0); // ratio of f_ams(2GV)/f_ace(2GV) 
		//printf("%s s_ratio = %0.4f \n", AMS_Element2_Cap[i], s_ratio);

		TH1 *h_scale[6];  

		h_scale[0] = (TH1 *) h_ace_C->Clone("h_scale");

		h_scale[0]->Scale(s_ratio); 			

		HistTools::SetStyle(h_scale[0], HistTools::GetColorPalette(i, 4), kFullCircle, 1.5, 1, 1); 

		// initializes the X and Y position of the spline nodes
		double *xnodes = HistTools::BuildLogBins(R1, R2, nnodes-1); // xnodes will be an array of nnodes+1 items
		xnodes[1] = 1;
		// xnodes[8] = 2.5e3;  

		double *ynodes = new double[nnodes+1];
		UShort_t inode;
		for (inode = 0; inode < nnodes+1; ++inode)
		{
   			if (xnodes[inode] > h_scale[0]->GetBinLowEdge(nacebins+1)) break;
   			ynodes[inode] = h_scale[0]->GetBinContent(h_scale[0]->FindBin(xnodes[inode]));
		}
		for (; inode < nnodes+1; ++inode)
		{
	   		double x = xnodes[inode] < h_ams->GetBinLowEdge(1) ? h_ams->GetBinLowEdge(1) : xnodes[inode];
	   		ynodes[inode] = h_ams->GetBinContent(h_ams->FindBin(x));
		}

		// create spline
		Spline *sp = new Spline("sp", nnodes, Spline::LogLog | Spline::PowerLaw, xnodes, ynodes);
		sp->SetSpectralIndices(3, -2.8, 0, 4, -3.5, -1); // set initial values and limits for the spectral index at first and last node
		TF1 *fit = sp->GetTF1Pointer();	
	
		// create array of data to be fitted together
		TObjArray data; 
		data.Add(h_scale[0]);
		data.Add(h_ams);

		// fit AMS and ACE data at the same time
		vector<double> rigmin, rigmax, chi2norm(2);
		ROOT::Fit::Fitter fitter;
		//HistTools::PrintFunction(fit);
		FitTools::SetCommonFitterOptions(fitter);
		FitTools::FitCombinedData(data, fit, "I", rigmin, rigmax, chi2norm, fitter, 3); 

		// data/fit
		TH1D *h_fitres[2];
 		
		for (UShort_t i = 0; i < 2; ++i) 
		{
  			TH1 *hist = HistTools::ToHist(data[i]);
   			h_fitres[i] = (TH1D *)HistTools::GetResiduals(hist, fit, "_fitres", false, true, true, 5, 1, 0.68, &fitter);
			HistTools::CopyStyle(hist, h_fitres[i]);

   			UShort_t ndf  = hist->GetNbinsX();
   			Double_t chi2 = chi2norm[i]*ndf;
   			printf("old %d nodes ### %-34s   chi2/ndf=%6.2f/%-2u   chi2norm=%5.2f   prob.=%5.2f%%\n", nnodes, hist->GetTitle(), chi2, ndf, chi2norm[i], TMath::Prob(chi2, ndf)*1e2);
		
		} 

		HistTools::SetMarkerStyle(h_fitres[0], kBlue, kFullCircle, 0.9);  
		HistTools::SetMarkerStyle(h_fitres[1], kBlue, kFullCircle, 0.9);  

		g_chi2_rescale->SetPoint(0, 1, fit->GetChisquare()); 
		g_residual_rescale->SetPoint(0, 1, h_fitres[1]->GetBinContent(1)); 

		int jmax = 0; // rescaling index 

		if (temp == 0) jmax = 6; 
		if (temp == 1) jmax = 2;
		if (temp == 2) jmax = 6; 
		if (temp == 3) jmax = 2; 

		for (int j=1; j<jmax; j++){ 

			// data/fit
			TH1D *h_fitres2[2];

			h_fitres2[1] = h_fitres[1]; 

			h_scale[j] = (TH1 *) h_scale[j-1]->Clone("h_scale");			
			HistTools::SetStyle(h_scale[j], HistTools::GetColorPalette(i, 4), kFullCircle, 1.5, 1, 1); 

			h_scale[j]->Scale(1.0 + h_fitres2[1]->GetBinContent(1)); 

			// Fit again !! 

			// create spline
			Spline *sp2 = new Spline("sp2", nnodes, Spline::LogLog | Spline::PowerLaw, xnodes, ynodes);
			sp2->SetSpectralIndices(3, -2.8, 0, 4, -3.5, -1); // set initial values and limits for the spectral index at first and last node
			TF1 *fit2 = sp2->GetTF1Pointer();

			// sp2->Print();  

			// create array of data to be fitted together
			TObjArray data2; 
			data2.Add(h_scale[j]);
			data2.Add(h_ams);

			// fit AMS and ACE data at the same time
			FitTools::FitCombinedData(data2, fit2, "I", rigmin, rigmax, chi2norm, fitter, 3); 

			for (UShort_t i = 0; i < 2; ++i)
			{
  				TH1 *hist = HistTools::ToHist(data2[i]);
   				h_fitres2[i] = (TH1D *)HistTools::GetResiduals(hist, fit2, "_fitres", false, true, true, 5, 1, 0.68, &fitter);
				HistTools::CopyStyle(hist, h_fitres2[i]);

   				UShort_t ndf  = hist->GetNbinsX();
   				Double_t chi2 = chi2norm[i]*ndf;
   				printf("new %d nodes ### %-34s   chi2/ndf=%6.2f/%-2u   chi2norm=%5.2f   prob.=%5.2f%%\n", nnodes, hist->GetTitle(), chi2, ndf, chi2norm[i], TMath::Prob(chi2, ndf)*1e2);
		
			}

			HistTools::SetMarkerStyle(h_fitres2[0], kBlue, kFullCircle, 0.9);
			HistTools::SetMarkerStyle(h_fitres2[1], kBlue, kFullCircle, 0.9);

			TLegend *legend = new TLegend(0.62,0.8,0.9,0.9); // left, down, right, top
			legend->AddEntry(h_scale[j], Form("Estimated AMS %s Low Energy Flux", AMS_Element2_Cap[i]), "p");
			legend->AddEntry(h_ams, Form("AMS Integrated %s Flux", AMS_Element2_Cap[i]), "p");
			legend->AddEntry(fit2, Form("Estimated+Actual AMS Combined %s Flux Reconstruction", AMS_Element2_Cap[i]), "l"); 

			// PRINT_HIST(h_scale); 

			c1->cd(1); 
			gPad->SetLogx();
			gPad->SetLogy();
			gPad->SetGrid();		

			ha->SetTitle(Form("Estimated AMS %s Flux vs. Rigidity ( %s Template, chi2/ndf=%6.2f/%d )", AMS_Element2_Cap[i], ACE_Element[temp], fit2->GetChisquare(), fit2->GetNDF()));

			ha->Draw("E1X0");
			fit2->Draw("SAME"); 
			h_ams->Draw("E1X0 SAME"); 
			h_scale[j]->Draw("E1X0 SAME");
			legend->Draw("SAME"); 

			c1->cd(2); 
			gPad->SetLogx();
			gPad->SetGrid();
	
			TH1 *h_resaxis = HistTools::CreateAxis("h_resaxis", " ; ; Data / Fit - 1", 0.1, 2500., 7, -0.6, 0.6, false);
			HistTools::CopyStyle(ha, h_resaxis); 
			h_resaxis->SetXTitle(Unit::GetEnergyLabel("GV")); 

			h_resaxis->Draw("E1X0"); 
			h_fitres2[0]->Draw("E1X0 SAME");
			h_fitres2[1]->Draw("E1X0 SAME"); 

			// PRINT_HIST(h_fitres[1]) 
			// PRINT_HIST(h_fitres2[1]) 
	
			g_chi2_rescale->SetPoint(j, 1.0+h_fitres2[1]->GetBinContent(1), fit2->GetChisquare()); 
			g_residual_rescale->SetPoint(j, 1.0+h_fitres2[1]->GetBinContent(1), h_fitres2[1]->GetBinContent(1)); 

			printf(" %s rescaled by %d times fit residual = %10.4f \n", AMS_Element2_Cap[i], j, h_fitres2[1]->GetBinContent(1));  
			printf(" %s rescaled by %d times global fit chi2/ndf = %10.4f/%d \n", AMS_Element2_Cap[i], j, fit2->GetChisquare(), fit2->GetNDF());

			TFile file0(Form("data/ACE/extend2/fit_%s_temp_%s_%dnodes.root", ACE_Element[temp], AMS_Element2[i], nnodes), "RECREATE"); 
		
			fit2->Write("fsp_li");
			h_fitres2[0]->Write("h_fitres_ace");
			h_fitres2[1]->Write("h_fitres_ams"); 

			file0.Write();
			file0.Close();  	

			c1->Print(Form("./data/ACE/extend2/%s_temp/extend_low_%s_rescale_%dtimes.png", ACE_Element[temp], AMS_Element2_Cap[i], j));			

		}

		// PRINT_GRAPH(g_chi2_rescale)
		// PRINT_GRAPH(g_residual_rescale)
		
		c2->cd(1);
		gPad->SetGrid(); 

		HistTools::SetStyle(g_chi2_rescale, kRed, kFullCircle, 1.5, 1, 1); 
		g_chi2_rescale->SetTitle(Form("%s Global Chi2 vs. Rescaling Factor; Rescaling Factor; Chi2", AMS_Element2_Cap[i]));  
		g_chi2_rescale->Draw("AP"); 

		c2->cd(2);
		gPad->SetGrid();
		HistTools::SetStyle(g_residual_rescale, kRed, kFullCircle, 1.5, 1, 1); 
		g_residual_rescale->SetTitle(Form("%s Fitting Residual of 1st Bin of AMS Range vs. Rescaling Factor; Rescaling Factor; Residual", AMS_Element2_Cap[i]));
		g_residual_rescale->Draw("AP");  

		c2->Print(Form("./data/ACE/extend2/%s_temp/extend_chi2_residual_vs_rescaling_factor_%s.png", ACE_Element[temp], AMS_Element2_Cap[i])); 

	   } else if (i==3){ 

		// Be, C template 
		int nnodes = 9; 

		int nnodes_ams = 6; 
		int nnodes_ace = 7; 

		TH1 *ha = HistTools::CreateAxis("ha", "haxis1", 0.1, 2500., 7, 1e-10, 1e3, false); 
		ha->SetXTitle(Unit::GetEnergyLabel("GV"));
  		ha->SetYTitle(Unit::GetDifferentialFluxLabel("GV m")); 
		ha->SetTitle(Form("Estimated AMS %s Flux vs. Rigidity", AMS_Element2_Cap[i]));

		TH1 *h_ams = Experiments::GetMeasurementHistogram(Experiments::AMS02, data_value[i], 0); // load AMS data
		HistTools::SetStyle(h_ams, kBlue, kFullCircle, 1.5, 1, 1); 
		TFile file1(Form("data/amsfit/fit_result_node%d.root", nnodes_ams)); // load AMS fit  

		TH1 *h_ams_i = (TH1 *) file1.Get(Form("h_%s", AMS_Element2[i])); // AMS p, He, Li, Be
		HistTools::SetStyle(h_ams_i, HistTools::GetColorPalette(i+1, 4), kFullCircle, 1.5, 1, 1);  

		Spline *sp_ams = new Spline("sp_ams", nnodes_ams, Spline::LogLog | Spline::PowerLaw);
		TF1 *fsp_ams = sp_ams->GetTF1Pointer(); 
		TF1 *fit_ams = (TF1*) file1.Get(Form("fsp_%s", AMS_Element2[i])); // the fsp in the root file is actually fit 

		HistTools::CopyParameters(fit_ams, fsp_ams); 
		double x1, x2;
		fit_ams->GetRange(x1,x2);
		fsp_ams->SetRange(x1,x2); 

		TFile file2(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[temp], nnodes_ace)); // load ACE Combined Template fit
	
		Spline *sp_comb_temp = new Spline("sp_comb_temp", nnodes_ace, Spline::LogLog | Spline::PowerLaw); 
		TF1 *fsp_comb_temp = sp_comb_temp->GetTF1Pointer();  // real function 
		TF1 *fit_comb_temp = (TF1*) file2.Get("fit_both");
	
		HistTools::CopyParameters(fit_comb_temp, fsp_comb_temp); 
		double x1_temp, x2_temp;
		fit_comb_temp->GetRange(x1_temp,x2_temp);
		fsp_comb_temp->SetRange(x1_temp,x2_temp);
		
		double s_ratio = fsp_ams->Eval(2.0)/fsp_comb_temp->Eval(2.0); // ratio of f_ams(2GV)/f_ace(2GV) 
		//printf("%s s_ratio = %0.4f \n", AMS_Element2_Cap[i], s_ratio);

		TH1 *h_scale[6];  

		h_scale[0] = (TH1 *) h_ace_C->Clone("h_scale");

		h_scale[0]->Scale(s_ratio); 			

		HistTools::SetStyle(h_scale[0], HistTools::GetColorPalette(i, 4), kFullCircle, 1.5, 1, 1); 

		// initializes the X and Y position of the spline nodes
		double *xnodes = HistTools::BuildLogBins(R1, R2, nnodes-1); // xnodes will be an array of nnodes+1 items
		xnodes[1] = 1; 
		// xnodes[8] = 2.5e3; 

		double *ynodes = new double[nnodes+1];
		UShort_t inode;
		for (inode = 0; inode < nnodes+1; ++inode)
		{
   			if (xnodes[inode] > h_scale[0]->GetBinLowEdge(nacebins+1)) break;
   			ynodes[inode] = h_scale[0]->GetBinContent(h_scale[0]->FindBin(xnodes[inode]));
		}
		for (; inode < nnodes+1; ++inode)
		{
	   		double x = xnodes[inode] < h_ams->GetBinLowEdge(1) ? h_ams->GetBinLowEdge(1) : xnodes[inode];
	   		ynodes[inode] = h_ams->GetBinContent(h_ams->FindBin(x));
		}

		// create spline
		Spline *sp = new Spline("sp", nnodes, Spline::LogLog | Spline::PowerLaw, xnodes, ynodes);
		sp->SetSpectralIndices(3, -2.8, 0, 4, -3.5, -1); // set initial values and limits for the spectral index at first and last node
		TF1 *fit = sp->GetTF1Pointer();	
	
		// create array of data to be fitted together
		TObjArray data; 
		data.Add(h_scale[0]);
		data.Add(h_ams);

		// fit AMS and ACE data at the same time
		vector<double> rigmin, rigmax, chi2norm(2);
		ROOT::Fit::Fitter fitter;
		//HistTools::PrintFunction(fit);
		FitTools::SetCommonFitterOptions(fitter);
		FitTools::FitCombinedData(data, fit, "I", rigmin, rigmax, chi2norm, fitter, 3); 

		// data/fit
		TH1D *h_fitres[2];
 		
		for (UShort_t i = 0; i < 2; ++i) 
		{
  			TH1 *hist = HistTools::ToHist(data[i]);
   			h_fitres[i] = (TH1D *)HistTools::GetResiduals(hist, fit, "_fitres", false, true, true, 5, 1, 0.68, &fitter);
			HistTools::CopyStyle(hist, h_fitres[i]);

   			UShort_t ndf  = hist->GetNbinsX();
   			Double_t chi2 = chi2norm[i]*ndf;
   			printf("old %d nodes ### %-34s   chi2/ndf=%6.2f/%-2u   chi2norm=%5.2f   prob.=%5.2f%%\n", nnodes, hist->GetTitle(), chi2, ndf, chi2norm[i], TMath::Prob(chi2, ndf)*1e2);
		
		} 

		HistTools::SetMarkerStyle(h_fitres[0], kBlue, kFullCircle, 0.9);  
		HistTools::SetMarkerStyle(h_fitres[1], kBlue, kFullCircle, 0.9);  

		g_chi2_rescale->SetPoint(0, 1, fit->GetChisquare()); 
		g_residual_rescale->SetPoint(0, 1, h_fitres[1]->GetBinContent(1)); 

		int jmax = 7; // rescaling index 

		if (temp == 0) jmax = 7;
		if (temp == 1) jmax = 2;
		if (temp == 2) jmax = 7;
		if (temp == 3) jmax = 6; 

		for (int j=1; j<jmax; j++){ 

			// data/fit
			TH1D *h_fitres2[2];

			h_fitres2[1] = h_fitres[1]; 

			h_scale[j] = (TH1 *) h_scale[j-1]->Clone("h_scale");			
			HistTools::SetStyle(h_scale[j], HistTools::GetColorPalette(i, 4), kFullCircle, 1.5, 1, 1); 

			h_scale[j]->Scale(1.0 + h_fitres2[1]->GetBinContent(1)); 

			// Fit again !! 

			// create spline
			Spline *sp2 = new Spline("sp2", nnodes, Spline::LogLog | Spline::PowerLaw, xnodes, ynodes);
			sp2->SetSpectralIndices(3, -2.8, 0, 4, -3.5, -1); // set initial values and limits for the spectral index at first and last node
			TF1 *fit2 = sp2->GetTF1Pointer();

			// sp2->Print();  

			// create array of data to be fitted together
			TObjArray data2; 
			data2.Add(h_scale[j]);
			data2.Add(h_ams);

			// fit AMS and ACE data at the same time
			FitTools::FitCombinedData(data2, fit2, "I", rigmin, rigmax, chi2norm, fitter, 3); 

			for (UShort_t i = 0; i < 2; ++i)
			{
  				TH1 *hist = HistTools::ToHist(data2[i]);
   				h_fitres2[i] = (TH1D *)HistTools::GetResiduals(hist, fit2, "_fitres", false, true, true, 5, 1, 0.68, &fitter);
				HistTools::CopyStyle(hist, h_fitres2[i]);

   				UShort_t ndf  = hist->GetNbinsX();
   				Double_t chi2 = chi2norm[i]*ndf;
   				printf("new %d nodes ### %-34s   chi2/ndf=%6.2f/%-2u   chi2norm=%5.2f   prob.=%5.2f%%\n", nnodes, hist->GetTitle(), chi2, ndf, chi2norm[i], TMath::Prob(chi2, ndf)*1e2);
		
			}

			HistTools::SetMarkerStyle(h_fitres2[0], kBlue, kFullCircle, 0.9);
			HistTools::SetMarkerStyle(h_fitres2[1], kBlue, kFullCircle, 0.9);

			TLegend *legend = new TLegend(0.62,0.8,0.9,0.9); // left, down, right, top
			legend->AddEntry(h_scale[j], Form("Estimated AMS %s Low Energy Flux", AMS_Element2_Cap[i]), "p");
			legend->AddEntry(h_ams, Form("AMS Integrated %s Flux", AMS_Element2_Cap[i]), "p");
			legend->AddEntry(fit2, Form("Estimated+Actual AMS Combined %s Flux Reconstruction", AMS_Element2_Cap[i]), "l"); 

			// PRINT_HIST(h_scale); 

			c1->cd(1); 
			gPad->SetLogx();
			gPad->SetLogy();
			gPad->SetGrid();	

			ha->SetTitle(Form("Estimated AMS %s Flux vs. Rigidity ( %s Template, chi2/ndf=%6.2f/%d )", AMS_Element2_Cap[i], ACE_Element[temp], fit2->GetChisquare(), fit2->GetNDF()));	

			ha->Draw("E1X0");
			fit2->Draw("SAME"); 
			h_ams->Draw("E1X0 SAME"); 
			h_scale[j]->Draw("E1X0 SAME");
			legend->Draw("SAME"); 

			c1->cd(2); 
			gPad->SetLogx();
			gPad->SetGrid();
	
			TH1 *h_resaxis = HistTools::CreateAxis("h_resaxis", " ; ; Data / Fit - 1", 0.1, 2500., 7, -0.6, 0.6, false);
			HistTools::CopyStyle(ha, h_resaxis); 
			h_resaxis->SetXTitle(Unit::GetEnergyLabel("GV")); 

			h_resaxis->Draw("E1X0"); 
			h_fitres2[0]->Draw("E1X0 SAME");
			h_fitres2[1]->Draw("E1X0 SAME"); 

			// PRINT_HIST(h_fitres[1]) 
			// PRINT_HIST(h_fitres2[1]) 
	
			g_chi2_rescale->SetPoint(j, 1.0+h_fitres2[1]->GetBinContent(1), fit2->GetChisquare()); 
			g_residual_rescale->SetPoint(j, 1.0+h_fitres2[1]->GetBinContent(1), h_fitres2[1]->GetBinContent(1)); 

			printf(" %s rescaled by %d times fit residual = %10.4f \n", AMS_Element2_Cap[i], j, h_fitres2[1]->GetBinContent(1));  
			printf(" %s rescaled by %d times global fit chi2/ndf = %10.4f/%d \n", AMS_Element2_Cap[i], j, fit2->GetChisquare(), fit2->GetNDF());  

			TFile file0(Form("data/ACE/extend2/fit_%s_temp_%s_%dnodes.root", ACE_Element[temp], AMS_Element2[i], nnodes), "RECREATE"); 
		
			fit2->Write("fsp_be");
			h_fitres2[0]->Write("h_fitres_ace");
			h_fitres2[1]->Write("h_fitres_ams"); 

			file0.Write();
			file0.Close();	

			c1->Print(Form("./data/ACE/extend2/%s_temp/extend_low_%s_rescale_%dtimes.png", ACE_Element[temp], AMS_Element2_Cap[i], j));			

		}

		// PRINT_GRAPH(g_chi2_rescale)
		// PRINT_GRAPH(g_residual_rescale)
		
		c2->cd(1);
		gPad->SetGrid(); 

		HistTools::SetStyle(g_chi2_rescale, kRed, kFullCircle, 1.5, 1, 1); 
		g_chi2_rescale->SetTitle(Form("%s Global Chi2 vs. Rescaling Factor; Rescaling Factor; Chi2", AMS_Element2_Cap[i]));  
		g_chi2_rescale->Draw("AP"); 

		c2->cd(2);
		gPad->SetGrid();
		HistTools::SetStyle(g_residual_rescale, kRed, kFullCircle, 1.5, 1, 1); 
		g_residual_rescale->SetTitle(Form("%s Fitting Residual of 1st Bin of AMS Range vs. Rescaling Factor; Rescaling Factor; Residual", AMS_Element2_Cap[i]));
		g_residual_rescale->Draw("AP");  

		c2->Print(Form("./data/ACE/extend2/%s_temp/extend_chi2_residual_vs_rescaling_factor_%s.png", ACE_Element[temp], AMS_Element2_Cap[i])); 

	   } else {

		continue;

	   }

	}

}

// plot the canvas for all contributions 
void ace_contribution(){

	Debug::Enable(Debug::ALL);

	TH1D *h_tot_flux[4];
	TH1D *h_tot_mw_flux[4];
	TH1D *h_tot_C_ratio[4];
	TH1D *h_tot_mw_C_ratio[4];

	TH1D *h_contribute[4][28]; 
	TH1D *h_contribute_mw[4][28]; 
	TH1D *h_contribute_ratio[4][28]; // ratio wrt C template 
	TH1D *h_contribute_mw_ratio[4][28]; 

	// plot total 
	TCanvas *c0 = new TCanvas("c0", "Total Contribution", 1800, 900); 
	c0->Divide(2, 2); 	 

	TLegend *legend1 = new TLegend(0.72,0.8,0.9,0.9); 
	TLegend *legend2 = new TLegend(0.72,0.8,0.9,0.9);  
	TLegend *legend3 = new TLegend(0.72,0.8,0.9,0.9); 
	TLegend *legend4 = new TLegend(0.72,0.8,0.9,0.9); 

	for(int i=0; i<4; i++){ 

		TFile *file = new TFile(Form("data/ACE/contribute/%s_temp/h_contribute.root", ACE_Element[i]));  

		h_tot_flux[i] = (TH1D *) file->Get("h_tot_cr_flux")->Clone(Form("h_tot_flux_%s_temp", ACE_Element[i]));
		HistTools::SetStyle(h_tot_flux[i], HistTools::GetColorPalette(i, 8) , 21, 1.5, 1, 1); 

		legend1->AddEntry(h_tot_flux[i], Form("%s template", ACE_Element[i]));

		h_tot_mw_flux[i] = (TH1D *) file->Get("h_tot_cr_flux_mw")->Clone(Form("h_tot_mw_flux_%s_temp", ACE_Element[i]));
		HistTools::SetStyle(h_tot_mw_flux[i], HistTools::GetColorPalette(i, 8) , 21, 1.5, 1, 1); 

		legend2->AddEntry(h_tot_mw_flux[i], Form("%s template", ACE_Element[i]));  

		c0->cd(1);
		gPad->SetLogx();
		gPad->SetLogy();
		gPad->SetGrid();
		gPad->SetMargin(0.11, 0.08, 0.08, 0.08);

		h_tot_mw_flux[i]->SetTitle("Normal Total Cosmic Ray Flux");

		if (i==0) h_tot_flux[i]->Draw("HIST");
		else h_tot_flux[i]->Draw("HIST SAME");  

		c0->cd(2);
		gPad->SetLogx();
		gPad->SetLogy();
		gPad->SetGrid();
		gPad->SetMargin(0.11, 0.08, 0.08, 0.08);

		h_tot_mw_flux[i]->SetTitle("Mass-Weighted Total Cosmic Ray Flux");

		if (i==0) h_tot_mw_flux[i]->Draw("HIST");
		else h_tot_mw_flux[i]->Draw("HIST SAME");  

	} 

	for(int i=0; i<4; i++){

		c0->cd(3);
		gPad->SetLogx();
		//gPad->SetLogy();
		gPad->SetGrid();
		gPad->SetMargin(0.11, 0.08, 0.08, 0.08);

		h_tot_C_ratio[i] = (TH1D *) h_tot_flux[i]->Clone(Form("h_tot_C_ratio_%s_temp", ACE_Element[i])); 
		h_tot_C_ratio[i]->SetTitle("Total Flux Templates wrt C");
		h_tot_C_ratio[i]->SetYTitle("Ratio");
		//h_tot_C_ratio[i]->Print("range"); 
		//h_tot_flux[1]->Print("range"); 

		h_tot_C_ratio[i]->Divide(h_tot_flux[1]); 
	
		if (i!=1) legend3->AddEntry(h_tot_C_ratio[i], Form("%s template", ACE_Element[i]));

		if (i==0) {

			h_tot_C_ratio[i]->GetYaxis()->SetRangeUser(0.975, 1.2); 
			h_tot_C_ratio[i]->Draw("HIST"); 

		} else if (i>1) h_tot_C_ratio[i]->Draw("HIST SAME");  

		c0->cd(4);
		gPad->SetLogx();
		//gPad->SetLogy();
		gPad->SetGrid();
		gPad->SetMargin(0.11, 0.08, 0.08, 0.08);

		h_tot_mw_C_ratio[i] = (TH1D *) h_tot_mw_flux[i]->Clone(Form("h_tot_mw_C_ratio_%s_temp", ACE_Element[i]));
		h_tot_mw_C_ratio[i]->SetTitle("Mass-weighted Total Flux Templates wrt C");
		h_tot_mw_C_ratio[i]->SetYTitle("Ratio");   

		//h_tot_mw_C_ratio[i]->Print("range"); 
		//h_tot_mw_flux[1]->Print("range"); 

		h_tot_mw_C_ratio[i]->Divide(h_tot_mw_flux[1]); 

		if (i!=1) legend4->AddEntry(h_tot_mw_C_ratio[i], Form("%s template", ACE_Element[i]));

		if (i==0) h_tot_mw_C_ratio[i]->Draw("HIST");
		else if (i>1) h_tot_mw_C_ratio[i]->Draw("HIST SAME"); 

	}

	c0->cd(1);
	legend1->Draw("SAME");

	c0->cd(2);
	legend2->Draw("SAME"); 

	c0->cd(3);
	legend3->Draw("SAME");
	
	c0->cd(4); 
	legend4->Draw("SAME"); 

	c0->Print("./data/ACE/contribute/contribution_analysis_BCNO.png");  

	for(int i=0; i<4; i++){
	
		TFile *file = new TFile(Form("data/ACE/contribute/%s_temp/h_contribute.root", ACE_Element[i])); 

		for (int j=0; j<n_total; j++){

			h_contribute[i][j] = (TH1D *) file->Get(Form("h_contribute_ratio_%s", Element[j]))->Clone(Form("h_contribute_%s_temp_%s", ACE_Element[i], Element[j])); 
			h_contribute[i][j]->SetTitle(Form("%s Contribution", Element[j])); 
			h_contribute_mw[i][j] = (TH1D *) file->Get(Form("h_contribute_mw_ratio_%s", Element[j]))->Clone(Form("h_contribute_mw_%s_temp_%s", ACE_Element[i], Element[j]));
			h_contribute_mw[i][j]->SetTitle(Form("%s MW Contribution", Element[j]));  
			HistTools::SetStyle(h_contribute[i][j], HistTools::GetColorPalette(i, 4) , kFullCircle, 1.1, 1, 1); 
			h_contribute[i][j]->SetFillStyle(4050);
			HistTools::SetStyle(h_contribute_mw[i][j], HistTools::GetColorPalette(i, 4) , kFullCircle, 1.1, 1, 1);
			h_contribute_mw[i][j]->SetFillStyle(4050);

		}

	}

	for(int i=0; i<4; i++){
	
		TFile *file = new TFile(Form("data/ACE/contribute/%s_temp/h_contribute.root", ACE_Element[i])); 

		for (int j=0; j<n_total; j++){

			if (i!=1){
				h_contribute_ratio[i][j] = (TH1D*) h_contribute[i][j]->Clone(Form("h_contribute_raio_%s_temp_%s", ACE_Element[i], Element[j])); 
				h_contribute_ratio[i][j]->SetTitle(Form("Ratio of %s Contribution wrt C Template; Rigidity [GV]; Ratio", Element[j]));  
				h_contribute_ratio[i][j]->Divide(h_contribute[1][j]); 
				h_contribute_mw_ratio[i][j] = (TH1D*) h_contribute_mw[i][j]->Clone(Form("h_contribute_mw_raio_%s_temp_%s", ACE_Element[i], Element[j]));
				h_contribute_mw_ratio[i][j]->SetTitle(Form("Ratio of %s MW Contribution wrt C Template; Rigidity [GV]; Ratio", Element[j]));
				h_contribute_mw_ratio[i][j]->Divide(h_contribute_mw[1][j]);  
				HistTools::SetStyle(h_contribute_ratio[i][j], HistTools::GetColorPalette(i, 4) , kFullCircle, 1.1, 1, 1); 
				h_contribute_ratio[i][j]->SetFillStyle(4050);
				HistTools::SetStyle(h_contribute_mw_ratio[i][j], HistTools::GetColorPalette(i, 4) , kFullCircle, 1.1, 1, 1); 
				h_contribute_mw_ratio[i][j]->SetFillStyle(4050);
			}

		}

	}

	for (int j=0; j<n_total; j++){ 

		TCanvas *c1 = new TCanvas("c1", "Individual Contribution", 1800, 900); 
		c1->Divide(2, 2);

		TLegend *legend1 = new TLegend(0.72,0.8,0.9,0.9); 
		TLegend *legend2 = new TLegend(0.72,0.8,0.9,0.9);  
		TLegend *legend3 = new TLegend(0.72,0.8,0.9,0.9); 
		TLegend *legend4 = new TLegend(0.72,0.8,0.9,0.9); 

		c1->cd(1);
		gPad->SetLogx();
		gPad->SetLogy();
		gPad->SetGrid(); 

		for(int i=0; i<4; i++){
		
			if (i==0) h_contribute[i][j]->Draw("HIST"); 
			else h_contribute[i][j]->Draw("HIST SAME"); 
			legend1->AddEntry(h_contribute[i][j], Form("%s template", ACE_Element[i]));

		}
		legend1->Draw("SAME");

		c1->cd(2);
		gPad->SetLogx();
		gPad->SetLogy();
		gPad->SetGrid(); 

		for(int i=0; i<4; i++){
		
			if (i==0) h_contribute_mw[i][j]->Draw("HIST");
			else h_contribute_mw[i][j]->Draw("HIST SAME"); 
			legend2->AddEntry(h_contribute[i][j], Form("%s template", ACE_Element[i]));
		
		}
		legend2->Draw("SAME");

		c1->cd(3);
		gPad->SetLogx();
		//gPad->SetLogy();
		gPad->SetGrid();

		for(int i=0; i<4; i++){
		
			if (i==0) h_contribute_ratio[i][j]->Draw("HIST"); 
			if (i!=1) {
				h_contribute_ratio[i][j]->Draw("HIST SAME"); 
				legend3->AddEntry(h_contribute_ratio[i][j], Form("%s template", ACE_Element[i]));
			}
		}
		legend3->Draw("SAME");

		c1->cd(4);
		gPad->SetLogx();
		//gPad->SetLogy();
		gPad->SetGrid();

		for(int i=0; i<4; i++){
		
			if (i==0) h_contribute_mw_ratio[i][j]->Draw("HIST"); 
			if (i!=1) {
				h_contribute_mw_ratio[i][j]->Draw("HIST SAME"); 
				legend4->AddEntry(h_contribute_mw_ratio[i][j], Form("%s template", ACE_Element[i]));
			}
		}
		legend4->Draw("SAME");

		c1->Print(Form("./data/ACE/contribute/contribution_analysis_BCNO_%s.png", Element[j]));  

	}

	// create the element contribution vs. rigidity table for BCNO templates

	const int n_bins = 11; 
	Double_t R_bins[n_bins+1] = {0.8, 1, 2, 5, 10, 20, 50, 100, 300, 1000, 2000};

	for (int i=0; i<4; i++){

		cout << ACE_Element[i] << " Template" << endl; 

		for (int j=0; j<n_total; j++){

			// cout << "Normal " << Element[j] << " Relative Contribution "; 
			
			for (int ibin=0; ibin<n_bins; ibin++){

				int bin = h_contribute[i][j]->FindBin(R_bins[ibin]); 
				double rel_con = h_contribute[i][j]->GetBinContent(bin); // relative contribution
				
				printf("%10.8f ", rel_con); 
	
			}

			cout << endl; 

		}

		cout << " " << endl; 

		cout << ACE_Element[i] << " Template" << endl; 

		for (int j=0; j<n_total; j++){

			// cout << "Mass-weighted " << Element[j] << " Relative Contribution "; 
			
			for (int ibin=0; ibin<n_bins; ibin++){

				int bin = h_contribute_mw[i][j]->FindBin(R_bins[ibin]); 
				double rel_con = h_contribute_mw[i][j]->GetBinContent(bin); // mass-weighted relative contribution
				
				printf("%10.8f ", rel_con); 
	
			}

			cout << endl;
		}

		cout << " " << endl; 
	} 



}

// contribution v2, the truth table 
void ace_contribution2(){

	Debug::Enable(Debug::ALL);

	TH1D *h_tot_flux;
	TH1D *h_tot_mw_flux;
	TH1D *h_tot_C_ratio;
	TH1D *h_tot_mw_C_ratio;

	TH1D *h_contribute[28]; 
	TH1D *h_contribute_mw[28]; 
	TH1D *h_contribute_ratio[28]; // ratio wrt C template 
	TH1D *h_contribute_mw_ratio[28]; 

	TFile *file = new TFile("data/ACE/contribute2/h_contribute.root");  

	h_tot_flux = (TH1D *) file->Get("h_tot_cr_flux")->Clone("h_tot_flux");
	HistTools::SetStyle(h_tot_flux, kBlue, 21, 1.5, 1, 1); 

	h_tot_mw_flux = (TH1D *) file->Get("h_tot_cr_flux_mw")->Clone("h_tot_mw_flux");
	HistTools::SetStyle(h_tot_mw_flux, kRed, 21, 1.5, 1, 1); 

	for (int j=0; j<n_total; j++){

		h_contribute[j] = (TH1D *) file->Get(Form("h_contribute_ratio_%s", Element[j]))->Clone(Form("h_contribute_%s", Element[j])); 
		h_contribute[j]->SetTitle(Form("%s Contribution", Element[j])); 
		h_contribute_mw[j] = (TH1D *) file->Get(Form("h_contribute_mw_ratio_%s", Element[j]))->Clone(Form("h_contribute_mw_%s", Element[j]));
		h_contribute_mw[j]->SetTitle(Form("%s MW Contribution", Element[j]));  
		HistTools::SetStyle(h_contribute[j], kBlue , kFullCircle, 1.1, 1, 1); 
		h_contribute[j]->SetFillStyle(4050);
		HistTools::SetStyle(h_contribute_mw[j], kRed, kFullCircle, 1.1, 1, 1);
		h_contribute_mw[j]->SetFillStyle(4050);

	}

	// create the element contribution vs. rigidity table aka the truth table 

	const int n_bins = 11; 
	Double_t R_bins[n_bins+1] = {0.8, 1, 2, 5, 10, 20, 50, 100, 300, 1000, 2000};

	for (int j=0; j<n_total; j++){

		// cout << "Normal " << Element[j] << " Relative Contribution "; 
			
		for (int ibin=0; ibin<n_bins; ibin++){

			int bin = h_contribute[j]->FindBin(R_bins[ibin]); 
			double rel_con = h_contribute[j]->GetBinContent(bin); // relative contribution
				
			printf("%10.8f ", rel_con); 
	
		}

		cout << endl; 

	}

	cout << " " << endl; 

	for (int j=0; j<n_total; j++){

		// cout << "Mass-weighted " << Element[j] << " Relative Contribution "; 
			
		for (int ibin=0; ibin<n_bins; ibin++){

			int bin = h_contribute_mw[j]->FindBin(R_bins[ibin]); 
			double rel_con = h_contribute_mw[j]->GetBinContent(bin); // mass-weighted relative contribution
				
			printf("%10.8f ", rel_con); 
	
		}

		cout << endl;
	}


}

// compute the contribution to the total cosmic ray flux with real AMS data + extended ACE data 
// option A, total; option B, mass-weighted 
TH1 *ace_contribute(int temp, const char *option){

	Debug::Enable(Debug::ALL); 

	Experiments::DataPath = "data";

	gSystem->mkdir(Form("data/ACE/contribute/%s_temp", ACE_Element[temp]), true);

	int data_value[n_ele] = {0, 18, 23, 27, 31, 20, 43, 22}; 

	const char *AMS_Element2[8] = { "p", "he", "li", "be", "b", "c", "n", "o" }; 
	const char *AMS_Element2_Cap[8] = { "Proton", "He", "Li", "Be", "B", "C", "N", "O" }; 

	const UInt_t FirstACEBR = 2240;
   	vector<UInt_t> BRs;
  	// we stop at BR 2493, which ends on 2016/05/23, just 3 days before the end of the data taking period for AMS nuclei
   	for (UInt_t br=2426; br<=2493; ++br) { 
	   	if (br != 2472 && br != 2473) BRs.push_back(br-FirstACEBR); 
	}

	TH1 *h_ams_C = Experiments::GetMeasurementHistogram(Experiments::AMS02, data_value[1], 0); // load AMS Template data in order to determine Rmax

	TH1 *h_ene_C = HistTools::GraphToHist(get_ace_average_graph( ACE_Element[1] , &BRs[0], BRs.size() ), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
	TH1 *h_ace_C = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene_C, ACE_Isotope[1], "MeV/n cm", "GV m", "_rig"); // load averaged ACE data for the same element in rigidity

	UShort_t namsbins = h_ams_C->GetNbinsX(); 
	UShort_t nacebins = h_ace_C->GetNbinsX(); 
	
	double R1 = h_ace_C->GetBinLowEdge(1);
	double R2 = h_ams_C->GetBinLowEdge(namsbins+1); 

	//printf("R1 = %f, R2 = %f \n", R1, R2); 

	TCanvas *c1 = new TCanvas("c1", "Extend Fit", 1800, 900); 
	c1->Divide(1, 2); 

	gStyle->SetOptStat(0); 

	TF1 *f_fit[n_ele+4]; 

	// Sum up the Total Cosmic Ray Flux 
	// p, He, Li, Be
	for (int i=0;i<4;++i){
		
		f_fit[i] = new TF1(); 

		TH1 *ha = HistTools::CreateAxis("ha", "haxis1", 0.1, 2500., 7, 1e-10, 1e3, false); 
		ha->SetXTitle(Unit::GetEnergyLabel("GV"));
  		ha->SetYTitle(Unit::GetDifferentialFluxLabel("GV m")); 
		ha->SetTitle(Form("Integrated AMS %s Flux vs. Rigidity", AMS_Element2_Cap[i]));

		TH1 *h_ams = Experiments::GetMeasurementHistogram(Experiments::AMS02, data_value[temp+4], 0); // AMS Template
		HistTools::SetStyle(h_ams, kBlue, kFullCircle, 1.5, 1, 1); 

		int nnodes = 9;  
		int nnodes_ams = 6; 
		
		TFile file0(Form("data/ACE/extend2/fit_%s_temp_%s_%dnodes.root", ACE_Element[temp], AMS_Element2[i], nnodes)); 

		TFile file1(Form("data/amsfit/fit_result_node%d.root", nnodes_ams));  

		TH1 *h_ams_i = (TH1 *) file1.Get(Form("h_%s", AMS_Element2[i])); // AMS p, He, Li, Be 
		HistTools::SetStyle(h_ams_i, kBlue, kFullCircle, 1.5, 1, 1);  

		Spline *sp_ams = new Spline("sp_ams", nnodes, Spline::LogLog | Spline::PowerLaw);
		f_fit[i] = sp_ams->GetTF1Pointer(); 
		TF1 *fit_ams = (TF1*) file0.Get(Form("fsp_%s", AMS_Element2[i])); // load the AMS fluxes, the fsp in the root file is actually fit 

		HistTools::CopyParameters(fit_ams, f_fit[i]);
		double x1, x2;
		fit_ams->GetRange(x1,x2);
		f_fit[i]->SetRange(x1,x2); 

		TH1 *h_fitres = (TH1 *) file1.Get(Form("h_%s_fitres", AMS_Element2[i])); 
		TH1 *h_fiterr = (TH1 *) file1.Get(Form("h_%s_fiterr", AMS_Element2[i]));  
	
		TLegend *legend = new TLegend(0.62,0.8,0.9,0.9); // left, down, right, top
		legend->AddEntry(h_ams_i, Form("AMS Integrated %s Flux", AMS_Element2_Cap[i]), "p");
		//legend->AddEntry(h_ams, Form("AMS Integrated %s Flux", ACE_Element[1]), "p"); 
		legend->AddEntry(f_fit[i], Form("AMS Integrated %s Flux Fit", AMS_Element2_Cap[i]), "l"); 
		
		c1->cd(1);
		gPad->SetLogx();
		gPad->SetLogy();
		gPad->SetGrid();

		ha->Draw("E1X0"); 
		h_ams_i->Draw("E1X0 SAME"); 
		// h_ams->Draw("E1X0 SAME"); 
		f_fit[i]->Draw("LSAME");
		legend->Draw("SAME");

		c1->cd(2);
		gPad->SetLogx();
		gPad->SetGrid(); 

		TH1 *h_resaxis = HistTools::CreateAxis("h_resaxis", " ; ; Data / Fit - 1", 0.1, 2500., 7, 0.9, 1.1, false);
		h_resaxis->SetXTitle(Unit::GetEnergyLabel("GV")); 

		h_resaxis->Draw("E1X0");
		h_fiterr->Draw("E3 SAME");
		h_fitres->Draw("E1X0 SAME"); 		 	

		c1->Print(Form("./data/ACE/contribute/%s_temp/extend_%d_%s.png", ACE_Element[temp], i, AMS_Element2_Cap[i]));  

	} 

	// B, C, N, O 
	for (int i=0;i<4;++i){ 

		int j = i+4; 		

		int nnodes = 9;
	
		f_fit[j] = new TF1(); 

		TH1 *h_ams_i[4];

		TH1 *ha = HistTools::CreateAxis("ha", "haxis1", 0.1, 2500., 7, 1e-10, 1e3, false);
		ha->SetXTitle(Unit::GetEnergyLabel("GV"));
  		ha->SetYTitle(Unit::GetDifferentialFluxLabel("GV m")); 
		ha->SetTitle(Form("Averaged ACE and Integrated AMS %s Flux vs. Rigidity", ACE_Element[i]));

		h_ams_i[i] = Experiments::GetMeasurementHistogram(Experiments::AMS02, data_value[j], 0); 
		HistTools::SetStyle(h_ams_i[i], kBlue , kFullCircle, 1.5, 1, 1);  

		TH1 *h_ene = HistTools::GraphToHist(get_ace_average_graph( ACE_Element[i] , &BRs[0], BRs.size() ), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
		TH1 *h_ace = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene, ACE_Isotope[i], "MeV/n cm", "GV m", "_rig"); // load averaged ACE data for the same element in rigidity 
		HistTools::SetStyle(h_ace, HistTools::GetColorPalette(i, n_total) , kFullCircle, 1.5, 1, 1); 

		TFile file2(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[i], nnodes)); // load ACE combined fit 

		Spline *sp_comb = new Spline("sp_comb", nnodes, Spline::LogLog | Spline::PowerLaw); 
		f_fit[j] = sp_comb->GetTF1Pointer();  // real function 
		TF1 *fit_comb = (TF1*) file2.Get("fit_both"); 

		HistTools::CopyParameters(fit_comb, f_fit[j]); 
		double x1, x2; 
		fit_comb->GetRange(x1,x2); 
		f_fit[j]->SetRange(x1,x2); 

		TH1 *h_fitres[2];

		h_fitres[0] = (TH1 *) file2.Get("h_fitres_ace");
		HistTools::CopyStyle(h_ace, h_fitres[0]); 
		h_fitres[1] = (TH1 *) file2.Get("h_fitres_ams"); 
		HistTools::CopyStyle(h_ams_i[i], h_fitres[1]);  

		TLegend *legend = new TLegend(0.1,0.8,0.28,0.9); // xlow, ylow, xhigh, yhigh
		legend->AddEntry(h_ace, Form("ACE %s Flux", ACE_Element[i]), "p");
		legend->AddEntry(h_ams_i[i], Form("AMS Integrated %s Flux", AMS_Element2_Cap[j]), "p");
		legend->AddEntry(f_fit[j], Form("ACE %s Flux Reconstruction", ACE_Element[i]), "l"); 

		c1->cd(1);
		gPad->SetLogx();
		gPad->SetLogy();
		gPad->SetGrid();
	
		ha->Draw("E1X0");
		h_ace->Draw("E1X0 SAME");  
		h_ams_i[i]->Draw("E1X0 SAME"); 
		f_fit[j]->Draw("LSAME");
		legend->Draw("SAME");

		c1->cd(2);
		gPad->SetLogx();
		gPad->SetGrid();

		TH1 *h_resaxis = HistTools::CreateAxis("h_resaxis", " ; ; Data / Fit - 1", 0.1, 2500., 7, -0.5, 0.5, false);
		h_resaxis->SetXTitle(Unit::GetEnergyLabel("GV")); 

		h_resaxis->Draw("E1X0");
		h_fitres[0]->Draw("E1X0 SAME"); 
		h_fitres[1]->Draw("E1X0 SAME");

		c1->Print(Form("./data/ACE/contribute/%s_temp/extend_%d_%s.png", ACE_Element[temp], j, ACE_Element[i])); 
	}

	// F, Ne, Na ... Va, Ni
	for (int i=4;i<24;++i){ 

		int j = i+4; 		

		int nnodes = 7;
	
		f_fit[j] = new TF1(); 

		TH1 *h_ams_i[4];

		TH1 *ha = HistTools::CreateAxis("ha", "haxis1", 0.1, 2500., 7, 1e-10, 1e3, false); 
		ha->SetXTitle(Unit::GetEnergyLabel("GV"));
  		ha->SetYTitle(Unit::GetDifferentialFluxLabel("GV m")); 
		ha->SetTitle(Form("Averaged ACE and Integrated AMS %s Flux vs. Rigidity", ACE_Element[i])); 
		 
		TH1 *h_ams_C = Experiments::GetMeasurementHistogram(Experiments::AMS02, data_value[temp+4], 0); 
		HistTools::SetStyle(h_ams_C, kBlue, kFullCircle, 1.5, 1, 1); 

		TH1 *h_ene = HistTools::GraphToHist(get_ace_average_graph( ACE_Element[i] , &BRs[0], BRs.size() ), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
		TH1 *h_ace = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene, ACE_Isotope[i], "MeV/n cm", "GV m", "_rig"); // load averaged ACE data for the same element in rigidity 
		HistTools::SetStyle(h_ace, HistTools::GetColorPalette(i, n_total) , kFullCircle, 1.5, 1, 1); 

		TFile file2(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[temp], nnodes)); // load ACE Template combined fit 

		Spline *sp_comb_C = new Spline("sp_comb_C", nnodes, Spline::LogLog | Spline::PowerLaw); 
		f_fit[j] = sp_comb_C->GetTF1Pointer();  // real function 
		TF1 *fit_comb_C = (TF1*) file2.Get("fit_both");
	
		HistTools::CopyParameters(fit_comb_C, f_fit[j]); 
		double x1_C, x2_C;
		fit_comb_C->GetRange(x1_C,x2_C);
		f_fit[j]->SetRange(x1_C,x2_C);
		
		TH1 *h_ratio = (TH1 *) h_ace->Clone("h_ratio");
	
		h_ratio->Divide(f_fit[j]); // WARNING! Doing this will change the # of entires of the histogram. 
		HistTools::SetStyle(h_ratio, kRed, kFullCircle, 0.9, 1, 1);
	
		double ratio_sum=0; // compute average of h_ratio manually  
		for(int k=0;k<14;k++){
			ratio_sum += h_ratio->GetBinContent(k); 
		}
		double ratio_ave = ratio_sum/7; 
	
		f_fit[j] = HistTools::CombineTF1Const(f_fit[j], ratio_ave, HistTools::MultiplyConst, "f_fit", R1, R2); 

		TH1 *h_fitres = HistTools::GetResiduals(h_ace, f_fit[j], "_fitratio", false, true, true, 4, 1); 
		HistTools::CopyStyle(h_ace, h_fitres);
	
		TLegend *legend = new TLegend(0.1,0.8,0.28,0.9); // left, down, right, top 
		legend->AddEntry(h_ace, Form("ACE %s Flux", ACE_Element[i]), "p");
		//legend->AddEntry(h_ams_C, Form("AMS Integrated %s Flux", ACE_Element[1]), "p");
		legend->AddEntry(f_fit[j], Form("ACE %s Flux Reconstruction", ACE_Element[i]), "l"); 
			
		c1->cd(1);
		gPad->SetLogx();
		gPad->SetLogy();
		gPad->SetGrid(); 

		ha->Draw("E1X0");
		h_ace->Draw("E1X0 SAME");  
		//h_ams_C->Draw("E1X0 SAME"); 
		f_fit[j]->Draw("LSAME");
		legend->Draw("SAME");

		c1->cd(2); 
		gPad->SetLogx();
		gPad->SetGrid();

		TH1 *h_resaxis = HistTools::CreateAxis("h_resaxis", " ; ; Data / Fit - 1", 0.1, 2500., 7, 0.2, 1.6, false);
		h_resaxis->SetXTitle(Unit::GetEnergyLabel("GV")); 

		h_resaxis->Draw("E1X0");
		//h_fiterr->Draw("E3 SAME");
		
		// PRINT_HIST(h_fitres)
		h_fitres->Draw("E1X0 SAME"); 

		// low rigidity outside temp folders, high rigidity inside temp folders
		TFile file0(Form("data/ACE/extend2/%s_temp/fit_%s_%dnodes.root", ACE_Element[temp], ACE_Element[i], nnodes), "RECREATE"); 
		
		f_fit[j]->Write(Form("fit_%s", Element[j]));

		file0.Write();
		file0.Close();

		c1->Print(Form("./data/ACE/contribute/%s_temp/extend_%d_%s.png", ACE_Element[temp], j, ACE_Element[i])); 
	}

	Int_t nbins = 200;

	Double_t *bins = HistTools::BuildLogBins(R1, R2, nbins); // Rmin is the minimum rigidity of the ACE fluxes, Rmax the maximum rigidity of the AMS fluxes
	TH1D *h_tot_cr_flux = new TH1D("h_tot_cr_flux", "Total Cosmic Ray Flux;Rigidity [GV];Flux [1/(m^{2} sr s GV)]", nbins, bins);
	TH1D *h_tot_cr_flux_mw = new TH1D("h_tot_cr_flux_mw", "Mass-weighted Total Cosmic Ray Flux;Rigidity [GV];Flux [1/(m^{2} sr s GV)]", nbins, bins); // Contributed by each flux multiplied by their mass number 
	HistTools::SetStyle(h_tot_cr_flux, kRed, kFullCircle, 0.9, 1, 1);
	HistTools::SetStyle(h_tot_cr_flux_mw, kBlue, kFullCircle, 0.9, 1, 1);
	
	for (int i = 0; i < n_total; ++i) 
	{
	
   		for (int bin = 1; bin <= nbins; ++bin)
   		{
      			Double_t R = h_tot_cr_flux->GetBinLowEdge(bin);
     			Double_t w = h_tot_cr_flux->GetBinWidth(bin);
      			Double_t flux = f_fit[i]->Integral(R, R+w)/w; 
      			h_tot_cr_flux->AddBinContent(bin, flux);
			h_tot_cr_flux_mw->AddBinContent(bin, flux*A[i]); 

			// printf("i=%d, bin=%d, flux=%f \n", i, bin, flux); 
   		}

		//break; 
	}

	TCanvas *c2 = new TCanvas("c2", "Total Cosmic Ray Flux", 1800, 900);
	c2->Divide(1, 2); 
	c2->cd(1); 
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGrid();

	TLegend *legend = new TLegend(0.62,0.8,0.9,0.9); 
	legend->AddEntry(h_tot_cr_flux, "normal");
	legend->AddEntry(h_tot_cr_flux_mw, "mass-weighted"); 

	h_tot_cr_flux->Draw("HIST");
	h_tot_cr_flux_mw->Draw("HIST SAME"); 
	legend->Draw("SAME");

	// PRINT_HIST(h_tot_cr_flux) 

	c2->cd(2); 
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGrid();

	TH1 *h_tot_cr_ratio = (TH1 *) h_tot_cr_flux_mw->Clone("h_tot_cr_ratio");
	h_tot_cr_ratio->Divide(h_tot_cr_flux);

	h_tot_cr_ratio->SetTitle("Total Cosmic Ray Flux / Mass-weighted Total Cosmic Ray Flux"); 
	h_tot_cr_ratio->SetXTitle(Unit::GetEnergyLabel("GV"));
	h_tot_cr_ratio->SetYTitle("Total Flux / Mass-weighted Total Flux");  
	h_tot_cr_ratio->Draw("HIST"); 

	// PRINT_HIST(h_tot_cr_flux_mw) 

	c2->Print(Form("./data/ACE/contribute/%s_temp/total_cr_flux.png", ACE_Element[temp])); 

	TCanvas *c3 = new TCanvas("c3", "Contribution", 1800, 900); 
	c3->Divide(1, 2); 

	TCanvas *c4 = new TCanvas("c4", "Comparison of Relative Contributions", 1800, 900);
	c4->Divide(1, 2);  

	THStack *hs = new THStack("hs", "Comparison of Relative Contributions");
	THStack *hs_mw = new THStack("hs_mw", "Comparison of Mass-Weighted Relative Contributions");

	// Now compute the ratio of h_element vs. h_tot_cr_flux 

	TH1 *h_contribute[n_total]; // flux
	TH1 *h_contribute_ratio[n_total]; // flux / total
	TH1 *h_contribute_mw[n_total]; // mass-weighted flux
	TH1 *h_contribute_mw_ratio[n_total]; // mass-weighted flux / total 

	// p, He, Li, Be
	for (int i=0;i<4;++i){

		int nnodes = 6;  

		h_contribute[i] = (TH1 *) h_tot_cr_flux->Clone("h_contribute");  
		h_contribute_mw[i] = (TH1 *) h_tot_cr_flux->Clone("h_contribute_mw"); 
		HistTools::SetStyle(h_contribute_mw[i], kBlue, kFullCircle, 0.9, 1, 1); 

		for (int k=1;k<=nbins;++k){ 

			Double_t R = h_tot_cr_flux->GetBinLowEdge(k);  
     			Double_t w = h_tot_cr_flux->GetBinWidth(k);  
      			Double_t flux = f_fit[i]->Integral(R, R+w)/w;  

			h_contribute[i]->SetBinContent(k, flux);
			h_contribute_mw[i]->SetBinContent(k, A[i]*flux);
			h_contribute[i]->SetBinError(k, 0); 
			h_contribute_mw[i]->SetBinError(k, 0); 

		} 

		// PRINT_HIST(h_contribute[i])

		c3->cd(1);
		gPad->SetLogx();
		gPad->SetLogy();
		gPad->SetGrid(); 

		// relative contribution 

		h_contribute_ratio[i] = (TH1 *) h_contribute[i]->Clone("h_contribute_ratio"); 
		h_contribute_ratio[i]->SetXTitle(Unit::GetEnergyLabel("GV")); 
   		h_contribute_ratio[i]->SetYTitle("Flux / Total Flux"); 
		h_contribute_ratio[i]->Divide(h_tot_cr_flux);   

		double R_1 = h_tot_cr_flux->GetBinLowEdge(1);
		double R_2 = h_tot_cr_flux->GetBinLowEdge(nbins+1);
		
		double contribute_rel = h_contribute[i]->Integral(R_1, R_2)/h_tot_cr_flux->Integral(R_1, R_2); 	

		printf("%s contribution = %f \n", AMS_Element2_Cap[i], contribute_rel); 

		// mass-weighted relative contribution 

		h_contribute_mw_ratio[i] = (TH1 *) h_contribute_mw[i]->Clone("h_contribute_mw_ratio"); 
		h_contribute_mw_ratio[i]->Divide(h_tot_cr_flux_mw); 
		HistTools::SetStyle(h_contribute_mw_ratio[i], kBlue, kFullCircle, 0.9, 1, 1);
		
		double contribute_rel_mw = h_contribute_mw[i]->Integral(R_1, R_2)/h_tot_cr_flux_mw->Integral(R_1, R_2); 			

		printf("%s mass-weighted contribution = %f \n", AMS_Element2_Cap[i], contribute_rel_mw);

		TLegend *legend2 = new TLegend(0.62,0.8,0.9,0.9); 
		legend2->AddEntry(h_contribute_ratio[i], "normal contribution");
		legend2->AddEntry(h_contribute_mw_ratio[i], "mass-weighted contribution");

		h_contribute_mw_ratio[i]->SetStats(0);
		h_contribute_ratio[i]->SetTitle(Form("%s Contribution = %10.6f", AMS_Element2_Cap[i], contribute_rel));
		h_contribute_mw_ratio[i]->SetTitle(Form("%s MW Contribution = %10.6f", AMS_Element2_Cap[i], contribute_rel_mw));
		h_contribute_mw_ratio[i]->SetXTitle(Unit::GetEnergyLabel("GV")); 
   		h_contribute_mw_ratio[i]->SetYTitle("Flux / Total Flux");
		h_contribute_mw_ratio[i]->Draw("HIST");
		h_contribute_ratio[i]->Draw("HIST SAME");
		legend2->Draw("SAME");

		c3->cd(2);
		gPad->SetLogx();
		gPad->SetLogy();
		gPad->SetGrid();

		TH1 *h_ratio = (TH1 *) h_contribute_mw_ratio[i]->Clone("h_ratio"); 
		h_ratio->Divide(h_contribute_ratio[i]); 

		TLegend *legend1 = new TLegend(0.62,0.8,0.9,0.9); 
		// legend1->AddEntry(h_contribute[i], "normal flux");
		// legend1->AddEntry(h_contribute_mw[i], "mass-weighted flux");

		h_ratio->SetStats(0);
		h_ratio->SetTitle(Form("%s Contribution / Mass-weighted Contribution", AMS_Element2_Cap[i]));
		h_ratio->SetXTitle(Unit::GetEnergyLabel("GV"));
   		h_ratio->SetYTitle("Ratio");
		h_ratio->Draw("HIST");
		// legend1->Draw("SAME"); 

		c3->Print(Form("./data/ACE/contribute/%s_temp/contribution_%d_%s.png", ACE_Element[temp], i, AMS_Element2_Cap[i])); 

	}

	// B, C, N, O
	for (int i=0;i<4;++i){ 

		int j = i+4; 		

		int nnodes = 6; 

		h_contribute[j] = (TH1 *) h_tot_cr_flux->Clone("h_contribute");  
		h_contribute_mw[j] = (TH1 *) h_tot_cr_flux->Clone("h_contribute_mw"); 
		HistTools::SetStyle(h_contribute_mw[j], kBlue, kFullCircle, 0.9, 1, 1); 

		for (int k=1;k<=nbins;++k){ 

			Double_t x = h_tot_cr_flux->GetBinCenter(k), y = h_tot_cr_flux->GetBinContent(k);  

			// printf("k=%d, x=%f, y=%f \n", k, x, y); 

			h_contribute[j]->SetBinContent(k, f_fit[j]->Eval(x));
			h_contribute_mw[j]->SetBinContent(k, A[j]*f_fit[j]->Eval(x));
			h_contribute[j]->SetBinError(k, 0); 
			h_contribute_mw[j]->SetBinError(k, 0); 

		} 

		// PRINT_HIST(h_contribute[j])

		c3->cd(1);
		gPad->SetLogx();
		gPad->SetLogy();
		gPad->SetGrid(); 

		// relative contribution

		h_contribute_ratio[j] = (TH1 *) h_contribute[j]->Clone("h_contribute_ratio");
		h_contribute_ratio[j]->SetXTitle(Unit::GetEnergyLabel("GV")); 
   		h_contribute_ratio[j]->SetYTitle("Flux / Total Flux"); 
		h_contribute_ratio[j]->Divide(h_tot_cr_flux); 

		double R_1 = h_tot_cr_flux->GetBinLowEdge(1);
		double R_2 = h_tot_cr_flux->GetBinLowEdge(nbins+1);

		double contribute_rel = h_contribute[j]->Integral(R_1, R_2)/h_tot_cr_flux->Integral(R_1, R_2); 			

		printf("%s contribution = %f \n", ACE_Element[i], contribute_rel); 

		// mass-weighted ratio

		h_contribute_mw_ratio[j] = (TH1 *) h_contribute_mw[j]->Clone("h_contribute_mw_ratio"); 
		h_contribute_mw_ratio[j]->Divide(h_tot_cr_flux_mw); 
		HistTools::SetStyle(h_contribute_mw_ratio[j], kBlue, kFullCircle, 0.9, 1, 1); 

		double contribute_rel_mw = h_contribute_mw[j]->Integral(R_1, R_2)/h_tot_cr_flux_mw->Integral(R_1, R_2); 			

		printf("%s mass-weighted contribution = %f \n", ACE_Element[i], contribute_rel_mw); 

		TLegend *legend2 = new TLegend(0.62,0.8,0.9,0.9); 
		legend2->AddEntry(h_contribute_ratio[j], "normal contribution");
		legend2->AddEntry(h_contribute_mw_ratio[j], "mass-weighted contribution");

		h_contribute_mw_ratio[j]->SetStats(0);
		h_contribute_ratio[j]->SetTitle(Form("%s Contribution = %10.6f", ACE_Element[i], contribute_rel));
		h_contribute_mw_ratio[j]->SetTitle(Form("%s MW Contribution = %10.6f", ACE_Element[i], contribute_rel_mw));
		h_contribute_mw_ratio[j]->SetXTitle(Unit::GetEnergyLabel("GV")); 
   		h_contribute_mw_ratio[j]->SetYTitle("Flux / Total Flux");
		// h_contribute_mw_ratio[j]->GetYaxis()->SetRangeUser(10e-4, 10e-2);
		h_contribute_mw_ratio[j]->Draw("HIST");
		h_contribute_ratio[j]->Draw("HIST SAME");
		legend2->Draw("SAME"); 

		c3->cd(2);
		gPad->SetLogx();
		gPad->SetLogy();
		gPad->SetGrid();

		TH1 *h_ratio = (TH1 *) h_contribute_mw_ratio[j]->Clone("h_ratio"); 
		h_ratio->Divide(h_contribute_ratio[j]); 

		TLegend *legend1 = new TLegend(0.62,0.8,0.9,0.9); 
		// legend1->AddEntry(h_contribute[j], "normal flux");
		// legend1->AddEntry(h_contribute_mw[j], "mass-weighted flux");

		h_ratio->SetStats(0);
		h_ratio->SetTitle(Form("%s Contribution / Mass-weighted Contribution", ACE_Element[i]));
		h_ratio->SetXTitle(Unit::GetEnergyLabel("GV"));
   		h_ratio->SetYTitle("Ratio");
		h_ratio->Draw("HIST");
		// legend1->Draw("SAME"); 

		c3->Print(Form("./data/ACE/contribute/%s_temp/contribution_%d_%s.png", ACE_Element[temp], j, ACE_Element[i])); 
 
	}

	// F, Ne, Na ... Va, Ni
	for (int i=4;i<24;++i){ 

		int j = i+4; 		

		h_contribute[j] = (TH1 *) h_tot_cr_flux->Clone("h_contribute");  
		h_contribute_mw[j] = (TH1 *) h_tot_cr_flux->Clone("h_contribute_mw"); 
		HistTools::SetStyle(h_contribute_mw[j], kBlue, kFullCircle, 0.9, 1, 1); 

		for (int k=1;k<=nbins;k++){ 

			Double_t x = h_tot_cr_flux->GetBinCenter(k), y = h_tot_cr_flux->GetBinContent(k);  

			// printf("k=%d, x=%f, y=%f \n", k, x, y); 

			h_contribute[j]->SetBinContent(k, f_fit[j]->Eval(x));
			h_contribute_mw[j]->SetBinContent(k, A[j]*f_fit[j]->Eval(x));
			h_contribute[j]->SetBinError(k, 0); 
			h_contribute_mw[j]->SetBinError(k, 0); 

		} 

		// PRINT_HIST(h_contribute[j])

		c3->cd(1);
		gPad->SetLogx();
		gPad->SetLogy();
		gPad->SetGrid(); 

		// average regular ratio 

		h_contribute_ratio[j] = (TH1 *) h_contribute[j]->Clone("h_contribute_ratio"); 
		h_contribute_ratio[j]->SetXTitle(Unit::GetEnergyLabel("GV")); 
   		h_contribute_ratio[j]->SetYTitle("Flux / Total Flux");
		h_contribute_ratio[j]->Divide(h_tot_cr_flux); 

		double R_1 = h_tot_cr_flux->GetBinLowEdge(1);
		double R_2 = h_tot_cr_flux->GetBinLowEdge(nbins+1); 

		// relative contribution 

		double contribute_rel = h_contribute[j]->Integral(R_1, R_2)/h_tot_cr_flux->Integral(R_1, R_2); 			

		printf("%s contribution = %f \n", ACE_Element[i], contribute_rel); 

		// mass-weighted contribution 

		h_contribute_mw_ratio[j] = (TH1 *) h_contribute_mw[j]->Clone("h_contribute_mw_ratio"); 
		h_contribute_mw_ratio[j]->Divide(h_tot_cr_flux_mw); 
		HistTools::SetStyle(h_contribute_mw_ratio[j], kBlue, kFullCircle, 0.9, 1, 1); 

		double contribute_rel_mw = h_contribute_mw_ratio[j]->Integral(R_1, R_2)/h_tot_cr_flux_mw->Integral(R_1, R_2); 			

		printf("%s mass-weighted contribution = %f \n", ACE_Element[i], contribute_rel_mw);

		TLegend *legend2 = new TLegend(0.62,0.8,0.9,0.9); 
		legend2->AddEntry(h_contribute_ratio[j], "normal contribution");
		legend2->AddEntry(h_contribute_mw_ratio[j], "mass-weighted contribution");

		h_contribute_mw_ratio[j]->SetStats(0);
		h_contribute_ratio[j]->SetTitle(Form("%s Contribution = %10.6f", ACE_Element[i], contribute_rel));
		h_contribute_mw_ratio[j]->SetTitle(Form("%s MW Contribution = %10.6f", ACE_Element[i], contribute_rel_mw));
		h_contribute_mw_ratio[j]->SetXTitle(Unit::GetEnergyLabel("GV")); 
   		h_contribute_mw_ratio[j]->SetYTitle("Flux / Total Flux");
		// h_contribute_mw_ratio[j]->GetYaxis()->SetRangeUser(10e-4, 10e-2);
		h_contribute_mw_ratio[j]->Draw("HIST");
		h_contribute_ratio[j]->Draw("HIST SAME");
		legend2->Draw("SAME"); 

		c3->cd(2);
		gPad->SetLogx();
		gPad->SetLogy();
		gPad->SetGrid();

		TH1 *h_ratio = (TH1 *) h_contribute_mw_ratio[j]->Clone("h_ratio"); 
		h_ratio->Divide(h_contribute_ratio[j]); 

		TLegend *legend1 = new TLegend(0.62,0.8,0.9,0.9); 
		// legend1->AddEntry(h_contribute[j], "normal flux");
		// legend1->AddEntry(h_contribute_mw[j], "mass-weighted flux");

		h_ratio->SetStats(0);
		h_ratio->SetTitle(Form("%s Contribution / Mass-weighted Contribution", ACE_Element[i]));
		h_ratio->SetXTitle(Unit::GetEnergyLabel("GV"));
   		h_ratio->SetYTitle("Ratio");
		h_ratio->Draw("HIST");
		// legend1->Draw("SAME"); 

		c3->Print(Form("./data/ACE/contribute/%s_temp/contribution_%d_%s.png", ACE_Element[temp], j, ACE_Element[i]));

	} 

	TLegend *legend_hs = new TLegend(0.85,0.4,0.9,0.9);
	legend_hs->SetNColumns(2); 
	TLegend *legend_hsmw = new TLegend(0.85,0.4,0.9,0.9);
	legend_hsmw->SetNColumns(2); 

	TFile file(Form("data/ACE/contribute/%s_temp/h_contribute.root", ACE_Element[temp]), "RECREATE"); 
	
	// input THStack for relative contribution comparison
	// save the h_contribute_ratio histograms
	for (int i = 0; i<28; ++i){
	
		HistTools::SetStyle(h_contribute_ratio[i], HistTools::GetColorPalette(i, n_total), 21, 0.9, 1, 1);
		h_contribute_ratio[i]->SetFillColor(HistTools::GetColorPalette(i, n_total));
		legend_hs->AddEntry(h_contribute_ratio[i], Form("%s", Element[i])); 
		hs->Add(h_contribute_ratio[i]);

		HistTools::SetStyle(h_contribute_mw_ratio[i], HistTools::GetColorPalette(i, n_total), 21, 0.9, 1, 1); 
		h_contribute_mw_ratio[i]->SetFillColor(HistTools::GetColorPalette(i, n_total)); 
		legend_hsmw->AddEntry(h_contribute_mw_ratio[i], Form("%s", Element[i])); 
		hs_mw->Add(h_contribute_mw_ratio[i]);

		h_contribute_ratio[i]->Write(Form("h_contribute_ratio_%s", Element[i]));
		h_contribute_mw_ratio[i]->Write(Form("h_contribute_mw_ratio_%s", Element[i]));   

		f_fit[i]->Write(Form("fit_%s", Element[i])); 
		
	} 

	h_tot_cr_flux->Write("h_tot_cr_flux");
	h_tot_cr_flux_mw->Write("h_tot_cr_flux_mw"); 

	file.Write();
	file.Close();

	c4->cd(1);
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGrid();

	//TH1 *h_a1 = HistTools::CreateAxis("h_axis", "Comparison of Relative Contributions; ; Flux / Total Flux", 0, 3000., 7, 0, 1, false);
	//h_a1->SetXTitle(Unit::GetEnergyLabel("GV")); 

	//h_a1->Draw("E1X0");

	hs->Draw("HIST"); // this will draw the relative contribution histograms stacked, i.e. one on top of the other, so that you can easily see which element is the dominant one in every rigidity bin. “pfc” will automatically choose a different color for each histogram from the current color palette.
	hs->GetYaxis()->SetTitle("Relative Contribution"); 
	hs->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV")); 
	legend_hs->Draw("SAME");

	c4->cd(2);
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGrid();

	//TH1 *h_a2 = HistTools::CreateAxis("h_axis", "Comparison of Mass-weighted Contributions; ; Mass-weighted Flux / Mass-weighted Total Flux", 0, 3000., 7, 0, 1, false);
	//h_a2->SetXTitle(Unit::GetEnergyLabel("GV")); 

	//h_a2->Draw("E1X0");

	hs_mw->Draw("HIST");
	hs_mw->GetYaxis()->SetTitle("Mass-weighted Relative Contribution"); 
	hs_mw->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV")); 
	legend_hsmw->Draw("SAME");

	c4->Print(Form("./data/ACE/contribute/%s_temp/contribution_stack.png", ACE_Element[temp])); 
	

	if (strcmp(option, "n") == 0) return h_tot_cr_flux; 
	if (strcmp(option, "mw") == 0) return h_tot_cr_flux_mw; 

}

// create the element contribution vs. rigidity table for BCNO templates 
// f_fit[i] 9-28 has not being saved into root files 
// understood the histogram bad, discarding this section 

void ace_contribute2(){	

	Debug::Enable(Debug::ALL); 

	Experiments::DataPath = "data"; 

	gSystem->mkdir("data/ACE/contribute2", true);

	for (int i=0; i<4; i++){ 

		gSystem->mkdir(Form("data/ACE/contribute2/%s_temp", ACE_Element[i]), true);

	} 

	int data_value[n_ele] = {0, 18, 23, 27, 31, 20, 43, 22}; 

	const char *AMS_Element2[8] = { "p", "he", "li", "be", "b", "c", "n", "o" }; 
	const char *AMS_Element2_Cap[8] = { "Proton", "He", "Li", "Be", "B", "C", "N", "O" }; 


	const UInt_t FirstACEBR = 2240;
   	vector<UInt_t> BRs;
  	// we stop at BR 2493, which ends on 2016/05/23, just 3 days before the end of the data taking period for AMS nuclei
   	for (UInt_t br=2426; br<=2493; ++br) { 
	   	if (br != 2472 && br != 2473) BRs.push_back(br-FirstACEBR); 
	}

	TH1 *h_ams_C = Experiments::GetMeasurementHistogram(Experiments::AMS02, data_value[1], 0); // load AMS Template data in order to determine Rmax

	TH1 *h_ene_C = HistTools::GraphToHist(get_ace_average_graph( ACE_Element[1] , &BRs[0], BRs.size() ), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
	TH1 *h_ace_C = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene_C, ACE_Isotope[1], "MeV/n cm", "GV m", "_rig"); // load averaged ACE data for the same element in rigidity

	UShort_t namsbins = h_ams_C->GetNbinsX(); 
	UShort_t nacebins = h_ace_C->GetNbinsX(); 
	
	double R1 = h_ace_C->GetBinLowEdge(1);
	double R2 = h_ams_C->GetBinLowEdge(namsbins+1); 

	//printf("R1 = %f, R2 = %f \n", R1, R2); 

	gStyle->SetOptStat(0);
	//gStyle->SetPalette(109); 

	TF1 *f_fit[n_total];
	
	//gPad->SetLogx();
	//gPad->SetLogy();
	//gPad->SetGrid();

	// Sum up the Total Cosmic Ray Flux 
	// p, He, Li, Be
	for (int i=0;i<4;i++){
		
		f_fit[i] = new TF1(); 

		int nnodes = 9;  
		int nnodes_ams = 6; 
		
		TFile file0(Form("data/ACE/extend2/fit_%s_temp_%s_%dnodes.root", ACE_Element[i_temp[i]], AMS_Element2[i], nnodes)); 
		cout << AMS_Element2_Cap[i] <<" "<< ACE_Element[i_temp[i]] << endl;

		Spline *sp_ams = new Spline("sp_ams", nnodes, Spline::LogLog | Spline::PowerLaw);
		f_fit[i] = sp_ams->GetTF1Pointer(); 
		TF1 *fit_ams = (TF1*) file0.Get(Form("fsp_%s", AMS_Element2[i])); // load the fluxes, the fsp in the root file is actually fit 

		HistTools::CopyParameters(fit_ams, f_fit[i]);
		double x1, x2;
		fit_ams->GetRange(x1,x2); 
		f_fit[i]->SetRange(x1,x2);   
			 	  
	} 

	// B, C, N, O 
	for (int i=0;i<4;i++){ 

		int j = i+4; 		

		int nnodes = 7;
	
		f_fit[j] = new TF1(); 

		TFile file0(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[i], nnodes)); // load ACE combined fit 

		cout << ACE_Element[i] << endl; 

		Spline *sp_comb = new Spline("sp_comb", nnodes, Spline::LogLog | Spline::PowerLaw); 
		f_fit[j] = sp_comb->GetTF1Pointer();  // real function 
		TF1 *fit_comb = (TF1*) file0.Get("fit_both"); 

		HistTools::CopyParameters(fit_comb, f_fit[j]); 
		double x1, x2; 
		fit_comb->GetRange(x1,x2); 
		f_fit[j]->SetRange(x1,x2);
		
	}

	// F, Ne, Na ... Va, Ni
	for (int i=4;i<24;i++){ 

		int j = i+4; 		

		int nnodes = 7; 
	
		f_fit[j] = new TF1();  

		TH1 *h_ene = HistTools::GraphToHist(get_ace_average_graph( ACE_Element[i] , &BRs[0], BRs.size() ), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
		TH1 *h_ace = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene, ACE_Isotope[i], "MeV/n cm", "GV m", "_rig"); // load averaged ACE data for the same element in rigidity 
		HistTools::SetStyle(h_ace, HistTools::GetColorPalette(j, n_total) , kFullCircle, 1.5, 1, 1); 

		TFile file0(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[i_temp[j]], nnodes)); // load ACE Template combined fit 

		Spline *sp_comb_C = new Spline("sp_comb_C", nnodes, Spline::LogLog | Spline::PowerLaw); 
		f_fit[j] = sp_comb_C->GetTF1Pointer();  // real function 
		TF1 *fit_comb_C = (TF1*) file0.Get("fit_both"); 
	
		HistTools::CopyParameters(fit_comb_C, f_fit[j]); 
		double x1_C, x2_C;
		fit_comb_C->GetRange(x1_C,x2_C); 
		f_fit[j]->SetRange(x1_C,x2_C); 
		
		TH1 *h_ratio = (TH1 *) h_ace->Clone("h_ratio");
	
		h_ratio->Divide(f_fit[j]); // WARNING! Doing this will change the # of entires of the histogram. 
		HistTools::SetStyle(h_ratio, kRed, kFullCircle, 0.9, 1, 1);
	
		double ratio_sum=0; // compute average of h_ratio manually  
		for(int k=0;k<14;k++){ 
			ratio_sum += h_ratio->GetBinContent(k); 
		}
		double ratio_ave = ratio_sum/7; 
	
		f_fit[j] = HistTools::CombineTF1Const(f_fit[j], ratio_ave, HistTools::MultiplyConst, "f_fit", R1, R2); 

// ####

	}

	Int_t nbins = 200;

	Double_t *bins = HistTools::BuildLogBins(R1, R2, nbins); // Rmin is the minimum rigidity of the ACE fluxes, Rmax the maximum rigidity of the AMS fluxes
	TH1D *h_tot_cr_flux = new TH1D("h_tot_cr_flux", "Total Cosmic Ray Flux;Rigidity [GV];Flux [1/(m^{2} sr s GV)]", nbins, bins);
	TH1D *h_tot_cr_flux_mw = new TH1D("h_tot_cr_flux_mw", "Mass-weighted Total Cosmic Ray Flux;Rigidity [GV];Flux [1/(m^{2} sr s GV)]", nbins, bins); // Contributed by each flux multiplied by their mass number 
	HistTools::SetStyle(h_tot_cr_flux, kRed, kFullCircle, 0.9, 1, 1); 
	HistTools::SetStyle(h_tot_cr_flux_mw, kBlue, kFullCircle, 0.9, 1, 1); 
	
	for (int i = 0; i < n_total; ++i) 
	{
	
   		for (int bin = 1; bin <= nbins; ++bin)
   		{
      			Double_t R = h_tot_cr_flux->GetBinLowEdge(bin); 
     			Double_t w = h_tot_cr_flux->GetBinWidth(bin); 
      			Double_t flux = f_fit[i]->Integral(R, R+w)/w; 
      			h_tot_cr_flux->AddBinContent(bin, flux);
			h_tot_cr_flux_mw->AddBinContent(bin, flux*A[i]); 

			// printf("i=%d, bin=%d, flux=%f \n", i, bin, flux); 
   		}

		//break; 
	}

	TCanvas *c2 = new TCanvas("c2", "Total Cosmic Ray Flux", 1800, 900); 
	c2->Divide(1, 2); 
	c2->cd(1); 
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGrid();

	TLegend *legend = new TLegend(0.62,0.8,0.9,0.9); 
	legend->AddEntry(h_tot_cr_flux, "normal");
	legend->AddEntry(h_tot_cr_flux_mw, "mass-weighted"); 

	h_tot_cr_flux->Draw("HIST");
	h_tot_cr_flux_mw->Draw("HIST SAME"); 
	legend->Draw("SAME"); 

	// PRINT_HIST(h_tot_cr_flux) 

	c2->cd(2); 
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGrid();

	TH1 *h_tot_cr_ratio = (TH1 *) h_tot_cr_flux_mw->Clone("h_tot_cr_ratio");
	h_tot_cr_ratio->Divide(h_tot_cr_flux);

	h_tot_cr_ratio->SetTitle("Total Cosmic Ray Flux / Mass-weighted Total Cosmic Ray Flux"); 
	h_tot_cr_ratio->SetXTitle(Unit::GetEnergyLabel("GV"));
	h_tot_cr_ratio->SetYTitle("Total Flux / Mass-weighted Total Flux");  
	h_tot_cr_ratio->Draw("HIST"); 

	// PRINT_HIST(h_tot_cr_flux_mw) 

	c2->Print("./data/ACE/contribute2/total_cr_flux.png"); 

	TCanvas *c3 = new TCanvas("c3", "Contribution", 1800, 900); 
	c3->Divide(1, 2); 

	TCanvas *c4 = new TCanvas("c4", "Comparison of Relative Contributions", 1800, 900);
	c4->Divide(1, 2);  

	THStack *hs = new THStack("hs", "Comparison of Relative Contributions");
	THStack *hs_mw = new THStack("hs_mw", "Comparison of Mass-Weighted Relative Contributions");

	// Now compute the ratio of h_element vs. h_tot_cr_flux 

	TH1 *h_contribute[n_total]; // flux
	TH1 *h_contribute_ratio[n_total]; // flux / total
	TH1 *h_contribute_mw[n_total]; // mass-weighted flux
	TH1 *h_contribute_mw_ratio[n_total]; // mass-weighted flux / total 

	// p, He, Li, Be
	for (int i=0;i<4;i++){

		int nnodes = 6;  

		h_contribute[i] = (TH1 *) h_tot_cr_flux->Clone("h_contribute");  
		h_contribute_mw[i] = (TH1 *) h_tot_cr_flux->Clone("h_contribute_mw"); 
		HistTools::SetStyle(h_contribute_mw[i], kBlue, kFullCircle, 0.9, 1, 1); 

		for (int k=1;k<=nbins;k++){ 

			Double_t R = h_tot_cr_flux->GetBinLowEdge(k);  
     			Double_t w = h_tot_cr_flux->GetBinWidth(k);  
      			Double_t flux = f_fit[i]->Integral(R, R+w)/w;  

			h_contribute[i]->SetBinContent(k, flux);
			h_contribute_mw[i]->SetBinContent(k, A[i]*flux);
			h_contribute[i]->SetBinError(k, 0); 
			h_contribute_mw[i]->SetBinError(k, 0); 

		} 

		// PRINT_HIST(h_contribute[i])

		c3->cd(1);
		gPad->SetLogx();
		gPad->SetLogy();
		gPad->SetGrid(); 

		// relative contribution 

		h_contribute_ratio[i] = (TH1 *) h_contribute[i]->Clone("h_contribute_ratio"); 
		h_contribute_ratio[i]->SetXTitle(Unit::GetEnergyLabel("GV")); 
   		h_contribute_ratio[i]->SetYTitle("Flux / Total Flux"); 
		h_contribute_ratio[i]->Divide(h_tot_cr_flux);   

		double R_1 = h_tot_cr_flux->GetBinLowEdge(1);
		double R_2 = h_tot_cr_flux->GetBinLowEdge(nbins+1);
		
		double contribute_rel = h_contribute[i]->Integral(R_1, R_2)/h_tot_cr_flux->Integral(R_1, R_2); 	

		printf("%s contribution = %f \n", AMS_Element2_Cap[i], contribute_rel); 

		// mass-weighted relative contribution 

		h_contribute_mw_ratio[i] = (TH1 *) h_contribute_mw[i]->Clone("h_contribute_mw_ratio"); 
		h_contribute_mw_ratio[i]->Divide(h_tot_cr_flux_mw); 
		HistTools::SetStyle(h_contribute_mw_ratio[i], kBlue, kFullCircle, 0.9, 1, 1);
		
		double contribute_rel_mw = h_contribute_mw[i]->Integral(R_1, R_2)/h_tot_cr_flux_mw->Integral(R_1, R_2); 			

		printf("%s mass-weighted contribution = %f \n", AMS_Element2_Cap[i], contribute_rel_mw);

		TLegend *legend2 = new TLegend(0.62,0.8,0.9,0.9); 
		legend2->AddEntry(h_contribute_ratio[i], "normal contribution");
		legend2->AddEntry(h_contribute_mw_ratio[i], "mass-weighted contribution");

		h_contribute_mw_ratio[i]->SetStats(0);
		h_contribute_ratio[i]->SetTitle(Form("%s Contribution = %10.6f", AMS_Element2_Cap[i], contribute_rel));
		h_contribute_mw_ratio[i]->SetTitle(Form("%s MW Contribution = %10.6f", AMS_Element2_Cap[i], contribute_rel_mw));
		h_contribute_mw_ratio[i]->SetXTitle(Unit::GetEnergyLabel("GV")); 
   		h_contribute_mw_ratio[i]->SetYTitle("Flux / Total Flux");
		h_contribute_mw_ratio[i]->Draw("HIST");
		h_contribute_ratio[i]->Draw("HIST SAME");
		legend2->Draw("SAME");

		c3->cd(2);
		gPad->SetLogx();
		gPad->SetLogy();
		gPad->SetGrid();

		TH1 *h_ratio = (TH1 *) h_contribute_mw_ratio[i]->Clone("h_ratio"); 
		h_ratio->Divide(h_contribute_ratio[i]); 

		TLegend *legend1 = new TLegend(0.62,0.8,0.9,0.9); 
		// legend1->AddEntry(h_contribute[i], "normal flux");
		// legend1->AddEntry(h_contribute_mw[i], "mass-weighted flux");

		h_ratio->SetStats(0);
		h_ratio->SetTitle(Form("%s Contribution / Mass-weighted Contribution", AMS_Element2_Cap[i]));
		h_ratio->SetXTitle(Unit::GetEnergyLabel("GV"));
   		h_ratio->SetYTitle("Ratio");
		h_ratio->Draw("HIST");
		// legend1->Draw("SAME"); 

		c3->Print(Form("./data/ACE/contribute2/%s_temp/contribution_%d_%s.png", ACE_Element[i_temp[i]], i, AMS_Element2_Cap[i])); 

	}

	// B, C, N, O
	for (int i=0;i<4;i++){ 

		int j = i+4; 		

		int nnodes = 6; 

		h_contribute[j] = (TH1 *) h_tot_cr_flux->Clone("h_contribute");  
		h_contribute_mw[j] = (TH1 *) h_tot_cr_flux->Clone("h_contribute_mw"); 
		HistTools::SetStyle(h_contribute_mw[j], kBlue, kFullCircle, 0.9, 1, 1); 

		for (int k=1;k<=nbins;k++){ 

			Double_t x = h_tot_cr_flux->GetBinCenter(k), y = h_tot_cr_flux->GetBinContent(k);  

			// printf("k=%d, x=%f, y=%f \n", k, x, y); 

			h_contribute[j]->SetBinContent(k, f_fit[j]->Eval(x));
			h_contribute_mw[j]->SetBinContent(k, A[j]*f_fit[j]->Eval(x));
			h_contribute[j]->SetBinError(k, 0); 
			h_contribute_mw[j]->SetBinError(k, 0); 

		} 

		// PRINT_HIST(h_contribute[j])

		c3->cd(1);
		gPad->SetLogx();
		gPad->SetLogy();
		gPad->SetGrid(); 

		// relative contribution

		h_contribute_ratio[j] = (TH1 *) h_contribute[j]->Clone("h_contribute_ratio");
		h_contribute_ratio[j]->SetXTitle(Unit::GetEnergyLabel("GV")); 
   		h_contribute_ratio[j]->SetYTitle("Flux / Total Flux"); 
		h_contribute_ratio[j]->Divide(h_tot_cr_flux); 

		double R_1 = h_tot_cr_flux->GetBinLowEdge(1);
		double R_2 = h_tot_cr_flux->GetBinLowEdge(nbins+1);

		double contribute_rel = h_contribute[j]->Integral(R_1, R_2)/h_tot_cr_flux->Integral(R_1, R_2); 			

		printf("%s contribution = %f \n", ACE_Element[i], contribute_rel); 

		// mass-weighted ratio

		h_contribute_mw_ratio[j] = (TH1 *) h_contribute_mw[j]->Clone("h_contribute_mw_ratio"); 
		h_contribute_mw_ratio[j]->Divide(h_tot_cr_flux_mw); 
		HistTools::SetStyle(h_contribute_mw_ratio[j], kBlue, kFullCircle, 0.9, 1, 1); 

		double contribute_rel_mw = h_contribute_mw[j]->Integral(R_1, R_2)/h_tot_cr_flux_mw->Integral(R_1, R_2); 			

		printf("%s mass-weighted contribution = %f \n", ACE_Element[i], contribute_rel_mw); 

		TLegend *legend2 = new TLegend(0.62,0.8,0.9,0.9); 
		legend2->AddEntry(h_contribute_ratio[j], "normal contribution");
		legend2->AddEntry(h_contribute_mw_ratio[j], "mass-weighted contribution");

		h_contribute_mw_ratio[j]->SetStats(0);
		h_contribute_ratio[j]->SetTitle(Form("%s Contribution = %10.6f", ACE_Element[i], contribute_rel));
		h_contribute_mw_ratio[j]->SetTitle(Form("%s MW Contribution = %10.6f", ACE_Element[i], contribute_rel_mw));
		h_contribute_mw_ratio[j]->SetXTitle(Unit::GetEnergyLabel("GV")); 
   		h_contribute_mw_ratio[j]->SetYTitle("Flux / Total Flux");
		// h_contribute_mw_ratio[j]->GetYaxis()->SetRangeUser(10e-4, 10e-2);
		h_contribute_mw_ratio[j]->Draw("HIST");
		h_contribute_ratio[j]->Draw("HIST SAME");
		legend2->Draw("SAME"); 

		c3->cd(2);
		gPad->SetLogx();
		gPad->SetLogy();
		gPad->SetGrid();

		TH1 *h_ratio = (TH1 *) h_contribute_mw_ratio[j]->Clone("h_ratio"); 
		h_ratio->Divide(h_contribute_ratio[j]); 

		TLegend *legend1 = new TLegend(0.62,0.8,0.9,0.9); 
		// legend1->AddEntry(h_contribute[j], "normal flux");
		// legend1->AddEntry(h_contribute_mw[j], "mass-weighted flux");

		h_ratio->SetStats(0);
		h_ratio->SetTitle(Form("%s Contribution / Mass-weighted Contribution", ACE_Element[i]));
		h_ratio->SetXTitle(Unit::GetEnergyLabel("GV"));
   		h_ratio->SetYTitle("Ratio");
		h_ratio->Draw("HIST");
		// legend1->Draw("SAME"); 

		c3->Print(Form("./data/ACE/contribute2/%s_temp/contribution_%d_%s.png", ACE_Element[i_temp[j]], j, ACE_Element[i])); 
 
	}

	// F, Ne, Na ... Va, Ni
	for (int i=4;i<24;i++){ 

		int j = i+4; 		

		h_contribute[j] = (TH1 *) h_tot_cr_flux->Clone("h_contribute");  
		h_contribute_mw[j] = (TH1 *) h_tot_cr_flux->Clone("h_contribute_mw"); 
		HistTools::SetStyle(h_contribute_mw[j], kBlue, kFullCircle, 0.9, 1, 1); 

		for (int k=1;k<=nbins;k++){ 

			Double_t x = h_tot_cr_flux->GetBinCenter(k), y = h_tot_cr_flux->GetBinContent(k);  

			// printf("k=%d, x=%f, y=%f \n", k, x, y); 

			h_contribute[j]->SetBinContent(k, f_fit[j]->Eval(x));
			h_contribute_mw[j]->SetBinContent(k, A[j]*f_fit[j]->Eval(x));
			h_contribute[j]->SetBinError(k, 0); 
			h_contribute_mw[j]->SetBinError(k, 0); 

		} 

		// PRINT_HIST(h_contribute[j])

		c3->cd(1);
		gPad->SetLogx();
		gPad->SetLogy();
		gPad->SetGrid(); 

		// average regular ratio 

		h_contribute_ratio[j] = (TH1 *) h_contribute[j]->Clone("h_contribute_ratio"); 
		h_contribute_ratio[j]->SetXTitle(Unit::GetEnergyLabel("GV")); 
   		h_contribute_ratio[j]->SetYTitle("Flux / Total Flux");
		h_contribute_ratio[j]->Divide(h_tot_cr_flux); 

		double R_1 = h_tot_cr_flux->GetBinLowEdge(1);
		double R_2 = h_tot_cr_flux->GetBinLowEdge(nbins+1); 

		// relative contribution 

		double contribute_rel = h_contribute[j]->Integral(R_1, R_2)/h_tot_cr_flux->Integral(R_1, R_2); 			

		printf("%s contribution = %f \n", ACE_Element[i], contribute_rel); 

		// mass-weighted contribution 

		h_contribute_mw_ratio[j] = (TH1 *) h_contribute_mw[j]->Clone("h_contribute_mw_ratio"); 
		h_contribute_mw_ratio[j]->Divide(h_tot_cr_flux_mw); 
		HistTools::SetStyle(h_contribute_mw_ratio[j], kBlue, kFullCircle, 0.9, 1, 1); 

		double contribute_rel_mw = h_contribute_mw_ratio[j]->Integral(R_1, R_2)/h_tot_cr_flux_mw->Integral(R_1, R_2); 			

		printf("%s mass-weighted contribution = %f \n", ACE_Element[i], contribute_rel_mw);

		TLegend *legend2 = new TLegend(0.62,0.8,0.9,0.9); 
		legend2->AddEntry(h_contribute_ratio[j], "normal contribution");
		legend2->AddEntry(h_contribute_mw_ratio[j], "mass-weighted contribution");

		h_contribute_mw_ratio[j]->SetStats(0);
		h_contribute_ratio[j]->SetTitle(Form("%s Contribution = %10.6f", ACE_Element[i], contribute_rel));
		h_contribute_mw_ratio[j]->SetTitle(Form("%s MW Contribution = %10.6f", ACE_Element[i], contribute_rel_mw));
		h_contribute_mw_ratio[j]->SetXTitle(Unit::GetEnergyLabel("GV")); 
   		h_contribute_mw_ratio[j]->SetYTitle("Flux / Total Flux");
		// h_contribute_mw_ratio[j]->GetYaxis()->SetRangeUser(10e-4, 10e-2);
		h_contribute_mw_ratio[j]->Draw("HIST");
		h_contribute_ratio[j]->Draw("HIST SAME");
		legend2->Draw("SAME"); 

		c3->cd(2);
		gPad->SetLogx();
		gPad->SetLogy();
		gPad->SetGrid();

		TH1 *h_ratio = (TH1 *) h_contribute_mw_ratio[j]->Clone("h_ratio"); 
		h_ratio->Divide(h_contribute_ratio[j]); 

		TLegend *legend1 = new TLegend(0.62,0.8,0.9,0.9); 
		// legend1->AddEntry(h_contribute[j], "normal flux");
		// legend1->AddEntry(h_contribute_mw[j], "mass-weighted flux");

		h_ratio->SetStats(0);
		h_ratio->SetTitle(Form("%s Contribution / Mass-weighted Contribution", ACE_Element[i]));
		h_ratio->SetXTitle(Unit::GetEnergyLabel("GV"));
   		h_ratio->SetYTitle("Ratio");
		h_ratio->Draw("HIST");
		// legend1->Draw("SAME"); 

		c3->Print(Form("./data/ACE/contribute2/%s_temp/contribution_%d_%s.png", ACE_Element[i_temp[j]], j, ACE_Element[i]));

	} 

	TLegend *legend_hs = new TLegend(0.85,0.4,0.9,0.9);
	legend_hs->SetNColumns(2); 
	TLegend *legend_hsmw = new TLegend(0.85,0.4,0.9,0.9);
	legend_hsmw->SetNColumns(2); 

	TFile file("data/ACE/contribute2/h_contribute.root", "RECREATE"); 
	
	// input THStack for relative contribution comparison
	// save the h_contribute_ratio histograms
	for (int i = 0; i<28; i++){
	
		HistTools::SetStyle(h_contribute_ratio[i], HistTools::GetColorPalette(i, n_total), 21, 0.9, 1, 1);
		h_contribute_ratio[i]->SetFillColor(HistTools::GetColorPalette(i, n_total));
		legend_hs->AddEntry(h_contribute_ratio[i], Form("%s", Element[i])); 
		hs->Add(h_contribute_ratio[i]);

		HistTools::SetStyle(h_contribute_mw_ratio[i], HistTools::GetColorPalette(i, n_total), 21, 0.9, 1, 1); 
		h_contribute_mw_ratio[i]->SetFillColor(HistTools::GetColorPalette(i, n_total)); 
		legend_hsmw->AddEntry(h_contribute_mw_ratio[i], Form("%s", Element[i])); 
		hs_mw->Add(h_contribute_mw_ratio[i]);

		h_contribute_ratio[i]->Write(Form("h_contribute_ratio_%s", Element[i]));
		h_contribute_mw_ratio[i]->Write(Form("h_contribute_mw_ratio_%s", Element[i]));   

		f_fit[i]->Write(Form("fit_%s", Element[i])); 
		
	} 

	h_tot_cr_flux->Write("h_tot_cr_flux");
	h_tot_cr_flux_mw->Write("h_tot_cr_flux_mw"); 

	file.Write();
	file.Close();

	c4->cd(1);
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGrid();

	//TH1 *h_a1 = HistTools::CreateAxis("h_axis", "Comparison of Relative Contributions; ; Flux / Total Flux", 0, 3000., 7, 0, 1, false);
	//h_a1->SetXTitle(Unit::GetEnergyLabel("GV")); 

	//h_a1->Draw("E1X0");

	hs->Draw("HIST"); // this will draw the relative contribution histograms stacked, i.e. one on top of the other, so that you can easily see which element is the dominant one in every rigidity bin. “pfc” will automatically choose a different color for each histogram from the current color palette.
	hs->GetYaxis()->SetTitle("Relative Contribution"); 
	hs->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV")); 
	legend_hs->Draw("SAME");

	c4->cd(2);
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGrid();

	//TH1 *h_a2 = HistTools::CreateAxis("h_axis", "Comparison of Mass-weighted Contributions; ; Mass-weighted Flux / Mass-weighted Total Flux", 0, 3000., 7, 0, 1, false);
	//h_a2->SetXTitle(Unit::GetEnergyLabel("GV")); 

	//h_a2->Draw("E1X0");

	hs_mw->Draw("HIST");
	hs_mw->GetYaxis()->SetTitle("Mass-weighted Relative Contribution"); 
	hs_mw->GetXaxis()->SetTitle(Unit::GetEnergyLabel("GV")); 
	legend_hsmw->Draw("SAME");

	c4->Print("./data/ACE/contribute2/contribution_stack.png"); 

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

TH1* ave_hist( TH1D **h_set, int nBRs ){

	TH1 *h_ave = (TH1*) h_set[0]->Clone("h_ave"); 

	for (int bin=1; bin<=h_set[0]->GetNbinsX(); ++bin){	

		double ave=0, sum=0; 	

		for (int iBR=0; iBR<nBRs; ++iBR){
			sum += h_set[iBR]->GetBinContent(bin); 
		}

		ave = sum/nBRs;  
		h_ave->SetBinContent(bin, ave); 
	}

	return h_ave; 

} 

TGraphErrors *ave_grapherrors( TGraphErrors *g_set[], int nset ){

	TGraphErrors *g_ave = (TGraphErrors*) g_set[0]->Clone("g_ave"); 

	for (int point=0; point<g_set[0]->GetN(); ++point){	

		double ave=0, sum=0, std=0, stdsum=0; 	

		for (int iset=0; iset<nset; ++iset){
			sum += g_set[iset]->GetY()[point]; 
		}

		ave = sum/nset; 

		for (int iset=0; iset<nset; ++iset){
			stdsum += pow(g_set[iset]->GetY()[point]-ave, 2); 
		}
 
		std = sqrt(stdsum/(nset-1)); 
		g_ave->SetPoint(point, g_set[0]->GetX()[point], ave); 
		g_ave->SetPointError(point, 0, std); 
	}

	HistTools::CopyStyle(g_set[0], g_ave); 
	return g_ave; 
} 

const char *get_template(const char *element){

	for (int i=0; i<n_total; i++){
		if (!strcmp(element, Form("%s", Element[i]))) return ACE_Element[i_temp[i]];
	}
}

int get_ams_data_value(const char *element){

	int data_value[8] = {0, 18, 23, 27, 31, 20, 43, 22}; // p, He, Li, Be, B, C, N, O 

	for (int i=0; i<8; i++){
		if (!strcmp(element, Form("%s", Element[i]))) return data_value[i];
	}
	
	for (int i=8; i<28; i++){
		if (!strcmp(element, Form("%s", Element[i]))) return data_value[4+i_temp[i]]; 
	}
}

TGraph *get_norm_graph(TGraph *g){

	double ave=0, sum=0; 

	for (int i=0; i<g->GetN(); i++){

		Double_t x=0, y=0; 
		g->GetPoint(i, x, y); 
		sum += y; 
	} 

	ave = sum/g->GetN(); 

	TGraph *g_norm = new TGraph(); 

	for (int i=0; i<g->GetN(); i++){

		Double_t x=0, y=0; 
		g->GetPoint(i, x, y); 

		g_norm->SetPoint(i, x, y/ave); 

	} 

	HistTools::CopyStyle(g, g_norm); 

	return g_norm; 
}

TGraphErrors *get_norm_graph(TGraphErrors *g){

	double ave=0, sum=0; 
	double ex, ey, dsum=0; 

	for (int i=0; i<g->GetN(); i++){

		Double_t x=0, y=0;
		Double_t dy=0;  
		g->GetPoint(i, x, y);
		dy = g->GetErrorY(i);  
		sum += y; 
		dsum += dy*dy; 
	} 

	ave = sum/g->GetN(); 
	ey = sqrt(dsum)/g->GetN(); 

	TGraphErrors *g_norm = new TGraphErrors(); 

	for (int i=0; i<g->GetN(); i++){

		Double_t x=0, y=0; 
		Double_t dx=0; 
		g->GetPoint(i, x, y); 

		dx = g->GetErrorX(i); 
		ex = dx; 

		g_norm->SetPoint(i, x, y/ave); 
		g_norm->SetPointError(i, ex, ey); 
	} 

	HistTools::CopyStyle(g, g_norm); 

	return g_norm; 
}


/* double *getmean_ace( TGraphErrors *g, const char *option ){ 

	double ave=0, sum=0, std=0, stdsum=0; 

	for (int point=0; point<g->GetN(); ++point){	
		sum += g_set[iset]->GetY()[point]; 
	}

	ave = sum/g->GetN(); 

	for (int point=0; point<g->GetN(); ++point){
		stdsum += pow(g->GetY()[point]-ave, 2); 
	}
 
	std = sqrt(stdsum/g->GetN()-1)); 

	if (!strcmp(option, "ave")) return ave; 
	if (!strcmp(option, "std")) return std; 
	else return 0; 	
} */


