// Manipulate ACE/CRIS data by Lingqiang He, updated on 04/17/2019
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

const int n_ele = 24; // number of elements  

const char *ACE_Element[n_ele] = { "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "Va", "Cr", "Mn", "Fe", "Co", "Ni" };
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
void ace_fill(const char *element, Particle::Type isotope); // fill ACE/CRIS data into root histograms then save 
void ace_all_average(); // plot averaged flux over energy bins for all elements
void ace_convert(const char *element, Particle::Type isotope); // convert h_ene into h_rig and also plot flux and normalized flux over time
void ace_fitboth(int nnodes); // fit ACE&AMS combined data 
void compare_nodes(int k); // compare among spline fit results for 5-9 nodes 
double compare_sig(TF1 *fit1, TF1 *fit2); // obtain sigma value between two fits
void ace_extend(); // divide remaining ACE data that is not measured by AMS by combined fit to find a flat residual for ACE assumption at high energy region
		   // divide elements with Z>8 ACE data by combinced C cit 
void ace_extend2(); // fit_comb of BNO / C
void ace_extend3(); // group ratios from extend() part II and check isotope assumptions 

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

	gSystem->Setenv("TZ", "UCT"); 
	gStyle->SetTimeOffset(0);
	gStyle->SetOptStat(0); 
	gStyle->SetNumberContours(99); 
	gStyle->SetPalette(55); 

	const int nAce = 28-5+1; 
	// const char *elements[nAce] = { "B11", "C12", "N15", "O16", "F19", "Ne20", "Na23", "Mg24", "Al27", "Si28", "P31", "S32", "Cl35", "Ar36", "K41", "Ca40", "Sc45", "Ti46", "Va51", "Cr52", "Mn55", "Fe56", "Co59", "Ni60" }; list of chosen isotopes
	
	if (strcmp(operation, "fill") == 0){

		gROOT->ProcessLine(".> data/ACE/fill/fill_all.txt"); 
	
		ace_fill( "B", Particle::BORON11 );
		ace_fill( "C", Particle::CARBON12 );
		ace_fill( "N", Particle::NITROGEN15 );
		ace_fill( "O", Particle::OXYGEN16 );
		ace_fill( "F", Particle::FLUORINE19 );
		ace_fill( "Ne", Particle::NEON20 );
		ace_fill( "Na", Particle::SODIUM23 );
		ace_fill( "Mg", Particle::MAGNESIUM24 );
		ace_fill( "Al", Particle::ALUMINUM27 );
		ace_fill( "Si", Particle::SILICON28 );
		ace_fill( "P", Particle::PHOSPHORUS31 );
		ace_fill( "S", Particle::SULFUR32 );
		ace_fill( "Cl", Particle::CHLORINE35 );
		ace_fill( "Ar", Particle::ARGON36 );
		ace_fill( "K", Particle::POTASSIUM41 );
		ace_fill( "Ca", Particle::CALCIUM40 );
		ace_fill( "Sc", Particle::SCANDIUM45 );
		ace_fill( "Ti", Particle::TITANIUM46 );
		ace_fill( "Va", Particle::VANADIUM51 );
		ace_fill( "Cr", Particle::CHROMIUM52 );
		ace_fill( "Mn", Particle::MANGANESE55 );
		ace_fill( "Fe", Particle::IRON56 );
		ace_fill( "Co", Particle::COBALT59 );
		ace_fill( "Ni", Particle::NICKEL60 );

		gROOT->ProcessLine(".> ");

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
	}

} 

//Fill CRIS Data into Bins
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
	for (int k=0; k<ace->GetEntries(); k++){

		TH1F *h = new TH1F("h","", 13, kin_bins);
			
		ace->GetEntry(k); //Get the nth entry of TTree!! 						

		for (int i=0; i<h->GetNbinsX(); ++i) { 
			if (i%2==0) { 
				double sys_err = F[i/2] * sqrt(8e-4 + SpallCorrUnc[i/2]*SpallCorrUnc[i/2]/SpallCorr[i/2]/SpallCorr[i/2]);
				double stat_err = F[i/2]/sqrt(C[i/2]); 
				double tot_err = sqrt(stat_err*stat_err + sys_err*sys_err);
				h->SetBinContent(i+1, F[i/2]); 

				if (C[i/2]!=0){
					h->SetBinError(i+1, tot_err); 
				} else if (C[i/2]==0){
					continue;
				}
			} 
		}

		h->SetMarkerStyle(kFullCircle);
		//HistTools::SetColors(h, 290, kFullCircle, 1.4);
		gPad->SetGrid(); 
		h->SetTitle(Form("%s BR-%d Energy Spectrum; Energy (MeV/nuc); Flux (/(cm^2 sr s)(MeV/nuc)", element, UTimeToBR(utime)));
		//h->Draw("PSAME"); 

		h->Write(Form("h_kin_%s_BR%d", element, UTimeToBR(utime)));
	}

	file1.Write();
	file1.Close();

	_file0->Close();

}

// Plot All Element Averaged Flux over Energy Bins
void ace_all_average(){

   gSystem->mkdir("data/ACE/average", true);

   TCanvas *c1 = new TCanvas("c1","", 1600, 1200);
   c1->Divide(2, 1);

   TLegend *legend1 = new TLegend(0.1,0.6,0.28,0.9); // left, down, right, top
   TLegend *legend2 = new TLegend(0.1,0.6,0.28,0.9); // left, down, right, top

   // create axis histogram
   TH1 *h1 = HistTools::CreateAxis("h1", "haxis1", 0., 500., 7, 1e-10, 0.2e-6, false); 
   TH1 *h2 = HistTools::CreateAxis("h2", "haxis2", 0., 2.5, 7, 0.001, 1., false);

   c1->cd(1);
   gPad->SetLogy();
   h1->Draw();
   h1->SetTitle("All Element Average Kinetic Energy Spectrum");
   h1->SetXTitle(Unit::GetEnergyLabel("MeV/n"));
   h1->SetYTitle(Unit::GetDifferentialFluxLabel("MeV/n cm"));
   c1->cd(2);
   gPad->SetLogy();
   h2->Draw();
   h2->SetTitle("All Element Average Energy Spectrum");
   h2->SetXTitle(Unit::GetEnergyLabel("GV"));
   h2->SetYTitle(Unit::GetDifferentialFluxLabel("GV m"));

   for (int i=0; i<n_ele; i++){

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
	legend1->AddEntry(h_ene, Form("%s", ACE_Element[i]), "p");
	legend1->SetNColumns(2);
	h_ene->Draw("E1X0 SAME"); 
	legend1->Draw("SAME");	

	HistTools::SetMarkerStyle(h_rig, HistTools::GetColorPalette(i, n_ele), kFullCircle, 0.9);
	c1->cd(2);
	gPad->SetGrid();
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

			//h_rig->Write(Form("h_rig_%s_BR%d", element, k)); 
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
			h_rig->SetTitle(Form("%s All BR Energy Spectrum; Rigidity; Differential Flux", element));
			c2->cd(2);
			gPad->SetMargin(0.12, 0.08, 0.08, 0.08);
			h_rig->Draw("E1X0 SAME"); 

			//h_rig->Write(Form("h_rig_%s_BR%d", element, k)); 
		}
			
	}

	c2->Print(Form("./data/ACE/convert/fluxrigidity/h_rig_%s_all.png", element));

	// plot flux_time 
	TCanvas *c3 = new TCanvas("c3","",800,600);
	TLegend *legend3 = new TLegend(0.1,0.7,0.28,0.9); // left, down, right, top 
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
			g->SetTitle(Form("%s All Energy Bin Flux Time Series; ;Flux [1/(m^2 sr s GV)]", element));
			//g->Print();
			c3->cd(1);
			
			legend3->AddEntry(g, Form("%d-th Bin Flux", i/2+1), "l");

			if (i/2==0){
				g->Draw("ALP");
				legend3->Draw("SAME");
			} else {
				g->Draw("LPSAME");
				legend3->Draw("SAME");
			}
		}
	} 

	// plot normalized flux_time 
	TCanvas *c4 = new TCanvas("c4","",800,600);
	TLegend *legend4 = new TLegend(0.1,0.7,0.28,0.9); // left, down, right, top 

	const UInt_t FirstACEBR = 2240;
	vector<UInt_t> BRs;
	// we stop at BR 2493, which ends on 2016/05/23, just 3 days before the end of the data taking period for AMS nuclei
	for (UInt_t br=2426; br<=2493; ++br) { 
		if (br != 2472 && br != 2473) BRs.push_back(br-FirstACEBR); 
	}

	TH1D *h_ene = HistTools::GraphToHist( get_ace_average_graph( element, &BRs[0], BRs.size()), DBL_MIN, -DBL_MAX, true, 0.5, 0.);
	TH1 *h_ave = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene, isotope, "MeV/n cm", "GV m", "_rig");

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
			g_norm->SetTitle(Form("%s All Energy Bin Flux Time Series (Normalized)", element));
			//g_norm->Print();
			c4->cd(1);
			
			legend4->AddEntry(g_norm, Form("%d-th Bin Flux", i/2+1), "l");

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
   int data_value[n_ele] = {0, 18, 23, 27, 31, 20, 43, 22};

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
	   ynodes[inode] = h_ams->GetBinContent(h_ams->FindBin(x));
	}

	// create spline
	Spline *sp = new Spline("sp", nnodes, Spline::LogLog | Spline::PowerLaw, xnodes, ynodes);
	sp->SetSpectralIndices(3, -2.8, 0, 4, -3.5, -1); // set initial values and limits for the spectral index at first and last node
	TF1 *fit = sp->GetTF1Pointer();	
	
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

			TH1 *ha_ratio = HistTools::CreateAxis("ha_ratio", Form("Norm %s/%s;Rigidity [GV];", ACE_Element[i], ACE_Element[j]), 0.8, 2.5, 7, 0.6, 1.2, false);
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

		break; 
		c1->Print(Form("./data/ACE/extend/ACE_extend_residuals_byBCNO_%s.png", ACE_Element[i]));
	} 

	TCanvas *c2 = new TCanvas("c2","f_ratio residuals for remaining elements", 2400, 900);
	c2->Divide(2, 2);

	int nnodes = 7;
	TFile file2(Form("data/ACE/compare/fit_%s_%dnodes.root", ACE_Element[1], nnodes)); // load C combined fit

	// ACE Z>8 element data over C Combined Fit Normalized Ratio for varies isotopes 
	for (int i=4; i<24; i++){
	
		c2->cd(1);

		gPad->SetGrid();
		gPad->SetLogy();
		gPad->SetLogx();

		TH1 *ha = HistTools::CreateAxis("ha", Form("ACE %s Flux & C Combined Fit", ACE_Element[i]), 0.1, 2500., 7, 1e-10, 1e2, false);
		ha->Draw("E1X0");

		c2->cd(2);

		TH1 *ha_ratio = HistTools::CreateAxis("ha_ratio", Form("Normalized %s/%s;Rigidity [GV];", ACE_Element[i], ACE_Element[1]), 0.8, 2.4, 7, 0.7, 1.3, false);
		ha_ratio->Draw("E1X0");

		TLegend *legend1 = new TLegend(0.1,0.7,0.2,0.9); // left, down, right, top
		TLegend *legend2 = new TLegend(0.1,0.7,0.48,0.9); // left, down, right, top
		//legend1->SetNColumns(2);
		legend2->SetNColumns(2);

		TH1D *h_chi2 = new TH1D();  
		TH1D *h_var = new TH1D();

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
			double ratio_ave = ratio_sum/h_res->GetEntries();
			
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

			TF1 *f_flat = new TF1("f_flat","[0]"); // create flat ratio line   
			f_flat->SetParameter(0, 1);  	

			h_ratio->Fit(f_flat, "NQ"); 

			double chi2_flat=0; 

			for(int k=0;k<14;k++){
				if(k%2==0) chi2_flat += pow(h_ratio->GetBinContent(k+1)-1, 2);
				//printf("%0.7f \n", h_ratio->GetBinError(k+1)); 
			}

			//double chi2_flat = f_flat->GetChisquare(); 

			h_chi2->SetBinContent(j+1, chi2_flat); 

			//h_chi2->Print("range");

			h_chi2->GetXaxis()->SetBinLabel(j+1, name_isotope[i][j].c_str());  
			HistTools::SetStyle(h_chi2, kRed, kFullCircle, 0.9, 1, 1);
			h_chi2->SetTitle(Form("%s/%s Chi-2 Test;;(data-1)^2/error", ACE_Element[i], ACE_Element[1])); 
			h_chi2->Draw("HIST P"); 

			c2->cd(4);
			gPad->SetGrid();

			double var_flat=0;

			for(int k=0;k<14;k++){
				if(k%2==0) var_flat += abs(h_ratio->GetBinContent(k+1)-1)/7;
			}

			h_var->SetBinContent(j+1, var_flat); 

			h_var->GetXaxis()->SetBinLabel(j+1, name_isotope[i][j].c_str());  
			HistTools::SetStyle(h_var, kRed, kFullCircle, 0.9, 1, 1);
			h_var->SetTitle(Form("%s/%s Variation Test;;abs(data-1)", ACE_Element[i], ACE_Element[1])); 
			h_var->Draw("HIST P"); 

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

	const UInt_t FirstACEBR = 2240;
   	vector<UInt_t> BRs;
  	// we stop at BR 2493, which ends on 2016/05/23, just 3 days before the end of the data taking period for AMS nuclei
   	for (UInt_t br=2426; br<=2493; ++br) { 
	   	if (br != 2472 && br != 2473) BRs.push_back(br-FirstACEBR); 
	}

	string group[5][6] = {  { "Al", "Ar", "Ca", "Cr", "F", "Na" },
				{ "Cl", "K", "P" }, 
				{ "Co", "Mn", "Ni", "S" }, 
				{ "Mg", "Ne", "Si" }, 
				{ "Fe", "Sc", "Ti", "Va" } };

	Particle::Type group_iso[5][6] = { { Particle::ALUMINUM27, Particle::ARGON36, Particle::CALCIUM40, Particle::CHROMIUM52, Particle::FLUORINE19, Particle::SODIUM23 }, 
					   { Particle::CHLORINE35, Particle::POTASSIUM41, Particle::PHOSPHORUS31 }, 
					   { Particle::COBALT59, Particle::MANGANESE55, Particle::NICKEL60, Particle::SULFUR32 }, 
					   { Particle::MAGNESIUM24, Particle::NEON20, Particle::SILICON28 }, 
					   { Particle::IRON56, Particle::SCANDIUM45, Particle::TITANIUM46, Particle::VANADIUM51 } }; 
	
	int size[5] = { 6, 3, 4, 3, 4 }; // update every time when the line above is changed 

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

			TH1 *h_ratio = (TH1D *)HistTools::GetResiduals(h_ace, fit_comb, "_ratio", false, true, true, 4, 1);
			HistTools::SetStyle(h_ratio, HistTools::GetColorPalette(j, 4), kFullCircle, 0.9, 1, 1);

			double ratio_sum=0; // compute average of h_ratio manually  
			for(int k=0;k<14;k++){
				ratio_sum += h_ratio->GetBinContent(k);
				//printf("ratio_sum = %0.6f \n", ratio_sum);
			}
			double ratio_ave = ratio_sum/h_ratio->GetEntries();
			
			double scale = 1./ratio_ave;
			h_ratio->Scale(scale); 

			legend1->AddEntry(h_ratio, Form("%s/C", group[i][j].c_str() )); 

			gPad->SetGrid();
			h_ratio->Draw("E1X0 SAME"); 
						
			file1.Close();
		}
		//break; 
		legend1->Draw("SAME");
	} 

	c1->Print("./data/ACE/extend/ACE_extend_group_ratio_dividebyC.png");

}

UInt_t UTimeToBR(Long64_t utime){
	Long64_t first_BR = -4351622400; // Unix time corresponding to the first Bartels rotation: Feb. 8th, 1832
	return (utime - first_BR) / (27*86400) + 1;
}

UInt_t UBRToTime(int BR){
	Long64_t first_BR = -4351622400; // Unix time corresponding to the first Bartels rotation: Feb. 8th, 1832
	return (BR - 1) * (27*86400) + first_BR;
}

// over time
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

// over energy bins
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
