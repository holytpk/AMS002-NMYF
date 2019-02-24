// Fill ACE data
// Extract energy band info: sed -r -e'1,5d' -e's/^ +//' -e's/ +/ /g' cris_energy_bands.txt | egrep "^6" | cut -d' ' -f3,4 | paste -sd' ' | sed -r 's/ /, /g'

#include "TFile.h"
#include "TTree.h"

void aceauto(const char *element){

	//TFile *file = TFile::Open("data/ACE/C_97226_18359.root");

	float F[7];

	//TTreeReader *ace = new TTreeReader("ace", &file); 

	//gROOT->ProcessLine("ace->Print();");	
	ace->Print(); 

	gSystem->Setenv("TZ", "UCT"); 
	gStyle->SetTimeOffset(0);
	gStyle->SetOptStat(0); 
	gStyle->SetNumberContours(99); 
	gStyle->SetPalette(55); 

	ace->Draw("F[0]:start_utime", "", "L");
	ace->Draw("F:Iteration$", "", "L same", 1, 1);

	TH2F *h2 = new TH2F("h2", Form("%s Flux - ACE/CRIS;Time;Energy bin", element), 290, 871552800,1548064800,7,0,7);
	
	ace->Draw("Iteration$:start_utime>>h2", "F", "zcol");
	//ace->Draw("Iteration$:start_utime>>h2", "C", "zcol");
	//h2->Draw("ZCOL"); 

	double kin_bins[] = { 59.0, 75.6, 77.2, 103.8, 105.1, 127.3, 128.3, 147.9, 148.8, 166.6, 167.4, 183.9, 184.8, 200.4 };
	TH1F *h = new TH1F("h", "", 13, kin_bins);  

	ace->SetBranchAddress("F", F);
	ace->GetEntry(0);
	//F[0]; 
	ace->Scan("F", "", "col=10.4e", 1, 0);

	for (int i=0; i<h->GetNbinsX(); ++i) { 
		if (i%2==0) { 
			h->SetBinContent(i+1, F[i/2]); 
		} 
	}
	for (int i=0; i<h->GetNbinsX(); ++i) { 
		printf("[%02u] %.2f-%.2f   %10.4e\n", i, h->GetBinLowEdge(i+1), h->GetBinLowEdge(i+2), h->GetBinContent(i+1)); 
	}

	ace->Scan("F", "", "col=10.4e", 1, 0); 

	//file->Close();

}
