// Fill ACE data, windows backup version 
// Extract energy band info: sed -r -e'1,5d' -e's/^ +//' -e's/ +/ /g' cris_energy_bands.txt | egrep "^6" | cut -d' ' -f3,4 | paste -sd' ' | sed -r 's/ /, /g'
// for i in {5..28}; do sed -r -e'1,5d' -e's/^ +//' -e's/ +/ /g' cris_energy_bands.txt | egrep "^$((i))" | cut -d' ' -f3,4 | paste -sd' ' | sed -r 's/ /, /g' >> ace_energy_band.dat; done
// Energy band txt file saved as ace_energy_band.dat

// - replace egrep "^6" to obtain the desired element data  

#include "TFile.h"
#include "TTree.h"

#include "commonlib/include/Experiments.hh"
#include "commonlib/include/HistTools.hh"
#include "commonlib/include/Spline.hh"
#include "commonlib/include/DateTimeTools.hh"

UInt_t UTimeToBR(Long64_t utime);
void ace_convert(const char *element, Particle::Type isotope);

void aceauto_win(){

	gROOT->ProcessLine(".L aceauto.C");

	const int nAce = 28-5+1; 
	const char *elements[nAce] = { "B11", "C12", "N15", "O16", "F19", "Ne20", "Na23", "Mg24", "Al27", "Si28", "P31", "S32", "Cl35", "Ar36", "K41", "Ca40", "Sc45", "Ti46", "Va51", "Cr52", "Mn55", "Fe56", "Co59", "Ni60" };
	
	acefill( "B", Particle::BORON11 );
	acefill( "C", Particle::CARBON12 );
	acefill( "N", Particle::NITROGEN15 );
	acefill( "O", Particle::OXYGEN16 );
	acefill( "F", Particle::FLUORINE19 );
	acefill( "Ne", Particle::NEON20 );
	acefill( "Na", Particle::SODIUM23 );
	acefill( "Mg", Particle::MAGNESIUM24 );
	acefill( "Al", Particle::ALUMINUM27 );
	acefill( "Si", Particle::SILICON28 );
	acefill( "P", Particle::PHOSPHORUS31 );
	acefill( "S", Particle::SULFUR32 );
	acefill( "Cl", Particle::CHLORINE35 );
	acefill( "Ar", Particle::ARGON36 );
	acefill( "K", Particle::POTASSIUM41 );
	acefill( "Ca", Particle::CALCIUM40 );
	acefill( "Sc", Particle::SCANDIUM45 );
	acefill( "Ti", Particle::TITANIUM46 );
	acefill( "Va", Particle::VANADIUM51 );
	acefill( "Cr", Particle::CHROMIUM52 );
	acefill( "Mn", Particle::MANGANESE55 );
	acefill( "Fe", Particle::IRON56 );
	acefill( "Co", Particle::COBALT59 );
	acefill( "Ni", Particle::NICKEL60 );
	
	
		
} 

//Fill & Convert CRIS Data into AMS Structure
void ace_convert(const char *element, Particle::Type isotope){

	gROOT->SetBatch(); 
	
	TFile *_file0 = new TFile(Form("data/ACE/%s_97226_18359.root", element));
	TTree *ace=(TTree*)_file0->Get("ace");
	
	ace->SetBranchAddress("F", F);
	ace->SetBranchAddress("start_utime", &utime);
	
	float F[7]; // must be initialized 
	Long64_t utime;

	ace->Print(); 

	gSystem->Setenv("TZ", "UCT"); 
	gStyle->SetTimeOffset(0);
	gStyle->SetOptStat(0); 
	gStyle->SetNumberContours(99); 
	gStyle->SetPalette(55); 

	//ace->Draw("F[0]:start_utime", "", "L");
	//ace->Draw("F:Iteration$", "", "L same", 1, 1);

	//TH2F *h2 = new TH2F("h2", Form("%s Flux - ACE/CRIS;Time;Energy bin", element), 290, 871552800,1548064800,7,0,7);
	
	//ace->Draw("Iteration$:start_utime>>h2", "F", "zcol");
	//ace->Draw("Iteration$:start_utime>>h2", "C", "zcol");
	//h2->Draw("ZCOL"); 

	const int nBins = 14;
	double kin_bins_B[nBins] = { 51.4, 65.8, 67.3, 90.4, 91.5, 110.6, 111.6, 128.4, 129.2, 144.6, 145.3, 159.4, 160.2, 173.7 };
	double kin_bins_C[nBins] = { 59.0, 75.6, 77.2, 103.8, 105.1, 127.3, 128.3, 147.9, 148.8, 166.6, 167.4, 183.9, 184.8, 200.4 };
	double kin_bins_N[nBins] = { 63.2, 81.0, 82.8, 111.4, 112.7, 136.6, 137.7, 158.8, 159.8, 179.0, 179.9, 197.6, 198.6, 215.5 };
	double kin_bins_O[nBins] = { 69.4, 89.0, 91.0, 122.5, 124.0, 150.3, 151.6, 174.9, 176.0, 197.3, 198.3, 218.0, 219.1, 237.9 };
	double kin_bins_F[nBins] = { 72.0, 92.4, 94.4, 127.2, 128.8, 156.2, 157.5, 181.8, 182.9, 205.1, 206.2, 226.8, 227.9, 247.5 };
	double kin_bins_Ne[nBins] = { 77.1, 99.0, 101.2, 136.4, 138.1, 167.6, 169.1, 195.2, 196.5, 220.4, 221.6, 243.8, 245.0, 266.3 };
	double kin_bins_Na[nBins] = { 81.0, 104.1, 106.4, 143.5, 145.3, 176.5, 178.0, 205.6, 206.9, 232.3, 233.5, 257.0, 258.4, 280.9 };
	double kin_bins_Mg[nBins] = { 86.3, 110.9, 113.4, 153.1, 155.0, 188.4, 190.1, 219.7, 221.1, 248.4, 249.7, 275.0, 276.4, 300.7 };
	double kin_bins_Al[nBins] = { 89.4, 115.0, 117.5, 158.8, 160.8, 195.5, 197.2, 228.1, 229.5, 257.9, 259.2, 285.6, 287.1, 312.4 };
	double kin_bins_Si[nBins] = { 94.8, 121.9, 124.7, 168.6, 170.7, 207.8, 209.6, 242.6, 244.1, 274.5, 275.9, 304.2, 305.7, 332.9 };
	double kin_bins_P[nBins] = { 97.1, 124.9, 127.7, 172.8, 175.0, 213.0, 214.9, 248.8, 250.4, 281.6, 283.1, 312.2, 313.8, 341.7 };
	double kin_bins_S[nBins] = { 101.6, 130.8, 133.7, 181.2, 183.4, 223.5, 225.4, 261.2, 262.8, 295.8, 297.3, 328.1, 329.8, 359.3 };
	double kin_bins_Cl[nBins] = { 103.2, 133.0, 136.0, 184.3, 186.6, 227.4, 229.4, 265.8, 267.5, 301.1, 302.7, 334.0, 335.8, 365.9 };
	double kin_bins_Ar[nBins] = { 107.7, 138.8, 141.9, 192.4, 194.9, 237.7, 239.7, 278.0, 279.8, 315.1, 316.7, 349.7, 351.6, 383.3 };
	double kin_bins_K[nBins] = { 110.0, 141.8, 145.0, 196.8, 199.3, 243.1, 245.3, 284.5, 286.3, 322.5, 324.2, 358.1, 360.0, 392.6 };
	double kin_bins_Ca[nBins] = { 113.3, 146.2, 149.5, 203.1, 205.6, 251.0, 253.2, 293.9, 295.8, 333.3, 335.1, 370.2, 372.2, 406.0 };
	double kin_bins_Sc[nBins] = { 114.7, 148.1, 151.4, 205.8, 208.4, 254.4, 256.7, 297.9, 299.9, 338.0, 339.8, 375.5, 377.5, 411.9 };
	double kin_bins_Ti[nBins] = { 117.8, 152.1, 155.5, 211.5, 214.1, 261.6, 263.9, 306.5, 308.5, 347.9, 349.7, 386.6, 388.7, 424.2 };
	double kin_bins_V[nBins] = { 120.2, 155.3, 158.8, 216.0, 218.8, 267.4, 269.8, 313.4, 315.4, 355.8, 357.7, 395.5, 397.6, 434.1 };
	double kin_bins_Cr[nBins] = { 123.4, 159.6, 163.2, 222.2, 225.0, 275.2, 277.6, 322.6, 324.8, 366.5, 368.4, 407.6, 409.8, 447.5 };
	double kin_bins_Mn[nBins] = { 126.0, 163.0, 166.7, 227.0, 230.0, 281.3, 283.9, 330.0, 332.2, 375.0, 377.0, 417.2, 419.4, 458.2 };
	double kin_bins_Fe[nBins] = { 129.1, 167.0, 170.8, 232.9, 235.9, 288.7, 291.3, 338.9, 341.1, 385.2, 387.3, 428.7, 431.0, 471.0 };
	double kin_bins_Co[nBins] = { 131.7, 170.6, 174.5, 238.0, 241.0, 295.2, 297.9, 346.6, 348.9, 394.1, 396.3, 438.8, 441.1, 482.2 };
	double kin_bins_Ni[nBins] = { 136.2, 176.5, 180.5, 246.4, 249.6, 305.9, 308.7, 359.4, 361.8, 408.9, 411.2, 455.5, 458.0, 500.8 };

	double *kin_bins;
	if (strcmp(element, "B") == 0) { 
		kin_bins = kin_bins_B; 
	} else if (strcmp(element, "C") == 0) { 
		kin_bins = kin_bins_C; 
	} else if (strcmp(element, "N") == 0) { 
		kin_bins = kin_bins_N; 
	} else if (strcmp(element, "O") == 0) { 
		kin_bins = kin_bins_O; 
	} else if (strcmp(element, "F") == 0) { 
		kin_bins = kin_bins_F; 
	} else if (strcmp(element, "Ne") == 0) { 
		kin_bins = kin_bins_Ne; 
	} else if (strcmp(element, "Na") == 0) { 
		kin_bins = kin_bins_Na; 
	} else if (strcmp(element, "Mg") == 0) { 
		kin_bins = kin_bins_Mg; 
	} else if (strcmp(element, "Al") == 0) { 
		kin_bins = kin_bins_Al; 
	} else if (strcmp(element, "Si") == 0) { 
		kin_bins = kin_bins_Si; 
	} else if (strcmp(element, "P") == 0) { 
		kin_bins = kin_bins_P; 
	} else if (strcmp(element, "S") == 0) { 
		kin_bins = kin_bins_S; 
	} else if (strcmp(element, "Cl") == 0) { 
		kin_bins = kin_bins_Cl; 
	} else if (strcmp(element, "Ar") == 0) { 
		kin_bins = kin_bins_Ar; 
	} else if (strcmp(element, "K") == 0) { 
		kin_bins = kin_bins_K; 
	} else if (strcmp(element, "Ca") == 0) { 
		kin_bins = kin_bins_Ca; 
	} else if (strcmp(element, "Sc") == 0) { 
		kin_bins = kin_bins_Sc; 
	} else if (strcmp(element, "Ti") == 0) { 
		kin_bins = kin_bins_Ti; 
	} else if (strcmp(element, "V") == 0) { 
		kin_bins = kin_bins_V; 
	} else if (strcmp(element, "Cr") == 0) { 
		kin_bins = kin_bins_Cr; 
	} else if (strcmp(element, "Mn") == 0) { 
		kin_bins = kin_bins_Mn; 
	} else if (strcmp(element, "Fe") == 0) { 
		kin_bins = kin_bins_Fe; 
	} else if (strcmp(element, "Co") == 0) { 
		kin_bins = kin_bins_Co; 
	} else if (strcmp(element, "Ni") == 0) { 
		kin_bins = kin_bins_Ni; 
	}
	
	TCanvas *c1 = new TCanvas("c1","",1800,1200);
	c1->cd(1);		

	gSystem->mkdir("data/ACE/fill", true);
	TFile _file1(Form("data/ACE/fill/%s_fill.root", element), "RECREATE");

	for (int k=0; k<ace->GetEntries(); k++){

			TH1F *h = new TH1F("h","", 13, kin_bins);
			
			ace->GetEntry(k); //Get the nth entry of TTree!! 			
			ace->Scan("F", "", "col=10.4e", 1, k); 

			ace->GetEntry(k);
			ace->Scan("start_utime", "", "col=10d", 1, k);			

			for (int i=0; i<h->GetNbinsX(); ++i) { 
				if (i%2==0) { 
					h->SetBinContent(i+1, F[i/2]); 
				} 
			}
			
			for (int i=0; i<h->GetNbinsX(); ++i) { 
				printf("[%02u] %.2f-%.2f   %10.4e\n", i, h->GetBinLowEdge(i+1), h->GetBinLowEdge(i+2), h->GetBinContent(i+1)); 
			}
			
			TH1 *h_rig = HistTools::TransformEnergyAndDifferentialFluxNew(h, isotope, "MeV/n cm", "GV m", Form("_rig_%s_BR%d", element, UTimeToBR(utime))); // (TH1 *hist, Particle::Type particle, const Char_t *from_flux_unit, const Char_t *to_flux_unit, const Char_t *suffix)
			
			h->SetMarkerStyle(kFullCircle);
			//HistTools::SetColors(h, 290, kFullCircle, 1.4);
			gPad->SetGrid(); 
			gPad->SetLogx(); 
			gPad->SetLogy();
			h->SetTitle(Form("%s BR-%d Energy Spectrum; Energy (MeV/nuc); Flux (/(cm^2 sr s)(MeV/nuc)", element, UTimeToBR(utime)));
			h_rig->SetTitle(Form("%s BR-%d Energy Spectrum; Rigidity (GeV); Flux (/(m^2 sr s)(GeV)", element, UTimeToBR(utime)));
			h->Draw("PSAME"); 

			h->Write(Form("h_kin_%s_BR%d", element, UTimeToBR(utime)));
	}
	
	_file1.Write();
	_file1.Close();

	_file0->Close();

}

UInt_t UTimeToBR(Long64_t utime)
{
    const first_BR = -4351622400; // Unix time corresponding to the 
first Bartels rotation: Feb. 8th, 1832
    return (utime - first_BR) / (27*86400) + 1;
}
