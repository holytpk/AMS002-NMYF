// Fill ACE data
// for i in {5..28}; do sed -r -e'1,5d' -e's/^ +//' -e's/ +/ /g' cris_energy_bands.txt | egrep "^$((i))" | cut -d' ' -f3,4 | paste -sd' ' | sed -r 's/ /, /g' >> ace_energy_band.dat; done

//CRIS: BR 2240-2529

#include "TFile.h"
#include "TTree.h"

#include "commonlib/include/Experiments.hh"
#include "commonlib/include/HistTools.hh"
#include "commonlib/include/Spline.hh"
#include "commonlib/include/DateTimeTools.hh"

UInt_t UTimeToBR(Long64_t utime);
UInt_t UBRToTime(int BR);
void ace_fill(const char *element);
void ace_convert(const char *element, Particle::Type isotope);
void ace_fluxtime();

void ace_auto(const char *operation){
 
	Experiments::DataPath = "data";

	gROOT->SetBatch();
	gROOT->ProcessLine(".L ace_auto.C");

	gSystem->mkdir("data/ACE/fill", true);
	gSystem->mkdir("data/ACE/convert", true);

	gSystem->Setenv("TZ", "UCT"); 
	gStyle->SetTimeOffset(0);
	gStyle->SetOptStat(0); 
	gStyle->SetNumberContours(99); 
	gStyle->SetPalette(55); 

	const int nAce = 28-5+1; 
	const char *elements[nAce] = { "B11", "C12", "N15", "O16", "F19", "Ne20", "Na23", "Mg24", "Al27", "Si28", "P31", "S32", "Cl35", "Ar36", "K41", "Ca40", "Sc45", "Ti46", "Va51", "Cr52", "Mn55", "Fe56", "Co59", "Ni60" };
	
	if (strcmp(operation, "fill") == 0){

		gROOT->ProcessLine(".> data/ACE/fill/fill_all.txt"); 
	
		ace_fill( "B" );
		ace_fill( "C" );
		ace_fill( "N" );
		ace_fill( "O" );
		ace_fill( "F" );
		ace_fill( "Ne" );
		ace_fill( "Na" );
		ace_fill( "Mg" );
		ace_fill( "Al" );
		ace_fill( "Si" );
		ace_fill( "P" );
		ace_fill( "S" );
		ace_fill( "Cl" );
		ace_fill( "Ar" );
		ace_fill( "K" );
		ace_fill( "Ca" );
		ace_fill( "Sc" );
		ace_fill( "Ti" );
		ace_fill( "Va" );
		ace_fill( "Cr" );
		ace_fill( "Mn" );
		ace_fill( "Fe" );
		ace_fill( "Co" );
		ace_fill( "Ni" );

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
	}

} 

//Fill CRIS Data into Bins
void ace_fill(const char *element){

	gSystem->mkdir("data/ACE/convert/fluxenergy", true);
		
	TFile *_file0 = new TFile(Form("data/ACE/%s_97226_18359.root", element));
	TTree *ace=(TTree*)_file0->Get("ace");
	
	float F[7]; // must be initialized 
	Long64_t utime;
	
	ace->SetBranchAddress("F", F);
	ace->SetBranchAddress("start_utime", &utime);	

	//ace->Print(); 

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

	double SpallCorr_B[nBins/2] = { 0.955, 0.929, 0.895, 0.863, 0.832, 0.802, 0.773 };
	double SpallCorr_C[nBins/2] = { 0.953, 0.926, 0.891, 0.857, 0.825, 0.794, 0.764 };
	double SpallCorr_N[nBins/2] = { 0.949, 0.920, 0.883, 0.847, 0.813, 0.781, 0.750 };
	double SpallCorr_O[nBins/2] = { 0.947, 0.917, 0.878, 0.842, 0.807, 0.773, 0.741 };
	double SpallCorr_F[nBins/2] = { 0.943, 0.911, 0.870, 0.832, 0.795, 0.760, 0.726 };
	double SpallCorr_Ne[nBins/2] = { 0.940, 0.908, 0.866, 0.826, 0.788, 0.752, 0.718 };
	double SpallCorr_Na[nBins/2] = { 0.938, 0.904, 0.860, 0.819, 0.780, 0.743, 0.708 };
	double SpallCorr_Mg[nBins/2] = { 0.936, 0.901, 0.857, 0.815, 0.775, 0.737, 0.701 };
	double SpallCorr_Al[nBins/2] = { 0.933, 0.897, 0.851, 0.808, 0.767, 0.728, 0.691 };
	double SpallCorr_Si[nBins/2] = { 0.932, 0.895, 0.848, 0.804, 0.763, 0.724, 0.686 };
	double SpallCorr_P[nBins/2] = { 0.929, 0.891, 0.843, 0.797, 0.754, 0.714, 0.676 };
	double SpallCorr_S[nBins/2] = { 0.928, 0.889, 0.840, 0.794, 0.751, 0.710, 0.671 };
	double SpallCorr_Cl[nBins/2] = { 0.925, 0.885, 0.834, 0.787, 0.742, 0.700, 0.660 };
	double SpallCorr_Ar[nBins/2] = { 0.922, 0.881, 0.830, 0.782, 0.737, 0.694, 0.654 };
	double SpallCorr_K[nBins/2] = { 0.920, 0.878, 0.826, 0.776, 0.730, 0.687, 0.646 };
	double SpallCorr_Ca[nBins/2] = { 0.917, 0.874, 0.821, 0.771, 0.724, 0.681, 0.639 };
	double SpallCorr_Sc[nBins/2] = { 0.915, 0.871, 0.817, 0.765, 0.718, 0.673, 0.631 };
	double SpallCorr_Ti[nBins/2] = { 0.913, 0.869, 0.813, 0.761, 0.713, 0.667, 0.625 };
	double SpallCorr_V[nBins/2] = { 0.911, 0.866, 0.809, 0.756, 0.707, 0.661, 0.618 };
	double SpallCorr_Cr[nBins/2] = { 0.910, 0.864, 0.807, 0.753, 0.703, 0.657, 0.614 };
	double SpallCorr_Mn[nBins/2] = { 0.907, 0.861, 0.802, 0.748, 0.698, 0.651, 0.607 };
	double SpallCorr_Fe[nBins/2] = { 0.906, 0.859, 0.799, 0.744, 0.693, 0.646, 0.602 };
	double SpallCorr_Co[nBins/2] = { 0.904, 0.856, 0.796, 0.740, 0.689, 0.641, 0.596 };
	double SpallCorr_Ni[nBins/2] = { 0.904, 0.855, 0.795, 0.739, 0.687, 0.639, 0.595 };

	double SpallCorrUnc_B[nBins/2] = { 0.005, 0.007, 0.011, 0.015, 0.019, 0.022, 0.026 };
	double SpallCorrUnc_C[nBins/2] = { 0.005, 0.008, 0.012, 0.016, 0.019, 0.023, 0.027 };
	double SpallCorrUnc_N[nBins/2] = { 0.005, 0.008, 0.013, 0.017, 0.021, 0.025, 0.029 };
	double SpallCorrUnc_O[nBins/2] = { 0.005, 0.009, 0.013, 0.017, 0.022, 0.026, 0.030 };
	double SpallCorrUnc_F[nBins/2] = { 0.006, 0.009, 0.014, 0.019, 0.023, 0.028, 0.033 };
	double SpallCorrUnc_Ne[nBins/2] = { 0.006, 0.010, 0.015, 0.019, 0.024, 0.029, 0.034 };
	double SpallCorrUnc_Na[nBins/2] = { 0.006, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035 };
	double SpallCorrUnc_Mg[nBins/2] = { 0.007, 0.010, 0.016, 0.021, 0.026, 0.031, 0.036 };
	double SpallCorrUnc_Al[nBins/2] = { 0.007, 0.011, 0.016, 0.022, 0.027, 0.032, 0.038 };
	double SpallCorrUnc_Si[nBins/2] = { 0.007, 0.011, 0.017, 0.022, 0.027, 0.033, 0.038 };
	double SpallCorrUnc_P[nBins/2] = { 0.007, 0.012, 0.017, 0.023, 0.029, 0.034, 0.040 };
	double SpallCorrUnc_S[nBins/2] = { 0.008, 0.012, 0.018, 0.023, 0.029, 0.035, 0.041 };
	double SpallCorrUnc_Cl[nBins/2] = { 0.008, 0.012, 0.018, 0.024, 0.030, 0.036, 0.042 };
	double SpallCorrUnc_Ar[nBins/2] = { 0.008, 0.013, 0.019, 0.025, 0.031, 0.037, 0.043 };
	double SpallCorrUnc_K[nBins/2] = { 0.008, 0.013, 0.019, 0.026, 0.032, 0.038, 0.045 };
	double SpallCorrUnc_Ca[nBins/2] = { 0.009, 0.013, 0.020, 0.026, 0.033, 0.039, 0.046 };
	double SpallCorrUnc_Sc[nBins/2] = { 0.009, 0.014, 0.020, 0.027, 0.034, 0.040, 0.047 };
	double SpallCorrUnc_Ti[nBins/2] = { 0.009, 0.014, 0.021, 0.028, 0.034, 0.041, 0.048 };
	double SpallCorrUnc_V[nBins/2] = { 0.009, 0.015, 0.021, 0.028, 0.035, 0.042, 0.049 };
	double SpallCorrUnc_Cr[nBins/2] = { 0.010, 0.015, 0.022, 0.029, 0.036, 0.043, 0.050 };
	double SpallCorrUnc_Mn[nBins/2] = { 0.010, 0.015, 0.022, 0.029, 0.037, 0.044, 0.051 };
	double SpallCorrUnc_Fe[nBins/2] = { 0.010, 0.015, 0.023, 0.030, 0.037, 0.045, 0.052 };
	double SpallCorrUnc_Co[nBins/2] = { 0.010, 0.016, 0.023, 0.031, 0.038, 0.045, 0.053 };
	double SpallCorrUnc_Ni[nBins/2] = { 0.010, 0.016, 0.023, 0.031, 0.038, 0.046, 0.053 };

	double *kin_bins, *spallcorr, *spallcorrunc;
	if (strcmp(element, "B") == 0) { 
		kin_bins = kin_bins_B;
		spallcorr = SpallCorr_B; 
		spallcorrunc = SpallCorrUnc_B; 
	} else if (strcmp(element, "C") == 0) { 
		kin_bins = kin_bins_C; 
		spallcorr = SpallCorr_C; 
		spallcorrunc = SpallCorrUnc_C;
	} else if (strcmp(element, "N") == 0) { 
		kin_bins = kin_bins_N; 
		spallcorr = SpallCorr_N; 
		spallcorrunc = SpallCorrUnc_N;
	} else if (strcmp(element, "O") == 0) { 
		kin_bins = kin_bins_O; 
		spallcorr = SpallCorr_O; 
		spallcorrunc = SpallCorrUnc_O;
	} else if (strcmp(element, "F") == 0) { 
		kin_bins = kin_bins_F; 
		spallcorr = SpallCorr_F; 
		spallcorrunc = SpallCorrUnc_F;
	} else if (strcmp(element, "Ne") == 0) { 
		kin_bins = kin_bins_Ne; 
		spallcorr = SpallCorr_Ne; 
		spallcorrunc = SpallCorrUnc_Ne;
	} else if (strcmp(element, "Na") == 0) { 
		kin_bins = kin_bins_Na; 
		spallcorr = SpallCorr_Na; 
		spallcorrunc = SpallCorrUnc_Na;
	} else if (strcmp(element, "Mg") == 0) { 
		kin_bins = kin_bins_Mg;
		spallcorr = SpallCorr_Mg; 
		spallcorrunc = SpallCorrUnc_Mg; 
	} else if (strcmp(element, "Al") == 0) { 
		kin_bins = kin_bins_Al;
		spallcorr = SpallCorr_Al; 
		spallcorrunc = SpallCorrUnc_Al; 
	} else if (strcmp(element, "Si") == 0) { 
		kin_bins = kin_bins_Si;
		spallcorr = SpallCorr_Si; 
		spallcorrunc = SpallCorrUnc_Si; 
	} else if (strcmp(element, "P") == 0) { 
		kin_bins = kin_bins_P;
		spallcorr = SpallCorr_P; 
		spallcorrunc = SpallCorrUnc_P; 
	} else if (strcmp(element, "S") == 0) { 
		kin_bins = kin_bins_S;
		spallcorr = SpallCorr_S; 
		spallcorrunc = SpallCorrUnc_S; 
	} else if (strcmp(element, "Cl") == 0) { 
		kin_bins = kin_bins_Cl;
		spallcorr = SpallCorr_Cl; 
		spallcorrunc = SpallCorrUnc_Cl; 
	} else if (strcmp(element, "Ar") == 0) { 
		kin_bins = kin_bins_Ar;
		spallcorr = SpallCorr_Ar; 
		spallcorrunc = SpallCorrUnc_Ar; 
	} else if (strcmp(element, "K") == 0) { 
		kin_bins = kin_bins_K;
		spallcorr = SpallCorr_K; 
		spallcorrunc = SpallCorrUnc_K; 
	} else if (strcmp(element, "Ca") == 0) { 
		kin_bins = kin_bins_Ca;
		spallcorr = SpallCorr_Ca; 
		spallcorrunc = SpallCorrUnc_Ca; 
	} else if (strcmp(element, "Sc") == 0) { 
		kin_bins = kin_bins_Sc;
		spallcorr = SpallCorr_Sc; 
		spallcorrunc = SpallCorrUnc_Sc; 
	} else if (strcmp(element, "Ti") == 0) { 
		kin_bins = kin_bins_Ti;
		spallcorr = SpallCorr_Ti; 
		spallcorrunc = SpallCorrUnc_Ti; 
	} else if (strcmp(element, "V") == 0) { 
		kin_bins = kin_bins_V; 
		spallcorr = SpallCorr_V; 
		spallcorrunc = SpallCorrUnc_V;
	} else if (strcmp(element, "Cr") == 0) { 
		kin_bins = kin_bins_Cr;
		spallcorr = SpallCorr_Cr; 
		spallcorrunc = SpallCorrUnc_Cr; 
	} else if (strcmp(element, "Mn") == 0) { 
		kin_bins = kin_bins_Mn;
		spallcorr = SpallCorr_Mn; 
		spallcorrunc = SpallCorrUnc_Mn; 
	} else if (strcmp(element, "Fe") == 0) { 
		kin_bins = kin_bins_Fe;
		spallcorr = SpallCorr_Fe; 
		spallcorrunc = SpallCorrUnc_Fe; 
	} else if (strcmp(element, "Co") == 0) { 
		kin_bins = kin_bins_Co;
		spallcorr = SpallCorr_Co; 
		spallcorrunc = SpallCorrUnc_Co; 
	} else if (strcmp(element, "Ni") == 0) { 
		kin_bins = kin_bins_Ni;
		spallcorr = SpallCorr_Ni; 
		spallcorrunc = SpallCorrUnc_Ni; 
	}
	
	TCanvas *c1 = new TCanvas("c1","",800,600);
	c1->cd(1);		

	TFile _file1(Form("data/ACE/fill/%s_fill.root", element), "RECREATE");

	for (int k=0; k<ace->GetEntries(); k++){

		ace->GetEntry(k); //Get the nth entry of TTree!! 

			TH1F *h = new TH1F("h","", 13, kin_bins);
			
			//ace->Scan("F", "", "col=10.4e", 1, k); 
			//ace->Scan("start_utime", "", "col=10d", 1, k);			

			for (int i=0; i<h->GetNbinsX(); ++i) { 
				if (i%2==0) { 
					h->SetBinContent(i+1, F[i/2]); 
				} 
			}
			
			for (int i=0; i<h->GetNbinsX(); ++i) { 
				printf("%s BR=%d [%02u] %.2f-%.2f   %10.4e\n", element, k, i, h->GetBinLowEdge(i+1), h->GetBinLowEdge(i+2), h->GetBinContent(i+1)); 
			}
			
			h->SetMarkerStyle(kFullCircle);
			h->SetMarkerColor(kBlue);
			h->SetMarkerSize(1.1);
			//HistTools::SetColors(h, 290, kFullCircle, 1.1);
			gPad->SetGrid(); 
			gPad->SetLogx(); 
			gPad->SetLogy();
			h->SetTitle(Form("%s BR-%d Kinetic Energy Spectrum; Energy (MeV/nuc); Flux (/(cm^2 sr s)(MeV/nuc)", element, UTimeToBR(utime)));
			h->Draw("E1X0 SAME"); 

			h->Write(Form("h_kin_%s_BR%d", element, UTimeToBR(utime)));	
	}

	c1->Print(Form("./data/ACE/convert/fluxenergy/h_kin_%s_all.png", element));
	
	_file1.Write();
	_file1.Close();

	_file0->Close();

}

// Convert CRIS Data into AMS Structure
void ace_convert(const char *element, Particle::Type isotope){

	gSystem->mkdir("data/ACE/convert/fluxtime", true);
	gSystem->mkdir("data/ACE/convert/fluxrigidity", true);	
	
	const int nBins = 14;
	
	TCanvas *c2 = new TCanvas("c2","",800,600);
	c2->cd(1);		

	TFile fin(Form("data/ACE/fill/%s_fill.root", element));
	TFile fout(Form("data/ACE/convert/%s_convert.root", element), "RECREATE");

	int utime_0, utime_1; // ams time range  

	time_t *tran_ams = Experiments::GetMeasurementTimeRange(Experiments::AMS02, 1, 0); 
		utime_0 = tran_ams[0];
		tran_ams = Experiments::GetMeasurementTimeRange(Experiments::AMS02, 1, 78);
		utime_1 = tran_ams[1];

	// convert to rigidity 
	for (int k=2240 ; k<=2529; k++){

		// utime (ace), utime_0 (ams), utime_1 (ams)
		// ace t_range = 2426 ~ 2506 
		if ( UBRToTime(k) <= utime_1 && UBRToTime(k) >= utime_0 ){

			TH1F *h = (TH1F*) fin.Get(Form("h_kin_%s_BR%d", element, k));			
			TH1 *h_rig = HistTools::TransformEnergyAndDifferentialFluxNew(h, isotope, "MeV/n cm", "GV m", Form("_rig_%s_BR%d", element, k)); // (TH1 *hist, Particle::Type particle, const Char_t *from_flux_unit, const Char_t *to_flux_unit, const Char_t *suffix)

			for (int i=0; i<h_rig->GetNbinsX(); ++i) { 
				//printf("%s BR=%d [%02u] %.2f-%.2f   %10.4e\n", element, k, i, h_rig->GetBinLowEdge(i+1), h_rig->GetBinLowEdge(i+2), h_rig->GetBinContent(i+1)); 
			}

			h_rig->SetMarkerStyle(kFullCircle);
			h_rig->SetMarkerColor(kBlue);
			h_rig->SetMarkerSize(1.1);
			//HistTools::SetColors(h_rig, 290, kFullCircle, 1.1);
			gPad->SetGrid(); 
			gPad->SetLogx(); 
			gPad->SetLogy();
			h_rig->SetTitle(Form("%s BR-%d Energy Spectrum; Rigidity (GeV); Flux (/(m^2 sr s)(GeV)", element, k));
			h_rig->Draw("E1X0 SAME"); 

			//h_rig->Write(Form("h_rig_%s_BR%d", element, k)); 
		}
			
	}

	c2->Print(Form("./data/ACE/convert/fluxrigidity/h_rig_%s_all.png", element));

	// plot flux_time
	for (int i=0; i<nBins; ++i){	
	
		if (i%2==0){	
			TCanvas *c3 = new TCanvas("c3","",800,600);
			c3->cd(1);		
			TGraph *g = new TGraph(nBins);			
			for (int k=2240; k<=2529; k++){	
				if ( UBRToTime(k) <= utime_1 && UBRToTime(k) >= utime_0 ){
					TH1 *h_rig = (TH1*) fout.Get(Form("h_rig_%s_BR%d", element, k));			
					g->SetPoint(k, UBRToTime(k), h_rig->GetBinContent(i+1)); 
					double x, y; 
					g->GetPoint(k, x, y); 
					printf("%s [%02u] BR=%d x=%0.0f y=%f \n", element, i, k, x, y); 
				}
			}

			g->GetXaxis()->SetTimeDisplay(1);
			g->GetXaxis()->SetTimeFormat("%m-%y");
			g->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00");
			g->GetXaxis()->SetTitleSize(0.7);			

			g->SetMarkerStyle(kFullCircle);
			g->SetMarkerColor(kRed);
			g->SetMarkerSize(1.1);

			//HistTools::SetColors(h, 290, kFullCircle, 1.4);
			gPad->SetGrid();  
			//gPad->SetLogy();
			g->GetXaxis()->SetRangeUser(UBRToTime(2416), UBRToTime(2516));
			g->SetTitle(Form("%s %d-th Energy Bin Flux Time Series; unix time (s); Flux (/(m^2 sr s)(GeV)", element, i));
			g->Draw("AP"); 
			c3->Print(Form("./data/ACE/convert/fluxtime/h_rig_%s_%dth.png", element, i));

		}
	}
	
	fout.Write();
	fout.Close();

	fin.Close();
}

void ace_fluxtime(){

	//gSystem->mkdir("data/ACE/fluxtime", true);
	
}

UInt_t UTimeToBR(Long64_t utime){
	Long64_t first_BR = -4351622400; // Unix time corresponding to the first Bartels rotation: Feb. 8th, 1832
	return (utime - first_BR) / (27*86400) + 1;
}

UInt_t UBRToTime(int BR){
	Long64_t first_BR = -4351622400; // Unix time corresponding to the first Bartels rotation: Feb. 8th, 1832
	return (BR - 1) * (27*86400) + first_BR;
}
