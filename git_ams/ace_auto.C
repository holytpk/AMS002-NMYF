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
void ace_fill(const char *element, Particle::Type isotope);
void ace_convert(const char *element, Particle::Type isotope);
void ace_fluxtime();

TGraphAsymmErrors *get_ace_average_graph(const char *element, UInt_t *BRs, UInt_t nBRs); 

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

	double EMed_B[nBins/2] = {59.6, 79.7, 102.0, 121.1, 138.2, 154.0, 168.6}; 
	double EMed_C[nBins/2] = {68.3, 91.5, 117.3, 139.3, 159.1, 177.4, 194.5};
	double EMed_N[nBins/2] = {73.3, 98.1, 125.9, 149.6, 171.0, 190.7, 209.2};
	double EMed_O[nBins/2] = {80.4, 107.8, 138.4, 164.7, 188.4, 210.3, 230.8};
	double EMed_F[nBins/2] = {83.5, 112.0, 143.8, 171.1, 195.9, 218.7, 240.0};
	double EMed_Ne[nBins/2] = {89.5, 120.1, 154.4, 183.9, 210.6, 235.3, 258.4};
	double EMed_Na[nBins/2] = {94.0, 126.2, 162.4, 193.5, 221.7, 247.8, 272.3};
	double EMed_Mg[nBins/2] = {100.2, 134.7, 173.4, 206.8, 237.1, 265.2, 291.5};
	double EMed_Al[nBins/2] = {103.8, 139.6, 179.8, 214.5, 246.1, 275.3, 302.8};
	double EMed_Si[nBins/2] = {110.1, 148.2, 191.1, 228.1, 261.8, 293.1, 322.6};
	double EMed_P[nBins/2] = {112.7, 151.8, 195.9, 233.9, 268.6, 300.8, 331.1};
	double EMed_S[nBins/2] = {118.2, 159.4, 205.8, 245.9, 282.5, 316.6, 348.7};
	double EMed_Cl[nBins/2] = {120.2, 162.1, 209.4, 250.3, 287.7, 322.4, 355.1};
	double EMed_Ar[nBins/2] = {125.0, 168.8, 218.1, 260.9, 300.0, 336.4, 370.8};
	double EMed_K[nBins/2] = {127.9, 172.8, 223.4, 267.4, 307.5, 344.9, 380.3};
	double EMed_Ca[nBins/2] = {131.6, 177.9, 230.1, 275.6, 317.1, 355.9, 392.4};
	double EMed_Sc[nBins/2] = {133.5, 180.5, 233.7, 279.9, 322.2, 361.6, 398.8};
	double EMed_Ti[nBins/2] = {137.1, 185.5, 240.3, 287.9, 331.6, 372.3, 410.8};
	double EMed_V[nBins/2] = {139.9, 189.5, 245.5, 294.3, 339.1, 380.8, 420.3};
	double EMed_Cr[nBins/2] = {144.0, 195.1, 253.0, 303.5, 349.8, 393.0, 434.0};
	double EMed_Mn[nBins/2] = {146.8, 199.1, 258.3, 309.9, 357.3, 401.6, 443.5};
	double EMed_Fe[nBins/2] = {150.4, 204.1, 265.0, 318.1, 366.9, 412.6, 455.9};
	double EMed_Co[nBins/2] = {153.6, 208.5, 270.9, 325.3, 375.4, 422.3, 466.7};
	double EMed_Ni[nBins/2] = {158.9, 215.9, 280.7, 337.3, 389.5, 438.4, 484.7};

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

	double *kin_bins, *SpallCorr, *SpallCorrUnc, *EMed;

void ace_auto(const char *operation){
 
	Experiments::DataPath = "data";

	//gROOT->SetBatch();
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
	}

} 

//Fill CRIS Data into Bins
void ace_fill(const char *element, Particle::Type isotope){

	if (strcmp(element, "B") == 0) { 
		kin_bins = kin_bins_B;
		SpallCorr = SpallCorr_B; 
		SpallCorrUnc = SpallCorrUnc_B; 
		EMed = EMed_B;
	} else if (strcmp(element, "C") == 0) { 
		kin_bins = kin_bins_C; 
		SpallCorr = SpallCorr_C; 
		SpallCorrUnc = SpallCorrUnc_C;
		EMed = EMed_C;
	} else if (strcmp(element, "N") == 0) { 
		kin_bins = kin_bins_N; 
		SpallCorr = SpallCorr_N; 
		SpallCorrUnc = SpallCorrUnc_N;
		EMed = EMed_N;
	} else if (strcmp(element, "O") == 0) { 
		kin_bins = kin_bins_O; 
		SpallCorr = SpallCorr_O; 
		SpallCorrUnc = SpallCorrUnc_O;
		EMed = EMed_O;
	} else if (strcmp(element, "F") == 0) { 
		kin_bins = kin_bins_F; 
		SpallCorr = SpallCorr_F; 
		SpallCorrUnc = SpallCorrUnc_F;
		EMed = EMed_F;
	} else if (strcmp(element, "Ne") == 0) { 
		kin_bins = kin_bins_Ne; 
		SpallCorr = SpallCorr_Ne; 
		SpallCorrUnc = SpallCorrUnc_Ne;
		EMed = EMed_Ne;
	} else if (strcmp(element, "Na") == 0) { 
		kin_bins = kin_bins_Na; 
		SpallCorr = SpallCorr_Na; 
		SpallCorrUnc = SpallCorrUnc_Na;
		EMed = EMed_Na;
	} else if (strcmp(element, "Mg") == 0) { 
		kin_bins = kin_bins_Mg;
		SpallCorr = SpallCorr_Mg; 
		SpallCorrUnc = SpallCorrUnc_Mg; 
		EMed = EMed_Mg;
	} else if (strcmp(element, "Al") == 0) { 
		kin_bins = kin_bins_Al;
		SpallCorr = SpallCorr_Al; 
		SpallCorrUnc = SpallCorrUnc_Al; 
		EMed = EMed_Al;
	} else if (strcmp(element, "Si") == 0) { 
		kin_bins = kin_bins_Si;
		SpallCorr = SpallCorr_Si; 
		SpallCorrUnc = SpallCorrUnc_Si; 
		EMed = EMed_Si;
	} else if (strcmp(element, "P") == 0) { 
		kin_bins = kin_bins_P;
		SpallCorr = SpallCorr_P; 
		SpallCorrUnc = SpallCorrUnc_P;
		EMed = EMed_P; 
	} else if (strcmp(element, "S") == 0) { 
		kin_bins = kin_bins_S;
		SpallCorr = SpallCorr_S; 
		SpallCorrUnc = SpallCorrUnc_S; 
		EMed = EMed_S;
	} else if (strcmp(element, "Cl") == 0) { 
		kin_bins = kin_bins_Cl;
		SpallCorr = SpallCorr_Cl; 
		SpallCorrUnc = SpallCorrUnc_Cl; 
		EMed = EMed_Cl;
	} else if (strcmp(element, "Ar") == 0) { 
		kin_bins = kin_bins_Ar;
		SpallCorr = SpallCorr_Ar; 
		SpallCorrUnc = SpallCorrUnc_Ar; 
		EMed = EMed_Ar;
	} else if (strcmp(element, "K") == 0) { 
		kin_bins = kin_bins_K;
		SpallCorr = SpallCorr_K; 
		SpallCorrUnc = SpallCorrUnc_K;
		EMed = EMed_K; 
	} else if (strcmp(element, "Ca") == 0) { 
		kin_bins = kin_bins_Ca;
		SpallCorr = SpallCorr_Ca; 
		SpallCorrUnc = SpallCorrUnc_Ca;
		EMed = EMed_Ca; 
	} else if (strcmp(element, "Sc") == 0) { 
		kin_bins = kin_bins_Sc;
		SpallCorr = SpallCorr_Sc; 
		SpallCorrUnc = SpallCorrUnc_Sc; 
		EMed = EMed_Sc;
	} else if (strcmp(element, "Ti") == 0) { 
		kin_bins = kin_bins_Ti;
		SpallCorr = SpallCorr_Ti; 
		SpallCorrUnc = SpallCorrUnc_Ti;
		EMed = EMed_Ti; 
	} else if (strcmp(element, "V") == 0) { 
		kin_bins = kin_bins_V; 
		SpallCorr = SpallCorr_V; 
		SpallCorrUnc = SpallCorrUnc_V;
		EMed = EMed_V;
	} else if (strcmp(element, "Cr") == 0) { 
		kin_bins = kin_bins_Cr;
		SpallCorr = SpallCorr_Cr; 
		SpallCorrUnc = SpallCorrUnc_Cr; 
		EMed = EMed_Cr;
	} else if (strcmp(element, "Mn") == 0) { 
		kin_bins = kin_bins_Mn;
		SpallCorr = SpallCorr_Mn; 
		SpallCorrUnc = SpallCorrUnc_Mn; 
		EMed = EMed_Mn;
	} else if (strcmp(element, "Fe") == 0) { 
		kin_bins = kin_bins_Fe;
		SpallCorr = SpallCorr_Fe; 
		SpallCorrUnc = SpallCorrUnc_Fe; 
		EMed = EMed_Fe;
	} else if (strcmp(element, "Co") == 0) { 
		kin_bins = kin_bins_Co;
		SpallCorr = SpallCorr_Co; 
		SpallCorrUnc = SpallCorrUnc_Co; 
		EMed = EMed_Co;
	} else if (strcmp(element, "Ni") == 0) { 
		kin_bins = kin_bins_Ni;
		SpallCorr = SpallCorr_Ni; 
		SpallCorrUnc = SpallCorrUnc_Ni;
		EMed = EMed_Ni; 
	}

	gSystem->mkdir("data/ACE/convert/fluxenergy", true);
		
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
		h->Draw("PSAME"); 

		h->Write(Form("h_kin_%s_BR%d", element, UTimeToBR(utime)));
	}

	file1.Write();
	file1.Close();

	TCanvas *c1 = new TCanvas("c1","",800,600);
	c1->Divide(2, 1);

	const UInt_t FirstACEBR = 2240;
	vector<UInt_t> BRs;
	// we stop at BR 2493, which ends on 2016/05/23, just 3 days before the end of the data taking period for AMS nuclei
	for (UInt_t br=2426; br<=2493; ++br) { 
		if (br != 2472 && br != 2473) BRs.push_back(br-FirstACEBR); 
	}
	TH1D *h_ene = HistTools::GraphToHist(get_ace_average_graph(element, &BRs[0], BRs.size()));
	h_ene->SetXTitle(Unit::GetEnergyLabel("MeV/n"));
	h_ene->SetYTitle(Unit::GetDifferentialFluxLabel("MeV/n cm"));
	TH1 *h_rig = HistTools::TransformEnergyAndDifferentialFluxNew(h_ene, isotope, "MeV/n cm", "GV m", "_rig");
	h_rig->SetXTitle(Unit::GetEnergyLabel("GV"));
	h_rig->SetYTitle(Unit::GetDifferentialFluxLabel("GV m"));			
			
	h_ene->SetMarkerStyle(kFullCircle);
	h_ene->SetMarkerColor(kBlue);
	h_ene->SetMarkerSize(1.1);
	//HistTools::SetColors(h, 290, kFullCircle, 1.1);
	h_ene->SetTitle(Form("%s BR-averaged Kinetic Energy Spectrum", element));
	c1->cd(1);
	gPad->SetGrid(); 
	h_ene->Draw("E1X0"); 

	h_rig->SetMarkerStyle(kFullCircle);
	h_rig->SetMarkerColor(kBlue);
	h_rig->SetMarkerSize(1.1);
	//HistTools::SetColors(h, 290, kFullCircle, 1.1); 
	h_rig->SetTitle(Form("%s BR-averaged Energy Spectrum", element));
	c1->cd(2);
	gPad->SetGrid();
	h_rig->Draw("E1X0"); 

	c1->Print(Form("./data/ACE/convert/fluxenergy/h_kin_%s_averaged.png", element));

	TFile file2(Form("data/ACE/fill/%s_averaged.root", element), "RECREATE");

	h_ene->Write(Form("h_kin_%s_BR_averaged", element));	
	h_rig->Write(Form("h_rig_%s_BR_averaged", element));
	
	file2.Write();
	file2.Close();

	_file0->Close();

}

// Convert CRIS Data into AMS Structure
void ace_convert(const char *element, Particle::Type isotope){

	if (strcmp(element, "B") == 0) { 
		kin_bins = kin_bins_B;
		SpallCorr = SpallCorr_B; 
		SpallCorrUnc = SpallCorrUnc_B; 
		EMed = EMed_B;
	} else if (strcmp(element, "C") == 0) { 
		kin_bins = kin_bins_C; 
		SpallCorr = SpallCorr_C; 
		SpallCorrUnc = SpallCorrUnc_C;
		EMed = EMed_C;
	} else if (strcmp(element, "N") == 0) { 
		kin_bins = kin_bins_N; 
		SpallCorr = SpallCorr_N; 
		SpallCorrUnc = SpallCorrUnc_N;
		EMed = EMed_N;
	} else if (strcmp(element, "O") == 0) { 
		kin_bins = kin_bins_O; 
		SpallCorr = SpallCorr_O; 
		SpallCorrUnc = SpallCorrUnc_O;
		EMed = EMed_O;
	} else if (strcmp(element, "F") == 0) { 
		kin_bins = kin_bins_F; 
		SpallCorr = SpallCorr_F; 
		SpallCorrUnc = SpallCorrUnc_F;
		EMed = EMed_F;
	} else if (strcmp(element, "Ne") == 0) { 
		kin_bins = kin_bins_Ne; 
		SpallCorr = SpallCorr_Ne; 
		SpallCorrUnc = SpallCorrUnc_Ne;
		EMed = EMed_Ne;
	} else if (strcmp(element, "Na") == 0) { 
		kin_bins = kin_bins_Na; 
		SpallCorr = SpallCorr_Na; 
		SpallCorrUnc = SpallCorrUnc_Na;
		EMed = EMed_Na;
	} else if (strcmp(element, "Mg") == 0) { 
		kin_bins = kin_bins_Mg;
		SpallCorr = SpallCorr_Mg; 
		SpallCorrUnc = SpallCorrUnc_Mg; 
		EMed = EMed_Mg;
	} else if (strcmp(element, "Al") == 0) { 
		kin_bins = kin_bins_Al;
		SpallCorr = SpallCorr_Al; 
		SpallCorrUnc = SpallCorrUnc_Al; 
		EMed = EMed_Al;
	} else if (strcmp(element, "Si") == 0) { 
		kin_bins = kin_bins_Si;
		SpallCorr = SpallCorr_Si; 
		SpallCorrUnc = SpallCorrUnc_Si; 
		EMed = EMed_Si;
	} else if (strcmp(element, "P") == 0) { 
		kin_bins = kin_bins_P;
		SpallCorr = SpallCorr_P; 
		SpallCorrUnc = SpallCorrUnc_P;
		EMed = EMed_P; 
	} else if (strcmp(element, "S") == 0) { 
		kin_bins = kin_bins_S;
		SpallCorr = SpallCorr_S; 
		SpallCorrUnc = SpallCorrUnc_S; 
		EMed = EMed_S;
	} else if (strcmp(element, "Cl") == 0) { 
		kin_bins = kin_bins_Cl;
		SpallCorr = SpallCorr_Cl; 
		SpallCorrUnc = SpallCorrUnc_Cl; 
		EMed = EMed_Cl;
	} else if (strcmp(element, "Ar") == 0) { 
		kin_bins = kin_bins_Ar;
		SpallCorr = SpallCorr_Ar; 
		SpallCorrUnc = SpallCorrUnc_Ar; 
		EMed = EMed_Ar;
	} else if (strcmp(element, "K") == 0) { 
		kin_bins = kin_bins_K;
		SpallCorr = SpallCorr_K; 
		SpallCorrUnc = SpallCorrUnc_K;
		EMed = EMed_K; 
	} else if (strcmp(element, "Ca") == 0) { 
		kin_bins = kin_bins_Ca;
		SpallCorr = SpallCorr_Ca; 
		SpallCorrUnc = SpallCorrUnc_Ca;
		EMed = EMed_Ca; 
	} else if (strcmp(element, "Sc") == 0) { 
		kin_bins = kin_bins_Sc;
		SpallCorr = SpallCorr_Sc; 
		SpallCorrUnc = SpallCorrUnc_Sc; 
		EMed = EMed_Sc;
	} else if (strcmp(element, "Ti") == 0) { 
		kin_bins = kin_bins_Ti;
		SpallCorr = SpallCorr_Ti; 
		SpallCorrUnc = SpallCorrUnc_Ti;
		EMed = EMed_Ti; 
	} else if (strcmp(element, "V") == 0) { 
		kin_bins = kin_bins_V; 
		SpallCorr = SpallCorr_V; 
		SpallCorrUnc = SpallCorrUnc_V;
		EMed = EMed_V;
	} else if (strcmp(element, "Cr") == 0) { 
		kin_bins = kin_bins_Cr;
		SpallCorr = SpallCorr_Cr; 
		SpallCorrUnc = SpallCorrUnc_Cr; 
		EMed = EMed_Cr;
	} else if (strcmp(element, "Mn") == 0) { 
		kin_bins = kin_bins_Mn;
		SpallCorr = SpallCorr_Mn; 
		SpallCorrUnc = SpallCorrUnc_Mn; 
		EMed = EMed_Mn;
	} else if (strcmp(element, "Fe") == 0) { 
		kin_bins = kin_bins_Fe;
		SpallCorr = SpallCorr_Fe; 
		SpallCorrUnc = SpallCorrUnc_Fe; 
		EMed = EMed_Fe;
	} else if (strcmp(element, "Co") == 0) { 
		kin_bins = kin_bins_Co;
		SpallCorr = SpallCorr_Co; 
		SpallCorrUnc = SpallCorrUnc_Co; 
		EMed = EMed_Co;
	} else if (strcmp(element, "Ni") == 0) { 
		kin_bins = kin_bins_Ni;
		SpallCorr = SpallCorr_Ni; 
		SpallCorrUnc = SpallCorrUnc_Ni;
		EMed = EMed_Ni; 
	}

	gSystem->mkdir("data/ACE/convert/fluxtime", true);
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

	// plot averaged flux (MeV/n cm) vs. kinetic energy within AMS range
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
				printf("%s BR=%d [%02u] %.2f-%.2f y=%10.4e dy=%10.4e  \n", element, k, i, h_rig->GetBinLowEdge(i+1), h_rig->GetBinLowEdge(i+2), h_rig->GetBinContent(i+1), h_rig->GetBinError(i+1)); 
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
	TLegend *legend = new TLegend(0.1,0.7,0.28,0.9); // left, down, right, top 
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
			g->SetMarkerColor(kRed-3*i);
			g->SetLineColor(kRed-3*i);
			g->SetLineWidth(1);
			g->SetMarkerSize(0.7);

			double x, y;	
			g->GetPoint(0, x, y);

   			//double *yaxis = g->GetY();
  		 	//double gmax = yaxis[TMath::LocMax(g->GetN(), yaxis);
			gPad->SetGrid();  
			//gPad->SetLogy();
			g->GetYaxis()->SetRangeUser(y-y*0.99, y+y*2.5);
			g->GetXaxis()->SetRangeUser(UBRToTime(2425), UBRToTime(2507));
			g->SetTitle(Form("%s All Energy Bin Flux Time Series", element));
			//g->Print();
			c3->cd(1);
			
			legend->AddEntry(g, Form("%d-th Bin Flux", i/2+1), "l");

			if (i/2==0){
				g->Draw("ALP");
				legend->Draw("SAME");
			} else {
				g->Draw("LPSAME");
				legend->Draw("SAME");
			}
		}
	} 

	c3->Print(Form("./data/ACE/convert/fluxtime/h_fluxtime_%s_all.png", element));	
	
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

TGraphAsymmErrors *get_ace_average_graph(const char *element, UInt_t *BRs, UInt_t nBRs){
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
