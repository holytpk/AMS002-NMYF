// Fit AMS p, He Fluxes by Lingqiang He and Claudio Corti, 01/20/2019
// Compare statistical fitting goodness w/ diferent node #s 

// To check dataset information, use "Experiments::Info[Experiments::AMS02].Print();"

#include "commonlib/include/Experiments.hh"
#include "commonlib/include/HistTools.hh"
#include "commonlib/include/Spline.hh"

void node(){

	gROOT->SetBatch(); // disable canvas loop	

	gStyle->SetOptStat(0);
	Experiments::DataPath = "data";
	// Experiments::Info[Experiments::AMS02].Print();

	// Define variables
	double sigma = 0, sigma2 = 0;
	int iBR = 0; // BR index
	
	double iBR_D = (double)iBR;
	int nnodes = 0; 
	int nnodes1 = 1;
	int nnodes2 = 2;

	// Load the AMS monthly p and He fluxes:
	const int nBRs = Experiments::Info[Experiments::AMS02].Dataset[1].nMeasurements; // read the number of BRs
	//int iBR=0;
	TH1D **h_BR_p = Experiments::GetDatasetHistograms(Experiments::AMS02, 1);
	TH1D **h_BR_he = Experiments::GetDatasetHistograms(Experiments::AMS02, 4);

	// Load the AMS integrated p and He fluxes: 
	TH1D *h_p = Experiments::GetMeasurementHistogram(Experiments::AMS02, 0, 0);
	TH1D *h_he = Experiments::GetMeasurementHistogram(Experiments::AMS02, 18, 0);

	// Create a new set of histograms which uses the monthly fluxes at low rigidities and the integrated fluxes at high rigidities:
	TH1D **h_BR_p_all = new TH1D *[nBRs]; // this creates a new array of pointers with nBRs items
	TH1D **h_BR_he_all = new TH1D *[nBRs]; // this creates a new array of pointers with nBRs items	

	// Create a Sigma test array
	TGraph *h_p_sig = new TGraph(nBRs); 
	TGraph *h_he_sig = new TGraph(nBRs); 

	// Open TFile
	TFile file_fitresult("fit_result.root", "RECREATE");

   for (nnodes1=5; nnodes1<=9; nnodes1++)
   {
   
	for(nnodes2=nnodes1+1; nnodes2<=9; nnodes2++)
	{
		for (iBR = 0; iBR < nBRs; iBR++)
		{
    			h_BR_p_all[iBR] = (TH1D *)h_p->Clone(Form("%s_all", h_BR_p[iBR]->GetName())); // set the histogram style as you prefer, with 
			HistTools::SetStyle(h_BR_p_all[iBR], kBlue, kFullCircle, 1.1, 1, 1);

			h_BR_he_all[iBR] = (TH1D *)h_he->Clone(Form("%s_all", h_BR_he[iBR]->GetName())); // set the histogram style as you prefer, with 
			HistTools::SetStyle(h_BR_he_all[iBR], kBlue, kFullCircle, 1.1, 1, 1);

    			for (int bin = 1; bin <= h_BR_p[0]->GetNbinsX(); ++bin)
    			{
       				h_BR_p_all[iBR]->SetBinContent(bin, h_BR_p[iBR]->GetBinContent(bin));
       				h_BR_p_all[iBR]->SetBinError(bin, h_BR_p[iBR]->GetBinError(bin));

    			}

			for (int bin = 1; bin <= h_BR_he[0]->GetNbinsX(); ++bin)
    			{
       				h_BR_he_all[iBR]->SetBinContent(bin, h_BR_he[iBR]->GetBinContent(bin));
       				h_BR_he_all[iBR]->SetBinError(bin, h_BR_he[iBR]->GetBinError(bin));

    			}

			// Fit the proton fluxes with a spline w/ nnodes1 nodes:
			nnodes = nnodes1;
			Spline *sp_p1 = Spline::BuildFromHistogram(h_BR_p_all[iBR], Form("fsp_BR_p_%02d", iBR), nnodes, Spline::LogLog | Spline::PowerLaw);
			TF1 *fsp_p1 = sp_p1->GetTF1Pointer();
			h_BR_p_all[iBR]->Fit(fsp_p1, "NQ"); // do not draw and print the fit
			h_BR_p_all[iBR]->Fit(fsp_p1, "NIQ"); // do not draw the fit, print the fit results, minimize the integral of the function

			// Compute the fit error and residuals:
			TH1 *h_fitres_p1 = HistTools::GetResiduals(h_BR_p_all[iBR], fsp_p1, "_fitratio", false, true, true, 4, 1); // this will compute data/fit
			HistTools::CopyStyle(h_BR_p_all[iBR], h_fitres_p1);
			TH1 *h_fiterr_p1 = HistTools::GetFitError(h_BR_p_all[iBR], fsp_p1, "_fiterr", true, false, true, 11, 1.); // this will compute the relative error of the fit, centered in 1 (so 1.01 = 1% upward error; 0.98 = 2% downward error, etc)
			HistTools::SetFillStyle(h_fiterr_p1, kRed-7, 1001);

			// Plot the proton fluxes fit results:
			// Create a new canvas 
			auto c1 = new TCanvas("c1","AMS Monthly Proton & Helium Fluxes at Low R & Integrated Flux at High R", 1600, 800);
			c1->Divide(2,2);
			c1->cd(1);
			gPad->SetLogx();
			gPad->SetLogy();
			gPad->SetGrid();
			h_BR_p_all[iBR]->SetTitle(Form("AMS-02 Proton Flux BR-%02d", iBR));
			h_BR_p_all[iBR]->Draw("E1X0");
			fsp_p1->Draw("SAME");

			c1->cd(3);
   			h_fiterr_p1->SetTitle("");
   			h_fiterr_p1->SetXTitle("Rigidity [GV]");
   			h_fiterr_p1->SetYTitle("Data / Fit");
   			h_fiterr_p1->Draw("E3");
   			h_fiterr_p1->GetYaxis()->SetRangeUser(0.9, 1.1);
   			gPad->SetLogx();
   			gPad->SetGrid();
   			h_fitres_p1->Draw("E1X0 SAME");

			// Fit the Helium fluxes with a spline w/ nnodes1 nodes:
			nnodes = nnodes1;
			Spline *sp_he1 = Spline::BuildFromHistogram(h_BR_he_all[iBR], Form("fsp_BR_he_%02d", iBR), nnodes, Spline::LogLog | Spline::PowerLaw);
			TF1 *fsp_he1 = sp_he1->GetTF1Pointer();
			h_BR_he_all[iBR]->Fit(fsp_he1, "NQ"); // do not draw and print the fit
			h_BR_he_all[iBR]->Fit(fsp_he1, "NIQ"); // do not draw the fit, print the fit results, minimize the integral of the function

			// Compute the Helium fluxes fit error and residuals:
			TH1 *h_fitres_he1 = HistTools::GetResiduals(h_BR_he_all[iBR], fsp_he1, "_fitratio", false, true, true, 4, 1); // this will compute data/fit
			HistTools::CopyStyle(h_BR_he_all[iBR], h_fitres_he1);
			TH1 *h_fiterr_he1 = HistTools::GetFitError(h_BR_he_all[iBR], fsp_he1, "_fiterr", true, false, true, 11, 1.); // this will compute the relative error of the fit, centered in 1 (so 1.01 = 1% upward error; 0.98 = 2% downward error, etc)
			HistTools::SetFillStyle(h_fiterr_he1, kRed-7, 1001);

			// Plot the Helium fluxes fit results:
			c1->cd(2);
			gPad->SetLogx();
			gPad->SetLogy();
			gPad->SetGrid();
			h_BR_he_all[iBR]->SetTitle(Form("AMS-02 Helium Flux BR-%02d", iBR));
			h_BR_he_all[iBR]->Draw("E1X0");
			fsp_he1->Draw("SAME");

			c1->cd(4);
   			h_fiterr_he1->SetTitle("");
   			h_fiterr_he1->SetXTitle("Rigidity [GV]");
   			h_fiterr_he1->SetYTitle("Data / Fit");
   			h_fiterr_he1->Draw("E3");
   			h_fiterr_he1->GetYaxis()->SetRangeUser(0.9, 1.1);
   			gPad->SetLogx();
   			gPad->SetGrid();
   			h_fitres_he1->Draw("E1X0 SAME");

			// Save the canvas		
			c1->Print(Form("./data/amsfit/BR-%dvs%d-%02d.png", nnodes1, nnodes2, iBR));

			h_BR_p_all[iBR]->Write(Form("h_BR_p_%02d", iBR));
			fsp_p1->Write();
			h_fitres_p1->Write(Form("h_BR_p_fitres_%02d", iBR));
			h_fiterr_p1->Write(Form("h_BR_p_fiterr_%02d", iBR));

			h_BR_he_all[iBR]->Write(Form("h_BR_he_%02d", iBR));
			fsp_he1->Write();
			h_fitres_he1->Write(Form("h_BR_he_fitres_%02d", iBR));
			h_fiterr_he1->Write(Form("h_BR_he_fiterr_%02d", iBR));


			// Fit the proton fluxes with a spline w/ nnodes2 nodes:
			nnodes = nnodes2;
			Spline *sp_p2 = Spline::BuildFromHistogram(h_BR_p_all[iBR], Form("fsp_BR_p_%02d", iBR), nnodes, Spline::LogLog | Spline::PowerLaw);
			TF1 *fsp_p2 = sp_p2->GetTF1Pointer();
			h_BR_p_all[iBR]->Fit(fsp_p2, "NQ"); // do not draw and print the fit
			h_BR_p_all[iBR]->Fit(fsp_p2, "NIQ"); // do not draw the fit, print the fit results, minimize the integral of the function

			// Compute the fit error and residuals:
			TH1 *h_fitres_p2 = HistTools::GetResiduals(h_BR_p_all[iBR], fsp_p2, "_fitratio", false, true, true, 4, 1); // this will compute data/fit
			HistTools::CopyStyle(h_BR_p_all[iBR], h_fitres_p2);
			TH1 *h_fiterr_p2 = HistTools::GetFitError(h_BR_p_all[iBR], fsp_p2, "_fiterr", true, false, true, 11, 1.); // this will compute the relative error of the fit, centered in 1 (so 1.01 = 1% upward error; 0.98 = 2% downward error, etc)
			HistTools::SetFillStyle(h_fiterr_p2, kRed-7, 1001);

			// Fit the Helium fluxes with a spline w/ nnodes2 nodes:
			nnodes = nnodes2;
			Spline *sp_he2 = Spline::BuildFromHistogram(h_BR_he_all[iBR], Form("fsp_BR_he_%02d", iBR), nnodes, Spline::LogLog | Spline::PowerLaw);
			TF1 *fsp_he2 = sp_he2->GetTF1Pointer();
			h_BR_he_all[iBR]->Fit(fsp_he2, "NQ"); // do not draw and print the fit
			h_BR_he_all[iBR]->Fit(fsp_he2, "NIQ"); // do not draw the fit, print the fit results, minimize the integral of the function

			// Compute the fit error and residuals:
			TH1 *h_fitres_he2 = HistTools::GetResiduals(h_BR_he_all[iBR], fsp_he2, "_fitratio", false, true, true, 4, 1); // this will compute data/fit
			HistTools::CopyStyle(h_BR_p_all[iBR], h_fitres_he2);
			TH1 *h_fiterr_he2 = HistTools::GetFitError(h_BR_he_all[iBR], fsp_he2, "_fiterr", true, false, true, 11, 1.); // this will compute the relative error of the fit, centered in 1 (so 1.01 = 1% upward error; 0.98 = 2% downward error, etc)
			HistTools::SetFillStyle(h_fiterr_he2, kRed-7, 1001);
		
			// Chi-2 Test for Different Spline Node #

			// p, nnodes1 
			double chi2_p1 = fsp_p1->GetChisquare(); // chi-squared = sum of the square of the residuals, where residual = (data - fit)/data_error
			double dof_p1 = fsp_p1->GetNDF(); // number of degrees of freedom = number of data points - number of free parameters in the fit: note that this number should be the same for all BRs
			double norm_chi2_p1 = chi2_p1/dof_p1; // normalized chi-squared

			// he, nnodes1 
			double chi2_he1 = fsp_he1->GetChisquare(); // chi-squared = sum of the square of the residuals, where residual = (data - fit)/data_error
			double dof_he1 = fsp_he1->GetNDF(); // number of degrees of freedom = number of data points - number of free parameters in the fit: note that this number should be the same for all BRs
			double norm_chi2_he1 = chi2_he1/dof_he1; // normalized chi-squared

			// p, nnodes2 
			double chi2_p2 = fsp_p2->GetChisquare(); // chi-squared = sum of the square of the residuals, where residual = (data - fit)/data_error
			double dof_p2 = fsp_p2->GetNDF(); // number of degrees of freedom = number of data points - number of free parameters in the fit: note that this number should be the same for all BRs
			double norm_chi2_p2 = chi2_p2/dof_p2; // normalized chi-squared

			// he, nnodes2 
			double chi2_he2 = fsp_he2->GetChisquare(); // chi-squared = sum of the square of the residuals, where residual = (data - fit)/data_error
			double dof_he2 = fsp_he2->GetNDF(); // number of degrees of freedom = number of data points - number of free parameters in the fit: note that this number should be the same for all BRs
			double norm_chi2_he2 = chi2_he2/dof_he2; // normalized chi-squared

			double pvalue = TMath::Prob(chi2_p1 - chi2_p2, dof_p1 - dof_p2); // proton
			double sigma = TMath::Sqrt2()*TMath::ErfcInverse(pvalue);

			double pvalue2 = TMath::Prob(chi2_he1 - chi2_he2, dof_he1 - dof_he2); // helium
			double sigma2 = TMath::Sqrt2()*TMath::ErfcInverse(pvalue2);

			// cout << iBR <<" " << chi2_p1 <<" "<< chi2_p2 <<" "<< norm_chi2_p1 <<" "<< norm_chi2_p2 <<" "<< pvalue <<" "<< sigma << endl;
	
			// Store Sigma Test Data
			h_p_sig->SetPoint(iBR+1, iBR+1, sigma);
			h_he_sig->SetPoint(iBR+1, iBR+1, sigma2);

		}

		// Plot chi-2 sigma vs. iBR 

		for (iBR = 0; iBR < nBRs; ++iBR)
		{
			h_p_sig->GetPoint(iBR, iBR_D, sigma);
			printf("iBR=%d, sigma=%f \n", iBR, sigma); 

		}

		cout << endl;

		for (iBR = 0; iBR < nBRs; ++iBR)
		{
			h_he_sig->GetPoint(iBR, iBR_D, sigma2);
			printf("iBR=%d, sigma2=%f \n", iBR, sigma2); 

		}
		
		TCanvas *c2 = new TCanvas("c2","Fitting Sigma Comparison",1800,900);
		c2->Divide(2,1);

		c2->cd(1);
		h_p_sig->SetTitle(Form("Proton Fitting Sigma (# of Nodes %d vs. %d) vs. BR index;BR index;Sigma", nnodes1, nnodes2));
		h_p_sig->GetXaxis()->SetRangeUser(0., 79.);
		h_p_sig->SetMarkerStyle(20);
   		h_p_sig->SetMarkerColor(kRed);
   		h_p_sig->SetMarkerSize(1.1);
		h_p_sig->Draw("AP");

		c2->cd(2);
		h_he_sig->SetTitle(Form("Helium Fitting Sigma (# of Nodes %d vs. %d) vs. BR index;BR index;Sigma", nnodes1, nnodes2));
		h_he_sig->GetXaxis()->SetRangeUser(0., 79.);
		h_he_sig->SetMarkerStyle(20);
   		h_he_sig->SetMarkerColor(kBlue);
   		h_he_sig->SetMarkerSize(1.1);
		h_he_sig->Draw("AP");

		c2->Print(Form("./data/amsfit/sigmatest%dvs%d.png", nnodes1, nnodes2));

		//h_BR_p_all[0]->Print("range");
		//h_BR_he_all[0]->Print("range");

		// Write in Results
		file_fitresult.Write();
		file_fitresult.Close();

	} // close nnodes2 loop

   } // close nnodes1 loop

	// return 0;
}

//continue with He flux 
