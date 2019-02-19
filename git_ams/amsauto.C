// Fit AMS p, He Fluxes by Lingqiang He and Claudio Corti, 01/20/2019
// Compare statistical fitting goodness w/ diferent node #s by amsnode()
// node = 6 for the best fit stats 

// To check dataset information, use "Experiments::Info[Experiments::AMS02].Print();"
// amsauto("element (lower case)", "redo fit? (1 or 0)", node_min, node max);

#include "commonlib/include/Experiments.hh"
#include "commonlib/include/HistTools.hh"
#include "commonlib/include/Spline.hh"

void amsfit(int nnodes);
void amsnode(char const *item, int i, int j);
void amsnodeBR(char const *item, int i, int j); 

void amsauto(char const *item, int redo, int node_min, int node_max){

	gROOT->ProcessLine(".L amsauto.C");

	for (int i=node_min; i<=node_max; i++){

		
		if (redo == 1){ 
			amsfit(i);
		}else if (redo == 0){	
			continue;
		}

	} // close i loop 

	for (int i=node_min; i<=node_max; i++){

		for (int j=i+1; j<=node_max; j++){

			if ( strcmp(item, "li") == 0 || strcmp(item, "be") == 0 || strcmp(item, "b") == 0 || strcmp(item, "c") == 0 || strcmp(item, "o") == 0){
				amsnode(item, i, j);
			}else if ( strcmp(item, "p") == 0 || strcmp(item, "he") == 0){
				amsnodeBR(item, i, j);
			}
	
		} // close j loop
 
	} // close i loop

}

void amsfit(int nnodes){

	gROOT->SetBatch(); // disable canvas loop	

	gStyle->SetOptStat(0);
	Experiments::DataPath = "data";
	// Experiments::Info[Experiments::AMS02].Print();

	// Define variables
	double sigma = 0, sigma2 = 0;
	int iBR = 0; // BR index
	
	double iBR_D = (double)iBR;

	// Load the AMS monthly p and He fluxes:
	const int nBRs = Experiments::Info[Experiments::AMS02].Dataset[1].nMeasurements; // read the number of BRs

	TH1D **h_BR_p = Experiments::GetDatasetHistograms(Experiments::AMS02, 1);
	TH1D **h_BR_he = Experiments::GetDatasetHistograms(Experiments::AMS02, 4);

	// Load the AMS integrated fluxes: 
	TH1D *h_p = Experiments::GetMeasurementHistogram(Experiments::AMS02, 0, 0);  // proton
	TH1D *h_he = Experiments::GetMeasurementHistogram(Experiments::AMS02, 18, 0);  // helium 
	TH1D *h_li = Experiments::GetMeasurementHistogram(Experiments::AMS02, 23, 0);  // lithium
	TH1D *h_be = Experiments::GetMeasurementHistogram(Experiments::AMS02, 27, 0);  // beryllium
	TH1D *h_b = Experiments::GetMeasurementHistogram(Experiments::AMS02, 31, 0);  // boron
	TH1D *h_c = Experiments::GetMeasurementHistogram(Experiments::AMS02, 20, 0);  // carbon
	TH1D *h_o = Experiments::GetMeasurementHistogram(Experiments::AMS02, 22, 0);  // oxygen

	// Create a new set of histograms which uses the monthly fluxes at low rigidities and the integrated fluxes at high rigidities:
	TH1D **h_BR_p_all = new TH1D *[nBRs]; // this creates a new array of pointers with nBRs items
	TH1D **h_BR_he_all = new TH1D *[nBRs]; // this creates a new array of pointers with nBRs items	

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

	} // acquiring the data

	// Open TFile
	TFile file_fitresult(Form("data/amsfit/fit_result_node%d.root", nnodes), "RECREATE");

	// Fit Li, Be, B, C, O

			Spline *sp_li = Spline::BuildFromHistogram(h_li, "fsp_li", nnodes, Spline::LogLog | Spline::PowerLaw);
			TF1 *fsp_li = sp_li->GetTF1Pointer(); 
			h_li->Fit(fsp_li, "NQ"); // do not draw and print the fit
			h_li->Fit(fsp_li, "NIQ"); // do not draw the fit, print the fit results, minimize the integral of the function 

			TH1 *h_fitres_li = HistTools::GetResiduals(h_li, fsp_li, "_fitratio", false, true, true, 4, 1); // this will compute data/fit
			HistTools::CopyStyle(h_li, h_fitres_li);
			TH1 *h_fiterr_li = HistTools::GetFitError(h_li, fsp_li, "_fiterr", true, false, true, 11, 1.); // this will compute the relative error of the fit, centered in 1 (so 1.01 = 1% upward error; 0.98 = 2% downward error, etc)
			HistTools::SetFillStyle(h_fiterr_li, kRed-7, 1001);

			Spline *sp_be = Spline::BuildFromHistogram(h_be, "fsp_be", nnodes, Spline::LogLog | Spline::PowerLaw);
			TF1 *fsp_be = sp_be->GetTF1Pointer(); 
			h_be->Fit(fsp_be, "NQ"); // do not draw and print the fit
			h_be->Fit(fsp_be, "NIQ"); // do not draw the fit, print the fit results, minimize the integral of the function 

			TH1 *h_fitres_be = HistTools::GetResiduals(h_be, fsp_be, "_fitratio", false, true, true, 4, 1); // this will compute data/fit
			HistTools::CopyStyle(h_be, h_fitres_be);
			TH1 *h_fiterr_be = HistTools::GetFitError(h_be, fsp_be, "_fiterr", true, false, true, 11, 1.); // this will compute the relative error of the fit, centered in 1 (so 1.01 = 1% upward error; 0.98 = 2% downward error, etc)
			HistTools::SetFillStyle(h_fiterr_be, kRed-7, 1001);

			Spline *sp_b = Spline::BuildFromHistogram(h_b, "fsp_b", nnodes, Spline::LogLog | Spline::PowerLaw);
			TF1 *fsp_b = sp_b->GetTF1Pointer(); 
			h_b->Fit(fsp_b, "NQ"); // do not draw and print the fit
			h_b->Fit(fsp_b, "NIQ"); // do not draw the fit, print the fit results, minimize the integral of the function 

			TH1 *h_fitres_b = HistTools::GetResiduals(h_b, fsp_b, "_fitratio", false, true, true, 4, 1); // this will compute data/fit
			HistTools::CopyStyle(h_b, h_fitres_b);
			TH1 *h_fiterr_b = HistTools::GetFitError(h_b, fsp_b, "_fiterr", true, false, true, 11, 1.); // this will compute the relative error of the fit, centered in 1 (so 1.01 = 1% upward error; 0.98 = 2% downward error, etc)
			HistTools::SetFillStyle(h_fiterr_b, kRed-7, 1001);

			Spline *sp_c = Spline::BuildFromHistogram(h_c, "fsp_c", nnodes, Spline::LogLog | Spline::PowerLaw);
			TF1 *fsp_c = sp_c->GetTF1Pointer(); 
			h_c->Fit(fsp_c, "NQ"); // do not draw and print the fit
			h_c->Fit(fsp_c, "NIQ"); // do not draw the fit, print the fit results, minimize the integral of the function 

			TH1 *h_fitres_c = HistTools::GetResiduals(h_c, fsp_c, "_fitratio", false, true, true, 4, 1); // this will compute data/fit
			HistTools::CopyStyle(h_c, h_fitres_c);
			TH1 *h_fiterr_c = HistTools::GetFitError(h_c, fsp_c, "_fiterr", true, false, true, 11, 1.); // this will compute the relative error of the fit, centered in 1 (so 1.01 = 1% upward error; 0.98 = 2% downward error, etc)
			HistTools::SetFillStyle(h_fiterr_c, kRed-7, 1001);

			Spline *sp_o = Spline::BuildFromHistogram(h_o, "fsp_o", nnodes, Spline::LogLog | Spline::PowerLaw);
			TF1 *fsp_o = sp_o->GetTF1Pointer(); 
			h_o->Fit(fsp_o, "NQ"); // do not draw and print the fit
			h_o->Fit(fsp_o, "NIQ"); // do not draw the fit, print the fit results, minimize the integral of the function 

			TH1 *h_fitres_o = HistTools::GetResiduals(h_o, fsp_o, "_fitratio", false, true, true, 4, 1); // this will compute data/fit
			HistTools::CopyStyle(h_o, h_fitres_o);
			TH1 *h_fiterr_o = HistTools::GetFitError(h_o, fsp_o, "_fiterr", true, false, true, 11, 1.); // this will compute the relative error of the fit, centered in 1 (so 1.01 = 1% upward error; 0.98 = 2% downward error, etc)
			HistTools::SetFillStyle(h_fiterr_o, kRed-7, 1001);

			// Draw fit results on TCanvas
			TCanvas *c0 = new TCanvas("c0","AMS-002 Li, Be, B, C, O Integrated Flux", 2400, 1600);
			c0->Divide(3, 4);
			
			c0->cd(1);
			gPad->SetLogx();
			gPad->SetLogy();
			gPad->SetGrid();
			h_li->SetTitle("AMS-02 Lithium Integrated Flux");
			h_li->SetMarkerStyle(20);
   			h_li->SetMarkerColor(kBlue);
   			h_li->SetMarkerSize(1.1);
			h_li->Draw("E1X0");
			fsp_li->Draw("SAME");

			c0->cd(4);
   			h_li->SetTitle("");
   			h_li->SetXTitle("Rigidity [GV]");
   			h_li->SetYTitle("Data / Fit");
   			h_fiterr_li->Draw("E3");
   			h_fiterr_li->GetYaxis()->SetRangeUser(0.9, 1.1);
   			gPad->SetLogx();
   			gPad->SetGrid();
   			h_fitres_li->Draw("E1X0 SAME");

			c0->cd(2);
			gPad->SetLogx();
			gPad->SetLogy();
			gPad->SetGrid();
			h_be->SetTitle("AMS-02 Beryllium Integrated Flux");
			h_be->SetMarkerStyle(20);
   			h_be->SetMarkerColor(kBlue);
   			h_be->SetMarkerSize(1.1);
			h_be->Draw("E1X0");
			fsp_be->Draw("SAME");

			c0->cd(5);
   			h_be->SetTitle("");
   			h_be->SetXTitle("Rigidity [GV]");
   			h_be->SetYTitle("Data / Fit");
   			h_fiterr_be->Draw("E3");
   			h_fiterr_be->GetYaxis()->SetRangeUser(0.9, 1.1);
   			gPad->SetLogx();
   			gPad->SetGrid();
   			h_fitres_be->Draw("E1X0 SAME");

			c0->cd(3);
			gPad->SetLogx();
			gPad->SetLogy();
			gPad->SetGrid();
			h_b->SetTitle("AMS-02 Boron Integrated Flux");
			h_b->SetMarkerStyle(20);
   			h_b->SetMarkerColor(kBlue);
   			h_b->SetMarkerSize(1.1);
			h_b->Draw("E1X0");
			fsp_b->Draw("SAME");

			c0->cd(6);
   			h_b->SetTitle("");
   			h_b->SetXTitle("Rigidity [GV]");
   			h_b->SetYTitle("Data / Fit");
   			h_fiterr_b->Draw("E3");
   			h_fiterr_b->GetYaxis()->SetRangeUser(0.9, 1.1);
   			gPad->SetLogx();
   			gPad->SetGrid();
   			h_fitres_b->Draw("E1X0 SAME");
	
			c0->cd(7);
			gPad->SetLogx();
			gPad->SetLogy();
			gPad->SetGrid();
			h_c->SetTitle("AMS-02 Carbon Integrated Flux");
			h_c->SetMarkerStyle(20);
   			h_c->SetMarkerColor(kBlue);
   			h_c->SetMarkerSize(1.1);
			h_c->Draw("E1X0");
			fsp_c->Draw("SAME");

			c0->cd(10);
   			h_c->SetTitle("");
   			h_c->SetXTitle("Rigidity [GV]");
   			h_c->SetYTitle("Data / Fit");
   			h_fiterr_c->Draw("E3");
   			h_fiterr_c->GetYaxis()->SetRangeUser(0.9, 1.1);
   			gPad->SetLogx();
   			gPad->SetGrid();
   			h_fitres_c->Draw("E1X0 SAME");

			c0->cd(8);
			gPad->SetLogx();
			gPad->SetLogy();
			gPad->SetGrid();
			h_o->SetTitle("AMS-02 Oxygen Integrated Flux");
			h_o->SetMarkerStyle(20);
   			h_o->SetMarkerColor(kBlue);
   			h_o->SetMarkerSize(1.1);
			h_o->Draw("E1X0");
			fsp_o->Draw("SAME");

			c0->cd(11);
   			h_o->SetTitle("");
   			h_o->SetXTitle("Rigidity [GV]");
   			h_o->SetYTitle("Data / Fit");
   			h_fiterr_o->Draw("E3");
   			h_fiterr_o->GetYaxis()->SetRangeUser(0.9, 1.1);
   			gPad->SetLogx();
   			gPad->SetGrid();
   			h_fitres_o->Draw("E1X0 SAME");

			c0->Print(Form("./data/amsfit/nodes%d/node%d_rest.png", nnodes, nnodes));

			// Store fit results
			h_li->Write("h_li");
			fsp_li->Write();
			h_fitres_li->Write("h_li_fitres");
			h_fiterr_li->Write("h_li_fiterr");

			h_be->Write("h_be");
			fsp_be->Write();
			h_fitres_be->Write("h_be_fitres");
			h_fiterr_be->Write("h_be_fiterr");

			h_b->Write("h_b");
			fsp_b->Write();
			h_fitres_b->Write("h_b_fitres");
			h_fiterr_b->Write("h_b_fiterr");

			h_c->Write("h_c");
			fsp_c->Write();
			h_fitres_c->Write("h_c_fitres");
			h_fiterr_c->Write("h_c_fiterr");

			h_o->Write("h_o");
			fsp_o->Write();
			h_fitres_o->Write("h_o_fitres");
			h_fiterr_o->Write("h_o_fiterr");

	// Fit p, He

	for (iBR = 0; iBR < nBRs; iBR++)
	{
			gSystem->mkdir(Form("data/amsfit/nodes%d", nnodes), true);

			// Fit the proton fluxes with a spline w/ nnodes nodes:
			Spline *sp_p = Spline::BuildFromHistogram(h_BR_p_all[iBR], Form("fsp_BR_p_%02d", iBR), nnodes, Spline::LogLog | Spline::PowerLaw);
			TF1 *fsp_p = sp_p->GetTF1Pointer(); 
			h_BR_p_all[iBR]->Fit(fsp_p, "NQ"); // do not draw and print the fit
			h_BR_p_all[iBR]->Fit(fsp_p, "NIQ"); // do not draw the fit, print the fit results, minimize the integral of the function
			// printf("BR=%d chi2=%17.10e\n", iBR, fsp_p->GetChisquare());

			// Compute the fit error and residuals:
			TH1 *h_fitres_p = HistTools::GetResiduals(h_BR_p_all[iBR], fsp_p, "_fitratio", false, true, true, 4, 1); // this will compute data/fit
			HistTools::CopyStyle(h_BR_p_all[iBR], h_fitres_p);
			TH1 *h_fiterr_p = HistTools::GetFitError(h_BR_p_all[iBR], fsp_p, "_fiterr", true, false, true, 11, 1.); // this will compute the relative error of the fit, centered in 1 (so 1.01 = 1% upward error; 0.98 = 2% downward error, etc)
			HistTools::SetFillStyle(h_fiterr_p, kRed-7, 1001);

			// Plot the proton fluxes fit results:
			// Create a new canvas 
			TCanvas *c1 = new TCanvas("c1","AMS Monthly Proton & Helium Fluxes at Low R & Integrated Flux at High R", 1600, 800);
			c1->Divide(2,2);
			c1->cd(1);
			gPad->SetLogx();
			gPad->SetLogy();
			gPad->SetGrid();
			h_BR_p_all[iBR]->SetTitle(Form("AMS-02 Proton Flux BR-%02d", iBR));
			h_BR_p_all[iBR]->Draw("E1X0");
			fsp_p->Draw("SAME");

			c1->cd(3);
   			h_fiterr_p->SetTitle("");
   			h_fiterr_p->SetXTitle("Rigidity [GV]");
   			h_fiterr_p->SetYTitle("Data / Fit");
   			h_fiterr_p->Draw("E3");
   			h_fiterr_p->GetYaxis()->SetRangeUser(0.9, 1.1);
   			gPad->SetLogx();
   			gPad->SetGrid();
   			h_fitres_p->Draw("E1X0 SAME");

			// Fit the Helium fluxes with a spline w/ nnodes1 nodes:
			// nnodes = nnodes1;
			Spline *sp_he = Spline::BuildFromHistogram(h_BR_he_all[iBR], Form("fsp_BR_he_%02d", iBR), nnodes, Spline::LogLog | Spline::PowerLaw);
			TF1 *fsp_he = sp_he->GetTF1Pointer();
			h_BR_he_all[iBR]->Fit(fsp_he, "NQ"); // do not draw and print the fit
			h_BR_he_all[iBR]->Fit(fsp_he, "NIQ"); // do not draw the fit, print the fit results, minimize the integral of the function

			// Compute the Helium fluxes fit error and residuals:
			TH1 *h_fitres_he = HistTools::GetResiduals(h_BR_he_all[iBR], fsp_he, "_fitratio", false, true, true, 4, 1); // this will compute data/fit
			HistTools::CopyStyle(h_BR_he_all[iBR], h_fitres_he);
			TH1 *h_fiterr_he = HistTools::GetFitError(h_BR_he_all[iBR], fsp_he, "_fiterr", true, false, true, 11, 1.); // this will compute the relative error of the fit, centered in 1 (so 1.01 = 1% upward error; 0.98 = 2% downward error, etc)
			HistTools::SetFillStyle(h_fiterr_he, kRed-7, 1001);

			// Plot the Helium fluxes fit results:
			c1->cd(2);
			gPad->SetLogx();
			gPad->SetLogy();
			gPad->SetGrid();
			h_BR_he_all[iBR]->SetTitle(Form("AMS-02 Helium Flux BR-%02d", iBR));
			h_BR_he_all[iBR]->Draw("E1X0");
			fsp_he->Draw("SAME");

			c1->cd(4);
   			h_fiterr_he->SetTitle("");
   			h_fiterr_he->SetXTitle("Rigidity [GV]");
   			h_fiterr_he->SetYTitle("Data / Fit");
   			h_fiterr_he->Draw("E3");
   			h_fiterr_he->GetYaxis()->SetRangeUser(0.9, 1.1);
   			gPad->SetLogx();
   			gPad->SetGrid();
   			h_fitres_he->Draw("E1X0 SAME");

			// Save the canvas		
			c1->Print(Form("./data/amsfit/nodes%d/node%d_BR%02d.png", nnodes, nnodes, iBR));

			h_BR_p_all[iBR]->Write(Form("h_BR_p_%02d", iBR));
			fsp_p->Write();
			h_fitres_p->Write(Form("h_BR_p_fitres_%02d", iBR));
			h_fiterr_p->Write(Form("h_BR_p_fiterr_%02d", iBR));

			h_BR_he_all[iBR]->Write(Form("h_BR_he_%02d", iBR));
			fsp_he->Write();
			h_fitres_he->Write(Form("h_BR_he_fitres_%02d", iBR));
			h_fiterr_he->Write(Form("h_BR_he_fiterr_%02d", iBR));

			// cout << iBR <<" " << chi2_p <<" "<< chi2_p <<" "<< norm_chi2_p <<" "<< norm_chi2_p <<" "<< pvalue <<" "<< sigma << endl;
			
		}

		// Plot chi-2 sigma vs. iBR 

		//h_BR_p_all[0]->Print("range");
		//h_BR_he_all[0]->Print("range");

		// Write in Results
		file_fitresult.Write();
		file_fitresult.Close();


	// return 0;
}

void amsnode(char const *item, int i, int j){

			gROOT->SetBatch(); 		

			gStyle->SetOptStat(0);
			Experiments::DataPath = "data";

			const int nBRs = Experiments::Info[Experiments::AMS02].Dataset[1].nMeasurements; 

			TFile file1(Form("data/amsfit/fit_result_node%d.root", i)); 
			TFile file2(Form("data/amsfit/fit_result_node%d.root", j)); 

			TF1 *fsp_1 = (TF1*)file1.Get(Form("fsp_%s", item)); 

			TF1 *fsp_2 = (TF1*)file2.Get(Form("fsp_%s", item)); 

			// w/ nnodes = i 
			double chi2_1 = fsp_1->GetChisquare(); // chi-squared = sum of the square of the residuals, where residual = (data - fit)/data_error
			double dof_1 = fsp_1->GetNDF(); // number of degrees of freedom = number of data points - number of free parameters in the fit
			double norm_chi2_1 = chi2_1/dof_1; // normalized chi-squared

			// w/ nnodes = j 
			double chi2_2 = fsp_2->GetChisquare(); // chi-squared = sum of the square of the residuals, where residual = (data - fit)/data_error
			double dof_2 = fsp_2->GetNDF(); // number of degrees of freedom = number of data points - number of free parameters in the fit
			double norm_chi2_2 = chi2_2/dof_2; // normalized chi-squared

			double pvalue = TMath::Prob(chi2_1 - chi2_2, dof_1 - dof_2);
			double sigma = TMath::Sqrt2()*TMath::ErfcInverse(pvalue);	
			
			printf("%s fit, node %dvs%d, sigma = %f \n", item, i, j, sigma);

}

void amsnodeBR(char const *item, int i, int j){ 

			gROOT->SetBatch(); 		

			gStyle->SetOptStat(0);
			Experiments::DataPath = "data";

			const int nBRs = Experiments::Info[Experiments::AMS02].Dataset[1].nMeasurements; 

			TGraph *g_sig = new TGraph();

			TFile file1(Form("data/amsfit/fit_result_node%d.root", i)); 
			TFile file2(Form("data/amsfit/fit_result_node%d.root", j)); 

			for (int iBR = 0; iBR < nBRs; ++iBR) 
			{

				TF1 *fsp_1 = (TF1*)file1.Get(Form("fsp_BR_%s_%02d", item, iBR)); 

				TF1 *fsp_2 = (TF1*)file2.Get(Form("fsp_BR_%s_%02d", item, iBR)); 

				// w/ nnodes = i 
				double chi2_1 = fsp_1->GetChisquare(); // chi-squared = sum of the square of the residuals, where residual = (data - fit)/data_error
				double dof_1 = fsp_1->GetNDF(); // number of degrees of freedom = number of data points - number of free parameters in the fit
				double norm_chi2_1 = chi2_1/dof_1; // normalized chi-squared

				// w/ nnodes = j 
				double chi2_2 = fsp_2->GetChisquare(); // chi-squared = sum of the square of the residuals, where residual = (data - fit)/data_error
				double dof_2 = fsp_2->GetNDF(); // number of degrees of freedom = number of data points - number of free parameters in the fit
				double norm_chi2_2 = chi2_2/dof_2; // normalized chi-squared

				double pvalue = TMath::Prob(chi2_1 - chi2_2, dof_1 - dof_2);
				double sigma = TMath::Sqrt2()*TMath::ErfcInverse(pvalue);	
				
				g_sig->SetPoint(iBR, iBR, sigma);

				// printf("iBR=%d, sigma=%f \n", iBR, sigma); 

				
			}

	TCanvas *c2 = new TCanvas("c2","Fitting Sigma Comparison",1800,900);
	//c2->Divide(2,1);

	c2->cd(1);
	g_sig->SetTitle(Form("%s Fitting Sigma (# of Nodes %d vs. %d) vs. BR index;BR index;Sigma", item, i, j));
	g_sig->GetXaxis()->SetRangeUser(0., 79.);
	g_sig->SetMarkerStyle(20);
   	g_sig->SetMarkerColor(kRed);
   	g_sig->SetMarkerSize(1.1);
	g_sig->Draw("AP");

	c2->Print(Form("./data/amsfit/%s_sigma_test_%dvs%d.png", item, i, j));

}


