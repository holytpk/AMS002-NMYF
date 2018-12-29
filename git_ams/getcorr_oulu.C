//Get correlation then fit a linear regression for selected neutron monitor stations
//Also graph fitting slopes of correlation of OULU vs. all other stations.   
//Consider the fact of different rigidity cutoff

#include <TGraph.h>
#include <TFile.h>
#include <cstdio>
#include <TString.h>
#include "commonlib/include/Experiments.hh"
#include "commonlib/include/HistTools.hh"
#include "commonlib/include/Spline.hh"
#include "commonlib/include/DateTimeTools.hh"

using namespace std;

TGraph *GetCorrelation2(TGraph *g1, TGraph *g2);

void getcorr_oulu(){

   Experiments::DataPath = "data";
   const int nNMs = 46; //# of NM stations
   const char *NM_name[nNMs] = { "PSNM", "TIBT", "DJON", "TSMB", "ATHN", "MXCO", "ARNM", "NANM", "PTFM", "CALM", "AATB", "ROME", "BKSN", "HRMS", "JUNG", "JUNG1", "LMKS", "IRK2", "IRK3", "IRKT", "DRBS", "MCRL", "MOSC", "NEWK", "KIEL", "KIEL2", "MGDN", "KERG", "OULU", "SANB", "SNAE", "APTY", "NRLK", "FSMT", "INVK", "JBGO", "MCMU", "NAIN", "PWNK", "THUL", "NEU3", "SOPB", "SOPO", "DOMB", "DOMC", "TERA" };

   const int nNMs_useful = 15;
   const char *NM_useful[nNMs_useful+1] = { "OULU", "PSNM", "MXCO", "HRMS", "JUNG", "JUNG1", "NEWK", "KIEL2", "APTY", "FSMT", "NAIN", "PWNK", "THUL", "SOPB", "SOPO" };

   const double Rcut[nNMs_useful+1] = { 0.81, 16.80, 8.28, 4.58, 4.49, 4.49, 2.40, 2.36, 0.65, 0.30, 0.30, 0.30, 0.30, 0.10, 0.10 }; // Rigidity Cutoff in GV
   const double alt[nNMs_useful+1] = { 15.0, 2565.0, 2274.0, 26.0, 3570.0, 3475.0, 50.0, 54.0, 181.0, 180.0, 46.0, 53.0, 26.0, 2820.0, 2820.0 }; // NM altitude

   auto c1 = new TCanvas("c1","NM BR-averaged Fluxes & Correlation",1600,900);
   c1->Divide(2,1);

   TH2D *mfit = new TH2D("mfit","Correlation Fit (slope);NM Stations;NM Stations", nNMs_useful, 0, nNMs_useful, nNMs_useful, 0, nNMs_useful );
   TH2D *qfit = new TH2D("qfit","Correlation Fit (constant);NM Stations;NM Stations", nNMs_useful, 0, nNMs_useful, nNMs_useful, 0, nNMs_useful );

   TGraphErrors *oulufit = new TGraphErrors("oulufit", "Oulu-Referenced Correlation Fit Slope vs. Rigidity Cutoff");
   TGraphAsymmErrors *oulucorr = new TGraphAsymmErrors("oulucorr", "Oulu-Reference Correlation Factor vs. Rigidity Cutoff");

   TGraphErrors *oulufit2 = new TGraphErrors("oulufit2", "Oulu-Referenced Correlation Fit Slope vs. Difference in Rigidity Cutoff");
   TGraphAsymmErrors *oulucorr2 = new TGraphAsymmErrors("oulucorr2", "Oulu-Reference Correlation Factor vs. Difference in Rigidity Cutoff");

   for(int i=0; i<=nNMs_useful; i++){
   	
	//printf("NM_useful: %s, rigidity cutoff: %f", NM_useful[i], Rcut[i]); 

   	for(int j=i+1; j<=nNMs_useful-1; j++){
   		   
		c1->cd(1);
		gPad->SetMargin(0.12, 0.08, 0.08, 0.08);
		//read in first file
		TFile fin1(Form("./data/nm/NM-%s.root", NM_useful[i]));
		TGraph *g1=(TGraph*) fin1.Get("g");
		TGraph *g1_ave=(TGraph*) fin1.Get("g_ave");
		if(g1_ave==NULL){
		      	printf("TGraph g1_ave not found \n");
		}else if(g1==NULL){
			printf("TGraph g1 not found \n");
		}else{
			g1_ave->SetTitle(Form("%s Flux & %s Flux;Time;Count Rate [Hz]", NM_useful[i], NM_useful[j]));
			g1_ave->GetXaxis()->SetTimeDisplay(1);
 			g1_ave->GetXaxis()->SetTimeFormat("%m-%y");
 			g1_ave->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00");
			g1_ave->GetYaxis()->SetRangeUser(0, 800);
			g1_ave->SetMarkerColor(kBlue);
			g1_ave->Draw("AP");
		}

		//read in second file
		TFile fin2(Form("./data/nm/NM-%s.root", NM_useful[j]));
		TGraph *g2=(TGraph*) fin2.Get("g");
		TGraph *g2_ave=(TGraph*) fin2.Get("g_ave");
		if(g2_ave==NULL){
			printf("TGraph g2_ave not found \n");
		}else if(g2==NULL){
			printf("TGraph g2 not found \n");
		}else{
			g2_ave->SetTitle(Form("%s Flux & %s Flux;Time;Count Rate [Hz]", NM_useful[i], NM_useful[j]));
			g2_ave->GetXaxis()->SetTimeDisplay(1);
	 		g2_ave->GetXaxis()->SetTimeFormat("%m-%y");
	 	 	g2_ave->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00");
			//g2_ave->GetYaxis()->SetRangeUser(0, 400);

			auto leg = new TLegend(0.1,0.7,0.5,0.9,"");
			leg->SetFillColor(0);
			leg->AddEntry(g1, Form("%s Daily Flux", NM_useful[i]),"l");
			leg->AddEntry(g1_ave,Form("%s BR-averaged Flux", NM_useful[i]),"p");
			leg->AddEntry(g2, Form("%s Daily Flux", NM_useful[j]),"l");
			leg->AddEntry(g2_ave,Form("%s BR-averaged Flux", NM_useful[j]),"p");
			leg->DrawClone("Same");

			g2_ave->SetMarkerColor(kRed);
			g2_ave->Draw("PSame");
		}

		//g1_ave->Print();
		//g2_ave->Print();

		g1->Draw("LSame");
		g2->Draw("LSame");

		c1->cd(2);
		gPad->SetMargin(0.12, 0.08, 0.08, 0.08);

		TGraph *gcor = GetCorrelation2(g1_ave, g2_ave);
		gcor->SetMarkerColor(kBlue);
		gcor->SetMarkerStyle(kFullCircle);
		gcor->SetTitle(Form("Correlation w/ Linear Regression;%s Count Rate [Hz];%s Count Rate [Hz]", NM_useful[i], NM_useful[j]));
		gcor->Draw("AP");

		double corr = gcor->GetCorrelationFactor();

		double corr_low=0, corr_up=0; 
		// these two variables will contain the lower and upper value of the correlation, accounting for statistical uncertainty at a given confidence level defined below
	
		double cl = TMath::Erf(1/sqrt(2)); 
		// confidence level corresponding to 1 sigma, roughly 0.68, which means ~68% probability

		double corr_pvalue = HistTools::CorrelationZTest(corr, gcor->GetN(), corr_low, corr_up, cl); 
		// the p-value is a number that tells you how much is probable that the correlation factor you computed is different from zero simply by chance, and not by some physical reason

		double corr_sign = TMath::ErfcInverse(corr_pvalue); 
		// get the statistical significance of the correlation being different from zero: the significance is given in unit of sigma. Any number above 3 is very significant, i.e. we can claim that most probably the computed correlation factor is different from zero not simply by chance.

		printf("%s vs. %s: correlation factor=%f, corr_low=%f, corr_up=%f \n", NM_useful[i], NM_useful[j], corr, corr_low, corr_up);
		//printf("%s vs. %s: correlation factor=%f, confidence level=%f, corr_low=%f, corr_up=%f, pvalue=%.8f, sign=%f \n", NM_useful[i], NM_useful[j], corr, cl, corr_low, corr_up, corr_pvalue, corr_sign);   		
		
		TF1 *f = new TF1("f", "[0]*x + [1]"); // y = m*x + q
		gcor->Fit(f);
		f->Draw("SAME");

		cout << NM_useful[i] <<" vs. "<< NM_useful[j] << endl;
		cout << i << "   vs.   " << j << endl;  

		if(i<1){		
			oulufit->SetPoint(j-1, Rcut[j], f->GetParameter(0) );
			oulufit->SetPointError(j-1, 0, f->GetParError(0) );

			oulucorr->SetPoint(j-1, Rcut[j], corr);
			oulucorr->SetPointError(j-1, 0, 0, corr-corr_low, corr_up-corr);

			oulufit2->SetPoint(j-1, Rcut[j]-Rcut[0], f->GetParameter(0) );
			oulufit2->SetPointError(j-1, 0, f->GetParError(0) );

			oulucorr2->SetPoint(j-1, Rcut[j]-Rcut[0], corr);
			oulucorr2->SetPointError(j-1, 0, 0, corr-corr_low, corr_up-corr);
		}

		mfit->SetBinContent(i+1, j+1, f->GetParameter(0) ); // get m		
		qfit->SetBinContent(i+1, j+1, f->GetParameter(1) ); // get q

		//c1->SaveAs(Form("./data/corr/Fit_%s_%s.png", NM_useful[i], NM_useful[j])); 

   	}

   } // end of selected NM combinations 

   for(int i=0; i<=nNMs_useful; i++){

	if (i<=nNMs_useful-1){
		printf("NM_useful: %s, rigidity cutoff: %f \n", NM_useful[i], Rcut[i]);
	}
	
	if (i>0){
		//oulufit->GetXaxis()->SetBinLabel(oulufit->GetXaxis()->FindBin(i), NM_useful[i]);
		TLatex *t1 = new TLatex(oulufit->GetX()[i-1], oulufit->GetY()[i-1], NM_useful[i]); 
		t1->SetTextSize(0.028);
		oulufit->GetListOfFunctions()->Add(t1); 
		oulufit->Draw("samelp");	

		TLatex *t2 = new TLatex(oulucorr->GetX()[i-1], oulucorr->GetY()[i-1], NM_useful[i]); 
		t2->SetTextSize(0.028);
		oulucorr->GetListOfFunctions()->Add(t2); 
		oulucorr->Draw("samelp");

		TLatex *t3 = new TLatex(oulufit2->GetX()[i-1], oulufit2->GetY()[i-1], NM_useful[i]); 
		t3->SetTextSize(0.028);
		oulufit2->GetListOfFunctions()->Add(t3); 
		oulufit2->Draw("samelp");	

		TLatex *t4 = new TLatex(oulucorr2->GetX()[i-1], oulucorr2->GetY()[i-1], NM_useful[i]); 
		t4->SetTextSize(0.028);
		oulucorr2->GetListOfFunctions()->Add(t4); 
		oulucorr2->Draw("samelp");
	}
   	
   	for(int j=i; j<=nNMs_useful-1; j++){

		//cout << i <<" "<< j << endl;

		mfit->GetXaxis()->SetBinLabel(j+1, NM_useful[j]);
		mfit->GetYaxis()->SetBinLabel(i+1, NM_useful[i]);
		qfit->GetXaxis()->SetBinLabel(j+1, NM_useful[j]);
		qfit->GetYaxis()->SetBinLabel(i+1, NM_useful[i]);

	}
   } // end of set label loop

   auto c2 = new TCanvas("c2","Correlation Fit Table",1600,900); 
   //c2->Divide(2,2);
   c2->cd(1);  
   gPad->SetMargin(0.12, 0.1, 0.08, 0.08);
   mfit->Draw("Colz");
   
   auto c3 = new TCanvas("c3","Oulu-Referenced Correlation Fit/Correlation Factor vs. Rigidity Cutoff", 1600, 900); 
   c3->Divide(2,2);
   c3->cd(1);
   gPad->SetLogx();
   gPad->SetMargin(0.11, 0.08, 0.08, 0.08);
   //qfit->Draw("Colz");
   oulufit->SetTitle("Oulu-Referenced Correlation Fit Slope vs. Rigidity Cutoff; NM Station Rigidity Cutoff (GV); Fit Slope Value");
   oulufit->SetMarkerStyle(20);
   oulufit->SetMarkerColor(kRed);
   oulufit->SetMarkerSize(0.7);
   //oulufit->GetXaxis()->SetTitleSize(0.4);
   oulufit->Draw("AP");
   c3->cd(2);  
   gPad->SetLogx();
   gPad->SetMargin(0.12, 0.1, 0.08, 0.08);
   oulucorr->SetTitle("Oulu-Referenced Correlation Factor vs. Rigidity Cutoff; NM Station Rigidity Cutoff (GV); Correlation Factor");
   oulucorr->SetMarkerStyle(20);
   oulucorr->SetMarkerColor(kRed);
   oulucorr->SetMarkerSize(0.7);
   //oulucorr->GetXaxis()->SetTitleSize(0.4);
   //oulucorr->GetYaxis()->SetRangeUser(0, 1.5);
   oulucorr->Draw("AP");
   c3->cd(3);
   gPad->SetLogx();
   gPad->SetMargin(0.11, 0.08, 0.08, 0.08);
   //qfit->Draw("Colz");
   oulufit2->SetTitle("Oulu-Referenced Correlation Fit Slope vs. Rigidity Cutoff Difference Respect to Oulu; NM Station Rigidity Cutoff Difference Respect to Oulu (GV); Fit Slope Value");
   oulufit2->SetMarkerStyle(20);
   oulufit2->SetMarkerColor(kRed);
   oulufit2->SetMarkerSize(0.7);
   //oulufit->GetXaxis()->SetTitleSize(0.4);
   oulufit2->Draw("AP");
   c3->cd(4);
   gPad->SetLogx();  
   gPad->SetMargin(0.12, 0.1, 0.08, 0.08);
   oulucorr2->SetTitle("Oulu-Referenced Correlation Factor vs. Rigidity Cutoff Difference Respect to Oulu; NM Station Rigidity Cutoff Difference Respect to Oulu (GV); Correlation Factor");
   oulucorr2->SetMarkerStyle(20);
   oulucorr2->SetMarkerColor(kRed);
   oulucorr2->SetMarkerSize(0.7);
   //oulucorr->GetXaxis()->SetTitleSize(0.4);
   //oulucorr2->GetYaxis()->SetRangeUser(0, 1.5);
   oulucorr2->Draw("AP");

   oulufit->Print();
   oulucorr->Print();
   oulufit2->Print();
   oulucorr2->Print();

   c3->SaveAs("./data/oulu-referenced.png"); 

}

//GetCorrelation ver. 2 for unmatched data points by Lingqiang He 11/29/2018
TGraph *GetCorrelation2(TGraph *g1, TGraph *g2){
   TGraph *g;

   UInt_t npoints = (g1->GetN() > g2->GetN()) ? g1->GetN() : g2->GetN() ;

     UInt_t ipoint = 0, ipoint1 = 0, ipoint2 = 0;   

   if (g1->GetEYlow() != NULL || g2->GetEYlow() != NULL)
   {
      g = new TGraphAsymmErrors(npoints);
   }
   else if (g1->GetEY() != NULL || g2->GetEY() != NULL) 
   {
      g = new TGraphErrors(npoints);
   }
   else
   {
      g = new TGraph(npoints);
   }

   g->SetTitle(Form("Correlation;%s;%s", g1->GetYaxis()->GetTitle(), g2->GetYaxis()->GetTitle()));

   while (true){

      Double_t y1   = g1->GetY()[ipoint];
      Double_t dy1l = g1->GetEYlow() != NULL ? g1->GetEYlow()[ipoint] : (g1->GetEY() != NULL ? g1->GetEY()[ipoint] : 0.);
      Double_t dy1u = g1->GetEYhigh() != NULL ? g1->GetEYhigh()[ipoint] : dy1l;

      Double_t y2   = g2->GetY()[ipoint];
      Double_t dy2l = g2->GetEYlow() != NULL ? g2->GetEYlow()[ipoint] : (g2->GetEY() != NULL ? g2->GetEY()[ipoint] : 0.);
      Double_t dy2u = g2->GetEYhigh() != NULL ? g2->GetEYhigh()[ipoint] : dy2l;

      if( g1->GetX()[ipoint1] == g2->GetX()[ipoint2] ){ 
         ipoint1 += 1;
         ipoint2 += 1;
         g->SetPoint(ipoint, y1, y2);
         ipoint += 1; 
      }else if( g1->GetX()[ipoint1] < g2->GetX()[ipoint2] ){
         ipoint1 += 1;
      }else if( g1->GetX()[ipoint1] > g2->GetX()[ipoint2] ){  
         ipoint2 += 1; 
      }

      if(g->GetEY() != NULL){
         g->GetEX()[ipoint] = dy1l;
         g->GetEY()[ipoint] = dy2l;
      }else if (g->GetEYlow() != NULL){
         g->GetEXlow()[ipoint]  = dy1l;
         g->GetEXhigh()[ipoint] = dy1u;
         g->GetEYlow()[ipoint]  = dy2l;
         g->GetEYhigh()[ipoint] = dy2u;
      }

      if( ipoint1 == g1->GetN() && ipoint2 == g2->GetN() ) break;   

   }

   g->Set(ipoint);

   return g;
}
