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

void getcorr(const char*file1, const char*file2){

  	Experiments::DataPath = "data";
  	const int nNMs = 46; //# of NM stations
  	const char *NM_name[nNMs] = { "PSNM", "TIBT", "DJON", "TSMB", "ATHN", "MXCO", "ARNM", "NANM", "PTFM", "CALM", "AATB", "ROME", "BKSN", "HRMS", "JUNG", "JUNG1", "LMKS", "IRK2", "IRK3", "IRKT", "DRBS", "MCRL", "MOSC", "NEWK", "KIEL", "KIEL2", "MGDN", "KERG", "OULU", "SANB", "SNAE", "APTY", "NRLK", "FSMT", "INVK", "JBGO", "MCMU", "NAIN", "PWNK", "THUL", "NEU3", "SOPB", "SOPO", "DOMB", "DOMC", "TERA" };

	const int nNMs_useful = 15;
	const char *NM_useful[nNMs_useful] = { "PSNM", "MXCO", "HRMS", "JUNG", "JUNG1", "NEWK", "KIEL2", "OULU", "APTY", "FSMT", "NAIN", "PWNK", "THUL", "SOPB", "SOPO" };

	auto c1 = new TCanvas("c1","NM BR-averaged Fluxes & Correlation",1600,900);
  	c1->Divide(2,1);

	c1->cd(1);

	//read in first file
	TFile fin1(Form("./data/nm/%s.root",file1));
	TGraph *g1=(TGraph*) fin1.Get("g");
	TGraph *g1_ave=(TGraph*) fin1.Get("g_ave");
	if(g1_ave==NULL){
		printf("TGraph g1_ave not found \n");
	}else if(g1==NULL){
		printf("TGraph g1 not found \n");
	}else{
		g1_ave->SetTitle(Form("%s Flux & %s Flux;Time;Flux [1/(m^2 sr s GV)]", file1,file2));
		g1_ave->GetXaxis()->SetTimeDisplay(1);
 		g1_ave->GetXaxis()->SetTimeFormat("%m-%y");
 		g1_ave->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00");
		g1_ave->SetMarkerColor(kBlue);
		g1_ave->Draw("AP");
	}

	//read in second file
	TFile fin2(Form("./data/nm/%s.root",file2));
	TGraph *g2=(TGraph*) fin2.Get("g");
	TGraph *g2_ave=(TGraph*) fin2.Get("g_ave");
	if(g2_ave==NULL){
		printf("TGraph g2_ave not found \n");
	}else if(g2==NULL){
		printf("TGraph g2 not found \n");
	}else{
		g2_ave->SetTitle(Form("%s Flux & %s Flux;Time;Flux [1/(m^2 sr s GV)]", file1,file2));
		g2_ave->GetXaxis()->SetTimeDisplay(1);
 		g2_ave->GetXaxis()->SetTimeFormat("%m-%y");
 		g2_ave->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00");
		g2_ave->SetMarkerColor(kRed);
		g2_ave->Draw("PSame");
	}

	//for (int i=0;

	g1_ave->Print();
	g2_ave->Print();

	g1->Draw("LSame");
	g2->Draw("LSame");

	c1->cd(2);

	TGraph *cor = GetCorrelation2(g1_ave, g2_ave);
	//TGraph *cor = GetCorrelation2(g1_ave, g2_ave, 0, (g1_ave->GetN() > g2_ave->GetN()) ? g1_ave->GetN() : g2_ave->GetN() );
	cor->SetMarkerColor(kBlue);
	cor->SetMarkerStyle(kFullCircle);
	cor->Print();
	cor->Draw("AP");
	printf("lastpoint in g1_ave or g2_ave = %d \n", (g1_ave->GetN() > g2_ave->GetN()) ? g1_ave->GetN() : g2_ave->GetN() );

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
