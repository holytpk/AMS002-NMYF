#include <TGraph.h>
#include <TFile.h>
#include <cstdio>
#include <TString.h>
#include "commonlib/include/Experiments.hh"
#include "commonlib/include/HistTools.hh"
#include "commonlib/include/Spline.hh"
#include "commonlib/include/DateTimeTools.hh"

using namespace std;

void getcorr(const char*file1,const char*file2){

  Experiments::DataPath = "data";
  const int nNMs = 46; //# of NM stations
  	const char *NM_name[nNMs] = { "PSNM", "TIBT", "DJON", "TSMB", "ATHN", "MXCO", "ARNM", "NANM", "PTFM", "CALM", "AATB", "ROME", "BKSN", "HRMS", "JUNG", "JUNG1", "LMKS", "IRK2", "IRK3", "IRKT", "DRBS", "MCRL", "MOSC", "NEWK", "KIEL", "KIEL2", "MGDN", "KERG", "OULU", "SANB", "SNAE", "APTY", "NRLK", "FSMT", "INVK", "JBGO", "MCMU", "NAIN", "PWNK", "THUL", "NEU3", "SOPB", "SOPO", "DOMB", "DOMC", "TERA" };

	auto c1 = new TCanvas("c1","NM BR-averaged Fluxes & Correlation",1600,900);
  	c1->Divide(2,1);

	c1->cd(1);

	//read in first file
	TFile fin1(Form("./data/nm/%s.root",file1));
	TGraph*g1=(TGraph*) fin1.Get("g");
	TGraph*g1_ave=(TGraph*) fin1.Get("g_ave");
	if(g1_ave==NULL){
		printf("TGraph g1_ave not found \n");
	}else if(g1==NULL){
		printf("TGraph g1 not found \n");
	}else{
		g1_ave->SetTitle(Form("%s Flux & AMS He Flux;Time;Flux [1/(m^2 sr s GV)]", file1));
		g1_ave->GetXaxis()->SetTimeDisplay(1);
 		g1_ave->GetXaxis()->SetTimeFormat("%m-%y");
 		g1_ave->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00");
		g1_ave->SetMarkerColor(kBlue);
		g1_ave->Draw("AP");
	}

	//read in second file
	TFile fin2(Form("./data/nm/%s.root",file2));
	TGraph*g2=(TGraph*) fin2.Get("g");
	TGraph*g2_ave=(TGraph*) fin2.Get("g_ave");
	if(g2_ave==NULL){
		printf("TGraph g2_ave not found \n");
	}else if(g2==NULL){
		printf("TGraph g2 not found \n");
	}else{
		g2_ave->SetTitle(Form("%s Flux & AMS He Flux;Time;Flux [1/(m^2 sr s GV)]", file1));
		g2_ave->GetXaxis()->SetTimeDisplay(1);
 		g2_ave->GetXaxis()->SetTimeFormat("%m-%y");
 		g2_ave->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00");
		g2_ave->SetMarkerColor(kRed);
		g2_ave->Draw("PSame");
	}

	g1->Draw("LSame");
	g2->Draw("LSame");

	c1->cd(2);

	TGraph *cor = HistTools::GetCorrelation(g1_ave, g2_ave, 0, (g1_ave->GetN() > (g2_ave->GetN()) ? g1_ave->GetN() : g2_ave->GetN()) );
	cor->Print();
	cor->SetMarkerColor(kBlue);
	cor->SetMarkerStyle(kFullCircle);

}
