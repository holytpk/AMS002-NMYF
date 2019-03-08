#include <TGraph.h>
#include <TFile.h>
#include <cstdio>
#include <TString.h>

#include "commonlib/include/Experiments.hh"
#include "commonlib/include/HistTools.hh"
#include "commonlib/include/Spline.hh"
#include "commonlib/include/DateTimeTools.hh"

using namespace std;

void readmean(const char*file1, const char*file2){

	const int nNMs = 46; //# of NM stations
  	const char *NM_name[nNMs] = { "PSNM", "TIBT", "DJON", "TSMB", "ATHN", "MXCO", "ARNM", "NANM", "PTFM", "CALM", "AATB", "ROME", "BKSN", "HRMS", "JUNG", "JUNG1", "LMKS", "IRK2", "IRK3", "IRKT", "DRBS", "MCRL", "MOSC", "NEWK", "KIEL", "KIEL2", "MGDN", "KERG", "OULU", "SANB", "SNAE", "APTY", "NRLK", "FSMT", "INVK", "JBGO", "MCMU", "NAIN", "PWNK", "THUL", "NEU3", "SOPB", "SOPO", "DOMB", "DOMC", "TERA" };

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
		g2_ave->Draw("PSame");
	}

	Experiments::DataPath = "data";
  	TH1D **hbr = Experiments::GetDatasetHistograms(Experiments::AMS02, 1);
  	int nfluxes = Experiments::Info[Experiments::AMS02].Dataset[1].nMeasurements;
  	TGraphErrors *gft = new TGraphErrors(nfluxes);
  	int bin = 1;

  	for (int iflux = 0; iflux < nfluxes; ++iflux) { 

		double flux = hbr[iflux]->GetBinContent(bin); 
		double err = hbr[iflux]->GetBinError(bin); 
		time_t *time_range = Experiments::GetMeasurementTimeRange(Experiments::AMS02, 1, iflux); 
		double t = 0.5*(time_range[0] + time_range[1]);
		double dt = time_range[1]-time_range[0];
		gft->SetPoint(iflux, t, flux);	
		gft->SetPointError(iflux, 0., err); 
		//printf("iflux=%d t=%10.0f flux=%f dt=%f \n",iflux,t,flux,dt); 
  
  	}

	g1->Draw("LSame");
	g2->Draw("LSame");

  	gft->SetTitle("AMS Flux in Time Domain;Time;Flux [1/(m^2 sr s GV)]");
  	gft->GetXaxis()->SetTimeDisplay(1);
  	gft->GetXaxis()->SetTimeFormat("%m-%y");
	gft->SetMarkerColor(kBlue);
	gft->SetMarkerStyle(kFullCircle);
  	gft->Draw("PSame");
  	//gft->Print();

	auto leg = new TLegend(0.1,0.7,0.3,0.9,"");
	leg->SetFillColor(0);
	leg->AddEntry(g1, Form("%s Flux",file1),"l");
	leg->AddEntry(g1_ave,Form("%s BR-averaged Flux",file1),"p");
	leg->AddEntry(gft, "AMS He Flux","p");
	leg->DrawClone("Same");
	
	//calculate correlation
	TGraph *cor = HistTools::GetCorrelation(gft, g1_ave, 0, g1_ave->GetN());
	cor->Print();
	cor->SetMarkerColor(kBlue);
	cor->SetMarkerStyle(kFullCircle);
	//cor->Draw("AP");

}



