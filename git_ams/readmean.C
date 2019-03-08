//Draft

#include <TGraph.h>
#include <TFile.h>
#include <cstdio>
#include <TString.h>

#include "commonlib/include/Experiments.hh"
#include "commonlib/include/HistTools.hh"
#include "commonlib/include/Spline.hh"
#include "commonlib/include/DateTimeTools.hh"

using namespace std;

void readmean(const char*filename){

	const int nNMs = 46; //# of NM stations
  	const char *NM_name[nNMs] = { "PSNM", "TIBT", "DJON", "TSMB", "ATHN", "MXCO", "ARNM", "NANM", "PTFM", "CALM", "AATB", "ROME", "BKSN", "HRMS", "JUNG", "JUNG1", "LMKS", "IRK2", "IRK3", "IRKT", "DRBS", "MCRL", "MOSC", "NEWK", "KIEL", "KIEL2", "MGDN", "KERG", "OULU", "SANB", "SNAE", "APTY", "NRLK", "FSMT", "INVK", "JBGO", "MCMU", "NAIN", "PWNK", "THUL", "NEU3", "SOPB", "SOPO", "DOMB", "DOMC", "TERA" };

	TFile fin(Form("./data/nm/%s.root",filename));
	TGraph*g=(TGraph*) fin.Get("g");
	TGraph*g_ave=(TGraph*) fin.Get("g_ave");
	if(g_ave==NULL){
		printf("TGraph g_ave not found \n");
	}else if(g==NULL){
		printf("TGraph g not found \n");
	}else{
		g_ave->SetTitle(Form("%s Flux & AMS He Flux;Time;Flux [1/(m^2 sr s GV)]", filename));
		g_ave->GetXaxis()->SetTimeDisplay(1);
 		g_ave->GetXaxis()->SetTimeFormat("%m-%y");
 		g_ave->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00");
		g_ave->Draw("AP");
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

	g->Draw("LSame");

  	gft->SetTitle("AMS Flux in Time Domain;Time;Flux [1/(m^2 sr s GV)]");
  	gft->GetXaxis()->SetTimeDisplay(1);
  	gft->GetXaxis()->SetTimeFormat("%m-%y");
	gft->SetMarkerColor(kBlue);
	gft->SetMarkerStyle(kFullCircle);
  	gft->Draw("PSame");
  	//gft->Print();

	auto leg = new TLegend(0.1,0.7,0.3,0.9,"");
	leg->SetFillColor(0);
	leg->AddEntry(g, Form("%s Flux",filename),"l");
	leg->AddEntry(g_ave,Form("%s BR-averaged Flux",filename),"p");
	leg->AddEntry(gft, "AMS He Flux","p");
	leg->DrawClone("Same");

	TGraph *cor = HistTools::GetCorrelation(gft, g_ave, 0, g_ave->GetN());
	cor->Print();
	cor->SetMarkerColor(kBlue);
	cor->SetMarkerStyle(kFullCircle);
	cor->Draw("AP");

}



