//Display the effective day coverage of each NM BR-averaged bin

#include <TGraph.h>
#include <TFile.h>
#include <cstdio>
#include <TString.h>
#include "commonlib/include/Experiments.hh"
#include "commonlib/include/HistTools.hh"
#include "commonlib/include/Spline.hh"
#include "commonlib/include/DateTimeTools.hh"

using namespace std;

void useful(){

  //Count the number of effective 27-day bin 
  const int nNMs = 46; //# of NM stations
  const char *NM_name[nNMs] = { "PSNM", "TIBT", "DJON", "TSMB", "ATHN", "MXCO", "ARNM", "NANM", "PTFM", "CALM", "AATB", "ROME", "BKSN", "HRMS", "JUNG", "JUNG1", "LMKS", "IRK2", "IRK3", "IRKT", "DRBS", "MCRL", "MOSC", "NEWK", "KIEL", "KIEL2", "MGDN", "KERG", "OULU", "SANB", "SNAE", "APTY", "NRLK", "FSMT", "INVK", "JBGO", "MCMU", "NAIN", "PWNK", "THUL", "NEU3", "SOPB", "SOPO", "DOMB", "DOMC", "TERA" };
  int NM_eff27days_counts[nNMs] = { 73,42,25,56,67,70,10,63,18,56,64,49,69,74,74,74,49,10,17,10,39,63,60,72,32,70,34,69,75,20,59,75,9,73,65,16,68,73,70,70,28,75,74,14,14,67 };
  TH1F *h_eff27days_counts = new TH1F("h_eff27days_counts", "Effective 27-days Counts", nNMs, 0, nNMs); //create the useful bin count histogram

  //Fill the histogram, in a loop
  for(int i=0;i<=nNMs-1;i++){
  	h_eff27days_counts->SetBinContent(i+1, NM_eff27days_counts[i]);
  	h_eff27days_counts->GetXaxis()->SetBinLabel(i+1, NM_name[i]);
	printf("NM-%s %d \n", NM_name[i], NM_eff27days_counts[i]);
	//printf(".x getmean2.C(\"NM-%s\")\n", NM_name[i]);
  }

  h_eff27days_counts->GetXaxis()->LabelsOption("v"); // draw horizontal labels: if they are too long, they will be tilted
  h_eff27days_counts->SetFillColor(kBlue);
  h_eff27days_counts->Draw("HIST][");

}

