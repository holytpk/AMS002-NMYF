//Calculate the average count rate for each Bartel rotation
//Match NM data point with AMS data point
//Count good quality NM bins to ensure the usefulness of data
//Updated 11/15/2018

#include <TGraph.h>
#include <TFile.h>
#include <cstdio>
#include <TString.h>
#include "commonlib/include/Experiments.hh"
#include "commonlib/include/HistTools.hh"
#include "commonlib/include/Spline.hh"
#include "commonlib/include/DateTimeTools.hh"

using namespace std;

void getmean2(const char*filename){

  Experiments::DataPath = "data";
  TGraph *g = new TGraph(Form("./data/nm/%s.dat", filename));
  TGraph *g_ave = new TGraph(g->GetN()/27+(g->GetN()%27>0)); //use less memory 
  double x, y, x_ave=0, y_ave=0, t_0=0, t_1=0, t;
  int j=0, k=0, count_x=0; //initiate counters, j is BR index
  time_t *tran = Experiments::GetMeasurementTimeRange(Experiments::AMS02, 1, j); 

  TH1F *h_useful = new TH1F("h_useful", "Number of good days per BR", 27, 0, 27); //count useful quality NM bins

  time_t GLE71_date_1 = DateTimeTools::UIntToUTCTime(20120517); //convert a date given in yyyymmdd format into a Unix time
  time_t GLE71_date_2 = DateTimeTools::UIntToUTCTime(20120518);
  time_t SUBGLE1_date = DateTimeTools::UIntToUTCTime(20120127);
  time_t SUBGLE2_date = DateTimeTools::UIntToUTCTime(20120307);
  time_t SUBGLE3_date = DateTimeTools::UIntToUTCTime(20140106);

  const int nNMs = 46; //# of NM stations
  const char *NM_name[nNMs] = { "PSNM", "TIBT", "DJON", "TSMB", "ATHN", "MXCO", "ARNM", "NANM", "PTFM", "CALM", "AATB", "ROME", "BKSN", "HRMS", "JUNG", "JUNG1", "LMKS", "IRK2", "IRK3", "IRKT", "DRBS", "MCRL", "MOSC", "NEWK", "KIEL", "KIEL2", "MGDN", "KERG", "OULU", "SANB", "SNAE", "APTY", "NRLK", "FSMT", "INVK", "JBGO", "MCMU", "NAIN", "PWNK", "THUL", "NEU3", "SOPB", "SOPO", "DOMB", "DOMC", "TERA" };

  //cout <<"GLEs " << DateTimeTools::UIntToUTCTime(20120517) <<" "<< DateTimeTools::UIntToUTCTime(20120127) << endl;

  auto c1 = new TCanvas("c1","NM BR-averaged Flux & Useful NM Bin Histogram",1600,900);
  c1->Divide(2,1);

  	for(int i=0;i<g->GetN();i++){ 
       
		//i is NM data point iterations, k is AMS BR period iteration
		g->GetPoint(i,x,y); 
		
		tran = Experiments::GetMeasurementTimeRange(Experiments::AMS02, 1, k);
 		t_0 = tran[0];
  		t_1 = tran[1];
		t   = 0.5*(t_1-t_0); //half of AMS time range

		//printf("x=%10.0f t_0=%f t_1=%f k=%d \n", x, t_0, t_1, k);

		if (x<t_0){
			
			continue; //skip data outside AMS time range			

		}else if (x>=t_0 && x<t_1){
			//check if inside AMS time range 

			

			if (y>0){

				if (x == GLE71_date_1 || x == GLE71_date_2){
	
					continue; //exclude ground level enhancement from data

				}else if (x == SUBGLE1_date || x == SUBGLE2_date || x == SUBGLE3_date){
	
					continue; //exclude sub ground level enhancement from data

				}	

				x_ave += x; 
				y_ave += y;
				count_x++;
			}	

			//printf("i=%d x=%10.0f y=%f y_ave=%f t_0=%f t_1=%f \n", i, x, y, y_ave,t_0,t_1);
  		
		}else if (x==t_1 && y_ave>0){

				//printf("x=%10.0f t_0=%f t_1=%f \n", x, t_0, t_1);

				y_ave /= count_x;
				x_ave /= count_x; 
				//printf("x=t_1 | count_x=%d x=%10.0f x_ave=%10.0f y_ave=%f j=%d k=%d \n", count_x, x, x_ave, y_ave, j, k);

				g_ave->SetPoint(j++,x-t,y_ave); //Input average values 
				
				h_useful->Fill(count_x); //Input number of useful days in each bin into histogram

				x_ave=0;
				y_ave=0;
				count_x=0;
  				//end of the BR loop, calculate the mean sigal 

		}else if (x>t_1){
			if (k<=Experiments::Info[Experiments::AMS02].Dataset[1].nMeasurements-2){
				k++;
				i--;
				//printf("k++ | i=%d x=%10.0f y=%f y_ave=%f t_0=%f t_1=%f \n", i, x, y, y_ave,t_0,t_1);
				//cout << "k++ | j=" << j << " " << "k=" << k << endl;
			} //check if k exceed AMS data limit 
		}
			
  	} //close iteration of each data point

  	/*if (g->GetN()%27>0 && y_ave>0){
		y_ave /= count_x;
		x_ave /= g->GetN()%27;
		printf("count_x=%d x_ave=%10.0f y_ave=%f j=%d remaining days=%d \n", count_x, x_ave, y_ave, j, g->GetN()%27);
		g_ave->SetPoint(j++,x_ave,y_ave);
  	} */ //calculate mean of remaining points

  g_ave->Set(j);
  //gStyle->SetLabelSize(0.05);

  //Display NM + NM BR_Average 
  c1->cd(1);
  g->GetXaxis()->SetTimeDisplay(1);
  g->GetXaxis()->SetTimeFormat("%m-%y");
  g->GetXaxis()->SetTimeOffset(0,"1970-01-01 00:00:00");
  gPad->SetLeftMargin(0.11);  
  gPad->SetRightMargin(0.05);
  g->SetTitle(Form("%s Observation with BR-average;Time;Flux [1/(m^2 sr s GV)]", filename));
  g->GetYaxis()->SetTitleOffset(1.65);
  g->Draw("ALSAME");
  g_ave->SetMarkerColor(kRed);
  g_ave->SetMarkerStyle(kFullCircle);
  g_ave->Draw("PSAME");
  g_ave->Print();

  //Display Useful Bin Statistics
  c1->cd(2);
  gPad->SetLeftMargin(0.069);
  gPad->SetRightMargin(0.06);
  h_useful->GetYaxis()->SetTitleOffset(1.1);
  h_useful->SetTitle(Form("%s Useful Bin Statistics; # of Useful Days in a Bin; # of Bins", filename));
  h_useful->SetFillColor(kBlue);
  h_useful->Draw("HIST][");

  Int_t bin = h_useful->GetBinContent(27); //count effective bin numbers
  printf("%s Effective Bin Counts = %d \n", filename, bin);

  //Store plots in root file
  TFile fout(Form("./data/nm/%s.root", filename), "recreate");
  g->Write("g");
  g_ave->Write("g_ave");
  h_useful->Write("h_useful");
  fout.Close();

  c1->SaveAs(Form("./data/useful/%s.png",filename));

  return 0;
}



