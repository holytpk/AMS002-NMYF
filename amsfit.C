#include "commonlib/include/Experiments.hh"
#include "commonlib/include/HistTools.hh"
#include "commonlib/include/Spline.hh"

void amsfit(){

  auto c1 = new TCanvas("c1","AMS Data",1600,800);
  //c1->Divide(2,2);

  float widths[] = { 1., 1. };
  float heights[] = { 3., 1. };
  float margins[] = { 0.07, 0.02, 0.08, 0.08 };
  HistTools::DivideCanvas(c1, 2, 2, widths, heights, margins, 0.05, 0.02);
  
  c1->cd(1);
  Experiments::DataPath = "data";
  TH1 *h = Experiments::GetMeasurementHistogram(Experiments::AMS02, 0, 0);
  h->GetXaxis()->SetLabelFont(gStyle->GetLabelFont("x"));
  h->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("x"));
  h->GetYaxis()->SetLabelFont(gStyle->GetLabelFont("y"));
  h->GetYaxis()->SetLabelSize(gStyle->GetLabelSize("y"));
  h->GetXaxis()->SetTitleFont(gStyle->GetTitleFont("x"));
  h->GetXaxis()->SetTitleSize(gStyle->GetTitleSize("x"));
  h->GetYaxis()->SetTitleFont(gStyle->GetTitleFont("y"));
  h->GetYaxis()->SetTitleSize(gStyle->GetTitleSize("y"));
  HistTools::SetStyle(h, kBlue, kFullCircle, 1.7, 1, 1);
  h->Draw("E1X0"); gPad->SetLogx(); gPad->SetLogy(); gPad->SetGrid();
  h->GetYaxis()->SetTitleOffset(2);
  
  
  /*TF1 *f_fit = new TF1("f_fit", "[0]*pow(x,[1])", 1, 2000);
  f_fit->SetParameters(1e4, -3); 
  h->Fit(f_fit, "N", "", 100, 2000);
  f_fit->Draw("SAME");*/
  
  
  Spline *sp = Spline::BuildFromHistogram(h, "sp", 7, Spline::LogX | Spline::LogY | Spline::PowerLaw);
  TF1 *fsp = sp->GetTF1Pointer();
  //HistTools::PrintFunction(fsp);
  h->Fit(fsp, "N");
  fsp->Draw("SAME");  

  c1->cd(2);
  //Execute for all files
  TH1D **hbr = Experiments::GetDatasetHistograms(Experiments::AMS02, 1);
  HistTools::SetColors(hbr, 79, kFullCircle, 1.4);
  for (int i=1; i<79; ++i) { hbr[i]->Draw("E1X0 SAME"); gPad->SetLogx(); gPad->SetLogy(); gPad->SetGrid(); }

  //F(R) for another set of data
/*  TH1 *h_2 = Experiments::GetMeasurementHistogram(Experiments::AMS02, 1, 1);
  h_2->SetMarkerStyle(kFullCircle);
  h_2->SetMarkerSize(1.7);
  h_2->SetMarkerColor(kRed);
  h_2->Draw("E1X0"); gPad->SetLogx(); gPad->SetLogy(); gPad->SetGrid();

  Spline *sp_2 = Spline::BuildFromHistogram(h_2, "sp", 7, Spline::LogX | Spline::LogY | Spline::PowerLaw);
  TF1 *fsp_2 = sp_2->GetTF1Pointer();
  h_2->Fit(fsp_2, "N");
  fsp_2->Draw("SAME"); 
*/

  c1->cd(3);
  TH1 *h_ratio = HistTools::GetResiduals(h, fsp, "_ratio", false, 0, true, 5, 1);
  HistTools::CopyStyle(h,h_ratio);
  h_ratio->GetXaxis()->SetLabelFont(gStyle->GetLabelFont("x"));
  h_ratio->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("x"));
  h_ratio->GetYaxis()->SetLabelFont(gStyle->GetLabelFont("y"));
  h_ratio->GetYaxis()->SetLabelSize(gStyle->GetLabelSize("y"));
  h_ratio->GetXaxis()->SetTitleFont(gStyle->GetTitleFont("x"));
  h_ratio->GetXaxis()->SetTitleSize(gStyle->GetTitleSize("x"));
  h_ratio->GetYaxis()->SetTitleFont(gStyle->GetTitleFont("y"));
  h_ratio->GetYaxis()->SetTitleSize(gStyle->GetTitleSize("y"));
  h_ratio->Draw("E1X0"); gPad->SetLogx(); gPad->SetGrid();
  h_ratio->SetTitle("");
  h_ratio->GetYaxis()->SetRangeUser(-0.12,0.12);
  h_ratio->GetYaxis()->SetTitle("Fit/Data -1");
  h_ratio->GetYaxis()->SetTitleOffset(2);
  h_ratio->GetXaxis()->SetTitleOffset(3.5);
  h_ratio->GetXaxis()->SetNoExponent(true);

  c1->cd(4);
/*  TH1 *h_ratio_2 = HistTools::GetResiduals(h_2, fsp_2, "_ratio", false, 0, true, 5, 1);
  HistTools::CopyStyle(h_2,h_ratio_2);
  h_ratio_2->GetXaxis()->SetLabelFont(gStyle->GetLabelFont("x"));
  h_ratio_2->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("x"));
  h_ratio_2->GetYaxis()->SetLabelFont(gStyle->GetLabelFont("y"));
  h_ratio_2->GetYaxis()->SetLabelSize(gStyle->GetLabelSize("y"));
  h_ratio_2->GetXaxis()->SetTitleFont(gStyle->GetTitleFont("x"));
  h_ratio_2->GetXaxis()->SetTitleSize(gStyle->GetTitleSize("x"));
  h_ratio_2->GetYaxis()->SetTitleFont(gStyle->GetTitleFont("y"));
  h_ratio_2->GetYaxis()->SetTitleSize(gStyle->GetTitleSize("y"));
  h_ratio_2->Draw("E1X0"); gPad->SetLogx(); gPad->SetGrid();
  h_ratio_2->SetTitle("");
  h_ratio_2->GetYaxis()->SetRangeUser(-0.12,0.12);
  h_ratio_2->GetYaxis()->SetTitle("Fit/Data -1");
  h_ratio_2->GetYaxis()->SetTitleOffset(2);
  h_ratio_2->GetXaxis()->SetTitleOffset(3.5);
  h_ratio_2->GetXaxis()->SetNoExponent(true);
*/

  TFile fout("spfit.root","recreate");
  fsp->Write();
  fout.Close();

}


