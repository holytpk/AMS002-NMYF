int tutorial(){

//Load 
#include "commonlib/include/Experiments.hh"
#include "commonlib/include/HistTools.hh"
#include "commonlib/include/Spline.hh"
#include "commonlib/include/DateTimeTools.hh"
Experiments::DataPath = "data";

//Check experiment info
Experiments::Info[Experiments::AMS02].Print();
Experiments::Info[Experiments::AMS02].Dataset[1].nMeasurements;

//Draw F(R) vs. R
TH1 *h = Experiments::GetMeasurementHistogram(Experiments::AMS02, 0, 0);
HistTools::SetStyle(h, kBlue, kFullCircle, 1.7, 1, 1);
h->Draw("E1X0"); gPad->SetLogx(); gPad->SetLogy(); gPad->SetGrid();
TH1D **hbr = Experiments::GetDatasetHistograms(Experiments::AMS02, 1);
HistTools::SetColors(hbr, 79, kFullCircle, 1.4);
for (int i=1; i<79; ++i) { hbr[i]->Draw("E1X0 SAME"); gPad->SetLogx(); gPad->SetLogy(); gPad->SetGrid(); }

//Get Spectral Index from TCanvas
gPad->GetListOfPrimitives()->GetEntries()
gPad->GetListOfPrimitives()->At(0)->ClassName()
gPad->GetListOfPrimitives()->At(1)->ClassName()
gPad->GetListOfPrimitives()->At(2)->ClassName()
TH1 *h = (TH1 *)gPad->GetListOfPrimitives()->At(1)
h->GetListOfFunctions()
h->GetListOfFunctions()->GetEntries()

TH1 *h = (TH1 *)gPad->GetListOfPrimitives()->At(1)
TGraphAsymmErrors *gsp = HistTools::GetSpectralIndex(h)
new TCanvas
gsp->Draw("AP"); gPad->SetLogx(); gPad->SetGrid()
Spline *sp = Spline::BuildFromHistogram(h, "sp", 7, Spline::LogX | Spline::LogY | Spline::PowerLaw)
TF1 *fsp = sp->GetTF1Pointer()
HistTools::PrintFunction(fsp)
h->Fit(fsp, "N")
c1->cd(1)
fsp->Draw("same")

c1->cd(1)
TH1 *h = (TH1 *)gPad->GetListOfPrimitives()->At(1)
TGraphAsymmErrors *gsp = HistTools::GetSpectralIndex(h)
new TCanvas
gsp->Draw("AP"); gPad->SetLogx(); gPad->SetGrid()
Spline *sp = Spline::BuildFromHistogram(h, "sp", 7, Spline::LogX | Spline::LogY | Spline::PowerLaw)
TF1 *fsp = sp->GetTF1Pointer()
HistTools::PrintFunction(fsp)
h->Fit(fsp, "N")
c1->cd(1)
fsp->Draw("same")
TF1 *fsp_sp = HistTools::GetSpectralIndex(fsp, "fsp_sp", 1, 2000)
new TCanvas
gsp->Draw("E1X0")
gsp->Draw("AP"); gPad->SetLogx(); gPad->SetGrid()
fsp_sp->Draw("same")
HistTools::PrintFunction(fsp)
TMath::Prob(4.17, 63)
fsp->GetChisquare()
fsp->GetNDF()
fsp->GetParameter(0)
fsp->GetParError(0)
TH1 *h_ratio = HistTools::GetResiduals(h, fsp, "_ratio", false, 0, true, 5, 1)
TH1 *h_res = HistTools::GetResiduals(h, fsp, "_res", false, 0, true, 1, 0)
TCanvas *cfitres = new TCanvas("cfitres")
cfitres->Divide(1,2)
cfitres->cd(1)
h_ratio->Draw("E1x0")
gPad->SetLogx(); gPad->SetGrid()
cfitres->cd(2)
h_res->Draw("HIST")
gPad->SetLogx(); gPad->SetGrid()
TH1 *h_fiterr = HistTools::GetFitError(h, fsp, "_fiterr", true, 0, true, 11, 0.)
HistTools::SetFillStyle(h_fiterr, kRed, 3001)
cfitres->cd(1)
h_fiterr->Draw("E3SAME")
h_fiterr->SetLineAttributes()
h_fiterr->GetNbinsX()
TH1 *h_fiterr = HistTools::GetFitError(h, fsp, "_fiterr", false, 0, true, 11, 0.)
h_fiterr = HistTools::GetFitError(h, fsp, "_fiterr", false, 0, true, 11, 0.)
h_fiterr->GetBinContent(1)
h_fiterr->GetBinContent(2)
h_fiterr->GetNbinsX()
h_fiterr->GetBinError(1)
h_fiterr->GetBinError(10)

//Spine Tutorial
TH1 *h = Experiments::GetMeasurementHistogram(Experiments::AMS02, 0, 0);
.L test_spline.C 
test_spline(h, 7, 2.7)
test_spline(h, 7, 1.5)
test_spline(h, 7, 0.001)
gStyle->GetLabelFont("x")
gStyle->GetLabelFont("y")
gStyle->GetTitleFont("y")

//BR flux to time domain tutorial
TH1D **hbr = Experiments::GetDatasetHistograms(Experiments::AMS02, 1);
int nfluxes = Experiments::Info[Experiments::AMS02].Dataset[1].nMeasurements;
TGraphErrors *gft = new TGraphErrors(nfluxes);
int bin = 1;
for (int iflux = 0; iflux < nfluxes; ++iflux) { double flux = hbr[iflux]->GetBinContent(bin); double err = hbr[iflux]->GetBinError(bin); time_t *time_range = Experiments::GetMeasurementTimeRange(Experiments::AMS02, 1, iflux); double t = 0.5*(time_range[0] + time_range[1]); gft->SetPoint(iflux, t, flux); gft->SetPointError(iflux, 0., err); };
gft->GetXaxis()->SetTimeDisplay(1);
gft->GetXaxis()->SetTimeFormat("%m-%y");
gft->Draw("APSAME");
gft->Print();
HistTools::SetStyle(gft, kBlack, kFullCircle, 1.7, 1, 1)
bin = hbr[0]->FindBin(5); for (int iflux = 0; iflux < nfluxes; ++iflux) { double flux = hbr[iflux]->GetBinContent(bin); double err = hbr[iflux]->GetBinError(bin); time_t *time_range = Experiments::GetMeasurementTimeRange(Experiments::AMS02, 1, iflux); double t = 0.5*(time_range[0] + time_range[1]); gft->SetPoint(iflux, t, flux); gft->SetPointError(iflux, 0., err); }
time_t begin_date = DateTimeTools::UIntToUTCTime(20110101)
time_t end_date = DateTimeTools::UIntToUTCTime(20180101)
TH1 *ha = HistTools::CreateTimeAxis("ha", "AMS-02 monthly proton flux", begin_date, end_date, 86400, 100, 150)
ha->Draw("SAME")
ha->GetXaxis()->SetNdivisions(7, 12, 0, false)
gft->Draw("PSAME")
gPad->Modified()
ha->GetXaxis()->SetTimeFormat("%Y"); gPad->Modified()
ha->GetXaxis()->SetTimeFormat("%Y/%m"); gPad->Modified()
ha->GetXaxis()->SetTimeFormat("%Y/%m/%d"); gPad->Modified()

//Multiply Histograms
TF1 *f1 = new TF1("f1", "sin(x)", 0, 2*TMath::Pi())
TF1 *f2 = new TF1("f2", "x", 0, 2*TMath:Pi())
TF1 *f2 = new TF1("f2", "x", 0, 2*TMath::Pi())
f2->SetLineColor(kBlue)
f1->Draw()
f2->Draw("same")
TF1 *f3 = new TF1("f3", "f1*f2", 0, 2*TMath::Pi())
f3->SetLineColor(kGreen+2)
f3->Draw("same")
TF1 *f4 = HistTools::CombineTF1(f1, f2, HistTools::Product, "f4", 0, 2*TMath::Pi())
TF1 *f4 = HistTools::CombineTF1(f1, f2, HistTools::Multiply, "f4", 0, 2*TMath::Pi())
HistTools::SetLineStyle(f4, kMagenta, 7, 3)
f4->Draw("same")

//Spline Correction by Recreation
#include "commonlib/include/Spline.hh"
int N=(sp->GetNpar()-2)/2; Spline *s1 = new Spline("fsp", N, Spline::LogLog | Spline::PowerLaw, sp->GetParameters()+2, sp->GetParameters()+2+N, sp->GetParameter(0), sp->GetParameter(1));
TF1 *fsp = s1->GetTF1Pointer();
HistTools::SetLineStyle(fsp, kBlue, 7, 3);
fsp->Draw("same"); gPad->SetGrid();



}
