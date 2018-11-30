#include <Experiments.hh>
#include <HistTools.hh>
#include <Spline.hh>
#include <debug.hh>
#include "commonlib/include/Experiments.hh"
#include "commonlib/include/HistTools.hh"
#include "commonlib/include/Spline.hh"
#include "commonlib/include/DateTimeTools.hh"

void Setup()
{
   Debug::Enable(Debug::ALL);

   // time axis style
   gSystem->Setenv("TZ", "UTC");
   gStyle->SetTimeOffset(0);

   // palette style
   gStyle->SetPalette(55);
   gStyle->SetNumberContours(99);

   gStyle->SetOptStat(0);

   // text style
   gStyle->SetLabelFont(43, "xyz");
   gStyle->SetLabelSize(26, "xyz");
   gStyle->SetLabelOffset(0.003, "xyz");
   gStyle->SetTitleFont(43, "xyz");
   gStyle->SetTitleSize(26, "xyz");
   gStyle->SetTitleOffset(1.2, "x");
   gStyle->SetTitleOffset(1., "y");

   // Ignore errors lower than the ignore level. Possible values:
   // kPrint, kInfo, kWarning, kError, kBreak, kSysError and kFatal
   // gErrorIgnoreLevel = kPrint;

   // fit options
   TVirtualFitter::SetDefaultFitter("Minuit2");
   TVirtualFitter::SetMaxIterations(5e4);
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
   ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(5e4);
   ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(5e4);
   ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.005);
   TMinuit minuit;
   minuit.SetPrintLevel(-1);
}

void test_fit()
{
   Setup();
   Experiments::DataPath = "data";

   TH1D *h = Experiments::GetMeasurementHistogram(Experiments::AMS02, 0, 0);
   h->SetMarkerStyle(kFullCircle);
   h->SetMarkerSize(1.7);

   Spline *sp = Spline::BuildFromHistogram(h, "fsp", 7, Spline::LogLog | Spline::PowerLaw);
   TF1 *fsp = sp->GetTF1Pointer();
   h->Fit(fsp, "N");
   //h->Fit(fsp, "NI");

   TH1 *hr = HistTools::GetResiduals(h, fsp, "_fitratio", false, true, true, 4, 1);
   HistTools::CopyStyle(h, hr);

   TH1 *he = HistTools::GetFitError(h, fsp, "_fiterr", true, false, true, 11, 1.);
   HistTools::SetFillStyle(he, kRed-7, 1001);
   PRINT_HIST(he)

   TCanvas *c = new TCanvas();
   c->Divide(1, 2);

   c->cd(1);
   h->Draw("E1X0");
   fsp->Draw("SAME");
   gPad->SetLogx();
   gPad->SetLogy();
   gPad->SetGrid();

   c->cd(2);
   he->SetTitle("");
   he->SetXTitle("Rigidity [GV]");
   he->SetYTitle("Data / Fit");
   he->Draw("E3");
   he->GetYaxis()->SetRangeUser(0.9, 1.1);
   gPad->SetLogx();
   gPad->SetGrid();
   hr->Draw("E1X0 SAME");

   c->Print("c.png");
}
