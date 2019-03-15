//modified on March 15, 2019
//replaced "data" w/ "ace"
//defined EMed[] globally based on the element 
//replaced EMin[], EMax[] w/ h->GetBinLowEdge(ibin+1), h->GetBinLowEdge(ibin+2)
//added condition ibin%2==0 to match the bins
//replaced ibin w/ ibin/2 to match the 7-bin arrays w/ 14-bin arrays


TGraphAsymmErrors *get_ace_average_graph(UInt_t *BRs, UInt_t nBRs)
{
   const UShort_t nbins = 14;
   TGraphAsymmErrors *graph = new TGraphAsymmErrors(nbins); 

   for (UShort_t ibin = 0; ibin < nbins; ++ibin)
   {
	if (ibin%2==0){
      		graph->SetPoint(ibin, EMed[ibin/2], 0.);
      		graph->SetPointError(ibin, EMed[ibin/2] - h->GetBinLowEdge(ibin+1), h->GetBinLowEdge(ibin+2) - EMed[ibin/2], 0., 0.);
	}
   }

   Double_t dT = 0.;
   for (UInt_t ibr = 0; ibr < nBRs; ++ibr)
   {
      if (BRs[ibr] >= ace->GetEntries()) continue;

      ace->GetEntry(BRs[ibr]);
      dT += livetime;

      for (UShort_t ibin = 0; ibin < nbins; ++ibin)
      {
	if (ibin%2==0){
         	graph->GetY()[ibin] += livetime*F[ibin/2];
         	graph->GetEYlow()[ibin] += C[ibin/2];
	}
      }
   }
   for (UShort_t ibin = 0; ibin < nbins; ++ibin)
   {
	if (ibin%2==0){
      		Double_t flux   = graph->GetY()[ibin]/dT;
      		Double_t counts = graph->GetEYlow()[ibin];
      		Double_t syst   = flux * sqrt(8e-4 + SpallCorrUnc[ibin/2]*SpallCorrUnc[ibin/2]/SpallCorr[ibin/2]/SpallCorr[ibin/2]);
      		Double_t stat   = flux/sqrt(counts);
      		Double_t err    = sqrt(stat*stat + syst*syst);

      		graph->GetY()[ibin]      = flux;
      		graph->GetEYlow()[ibin]  = err;
      		graph->GetEYhigh()[ibin] = err;
	} 
   }

   return graph;
}
