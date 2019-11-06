// barEffSmearing.C
// Simulates the effect of time resolution on the efficiency correction used
const std::vector<double> timeRes  = {1, 1.2, 1.4, 1.6, 1.8, 2.}; // in ns
const std::vector<double> effSigma = {0.2, 0.25, 0.3, 0.33, 0.35, 0.38, 0.4}; // one sigma values for efficiency curves
const double timeToDist = 0.0768; // Factor to go from time to position in m

void barEffSmearing(const char* outFile, const int nEntries=10000000)
{
  TRandom3 *r = new TRandom3(666);

  std::vector<TH1D*> effPlots;
  for (int e=0; e<effSigma.size(); e++) {
    TH1D *hEff = new TH1D(Form("hEff%d", (int)(effSigma[e]*100)), Form("Efficiency without smearing along S4 bar, %2g m sigma; Bar position / m; Events", effSigma[e]), 20, 0, 1.4);
    effPlots.push_back(hEff);
  }
  std::vector< std::vector<TH1D*> > resPlots;
  for (int e=0; e<effPlots.size(); e++) {
    std::vector<TH1D*> resPlotsSmall;
    for (int h=0; h<timeRes.size(); h++) {
      TH1D *hTimeResEff = new TH1D(Form("hTimeResEff%dRes%d", (int)(effSigma[e]*100.), (int)(timeRes[h]*10.)), Form("Efficiency with %.2g ns smearing along S4 bar; Bar position / m; Events", timeRes[h]), 20, 0, 1.4);
      hTimeResEff->Sumw2();
      resPlotsSmall.push_back(hTimeResEff);
    }
    resPlots.push_back(resPlotsSmall);
  }
  // Number of histogram entries to do
  for (int i=0; i<nEntries; i++) {
    if (i % 100000 == 0) std::cout<<i<<" of "<<nEntries<<std::endl;
    for (int e=0; e<effPlots.size(); e++) {
      double pos = r->Gaus(0.7, effSigma[e]);
      effPlots[e]->Fill(pos); // Histogram with no smearing
      // Loop over resolutions
      for (int h=0; h<timeRes.size(); h++) {
	resPlots[e][h]->Fill(r->Gaus(pos, timeRes[h]*timeToDist));
      } // Loop over resolutions
    } // Loop over efficiency sigmas
  } // Loop over the entries

  TFile *fout = new TFile(outFile, "recreate");
  fout->cd();
  for (int e=0; e<effSigma.size(); e++) {
    effPlots[e]->Write();
    THStack *hsRatioEff = new THStack(Form("hsRatioEff%d", (int)(effSigma[e]*100)), Form("Smeared/unsmeared effs, efficiency sigma %.2g m; Bar position / m; Smeared/unsmeared", effSigma[e]));
    TLegend *legEff = new TLegend(0.3, 0.3, 0.7, 0.6);
    for (int h=0; h<timeRes.size(); h++) {
      resPlots[e][h]->SetLineColor(kBlack);
      resPlots[e][h]->SetLineWidth(2);
      resPlots[e][h]->Write();
      TH1D *hRatio = new TH1D(Form("hRatioEff%dRes%d", (int)(effSigma[e]*100.), (int)(timeRes[h]*10.)), Form("Efficiency ratio with %.2g ns smearing along S4 bar and efficiency sigma of %.2g m; Bar position / m; Smeared/unsmeared", timeRes[h], effSigma[e]), 20, 0, 1.4);
      hRatio->Divide(resPlots[e][h], effPlots[e]);
      hRatio->Write();
      hRatio->SetLineColor(52+h*4);
      hsRatioEff->Add(hRatio);
      legEff->AddEntry(hRatio, Form("Time res = %.2g", timeRes[h]), "l");

    }
    hsRatioEff->Write();
    if (e==0) legEff->Write("legEff");
  }

  for (int h=0; h<timeRes.size(); h++) {
    THStack *hsRatioRes = new THStack(Form("hsRatioRes%d", (int)(timeRes[h]*10)), Form("Smeared/unsmeared effs, time resolution %.2g ns; Bar position / m; Smeared/unsmeared", timeRes[h]));
    TLegend *legRes = new TLegend(0.3, 0.3, 0.7, 0.6);
    for (int e=0; e<effSigma.size(); e++) {
      TH1D *hRatio = new TH1D(Form("hRatioEff%dRes%d", (int)(effSigma[e]*100.), (int)(timeRes[h]*10.)), Form("Efficiency ratio with %.2g ns smearing along S4 bar and efficiency sigma of %.2g m; Bar position / m; Smeared/unsmeared", timeRes[h], effSigma[e]), 20, 0, 1.4);
      hRatio->Divide(resPlots[e][h], effPlots[e]);
      hRatio->SetLineColor(52+e*4);
      hsRatioRes->Add(hRatio);
      legRes->AddEntry(hRatio, Form("#sigma_{Eff} = %.2g", effSigma[e]), "l");
    }
    hsRatioRes->Write();
    if (h==0) legRes->Write("legRes");
  }
    
  fout->Close();
  delete fout;
} // barEffSmearing
