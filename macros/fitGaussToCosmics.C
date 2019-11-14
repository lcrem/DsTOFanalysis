// fitGaussToCosmics.C
// Attempt to fit gaussians to 1D projections of the bar by bar cosmics

void fitGaussToCosmics(const char* outFile, const char* s4Plots="~/work/hptpc/angularDistS4_newSample/includesUnweightedHists.root") 
{
  gROOT->SetBatch(kTRUE);
  TFile *fout = new TFile(outFile, "recreate");
  TFile *fin = new TFile(s4Plots, "read");
  TH1D *hGausSigma = new TH1D("hGausSigma", "#sigma for bar efficiencies; #sigma / cm; Events", 10, 26, 36);
  TH1D *hGausMean = new TH1D("hGausMean", "#mu for bar efficiencies; #mu / cm; Events", 20, 50., 90.);

  TH1D *hHalfGausSigma = new TH1D("hHalfGausSigma", "#sigma for bar efficiencies with half gaussians fitted; #sigma / cm; Events", 20, 4, 24);
  TH1D *hHalfGausMean1 = new TH1D("hHalfGausMean1", "#mu_{1} for bar efficiencies with half gaussians fitted; #mu / cm; Events", 20, 20., 50.);
  TH1D *hHalfGausMean2 = new TH1D("hHalfGausMean2", "#mu_{2} for bar efficiencies with half gaussians fitted; #mu / cm; Events", 20, 90., 120.); 
  // Loop over samples
  for (int sample=0; sample < 4; sample++) {
    cout<<"Sample "<<sample<<endl;
    TH2D *h2Cosmics = (TH2D*)fin->Get(Form("h2Cosmics%d", sample));
    for (int b=1; b<10; b++) {
      TF1 *f1 = new TF1( "f1", "[0]*exp(-1*pow((x-[1])/(2*[2]), 4))" );
      f1->SetParameter(0, 1);
      f1->SetParameter(1, 70);
      f1->SetParameter(2, 30);
      TH1D *hCosmic = h2Cosmics->ProjectionX(Form("block%dbar%d", sample, b), b, b);
      hCosmic->Fit(f1);
      hGausSigma->Fill(f1->GetParameter(2));
      hGausMean->Fill(f1->GetParameter(1));
      // Same as this but fit some half Gaussians to the edges of the distributions
      // TF1 *fHalf1 = new TF1("fHalf1", "gaus", 0, 40);
      // TF1 *fHalf2 = new TF1("fHalf2", "gaus", 100, 140);
      // hCosmic->Fit(fHalf1, "R");
      // hCosmic->Fit(fHalf2, "R");
      // hHalfGausSigma->Fill(fHalf1->GetParameter(2));
      // hHalfGausMean1->Fill(fHalf1->GetParameter(1));
      // hHalfGausSigma->Fill(fHalf2->GetParameter(2));
      // hHalfGausMean2->Fill(fHalf2->GetParameter(1));
      fout->cd();
      hCosmic->Write();
      f1->Write(Form("fitBlock%dBar%d", sample, b));
    }
    delete h2Cosmics;
  } // Loop over samples
  fout->cd();
  hGausSigma->Write();
  hGausMean->Write();
  hHalfGausSigma->Write();
  hHalfGausMean1->Write();
  hHalfGausMean2->Write();
  
  fin->Close();
  delete fin;
  fout->Close();
  delete fout;
} // fitGaussToCosmics
