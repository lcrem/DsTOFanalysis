// fitGaussToCosmics.C
// Attempt to fit gaussians to 1D projections of the bar by bar cosmics

void fitGaussToCosmics(const char* outFile, const char* s4Plots="~/work/hptpc/angularDistS4_newSample/newCosmicsHists.root") 
{
  gROOT->SetBatch(kTRUE);
  TFile *fout = new TFile(outFile, "recreate");
  TFile *fin = new TFile(s4Plots, "read");
  // Loop over samples
  int nBad  = 0;
  int nStillBad  = 0;
  int nFits = 0;
  for (int sample=0; sample < 4; sample++) {
    cout<<"Sample "<<sample<<endl;
    TH2D *h2Cosmics     = (TH2D*)fin->Get(Form("h2CosmicsEff%d", sample));
    TH2D *h2CosmicsBins = (TH2D*)fin->Get(Form("h2CosmicsEffBins%d", sample));
    for (int b=1; b<10; b++) {
      cout<<"Block "<<sample<<", bar "<<b<<endl;
      fout->cd();
      TH1D *hCosmic     = h2Cosmics->ProjectionX(Form("block%dbar%d", sample, b), b, b);
      TH1D *hCosmicBins = h2CosmicsBins->ProjectionX(Form("block%dbar%dbins", sample, b), b, b);
      hCosmic->SetTitle(Form("%d blocks, bar %d; Efficiency; x / cm", sample, b));
      hCosmicBins->SetTitle(Form("%d blocks, bar %d (finer bins); Efficiency; x / cm", sample, b));
      hCosmic->SetLineWidth(2);
      hCosmicBins->SetLineWidth(2);
      cout<<"Two sigmoids fitted"<<endl;
      // Fit the sigmoids individually  to one side then use the results to tune
      TF1 *f2_1 = new TF1(Form("f2_1block%dbar%d", sample, b), "[0]/(1+exp(-[1]*(x-[2])))", 0, 70);
      f2_1->SetParameter(0, 1.);
      f2_1->SetParameter(1, 0.2);
      f2_1->SetParameter(2, 20.);
      hCosmic->Fit(f2_1, "r");
      TF1 *f2_2 = new TF1(Form("f2_2block%dbar%d", sample, b), "[0]/(1+exp(-[1]*(x-[2])))", 70, 140);
      f2_2->SetParameter(0, 1.);
      f2_2->SetParameter(1, -0.2);
      f2_2->SetParameter(2, 120.);
      hCosmic->Fit(f2_2, "r");
      
      TF1 *f2 = new TF1(Form("f2block%dbar%d", sample, b), "([0]/(1+exp(-[1]*(x-[2]))))*(1/(1+exp(-[3]*(x-[4]))))", 0, 140 );
      f2->SetParameter(0, f2_1->GetParameter(0));
      f2->SetParameter(1, f2_1->GetParameter(1));
      f2->SetParameter(2, f2_1->GetParameter(2));
      f2->SetParameter(3, f2_2->GetParameter(1));
      f2->SetParameter(4, f2_2->GetParameter(2));
      TFitResultPtr r = hCosmic->Fit(f2);
      f2->Write();
      TCanvas *c2 = new TCanvas(Form("c2_%d_%d", sample, b), Form("c2_%d_%d", sample, b), 1300, 500);
      c2->Divide(2, 1);
      c2->cd(1);
      hCosmic->Draw();
      TText *t1 = new TText(40, .15, Form("#Chi^{2} = %.4g/%d", f2->GetChisquare(), f2->GetNDF()));
      t1->Draw();
      //f2->Draw("same");
      hCosmic->Write();
      Int_t fitStatus = r;
      if (fitStatus != 0) {
	nBad++;
	// Try with slightly different parameters
	f2->SetParameter(0, 2);
	f2->SetParameter(1, 0.5);
	f2->SetParameter(2, 10);
	f2->SetParameter(3, -0.5);
	f2->SetParameter(4, 130);
	TFitResultPtr r2 = hCosmic->Fit(f2);
	Int_t fitStatus2 = r2;
	if (fitStatus2 != 0) {
	  nStillBad++;
	}
	else {
	  hCosmic->Write();
	}
      }
      nFits++;
      c2->cd(2);
      hCosmicBins->Fit(f2);
      hCosmicBins->Draw();
      TText *t2 = new TText(40, .15, Form("#Chi^{2} = %.4g/%d", f2->GetChisquare(), f2->GetNDF()));
      t2->Draw();
      hCosmicBins->Write();
      c2->Print(Form("~/work/hptpc/fitGaussToCosmics/compBlock%dbar%d.pdf", sample, b));
      c2->Print(Form("~/work/hptpc/fitGaussToCosmics/compBlock%dbar%d.tex", sample, b));
    } // Loop over bars
    delete h2Cosmics;
  } // Loop over samples
  cout<<nBad<<" bad fits out of "<<nFits<<endl;
  cout<<nStillBad<<" still bad fits (after changing start points) out of "<<nFits<<endl;
  
  fin->Close();
  delete fin;
  fout->Close();
  delete fout;
} // fitGaussToCosmics
