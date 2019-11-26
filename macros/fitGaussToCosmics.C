// fitGaussToCosmics.C
// Attempt to fit gaussians to 1D projections of the bar by bar cosmics

void fitGaussToCosmics(const char* outFile, const char* s4Plots="~/work/hptpc/angularDistS4_newSample/newCosmicsHists.root") 
{
  gROOT->SetBatch(kTRUE);
  TFile *fout = new TFile(outFile, "recreate");
  TFile *fin = new TFile(s4Plots, "read");
  TH1D *hGausSigma = new TH1D("hGausSigma", "#sigma for bar efficiencies; #sigma / cm; Events", 10, 26, 36);
  TH1D *hGausMean = new TH1D("hGausMean", "#mu for bar efficiencies; #mu / cm; Events", 20, 50., 90.);
  // Loop over samples
  int nBad  = 0;
  int nStillBad  = 0;
  int nFits = 0;
  for (int sample=0; sample < 4; sample++) {
    cout<<"Sample "<<sample<<endl;
    TH2D *h2Cosmics = (TH2D*)fin->Get(Form("h2CosmicsEffBins%d", sample));
    for (int b=1; b<10; b++) {
      cout<<"Block "<<sample<<", bar "<<b<<endl;
      cout<<"Gaussian to 4th power"<<endl;
      TF1 *f1 = new TF1( Form("f1block%dbar%d", sample, b), "[0]*exp(-1*pow((x-[1])/(2*[2]), 4))" );
      f1->SetParameter(0, 1);
      f1->SetParameter(1, 70);
      f1->SetParameter(2, 30);
      TH1D *hCosmic = h2Cosmics->ProjectionX(Form("block%dbar%d", sample, b), b, b);
      hCosmic->Fit(f1,"q");
      hGausSigma->Fill(f1->GetParameter(2));
      hGausMean->Fill(f1->GetParameter(1));
      fout->cd();
      f1->Write();

      cout<<"One sigmoid fitted"<<endl;
      TF1 *ftmp = new TF1(Form("ftmpblock%dbar%d", sample, b), "([0]/(1+exp(-[1]*(x-[2]))))", 0, 40);
      ftmp->SetParameter(0, 2);
      ftmp->SetParameter(1, 0.15);
      ftmp->SetParameter(2, 20);
      hCosmic->Fit(ftmp, "q", "R");
      ftmp->Write();
      
      cout<<"Gaussian to 6th power"<<endl;
      TF1 *f3 = new TF1( Form("f3block%dbar%d", sample, b), "[0]*exp(-1*pow((x-[1])/(2*[2]), 6))" );
      f3->SetParameter(0, 1);
      f3->SetParameter(1, 70);
      f3->SetParameter(2, 25);
      hCosmic->Fit(f3, "q");
      f3->Write();

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
    }
    delete h2Cosmics;
  } // Loop over samples
  cout<<nBad<<" bad fits out of "<<nFits<<endl;
  cout<<nStillBad<<" still bad fits (after changing start points) out of "<<nFits<<endl;
  fout->cd();
  hGausSigma->Write();
  hGausMean->Write();
  
  fin->Close();
  delete fin;
  fout->Close();
  delete fout;
} // fitGaussToCosmics
