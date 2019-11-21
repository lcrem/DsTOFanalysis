// fitGaussToCosmics.C
// Attempt to fit gaussians to 1D projections of the bar by bar cosmics

void fitGaussToCosmics(const char* outFile, const char* s4Plots="~/work/hptpc/angularDistS4_newSample/includesUnweightedHists.root") 
{
  gROOT->SetBatch(kTRUE);
  TFile *fout = new TFile(outFile, "recreate");
  TFile *fin = new TFile(s4Plots, "read");
  TH1D *hGausSigma = new TH1D("hGausSigma", "#sigma for bar efficiencies; #sigma / cm; Events", 10, 26, 36);
  TH1D *hGausMean = new TH1D("hGausMean", "#mu for bar efficiencies; #mu / cm; Events", 20, 50., 90.);
  // Loop over samples
  for (int sample=0; sample < 4; sample++) {
    cout<<"Sample "<<sample<<endl;
    TH2D *h2Cosmics = (TH2D*)fin->Get(Form("h2Cosmics%d", sample));
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
      cout<<"Gauss4 - Gauss6 = "<<f1->GetChisquare() - f3->GetChisquare()<<endl;

      cout<<"Two sigmoids fitted"<<endl;
      TF1 *f2 = new TF1(Form("f2block%dbar%d", sample, b), "([0]/(1+exp(-[1]*(x-[2]))))*(1/(1+exp(-[3]*(x-[4]))))", 0, 140 );
      f2->SetParameter(0, 2);
      f2->SetParameter(1, 0.15);
      f2->SetParameter(2, 20);
      f2->SetParameter(3, -0.15);
      f2->SetParameter(4, 120);
      f2->SetParameter(5, 0.2);
      f2->SetParameter(6, 0.2);
      hCosmic->Fit(f2);
      f2->Write();

      hCosmic->Write();
    }
    delete h2Cosmics;
  } // Loop over samples
  fout->cd();
  hGausSigma->Write();
  hGausMean->Write();
  
  fin->Close();
  delete fin;
  fout->Close();
  delete fout;
} // fitGaussToCosmics
