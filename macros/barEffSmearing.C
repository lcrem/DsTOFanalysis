// barEffSmearing.C
// Simulates the effect of time resolution on the efficiency correction used
const std::vector<double> timeRes  = {1, 1.2, 1.4, 1.6, 1.8, 2.}; // in ns
const std::vector<double> effWidth = {0.2, 0.22, 0.24, 0.26}; // one sigma values for efficiency curves
const double timeToDist = 0.0768; // Factor to go from time to position in m

void setHistAttr(TH1D *h) 
{
  h->SetLineWidth(2);
  h->GetXaxis()->SetTitleSize(.05);
  h->GetYaxis()->SetTitleSize(.05);
  h->GetXaxis()->SetLabelSize(.05);
  h->GetYaxis()->SetLabelSize(.05);
  h->Sumw2();
}

void setHistAttr(TH2D *h2) 
{
  h2->GetXaxis()->SetTitleSize(.05);
  h2->GetYaxis()->SetTitleSize(.05);
  h2->GetZaxis()->SetTitleSize(.05);
  h2->GetXaxis()->SetLabelSize(.05);
  h2->GetYaxis()->SetLabelSize(.05);
  h2->GetZaxis()->SetLabelSize(.05);
  h2->Sumw2();
}

// Function for generating MC points
double gauss4(const double x, const double width, const double mean=.7)
{
  double val = TMath::Exp( -1*pow((x-mean)/(2*width), 4) );
  return val;
}

void barEffSmearing(const char* outFile, bool trueBars=false, const int nEntries=10000000, const char* s4Plots="/scratch0/sjones/plots/angularDistS4_newSample/withTree/includesUnweightedHists.root")
{
  gROOT->SetBatch(kTRUE);

  TFile *fout = new TFile(outFile, "recreate");
  TRandom3 *r = new TRandom3(666);

  if (!trueBars) {
    std::vector<TH1D*> effPlots;
    for (int e=0; e<effWidth.size(); e++) {
      TH1D *hEff = new TH1D(Form("hEff%d", (int)(effWidth[e]*100)), Form("Efficiency without smearing along S4 bar, %2g m sigma; Bar position / m; Events", effWidth[e]), 20, 0, 1.4);
      setHistAttr(hEff);
      effPlots.push_back(hEff);
    }
    std::vector< std::vector<TH1D*> > resPlots;
    for (int e=0; e<effPlots.size(); e++) {
      std::vector<TH1D*> resPlotsSmall;
      for (int h=0; h<timeRes.size(); h++) {
	TH1D *hTimeResEff = new TH1D(Form("hTimeResEff%dRes%d", (int)(effWidth[e]*100.), (int)(timeRes[h]*10.)), Form("Efficiency with %.2g ns smearing along S4 bar; Bar position / m; Events", timeRes[h]), 20, 0, 1.4);
	setHistAttr(hTimeResEff);
	resPlotsSmall.push_back(hTimeResEff);
      }
      resPlots.push_back(resPlotsSmall);
    }
    // Number of histogram entries to do
    for (int i=0; i<nEntries; i++) {
      if (i % 100000 == 0) std::cout<<i<<" of "<<nEntries<<std::endl;
      for (int e=0; e<effPlots.size(); e++) {
	double r1 = r->Rndm();
	double r2 = r->Rndm() * 1.4;
	if (r1 < gauss4(r2, effWidth.at(e))) {
	  effPlots[e]->Fill(r2); // Histogram with no smearing
	  // Loop over resolutions
	  for (int h=0; h<timeRes.size(); h++) {
	    resPlots[e][h]->Fill(r->Gaus(r2, timeRes[h]*timeToDist));
	  } // Loop over resolutions
	}
      } // Loop over efficiency sigmas
    } // Loop over the entries

    // Fit to smeared graphs
    for (int e=0; e<effPlots.size(); e++) {
      for (int h=0; h<timeRes.size(); h++) {
	fout->cd();
	TF1 *fGauss = new TF1(Form("fGaussEff%dRes%d", (int)(effWidth[e]*100.), (int)(timeRes[h]*10.)), "[0]*exp( -1*pow((x - [1])/(2*[2]), 4) )");
	fGauss->SetParameter(0, 1);
	fGauss->SetParameter(1, 0.7);
	fGauss->SetParameter(2, 0.25);
	resPlots[e][h]->Fit(fGauss, "q");
	fGauss->Write();
	cout<<"Efficiency sigma, time resolution "<<(int)(effWidth[e]*100.)<<"cm, "<<timeRes[h]<<"ns   ";
	cout<<"Measured mean, sigma "<<fGauss->GetParameter(1)*100.<<"cm, "<<fGauss->GetParameter(2)*100.<<"cm"<<endl;
      }
    }

    for (int e=0; e<effWidth.size(); e++) {
      effPlots[e]->Write();
      THStack *hsRatioEff = new THStack(Form("hsRatioEff%d", (int)(effWidth[e]*100)), Form("Smeared/unsmeared effs, efficiency sigma %.2g m; Bar position / m; Smeared/unsmeared", effWidth[e]));
      TLegend *legEff = new TLegend(0.3, 0.3, 0.7, 0.6);
      for (int h=0; h<timeRes.size(); h++) {
	resPlots[e][h]->SetLineColor(kBlack);
	resPlots[e][h]->SetLineWidth(2);
	resPlots[e][h]->Write();
	TH1D *hRatio = new TH1D(Form("hRatioEff%dRes%d", (int)(effWidth[e]*100.), (int)(timeRes[h]*10.)), Form("Efficiency ratio with %.2g ns smearing along S4 bar and efficiency sigma of %.2g m; Bar position / m; Smeared/unsmeared", timeRes[h], effWidth[e]), 20, 0, 1.4);
	setHistAttr(hRatio);
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
      for (int e=0; e<effWidth.size(); e++) {
	TH1D *hRatio = new TH1D(Form("hRatioEff%dRes%d", (int)(effWidth[e]*100.), (int)(timeRes[h]*10.)), Form("Efficiency ratio with %.2g ns smearing along S4 bar and efficiency sigma of %.2g m; Bar position / m; Smeared/unsmeared", timeRes[h], effWidth[e]), 20, 0, 1.4);
	setHistAttr(hRatio);
	hRatio->Divide(resPlots[e][h], effPlots[e]);
	hRatio->SetLineColor(52+e*4);
	hsRatioRes->Add(hRatio);
	legRes->AddEntry(hRatio, Form("#sigma_{Eff} = %.2g", effWidth[e]), "l");
      }
      hsRatioRes->Write();
      if (h==0) legRes->Write("legRes");
    }
  } // Generic bars

  else {
    TFile *fin = new TFile(s4Plots, "read");
    // Get the cosmic histograms
    for (int sample=0; sample < 4; sample++) {
      vector<double> meanVec;
      meanVec.resize(10, 0);
      cout<<"Sample "<<sample<<endl;
      TH2D *h2Cosmics = (TH2D*)fin->Get(Form("h2Cosmics%d", sample));
      for (int b=1; b<10; b++) {
	TH1D *hBarSmear = new TH1D(Form("hBlock%dBar%dSmear", sample, b), Form("Block %d, bar %d", sample, b), 20, 0., 140.);
	TH1D *hBarEff = new TH1D(Form("hBlock%dBar%dEff", sample, b), Form("Block %d, bar %d", sample, b), 20, 0., 140.);
	TH1D *hBarRatio = new TH1D(Form("hBlock%dBar%dRatio", sample, b), Form("Block %d, bar %d", sample, b), 20, 0., 140.);
	TF1 *f1 = new TF1( "f1", "[0]*exp(-1*pow((x-[1])/(2*[2]), 4))" );
	f1->SetParameter(0, 1.);
	f1->SetParameter(1, 70.);
	f1->SetParameter(2, 30.);
	TH1D *hCosmic = h2Cosmics->ProjectionX(Form("block%dbar%d", sample, b), b, b);
	hCosmic->Fit(f1);
	// We don't want to use the width, just the mean
	meanVec.at(b-1) = f1->GetParameter(1);
	// Loop over number of entries
	for (int n=0; n<nEntries; n++) {
	  if (n % 100000 == 0) std::cout<<n<<" of "<<nEntries<<std::endl;
	  double r1 = r->Rndm();
	  double r2 = r->Rndm() * 140;
	  if (r1 < gauss4(r2, 23., f1->GetParameter(1))) {
	    hBarEff->Fill(r2);
	    hBarSmear->Fill(r->Gaus(r2, 1.4*timeToDist*100.));
	  } 
	} // Loop over entries
	hBarRatio->Divide(hBarSmear, hBarEff, 1., 1., "B");
	fout->cd();
	hBarEff->Write();
	hBarSmear->Write();
	hBarRatio->Write();
      } // Loop over bars
    } // Loop over samples
    fin->Close();
    delete fin;
  }    

  fout->Close();
  delete fout;
} // barEffSmearing

