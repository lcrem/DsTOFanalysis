// dataS4Cosmics.C

void dataS4Cosmics(const char* saveDir, const int startRun, const int endRun, 
		   const char* dstofDir="/nfs/scratch0/dbrailsf/data_backup/dtof_backup/")
{
  gROOT->SetBatch(kTRUE);
  TFile *fout = new TFile(saveDir, "recreate");

  TH1D *hCosmicRate = new TH1D("hCosmicRate", "Cosmic rate in S4; Bar; Rate / s^{-1}", 10, 0.5, 10.5);
  hCosmicRate->GetXaxis()->SetTitleSize(.05);
  hCosmicRate->GetXaxis()->SetLabelSize(.05);
  hCosmicRate->GetYaxis()->SetTitleSize(.05);
  hCosmicRate->GetYaxis()->SetLabelSize(.05);
  hCosmicRate->GetZaxis()->SetTitleSize(.05);
  hCosmicRate->GetZaxis()->SetLabelSize(.05);
  hCosmicRate->Sumw2();
  TH2D *h2CosmicRate = new TH2D("h2CosmicRate", "Cosmic rate in S4; Bar position / cm; Bar; Rate / s^{-1}", 20, 0., 140., 10, 0.5, 10.5);
  h2CosmicRate->GetXaxis()->SetTitleSize(.05);
  h2CosmicRate->GetXaxis()->SetLabelSize(.05);
  h2CosmicRate->GetYaxis()->SetTitleSize(.05);
  h2CosmicRate->GetYaxis()->SetLabelSize(.05);
  h2CosmicRate->GetZaxis()->SetTitleSize(.05);
  h2CosmicRate->GetZaxis()->SetLabelSize(.05);
  h2CosmicRate->Sumw2();


  int nSteps = 10;
  for (int n = 0; n < nSteps; n++) {
    int nSpills = 0;
    double deadtimeCut = (700./nSteps) * n;
    vector<double> deadtimeVec;
    deadtimeVec.resize(10, 0.);
    TH2D *hCutCosmics = new TH2D(Form("hCutCosmics_%d", (int)deadtimeCut), "Cosmics; xpos / cm; Bar", 20, 0., 140, 10., 0.5, 10.5);
    for (int irun = startRun; irun <= endRun; irun++) {
      cout<<"Run "<<irun<<endl;
      double lastSpill = 0.;
      for (int itdc = 0; itdc < 2; itdc++) {
	TFile *fin = new TFile(Form("%srun%d/DsTOFcoincidenceRun%d_tdc%d.root", dstofDir, irun, irun, itdc+1), "read");
	RawDsTofCoincidence *tofCoin = NULL;
	TTree *tree = (TTree*)fin->Get("tofCoinTree");
	tree->SetBranchAddress("tofCoin", &tofCoin);
	for (int e=0; e<tree->GetEntries(); e++) {
	  tree->GetEntry(e);

	  if (tofCoin->lastDelayedBeamSignal != lastSpill && itdc == 0) {
	    lastSpill = tofCoin->lastDelayedBeamSignal;
	    nSpills++;
	  }
	  if (tofCoin->inSpill) continue;
	  double dstofHitT = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (10.-TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1])/2.);
	  if (!tofCoin->inSpill && dstofHitT - deadtimeVec.at(tofCoin->bar-1) > deadtimeCut) {
	    hCosmicRate->Fill(tofCoin->bar);
	    deadtimevec.at(tofCoin->bar-1) = dstofHitT;
	    double xpos = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 70.));
	    hCutCosmics->Fill(xpos, tofCoin->bar);
	    h2CosmicRate->Fill(xpos, tofCoin->bar);
	  }
	} // Loop over entries
	delete tofCoin;
	fin->Close();
	delete fin;
      } // Loop over TDCs
    } // Loop over runs
    fout->cd();
    double totalTime = (endRun - startRun + 1.) * 3600. - nSpills;
    hCutCosmics->Scale(1./totalTime);
    hCutCosmics->Write();
  }
  //  hCosmicRate->Scale(1./totalTime);
  //  h2CosmicRate->Scale(1./totalTime);

  fout->cd();
  hCosmicRate->Write();
  h2CosmicRate->Write();

  fout->Close();
  delete fout;
} // dataS4Cosmics
