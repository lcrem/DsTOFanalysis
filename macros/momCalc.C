// momCalc.C
// The same analysis process as for bkgSub.C (efficiency calc and background subtraction)
// but now calculates the monetum of particles (of a given mass)

// Outputs momentum in GeV/c
double momFromTime(const double mass, const double baseline, const double time)
{
  double mom = 0.;
  mom = (mass / 3e8) * baseline * ( 1. / TMath::Sqrt( pow(time*1e-9, 2) - (pow((baseline / 9e16), 2)) ) );
  return mom;
}

void momCalc(const char* saveDir, const double mass,
	     const char* dstofDir="/scratch0/dbrailsf/temp/mylinktodtof/") 
{
  gSystem->Load("libdstof.so");
  gROOT->SetBatch(kTRUE);
  // Unix timestamps for variable block moves
  // 0.8GeV/c, 0 blocks
  const double start0Block = 1535713289; 
  const double end0Block   = 1535716132;
  // 0.8GeV/c, 1 block
  const double start1Block = 1535796057;
  const double end1Block   = 1535799112;
  // 0.8GeV/c, 2 blocks
  const double start2Block = 1535789157;
  const double end2Block   = 1535792026;
  // 0.8GeV/c, 3 block
  const double start3Block = 1535792404;
  const double end3Block   = 1535795300;
  // 0.8GeV/c, 4 block
  // Most runs were in this configuration so don't need to use necessarily
  const double start4Block = 1536537600; 
  const double end4Block   = 1536669600;
  // Timing cuts
  const double piLow  = 80.;
  const double piHi   = 95.;
  const double proLow = 106.;
  const double proHi  = 134.;
  // S1 -> S4 baseline length
  const double baselineS1S4 = 13.97;
  const double baselineS1S3 = 10.9;

  // ustof-dstof cable delay
  const double ustofDelay = 184.7;
  // Shift in ns required to to pion peak at speed of light
  const double dstofShift = 40.;

  TFile *fout = new TFile(Form("%s/bkgSubPlots_eff.root", saveDir), "recreate");

  for (int nBlocks = 0; nBlocks <= 3; nBlocks++) {
    double nP  = 0.;
    double nPi = 0.;
    // Define signal and background functions to be fitted
    // Signals are gaussians
    TF1 *sPro = new TF1("sPro", "gaus", 66, 100);
    TF1 *sPi  = new TF1("sPi", "gaus", 35, 55);
    // Exponential background
    TF1 *fBkgExp = new TF1("fBkgExp","expo", 30, 160);
    sPro->SetLineColor(kGreen+2);
    sPi->SetLineColor(kRed);
    TF1 *fSplusBExp = new TF1("signal+bkg exp", "gaus(0)+gaus(3)+expo(6)", 30, 160);
    fSplusBExp->SetParNames("const 1", "mean 1", "sigma 1",
			    "const 2", "mean 2", "sigma 2",
			    "bkgconst", "bkgdecay");
    fSplusBExp->SetLineColor(kBlack);

    int nSpills = 0;
    int nSpillsTrue = 0;
    double lastSpill = 0.;

    TH1D *hdtof1d = new TH1D(Form("hdtof1d_%d",nBlocks), Form("Time of flight, %d blocks; S4 - S1 / ns; Events / spill", nBlocks), 260, 30, 160);
    TH1D *hpromom = new TH1D(Form("hpromom_%d",nBlocks), Form("Proton momentum in S4, %d blocks; Proton momentum [GeV/c]; Events / spill", nBlocks), 300, 0., 1);

    // Find the correct dstof files
    Int_t runMin=-1;
    Int_t runMax=-1;

    double startTime = 0;
    double endTime   = 0;
    if (nBlocks == 0) {
      startTime = start0Block;
      endTime   = end0Block;
    }
    else if (nBlocks == 1) {
      startTime = start1Block;
      endTime   = end1Block;
    }
    else if (nBlocks == 2) {
      startTime = start2Block;
      endTime   = end2Block;
    }
    else if (nBlocks == 3) {
      startTime = start3Block;
      endTime   = end3Block;
    }
    else if (nBlocks == 4) {
      startTime = start4Block;
      endTime   = end4Block;
    }

    for (int irun=950; irun<1400; irun++){
      TFile *fin = new TFile(Form("%srun%d/DsTOFcoincidenceRun%d_tdc1.root", dstofDir, irun, irun), "read");
      RawDsTofCoincidence *tofCoinTemp = NULL;
      TTree *tree = (TTree*) fin->Get("tofCoinTree");
      tree->SetBranchAddress("tofCoin", &tofCoinTemp);
      tree->GetEntry(0);
      UInt_t firstTemp = tofCoinTemp->unixTime[0];
      tree->GetEntry(tree->GetEntries()-1);
      UInt_t lastTemp = tofCoinTemp->unixTime[0];
      
      fin->Close();
      delete fin;
      
      if (firstTemp>endTime){
	break;
      }
      
      if (firstTemp<startTime && lastTemp>startTime){
	runMin = irun;
      }
      
      if (firstTemp<endTime && lastTemp>endTime){
	runMax = irun;
      }   
    }
    
    cout << "Min and max runs are " << runMin << " " << runMax << endl;

    TH1D *hCoins = new TH1D(Form("hCoins_%d",nBlocks), Form("Bar coincidences + S_{1,2} coincidences, %d blocks; Bar; Events",nBlocks), 10, 0.5, 10.5);
    TH1D *hHits  = new TH1D(Form("hHits_%d",nBlocks), Form("PMT hits + S_{1,2} coincidences, %d blocks; Bar; Events",nBlocks), 10, 0.5, 10.5);
    TH1D *hEff   = new TH1D(Form("hEff_%d",nBlocks), Form("S4 bar efficiencies %d blocks; Bar; Events",nBlocks), 10, 0.5, 10.5);
    hHits->Sumw2();
    hCoins->Sumw2();

    // In this loop calculate the bar-by-bar efficiencies
    for (int itdc=0; itdc<2; itdc++) {
      // Need both the coincidence and the raw trees for this
      double tempUstof;
      double ustofNs;
      for (int irun=runMin; irun<runMax+1; irun++){
	// Load input files
	TFile *tofCoinFile = new TFile(Form("%srun%d/DsTOFcoincidenceRun%d_tdc%d.root", dstofDir, irun, irun, itdc+1));
	TFile *tofFile     = new TFile(Form("%srun%d/DsTOFtreeRun%d_tdc%d.root", dstofDir, irun, irun, itdc+1));
	TTree *tofCoinTree = (TTree*)tofCoinFile->Get("tofCoinTree");
	TTree *tofTree = (TTree*)tofFile->Get("tofTree");
	RawDsTofCoincidence *tofCoin = NULL;
	RawDsTofHeader *tof = NULL;
	tofCoinTree->SetBranchAddress("tofCoin", &tofCoin);
	tofTree->SetBranchAddress("tof", &tof);
	// Create new TTree for ustof signals
	TTree *ustofTree = new TTree("ustofTree", "ustof");;
	ustofTree->SetDirectory(0);
	ustofTree->Branch("ustofNs", &ustofNs, "ustofNs/D");
	tempUstof=0;
	// Build tree of all the ustof hit times
	for (int i=0; i<tofTree->GetEntries(); i++) {
	  tofTree->GetEntry(i);
	  if (tof->unixTime < startTime) continue;
	  if (tof->unixTime > endTime) break;

	  if (tof->channel == 13) {
	    ustofNs = tof->fakeTimeNs - ustofDelay;
	    // ustof signals shouldn't be coming closer than 500ns
	    if ( (ustofNs - tempUstof) < 500.) continue;
	    ustofTree->Fill();
	    tempUstof = ustofNs;
	  } // if (tof->channel == 13) 
	} // for (int i=0; i<tof->GetEntries(); i++)
	ustofTree->BuildIndex("ustofNs");
	// Now loop over coincidence file to find number of coincidences for each bar
	for (int t=0; t<tofCoinTree->GetEntries(); t++) {
	  tofCoinTree->GetEntry(t);
	  if (tofCoin->unixTime[0] < startTime) continue;
	  if (tofCoin->unixTime[0] > endTime) break;
	  double deltat = TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1]);
	  double dstofHitT = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (10. - TMath::Abs(deltat) / 2. );
	  double tofCalc = dstofHitT - tofCoin->usTofSignal;
	  if (tofCalc > 70. && tofCalc < 200. && tofCoin->bar != 10) {
	    hCoins->Fill(tofCoin->bar);
	  } // if (tofCalc > 70. && tofCalc < 200.)
	} // for (int t=0; t<tofTree->GetEntries(); t++) 
	// Loop over hit file to find number of hits in coincidence 
	for (int t=0; t<tofTree->GetEntries(); t++) {
	  tofTree->GetEntry(t);
	  if (tof->unixTime < startTime) continue;
	  if (tof->unixTime > endTime) break;
	  if (tof->channel > 10 || tof->channel == 0) continue;
	  int ustofEntry = ustofTree->GetEntryNumberWithBestIndex(tof->fakeTimeNs);
	  ustofTree->GetEntry(ustofEntry);
	  if (tof->fakeTimeNs - ustofNs < 260. && tof->fakeTimeNs - ustofNs > 130.) {
	    if (itdc == 0) { // TDC1
	      if (tof->channel % 2 == 1) { // PMT B
		hHits->Fill(((tof->channel + 1)/ -2) + 11);
	      } // if (tof->channel % 2 == 1)
	      else { // PMT A
		hHits->Fill((tof->channel / -2) + 11);
	      }
	    } // TDC 1
	    else { // TDC 2
	      if (tof->channel % 2 == 1) { // PMT B
		hHits->Fill(((tof->channel + 1)/ -2) + 6);
	      } // if (tof->channel % 2 == 1)
	      else { // PMT A
		hHits->Fill((tof->channel / -2) + 6);
	      }
	    } // TDC 2
	  } // if (tof->fakeTimeNs - ustofNs < 260. && tof->fakeTimeNs - ustofNs > 130.) 
	} // for (int t=0; t<tofTree->GetEntries(); t++)
	delete tof;
	delete tofCoin;
	delete tofCoinTree;
	delete ustofTree;
      } // for (int irun=runMin; irun<runMax+1; irun++)
    } // for (int itdc=0; itdc<2; itdc++) 

    fout->cd(0);
    TCanvas *cHits = new TCanvas(Form("%d_cHits",nBlocks));
    hHits->Draw("hist");
    hHits->Write();
    cHits->Print(Form("%s/%d_barHits.png",saveDir,nBlocks));
    cHits->Print(Form("%s/%d_barHits.pdf",saveDir,nBlocks));
    TCanvas *cCoins = new TCanvas(Form("%d_cCoins",nBlocks));
    hCoins->Draw("hist");
    hCoins->Write();
    cCoins->Print(Form("%s/%d_barCoins.png",saveDir,nBlocks));
    cCoins->Print(Form("%s/%d_barCoins.pdf",saveDir,nBlocks));
    TCanvas *cEff = new TCanvas(Form("%d_cEff",nBlocks));
    hHits->Add(hCoins, -1.);
    hEff->Divide(hCoins, hHits, 1., 1., "B");
    hEff->Draw("hist");
    hEff->Write();
    cEff->Print(Form("%s/%d_barEff.png",saveDir,nBlocks));
    cEff->Print(Form("%s/%d_barEff.pdf",saveDir,nBlocks));

    for (int itdc=0; itdc<2; itdc++) {
      TChain *tofCoinChain = new TChain("tofCoinTree");
      for (int irun=runMin; irun<runMax+1; irun++){
	tofCoinChain->Add(Form("%srun%d/DsTOFcoincidenceRun%d_tdc%d.root", dstofDir, irun, irun, itdc+1));
      } // for (int irun=runMin; irun<runMax+1; irun++)
      RawDsTofCoincidence *tofCoin = NULL;
      tofCoinChain->SetBranchAddress("tofCoin", &tofCoin);
      for (int h=0; h<tofCoinChain->GetEntries(); h++) {
	tofCoinChain->GetEntry(h);
	if (tofCoin->unixTime[0]<startTime) continue;
	if (tofCoin->unixTime[0]>endTime) break;

	if (tofCoin->lastDelayedBeamSignal != lastSpill && itdc == 0) {
	  lastSpill = tofCoin->lastDelayedBeamSignal;
	  nSpills++;
	  for (int sp=h; sp<tofCoinChain->GetEntries(); sp++) {
	    tofCoinChain->GetEntry(sp);
	    if (tofCoin->unixTime[0]<startTime) continue;
	    if (tofCoin->unixTime[0]>endTime) break;

	    if ((tofCoin->fakeTimeNs[0] - tofCoin->usTofSignal) < 200. &&
		(tofCoin->fakeTimeNs[0] - tofCoin->usTofSignal) > 70. &&
		(tofCoin->fakeTimeNs[0] - lastSpill) < 1e9 &&
		(tofCoin->fakeTimeNs[0] - lastSpill) > 0) {
	      nSpillsTrue++;
	      break;
	    }

	  } // for (int sp=h; sp<tofCoinChain->GetEntries(); sp++) 
	} // if (tofCoin->lastDelayedBeamSignal != lastSpill && itdc == 0)
	// Need to calculate total signal hits here
	// Weight these by the efficiency of calculated above
	double deltat = TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1]  );
	double dstofHitT = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (10. - TMath::Abs(deltat) / 2. );
	double tofCalc = dstofHitT - tofCoin->usTofSignal - dstofShift;
	if (tofCalc < 160. && tofCalc > 30. && tofCoin->bar != 10) {
	  hdtof1d->Fill(tofCalc, 1. / hEff->GetBinContent(tofCoin->bar));
	  if (tofCalc > proLow - dstofShift && tofCalc < proHi - dstofShift) {
	    hpromom->Fill(momFromTime(0.938, baselineS1S4, tofCalc), 1. / hEff->GetBinContent(tofCoin->bar));
	  }
	} // if (tofCalc < 200. && tofCalc > 70.) 
      } // for (int h=0; h<tofCoinChain->GetEntries(); h++) 
      delete tofCoin;
      delete tofCoinChain;
    } // for (int itdc=0; itdc<2; itdc++)

    fout->cd(0);
    TCanvas *c2d_exp = new TCanvas(Form("%d_c2d_exp",nBlocks));
    c2d_exp->SetLogy();
    hdtof1d->Fit(sPi, "R");
    hdtof1d->Fit(sPro, "R");
    hdtof1d->Fit(fBkgExp, "R");
    Double_t parExp[8];
    sPro->GetParameters(&parExp[0]);
    sPi->GetParameters(&parExp[3]);
    fBkgExp->GetParameters(&parExp[6]);
    fSplusBExp->SetParameters(parExp);
    hdtof1d->Fit(fSplusBExp, "R");
    hdtof1d->Draw("hist");
    fSplusBExp->Draw("same");
    hdtof1d->Write();

    TCanvas *cmom = new TCanvas(Form("%d_cmom",nBlocks));
    cmom->SetLogy();
    hpromom->Draw("hist");
    hpromom->Write();
    cmom->Print(Form("%s/%d_protonMomentum.png", saveDir, nBlocks));
    cmom->Print(Form("%s/%d_protonMomentum.pdf", saveDir, nBlocks));

    // Now we have the fit values, loop over again and subtract the background
    // For each bin, find the fraction of each particle type which are background and 
    // this fraction from the bin
    /*    
    TF1 *fSub = new TF1("fSub", "exp([0]+[1]*x)", 70, 200);
    fSub->SetParameter(0, fSplusBExp->GetParameter("bkgconst"));  
    fSub->SetParameter(1, fSplusBExp->GetParameter("bkgdecay"));
    // Background subtracted tof spectrum
    TH1D *hdtof1d_sub = (TH1D*)hdtof1d->Clone(Form("%d_hdtof1d_sub",nBlocks));
    hdtof1d_sub->Add(fSub, -1.);
    // If the bin content drops below 0, set to 0
    for (int b = 0; b <=  hdtof1d_sub->GetNbinsX(); b++) {
      if (hdtof1d_sub->GetBinContent(b) < 0.) {
	hdtof1d_sub->SetBinContent(b, 0);
      } // if (hdtof1d_sub->GetBinContent(b) < 0.)
    } // for (int b = 0; hdtof1d_sub->GetNbinsX(); b++)
    TCanvas *cdtof_sub = new TCanvas(Form("%d_cdtof_sub",nBlocks));
    cdtof_sub->SetLogy();
    hdtof1d_sub->Draw("hist");
    hdtof1d_sub->Write();
    */

  } // for (int nBlocks=firstBlock; nBlocks <= lastBlock; nBlocks++)

} // momCalc
