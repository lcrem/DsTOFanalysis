// bkgSub.C
// Performs background subtraction on the dtof spectrum

void bkgSub (const char* saveDir, 
	     const int firstBlock,
	     const int lastBlock,
	     const char* dstofDir="/scratch0/dbrailsf/temp/mylinktodtof/") 
{
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
  // ustof-dstof cable delay
  const double ustofDelay = 184.7;

  TFile *fout = new TFile(Form("%s/bkgSubPlots_eff.root", saveDir), "recreate");

  for (int nBlocks=firstBlock; nBlocks <= lastBlock; nBlocks++) {
    double nP  = 0.;
    double nPi = 0.;
    double nP_noeff  = 0.;
    double nPi_noeff = 0.;
    // Define signal and background functions to be fitted
    // Signals are gaussians
    TF1 *sPro = new TF1("sPro", "gaus", 106, 140);
    TF1 *sPi  = new TF1("sPi", "gaus", 75, 95);
    TF1 *sPro_noeff = new TF1("sPro_noeff", "gaus", 106, 140);
    TF1 *sPi_noeff  = new TF1("sPi_noeff", "gaus", 75, 95);
    // Background is flat
    TF1 *fBkg    = new TF1("fBkg", "[0]", 70, 180);
    TF1 *fBkgLin = new TF1("fBkgLin","[0] + [1]*x", 70, 200);
    TF1 *fBkgExp = new TF1("fBkgExp","expo", 70, 200);
    TF1 *fBkgExp_noeff = new TF1("fBkgExp_noeff","expo", 70, 200);
    fBkg->SetLineColor(kBlue);
    fBkgLin->SetLineColor(kBlue);
    sPro->SetLineColor(kGreen+2);
    sPi->SetLineColor(kRed);

    TF1 *fSplusB    = new TF1("signal+bkg", "gaus(0)+gaus(3)+[6]", 70, 180);
    TF1 *fSplusBLin = new TF1("signal+bkg lin", "gaus(0)+gaus(3)+[6]+[7]*x", 70, 200);
    TF1 *fSplusBExp = new TF1("signal+bkg exp", "gaus(0)+gaus(3)+expo(6)", 70, 200);
    TF1 *fSplusBExp_noeff = new TF1("signal+bkg exp no eff", "gaus(0)+gaus(3)+expo(6)", 70, 200);
    fSplusB->SetParNames("const 1", "mean 1", "sigma 1",
			 "const 2", "mean 2", "sigma 2",
			 "bkg");
    fSplusBLin->SetParNames("const 1", "mean 1", "sigma 1",
			    "const 2", "mean 2", "sigma 2",
			    "bkgconst", "bkgslope");
    fSplusBExp->SetParNames("const 1", "mean 1", "sigma 1",
			    "const 2", "mean 2", "sigma 2",
			    "bkgconst", "bkgdecay");
    fSplusBExp_noeff->SetParNames("const 1", "mean 1", "sigma 1",
				  "const 2", "mean 2", "sigma 2",
				  "bkgconst", "bkgdecay");
    fSplusB->SetLineColor(kBlack);
    fSplusBLin->SetLineColor(kBlack);
    fSplusBExp->SetLineColor(kBlack);
    fSplusBExp_noeff->SetLineColor(kBlack);

    int nSpills = 0;
    int nSpillsTrue = 0;
    double lastSpill = 0.;
    
    TH1D *hproHitsDstofVert = new TH1D(Form("hpiHitsDstofVert_%d",nBlocks), Form("Number of MIPs in S4, %d blocks; Bar; Protons",nBlocks), 10, 0.5, 10.5);
    TH1D *hpiHitsDstofVert = new TH1D(Form("hproHitsDstofVert_%d",nBlocks), Form("Number of MIPs, %d blocks; Bar; MIPs",nBlocks), 10, 0.5, 10.5);
    hproHitsDstofVert->Sumw2();
    hpiHitsDstofVert->Sumw2();
    TH1D *hproHitsDstofHorz = new TH1D(Form("hpiHitsDstofHorz_%d",nBlocks), Form("Number of MIPs in S4, %d blocks; x / cm; Protons",nBlocks), 20, 0., 140.);
    TH1D *hpiHitsDstofHorz = new TH1D(Form("hproHitsDstofHorz_%d",nBlocks), Form("Number of MIPs, %d blocks; x / cm; MIPs",nBlocks), 20, 0., 140.);
    hproHitsDstofHorz->Sumw2();
    hpiHitsDstofHorz->Sumw2();
    // No efficiency
    TH1D *hproHitsDstofVert_noeff = new TH1D(Form("hpiHitsDstofVert_noeff_%d",nBlocks), Form("Number of MIPs in S4, %d blocks (no efficiency); Bar; Protons",nBlocks), 10, 0.5, 10.5);
    TH1D *hpiHitsDstofVert_noeff = new TH1D(Form("hproHitsDstofVert_noeff_%d",nBlocks), Form("Number of MIPs, %d blocks (no efficiency); Bar; MIPs",nBlocks), 10, 0.5, 10.5);
    hproHitsDstofVert_noeff->Sumw2();
    hpiHitsDstofVert_noeff->Sumw2();
    TH1D *hproHitsDstofHorz_noeff = new TH1D(Form("hpiHitsDstofHorz_noeff_%d",nBlocks), Form("Number of MIPs in S4, %d blocks (no efficiency); x / cm; Protons",nBlocks), 20, 0., 140.);
    TH1D *hpiHitsDstofHorz_noeff = new TH1D(Form("hproHitsDstofHorz_noeff_%d",nBlocks), Form("Number of MIPs, %d blocks (no efficiency); x / cm; MIPs",nBlocks), 20, 0., 140.);
    hproHitsDstofHorz_noeff->Sumw2();
    hpiHitsDstofHorz_noeff->Sumw2();

    TH1D *hdtof1d = new TH1D(Form("hdtof1d_%d",nBlocks), Form("Time of flight, %d blocks; S4 - S1 / ns; Events / spill", nBlocks), 260, 70, 200);
    TH1D *hdtof1d_noeff = new TH1D(Form("hdtof1d_noeff_%d",nBlocks), Form("Time of flight, %d blocks (no efficiency); S4 - S1 / ns; Events / spill", nBlocks), 260, 70, 200);


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
    TH1D *htof1dHitsA = new TH1D(Form("htof1dHitsA_%d",nBlocks), Form("Time of flight, %d blocks; PMT_{A} - UsToF / ns; Events", nBlocks), 260, 130, 260);
    // In this loop calculate the bar-by-bar efficiencies
    for (int itdc=0; itdc<2; itdc++) {
      // Need both the coincidence and the raw trees for this
      // TChain *tofCoinChain = new TChain("tofCoinTree");
      // TChain *tofChain     = new TChain("tofTree");
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
	  double deltat = TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1]  );
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
	  htof1dHitsA->Fill(tof->fakeTimeNs - ustofNs);
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
      } // for (int irun=runMin; irun<runMax+1; irun++)
    } // for (int itdc=0; itdc<2; itdc++) 

    fout->cd(0);
    TCanvas *cHits = new TCanvas("cHits");
    hHits->Draw("hist");
    hHits->Write();
    cHits->Print(Form("%s/%d_barHits.png",saveDir,nBlocks));
    cHits->Print(Form("%s/%d_barHits.pdf",saveDir,nBlocks));
    TCanvas *cCoins = new TCanvas("cCoins");
    hCoins->Draw("hist");
    hCoins->Write();
    cCoins->Print(Form("%s/%d_barCoins.png",saveDir,nBlocks));
    cCoins->Print(Form("%s/%d_barCoins.pdf",saveDir,nBlocks));
    TCanvas *cEff = new TCanvas("cEff");
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

	  } // for (int sp=ientry; sp<tof1CoinChain->GetEntries(); sp++) 
	} // if (tof1Coin->lastDelayedBeamSignal != lastSpill && itdc == 0)
	// Need to calculate total signal hits here
	// Weight these by the efficiency of calculated above
	double deltat = TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1]  );
	double dstofHitT = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (10. - TMath::Abs(deltat) / 2. );
	double tofCalc = dstofHitT - tofCoin->usTofSignal;
	if (tofCalc < 200. && tofCalc > 70. && tofCoin->bar != 10) {
	  hdtof1d->Fill(tofCalc, 1. / hEff->GetBinContent(tofCoin->bar));
	  hdtof1d_noeff->Fill(tofCalc);
	} // if (tofCalc < 200. && tofCalc > 50.) 
	// Calculate position of the protons and pions with calculated efficiency
	if (tofCalc > piLow && tofCalc < piHi) {
	  if (tofCoin->bar != 10) {
	    hpiHitsDstofVert->Fill(tofCoin->bar, 1. / hEff->GetBinContent(tofCoin->bar));
	    hpiHitsDstofHorz->Fill((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.)+70., 1. / hEff->GetBinContent(tofCoin->bar));
	    hpiHitsDstofVert_noeff->Fill(tofCoin->bar);
	    hpiHitsDstofHorz_noeff->Fill((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.)+70.);
	    nPi += (1. / hEff->GetBinContent(tofCoin->bar));
	    nPi_noeff++;
	  } // if (tof1Coin->bar != 10)
	} // if (tofCalc > piLow && tofCalc < piHi) 
	else if (tofCalc > proLow && tofCalc < proHi) {
	  if (tofCoin->bar != 10) {
	    hproHitsDstofVert->Fill(tofCoin->bar, 1. / hEff->GetBinContent(tofCoin->bar));
	    hproHitsDstofHorz->Fill((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.)+70., 1. / hEff->GetBinContent(tofCoin->bar));
	    hproHitsDstofVert_noeff->Fill(tofCoin->bar);
	    hproHitsDstofHorz_noeff->Fill((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.)+70.);
	    nP += (1. / hEff->GetBinContent(tofCoin->bar));
	    nP_noeff++;
	  } // if (tof1Coin->bar != 10) 
	} // else if (tofCalc > proLow && tofCalc < proHi) 
      } // for (int h=0; h<tofCoinChain->GetEntries(); h++) 
    } // for (int itdc=0; itdc<2; itdc++)

    if (nBlocks==0) {
      hdtof1d->SetLineColor(kRed);
      hdtof1d_noeff->SetLineColor(kRed);
    }
    else if (nBlocks==1) {
      hdtof1d->SetLineColor(kBlue);
      hdtof1d_noeff->SetLineColor(kBlue);
    }
    else if (nBlocks==2) {
      hdtof1d->SetLineColor(kBlack);
      hdtof1d_noeff->SetLineColor(kBlack);
    }
    else if (nBlocks==3) {
      hdtof1d->SetLineColor(kGreen+2);
      hdtof1d_noeff->SetLineColor(kGreen+2);
    }
    else if (nBlocks==4) {
      hdtof1d->SetLineColor(kMagenta);
      hdtof1d_noeff->SetLineColor(kMagenta);
    }

    fout->cd(0);
    TCanvas *c2d = new TCanvas(Form("%d_c2d",nBlocks));
    c2d->SetLogy();
    hdtof1d->Fit(sPi, "R");
    hdtof1d->Fit(sPro, "R");
    hdtof1d->Fit(fBkg, "R");

    Double_t par[7];
    sPro->GetParameters(&par[0]);
    sPi->GetParameters(&par[3]);
    fBkg->GetParameters(&par[6]);
    fSplusB->SetParameters(par);
    hdtof1d->Fit(fSplusB, "R"); 

    hdtof1d->Draw("hist");
    sPro->Draw("same");
    sPi->Draw("same");
    fBkg->Draw("same");
    fSplusB->Draw("same");

    c2d->Print(Form("%s/%d_dtofFitted.png",saveDir,nBlocks));
    c2d->Print(Form("%s/%d_dtofFitted.pdf",saveDir,nBlocks));

    TCanvas *c2d_lin = new TCanvas(Form("%d_c2d_lin",nBlocks));
    c2d_lin->SetLogy();
    hdtof1d->Fit(fBkgLin, "R");
    Double_t parLin[8];
    sPro->GetParameters(&parLin[0]);
    sPi->GetParameters(&parLin[3]);
    fBkgLin->GetParameters(&parLin[6]);
    fSplusBLin->SetParameters(parLin);
    hdtof1d->Fit(fSplusBLin, "R");

    hdtof1d->Draw("hist");

    sPro->Draw("same");
    sPi->Draw("same");
    fBkgLin->Draw("same");
    fSplusBLin->Draw("same");

    c2d_lin->Print(Form("%s/%d_dtofFittedLinear.png",saveDir,nBlocks));
    c2d_lin->Print(Form("%s/%d_dtofFittedLinear.pdf",saveDir,nBlocks));

    TCanvas *c2d_exp = new TCanvas(Form("%d_c2d_exp",nBlocks));
    c2d_exp->SetLogy();
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
    c2d_exp->Print(Form("%s/%d_dtofFittedExp.png",saveDir,nBlocks));
    c2d_exp->Print(Form("%s/%d_dtofFittedExp.pdf",saveDir,nBlocks));

    TCanvas *c2d_exp_noeff = new TCanvas(Form("%d_c2d_exp_noeff",nBlocks));
    c2d_exp_noeff->SetLogy();
    hdtof1d_noeff->Fit(fBkgExp_noeff, "R");
    Double_t parExp_noeff[8];
    sPro->GetParameters(&parExp_noeff[0]);
    sPi->GetParameters(&parExp_noeff[3]);
    fBkgExp_noeff->GetParameters(&parExp_noeff[6]);
    fSplusBExp_noeff->SetParameters(parExp_noeff);
    hdtof1d_noeff->Fit(fSplusBExp_noeff, "R");
    hdtof1d_noeff->Draw("hist");
    fSplusBExp_noeff->Draw("same");
    hdtof1d_noeff->Write();
    c2d_exp_noeff->Print(Form("%s/%d_dtofFittedExp_noeff.png",saveDir,nBlocks));
    c2d_exp_noeff->Print(Form("%s/%d_dtofFittedExp_noeff.pdf",saveDir,nBlocks));

    // Now we have the fit values, loop over again and subtract the background
    // For each bin, find the fraction of each particle type which are background and 
    // this fraction from the bin
    TF1 *fSub = new TF1("fSub", "exp([0]+[1]*x)", 70, 200);
    fSub->SetParameter(0, fSplusBExp->GetParameter("bkgconst"));  
    fSub->SetParameter(1, fSplusBExp->GetParameter("bkgdecay"));
    // Background subtracted tof spectrum
    TH1D *hdtof1d_sub = (TH1D*)hdtof1d->Clone(Form("%d_hdtof1d_sub",nBlocks));
    hdtof1d_sub->Add(fSub, -1.);
    TCanvas *cdtof_sub = new TCanvas(Form("%d_cdtof_sub",nBlocks));
    hdtof1d_sub->Draw("hist");
    hdtof1d_sub->Write();    
    cdtof_sub->Print(Form("%s/%d_dtof1d_sub.png",saveDir,nBlocks));
    cdtof_sub->Print(Form("%s/%d_dtof1d_sub.pdf",saveDir,nBlocks));
    // No efficiency calculation
    TF1 *fSub_noeff = new TF1("fSub_noeff", "exp([0]+[1]*x)", 70, 200);
    fSub_noeff->SetParameter(0, fSplusBExp_noeff->GetParameter("bkgconst"));  
    fSub_noeff->SetParameter(1, fSplusBExp_noeff->GetParameter("bkgdecay"));
    // Background subtracted tof spectrum
    TH1D *hdtof1d_sub_noeff = (TH1D*)hdtof1d_noeff->Clone(Form("%d_hdtof1d_sub_noeff",nBlocks));
    hdtof1d_sub_noeff->Add(fSub_noeff, -1.);
    TCanvas *cdtof_sub_noeff = new TCanvas(Form("%d_cdtof_sub_noeff",nBlocks));
    hdtof1d_sub_noeff->Draw("hist");
    hdtof1d_sub_noeff->Write();    
    cdtof_sub_noeff->Print(Form("%s/%d_dtof1d_sub_noeff.png",saveDir,nBlocks));
    cdtof_sub_noeff->Print(Form("%s/%d_dtof1d_sub_noeff.pdf",saveDir,nBlocks));

    double proSubFrac = fSub->Integral(proLow, proHi) / nP;
    double piSubFrac  = fSub->Integral(piLow, piHi) / nPi;
    cout<<"Pro bkg frac "<<proSubFrac<<" Pi bkg frac "<<piSubFrac<<endl;
    for (int i=0; i <= hproHitsDstofVert->GetNbinsX(); i++) {
      hproHitsDstofVert->SetBinContent(i, hproHitsDstofVert->GetBinContent(i) * (1. - proSubFrac));
      hpiHitsDstofVert->SetBinContent(i, hpiHitsDstofVert->GetBinContent(i) * (1. - piSubFrac));
    } // for (int i=0; i <= hproHitsDstofVert->GetNBins(); i++) 
    for (int j=0; j <= hproHitsDstofHorz->GetNbinsX(); j++) {
      hproHitsDstofHorz->SetBinContent(j, hproHitsDstofHorz->GetBinContent(j) * (1. - proSubFrac));
      hpiHitsDstofHorz->SetBinContent(j, hpiHitsDstofHorz->GetBinContent(j) * (1. - piSubFrac));
    } // for (int j=0; j <= hproHitsDstofHorz->GetNBins(); j++) 
    // No efficiency
    double proSubFrac_noeff = fSub_noeff->Integral(proLow, proHi) / nP_noeff;
    double piSubFrac_noeff  = fSub_noeff->Integral(piLow, piHi) / nPi_noeff;
    cout<<"Pro bkg frac "<<proSubFrac<<" Pi bkg frac "<<piSubFrac<<endl;
    for (int i=0; i <= hproHitsDstofVert_noeff->GetNbinsX(); i++) {
      hproHitsDstofVert_noeff->SetBinContent(i, hproHitsDstofVert_noeff->GetBinContent(i) * (1. - proSubFrac_noeff));
      hpiHitsDstofVert_noeff->SetBinContent(i, hpiHitsDstofVert_noeff->GetBinContent(i) * (1. - piSubFrac_noeff));
    } // for (int i=0; i <= hproHitsDstofVert->GetNBins(); i++) 
    for (int j=0; j <= hproHitsDstofHorz_noeff->GetNbinsX(); j++) {
      hproHitsDstofHorz_noeff->SetBinContent(j, hproHitsDstofHorz_noeff->GetBinContent(j) * (1. - proSubFrac_noeff));
      hpiHitsDstofHorz_noeff->SetBinContent(j, hpiHitsDstofHorz_noeff->GetBinContent(j) * (1. - piSubFrac_noeff));
    } // for (int j=0; j <= hproHitsDstofHorz->GetNBins(); j++) 
    
    TH1D *hproPiRatioHorz = new TH1D(Form("hproPiRatioHorz_%d",nBlocks), Form("Proton/MIP, %d blocks; x / cm; Proton/MIP",nBlocks), 20, 0., 140.);
    TH1D *hproPiRatioVert = new TH1D(Form("hproPiRatioVert_%d",nBlocks), Form("Proton/MIP, %d blocks; Bar; Proton/MIP",nBlocks), 10, 0.5, 10.5);
    TH1D *hproPiRatioHorz_noeff = new TH1D(Form("hproPiRatioHorz_noeff_%d",nBlocks), Form("Proton/MIP, %d blocks (no efficiency); x / cm; Proton/MIP",nBlocks), 20, 0., 140.);
    TH1D *hproPiRatioVert_noeff = new TH1D(Form("hproPiRatioVert_noeff_%d",nBlocks), Form("Proton/MIP, %d blocks (no efficiency); Bar; Proton/MIP",nBlocks), 10, 0.5, 10.5);
    hproPiRatioHorz->Divide(hproHitsDstofHorz, hpiHitsDstofHorz, 1., 1., "B");
    hproPiRatioVert->Divide(hproHitsDstofVert, hpiHitsDstofVert, 1., 1., "B");
    hproPiRatioHorz_noeff->Divide(hproHitsDstofHorz_noeff, hpiHitsDstofHorz_noeff, 1., 1., "B");
    hproPiRatioVert_noeff->Divide(hproHitsDstofVert_noeff, hpiHitsDstofVert_noeff, 1., 1., "B");

    fout->cd(0);
    new TCanvas;
    hproPiRatioHorz->Draw("hist e");
    hproPiRatioHorz->Write();
    new TCanvas;
    hproPiRatioVert->Draw("hist e");
    hproPiRatioVert->Write();

    new TCanvas;
    hproPiRatioHorz_noeff->Draw("hist e");
    hproPiRatioHorz_noeff->Write();
    new TCanvas;
    hproPiRatioVert_noeff->Draw("hist e");
    hproPiRatioVert_noeff->Write();

    new TCanvas;
    hproHitsDstofHorz->Draw("hist e");
    hproHitsDstofHorz->Write();
    new TCanvas;
    hproHitsDstofVert->Draw("hist e");
    hproHitsDstofVert->Write();
    new TCanvas;
    hpiHitsDstofHorz->Draw("hist e");
    hpiHitsDstofHorz->Write();
    new TCanvas;
    hpiHitsDstofVert->Draw("hist e");
    hpiHitsDstofVert->Write();

    new TCanvas;
    hproHitsDstofHorz_noeff->Draw("hist e");
    hproHitsDstofHorz_noeff->Write();
    new TCanvas;
    hproHitsDstofVert_noeff->Draw("hist e");
    hproHitsDstofVert_noeff->Write();
    new TCanvas;
    hpiHitsDstofHorz_noeff->Draw("hist e");
    hpiHitsDstofHorz_noeff->Write();
    new TCanvas;
    hpiHitsDstofVert_noeff->Draw("hist e");
    hpiHitsDstofVert_noeff->Write();

  } // nBlocks

} // bkgSub


