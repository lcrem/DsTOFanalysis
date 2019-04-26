// proPiFrac.C
// Macro for checking proton/pion fraction across the length of the spill

// Takes 2 coincidence files
void proPiFrac(const char* dstofDir="/scratch0/dbrailsf/temp/mylinktodtof/") {

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
  const double end3Block   = 1535798437;
  // 0.8GeV/c, 4 block
  // 4 moderator blocks with -4cm bend
  const double start4Block = 1535836129;
  const double end4Block   = 1535879634;  

  // Timing cuts
  const double piLow  = 80.;
  const double piHi   = 95.;
  const double proLow = 106.;
  const double proHi  = 145.;

  THStack *hsproton = new THStack("hsproton", "Proton profile during beam spill in S4; Time since beam spill start / ns; Events / spill");
  THStack *hspion   = new THStack("hspion", "MIP profile during beam spill in S4; Time since beam spill start / ns; Events / spill"); 
  THStack *hsratio  = new THStack("hsratio", "Proton/MIP ratio during beam spill in S4; Time since beam spill start / ns; Protons/MIPs");

  for (int nBlocks=0; nBlocks <= 4; nBlocks++) {
    double nSpills = 0.;
    double nSpillsTrue = 0.;
    double lastSpill = 0.;

    TH1D *hproton = new TH1D(Form("hproton_%d", nBlocks), Form("%d blocks: Number of protons; Time since spill start / ns; Events / spill", nBlocks), 40, 0, 1e9);
    TH1D *hpion = new TH1D(Form("hpion_%d", nBlocks), Form("%d blocks: Number of pions; Time since spill start / ns; Events / spill", nBlocks), 40, 0, 1e9);
    TH1D *hratio = new TH1D(Form("hratio_%d", nBlocks), Form("%d blocks: Proton/MIP; Time since spill start / ns; Proton / MIP", nBlocks), 40, 0, 1e9);
    hproton->Sumw2();
    hpion->Sumw2();
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

    for (int irun=950; irun<1200; irun++){
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

    // Need these to calculate bar efficiencies
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
      for (int irun=runMin; irun<runMax+1; irun++) {
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
	TTree *ustofTree = new TTree("ustofTree", "ustof");
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
	} // for (int t=0; t<tofCoinTree->GetEntries(); t++) 
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
	delete tofTree;
	delete tofCoinTree;
	delete ustofTree;
	tofFile->Close();
	tofCoinFile->Close();
      } // for (int irun=runMin; irun<runMax+1; irun++) 
    } // for (int itdc=0; itdc<2; itdc++)
    fout->cd();
    TCanvas *cHits = new TCanvas(Form("%d_cHits",nBlocks));
    hHits->Draw("hist e");
    hHits->Write();
    cHits->Print(Form("%s/%d_barHits.png",saveDir,nBlocks));
    cHits->Print(Form("%s/%d_barHits.pdf",saveDir,nBlocks));
    TCanvas *cCoins = new TCanvas(Form("%d_cCoins",nBlocks));
    hCoins->Draw("hist e");
    hCoins->Write();
    cCoins->Print(Form("%s/%d_barCoins.png",saveDir,nBlocks));
    cCoins->Print(Form("%s/%d_barCoins.pdf",saveDir,nBlocks));
    TCanvas *cEff = new TCanvas(Form("%d_cEff",nBlocks));
    hHits->Add(hCoins, -1.);
    hEff->Divide(hCoins, hHits, 1., 1., "B");
    hEff->Draw("hist e");
    hEff->Write();
    cEff->Print(Form("%s/%d_barEff.png",saveDir,nBlocks));
    cEff->Print(Form("%s/%d_barEff.pdf",saveDir,nBlocks));

    // Now loop over again and make the actual plots
    for (int itdc=0; itdc<2; itdc++) {
      TChain *tofCoinChain = new TChain("tofCoinTree");
      for (int irun=runMin; irun<runMax+1; irun++){
	tofCoinChain->Add(Form("%srun%d/DsTOFcoincidenceRun%d_tdc%d.root", dstofDir, irun, irun, itdc+1));
      }
      RawDsTofCoincidence *tofCoin = NULL;
      tofCoinChain->SetBranchAddress("tofCoin", &tofCoin);

      double dstofHitT, deltat;

      for (int ientry=0; ientry<tofCoinChain->GetEntries(); ientry++){
	tofCoinChain->GetEntry(ientry);
	if (tofCoin->unixTime[0]<startTime) continue;
	if (tofCoin->unixTime[0]>endTime) break;
	
	if (tofCoin->lastDelayedBeamSignal != lastSpill && itdc==0) {
	  lastSpill = tofCoin->lastDelayedBeamSignal;
	  nSpills++;
	  for (int sp=ientry; sp < tofCoinChain->GetEntries(); sp++) {
	    tofCoinChain->GetEntry(sp);
	    if ((tofCoin->fakeTimeNs[0] - tofCoin->usTofSignal) < 200. && 
		(tofCoin->fakeTimeNs[0] - tofCoin->usTofSignal) > 70. &&
		(tofCoin->fakeTimeNs[0] - lastSpill) < 1e9 &&
		(tofCoin->fakeTimeNs[0] - lastSpill) > 0) {
	      nSpillsTrue++;
	      break;
	    }

	    if (tofCoin->unixTime[0]>endTime) break;
	  }
	}      

	deltat = TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1]  );
	dstofHitT = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (10. - TMath::Abs(deltat) / 2 );
	double tof = dstofHitT - tofCoin->usTofSignal;

	if (tof > piLow && tof < piHi && tofCoin->bar!=10) {
	  hpion->Fill(dstofHitT - tofCoin->lastDelayedBeamSignal, hEff->GetBinContent(tofCoin->bar));
	}
	else if (tof > proLow && tof < proHi && tofCoin->bar!=10) {
	  hproton->Fill(dstofHitT - tofCoin->lastDelayedBeamSignalm , hEff->GetBinContent(tofCoin->bar));
	}
      }
    }
    gStyle->SetOptStat(0);


    hratio->Divide(hproton, hpion, 1., 1., "B");
    hpion->Scale(1. / nSpillsTrue);
    hproton->Scale(1. / nSpillsTrue);
    if (nBlocks == 0) {
      hratio->SetLineColor(kBlue);
      hpion->SetLineColor(kBlue);
      hproton->SetLineColor(kBlue);
    }
    if (nBlocks == 1) {
      hratio->SetLineColor(kRed);
      hpion->SetLineColor(kRed);
      hproton->SetLineColor(kRed);
    }
    if (nBlocks == 2) {
      hratio->SetLineColor(kBlack);
      hpion->SetLineColor(kBlack);
      hproton->SetLineColor(kBlack);
    }
    if (nBlocks == 3) {
      hratio->SetLineColor(kGreen+2);
      hpion->SetLineColor(kGreen+2);
      hproton->SetLineColor(kGreen+2);
    }
    if (nBlocks == 4) {
      hratio->SetLineColor(kMagenta);
      hpion->SetLineColor(kMagenta);
      hproton->SetLineColor(kMagenta);
    }
    TCanvas *c1 = new TCanvas(Form("c1_%d", nBlocks));
    hratio->Draw("hist E");
    c1->Print(Form("%s/%dblocks_proPiSpillRatio.png", saveDir, nBlocks));
    c1->Print(Form("../nBlocksPlots/%dblocks_proPiSpillRatio.pdf", nBlocks));

  } 
}
