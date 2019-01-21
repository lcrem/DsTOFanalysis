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
  // Most runs were in this configuration so don't need to use necessarily
  const double start4Block = 1535608220;
  const double end4Block   = 1535617102;
  
  // Timing cuts
  const double piLow  = 80.;
  const double piHi   = 95.;
  const double proLow = 106.;
  const double proHi  = 134.;

  for (int nBlocks=0; nBlocks < 4; nBlocks++) {
    double nSpills = 0.;
    double nSpillsTrue = 0.;
    double lastSpill = 0.;

    TH1D *hproton = new TH1D(Form("hproton_%d", nBlocks), Form("%d blocks: Number of protons; Time since spill start / ns; P", nBlocks), 40, 0, 1e9);
    TH1D *hpion = new TH1D(Form("hpion_%d", nBlocks), Form("%d blocks: Number of pions; Time since spill start / ns; #pi", nBlocks), 40, 0, 1e9);
    TH1D *hratio = new TH1D(Form("hratio_%d", nBlocks), Form("%d blocks: Proton/(#pi+#mu); Time since spill start / ns; P / (#pi+#mu)", nBlocks), 40, 0, 1e9);
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

	if (tof > piLow && tof < piHi /*&& tofCoin->bar!=10*/) {
	  hpion->Fill(dstofHitT - tofCoin->lastDelayedBeamSignal);
	}
	else if (tof > proLow && tof < proHi /*&& tofCoin->bar!=10*/) {
	  hproton->Fill(dstofHitT - tofCoin->lastDelayedBeamSignal);
	}
      }
    }
    gStyle->SetOptStat(0);

    TCanvas *c1 = new TCanvas(Form("c1_%d", nBlocks));
    hratio->Divide(hproton, hpion, 1., 1., "B");
    if (nBlocks == 0) {
      hratio->SetBinContent(25, 0);
    }
    hratio->Draw("hist E");
    c1->Print(Form("../nBlocksPlots/%dblocks_proPiSpillRatio.png", nBlocks));
    c1->Print(Form("../nBlocksPlots/%dblocks_proPiSpillRatio.pdf", nBlocks));

  } 
}
