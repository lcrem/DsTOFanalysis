// angularDist.C
// Angular distribution of protons and pions for different moderator blocks
void angularDist (const int nBlocks) {
  // Unix timestamps for variable block moves
  // 0.8GeV/c, 0 blocks
  const double start0Block = 1535713289;
  const double end0Block   = 1535716132;
  // 0.8GeV/c, 1 block
  const double start1Block = 1535796057;
  const double end1Block   = 1535799112;
  // 0.8GeV/c, 2 blocks
  const double start2Block = 1535796057;
  const double end2Block   = 1535799112;
  // 0.8GeV/c, 3 block
  const double start3Block = 1535792404;
  const double end3Block   = 1535798437;
  // 0.8GeV/c, 4 block
  // Most runs were in this configuration so don't need to use necessarily
  const double start4Block = 1535608220;
  const double end4Block   = 1535617102;
  // Dstof directory
  const char* dsDir = "/scratch0/dbrailsf/temp/mylinktodtof/";
  // Timing cuts
  const double piLow  = 80.;
  const double piHi   = 95.;
  const double proLow = 116.;
  const double proHi  = 134.;

  TH1D *htof1d = new TH1D("htof1d", Form("Time of flight, %d blocks; DsToF - UsToF / ns; Events", nBlocks), 100, 50, 150);
  TH2D *piHitsDstof = new TH2D("piHitsDstof", Form("%d blocks, position of S4 #pi hits; x / cm; y / cm", nBlocks), 60, 0, 140, 10, 0, 77.5);
  TH2D *proHitsDstof = new TH2D("proHitsDstof", Form("%d blocks, position of S4 proton hits; x / cm; y / cm", nBlocks), 60, 0, 140, 10, 0, 77.5);
  TH2D *proPiDstof = new TH2D("proPiDstof", Form("%d blocks, proton/pion ratio in S4; x / cm; y / cm", nBlocks), 60, 0, 140, 10, 0, 77.5);
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
    TFile *fin = new TFile(Form("%srun%d/DsTOFcoincidenceRun%d_tdc1.root", dsDir, irun, irun), "read");
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
      tofCoinChain->Add(Form("%srun%d/DsTOFcoincidenceRun%d_tdc%d.root", dsDir, irun, irun, itdc+1));
    }

    RawDsTofCoincidence *tofCoin = NULL;
    tofCoinChain->SetBranchAddress("tofCoin", &tofCoin);

    double dstofHitT, deltat;

    for (int ientry=0; ientry<tofCoinChain->GetEntries(); ientry++){
      tofCoinChain->GetEntry(ientry);
      if (tofCoin->unixTime[0]<startTime) continue;
      if (tofCoin->unixTime[0]>endTime) break;
      
      
      deltat = TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1]  );
      dstofHitT = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (10. - TMath::Abs(deltat) / 2 );
      double tof = dstofHitT - tofCoin->usTofSignal;
      htof1d->Fill(dstofHitT - tofCoin->usTofSignal);
      if (tof > piLow && tof < piHi) {
	piHitsDstof->Fill((tofCoin->fakeTimeNs[0] - tofCoin->fakeTimeNs[1])*(7./2.)+70., (tofCoin->bar*7.5) - 2.5);
      }
      else if (tof > proLow && tof < proHi) {
	proHitsDstof->Fill((tofCoin->fakeTimeNs[0] - tofCoin->fakeTimeNs[1])*(7./2.)+70., (tofCoin->bar*7.5) - 2.5);
      }
    }
  }
  new TCanvas;
  proPiDstof->Divide(proHitsDstof, piHitsDstof);
  proPiDstof->Draw("colz");
  new TCanvas;
  proHitsDstof->Draw("colz");
  new TCanvas;
  piHitsDstof->Draw("colz");
  new TCanvas;
  htof1d->Draw("hist");
}
