// beamTree.C


void beamTree () {
  gSystem->Load("libdstof.so");

  TTree *beamTree = new TTree("beamTree", "HPTPC beam spills");
  beamTree->SetDirectory(0);

  int run;
  int unixRunStart;
  double nsTime; // Time since run start
  double fakeUnixTime; // Unix + some time in ns


  beamTree->Branch("run", &run);
  beamTree->Branch("unixRunStart", &unixRunStart);
  beamTree->Branch("nsTime", &nsTime);
  beamTree->Branch("fakeUnixTime", &fakeUnixTime);

  int nSpills = 0;
  int nSpillsTrue = 0;
  // Loop over all Dstof files
  for (int file = 1; file <= 1442; file++) {
    std::cout<<"Opening run "<<file<<std::endl;
    string filename = Form("/unix/dune/hptpctof/run%d/DsTOFcoincidenceRun%d_tdc1.root", file, file);
    if(!gSystem->AccessPathName(filename.c_str())) {

      TFile *inFile = TFile::Open(filename.c_str(), "read");
      
      TTree *coinTree = (TTree*)inFile->Get("tofCoinTree");
      RawDsTofCoincidence *tempcoin = NULL;
      coinTree->SetBranchAddress("tofCoin", &tempcoin);

      coinTree->GetEntry(0);
      int runstart = tempcoin->unixTime[0];
      
      double lastdelayed = 0.;
      double lastdelayed2=0.;
      
      std::cout<<"Run has "<<coinTree->GetEntries()<<" entries"<<std::endl;
      // Now loop over all events and find the correct beam spills
      for (int t = 0; t < coinTree->GetEntries(); t++) {
	coinTree->GetEntry(t);
	// Selects the right kind of beam signals
	if (tempcoin->lastDelayedBeamSignal > 0. && tempcoin->lastDelayedBeamSignal != lastdelayed && tempcoin->lastDelayedBeamSignal - tempcoin->lastRawBeamSignal > 0.5e9 && tempcoin->lastDelayedBeamSignal - tempcoin->lastRawBeamSignal < 1e9) {
	  nSpills++;
	  bool hasWritten = false;
	  // Require there to be at least one coincidence between S1 and S4 (i.e. beam stopper out) 
	  // in the next second
	  lastdelayed = tempcoin->lastDelayedBeamSignal;
	  for (int sp = t; sp < coinTree->GetEntries(); sp++) {
	    coinTree->GetEntry(sp);
	    if ((tempcoin->fakeTimeNs[0] - tempcoin->usTofSignal) < 200. &&  (tempcoin->fakeTimeNs[0] - tempcoin->usTofSignal) > 0. && (tempcoin->fakeTimeNs[0] - lastdelayed) < 1e9 && (tempcoin->fakeTimeNs[0] - lastdelayed) > 0) {//&& !hasWritten) {
	      nSpillsTrue++;
	      // Set new tree vars
	      run = tempcoin->run;
	      unixRunStart = runstart;
	      nsTime = tempcoin->lastDelayedBeamSignal;
	      fakeUnixTime = runstart + (tempcoin->lastDelayedBeamSignal / 1e9);
      
	      beamTree->Fill();
	      // Fill once and then stop
	      break;
	      //	      hasWritten = true;
	    }
	  }
	}
      }
      std::cout<<"Have found "<<nSpills<<" total"<<" and "<<nSpillsTrue<<" with stopper out"<<std::endl;
      delete tempcoin;
      delete coinTree;
      delete inFile;
    } 
  }
  TFile *fout = new TFile("/unix/dune/sjones/TOF/dtofBeamSpills.root", "recreate");
  beamTree->Write("dtofBeamTree");
  fout->Close();
}
