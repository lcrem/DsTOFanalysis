// beamTree.C



void beamTree () {
  gSystem->Load("libdstof.so");

  TTree *beamTree = new TTree("beamTree", "HPTPC beam spills");
  beamTree->SetDirectory(0);

  int run;
  int unixRunStart;
  double nsTime; // Time since run start
  double fakeUnixNs; // Unix + some time in ns


  beamTree->Branch("run", &run);
  beamTree->Branch("unixRunStart", &unixRunStart);
  beamTree->Branch("nsTime", &nsTime);
  beamTree->Branch("fakeUnixNs", &fakeUnixNs);

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
      
      std::cout<<"Run has "<<coinTree->GetEntries()<<" entries"<<std::endl;
      // Now loop over all events and find the correct beam spills
      for (int t = 0; t < coinTree->GetEntries(); t++) {
	coinTree->GetEntry(t);
	// Selects the right kind of beam signals
	if (tempcoin->lastDelayedBeamSignal > 0. && tempcoin->lastDelayedBeamSignal != lastdelayed && tempcoin->lastDelayedBeamSignal - tempcoin->lastRawBeamSignal > 9.00981e8 && tempcoin->lastDelayedBeamSignal - tempcoin->lastRawBeamSignal < 9.0098135e8) {
	  lastdelayed = tempcoin->lastDelayedBeamSignal;
	  //	std::cout<<"Beam signal"<<std::endl;
	  // Set new tree vars
	  run = tempcoin->run;
	  unixRunStart = runstart;
	  nsTime = tempcoin->lastDelayedBeamSignal;
	  fakeUnixNs = runstart + (tempcoin->lastDelayedBeamSignal / 1e9);

	  beamTree->Fill();
	}
      }

      delete tempcoin;
      delete coinTree;
      delete inFile;
    } 
  }
  TFile *fout = new TFile("/unix/dune/sjones/TOF/beamSpills.root", "recreate");
  beamTree->Write("beamTree");
  fout->Close();
}
