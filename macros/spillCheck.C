// spillCheck.C

// Checks a if an umatched spill in a given time frame is empty
void spillCheck (const int timeStart , const int timeEnd) {
  TFile *usDB = new TFile("/unix/dune/sjones/utofBeamSpills.root");
  TFile *dsDB = new TFile("/unix/dune/sjones/dtofBeamSpills.root");
  TTree *dsTree = (TTree*)dsDB->Get("beamTree");
  TTree *usTree = (TTree*)usDB->Get("utofBeamTree");

  double fakeUnixTimeUs;
  double fakeUnixTimeDs;
  double nsTimeUs;
  double nsTimeDs;
  int unixRunStartUs;
  int unixRunStartDs;
  int run;
  dsTree->SetBranchAddress("fakeUnixNs", &fakeUnixTimeDs);
  dsTree->SetBranchAddress("unixRunStart", &unixRunStartDs);
  dsTree->SetBranchAddress("nsTime", &nsTimeDs);
  dsTree->SetBranchAddress("run", &run);
  usTree->SetBranchAddress("fakeUnixTime", &fakeUnixTimeUs);
  usTree->SetBranchAddress("unixRunStart", &unixRunStartUs);
  usTree->SetBranchAddress("nsTime", &nsTimeUs);

  std::vector< pair<double, int> > unmatchedVec;
  
  int totalSpills = 0;
  int unmatched = 0;
  for (int us=0; us<usTree->GetEntries(); us++) {
    usTree->GetEntry(us);
    // Get spills in time range
    if (fakeUnixTimeUs >= timeStart && fakeUnixTimeUs <= timeEnd) {
      totalSpills++;
      double lowestDiff = 99999999.;
      int lowestRun = 0;
      double lowestNs = 0.;
      double lowestNsDs = 0.;
      double lowestNsUs = 0.;
      double lowestfakeUnixDs = 0.;
      double lowestfakeUnixUs = 0.;
      double lowestunixRunStartDs = 0.;
      for (int ds=0; ds<dsTree->GetEntries(); ds++) {
	dsTree->GetEntry(ds);
	if (fakeUnixTimeDs >= timeStart && fakeUnixTimeDs <=timeEnd) {
	  double diff = fakeUnixTimeDs - fakeUnixTimeUs;
	  if (abs(diff) < abs(lowestDiff)) {
	    lowestDiff = diff;
	    lowestRun  = run;
	    lowestNsUs = nsTimeUs;
	    //	    lowestNsDs = nsTimeDs;
	    lowestunixRunStartDs = unixRunStartDs;
	    lowestfakeUnixUs = fakeUnixTimeUs;
	    lowestfakeUnixDs = fakeUnixTimeDs;
	    // Take the unix time of the unmatched ustof spill and subtract the start of the dstof run within which it is contained
	    lowestNsDs = (fakeUnixTimeUs - unixRunStartDs) * 1e9;
	    // Work out when the unmatched ustof signal is in the dstof file
	    
	  }
	} 
      } // dstof spills
      // Is an umatched spill
      if (abs(lowestDiff) > 2) {
	unmatched++;
	cout.precision(10);
	cout<<"Unmatched spill in dstof run "<<lowestRun<<", unix time "<<fakeUnixTimeUs<<", "<<lowestNsDs<<" from run start, with difference "<<lowestDiff<<endl;
	unmatchedVec.push_back(make_pair(lowestNsDs, lowestRun));
      }
    }
  } // ustof spills
  double perc = (float(unmatched) / float(totalSpills)) * 100;
  std::cout<<unmatched<<" unmatched spills found out of "<<totalSpills<< " ("<<perc<<"%)"<<endl;
  std::cout<<unmatchedVec.size()<<std::endl;
  // Now look up if there are spills in the dtof files
  for (int i=0; i<unmatchedVec.size(); i++) {
    // Find and open file
    for (int file = 1; file <= 1442; file++) {
      if (file == unmatchedVec[i].second) {
	TH1D *htof = new TH1D("htof", Form("Dstof hits in missing spill %d; Time / ns", i), 250, 0, 1e9);
	TFile *dtoffile = new TFile(Form("/unix/dune/hptpctof/run%d/DsTOFcoincidenceRun%d_tdc1.root", file, file), "read");
	TTree *dtofTree = (TTree*)dtoffile->Get("tofCoinTree");
	RawDsTofCoincidence *tempcoin = NULL;
	dtofTree->SetBranchAddress("tofCoin", &tempcoin);

	// Plot the events that are in the missing spill
	for (int t=0; t<dtofTree->GetEntries(); t++) {
	  dtofTree->GetEntry(t);
	  if (tempcoin->fakeTimeNs[0] >= unmatchedVec[i].first && tempcoin->fakeTimeNs[0] <= unmatchedVec[i].first + 1e9) {
	    htof->Fill(tempcoin->fakeTimeNs[0] - unmatchedVec[i].first);
	  }
	}
	TCanvas *c = new TCanvas(Form("c%d", i));
	htof->Draw("hist");
	c->Print(Form("unmatchedspill%d.png", i));
	c->Print(Form("unmatchedspill%d.pdf", i));
	delete tempcoin;
	delete dtofTree;
	delete dtoffile;
      }
    }
  }
}
