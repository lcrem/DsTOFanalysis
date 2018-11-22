// spillCheck.C

// Checks a if an umatched spill in a given time frame is empty
void spillCheck (const int timeStart , const int timeEnd) {
  TFile *usDB = new TFile("utofBeamSpills.root");
  TFile *dsDB = new TFile("dtofBeamSpills.root");
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

  vector< pair<double, int> > unmatchedVec;
  
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

  // Now look up if there are spills in the dtof files
  for (int i=0; i<unmatchedVec.size(); i++) {

  }
}
