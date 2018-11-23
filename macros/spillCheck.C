// spillCheck.C

// Checks a if an umatched spill in a given time frame is empty
void spillCheck (const int timeStart , const int timeEnd) {

  assert(timeStart < timeEnd);

  TFile *usDB = new TFile("~/utofBeamTree/utofBeamSpills.root");
  TFile *dsDB = new TFile("~/utofBeamTree/dtofBeamSpills.root");
  TTree *dsTree = (TTree*)dsDB->Get("dtofBeamTree");
  TTree *usTree = (TTree*)usDB->Get("utofBeamTree");

  double fakeUnixTimeUs;
  double fakeUnixTimeDs;
  double nsTimeUs;
  double nsTimeDs;
  int unixRunStartUs;
  int unixRunStartDs;
  int run;
  dsTree->SetBranchAddress("fakeUnixTime", &fakeUnixTimeDs);
  dsTree->SetBranchAddress("unixRunStart", &unixRunStartDs);
  dsTree->SetBranchAddress("nsTime", &nsTimeDs);
  dsTree->SetBranchAddress("run", &run);
  usTree->SetBranchAddress("fakeUnixTime", &fakeUnixTimeUs);
  usTree->SetBranchAddress("unixRunStart", &unixRunStartUs);
  usTree->SetBranchAddress("nsTime", &nsTimeUs);

  vector< pair<double, int> > unmatchedVec;
  
  int totalSpills = 0;
  int unmatched = 0;
  for (int ds=0; ds<dsTree->GetEntries(); ds++) {
    dsTree->GetEntry(ds);
    // Get spills in time range
    if (fakeUnixTimeDs >= timeStart && fakeUnixTimeDs <= timeEnd) {
      totalSpills++;
      double lowestDiff = 99999999.;
      int lowestRun = 0;
      double lowestNs = 0.;
      double lowestNsDs = 0.;
      double lowestNsUs = 0.;
      double lowestfakeUnixDs = 0.;
      double lowestfakeUnixUs = 0.;
      double lowestunixRunStartDs = 0.;
      for (int us=0; us<usTree->GetEntries(); us++) {
	usTree->GetEntry(us);
	if (fakeUnixTimeUs >= timeStart && fakeUnixTimeUs <=timeEnd) {
	  double diff = fakeUnixTimeDs - fakeUnixTimeUs;
	  if (abs(diff) < abs(lowestDiff)) {
	    lowestDiff = diff;
	    lowestRun  = run;
	    lowestNsDs = nsTimeDs;
	    //	    lowestNsDs = nsTimeDs;
	    lowestunixRunStartUs = unixRunStartUs;
	    lowestfakeUnixDs = fakeUnixTimeDs;
	    lowestfakeUnixUs = fakeUnixTimeUs;
	    lowestNsUs = (fakeUnixTimeDs - unixRunStartUs) * 1e9;
	    // Work out when the unmatched dstof signal is in the ustof file
	    
	  }
	} 
      } // ustof spills
      // Is an umatched spill
      if (abs(lowestDiff) > 2) {
	unmatched++;
	cout.precision(10);
	cout<<"Unmatched spill in dstof run "<<lowestRun<<", unix time "<<fakeUnixTimeDs<<", "<<lowestNsUs<<" from run start, with difference "<<lowestDiff<<endl;
	unmatchedVec.push_back(make_pair(lowestNsUs, lowestRun));
      }
    }
  } // dstof spills
  double perc = (float(unmatched) / float(totalSpills)) * 100;
  std::cout<<unmatched<<" unmatched spills found out of "<<totalSpills<< " ("<<perc<<"%)"<<endl;

  // Now look up if there are spills in the utof files
  for (int i=0; i<unmatchedVec.size(); i++) {

  }
}
