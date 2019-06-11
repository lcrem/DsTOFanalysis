// deadtimeNoRain.C
void deadtimeTimestamp(const char* utofFile,
		       const char *outDir,	    
		       const char* dstofDir="/scratch0/dbrailsf/temp/mylinktodtof",
		       const char* ustofDir="/zfs_home/sjones/mylinktoutof") 
{
  gROOT->SetBatch(kTRUE);
  TFile *fout = new TFile(Form("%s/%s/%s_out.root", outDir, utofFile, utofFile), "recreate");

  // Ustof-dstof cable delay
  const double ustofDelay = 184.7;
  // Apparent cut off time from spill profiles
  const double cutOffT = 0.29;

  double startTime = 0;
  double endTime   = 0;

  // Go through utof files and count the number of hits
  TFile *futof = new TFile(Form("%s/%s.root", ustofDir, utofFile), "read");

  double tToF[50];
  double tTrig;
  double tS1;
  double tSoSd;
  int nhit;
  int nBar[50];

  TTree *tree = (TTree*)futof->Get("tree");

  tree->SetBranchAddress("nhit", &nhit);
  tree->SetBranchAddress("tS1", &tS1);
  tree->SetBranchAddress("tToF", tToF);
  tree->SetBranchAddress("tTrig", &tTrig);
  tree->SetBranchAddress("tSoSd", &tSoSd);
  tree->SetBranchAddress("nBar", nBar);

  double lastS1S2utof = 0.;

  double lastUtofSpill = 0.;
  int nUtofSpills = 0;
  int nDtofSpills = 0;

  std::vector<double> utofTimes;
  std::vector<double> dtofTimes;

  // Get the start and end time of the file
  TNamed *start = 0;
  TNamed *end   = 0;
  futof->GetObject("start_of_run", start);
  futof->GetObject("end_of_run", end);

  const char* startchar = start->GetTitle();
  std::string startstr(startchar);
  std::string unixstart = startstr.substr(25,10);
  startTime = stoi(unixstart);
  
  const char* endchar = end->GetTitle();
  std::string endstr(endchar);
  std::string unixend = endstr.substr(23,10);
  endTime = stoi(unixend);

  cout.precision(13);
  cout<<"Utof file start, end "<<startTime<<", "<<endTime<<endl;

  Int_t runMin=-1;
  Int_t runMax=-1;

  // Find the correct dtof files
  for (int irun=900; irun<1400; irun++) {
    TFile *fin = new TFile(Form("%s/run%d/DsTOFcoincidenceRun%d_tdc1.root", dstofDir, irun, irun), "read");
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
  } // for (int irun=950; irun<1400; irun++) 
    
  cout << "Min and max runs are " << runMin << " " << runMax << endl;

  TTree *hitSpillTreeD = new TTree("hitSpillTreeD", "Spill hits");
  TTree *hitSpillTreeU = new TTree("hitSpillTreeU", "Spill hits");
  int dtofHits;
  int utofHits;
  hitSpillTreeD->Branch("dtofHits", &dtofHits);
  hitSpillTreeD->Branch("utofHits", &utofHits);

  for (int irun = runMin; irun < runMax+1; irun++) {
    cout<<"Getting hits for dtof run "<<irun<<endl;

    TFile *dbFile = new TFile(Form("/scratch2/sjones/spillDB/spillDB_run%d_run%d.root", irun, irun), "read");
    std::vector<double> tempDtofTimes;
    TTree *spillTree = (TTree*)dbFile->Get("spillTree");
    double globalSpillTime;
    double ustofSpillTime;
    spillTree->SetBranchAddress("globalSpillTime", &globalSpillTime);
    spillTree->SetBranchAddress("ustofSpillTime", &ustofSpillTime);
    // Loop over these entries
    // Only add them to the vectors if they are between
    // correct unix start and end times
    for (int t = 0; t < spillTree->GetEntries(); t++) {
      spillTree->GetEntry(t);

      if (globalSpillTime < startTime) continue;
      if (globalSpillTime > endTime) break;

      utofTimes.push_back(ustofSpillTime);
      dtofTimes.push_back(globalSpillTime);
    } // for (int t = 0; t < spillTree->GetEntries(); t++)
    dbFile->Close();
    delete dbFile;
    // Loop over spills then find the number of utof a dtof hits 
    int u = 0;
    for (int spill = 0; spill < tempUtofTimes.size(); spill++) {
      utofHits = 0;
      dtofHits = 0;
      double tempUtofSpill = tempUtofTimes.at(spill);
      double tempDtofSpill = tempDtofTimes.at(spill);
      while(u < tree->GetEntries()) {
	tree->GetEntry(u);
	u++;
      }

      hitSpillTree->Fill();
    }
  }

  fout->cd();
  hitSpillTree->Write();
  fout->Close();
}
