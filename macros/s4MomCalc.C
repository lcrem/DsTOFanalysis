// s4MomCalc.C
// Calculates proton momentum spectrum in S4 by matching hits between S3 and S4

// Outputs momentum in GeV/c
double momFromTime(const double mass, const double baseline, const double time)
{
  double mom = 0.;
  mom = (mass / 3e8) * baseline * ( 1. / TMath::Sqrt( pow(time*1e-9, 2) - (pow((baseline / 9e16), 2)) ) );
  return mom;
}

void s4MomCalc(const char* saveDir,
	       const char* ustofDir="/zfs_home/sjones/mylinktoutof/",
	       const char* dstofDir="/scratch0/dbrailsf/temp/mylinktodtof/",
	       const char* spillDBDir="/zfs_home/sjones/spillDB/")
{
  gSystem->Load("libdstof.so");
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
  // 4 moderator blocks with -4cm bend
  const double start4Block = 1535836129;
  const double end4Block   = 1535879634;
  // Timing cuts
  const double piLow  = 80.;
  const double piHi   = 95.;
  const double proLow = 106.;
  const double proHi  = 140.;
  // S1 -> S4 baseline length
  const double baselineS1S4 = 13.97;
  const double baselineS1S3 = 10.9;
  // Define the runs to be used for varying number of blocks
  const char* str0Block = "Data_2018_8_31_b2_800MeV_0block.root";
  const char* str1Block = "Data_2018_9_1_b4_800MeV_1block_bend4cm.root";
  const char* str2Block = "Data_2018_9_1_b2_800MeV_2block_bend4cm.root";
  const char* str3Block = "Data_2018_9_1_b3_800MeV_3block_bend4cm.root";
  const char* str4Block = "Data_2018_9_1_b8_800MeV_4block_bend4cm.root";
  // ustof-dstof cable delay
  const double ustofDelay = 184.7;
  // Shift in ns required to to pion peak at speed of light
  const double dstofShift = 39.3;
  // S3 amplitude cut for protons
  // Apply to A1ToF and A2ToF
  // AK's standard cut
  const double ACut = 0.25;
  // SJ's bar by bar amplitude cut (by eye)
  // For A1 
  const std::vector<double> A1CutVec = {0.25, 0.25, 0.275, 0.2, 0.25, 0.25, 0.25, 0.275, 0.325, 
					0.3, 0.3, 0.2, 0.2, 0.25, 0.225, 0.25, 0.25, 0.3, 0.3, 
					0.25, 0.3, 0.3};
  // For A2
  const std::vector<double> A2CutVec = {0.3, 0.275, 0.25, 0.125, 0.3, 0.3, 0.15, 0.225, 0.25, 0.25,
				       0.225, 0.225, 0.2, 0.225, 0.225, 0.225, 0.2, 0.275, 0.25, 
				       0.225, 0.3, 0.3};

  THStack *hsMom = new THStack("hsMom", "Proton momentum as measured in S4; Proton momentum [GeV/c]; Events / spill");
  TLegend *leg = new TLegend(0.15, 0.5, 0.45, 0.88);
  
  TFile *fout = new TFile(Form("%s/s4MomCalcPlots_eff.root", saveDir), "recreate");

  for (int nBlocks = 0; nBlocks <= 4; nBlocks++) {
    cout<<nBlocks<<" block case"<<endl;
    // Find the correct dstof files
    Int_t runMin=-1;
    Int_t runMax=-1;
    double startTime = 0;
    double endTime   = 0;

    const char* nustof;
    if (nBlocks==0) {
      nustof = Form("%sData_2018_8_31_b2_800MeV_0block.root", ustofDir);
      startTime = start0Block;
      endTime   = end0Block;
    }
    else if (nBlocks==1) {
      nustof = Form("%sData_2018_9_1_b4_800MeV_1block_bend4cm.root", ustofDir);
      startTime = start1Block;
      endTime   = end1Block;
    }
    else if (nBlocks==2) {
      nustof = Form("%sData_2018_9_1_b2_800MeV_2block_bend4cm.root", ustofDir);
      startTime = start2Block;
      endTime   = end2Block;
    }
    else if (nBlocks==3) {
      nustof = Form("%sData_2018_9_1_b3_800MeV_3block_bend4cm.root", ustofDir);
      startTime = start3Block;
      endTime   = end3Block;
    }
    else if (nBlocks==4) {
      nustof = Form("%sData_2018_9_1_b8_800MeV_4block_bend4cm.root", ustofDir);
      startTime = start4Block;
      endTime   = end4Block;
    }

    // Open ustof file
    TFile *ustofIn = new TFile(nustof, "read");
    double tToF[50];
    float xToF[50];
    float yToF[50];
    float A1ToF[50];
    float A2ToF[50];
    double tTrig;
    double tS1;
    double tSoSd;
    int nhit;
    int nBar[50];

    TTree *tree = (TTree*)ustofIn->Get("tree");
    tree->SetBranchAddress("xToF", xToF);
    tree->SetBranchAddress("yToF", yToF);
    tree->SetBranchAddress("A1ToF", A1ToF);
    tree->SetBranchAddress("A2ToF", A2ToF);
    tree->SetBranchAddress("nhit", &nhit);
    tree->SetBranchAddress("tS1", &tS1);
    tree->SetBranchAddress("tToF", tToF);
    tree->SetBranchAddress("tTrig", &tTrig);
    tree->SetBranchAddress("tSoSd", &tSoSd);
    tree->SetBranchAddress("nBar", nBar);
    // Get unix start and end times for the ustof files
    TNamed *start = 0;
    TNamed *end = 0;
    ustofIn->GetObject("start_of_run", start);
    ustofIn->GetObject("end_of_run", end);

    const char* startchar = start->GetTitle();
    std::string startstr(startchar);
    std::string unixstart = startstr.substr(25,10);
    const int ustofStart = stoi(unixstart);
  
    const char* endchar = end->GetTitle();
    std::string endstr(endchar);
    std::string unixend = endstr.substr(23,10);
    const int ustofEnd = stoi(unixend);

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
    } // for (int irun=950; irun<1400; irun++)
    
    cout << "Min and max runs are " << runMin << " " << runMax << endl;

    // Use the spill DB files to match the spills
    for (int irun = runMin; irun < runMax+1; irun++) {
      // Open spill DB files
      TFile *dbFile = new TFile(Form("%s/spillDB_run%d_run%d.root", spillDBDir, irun, irun), "read");
      double globalSpillTime;
      double ustofSpillTime;
      TTree *spillTree = (TTree*)dbFile->Get("spillTree");
      spillTree->SetBranchAddress("globalSpillTime", &globalSpillTime);
      spillTree->SetBranchAddress("ustofSpillTime", &ustofSpillTime);
      // Open dstof files
      TFile *dstofFile1 = new TFile(Form("%s/run%d/DsTOFcoincidenceRun%d_tdc1.root", dstofDir, irun, irun), "read");
      TFile *dstofFile2 = new TFile(Form("%s/run%d/DsTOFcoincidenceRun%d_tdc2.root", dstofDir, irun, irun), "read");
      RawDsTofCoincidence *tofCoin1 = NULL;
      RawDsTofCoincidence *tofCoin2 = NULL;
      TTree* dstofTree1 = (TTree*)dstofFile1->Get("tofCoinTree");
      TTree* dstofTree2 = (TTree*)dstofFile2->Get("tofCoinTree");
      dstofTree1->SetBranchAddress("tofCoin", &tofCoin1);
      dstofTree2->SetBranchAddress("tofCoin", &tofCoin2);
      dstofTree1->GetEntry(0);
      UInt_t dstofStart = tofCoin1->unixTime[0];
      // Make vectors of spill times
      vector<double> ustofSpillTimes;
      vector<double> dstofSpillTimes;
      for (int db=0; db<spillTree->GetEntries(); db++) {
	spillTree->GetEntry(db);
	if (globalSpillTime < startTime) continue;
	if (globalSpillTime > endTime) break;
	dstofSpillTimes.push_back(globalSpillTime - dstofStart);
	ustofSpillTimes.push_back(ustofSpillTime - ustofStart);
      } // for (int db=0; db<spillTree->GetEntries(); db++)
      cout<<ustofSpillTimes.size()<<" spills in run "<<irun<<endl;
      // For each matched spill make a vector of all the hit times in each system
      //for (int spill = 0; spill < utofSpillTimes.size(); spill++) {
	// Calculate time drift for this spill by taking the 5 spills either side of it
      //} // for (int spill = 0; spill < utofSpillTimes.size(); spill++)

      // For second spill in the run make the vector of hits
      if (ustofSpillTimes.size() > 10) {
      TGraph *grTmp = new TGraph();
	for (int i=0; i<10; i++) {
	  grTmp->SetPoint(grTmp->GetN(), ustofSpillTimes[i], ustofSpillTimes[i] - dstofSpillTimes[i]);
	} // for (int i=0; i<ustofSpillTimes; i++) 
	// Fit straight line to the points
	TF1 *f1 = new TF1("f1", "pol1");
	grTmp->Fit(f1);
	double grad = f1->GetParameter(1);
	cout<<"Time drift is "<<grad<<" seconds / seconds"<<endl;

	// Get dstof hits in the spill in vectors
	vector<double> dstofHits1;
	vector<double> dstofHits2;
	vector<double> ustofHits;
	cout<<"Dtof spill time is "<<dstofSpillTimes[1]<<endl;
	
	for (int d = 0; d < dstofTree1->GetEntries(); d++) {
	  dstofTree1->GetEntry(d);
	  if ((tofCoin1->fakeTimeNs[0]/1e9 - dstofSpillTimes[1]) < 0.) continue;
	  if ((tofCoin1->fakeTimeNs[0]/1e9 - dstofSpillTimes[1]) > 1.) break;
	  double deltat = TMath::Abs(tofCoin1->fakeTimeNs[0]-tofCoin1->fakeTimeNs[1]);
	  double dstofHitT = min(tofCoin1->fakeTimeNs[0], tofCoin1->fakeTimeNs[1]) - (10. - TMath::Abs(deltat) / 2.);
	  dstofHits1.push_back(dstofHitT/1e9 - dstofSpillTimes[1]);
	} // for (int d = 0; d < dstofTree1->GetEntries(); d++)
	for (int d = 0; d < dstofTree2->GetEntries(); d++) {
	  dstofTree2->GetEntry(d);
	  if ((tofCoin2->fakeTimeNs[0]/1e9 - dstofSpillTimes[1]) < 0.) continue;
	  if ((tofCoin2->fakeTimeNs[0]/1e9 - dstofSpillTimes[1]) > 1.) break;
	  double deltat = TMath::Abs(tofCoin2->fakeTimeNs[0]-tofCoin2->fakeTimeNs[1]);
	  double dstofHitT = min(tofCoin2->fakeTimeNs[0], tofCoin2->fakeTimeNs[1]) - (10. - TMath::Abs(deltat) / 2.);
	  dstofHits2.push_back(dstofHitT/1e9 - dstofSpillTimes[1]);
	}
	// Get ustof hits in a vector (accounting for time drift)
	for (int u = 0; u < tree->GetEntries(); u++) {
	  tree->GetEntry(u);
	  if ((tS1/1e9 - ustofSpillTimes[1]) < 0.) continue;
	  if ((tS1/1e9 - ustofSpillTimes[1]) > 1.) break;
	  for (int n = 0; n < nhit; n++) {
	    ustofHits.push_back((tToF[n]/1e9 - ustofSpillTimes[1]) / (1. - grad));
	  } // for(int n = 0; n < nhit; n++) 
	} // for (int u = 0; u < tree->GetEntries(); u++)

	  cout<<"Have counted "<<dstofHits1.size()+dstofHits2.size()<<" dstof hits and "<<ustofHits.size()<<" ustof hits"<<endl;
      } // if (ustofSpillTimes.size() > 10) 


      dbFile->Close();
      delete tofCoin1; 
      delete tofCoin2;
      dstofFile1->Close();
      dstofFile2->Close();

    } // for (int irun = runMin; irun < runMax+1; irun++) 

    ustofIn->Close();
  } // for (int nBlocks = 0; nBlocks <= 4; nBlocks++) 

} // s4MomCalc
