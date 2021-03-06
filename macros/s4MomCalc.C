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
  // Timing cuts for S3
  const double proLowS3 = -45;
  const double proHiS3  = 10.;
  const double piLowS3  = -63;
  const double piHiS3   = -57;
  // Timing cuts for S4
  const double piLowS4  = 80.;
  const double piHiS4   = 95.;
  const double proLowS4 = 106.;
  const double proHiS4  = 140.;
  // S1 -> S4 baseline length
  const double baselineS1S4 = 13.97;
  // S1 -> S4 baseline length
  const double baselineS1S3 = 10.9;
  // S3 -> S4 baseline length
  const double baselineS3S4 = 3.3;
  // Shift required to put MIP peak at 11ns
  const double s3s4Shift = 115. - 11.;
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
  // Window for matching hits between utof and dtof
  const double matchWindow = 10000.;

  THStack *hsMom = new THStack("hsMom", "Proton momentum as measured between S3 & S4; Proton momentum [GeV/c]; Events / spill");
  TLegend *leg = new TLegend(0.15, 0.5, 0.45, 0.88);
  
  TFile *fout = new TFile(Form("%s/s4MomCalcPlots.root", saveDir), "recreate");
  
  for (int nBlocks = 0; nBlocks < 4; nBlocks++) {
    TTree *outTree = new TTree(Form("outTree%dBlocks", nBlocks), "Spill tree");
    TH1D *s3s4TofAll = new TH1D(Form("s3s4Tof_%dblocks", nBlocks), Form("S3 - S4 ToF, %d blocks; S4 - S3 / ns; Events", nBlocks), 440, 90, 200);
    TH1D *s3s4TofAllMatched = new TH1D(Form("s3s4TofMatched_%dblocks", nBlocks), Form("S3 - S4 ToF, %d blocks; S4 - S3 / ns; Events", nBlocks), 350, 70, 200);
    TH1D *s3s4TofAlls1 = new TH1D(Form("s3s4TofMatcheds1_%dblocks", nBlocks), Form("S3 - S4 ToF, %d blocks; S4 - S3 / ns; Events", nBlocks), 400, 70, 200);
    TH1D *s3s4Mom = new TH1D(Form("s3s4Mom_%dblocks", nBlocks), Form("S3 - S4 proton momentum spectrum, %d blocks; Momentum [GeV/c]; Events", nBlocks), 100, 0.2, 0.8);
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

    vector<double> tempdstofHitsS4;
    vector<double> tempdstofHitsS2;
    vector<double> tempustofHitsS2;
    vector<double> tempustofHitsS2Drift;
    vector<double> tempustofHitsS1;
    vector<double> tempustofHitsS1Drift;
    vector<double> tempustofHitsS3;
    vector<double> temps3s4T;
    double tempDstofSpill;
    double tempUstofSpill;

    outTree->Branch("dstofHitsS4", &tempdstofHitsS4);
    outTree->Branch("dstofHitsS2", &tempdstofHitsS2);
    outTree->Branch("ustofHitsS2", &tempustofHitsS2);
    outTree->Branch("ustofHitsS2Drift", &tempustofHitsS2Drift);
    outTree->Branch("ustofHitsS1", &tempustofHitsS1);
    outTree->Branch("ustofHitsS1Drift", &tempustofHitsS1Drift);
    outTree->Branch("ustofHitsS3", &tempustofHitsS3);	
    outTree->Branch("s3s4T", &temps3s4T);
    outTree->Branch("DstofSpillTime", &tempDstofSpill);
    outTree->Branch("UstofSpillTime", &tempUstofSpill);

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
      TGraph *grTDiff = new TGraph();
      for (int db=0; db<spillTree->GetEntries(); db++) {
	spillTree->GetEntry(db);
	if (globalSpillTime < startTime) continue;
	if (globalSpillTime > endTime) break;
	dstofSpillTimes.push_back((globalSpillTime - dstofStart) * 1e9);
	ustofSpillTimes.push_back((ustofSpillTime - ustofStart) * 1e9);
	grTDiff->SetPoint(grTDiff->GetN(), grTDiff->GetN(), ustofSpillTime - globalSpillTime);
      } // for (int db=0; db<spillTree->GetEntries(); db++)
      cout<<ustofSpillTimes.size()<<" spills in run "<<irun<<endl;
      fout->cd();
      grTDiff->Write(Form("grTDiffRun%d_%dblocks", irun, nBlocks));

      // For each spill in the run make the vector of hits
      if (ustofSpillTimes.size() > 4) {
	for (int spill = 0; spill < ustofSpillTimes.size(); spill++) {
	  if (spill % 25 == 0) {
	    cout<<"Spill "<<spill<<" of "<<ustofSpillTimes.size()<<" in run "<<irun<<" of "<<runMax<<endl;
	  }

	  TH1D *s3s4Tof = new TH1D(Form("s3s4TofRun%dspill%d_%dblocks", irun, spill, nBlocks), Form("S3 - S4 ToF, run %d, spill, %d, %d blocks", irun, spill, nBlocks), 100, -200, 200);
	  TH1D *s3s4TofShift = new TH1D(Form("s3s4TofShiftRun%dspill%d_%dblocks", irun, spill, nBlocks), Form("S3 - S4 ToF, run %d, spill, %d, %d blocks", irun, spill, nBlocks), 100, -200, 200);
	  TGraph *grs3s4= new TGraph();
	  grs3s4->SetTitle(Form("S3 - S4 vs. S3, run %d, spill %d, %d blocks", irun, spill, nBlocks));
	  TGraph *grTmp = new TGraph();
	  if (spill < 4) {
	    for (int i = 0; i < 5; i++) {
	      grTmp->SetPoint(grTmp->GetN(),ustofSpillTimes[i],ustofSpillTimes[i]-dstofSpillTimes[i]);
	    }
	  }
	  else if (spill > ustofSpillTimes.size() - 3) {
	    for (int i = ustofSpillTimes.size() - 5; i < ustofSpillTimes.size(); i++) {
	      grTmp->SetPoint(grTmp->GetN(),ustofSpillTimes[i],ustofSpillTimes[i]-dstofSpillTimes[i]);
	    }
	  }
	  else {
	    for (int i = spill-2; i < spill+2; i++) {
	      grTmp->SetPoint(grTmp->GetN(),ustofSpillTimes[i],ustofSpillTimes[i]-dstofSpillTimes[i]);
	    }
	  }
	  // Fit straight line to the points
	  TF1 *f1 = new TF1("f1", "pol1");
	  grTmp->Fit(f1, "Q");
	  double grad = f1->GetParameter(1);

	  // Get dstof hits in the spill in vectors
	  vector<double> dstofHitsS4_1;
	  vector<double> dstofHitsS4_2;
	  vector<double> dstofHitsS2_1;
	  vector<double> dstofHitsS2_2;
	  vector<double> ustofHitsS2;
	  vector<double> ustofHitsS2Drift;
	  vector<double> ustofHitsS1;
	  vector<double> ustofHitsS1Drift;
	  vector<double> ustofHitsS3;
	  
	  vector<double> s3s4TofVec;

	  tempDstofSpill = dstofSpillTimes[spill]/1e9;
	  tempUstofSpill = ustofSpillTimes[spill]/1e9;

	  // These are for the matched hits only 
	  // After differences have all been calculated
	  vector<double> dstofHitsS4M;
	  vector<double> dstofHitsS2M;
	  vector<double> s3s4TM;
	  vector<double> ustofHitsS2M;
	  vector<double> ustofHitsS2DriftM;
	  vector<double> ustofHitsS1M;
	  vector<double> ustofHitsS1DriftM;
	  vector<double> ustofHitsS3M;

	  for (int d = 0; d < dstofTree1->GetEntries(); d++) {
	    dstofTree1->GetEntry(d);
	    if ((tofCoin1->fakeTimeNs[0] - dstofSpillTimes[spill]) < 0.) continue;
	    if ((tofCoin1->fakeTimeNs[0] - dstofSpillTimes[spill]) > 1e9) break;
	    double deltat = TMath::Abs(tofCoin1->fakeTimeNs[0]-tofCoin1->fakeTimeNs[1]);
	    double dstofHitT = min(tofCoin1->fakeTimeNs[0], tofCoin1->fakeTimeNs[1]) - (10. - TMath::Abs(deltat) / 2.);
	    // Only put things that are signal events in the file
	    if ((dstofHitT - tofCoin1->usTofSignal) < 200. &&
		(dstofHitT - tofCoin1->usTofSignal) > 70.) {
	      dstofHitsS2_1.push_back(tofCoin1->usTofSignal - dstofSpillTimes[spill]);
	      dstofHitsS4_1.push_back(dstofHitT - dstofSpillTimes[spill]);
	    }
	  } // for (int d = 0; d < dstofTree1->GetEntries(); d++)
	  for (int d = 0; d < dstofTree2->GetEntries(); d++) {
	    dstofTree2->GetEntry(d);
	    if ((tofCoin2->fakeTimeNs[0] - dstofSpillTimes[spill]) < 0.) continue;
	    if ((tofCoin2->fakeTimeNs[0] - dstofSpillTimes[spill]) > 1e9) break;
	    double deltat = TMath::Abs(tofCoin2->fakeTimeNs[0]-tofCoin2->fakeTimeNs[1]);
	    double dstofHitT = min(tofCoin2->fakeTimeNs[0], tofCoin2->fakeTimeNs[1]) - (10. - TMath::Abs(deltat) / 2.);
	    // Only put things that are signal events in the file
	    if ((dstofHitT - tofCoin2->usTofSignal) < 200. &&
		(dstofHitT - tofCoin2->usTofSignal) > 70.) {
	      dstofHitsS2_2.push_back(tofCoin2->usTofSignal - dstofSpillTimes[spill]);
	      dstofHitsS4_2.push_back(dstofHitT - dstofSpillTimes[spill]);
	    }
	  }
	  // Get ustof hits in a vector (accounting for time drift)
	  for (int u = 0; u < tree->GetEntries(); u++) {
	    tree->GetEntry(u);
	    if ((tS1 - ustofSpillTimes[spill]) < 0.) continue;
	    if ((tS1 - ustofSpillTimes[spill]) > 1e9) break;
	    // We only want the events with an S1-S2 coincidence
	    if (tTrig !=0) {
	      // Vector of unshifted hit times and also one of drifted one
	      for (int n = 0; n < nhit; n++) {
		ustofHitsS3.push_back(tToF[n] - ustofSpillTimes[spill]);
		ustofHitsS2.push_back(tTrig - ustofSpillTimes[spill]);
		ustofHitsS1.push_back(tS1 - ustofSpillTimes[spill]);
		ustofHitsS2Drift.push_back((tTrig - ustofSpillTimes[spill]) / (1. + grad));
		ustofHitsS1Drift.push_back((tS1-ustofSpillTimes[spill]) / (1. + grad));
	      } // for(int n = 0; n < nhit; n++) 
	    } // if (tTrig !=0)
	  } // for (int u = 0; u < tree->GetEntries(); u++)

	  //cout<<"Have counted "<<dstofHitsS4_1.size()+dstofHitsS4_2.size()<<" dstof hits and "<<ustofHitsS3.size()<<" ustof hits"<<endl;
	  // Now we want to try and match these hits together
	  // Match S2 hits in one filesystem with another
	  // Try and match each of these with the corresponding closest dtof S2 hit
	  int lastdh1 = 0;
	  int lastdh2 = 0;
	  for (int uh = 0; uh < ustofHitsS2Drift.size(); uh++) {
	    double diffLow = 999999999999.;
	    double diff = 0.;
	    int hitLow = -1;
	    int tdc = 0;
	    for (int dh1=lastdh1; dh1<dstofHitsS2_1.size(); dh1++) {
	      diff = ustofHitsS2Drift[uh] - dstofHitsS2_1[dh1];
	      if (abs(diff) < abs(diffLow)) {
		diffLow = diff;
		hitLow = dh1;
		tdc = 1;
		lastdh1 = dh1;
	      } // if (abs(diff) < abs(diffLow))
	    } // for (int dh=0; dh<dstofHits1.size(); dh++) 
	    for (int dh2=lastdh2; dh2<dstofHitsS2_2.size(); dh2++) {
	      diff = ustofHitsS2Drift[uh] - dstofHitsS2_2[dh2];
	      if (abs(diff) < abs(diffLow)) {
		diffLow = diff;
		hitLow = dh2;
		tdc = 2;
		lastdh2 = dh2;
	      } // if (abs(diff) < abs(diffLow))
	    } // for (int dh=0; dh<dstofHits2.size(); dh++)

	    ustofHitsS1M.push_back(ustofHitsS1[uh]);
	    ustofHitsS2M.push_back(ustofHitsS2[uh]);
	    ustofHitsS3M.push_back(ustofHitsS3[uh]);

	    if (tdc == 1) {
	      // Make sure we are only plotting events where particles are ID'd
	      // as the same kind in both TOFs
	      if (((ustofHitsS3[uh] - ustofHitsS2[uh]) < piHiS3 &&
		   (ustofHitsS3[uh] - ustofHitsS2[uh]) > piLowS3 &&
		   (dstofHitsS4_1[hitLow] - dstofHitsS2_1[hitLow]) > piLowS4 &&
		   (dstofHitsS4_1[hitLow] - dstofHitsS2_1[hitLow]) < piHiS4) ||
		  ((ustofHitsS3[uh] - ustofHitsS2[uh]) < proHiS3 &&
		   (ustofHitsS3[uh] - ustofHitsS2[uh]) > proLowS3 &&
		   (dstofHitsS4_1[hitLow] - dstofHitsS2_1[hitLow]) > proLowS4 &&
		   (dstofHitsS4_1[hitLow] - dstofHitsS2_1[hitLow]) < proHiS4))
		{
		  dstofHitsS2M.push_back(dstofHitsS2_1[hitLow]);
		  dstofHitsS4M.push_back(dstofHitsS4_1[hitLow]);
		  s3s4TofAll->Fill((ustofHitsS2[uh] - ustofHitsS3[uh]) - (dstofHitsS2_1[hitLow] - dstofHitsS4_1[hitLow]));
		  s3s4Tof->Fill((ustofHitsS2[uh] - ustofHitsS3[uh]) - (dstofHitsS2_1[hitLow] - dstofHitsS4_1[hitLow]));
		  s3s4TM.push_back((ustofHitsS2[uh] - ustofHitsS3[uh]) - (dstofHitsS2_1[hitLow] - dstofHitsS4_1[hitLow]));
		  if (abs(ustofHitsS2Drift[uh] - dstofHitsS2_1[hitLow]) < matchWindow) {
		    s3s4TofAllMatched->Fill((ustofHitsS2[uh] - ustofHitsS3[uh]) - (dstofHitsS2_1[hitLow] - dstofHitsS4_1[hitLow]));
		  } // if (abs(ustofHitsS2Drift[uh]-dstofHitsS2_1[hitLow]) < matchWindow)
		} // Is same particle in each TOF
	      
	      // Do this but for matched S1 utof hits 
	      if (((ustofHitsS3[uh] - ustofHitsS1[uh] - 34.) < piHiS3 &&
		   (ustofHitsS3[uh] - ustofHitsS1[uh] - 34.) > piLowS3 &&
		   (dstofHitsS4_1[hitLow] - dstofHitsS2_1[hitLow]) > piLowS4 &&
		   (dstofHitsS4_1[hitLow] - dstofHitsS2_1[hitLow]) < piHiS4) ||
		  ((ustofHitsS3[uh] - ustofHitsS1[uh] - 34.) < proHiS3 &&
		   (ustofHitsS3[uh] - ustofHitsS1[uh] - 34.) > proLowS3 &&
		   (dstofHitsS4_1[hitLow] - dstofHitsS2_1[hitLow]) > proLowS4 &&
		   (dstofHitsS4_1[hitLow] - dstofHitsS2_1[hitLow]) < proHiS4))
		{

		  if (abs(ustofHitsS1Drift[uh] - dstofHitsS2_1[hitLow]) < matchWindow) {
		    s3s4TofAlls1->Fill((ustofHitsS1[uh] - ustofHitsS3[uh]) - (dstofHitsS2_1[hitLow] - dstofHitsS4_1[hitLow]));
		    if (((ustofHitsS1[uh] - ustofHitsS3[uh]) - (dstofHitsS2_1[hitLow] - dstofHitsS4_1[hitLow])) > 120) 
		      {
			s3s4Mom->Fill(momFromTime(0.938, baselineS3S4, (ustofHitsS1[uh] - ustofHitsS3[uh]) - (dstofHitsS2_1[hitLow] - dstofHitsS4_1[hitLow]) - s3s4Shift));
		      }
		  } // if (abs(ustofHitsS1Drift[uh]-dstofHitsS2_1[hitLow]) < matchWindow)
		} // Is same particle in each TOF
	      s3s4TofVec.push_back((ustofHitsS2[uh] - ustofHitsS3[uh]) - (dstofHitsS2_1[hitLow] - dstofHitsS4_1[hitLow]));
	    }
	    else if (tdc == 2) {
	      if (((ustofHitsS3[uh] - ustofHitsS2[uh]) < piHiS3 &&
		   (ustofHitsS3[uh] - ustofHitsS2[uh]) > piLowS3 &&
		   (dstofHitsS4_2[hitLow] - dstofHitsS2_2[hitLow]) > piLowS4 &&
		   (dstofHitsS4_2[hitLow] - dstofHitsS2_2[hitLow]) < piHiS4) ||
		  ((ustofHitsS3[uh] - ustofHitsS2[uh]) < proHiS3 &&
		   (ustofHitsS3[uh] - ustofHitsS2[uh]) > proLowS3 &&
		   (dstofHitsS4_2[hitLow] - dstofHitsS2_2[hitLow]) > proLowS4 &&
		   (dstofHitsS4_2[hitLow] - dstofHitsS2_2[hitLow]) < proHiS4))
		{
		  dstofHitsS2M.push_back(dstofHitsS2_2[hitLow]);
		  dstofHitsS4M.push_back(dstofHitsS4_2[hitLow]);
		  s3s4TofAll->Fill((ustofHitsS2[uh] - ustofHitsS3[uh]) - (dstofHitsS2_2[hitLow] - dstofHitsS4_2[hitLow]));
		  s3s4Tof->Fill((ustofHitsS2[uh] - ustofHitsS3[uh]) - (dstofHitsS2_2[hitLow] - dstofHitsS4_2[hitLow]));
		  s3s4TM.push_back((ustofHitsS2[uh] - ustofHitsS3[uh]) - (dstofHitsS2_2[hitLow] - dstofHitsS4_2[hitLow]));
		  if (abs(ustofHitsS2Drift[uh] - dstofHitsS2_2[hitLow]) < matchWindow) {
		    s3s4TofAllMatched->Fill((ustofHitsS2[uh] - ustofHitsS3[uh]) - (dstofHitsS2_2[hitLow] - dstofHitsS4_2[hitLow]));
		    
		  } // if (abs(ustofHitsS2Drift[uh]-dstofHitsS2_1[hitLow]) < matchWindow)
		} // Is same particle in each TOF
	      
	      // Do this but for matched S1 utof hits 
	      if (((ustofHitsS3[uh] - ustofHitsS1[uh] - 34.) < piHiS3 &&
		   (ustofHitsS3[uh] - ustofHitsS1[uh] - 34.) > piLowS3 &&
		   (dstofHitsS4_2[hitLow] - dstofHitsS2_2[hitLow]) > piLowS4 &&
		   (dstofHitsS4_2[hitLow] - dstofHitsS2_2[hitLow]) < piHiS4) ||
		  ((ustofHitsS3[uh] - ustofHitsS1[uh] - 34.) < proHiS3 &&
		   (ustofHitsS3[uh] - ustofHitsS1[uh] - 34.) > proLowS3 &&
		   (dstofHitsS4_2[hitLow] - dstofHitsS2_2[hitLow]) > proLowS4 &&
		   (dstofHitsS4_2[hitLow] - dstofHitsS2_2[hitLow]) < proHiS4))
		{

		  if (abs(ustofHitsS1Drift[uh] + 34. - dstofHitsS2_2[hitLow]) < matchWindow) {
		    s3s4TofAlls1->Fill((ustofHitsS1[uh] - ustofHitsS3[uh]) - (dstofHitsS2_2[hitLow] - dstofHitsS4_2[hitLow]));
		    if (((ustofHitsS1[uh] - ustofHitsS3[uh]) - (dstofHitsS2_2[hitLow] - dstofHitsS4_2[hitLow])) > 118) 
		      {
			s3s4Mom->Fill(momFromTime(0.938, baselineS3S4, (ustofHitsS1[uh] - ustofHitsS3[uh]) - (dstofHitsS2_2[hitLow] - dstofHitsS4_2[hitLow]) - s3s4Shift));
		      }
		  } // if (abs(ustofHitsS1Drift[uh]-dstofHitsS2_2[hitLow]) < matchWindow)
		} // Is same particle in each TOF
	      s3s4TofVec.push_back((ustofHitsS2[uh] - ustofHitsS3[uh]) - (dstofHitsS2_2[hitLow] - dstofHitsS4_2[hitLow]));
	    }
	    
	    grs3s4->SetPoint(grs3s4->GetN(), diffLow, ustofHitsS2Drift[uh]);
	  } // for (int uh = 0; uh < ustofHits.size(); uh++)
	  fout->cd();
	  grs3s4->SetTitle(Form("Time difference vs. time since spill, Run %d, %d blocks; S3-S4 / ns; Time since spill start / ns", irun, nBlocks));
	  grs3s4->Write(Form("grs3s4Run%dspill%d_%dblocks", irun, spill, nBlocks)); 
     
	  s3s4Tof->Write();
	  tempdstofHitsS4 = dstofHitsS4M;
	  tempdstofHitsS2 = dstofHitsS2M;
	  tempustofHitsS2 = ustofHitsS2M;
	  tempustofHitsS3 = ustofHitsS3M;
	  temps3s4T = s3s4TM;
	  tempustofHitsS2Drift = ustofHitsS2Drift;
	  outTree->Fill();
	} // Loop over spills
	// Signals are gaussians
	TCanvas *c1 = new TCanvas(Form("c1%d",nBlocks));
	c1->SetLogy();
	s3s4TofAll->Draw("hist");
	c1->Print(Form("%s/s3s4Tof%dblocks.png", saveDir, nBlocks));
	c1->Print(Form("%s/s3s4Tof%dblocks.pdf", saveDir, nBlocks));

	TCanvas *c2 = new TCanvas(Form("c2%d",nBlocks));
	s3s4TofAllMatched->Draw("hist");
	c2->Print(Form("%s/s3s4TofMatched%dblocks.png", saveDir, nBlocks));
	c2->Print(Form("%s/s3s4TofMatched%dblocks.pdf", saveDir, nBlocks));

	TCanvas *c3 = new TCanvas("c3");
	TF1 *sPro = new TF1("sPro", "gaus", 115, 140);
	TF1 *sPi  = new TF1("sPi", "gaus", 105, 120);
	TF1 *fSignal = new TF1("signal", "gaus(0)+gaus(3)", 105, 140);
	fSignal->SetParNames("const 1", "mean 1", "sigma 1",
			     "const 2", "mean 2", "sigma 2");
	fSignal->SetLineColor(kBlack);
	s3s4TofAlls1->Fit(sPi, "R");
	s3s4TofAlls1->Fit(sPro, "R");
	Double_t parExp[6];
	sPro->GetParameters(&parExp[0]);
	sPi->GetParameters(&parExp[3]);
	fSignal->SetParameters(parExp);
	s3s4TofAlls1->Fit(fSignal, "R");
	s3s4TofAlls1->Draw("hist");
	fSignal->Draw("same");
	c3->Print(Form("%s/s3s4Tofs1Matched%dblocks.png", saveDir, nBlocks));
	c3->Print(Form("%s/s3s4Tofs1Matched%dblocks.pdf", saveDir, nBlocks));

	TCanvas *c4 = new TCanvas(Form("c4%d",nBlocks));
	s3s4Mom->Draw("hist");
	c4->Print(Form("%s/s3s4Mom%dblocks.png", saveDir, nBlocks));
	c4->Print(Form("%s/s3s4Mom%dblocks.pdf", saveDir, nBlocks));

	s3s4Mom->SetLineWidth(2);
	if (nBlocks==0) s3s4Mom->SetLineColor(kBlack);
	else if (nBlocks==1) s3s4Mom->SetLineColor(kRed);
	else if (nBlocks==2) s3s4Mom->SetLineColor(kBlue);
	else if (nBlocks==3) s3s4Mom->SetLineColor(kCyan+1);
	else if (nBlocks==4) s3s4Mom->SetLineColor(kOrange+1);

	hsMom->Add(s3s4Mom);


	s3s4Mom->Write();
	s3s4TofAll->Write();
	s3s4TofAllMatched->Write();
	s3s4TofAlls1->Write();
      } // if (ustofSpillTimes.size() > 10) 

      dbFile->Close();
      delete tofCoin1; 
      delete tofCoin2;
      dstofFile1->Close();
      dstofFile2->Close();

    } // for (int irun = runMin; irun < runMax+1; irun++) 

    ustofIn->Close();
    fout->cd();
    outTree->Write();
  } // for (int nBlocks = 0; nBlocks <= 4; nBlocks++) 
  fout->cd();
  hsMom->Write();
  fout->Close();
} // s4MomCalc
