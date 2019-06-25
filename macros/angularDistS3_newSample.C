// angularDistS3.C
// Outputs momentum in GeV/c
double momFromTime(const double mass, const double baseline, const double time)
{
  double mom = 0.;
  mom = (mass / 3e8) * baseline * ( 1. / TMath::Sqrt( pow(time*1e-9, 2) - (pow((baseline / 9e16), 2)) ) );
  return mom;
}

void angularDistS3_newSample(const char* saveDir,
		   const char* ustofDir="/nfs/scratch0/dbrailsf/data_backup/utof_backup_firsthitpinnedtounixtime/Data_root_v3_wo_walk_corr/",
		   const char* dstofDir="/nfs/scratch0/dbrailsf/data_backup/dtof_backup/",
		   const char* spillDir="/scratch0/sjones/spillDB/") 
{
  gROOT->SetBatch(kTRUE);
  // Edges of S3 in beam coordinate system
  const double s3StartX = -0.5168;
  const double s3EndX   = 0.9970;
  const double s3s1StartY = 9.0569 + 1.77;
  const double s3s1EndY   = 8.9146 + 1.77;
  // z coordinates
  const double s3BarTop    = 62.; 
  const double s3BarBottom = -60.;
  // Define the runs to be used for varying number of blocks
  const char* str0Block = "Data_2018_8_31_b2_800MeV_0block.root";
  const char* str1Block = "Data_2018_9_1_b4_800MeV_1block_bend4cm.root";
  const char* str2Block = "Data_2018_9_1_b2_800MeV_2block_bend4cm.root";
  const char* str3Block = "Data_2018_9_1_b3_800MeV_3block_bend4cm.root";
  // New datasets to combine for the 4 block data
  const char* str4Block0 = "Data_2018_8_28_b5.root";
  const char* str4Block1 = "Data_2018_9_3_b2_800MeV_4block_bend4cm.root";
  const char* str4Block2 = "Data_2018_8_29_b4.root";
  const char* str4Block3 = "Data_2018_8_29_b1.root";
  std::vector<const char*> str4BlockVec = {str4Block0, str4Block1, str4Block2, str4Block3};
  // Deadtime corrections
  // Just use a constant ratio for the 0 block case
  const double block0Slope = -0.0003738;
  const double block0Const = 0.2332;
  const double block1Slope = -0.0002569;
  const double block1Const =  0.3965;
  const double block2Slope = -0.0001739;
  const double block2Const =  0.419;
  const double block3Slope = -0.0001498;
  const double block3Const =  0.4343;
  // 4 block data
  const double block4Slope0 = -0.0005476;
  const double block4Const0 = 0.983;
  const double block4Slope1 = -0.0001622;
  const double block4Const1 = 0.5246;
  const double block4Slope2 = -0.0003664;
  const double block4Const2 = 0.7991;
  const double block4Slope3 = -0.0005316;
  const double block4Const3 = 1.018;
  std::vector<double> block4SlopeVec = {block4Slope0, block4Slope1, 
					block4Slope2, block4Slope3};
  std::vector<double> block4ConstVec = {block4Const0, block4Const1, 
					block4Const2, block4Const3};
  // Unix timestamps for variable block moves
  // 0.8GeV/c, 0 blocks
  // 31/08/2018
  const double start0Block = 1535713289;
  const double end0Block   = 1535716132;
  // 0.8GeV/c, 1 block
  // 01/09/2018
  const double start1Block = 1535796057;
  const double end1Block   = 1535799112;
  // 0.8GeV/c, 2 blocks
  // 01/09/2018
  const double start2Block = 1535789157;
  const double end2Block   = 1535792026;
  // 0.8GeV/c, 3 block
  // 01/09/2018
  const double start3Block = 1535792404;
  const double end3Block   = 1535795300;
  // const double end3Block   = 1535798437;
  // 0.8GeV/c, 4 block
  // 4 moderator blocks with -4cm bend
  const double start4Block = 1535836129;
  const double end4Block   = 1535879634;
  // Time of flight cuts for S1 to S3
  // Particles travelling at c should cross the distance in between about 36.9ns
  // and 35.6ns
  const double tLight = 35.8; // Expected time in ns of light speed particles according to AK
  // This is before the shift is applied
  const double proLow  = 54.5;
  const double proHi   = 115.;
  const double piLow = 35.75;
  const double piHi  = 37.75;
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
  // Time cut for hits interacting in neighbouring S3 bars
  const double twoBarCut = 0.5; // ns

  THStack *hsThetaS1pro   = new THStack("hsThetaS1pro", "S1 #cap S3 angular distribution of proton hits; #theta / degrees; Events / spill / degree");
  THStack *hsThetaS1S2pro = new THStack("hsThetaS1S2pro", "S1 #cap S2 #cap S3 angular distribution of proton hits; #theta / degrees; Events / spill / degree");
  THStack *hsPhiS1pro     = new THStack("hsPhiS1pro", "S1 #cap S3 angular distribution of proton hits; #phi / degrees; Events / spill / degree");
  THStack *hsPhiS1S2pro   = new THStack("hsPhiS1S2pro", "S1 #cap S2 #cap S3 angular distribution of proton hits; #phi / degrees; Events / spill / degree");

  THStack *hsThetaS1pi   = new THStack("hsThetaS1pi", "S1 #cap S3 angular distribution of MIP hits; #theta / degrees; Events / spill / degree");
  THStack *hsThetaS1S2pi = new THStack("hsThetaS1S2pi", "S1 #cap S2 #cap S3 angular distribution of MIP hits; #theta / degrees; Events / spill / degree");
  THStack *hsPhiS1pi     = new THStack("hsPhiS1pi", "S1 #cap S3 angular distribution of MIP hits; #phi / degrees; Events / spill / degree");
  THStack *hsPhiS1S2pi   = new THStack("hsPhiS1S2pi", "S1 #cap S2 #cap S3 angular distribution of MIP hits; #phi / degrees; Events / spill / degree");

  THStack *hsPhiS1S2ratio   = new THStack("hsPhiS1S2ratio", "S1 #cap S2 #cap S3 angular distribution of proton/MIP ratio; #phi / degrees;  Protons/MIPs");
  THStack *hsThetaS1S2ratio = new THStack("hsThetaS1S2ratio", "S1 #cap S2 #cap S3 angular distribution of proton/MIP ratio; #theta / degrees;  Protons/MIPs");
  THStack *hsPhiS1ratio   = new THStack("hsPhiS1ratio", "S1 #cap S3 angular distribution of proton/MIP ratio; #phi / degrees;  Protons/MIPs");
  THStack *hsThetaS1ratio = new THStack("hsThetaS1ratio", "S1 #cap S3 angular distribution of proton/MIP ratio; #theta / degrees;  Protons/MIPs");

  THStack *hsMomS1S2 = new THStack("hsMomS1S2", "Proton momentum measured in S3 (S1 & S2 trigger); Proton momentum [GeV/c]; Events / spill");
  THStack *hsMomS1 = new THStack("hsMomS1", "Proton momentum measured in S3 (S1 trigger only); Proton momentum [GeV/c]; Events / spill");

  THStack *hsutof1dS1S2 = new THStack("hsutof1dS1S2", "Time of flight as measured in S3 (S1 & S2 trigger); Time of flight / ns; Events / spill");
  THStack *hsutof1dS1   = new THStack("hsutof1dS1", "Time of flight as measured in S3 (S1 trigger only); Time of flight / ns; Events / spill");

  TFile *fout = new TFile(Form("%s/angularDistS3plots.root", saveDir), "recreate");

  TLegend *leg = new TLegend(0.15, 0.55, 0.3, 0.85);
  TLegend *legTheta = new TLegend(0.23, 0.5, 0.38, 0.85);
  TLegend *legTof = new TLegend(0.71, 0.53, 0.88, 0.85);
  TLegend *legPhiRatio = new TLegend(0.38, 0.56, 0.55, 0.88);

  TH2D *hMom2D_0blkQ = new TH2D("hMom2D_0blkQ", "Quick peak 0 block; x / cm; Bar", 100, -10, 170, 22, 0.5, 22.5);
  TH2D *hMom2D_0blkS = new TH2D("hMom2D_0blkS", "Slow peak 0 block; x / cm; Bar", 100, -10, 170, 22, 0.5, 22.5);
  TH2D *hMom2D_1blkQ = new TH2D("hMom2D_1blkQ", "Quick peak 1 block; x / cm; Bar", 100, -10, 170, 22, 0.5, 22.5);
  TH2D *hMom2D_1blkS = new TH2D("hMom2D_1blkS", "Slow peak 1 block; x / cm; Bar", 100, -10, 170, 22, 0.5, 22.5);
  TH2D *hMom2D_2blkQ = new TH2D("hMom2D_2blkQ", "Quick peak 2 block; x / cm; Bar", 100, -10, 170, 22, 0.5, 22.5);
  TH2D *hMom2D_2blkS = new TH2D("hMom2D_2blkS", "Slow peak 2 block; x / cm; Bar", 100, -10, 170, 22, 0.5, 22.5);

  double binsTheta[] = {-3.8, -3.7, -3.6, -3.5, -3.4, -3.3, -3.2, -3.1, -3.,
			-2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, -2., 
			-1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.,
			-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 
			0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
			1.0, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 
			2.0, 2.2, 2.4, 2.6, 2.8, 
			3.0, 3.2, 3.4, 3.6, 3.8,  
			4.0, 4.25, 4.5, 4.75, 
			5.0, 5.25, 5.8,
		        6.2};
  int binnum = sizeof(binsTheta)/sizeof(double) - 1;

  for (int nBlocks = 0; nBlocks < 5; nBlocks++) {
    int nSpills = 0;
    TH1D *hThetaS1pro   = new TH1D(Form("hThetaS1pro%d", nBlocks), Form("Angular distribution of proton hits in S3 (S1 trigger only), %d blocks; #theta / degrees; Events / spill / degree", nBlocks), binnum, binsTheta);
    hThetaS1pro->Sumw2();
    TH1D *hThetaS1S2pro = new TH1D(Form("hThetaS1S2pro%d", nBlocks), Form("Angular distribution of proton hits in S3 (S1 & S2 triggers), %d blocks; #theta / degrees; Events / spill", nBlocks), binnum, binsTheta);
    hThetaS1S2pro->Sumw2();
    TH1D *hPhiS1pro   = new TH1D(Form("hPhiS1pro%d", nBlocks), Form("Angular distribution of proton hits in S3 (S1 trigger only), %d blocks; #phi / degrees; Events / spill / degree", nBlocks), 22, -3.15, 3.25);
    hPhiS1pro->Sumw2();
    TH1D *hPhiS1S2pro = new TH1D(Form("hPhiS1S2pro%d", nBlocks), Form("Angular distribution of proton hits in S3 (S1 & S2 triggers), %d blocks; #phi / degrees; Events / spill", nBlocks), 22, -3.15, 3.25);
    hPhiS1S2pro->Sumw2();

    TH1D *hThetaS1pi   = new TH1D(Form("hThetaS1pi%d", nBlocks), Form("Angular distribution of pion hits in S3 (S1 trigger only), %d blocks; #theta / degrees; Events / spill / degree", nBlocks), binnum, binsTheta);
    hThetaS1pi->Sumw2();
    TH1D *hThetaS1S2pi = new TH1D(Form("hThetaS1S2pi%d", nBlocks), Form("Angular distribution of pion hits in S3 (S1 & S2 triggers), %d blocks; #theta / degrees; Events / spill", nBlocks), binnum, binsTheta);
    hThetaS1S2pi->Sumw2();
    TH1D *hPhiS1pi   = new TH1D(Form("hPhiS1pi%d", nBlocks), Form("Angular distribution of pion hits in S3 (S1 trigger only), %d blocks; #phi / degrees; Events / spill / degree", nBlocks), 22, -3.15, 3.25);
    hPhiS1pi->Sumw2();
    TH1D *hPhiS1S2pi = new TH1D(Form("hPhiS1S2pi%d", nBlocks), Form("Angular distribution of pion hits in S3 (S1 & S2 triggers), %d blocks; #phi / degrees; Events / spill", nBlocks), 22, -3.15, 3.25);
    hPhiS1S2pi->Sumw2();

    TH1D *hThetaS1ratio = new TH1D(Form("hThetaS1ratio%d", nBlocks), Form("S1 #cap S3 angular distribution of proton/MIP ratio, %d blocks; #phi / degrees; Protons/MIPs", nBlocks), binnum, binsTheta);
    TH1D *hPhiS1ratio   = new TH1D(Form("hPhiS1ratio%d", nBlocks), Form("S1 #cap S3 angular distribution of proton/MIP, %d blocks; #phi / degrees; Protons/MIPs", nBlocks), 22, -3.15, 3.25);
    TH1D *hThetaS1S2ratio = new TH1D(Form("hThetaS1S2ratio%d", nBlocks), Form("S1 #cap S2 #cap S3 angular distribution of proton/MIP ratio, %d blocks; #phi / degrees; Protons/MIPs", nBlocks), binnum, binsTheta);
    TH1D *hPhiS1S2ratio   = new TH1D(Form("hPhiS1S2ratio%d", nBlocks), Form("S1 #cap S2 #cap S3 angular distribution of proton/MIP ratio, %d blocks; #phi / degrees; Protons/MIPs", nBlocks), 22, -3.15, 3.25);

    TH1D *hutof1dS1 = new TH1D(Form("hutof1dS1_%d",nBlocks), Form("Time of flight, %d blocks (S1 trigger only); S3 - S1 / ns; Events / spill", nBlocks), 250, 25, 125);
    hutof1dS1->Sumw2();
    TH1D *hutof1dS1S2 = new TH1D(Form("hutof1dS1S2_%d",nBlocks), Form("Time of flight, %d blocks (S1 & S2 trigger); S3 - S1 / ns; Events / spill", nBlocks), 250, 25, 125);
    hutof1dS1S2->Sumw2();

    TH1D *hMomS1S2 = new TH1D(Form("hMomS1S2_%d",nBlocks), Form("Proton momentum measured in S3, %d blocks; Proton momentum [GeV/c]; Events / spill", nBlocks), 100, 0.3, 0.7);
    hMomS1S2->Sumw2();
    TH1D *hMomS1 = new TH1D(Form("hMomS1_%d",nBlocks), Form("Proton momentum measured in S3, %d blocks; Proton momentum [GeV/c]; Events / spill", nBlocks), 100, 0.3, 0.7);
    hMomS1->Sumw2();

    // XY distributions in S3
    TH2D *hprotonXY = new TH2D(Form("hprotonXY%d", nBlocks), Form("S3 spatial distribution of proton hits, %d blocks; x / cm; y / cm; Events / spill", nBlocks), 100, -10., 170., 22, 0., 120.);
    TH2D *hpionXY   = new TH2D(Form("hpionXY%d", nBlocks), Form("S3 spatial distribution of MIP hits, %d blocks; x / cm; y / cm; Events / spill", nBlocks), 100, -10., 170., 22, 0., 120.);

    if (nBlocks != 4) {
      // Number of protons and number of MIPs
      int nP  = 0;
      int nPi = 0;
      // Define signal and background functions to be fitted
      TF1 *sPro = new TF1("sPro", "gaus", proLow, proHi);
      TF1 *sPi  = new TF1("sPi", "gaus", piLow, piHi);

      // Find the correct dstof files
      Int_t runMin=-1;
      Int_t runMax=-1;

      double startTime = 0;
      double endTime   = 0;

      const char* nustof;
      double slope;
      double constant;
      if (nBlocks==0) {
	nustof = str0Block;
	startTime = start0Block;
	endTime   = end0Block;
	slope = block0Slope;
	constant = block0Const;
      }
      else if (nBlocks==1) {
	nustof = str1Block;
	startTime = start1Block;
	endTime   = end1Block;
	slope = block1Slope;
	constant = block1Const;
      }
      else if (nBlocks==2) {
	nustof = str2Block;
	startTime = start2Block;
	endTime   = end2Block;
	slope = block2Slope;
	constant = block2Const;
      }
      else if (nBlocks==3) {
	nustof = str3Block;
	startTime = start3Block;
	endTime   = end3Block;
	slope    = block3Slope;
	constant = block3Const;
      }

      // Find dtof runs
      for (int irun=950; irun<1400; irun++) {
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
    
      cout << "Min and max dtof runs are " << runMin << " " << runMax << endl;

      std::vector<double> dtofTimes;
      std::vector<int> dtofS1S2Hits;
      std::vector<double> utofTimes;
      std::vector<int> utofS1S2Hits;

      // Open the appropriate spill DB files and get the spill times
      for (int irun = runMin; irun < runMax+1; irun++) {
	TFile *dbFile = new TFile(Form("%s/spillDB_run%d_run%d.root", spillDir, irun, irun), "read");
	TTree *spillTree = (TTree*)dbFile->Get("spillTree");
	double globalSpillTime;
	double ustofSpillTime;
	spillTree->SetBranchAddress("globalSpillTime", &globalSpillTime);
	spillTree->SetBranchAddress("ustofSpillTime", &ustofSpillTime);
	for (int t = 0; t < spillTree->GetEntries(); t++) {
	  spillTree->GetEntry(t);
	  if (globalSpillTime >= startTime && globalSpillTime <= endTime) {
	    dtofTimes.push_back(globalSpillTime);
	    utofTimes.push_back(ustofSpillTime);
	  }
	} // for (int t = 0; t < spillTree->GetEntries(); t++)

	dbFile->Close();
	delete dbFile;
      } // for (int irun = runMin; irun < runMax+1; irun++) 

      dtofS1S2Hits.resize(dtofTimes.size(), 0);
      utofS1S2Hits.resize(dtofTimes.size(), 0);

      cout<<"Finding number of dtof hits in each spill"<<endl;
      int lastt = 0;
      int lastrun = 0;
      for (int s=0; s<dtofTimes.size(); s++) {
	if (s % 25 == 0) {
	  cout<<"Spill "<<s<<" of "<<dtofTimes.size()<<endl;
	}
	// Loop over the all the files
	for (int irun = runMin; irun < runMax+1; irun++) {
	
	  TFile *dtofFile = new TFile(Form("%srun%d/DsTOFtreeRun%d_tdc1.root", dstofDir, irun, irun), "read");
	  RawDsTofHeader *tof = NULL;
	  TTree *tofTree = (TTree*)dtofFile->Get("tofTree");
	  tofTree->SetBranchAddress("tof", &tof);
	  tofTree->GetEntry(0);
	  double firstTime = tof->unixTime;
	  tofTree->GetEntry(tofTree->GetEntries()-1);
	  double lastTime = tof->unixTime;
	  double lastS1S2Dtof = 0.;
	  // Spill is in this file
	  if (firstTime <= dtofTimes[s] && lastTime >= dtofTimes[s]) {
	    if (irun != lastrun) {
	      lastt = 0;
	      lastrun = irun;
	    }
	    // Loop over all entries and count the number of S1 S2 hits within the spill
	    for (int t=lastt; t<tofTree->GetEntries(); t++) {
	      tofTree->GetEntry(t);
	      if ((tof->fakeTimeNs/1e9)+firstTime < dtofTimes[s]) continue;
	      if ((tof->fakeTimeNs/1e9)+firstTime > dtofTimes[s]+1.) break;
	      // Is within a spill
	      if ((tof->fakeTimeNs/1e9)+firstTime >= dtofTimes[s] &&
		  (tof->fakeTimeNs/1e9)+firstTime <= dtofTimes[s]+1. && 
		  tof->channel == 13 && (tof->fakeTimeNs - lastS1S2Dtof) > 500.) {
		dtofS1S2Hits[s]++;
		lastt = t;
		lastS1S2Dtof = tof->fakeTimeNs;
	      } // Is within the spill
	    } // for (int t=0; t<tofTree->GetEntries(); t++)
	  } // if (firstTime <= dtofTimes[s] && lastTime >=  dtofTimes[s])

	  delete tof;
	  dtofFile->Close();
	  delete dtofFile;
	} // for (int irun = runMin; irun < runMax+1; irun++) 
      } // for (int s=0; s<dtofTimes.size(); s++)

      TFile *futof = new TFile(Form("%s/%s",ustofDir,nustof), "read");

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

      TTree *tree = (TTree*)futof->Get("tree");

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

      double lastSpill = 0.; 

      TNamed *start = 0;
      TNamed *end   = 0;
      futof->GetObject("start_of_run", start);
      futof->GetObject("end_of_run", end);
    
      const char* startchar = start->GetTitle();
      std::string startstr(startchar);
      std::string unixstart = startstr.substr(25,10);
      int startTimeUtof = stoi(unixstart);

      int lastut = 0;
      // Loop over the spills and perform the adjustment for each spill
      for (int s = 0; s < utofTimes.size(); s++) {
	if (s % 25 == 0) cout<<"Getting hits from spill "<<s<<" of "<<utofTimes.size()<<endl;
	double deadtimeWeight = dtofS1S2Hits[s] * slope + constant;
	// Initial data quality loop
	bool isGood = false;
	for (int t=lastut; t<tree->GetEntries(); t++) {
	  tree->GetEntry(t);
	  if ((tS1/1e9) + startTimeUtof < utofTimes[s]) continue;
	  if ((tS1/1e9) + startTimeUtof > utofTimes[s] + 1.) break;

	  if ((tTrig/1e9) + startTimeUtof > utofTimes[s] + 0.47) {
	    isGood = true;
	    nSpills++;
	    break;
	  }
	} // Initial data quality loop
	// Only do this spill if it's good
	if (isGood) {
	  for (int t=lastut; t<tree->GetEntries(); t++) {
	    tree->GetEntry(t);
	    if ((tS1/1e9) + startTimeUtof < utofTimes[s]) continue;
	    if ((tS1/1e9) + startTimeUtof > utofTimes[s] + 1.) break;
	    // Count number of spills
	    // if (tSoSd != lastSpill && tSoSd -lastSpill > 1e9) {
	    // nSpills++;
	    // lastSpill = tSoSd;
	    // } // if (tSoSd != lastSpill && tSoSd -lastSpill > 1e9) 

	    for (int nh=0; nh<nhit; nh++) {
	      // Checks if the next S3 hit is in a neighbouring bar
	      // Only want to count one of them if it is
	      // Check from this hit onwards to see if there are any double hit candidates
	      bool isDouble = false;
	      if (nh < nhit-1) {
		for (int nh2=nh+1; nh2<nhit; nh2++) {
		  // Timing cut and checks if in the neighbouring bar
		  if (abs(tToF[nh] - tToF[nh2]) < twoBarCut && 
		      (nBar[nh] == nBar[nh2]-1 || nBar[nh] == nBar[nh2]+1)) {
		    isDouble = true;
		    break;
		  }
		} // for (int nh2=nh; nh2<nhit; nh2++)
	      } // if (nh < nh-1)
	      // If it's not a double hit then go ahead
	      if (!isDouble) {
		double tofCalc = tToF[nh] - tS1;
		// Calculate x, y z positions relative to S1
		double positionX = ((xToF[nh] - 4.) / 152.)*(s3EndX - s3StartX) + s3StartX;;
		double positionY = ((xToF[nh] - 4.) / 152.)*(s3s1EndY - s3s1StartY) + s3s1StartY;
		double positionZ = (yToF[nh] + s3BarBottom + 2.75) / 100.;
		double angleTheta = TMath::ATan(positionX / positionY) * (180./TMath::Pi());
		double anglePhi   = TMath::ATan(positionZ / positionY) * (180./TMath::Pi());
		// All triggers
		//	if (tTrig == 0) {
		hutof1dS1->Fill(tofCalc, 1./deadtimeWeight);
		// Separate protons and MIPs using timing and amplitude cuts
		// Is a MIP
		if ( tofCalc > piLow && tofCalc < piHi ) {
		  nPi++;
		  hThetaS1pi->Fill(angleTheta, 1./deadtimeWeight);
		  hPhiS1pi->Fill(anglePhi, 1./deadtimeWeight);
		  hpionXY->Fill(xToF[nh], yToF[nh]);
		  lastut = t;
		} // if ( tofCalc > (tLight - (piLow+piHi)/2.) + piLow && tofCalc < (tLight - (piLow+piHi)/2.) + piHi )
		// Is a proton
		else if ( tofCalc > proLow && tofCalc < proHi && A1ToF[nh] > A1CutVec[nBar[nh]] && A2ToF[nh] > A2CutVec[nBar[nh]]) {
		  nP++;
		  hThetaS1pro->Fill(angleTheta, 1./deadtimeWeight);
		  hPhiS1pro->Fill(anglePhi, 1./deadtimeWeight);
		  hprotonXY->Fill(xToF[nh], yToF[nh]);
		  hMomS1->Fill(momFromTime(0.938, 10.9, tofCalc), 1./deadtimeWeight);
		  // hMomZS1->Fill(momFromTime(0.938, 10.9, tofCalc), nBar[nh]);
		  // hMomYS1->Fill(momFromTime(0.938, 10.9, tofCalc), xToF[nh]);
		  lastut = t;
		  if (nBlocks == 0 && momFromTime(0.938, 10.9, tofCalc) > 0.595) hMom2D_0blkQ->Fill(xToF[nh], nBar[nh]);
		  else if (nBlocks == 0 && momFromTime(0.938, 10.9, tofCalc) < 0.595) hMom2D_0blkS->Fill(xToF[nh], nBar[nh]);
		  else if (nBlocks == 1 && momFromTime(0.938, 10.9, tofCalc) > 0.570) hMom2D_1blkQ->Fill(xToF[nh], nBar[nh]);
		  else if (nBlocks == 1 && momFromTime(0.938, 10.9, tofCalc) < 0.570) hMom2D_1blkS->Fill(xToF[nh], nBar[nh]);
		  else if (nBlocks == 2 && momFromTime(0.938, 10.9, tofCalc) > 0.525) hMom2D_2blkQ->Fill(xToF[0], nBar[nh]);
		  else if (nBlocks == 2 && momFromTime(0.938, 10.9, tofCalc) < 0.525) hMom2D_2blkS->Fill(xToF[nh], nBar[nh]);
	    
		} // else if ( tofCalc > (tLight - (piLow+piHi)/2.) + proLow && tofCalc < (tLight - (piLow+piHi)/2.) + proHi )
		//	}
		// S1 & S2 trigger only
		if (tTrig !=0) {
		  hutof1dS1S2->Fill(tofCalc, 1./deadtimeWeight);
		  // Separate protons and MIPs using timing and amplitude cuts
		  // Is a MIP
		  if ( tofCalc > piLow && tofCalc < piHi ) {
		    nPi++;
		    hThetaS1S2pi->Fill(angleTheta, 1./deadtimeWeight);
		    hPhiS1S2pi->Fill(anglePhi, 1./deadtimeWeight);
		  } // if ( tofCalc > (tLight - (piLow+piHi)/2.) + piLow && tofCalc < (tLight - (piLow+piHi)/2.) + piHi )
		  // Is a proton
		  else if ( tofCalc > proLow && tofCalc < proHi && A1ToF[nh] > A1CutVec[nBar[nh]] && A2ToF[nh] > A2CutVec[nBar[nh]]) {
		    nP++;
		    hThetaS1S2pro->Fill(angleTheta, 1./deadtimeWeight);
		    hPhiS1S2pro->Fill(anglePhi, 1./deadtimeWeight);
		    hMomS1S2->Fill(momFromTime(0.938, 10.9, tofCalc), 1./deadtimeWeight);
		    // hMomZS12->Fill(momFromTime(0.938, 10.9, tofCalc), nBar[nh]);
		    // hMomYS12->Fill(momFromTime(0.938, 10.9, tofCalc), xToF[nh]);
		  } // else if ( tofCalc > (tLight - (piLow+piHi)/2.) + proLow && tofCalc < (tLight - (piLow+piHi)/2.) + proHi )
		} // S1 + S2 trigger
	      }
	    } // Loop over nhits
	  } // for (int t=0; t<tree->GetEntries(); t++)
	} // if (isGood)
      } 

      cout<<"Utof spills "<<nSpills<<" vs. "<<utofTimes.size()<<" in spill DB"<<endl;
      fout->cd();

      hThetaS1S2ratio->Divide(hThetaS1S2pro, hThetaS1S2pi, 1., 1., "B");
      hPhiS1S2ratio->Divide(hPhiS1S2pro, hPhiS1S2pi, 1., 1., "B");
      hThetaS1ratio->Divide(hThetaS1pro, hThetaS1pi, 1., 1., "B");
      hPhiS1ratio->Divide(hPhiS1pro, hPhiS1pi, 1., 1., "B");
      hThetaS1S2ratio->Write();
      hPhiS1S2ratio->Write();
      hThetaS1ratio->Write();
      hPhiS1ratio->Write();

      hThetaS1S2pro->SetLineWidth(2);
      hThetaS1S2pi->SetLineWidth(2);
      hPhiS1S2pro->SetLineWidth(2);
      hPhiS1S2pi->SetLineWidth(2);
      hPhiS1S2ratio->SetLineWidth(2);
      hThetaS1S2ratio->SetLineWidth(2);
      hThetaS1pro->SetLineWidth(2);
      hThetaS1pi->SetLineWidth(2);
      hPhiS1pro->SetLineWidth(2);
      hPhiS1pi->SetLineWidth(2);
      hPhiS1ratio->SetLineWidth(2);
      hThetaS1ratio->SetLineWidth(2);
      hMomS1S2->SetLineWidth(2);
      hMomS1->SetLineWidth(2);
      hutof1dS1S2->SetLineWidth(2);
      hutof1dS1->SetLineWidth(2);
    
      if (nBlocks==0) {
	hThetaS1S2pro->SetLineColor(kBlack);
	hThetaS1S2pi->SetLineColor(kBlack);
	hPhiS1S2pro->SetLineColor(kBlack);
	hPhiS1S2pi->SetLineColor(kBlack);
	hPhiS1S2ratio->SetLineColor(kBlack);
	hThetaS1S2ratio->SetLineColor(kBlack);
	hThetaS1pro->SetLineColor(kBlack);
	hThetaS1pi->SetLineColor(kBlack);
	hPhiS1pro->SetLineColor(kBlack);
	hPhiS1pi->SetLineColor(kBlack);
	hPhiS1ratio->SetLineColor(kBlack);
	hThetaS1ratio->SetLineColor(kBlack);
	hMomS1S2->SetLineColor(kBlack);
	hMomS1->SetLineColor(kBlack);
	hutof1dS1S2->SetLineColor(kBlack);
	hutof1dS1->SetLineColor(kBlack);
      
	leg->AddEntry(hThetaS1S2pro, "0 blocks", "l");
	legTheta->AddEntry(hThetaS1ratio, "0 blocks", "l");
	legTof->AddEntry(hutof1dS1, "0 blocks", "l");
	legPhiRatio->AddEntry(hPhiS1ratio, "0 blocks", "l");

	hMom2D_0blkQ->Scale(1. / (double)nSpills);
	hMom2D_0blkS->Scale(1. / (double)nSpills);
	hMom2D_0blkQ->Write();
	hMom2D_0blkS->Write();
      }
      if (nBlocks==1) {
	hThetaS1S2pro->SetLineColor(kRed);
	hThetaS1S2pi->SetLineColor(kRed);
	hPhiS1S2pro->SetLineColor(kRed);
	hPhiS1S2pi->SetLineColor(kRed);
	hPhiS1S2ratio->SetLineColor(kRed);
	hThetaS1S2ratio->SetLineColor(kRed);
	hThetaS1pro->SetLineColor(kRed);
	hThetaS1pi->SetLineColor(kRed);
	hPhiS1pro->SetLineColor(kRed);
	hPhiS1pi->SetLineColor(kRed);
	hPhiS1ratio->SetLineColor(kRed);
	hThetaS1ratio->SetLineColor(kRed);
	hMomS1S2->SetLineColor(kRed);
	hMomS1->SetLineColor(kRed);
	hutof1dS1S2->SetLineColor(kRed);
	hutof1dS1->SetLineColor(kRed);
      
	leg->AddEntry(hThetaS1S2pro, "1 block", "l");
	legTheta->AddEntry(hThetaS1ratio, "1 block", "l");
	legTof->AddEntry(hutof1dS1, "1 block", "l");
	legPhiRatio->AddEntry(hPhiS1ratio, "1 block", "l");

	hMom2D_1blkQ->Scale(1. / (double)nSpills);
	hMom2D_1blkS->Scale(1. / (double)nSpills);
	hMom2D_1blkQ->Write();
	hMom2D_1blkS->Write();
      }
      if (nBlocks==2) {
	hThetaS1S2pro->SetLineColor(kBlue);
	hThetaS1S2pi->SetLineColor(kBlue);
	hPhiS1S2pro->SetLineColor(kBlue);
	hPhiS1S2pi->SetLineColor(kBlue);
	hPhiS1S2ratio->SetLineColor(kBlue);
	hThetaS1S2ratio->SetLineColor(kBlue);

	hThetaS1pro->SetLineColor(kBlue);
	hThetaS1pi->SetLineColor(kBlue);
	hPhiS1pro->SetLineColor(kBlue);
	hPhiS1pi->SetLineColor(kBlue);
	hPhiS1ratio->SetLineColor(kBlue);
	hThetaS1ratio->SetLineColor(kBlue);

	hMomS1S2->SetLineColor(kBlue);
	hMomS1->SetLineColor(kBlue);

	hutof1dS1S2->SetLineColor(kBlue);
	hutof1dS1->SetLineColor(kBlue);

	leg->AddEntry(hThetaS1S2pro, "2 blocks", "l");
	legTheta->AddEntry(hThetaS1ratio, "2 blocks", "l");
	legTof->AddEntry(hutof1dS1, "2 blocks", "l");
	legPhiRatio->AddEntry(hPhiS1ratio, "2 blocks", "l");

	hMom2D_2blkQ->Scale(1. / (double)nSpills);
	hMom2D_2blkS->Scale(1. / (double)nSpills);
	hMom2D_2blkQ->Write();
	hMom2D_2blkS->Write();
      }
      if (nBlocks==3) {
	hThetaS1S2pro->SetLineColor(kCyan+1);
	hThetaS1S2pi->SetLineColor(kCyan+1);
	hPhiS1S2pro->SetLineColor(kCyan+1);
	hPhiS1S2pi->SetLineColor(kCyan+1);
	hPhiS1S2ratio->SetLineColor(kCyan+1);
	hThetaS1S2ratio->SetLineColor(kCyan+1);

	hThetaS1pro->SetLineColor(kCyan+1);
	hThetaS1pi->SetLineColor(kCyan+1);
	hPhiS1pro->SetLineColor(kCyan+1);
	hPhiS1pi->SetLineColor(kCyan+1);
	hPhiS1ratio->SetLineColor(kCyan+1);
	hThetaS1ratio->SetLineColor(kCyan+1);

	hMomS1S2->SetLineColor(kCyan+1);
	hMomS1->SetLineColor(kCyan+1);

	hutof1dS1S2->SetLineColor(kCyan+1);
	hutof1dS1->SetLineColor(kCyan+1);

	leg->AddEntry(hThetaS1S2pro, "3 blocks", "l");
	legTheta->AddEntry(hThetaS1ratio, "3 blocks", "l");
	legTof->AddEntry(hutof1dS1, "3 blocks", "l");
	legPhiRatio->AddEntry(hPhiS1ratio, "3 blocks", "l");
      }
      if (nBlocks==4) {
	hThetaS1S2pro->SetLineColor(kOrange+1);
	hThetaS1S2pi->SetLineColor(kOrange+1);
	hPhiS1S2pro->SetLineColor(kOrange+1);
	hPhiS1S2pi->SetLineColor(kOrange+1);
	hPhiS1S2ratio->SetLineColor(kOrange+1);
	hThetaS1S2ratio->SetLineColor(kOrange+1);

	hThetaS1pro->SetLineColor(kOrange+1);
	hThetaS1pi->SetLineColor(kOrange+1);
	hPhiS1pro->SetLineColor(kOrange+1);
	hPhiS1pi->SetLineColor(kOrange+1);
	hPhiS1ratio->SetLineColor(kOrange+1);
	hThetaS1ratio->SetLineColor(kOrange+1);

	hMomS1S2->SetLineColor(kOrange+1);
	hMomS1->SetLineColor(kOrange+1);

	hutof1dS1S2->SetLineColor(kOrange+1);
	hutof1dS1->SetLineColor(kOrange+1);

	leg->AddEntry(hThetaS1S2pro, "4 blocks", "l");
	legTheta->AddEntry(hThetaS1ratio, "4 blocks", "l");
	legTof->AddEntry(hutof1dS1, "4 blocks", "l");
	legPhiRatio->AddEntry(hPhiS1ratio, "4 blocks", "l");
      }

      hPhiS1S2pro->Scale(1. / (double)nSpills);
      hPhiS1S2pi->Scale(1. / (double)nSpills);
      hThetaS1S2pro->Scale(1. / (double)nSpills);
      hThetaS1S2pi->Scale(1. / (double)nSpills);

      hPhiS1pro->Scale(1. / (double)nSpills);
      hPhiS1pi->Scale(1. / (double)nSpills);
      hPhiS1pro->Scale(22./6.4);
      hPhiS1pi->Scale(22./6.4);
      hPhiS1S2pro->Scale(22./6.4);
      hPhiS1S2pi->Scale(22./6.4);
      hThetaS1pro->Scale(1. / (double)nSpills);
      hThetaS1pi->Scale(1. / (double)nSpills);
      for (int i=1; i < hThetaS1pi->GetNbinsX(); i++) {
	double binWidth = (hThetaS1pi->GetXaxis()->GetBinUpEdge(i) - hThetaS1pi->GetXaxis()->GetBinLowEdge(i));
	hThetaS1pi->SetBinContent(i, hThetaS1pi->GetBinContent(i) / binWidth);
	hThetaS1pro->SetBinContent(i, hThetaS1pro->GetBinContent(i) / binWidth);
	hThetaS1S2pi->SetBinContent(i, hThetaS1S2pi->GetBinContent(i) / binWidth);
	hThetaS1S2pro->SetBinContent(i, hThetaS1S2pro->GetBinContent(i) / binWidth);
      }

      hMomS1S2->Scale(1. / (double)nSpills);
      hMomS1->Scale(1. / (double)nSpills);

      hutof1dS1S2->Scale(1. / (double)nSpills);
      hutof1dS1->Scale(1. / (double)nSpills);

      // hMomYS1->Scale(1. / (double)nSpills);
      // hMomYS12->Scale(1. / (double)nSpills);
      // hMomZS1->Scale(1. / (double)nSpills);
      // hMomZS12->Scale(1. / (double)nSpills);

      hprotonXY->Scale(1. / (double)nSpills);
      hpionXY->Scale(1. / (double)nSpills);
      hprotonXY->Write();
      hpionXY->Write();

      gStyle->SetPalette(55);
      gStyle->SetOptStat(0);
      TCanvas *cpionXY   = new TCanvas(Form("cpionXY%d", nBlocks));
      hpionXY->Draw("colz");
      cpionXY->SetRightMargin(0.15);
      cpionXY->SetLeftMargin(0.13);
      cpionXY->SetBottomMargin(0.13);
      hpionXY->GetXaxis()->SetLabelSize(0.05);
      hpionXY->GetYaxis()->SetLabelSize(0.05);
      hpionXY->GetXaxis()->SetTitleSize(0.05);
      hpionXY->GetYaxis()->SetTitleSize(0.05);
      hpionXY->GetZaxis()->SetTitleSize(0.05);
      hpionXY->GetZaxis()->SetLabelSize(0.05);
      cpionXY->Print(Form("%s/%d_pionXY.png", saveDir, nBlocks));
      cpionXY->Print(Form("%s/%d_pionXY.pdf", saveDir, nBlocks));
      cpionXY->Print(Form("%s/%d_pionXY.tex", saveDir, nBlocks));
      TCanvas *cprotonXY = new TCanvas(Form("cprotonXY%d", nBlocks));
      hprotonXY->Draw("colz");
      cprotonXY->SetRightMargin(0.15);
      cprotonXY->SetLeftMargin(0.13);
      cprotonXY->SetBottomMargin(0.13);
      hprotonXY->GetXaxis()->SetLabelSize(0.05);
      hprotonXY->GetYaxis()->SetLabelSize(0.05);
      hprotonXY->GetXaxis()->SetTitleSize(0.05);
      hprotonXY->GetYaxis()->SetTitleSize(0.05);
      hprotonXY->GetZaxis()->SetTitleSize(0.05);
      hprotonXY->GetZaxis()->SetLabelSize(0.05);
      cprotonXY->Print(Form("%s/%d_protonXY.png", saveDir, nBlocks));
      cprotonXY->Print(Form("%s/%d_protonXY.pdf", saveDir, nBlocks));
      cprotonXY->Print(Form("%s/%d_protonXY.tex", saveDir, nBlocks));

      hThetaS1S2pro->Write();
      hThetaS1S2pi->Write();
      hPhiS1S2pro->Write();
      hPhiS1S2pi->Write();
      hThetaS1pro->Write();
      hThetaS1pi->Write();
      hPhiS1pro->Write();
      hPhiS1pi->Write();

      hMomS1S2->Write();
      hMomS1->Write();

      hutof1dS1S2->Write();
      hutof1dS1->Write();

      hsThetaS1S2pro->Add(hThetaS1S2pro);
      hsThetaS1S2pi->Add(hThetaS1S2pi);
      hsPhiS1S2pro->Add(hPhiS1S2pro);
      hsPhiS1S2pi->Add(hPhiS1S2pi);
      hsThetaS1S2ratio->Add(hThetaS1S2ratio);
      hsPhiS1S2ratio->Add(hPhiS1S2ratio);

      hsThetaS1pro->Add(hThetaS1pro);
      hsThetaS1pi->Add(hThetaS1pi);
      hsPhiS1pro->Add(hPhiS1pro);
      hsPhiS1pi->Add(hPhiS1pi);
      hsThetaS1ratio->Add(hThetaS1ratio);
      hsPhiS1ratio->Add(hPhiS1ratio);

      hsutof1dS1S2->Add(hutof1dS1S2);
      hsutof1dS1->Add(hutof1dS1);

      hsMomS1S2->Add(hMomS1S2);
      hsMomS1->Add(hMomS1);
    } // if (nBlocks != 4)
    else {
      // Loop through 4 block data
      for (int j=0; j<str4BlockVec.size(); j++) {
	int nSpillsTmp = 0;
	std::cout<<"Analysing 4 block sample "<<j+1<<" of "<<str4BlockVec.size()<<std::endl;
	TH1D *tofTmp = new TH1D(Form("tofTmp%s", str4BlockVec.at(j)), Form("S3 ToF %s", str4BlockVec.at(j)), 250, 25, 125);
	TH1D *hThetaS1piTmp = new TH1D(Form("hThetaS1piTmp%s", str4BlockVec.at(j)), Form("Angular distribution of pion hits in S3 (S1 trigger only), %d blocks; #theta / degrees; Events / spill / degree", nBlocks), binnum, binsTheta);
	TH1D *hThetaS1proTmp = new TH1D(Form("hThetaS1proTmp%s", str4BlockVec.at(j)), Form("Angular distribution of proton hits in S3 (S1 trigger only), %d blocks; #theta / degrees; Events / spill / degree", nBlocks), binnum, binsTheta);
	TH1D *hThetaS1ratioTmp = new TH1D(Form("hThetaS1ratioTmp%s", str4BlockVec.at(j)), Form("Angular distribution of proton:MIP ratio in S3 (S1 trigger only), %d blocks; #theta / degrees; Ratio", nBlocks), binnum, binsTheta);
	// Number of protons and number of MIPs
	int nP  = 0;
	int nPi = 0;
	// Define signal and background functions to be fitted
	TF1 *sPro = new TF1("sPro", "gaus", proLow, proHi);
	TF1 *sPi  = new TF1("sPi", "gaus", piLow, piHi);

	// Find the correct dstof files
	Int_t runMin=-1;
	Int_t runMax=-1;

	double startTime = 0;
	double endTime   = 0;

	const char* nustof;
	double slope;
	double constant;

	nustof = str4BlockVec.at(j);
	slope    = block4SlopeVec.at(j);
	constant = block4ConstVec.at(j);

	TFile *futof = new TFile(Form("%s/%s",ustofDir,nustof), "read");
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
	TTree *tree = (TTree*)futof->Get("tree");
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
	TNamed *start = 0;
	TNamed *end   = 0;
	futof->GetObject("start_of_run", start);
	futof->GetObject("end_of_run", end);
	const char* startchar = start->GetTitle();
	std::string startstr(startchar);
	std::string unixstart = startstr.substr(25,10);
	int startTimeUtof = stoi(unixstart);
	const char* endchar = end->GetTitle();
	std::string endstr(endchar);
	std::string unixend = endstr.substr(23,10);
	int endTimeUtof = stoi(unixend);

	startTime = startTimeUtof;
	endTime = endTimeUtof;

	// Find dtof runs
	for (int irun=950; irun<1400; irun++) {
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
    
	cout << "Min and max dtof runs are " << runMin << " " << runMax << endl;

	std::vector<double> dtofTimes;
	std::vector<int> dtofS1S2Hits;
	std::vector<double> utofTimes;
	std::vector<int> utofS1S2Hits;

	// Open the appropriate spill DB files and get the spill times
	for (int irun = runMin; irun < runMax+1; irun++) {
	  TFile *dbFile = new TFile(Form("%s/spillDB_run%d_run%d.root", spillDir, irun, irun), "read");
	  TTree *spillTree = (TTree*)dbFile->Get("spillTree");
	  double globalSpillTime;
	  double ustofSpillTime;
	  spillTree->SetBranchAddress("globalSpillTime", &globalSpillTime);
	  spillTree->SetBranchAddress("ustofSpillTime", &ustofSpillTime);
	  for (int t = 0; t < spillTree->GetEntries(); t++) {
	    spillTree->GetEntry(t);
	    if (globalSpillTime >= startTime && globalSpillTime <= endTime) {
	      dtofTimes.push_back(globalSpillTime);
	      utofTimes.push_back(ustofSpillTime);
	    }
	  } // for (int t = 0; t < spillTree->GetEntries(); t++)

	  dbFile->Close();
	  delete dbFile;
	} // for (int irun = runMin; irun < runMax+1; irun++) 

	dtofS1S2Hits.resize(dtofTimes.size(), 0);
	utofS1S2Hits.resize(dtofTimes.size(), 0);

	cout<<"Finding number of dtof hits in each spill"<<endl;
	int lastt = 0;
	int lastrun = 0;
	for (int s=0; s<dtofTimes.size(); s++) {
	  if (s % 25 == 0) {
	    cout<<"Spill "<<s<<" of "<<dtofTimes.size()<<endl;
	  }
	  // Loop over the all the files
	  for (int irun = runMin; irun < runMax+1; irun++) {
	
	    TFile *dtofFile = new TFile(Form("%srun%d/DsTOFtreeRun%d_tdc1.root", dstofDir, irun, irun), "read");
	    RawDsTofHeader *tof = NULL;
	    TTree *tofTree = (TTree*)dtofFile->Get("tofTree");
	    tofTree->SetBranchAddress("tof", &tof);
	    tofTree->GetEntry(0);
	    double firstTime = tof->unixTime;
	    tofTree->GetEntry(tofTree->GetEntries()-1);
	    double lastTime = tof->unixTime;
	    double lastS1S2Dtof = 0.;
	    // Spill is in this file
	    if (firstTime <= dtofTimes[s] && lastTime >= dtofTimes[s]) {
	      if (irun != lastrun) {
		lastt = 0;
		lastrun = irun;
	      }
	      // Loop over all entries and count the number of S1 S2 hits within the spill
	      for (int t=lastt; t<tofTree->GetEntries(); t++) {
		tofTree->GetEntry(t);
		if ((tof->fakeTimeNs/1e9)+firstTime < dtofTimes[s]) continue;
		if ((tof->fakeTimeNs/1e9)+firstTime > dtofTimes[s]+1.) break;
		// Is within a spill
		if ((tof->fakeTimeNs/1e9)+firstTime >= dtofTimes[s] &&
		    (tof->fakeTimeNs/1e9)+firstTime <= dtofTimes[s]+1. && 
		    tof->channel == 13 && (tof->fakeTimeNs - lastS1S2Dtof) > 500.) {
		  dtofS1S2Hits[s]++;
		  lastt = t;
		  lastS1S2Dtof = tof->fakeTimeNs;
		} // Is within the spill
	      } // for (int t=0; t<tofTree->GetEntries(); t++)
	    } // if (firstTime <= dtofTimes[s] && lastTime >=  dtofTimes[s])

	    delete tof;
	    dtofFile->Close();
	    delete dtofFile;
	  } // for (int irun = runMin; irun < runMax+1; irun++) 
	} // for (int s=0; s<dtofTimes.size(); s++)
	double lastSpill = 0.; 
	int lastut = 0;
       	// Loop over the spills and perform the adjustment for each spill
	for (int s = 0; s < utofTimes.size(); s++) {
	  if (s % 25 == 0) cout<<"Getting hits from spill "<<s<<" of "<<utofTimes.size()<<endl;
	  double deadtimeWeight = dtofS1S2Hits[s] * slope + constant;
	  // Do initial loop to check data quality
	  bool isGood = false;
	  for (int t=lastut; t<tree->GetEntries(); t++) {
	    tree->GetEntry(t);
	    if ((tS1/1e9) + startTimeUtof < utofTimes[s]) continue;
	    if ((tS1/1e9) + startTimeUtof > utofTimes[s] + 1.) break;
	    // Quality for the 
	    if ((tTrig/1e9) + startTimeUtof > utofTimes[s] + 0.47) {
	      isGood = true;
	      nSpillsTmp++;
	      break;
	    }
	  } // Initial loop to check data quality
	  if (isGood) {
	    for (int t=lastut; t<tree->GetEntries(); t++) {
	      tree->GetEntry(t);
	      if ((tS1/1e9) + startTimeUtof < utofTimes[s]) continue;
	      if ((tS1/1e9) + startTimeUtof > utofTimes[s] + 1.) break;
	      // Count number of spills
	      // if (tSoSd != lastSpill && tSoSd -lastSpill > 1e9) {
	      // lastSpill = tSoSd;
	      // } // if (tSoSd != lastSpill && tSoSd -lastSpill > 1e9) 
	      for (int nh=0; nh < nhit; nh++) {
		// Checks if the next S3 hit is in a neighbouring bar
		// Only want to count one of them if it is
		// Check from this hit onwards to see if there are any double hit candidates
		bool isDouble = false;
		if (nh < nhit-1) {
		  for (int nh2=nh+1; nh2<nhit; nh2++) {
		    // Timing cut and checks if in the neighbouring bar
		    if (abs(tToF[nh] - tToF[nh2]) < twoBarCut && 
			(nBar[nh] == nBar[nh2]-1 || nBar[nh] == nBar[nh2]+1)) {
		      isDouble = true;
		      break;
		    }
		  } // for (int nh2=nh; nh2<nhit; nh2++)
		} // if (nh < nh-1)
		// If it's not a double hit then go ahead
		if (!isDouble) {
		  double tofCalc = tToF[nh] - tS1;
		  // Calculate x, y z positions relative to S1
		  double positionX = ((xToF[nh] - 4.) / 152.)*(s3EndX - s3StartX) + s3StartX;;
		  double positionY = ((xToF[nh] - 4.) / 152.)*(s3s1EndY - s3s1StartY) + s3s1StartY;
		  double positionZ = (yToF[nh] + s3BarBottom + 2.75) / 100.;
		  double angleTheta = TMath::ATan(positionX / positionY) * (180./TMath::Pi());
		  double anglePhi   = TMath::ATan(positionZ / positionY) * (180./TMath::Pi());
		  // All triggers
		  hutof1dS1->Fill(tofCalc, 1./deadtimeWeight);
		  tofTmp->Fill(tofCalc, 1./deadtimeWeight);
		  // Separate protons and MIPs using timing and amplitude cuts
		  // Is a MIP
		  if ( tofCalc > piLow && tofCalc < piHi ) {
		    nPi++;
		    hThetaS1pi->Fill(angleTheta, 1./deadtimeWeight);
		    hPhiS1pi->Fill(anglePhi, 1./deadtimeWeight);
		    hpionXY->Fill(xToF[nh], yToF[nh], 1./deadtimeWeight);
		    lastut = t;
		  } // if ( tofCalc > piLow && tofCalc < piHi )
		  // Is a proton
		  else if ( tofCalc > proLow && tofCalc < proHi && A1ToF[nh] > A1CutVec[nBar[nh]] && A2ToF[nh] > A2CutVec[nBar[nh]]) {
		    nP++;
		    hThetaS1pro->Fill(angleTheta, 1./deadtimeWeight);
		    hThetaS1proTmp->Fill(angleTheta, 1./deadtimeWeight);
		    hThetaS1piTmp->Fill(angleTheta, 1./deadtimeWeight);
		    hPhiS1pro->Fill(anglePhi, 1./deadtimeWeight);
		    hprotonXY->Fill(xToF[nh], yToF[nh], 1./deadtimeWeight);
		    hMomS1->Fill(momFromTime(0.938, 10.9, tofCalc), 1./deadtimeWeight);
		    // hMomZS1->Fill(momFromTime(0.938, 10.9, tofCalc), nBar[nh]);
		    // hMomYS1->Fill(momFromTime(0.938, 10.9, tofCalc), xToF[nh]);
		    lastut = t;
		    if (nBlocks == 0 && momFromTime(0.938, 10.9, tofCalc) > 0.595) hMom2D_0blkQ->Fill(xToF[nh], nBar[nh]);
		    else if (nBlocks == 0 && momFromTime(0.938, 10.9, tofCalc) < 0.595) hMom2D_0blkS->Fill(xToF[nh], nBar[nh]);
		    else if (nBlocks == 1 && momFromTime(0.938, 10.9, tofCalc) > 0.570) hMom2D_1blkQ->Fill(xToF[nh], nBar[nh]);
		    else if (nBlocks == 1 && momFromTime(0.938, 10.9, tofCalc) < 0.570) hMom2D_1blkS->Fill(xToF[nh], nBar[nh]);
		    else if (nBlocks == 2 && momFromTime(0.938, 10.75, tofCalc) > 0.525) hMom2D_2blkQ->Fill(xToF[nh], nBar[nh]);
		    else if (nBlocks == 2 && momFromTime(0.938, 10.75, tofCalc) < 0.525) hMom2D_2blkS->Fill(xToF[nh], nBar[nh]);
	    
		  } // else if ( tofCalc > (tLight - (piLow+piHi)/2.) + proLow && tofCalc < (tLight - (piLow+piHi)/2.) + proHi )
		  //	}
		  // S1 & S2 trigger only
		  if (tTrig !=0) {
		    hutof1dS1S2->Fill(tofCalc, 1./deadtimeWeight);
		    // Separate protons and MIPs using timing and amplitude cuts
		    // Is a MIP
		    if ( tofCalc > piLow && tofCalc < piHi ) {
		      nPi++;
		      hThetaS1S2pi->Fill(angleTheta, 1./deadtimeWeight);
		      hPhiS1S2pi->Fill(anglePhi, 1./deadtimeWeight);
		    } // if ( tofCalc > (tLight - (piLow+piHi)/2.) + piLow && tofCalc < (tLight - (piLow+piHi)/2.) + piHi )
		    // Is a proton
		    else if ( tofCalc > proLow && tofCalc < proHi && A1ToF[nh] > A1CutVec[nBar[nh]] && A2ToF[nh] > A2CutVec[nBar[nh]]) {
		      nP++;
		      hThetaS1S2pro->Fill(angleTheta, 1./deadtimeWeight);
		      hPhiS1S2pro->Fill(anglePhi, 1./deadtimeWeight);
		      hMomS1S2->Fill(momFromTime(0.938, 10.9, tofCalc), 1./deadtimeWeight);
		      // hMomZS12->Fill(momFromTime(0.938, 10.9, tofCalc), nBar[nh]);
		      // hMomYS12->Fill(momFromTime(0.938, 10.9, tofCalc), xToF[nh]);
		    } // else if ( tofCalc > (tLight - (piLow+piHi)/2.) + proLow && tofCalc < (tLight - (piLow+piHi)/2.) + proHi )
		  } // S1 + S2 trigger
		} // if (!isDouble)
	      } // Loop over nhits
	    } // for (int t=lastut; t<tree->GetEntries(); t++)
	  } // if (isGood)
	}
	nSpills += nSpillsTmp;
	fout->cd();
	for (int i=1; i < hThetaS1pi->GetNbinsX(); i++) {
	  double binWidth = (hThetaS1piTmp->GetXaxis()->GetBinUpEdge(i) - hThetaS1piTmp->GetXaxis()->GetBinLowEdge(i));
	  hThetaS1piTmp->SetBinContent(i, hThetaS1piTmp->GetBinContent(i) / binWidth);
	  hThetaS1proTmp->SetBinContent(i, hThetaS1proTmp->GetBinContent(i) / binWidth);
	}
	tofTmp->Scale(1. / (double)nSpillsTmp);
	hThetaS1proTmp->Scale(1. / (double)nSpillsTmp);
	hThetaS1piTmp->Scale(1. / (double)nSpillsTmp);
	hThetaS1ratioTmp->Divide(hThetaS1proTmp, hThetaS1piTmp, 1., 1., "B");
	tofTmp->Write();
	hThetaS1proTmp->Write();
	hThetaS1piTmp->Write();
	hThetaS1ratioTmp->Write();
	cout<<nSpillsTmp<<" good utof spills out of "<<utofTimes.size()<<endl;
      } // for (int j=0; j<str4BlockVec.size(); j++)
      fout->cd();
      hThetaS1S2ratio->Divide(hThetaS1S2pro, hThetaS1S2pi, 1., 1., "B");
      hPhiS1S2ratio->Divide(hPhiS1S2pro, hPhiS1S2pi, 1., 1., "B");
      hThetaS1ratio->Divide(hThetaS1pro, hThetaS1pi, 1., 1., "B");
      hPhiS1ratio->Divide(hPhiS1pro, hPhiS1pi, 1., 1., "B");
      hThetaS1S2ratio->Write();
      hPhiS1S2ratio->Write();
      hThetaS1ratio->Write();
      hPhiS1ratio->Write();

      hThetaS1S2pro->SetLineWidth(2);
      hThetaS1S2pi->SetLineWidth(2);
      hPhiS1S2pro->SetLineWidth(2);
      hPhiS1S2pi->SetLineWidth(2);
      hPhiS1S2ratio->SetLineWidth(2);
      hThetaS1S2ratio->SetLineWidth(2);
      hThetaS1pro->SetLineWidth(2);
      hThetaS1pi->SetLineWidth(2);
      hPhiS1pro->SetLineWidth(2);
      hPhiS1pi->SetLineWidth(2);
      hPhiS1ratio->SetLineWidth(2);
      hThetaS1ratio->SetLineWidth(2);
      hMomS1S2->SetLineWidth(2);
      hMomS1->SetLineWidth(2);
      hutof1dS1S2->SetLineWidth(2);
      hutof1dS1->SetLineWidth(2);
    
      hThetaS1S2pro->SetLineColor(kOrange+1);
      hThetaS1S2pi->SetLineColor(kOrange+1);
      hPhiS1S2pro->SetLineColor(kOrange+1);
      hPhiS1S2pi->SetLineColor(kOrange+1);
      hPhiS1S2ratio->SetLineColor(kOrange+1);
      hThetaS1S2ratio->SetLineColor(kOrange+1);

      hThetaS1pro->SetLineColor(kOrange+1);
      hThetaS1pi->SetLineColor(kOrange+1);
      hPhiS1pro->SetLineColor(kOrange+1);
      hPhiS1pi->SetLineColor(kOrange+1);
      hPhiS1ratio->SetLineColor(kOrange+1);
      hThetaS1ratio->SetLineColor(kOrange+1);

      hMomS1S2->SetLineColor(kOrange+1);
      hMomS1->SetLineColor(kOrange+1);

      hutof1dS1S2->SetLineColor(kOrange+1);
      hutof1dS1->SetLineColor(kOrange+1);

      leg->AddEntry(hThetaS1S2pro, "4 blocks", "l");
      legTheta->AddEntry(hThetaS1ratio, "4 blocks", "l");
      legTof->AddEntry(hutof1dS1, "4 blocks", "l");
      legPhiRatio->AddEntry(hPhiS1ratio, "4 blocks", "l");

      hPhiS1S2pro->Scale(1. / (double)nSpills);
      hPhiS1S2pi->Scale(1. / (double)nSpills);
      hThetaS1S2pro->Scale(1. / (double)nSpills);
      hThetaS1S2pi->Scale(1. / (double)nSpills);

      hPhiS1pro->Scale(1. / (double)nSpills);
      hPhiS1pi->Scale(1. / (double)nSpills);
      hThetaS1pro->Scale(1. / (double)nSpills);
      hThetaS1pi->Scale(1. / (double)nSpills);

      hPhiS1pro->Scale(22./6.4);
      hPhiS1pi->Scale(22./6.4);
      hPhiS1S2pro->Scale(22./6.4);
      hPhiS1S2pi->Scale(22./6.4);
      for (int i=1; i < hThetaS1pi->GetNbinsX(); i++) {
	double binWidth = (hThetaS1pi->GetXaxis()->GetBinUpEdge(i) - hThetaS1pi->GetXaxis()->GetBinLowEdge(i));
	hThetaS1pi->SetBinContent(i, hThetaS1pi->GetBinContent(i) / binWidth);
	hThetaS1pro->SetBinContent(i, hThetaS1pro->GetBinContent(i) / binWidth);
	hThetaS1S2pi->SetBinContent(i, hThetaS1S2pi->GetBinContent(i) / binWidth);
	hThetaS1S2pro->SetBinContent(i, hThetaS1S2pro->GetBinContent(i) / binWidth);
      }

      hMomS1S2->Scale(1. / (double)nSpills);
      hMomS1->Scale(1. / (double)nSpills);

      hutof1dS1S2->Scale(1. / (double)nSpills);
      hutof1dS1->Scale(1. / (double)nSpills);

      // hMomYS1->Scale(1. / (double)nSpills);
      // hMomYS12->Scale(1. / (double)nSpills);
      // hMomZS1->Scale(1. / (double)nSpills);
      // hMomZS12->Scale(1. / (double)nSpills);

      hprotonXY->Scale(1. / (double)nSpills);
      hpionXY->Scale(1. / (double)nSpills);
      hprotonXY->Write();
      hpionXY->Write();

      gStyle->SetPalette(55);
      gStyle->SetOptStat(0);
      TCanvas *cpionXY   = new TCanvas(Form("cpionXY%d", nBlocks));
      hpionXY->Draw("colz");
      cpionXY->SetRightMargin(0.15);
      cpionXY->SetLeftMargin(0.13);
      cpionXY->SetBottomMargin(0.13);
      hpionXY->GetXaxis()->SetLabelSize(0.05);
      hpionXY->GetYaxis()->SetLabelSize(0.05);
      hpionXY->GetXaxis()->SetTitleSize(0.05);
      hpionXY->GetYaxis()->SetTitleSize(0.05);
      hpionXY->GetZaxis()->SetTitleSize(0.05);
      hpionXY->GetZaxis()->SetLabelSize(0.05);
      cpionXY->Print(Form("%s/%d_pionXY.png", saveDir, nBlocks));
      cpionXY->Print(Form("%s/%d_pionXY.pdf", saveDir, nBlocks));
      cpionXY->Print(Form("%s/%d_pionXY.tex", saveDir, nBlocks));
      TCanvas *cprotonXY = new TCanvas(Form("cprotonXY%d", nBlocks));
      hprotonXY->Draw("colz");
      cprotonXY->SetRightMargin(0.15);
      cprotonXY->SetLeftMargin(0.13);
      cprotonXY->SetBottomMargin(0.13);
      hprotonXY->GetXaxis()->SetLabelSize(0.05);
      hprotonXY->GetYaxis()->SetLabelSize(0.05);
      hprotonXY->GetXaxis()->SetTitleSize(0.05);
      hprotonXY->GetYaxis()->SetTitleSize(0.05);
      hprotonXY->GetZaxis()->SetTitleSize(0.05);
      hprotonXY->GetZaxis()->SetLabelSize(0.05);
      cprotonXY->Print(Form("%s/%d_protonXY.png", saveDir, nBlocks));
      cprotonXY->Print(Form("%s/%d_protonXY.pdf", saveDir, nBlocks));
      cprotonXY->Print(Form("%s/%d_protonXY.tex", saveDir, nBlocks));

      hThetaS1S2pro->Write();
      hThetaS1S2pi->Write();
      hPhiS1S2pro->Write();
      hPhiS1S2pi->Write();
      hThetaS1pro->Write();
      hThetaS1pi->Write();
      hPhiS1pro->Write();
      hPhiS1pi->Write();

      hMomS1S2->Write();
      hMomS1->Write();

      hutof1dS1S2->Write();
      hutof1dS1->Write();

      hsThetaS1S2pro->Add(hThetaS1S2pro);
      hsThetaS1S2pi->Add(hThetaS1S2pi);
      hsPhiS1S2pro->Add(hPhiS1S2pro);
      hsPhiS1S2pi->Add(hPhiS1S2pi);
      hsThetaS1S2ratio->Add(hThetaS1S2ratio);
      hsPhiS1S2ratio->Add(hPhiS1S2ratio);

      hsThetaS1pro->Add(hThetaS1pro);
      hsThetaS1pi->Add(hThetaS1pi);
      hsPhiS1pro->Add(hPhiS1pro);
      hsPhiS1pi->Add(hPhiS1pi);
      hsThetaS1ratio->Add(hThetaS1ratio);
      hsPhiS1ratio->Add(hPhiS1ratio);

      hsutof1dS1S2->Add(hutof1dS1S2);
      hsutof1dS1->Add(hutof1dS1);

      hsMomS1S2->Add(hMomS1S2);
      hsMomS1->Add(hMomS1);
      cout<<"Completed this dataset"<<endl;
    } // 4 block data
  } // for (int nBlocks = 0; nBlocks <= 4; nBlocks++) 

  fout->cd();
  leg->Write("leg");
  hsThetaS1S2pro->Write();
  hsThetaS1S2pi->Write();
  hsPhiS1S2pro->Write();
  hsPhiS1S2pi->Write();
  hsThetaS1S2ratio->Write();
  hsPhiS1S2ratio->Write();

  hsutof1dS1S2->Write();
  hsutof1dS1->Write();

  hsThetaS1pro->Write();
  hsThetaS1pi->Write();
  hsPhiS1pro->Write();
  hsPhiS1pi->Write();
  hsThetaS1ratio->Write();
  hsPhiS1ratio->Write();

  hsMomS1S2->Write();
  hsMomS1->Write();

  TCanvas *c1_Log1 = new TCanvas("c1_Log1");
  c1_Log1->SetLogy();
  hsThetaS1S2pro->Draw("hist e nostack");
  hsThetaS1S2pro->GetXaxis()->SetLabelSize(0.05);
  hsThetaS1S2pro->GetYaxis()->SetLabelSize(0.05);
  hsThetaS1S2pro->GetXaxis()->SetTitleSize(0.05);
  hsThetaS1S2pro->GetYaxis()->SetTitleSize(0.05);
  leg->Draw();
  c1_Log1->Print(Form("%s/thetaS12proLog.png", saveDir));
  c1_Log1->Print(Form("%s/thetaS12proLog.pdf", saveDir));
  c1_Log1->Print(Form("%s/thetaS12proLog.tex", saveDir));
  TCanvas *c1_Log2 = new TCanvas("c1_Log2");
  c1_Log2->SetLogy();
  hsThetaS1S2pi->Draw("hist e nostack");
  hsThetaS1S2pi->GetXaxis()->SetLabelSize();
  hsThetaS1S2pi->GetYaxis()->SetLabelSize();
  hsThetaS1S2pi->GetXaxis()->SetTitleSize();
  hsThetaS1S2pi->GetYaxis()->SetTitleSize();
  leg->Draw();
  c1_Log2->Print(Form("%s/thetaS12piLog.png", saveDir));
  c1_Log2->Print(Form("%s/thetaS12piLog.pdf", saveDir));
  c1_Log2->Print(Form("%s/thetaS12piLog.tex", saveDir));
  TCanvas *c1_Log3 = new TCanvas("c1_Log3");
  c1_Log3->SetLogy();
  hsPhiS1S2pro->Draw("hist e nostack");
  leg->Draw();
  c1_Log3->Print(Form("%s/phiS12proLog.png", saveDir));
  c1_Log3->Print(Form("%s/phiS12proLog.pdf", saveDir));
  c1_Log3->Print(Form("%s/phiS12proLog.tex", saveDir));
  TCanvas *c1_Log4 = new TCanvas("c1_Log4");
  c1_Log4->SetLogy();
  hsPhiS1S2pi->Draw("hist e nostack");
  leg->Draw();
  c1_Log4->Print(Form("%s/phiS12piLog.png", saveDir));
  c1_Log4->Print(Form("%s/phiS12piLog.pdf", saveDir));
  c1_Log4->Print(Form("%s/phiS12piLog.tex", saveDir));
  TCanvas *c1_Log5 = new TCanvas("c1_Log5");
  c1_Log5->SetLogy();
  hsThetaS1S2ratio->Draw("hist e nostack");
  leg->Draw();
  c1_Log5->Print(Form("%s/thetaS12ratioLog.png", saveDir));
  c1_Log5->Print(Form("%s/thetaS12ratioLog.pdf", saveDir));
  c1_Log5->Print(Form("%s/thetaS12ratioLog.tex", saveDir));
  TCanvas *c1_Log6 = new TCanvas("c1_Log6");
  c1_Log6->SetLogy();
  hsPhiS1S2ratio->Draw("hist e nostack");
  leg->Draw();
  c1_Log6->Print(Form("%s/phiS12ratioLog.png", saveDir));
  c1_Log6->Print(Form("%s/phiS12ratioLog.pdf", saveDir));
  c1_Log6->Print(Form("%s/phiS12ratioLog.tex", saveDir));

  TCanvas *c1_1 = new TCanvas("c1_1");
  hsThetaS1S2pro->Draw("hist e nostack");
  c1_1->SetGridx();
  c1_1->SetGridy();
  leg->Draw();
  c1_1->Print(Form("%s/thetaS12pro.png", saveDir));
  c1_1->Print(Form("%s/thetaS12pro.pdf", saveDir));
  c1_1->Print(Form("%s/thetaS12pro.tex", saveDir));
  TCanvas *c1_2 = new TCanvas("c1_2");
  hsThetaS1S2pi->Draw("hist e nostack");
  c1_2->SetGridx();
  c1_2->SetGridy();
  leg->Draw();
  c1_2->Print(Form("%s/thetaS12pi.png", saveDir));
  c1_2->Print(Form("%s/thetaS12pi.pdf", saveDir));
  c1_2->Print(Form("%s/thetaS12pi.tex", saveDir));
  TCanvas *c1_3 = new TCanvas("c1_3");
  hsPhiS1S2pro->Draw("hist e nostack");
  c1_3->SetGridx();
  c1_3->SetGridy();
  leg->Draw();
  c1_3->Print(Form("%s/phiS12pro.png", saveDir));
  c1_3->Print(Form("%s/phiS12pro.pdf", saveDir));
  c1_3->Print(Form("%s/phiS12pro.tex", saveDir));
  TCanvas *c1_4 = new TCanvas("c1_4");
  hsPhiS1S2pi->Draw("hist e nostack");
  c1_4->SetGridx();
  c1_4->SetGridy();
  leg->Draw();
  c1_4->Print(Form("%s/phiS12pi.png", saveDir));
  c1_4->Print(Form("%s/phiS12pi.pdf", saveDir));
  c1_4->Print(Form("%s/phiS12pi.tex", saveDir));
  TCanvas *c1_5 = new TCanvas("c1_5");
  hsThetaS1S2ratio->Draw("hist e nostack");
  c1_5->SetGridx();
  c1_5->SetGridy();
  leg->Draw();
  c1_5->Print(Form("%s/thetaS12ratio.png", saveDir));
  c1_5->Print(Form("%s/thetaS12ratio.pdf", saveDir));
  c1_5->Print(Form("%s/thetaS12ratio.tex", saveDir));
  TCanvas *c1_6 = new TCanvas("c1_6");
  hsPhiS1S2ratio->Draw("hist e nostack");
  c1_6->SetGridx();
  c1_6->SetGridy();
  leg->Draw();
  c1_6->Print(Form("%s/phiS12ratio.png", saveDir));
  c1_6->Print(Form("%s/phiS12ratio.pdf", saveDir));
  c1_6->Print(Form("%s/phiS12ratio.tex", saveDir));

  TCanvas *c1_Logs11 = new TCanvas("c1_Logs11");
  c1_Logs11->SetLogy();
  hsThetaS1pro->Draw("hist e nostack");
  hsThetaS1pro->GetXaxis()->SetLabelSize(0.05);
  hsThetaS1pro->GetYaxis()->SetLabelSize(0.05);
  hsThetaS1pro->GetXaxis()->SetTitleSize(0.05);
  hsThetaS1pro->GetYaxis()->SetTitleSize(0.05);
  legTof->Draw();
  c1_Logs11->Print(Form("%s/thetaS1proLog.png", saveDir));
  c1_Logs11->Print(Form("%s/thetaS1proLog.pdf", saveDir));
  c1_Logs11->Print(Form("%s/thetaS1proLog.tex", saveDir));
  TCanvas *c1_Logs12 = new TCanvas("c1_Logs12");
  c1_Logs12->SetLogy();
  hsThetaS1pi->Draw("hist e nostack");
  hsThetaS1pi->GetXaxis()->SetLabelSize(0.05);
  hsThetaS1pi->GetYaxis()->SetLabelSize(0.05);
  hsThetaS1pi->GetXaxis()->SetTitleSize(0.05);
  hsThetaS1pi->GetYaxis()->SetTitleSize(0.05);
  legTof->Draw();
  c1_Logs12->Print(Form("%s/thetaS1piLog.png", saveDir));
  c1_Logs12->Print(Form("%s/thetaS1piLog.pdf", saveDir));
  c1_Logs12->Print(Form("%s/thetaS1piLog.tex", saveDir));
  TCanvas *c1_Logs13 = new TCanvas("c1_Logs13");
  c1_Logs13->SetLogy();
  hsPhiS1pro->Draw("hist e nostack");
  hsPhiS1pro->GetXaxis()->SetLabelSize(0.05); 
  hsPhiS1pro->GetYaxis()->SetLabelSize(0.05);
  hsPhiS1pro->GetXaxis()->SetTitleSize(0.05);
  hsPhiS1pro->GetYaxis()->SetTitleSize(0.05);
  legTof->Draw();
  c1_Logs13->Print(Form("%s/phiS1proLog.png", saveDir));
  c1_Logs13->Print(Form("%s/phiS1proLog.pdf", saveDir));
  c1_Logs13->Print(Form("%s/phiS1proLog.tex", saveDir));
  TCanvas *c1_Logs14 = new TCanvas("c1_Logs14");
  c1_Logs14->SetLogy();
  hsPhiS1pi->Draw("hist e nostack");
  hsPhiS1pi->GetXaxis()->SetLabelSize(0.05); 
  hsPhiS1pi->GetYaxis()->SetLabelSize(0.05);
  hsPhiS1pi->GetXaxis()->SetTitleSize(0.05);
  hsPhiS1pi->GetYaxis()->SetTitleSize(0.05);
  legTof->Draw();
  c1_Logs14->Print(Form("%s/phiS1piLog.png", saveDir));
  c1_Logs14->Print(Form("%s/phiS1piLog.pdf", saveDir));
  c1_Logs14->Print(Form("%s/phiS1piLog.tex", saveDir));
  TCanvas *c1_Logs15 = new TCanvas("c1_Logs15");
  c1_Logs15->SetLogy();
  hsThetaS1ratio->Draw("hist e nostack");
  hsThetaS1ratio->GetXaxis()->SetLabelSize(0.05);
  hsThetaS1ratio->GetYaxis()->SetLabelSize(0.05);
  hsThetaS1ratio->GetXaxis()->SetTitleSize(0.05);
  hsThetaS1ratio->GetYaxis()->SetTitleSize(0.05);
  leg->Draw();
  c1_Logs15->Print(Form("%s/thetaS1ratioLog.png", saveDir));
  c1_Logs15->Print(Form("%s/thetaS1ratioLog.pdf", saveDir));
  c1_Logs15->Print(Form("%s/thetaS1ratioLog.tex", saveDir));
  TCanvas *c1_Logs16 = new TCanvas("c1_Logs16");
  c1_Logs16->SetLogy();
  hsPhiS1ratio->Draw("hist e nostack");
  hsPhiS1ratio->GetXaxis()->SetLabelSize(0.05); 
  hsPhiS1ratio->GetYaxis()->SetLabelSize(0.05);
  hsPhiS1ratio->GetXaxis()->SetTitleSize(0.05);
  hsPhiS1ratio->GetYaxis()->SetTitleSize(0.05);
  legTof->Draw();
  c1_Logs16->Print(Form("%s/phiS1ratioLog.png", saveDir));
  c1_Logs16->Print(Form("%s/phiS1ratioLog.pdf", saveDir));
  c1_Logs16->Print(Form("%s/phiS1ratioLog.tex", saveDir));

  TCanvas *c1_s11 = new TCanvas("c1_s11");
  hsThetaS1pro->Draw("hist e nostack");
  hsThetaS1pro->GetXaxis()->SetLabelSize(0.06);
  hsThetaS1pro->GetYaxis()->SetLabelSize(0.06);
  hsThetaS1pro->GetXaxis()->SetTitleSize(0.06);
  hsThetaS1pro->GetYaxis()->SetTitleSize(0.06);
  c1_s11->SetLeftMargin(0.13);
  c1_s11->SetBottomMargin(0.13);
  c1_s11->SetGridx();
  c1_s11->SetGridy();
  legTof->Draw();
  c1_s11->Print(Form("%s/thetaS1pro.png", saveDir));
  c1_s11->Print(Form("%s/thetaS1pro.pdf", saveDir));
  c1_s11->Print(Form("%s/thetaS1pro.tex", saveDir));
  TCanvas *c1_s12 = new TCanvas("c1_s12");
  hsThetaS1pi->Draw("hist e nostack");
  hsThetaS1pi->GetXaxis()->SetLabelSize(0.06);
  hsThetaS1pi->GetYaxis()->SetLabelSize(0.06);
  hsThetaS1pi->GetXaxis()->SetTitleSize(0.06);
  hsThetaS1pi->GetYaxis()->SetTitleSize(0.06);
  c1_s12->SetLeftMargin(0.13);
  c1_s12->SetBottomMargin(0.13);
  c1_s12->SetGridx();
  c1_s12->SetGridy();
  legTof->Draw();
  c1_s12->Print(Form("%s/thetaS1pi.png", saveDir));
  c1_s12->Print(Form("%s/thetaS1pi.pdf", saveDir));
  c1_s12->Print(Form("%s/thetaS1pi.tex", saveDir));
  TCanvas *c1_s13 = new TCanvas("c1_s13");
  hsPhiS1pro->Draw("hist e nostack");
  hsPhiS1pro->GetXaxis()->SetLabelSize(0.06); 
  hsPhiS1pro->GetYaxis()->SetLabelSize(0.06);
  hsPhiS1pro->GetXaxis()->SetTitleSize(0.06);
  hsPhiS1pro->GetYaxis()->SetTitleSize(0.06);
  c1_s13->SetGridx();
  c1_s13->SetGridy();
  leg->Draw();
  c1_s13->Print(Form("%s/phiS1pro.png", saveDir));
  c1_s13->Print(Form("%s/phiS1pro.pdf", saveDir));
  c1_s13->Print(Form("%s/phiS1pro.tex", saveDir));
  TCanvas *c1_s14 = new TCanvas("c1_s14");
  hsPhiS1pi->Draw("hist e nostack");
  hsPhiS1pi->GetXaxis()->SetLabelSize(0.06); 
  hsPhiS1pi->GetYaxis()->SetLabelSize(0.06);
  hsPhiS1pi->GetXaxis()->SetTitleSize(0.06);
  hsPhiS1pi->GetYaxis()->SetTitleSize(0.06);
  c1_s14->SetLeftMargin(0.13);
  c1_s14->SetBottomMargin(0.13);
  c1_s14->SetGridx();
  c1_s14->SetGridy();
  leg->Draw();
  c1_s14->Print(Form("%s/phiS1pi.png", saveDir));
  c1_s14->Print(Form("%s/phiS1pi.pdf", saveDir));
  c1_s14->Print(Form("%s/phiS1pi.tex", saveDir));
  TCanvas *c1_s15 = new TCanvas("c1_s15");
  hsThetaS1ratio->Draw("hist e nostack");
  hsThetaS1ratio->GetXaxis()->SetLabelSize(0.06);
  hsThetaS1ratio->GetYaxis()->SetLabelSize(0.06);
  hsThetaS1ratio->GetXaxis()->SetTitleSize(0.06);
  hsThetaS1ratio->GetYaxis()->SetTitleSize(0.06);
  c1_s15->SetLeftMargin(0.13);
  c1_s15->SetBottomMargin(0.13);
  c1_s15->SetGridx();
  c1_s15->SetGridy();
  legTheta->Draw();
  c1_s15->Print(Form("%s/thetaS1ratio.png", saveDir));
  c1_s15->Print(Form("%s/thetaS1ratio.pdf", saveDir));
  c1_s15->Print(Form("%s/thetaS1ratio.tex", saveDir));
  TCanvas *c1_s16 = new TCanvas("c1_s16");
  hsPhiS1ratio->Draw("hist e nostack");
  hsPhiS1ratio->GetXaxis()->SetLabelSize(0.06); 
  hsPhiS1ratio->GetYaxis()->SetLabelSize(0.06);
  hsPhiS1ratio->GetXaxis()->SetTitleSize(0.06);
  hsPhiS1ratio->GetYaxis()->SetTitleSize(0.06);
  c1_s16->SetLeftMargin(0.13);
  c1_s16->SetBottomMargin(0.13);
  c1_s16->Update();
  c1_s16->SetGridx();
  c1_s16->SetGridy();
  legPhiRatio->Draw();
  c1_s16->Print(Form("%s/phiS1ratio.png", saveDir));
  c1_s16->Print(Form("%s/phiS1ratio.pdf", saveDir));
  c1_s16->Print(Form("%s/phiS1ratio.tex", saveDir));

  TCanvas *cMomS1S2 = new TCanvas("cMomS1S2");
  hsMomS1S2->Draw("hist e nostack");
  hsMomS1S2->GetXaxis()->SetLabelSize(0.05);
  hsMomS1S2->GetYaxis()->SetLabelSize(0.05);
  hsMomS1S2->GetXaxis()->SetTitleSize(0.05);
  hsMomS1S2->GetYaxis()->SetTitleSize(0.05);
  leg->Draw();
  cMomS1S2->Print(Form("%s/proMomS1S2.png", saveDir));
  cMomS1S2->Print(Form("%s/proMomS1S2.pdf", saveDir));
  cMomS1S2->Print(Form("%s/proMomS1S2.tex", saveDir));

  TCanvas *cMomS1 = new TCanvas("cMomS1");
  hsMomS1->Draw("hist e nostack");
  cMomS1->SetGridx();
  cMomS1->SetGridy();
  hsMomS1->GetXaxis()->SetLabelSize(0.05);
  hsMomS1->GetYaxis()->SetLabelSize(0.05);
  hsMomS1->GetXaxis()->SetTitleSize(0.05);
  hsMomS1->GetYaxis()->SetTitleSize(0.05);
  cMomS1->SetLeftMargin(0.13);
  cMomS1->SetBottomMargin(0.13);
  cMomS1->Update();
  leg->Draw();
  cMomS1->Print(Form("%s/proMomS1.png", saveDir));
  cMomS1->Print(Form("%s/proMomS1.pdf", saveDir));
  cMomS1->Print(Form("%s/proMomS1.tex", saveDir));

  TCanvas *cutofS1S2 = new TCanvas("cutofS1S2");
  hsutof1dS1S2->Draw("hist e nostack");
  hsutof1dS1S2->GetXaxis()->SetLabelSize(0.05);
  hsutof1dS1S2->GetYaxis()->SetLabelSize(0.05);
  hsutof1dS1S2->GetXaxis()->SetTitleSize(0.05);
  hsutof1dS1S2->GetYaxis()->SetTitleSize(0.05);
  legTof->Draw();
  cutofS1S2->Print(Form("%s/utof1dS1S2.png", saveDir));
  cutofS1S2->Print(Form("%s/utof1dS1S2.pdf", saveDir));
  cutofS1S2->Print(Form("%s/utof1dS1S2.tex", saveDir));
  TCanvas *cutofS1 = new TCanvas("cutofS1");
  hsutof1dS1->Draw("hist e nostack");
  hsutof1dS1->GetXaxis()->SetLabelSize(0.05);
  hsutof1dS1->GetYaxis()->SetLabelSize(0.05);
  hsutof1dS1->GetXaxis()->SetTitleSize(0.05);
  hsutof1dS1->GetYaxis()->SetTitleSize(0.05);
  legTof->Draw();
  legTof->Write("legTof");
  cutofS1->Print(Form("%s/utof1dS1.png", saveDir));
  cutofS1->Print(Form("%s/utof1dS1.pdf", saveDir));
  cutofS1->Print(Form("%s/utof1dS1.tex", saveDir));

  TCanvas *cutofS1S2Log = new TCanvas("cutofS1S2Log");
  cutofS1S2Log->SetLogy();
  hsutof1dS1S2->Draw("hist nostack");
  cutofS1S2Log->SetGridx();
  cutofS1S2Log->SetGridy();
  hsutof1dS1S2->GetXaxis()->SetLabelSize(0.05);
  hsutof1dS1S2->GetYaxis()->SetLabelSize(0.05);
  hsutof1dS1S2->GetXaxis()->SetTitleSize(0.05);
  hsutof1dS1S2->GetYaxis()->SetTitleSize(0.05);
  cutofS1S2Log->SetLeftMargin(0.13);
  cutofS1S2Log->SetBottomMargin(0.13);
  legTof->Draw();
  cutofS1S2Log->Print(Form("%s/utof1dS1S2Log.png", saveDir));
  cutofS1S2Log->Print(Form("%s/utof1dS1S2Log.pdf", saveDir));
  cutofS1S2Log->Print(Form("%s/utof1dS1S2Log.tex", saveDir));
  TCanvas *cutofS1Log = new TCanvas("cutofS1Log");
  cutofS1Log->SetLogy();
  hsutof1dS1->Draw("hist nostack");
  cutofS1Log->SetGridx();
  cutofS1Log->SetGridy();
  hsutof1dS1->GetXaxis()->SetLabelSize(0.05);
  hsutof1dS1->GetYaxis()->SetLabelSize(0.05);
  hsutof1dS1->GetXaxis()->SetTitleSize(0.05);
  hsutof1dS1->GetYaxis()->SetTitleSize(0.05);
  cutofS1Log->SetLeftMargin(0.13);
  cutofS1Log->SetBottomMargin(0.13);
  legTof->Draw();
  cutofS1Log->Print(Form("%s/utof1dS1Log.png", saveDir));
  cutofS1Log->Print(Form("%s/utof1dS1Log.pdf", saveDir));
  cutofS1Log->Print(Form("%s/utof1dS1Log.tex", saveDir));
 
  fout->Close();
  delete fout;
} // angularDistS3
