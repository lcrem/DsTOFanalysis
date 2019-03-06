// deadtime.C

void deadtime(const char* saveDir,
	      const char* dstofDir="/scratch0/dbrailsf/temp/mylinktodtof",
	      const char* ustofDir="/zfs_home/sjones/mylinktoutof") 
{
  gROOT->SetBatch(kTRUE);
  // Define the runs to be used for varying number of blocks
  const char* str0Block = "Data_2018_8_31_b2_800MeV_0block.root";
  const char* str1Block = "Data_2018_9_1_b4_800MeV_1block_bend4cm.root";
  const char* str2Block = "Data_2018_9_1_b2_800MeV_2block_bend4cm.root";
  const char* str3Block = "Data_2018_9_1_b3_800MeV_3block_bend4cm.root";
  const char* str4Block = "Data_2018_9_1_b8_800MeV_4block_bend4cm.root";
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
  // Ustof-dstof cable delay
  const double ustofDelay = 184.7;

  for (int nBlocks = 0; nBlocks <=4; nBlocks++) {
    cout<<nBlocks<<" blocks"<<endl;

    Int_t runMin=-1;
    Int_t runMax=-1;

    double startTime = 0;
    double endTime   = 0;

    const char* nustof;
    if (nBlocks == 0) {
      startTime = start0Block;
      endTime   = end0Block;
      nustof = Form("%s/%s", ustofDir, str0Block);
    }
    else if (nBlocks == 1) {
      startTime = start1Block;
      endTime   = end1Block;
      nustof = Form("%s/%s", ustofDir, str1Block);
    }
    else if (nBlocks == 2) {
      startTime = start2Block;
      endTime   = end2Block;
      nustof = Form("%s/%s", ustofDir, str2Block);
    }
    else if (nBlocks == 3) {
      startTime = start3Block;
      endTime   = end3Block;
      nustof = Form("%s/%s", ustofDir, str3Block);
    }
    else if (nBlocks == 4) {
      startTime = start4Block;
      endTime   = end4Block;
      nustof = Form("%s/%s", ustofDir, str4Block);
    }

    // Go through utof files and count the number of hits
    TFile *futof = new TFile(nustof, "read");

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

    int nS1S2utof = 0;
    int nS1S2dtof = 0;
    // Count number of S1 x S2 hits
    for (int t=0; t<tree->GetEntries(); t++ ) {
      tree->GetEntry(t);
      if (tTrig !=0) {
	nS1S2utof++;
      }
    } // for (int t=0; t<tree->GetEntries(); t++ )

    // Find the correct dtof files
    for (int irun=1000; irun<1100; irun++) {
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

    // Now go over dtof files and count the number of s1s2 coincidences
    for (int irun = runMin; irun < runMax+1; irun++) {
      // Utof and beam signals go into both TDCs so only need to run over one
      TFile *tofFile = new TFile(Form("%s/run%d/DsTOFtreeRun%d_tdc1.root", dstofDir, irun, irun));
      RawDsTofHeader *tof = NULL;
      TTree *tofTree = (TTree*)tofFile->Get("tofTree");
      tofTree->SetBranchAddress("tof", &tof);
      double tempBeam = 0.;
      double beamNs = 0.;
      double lastRawBeamSpillNs = 0.;
      TTree *beamTree = new TTree("beamTree", "beam");
      beamTree->SetDirectory(0);
      beamTree->Branch("beamNs", &beamNs, "beamNs/D");
      for (int i=0; i<tofTree->GetEntries(); i++) {
	tofTree->GetEntry(i);
	if (tof->unixTime < startTime) continue;
	if (tof->unixTime > endTime) break;
	
	if (tof->channel == 15) {
	  lastRawBeamSpillNs = tof->fakeTimeNs;
	  continue;
	} // if (tof->channel == 15) 

	if (tof->channel == 14) {
	  if (tof->fakeTimeNs>(lastRawBeamSpillNs+900.98e6) && tof->fakeTimeNs<(lastRawBeamSpillNs+900.982e6)) {
	    beamNs = tof->fakeTimeNs - ustofDelay;
	    tempBeam = beamNs;
	    beamTree->Fill();
	    cout.precision(12);

	  }
	}
      } // for (int i=0; i<tof->GetEntries(); i++)

      beamTree->BuildIndex("beamNs");
      cout<<"Beam spill tree has "<<beamTree->GetEntries()<<" entries"<<endl;
      // Now loop again and this time count all the S1S2 coincidences in the spill windows
      for (int i=0; i<tofTree->GetEntries(); i++) {
	tofTree->GetEntry(i);
	if (tof->unixTime < startTime) continue;
	if (tof->unixTime > endTime) break;
	if (tof->channel == 13) {
	  double ustofTemp = tof->fakeTimeNs;
	  int beamEntry = beamTree->GetEntryNumberWithBestIndex(tof->fakeTimeNs);
	  beamTree->GetEntry(beamEntry);
	  // Is in a spill
	  if ((ustofTemp - beamNs) > 0. && (ustofTemp - beamNs) < 1e9) {
	    nS1S2dtof++;
	  }
	}
      } // for (int i=0; i<tof->GetEntries(); i++)
      delete tof;
      tofFile->Close();
      //      delete ustofTree;
      delete beamTree;
    } // for (int irun = runMin; irun < runMax+1; irun++)

    cout<<"Dtof S1 S2 "<<nS1S2dtof<<", Utof S1 S2 "<<nS1S2utof<<endl;
  } // for (int nBlocks = 0; nBlocks <=4; nBlocks++)
   
  
} // deadtime
