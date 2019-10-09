// barCoins.C
// Big matrix of coincidences between various bars

void barCoins(const char* saveDir, 
	      const char* dstofDir="/nfs/scratch0/dbrailsf/data_backup/dtof_backup/",
	      const char* ustofDir="/nfs/scratch0/dbrailsf/data_backup/utof_backup_firsthitpinnedtounixtime/Data_root_v3_wo_walk_corr/")
{
  gROOT->SetBatch(kTRUE);
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
  // Define the runs to be used for varying number of blocks for ustof
  const char* str0Block = "Data_2018_8_31_b2_800MeV_0block.root";
  const char* str1Block = "Data_2018_9_1_b4_800MeV_1block_bend4cm.root";
  const char* str2Block = "Data_2018_9_1_b2_800MeV_2block_bend4cm.root";
  const char* str3Block = "Data_2018_9_1_b3_800MeV_3block_bend4cm.root";
  const char* str4Block = "Data_2018_9_1_b8_800MeV_4block_bend4cm.root";
  // New datasets to combine for the 4 block data
  const char* str4Block0 = "Data_2018_8_28_b5.root";
  const char* str4Block1 = "Data_2018_8_30_b1.root";
  const char* str4Block2 = "Data_2018_8_29_b4.root";
  const char* str4Block3 = "Data_2018_8_29_b1.root";
  std::vector<const char*> str4BlockVec = {/*str4Block0*/ str4Block1, str4Block2, str4Block3};

  const double startLongTime=1536919200;
  const double endLongTime  =1537005600;

  const double coinCut = 1.; // Cut for determining if a coincidence has occurred
  const double deadtimeCut = 185.;

  TFile *fout = new TFile(saveDir, "recreate");

  for (int sample=0; sample<6; sample++) {
    int totalTime = 0;
    string name;
    vector<double> startTimes;
    vector<double> endTimes;
    bool isMultiple = false;
    startTimes.clear();
    endTimes.clear();

    if (sample==0) {
      name="0";
      startTimes.push_back(start0Block);
      endTimes.push_back(end0Block);
    }
    else if (sample==1) {
      name="1";
      startTimes.push_back(start1Block);
      endTimes.push_back(end1Block);
    }
    else if (sample==2) {
      name="2";
      startTimes.push_back(start2Block);
      endTimes.push_back(end2Block);
    }
    else if (sample==3) {
      name="3";
      startTimes.push_back(start3Block);
      endTimes.push_back(end3Block);
    }
    else if (sample==4) {
      name="4";
      isMultiple = true;
      for (int b4=0; b4<str4BlockVec.size(); b4++) {
	TFile *futofTmp = new TFile(Form("%s/%s",ustofDir, str4BlockVec.at(b4)), "read");
	TTree *treeTmp = (TTree*)futofTmp->Get("tree");
	double tS1Tmp;
	treeTmp->SetBranchAddress("tS1", &tS1Tmp);
	TNamed *start = 0;
	TNamed *end   = 0;
	futofTmp->GetObject("start_of_run", start);
	const char* startchar = start->GetTitle();
	std::string startstr(startchar);
	std::string unixstart = startstr.substr(25,10);
	int startTime = stoi(unixstart);
	treeTmp->GetEntry(treeTmp->GetEntries() - 1);
	int endTime = startTime + (tS1Tmp/1e9);
	futofTmp->Close();
	delete futofTmp;
	startTimes.push_back(startTime);
	endTimes.push_back(endTime);
      }
    }
    else if (sample==5) {
      name="long";
      startTimes.push_back(startLongTime);
      endTimes.push_back(endLongTime);
    }

    TH2D *h2CosmicRate = new TH2D(Form("h2CosmicRate_%s", name.c_str()), "Cosmic rate in S4; Bar position / cm; Bar; Rate / s^{-1}", 20, 0., 140., 10, 0.5, 10.5);
    h2CosmicRate->GetXaxis()->SetTitleSize(.05);
    h2CosmicRate->GetXaxis()->SetLabelSize(.05);
    h2CosmicRate->GetYaxis()->SetTitleSize(.05);
    h2CosmicRate->GetYaxis()->SetLabelSize(.05);
    h2CosmicRate->GetZaxis()->SetTitleSize(.05);
    h2CosmicRate->GetZaxis()->SetLabelSize(.05);
    h2CosmicRate->Sumw2();
    TH2D *h2matrix = new TH2D(Form("matrix_%s", name.c_str()), "Bar coincidences; Bar 1; Bar 2; Rate / s^{-1}", 10, 0.5, 10.5, 10, 0.5, 10.5);
    h2matrix->GetXaxis()->SetTitleSize(.05);
    h2matrix->GetXaxis()->SetLabelSize(.05);
    h2matrix->GetYaxis()->SetTitleSize(.05);
    h2matrix->GetYaxis()->SetLabelSize(.05);
    h2matrix->GetZaxis()->SetTitleSize(.05);
    h2matrix->GetZaxis()->SetLabelSize(.05);

    // Loop over samples
    for (int n=0; n<startTimes.size(); n++) {
      int startTime = startTimes.at(n);
      int endTime   = endTimes.at(n);
      int nSpills = 0;
      double lastSpill = 0.;
      int runMin = -1;
      int runMax = -1;

      vector<double> deadtimeVec;
      deadtimeVec.resize(10, 0.);
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
      cout << "Min and max runs are " << runMin << " " << runMax << endl;

      TChain *tofCoinChain1 = new TChain("tofCoinTree");
      TChain *tofCoinChain2 = new TChain("tofCoinTree");
      for (int irun=runMin; irun<runMax+1; irun++) {
	// Load input files
	tofCoinChain1->Add(Form("%srun%d/DsTOFcoincidenceRun%d_tdc1.root", dstofDir, irun, irun));
	tofCoinChain2->Add(Form("%srun%d/DsTOFcoincidenceRun%d_tdc2.root", dstofDir, irun, irun));
      } // Loop over runs
      RawDsTofCoincidence *tofCoin1 = NULL;
      RawDsTofCoincidence *tofCoin2 = NULL;
      tofCoinChain1->SetBranchAddress("tofCoin", &tofCoin1);
      tofCoinChain2->SetBranchAddress("tofCoin", &tofCoin2);
      for (int h=0; h<tofCoinChain1->GetEntries(); h++) {
	tofCoinChain1->GetEntry(h);
	if (h % 100000 == 0) cout<<"Entry "<<h<<" of "<<tofCoinChain1->GetEntries()<<endl;

	if (tofCoin1->unixTime[0]<startTime) continue;
	if (tofCoin1->unixTime[0]>endTime) break;

	if (tofCoin1->lastDelayedBeamSignal != lastSpill) {
	  nSpills++;
	  lastSpill = tofCoin1->lastDelayedBeamSignal;
	}

	int bar1 = tofCoin1->bar;
	double dstofHitT1 = min(tofCoin1->fakeTimeNs[0], tofCoin1->fakeTimeNs[1]) - (10.-TMath::Abs(tofCoin1->fakeTimeNs[0]-tofCoin1->fakeTimeNs[1])/2.);
	if (!tofCoin1->inSpill && (dstofHitT1 - deadtimeVec.at(bar1-1)) > deadtimeCut) {
	  deadtimeVec.at(bar1-1) = dstofHitT1;
	  // Skip forward in this file to see if there any coincidences
	  for (int h2=h; h2<tofCoinChain1->GetEntries(); h2++) {
	    tofCoinChain1->GetEntry(h2);
	    double dstofHitT2 = min(tofCoin1->fakeTimeNs[0], tofCoin1->fakeTimeNs[1]) - (10.-TMath::Abs(tofCoin1->fakeTimeNs[0]-tofCoin1->fakeTimeNs[1])/2.);
	    if (dstofHitT2-dstofHitT1 > 20.) break;
	    int bar2 = tofCoin1->bar;
	    if (!tofCoin1->inSpill && (dstofHitT2 - dstofHitT1) < coinCut && 
		bar1 != bar2 && (dstofHitT2 - deadtimeVec.at(bar2-1)) > deadtimeCut) {
	      deadtimeVec.at(bar1-1) = dstofHitT2;
	      h2matrix->Fill(bar1, bar2);
	    } // Is not in a spill
	  } // Loop over TChain for a second time
	} // Is not in a spill
      } // Loop over entries

      deadtimeVec.clear();
      deadtimeVec.resize(10, 0.);

      for (int h=0; h<tofCoinChain2->GetEntries(); h++) {
	tofCoinChain2->GetEntry(h);
	if (h % 100000 == 0) cout<<"Entry "<<h<<" of "<<tofCoinChain2->GetEntries()<<endl;

	if (tofCoin2->unixTime[0]<startTime) continue;
	if (tofCoin2->unixTime[0]>endTime) break;

	int bar1 = tofCoin2->bar;
	double dstofHitT1 = min(tofCoin2->fakeTimeNs[0], tofCoin2->fakeTimeNs[1]) - (10.-TMath::Abs(tofCoin2->fakeTimeNs[0]-tofCoin2->fakeTimeNs[1])/2.);
	if (!tofCoin2->inSpill && (dstofHitT1 - deadtimeVec.at(bar1-1)) < deadtimeCut) {
	  deadtimeVec.at(bar1-1) = dstofHitT1;
	  // Skip forward in this file to see if there any coincidences
	  for (int h2=h; h2<tofCoinChain2->GetEntries(); h2++) {
	    tofCoinChain2->GetEntry(h2);
	    double dstofHitT2 = min(tofCoin2->fakeTimeNs[0], tofCoin2->fakeTimeNs[1]) - (10.-TMath::Abs(tofCoin2->fakeTimeNs[0]-tofCoin2->fakeTimeNs[1])/2.);
	    if (dstofHitT2-dstofHitT1 > 20.) break;
	    int bar2 = tofCoin2->bar;
	    if (!tofCoin2->inSpill && (dstofHitT2 - dstofHitT1) < coinCut && 
		bar1 != bar2 && (dstofHitT2 - deadtimeVec.at(bar2-1)) > deadtimeCut) {
	      deadtimeVec.at(bar1-1) = dstofHitT2;
	      h2matrix->Fill(bar1, bar2);
	    } // Is not in a spill
	  } // Loop over TChain for a second time
	} // Is not in a spill
      } // Loop over entries

      totalTime += endTime - startTime - nSpills;
    } // Loop over the sub samples within a sample

    h2matrix->Scale(1. / totalTime);
    fout->cd();
    h2matrix->Write();
  }

  fout->Close();
  delete fout;
} // barCoins
