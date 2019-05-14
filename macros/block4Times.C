// block4Times.C

void block4Times(const char* saveDir,
		 const char* ustofDir="/zfs_home/sjones/mylinktoutof/",
		 const char* dstofDir="/scratch0/dbrailsf/temp/mylinktodtof/",
		 const char* spillDir="/scratch2/sjones/spillDB/") 
{
  //  gROOT->SetBatch(kTRUE);
  // For each set of data get the start and end times and put in an individual graph
  const char* str4Block1 = "Data_2018_9_1_b8_800MeV_4block_bend4cm.root";
  const char* str4Block2 = "Data_2018_9_3_b2_800MeV_4block_bend4cm.root";
  const char* str4Block3 = "Data_2018_9_1_b5_800MeV_4block.root";
  const char* str4Block4 = "Data_2018_9_4_b1_800MeV_4block.root";
  const char* str4Block5 = "Data_2018_9_4_b2_800MeV_4block.root";
  const char* str4Block6 = "Data_2018_9_3_b4_800MeV_4block_1mUSB.root";
  const char* str0Block = "Data_2018_8_31_b2_800MeV_0block.root";
  const char* str1Block = "Data_2018_9_1_b4_800MeV_1block_bend4cm.root";
  const char* str2Block = "Data_2018_9_1_b2_800MeV_2block_bend4cm.root";
  const char* str3Block = "Data_2018_9_1_b3_800MeV_3block_bend4cm.root";
  std::vector<const char*> blockVec = {str4Block1, str4Block2, str4Block3,
				       str4Block4, str4Block5, str4Block6,
				       str0Block, str1Block, str2Block, str3Block};
  TMultiGraph *mg = new TMultiGraph();
  TGraph *grDtof = new TGraph();
  for (int b=0; b<blockVec.size(); b++) {
    TFile *futof = new TFile(Form("%s/%s",ustofDir,blockVec[b]), "read");
    TTree *tree = (TTree*)futof->Get("tree");
    double tS1;
    tree->SetBranchAddress("tS1", &tS1);
    TNamed *start = 0;
    TNamed *end   = 0;
    futof->GetObject("start_of_run", start);
    futof->GetObject("end_of_run", end);
    const char* startchar = start->GetTitle();
    std::string startstr(startchar);
    std::string unixstart = startstr.substr(25,10);
    int startTimeUtof = stoi(unixstart);
    tree->GetEntry(tree->GetEntries()-1);
    int endTimeUtof = (tS1/1e9) + startTimeUtof;

    TGraph *grUtof = new TGraph();
    if (b < 6) {
      grUtof->SetPoint(grUtof->GetN(), startTimeUtof, 3500);
      grUtof->SetPoint(grUtof->GetN(), endTimeUtof, 3500);
    }
    else {
      grUtof->SetPoint(grUtof->GetN(), startTimeUtof, 3000);
      grUtof->SetPoint(grUtof->GetN(), endTimeUtof, 3000);
    }

    if (b==0) {
      grUtof->SetMarkerColor(kBlack);
      grUtof->SetMarkerStyle(21);
    }
    else if (b==1) {
      grUtof->SetMarkerColor(kRed);
      grUtof->SetMarkerStyle(47);
    }
    else if (b==2) {
      grUtof->SetMarkerColor(kBlue);
      grUtof->SetMarkerStyle(20);
    }
    else if (b==3) {
      grUtof->SetMarkerColor(kCyan+1);
      grUtof->SetMarkerStyle(22);
    }
    else if (b==4) {
      grUtof->SetMarkerColor(kGreen+2);
      grUtof->SetMarkerStyle(34);
    }
    else if (b==5) {
      grUtof->SetMarkerColor(kOrange+1);
      grUtof->SetMarkerStyle(49);
    }
    else if (b==6) {
      grUtof->SetMarkerColor(kMagenta+1);
      grUtof->SetMarkerStyle(39);
    }
    else if (b==7) {
      grUtof->SetMarkerColor(kMagenta+1);
      grUtof->SetMarkerStyle(39);
    }
    else if (b==8) {
      grUtof->SetMarkerColor(kMagenta+1);
      grUtof->SetMarkerStyle(39);
    }
    else if (b==9) {
      grUtof->SetMarkerColor(kMagenta+1);
      grUtof->SetMarkerStyle(39);
    }
    //grUtof->SetMarkerStyle(20);
    grUtof->SetMarkerSize(2);
    mg->Add(grUtof);
    futof->Close();
    delete futof;

    // Find dtof runs
    int runMin = -1;
    int runMax = -1;
    int startTime = startTimeUtof;
    int endTime   = endTimeUtof;
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
    for (int irun = runMin; irun < runMax+1; irun++) {
      TFile *dbFile = new TFile(Form("%s/spillDB_run%d_run%d.root", spillDir, irun, irun), "read");
      TTree *spillTree = (TTree*)dbFile->Get("spillTree");
      double globalSpillTime;
      double ustofSpillTime;
      spillTree->SetBranchAddress("globalSpillTime", &globalSpillTime);
      spillTree->SetBranchAddress("ustofSpillTime", &ustofSpillTime);
      TFile *fdtof = new TFile(Form("%s/run%d/DsTOFtreeRun%d_tdc1.root", dstofDir, irun, irun), "read");
      RawDsTofHeader *tof = NULL;
      TTree *tofTree = (TTree*)fdtof->Get("tofTree");
      tofTree->SetBranchAddress("tof", &tof);
      int lasts = 0;
      tofTree->GetEntry(0);
      int firstTime = tof->unixTime;
      double lastS1S2Dtof = 0.;
      for (int t = 0; t < spillTree->GetEntries(); t++) {
	spillTree->GetEntry(t);
	int hits = 0;
	if (globalSpillTime >= startTime && globalSpillTime <= endTime) {
	  for (int s = lasts; s<tofTree->GetEntries(); s++) {
	    tofTree->GetEntry(s);
	    if ((tof->fakeTimeNs/1e9)+firstTime >= globalSpillTime+1.) break;
	    if ((tof->fakeTimeNs/1e9)+firstTime <= globalSpillTime) continue;

	    if ((tof->fakeTimeNs/1e9)+firstTime >= globalSpillTime &&
		(tof->fakeTimeNs/1e9)+firstTime <= globalSpillTime + 1. && 
		tof->channel == 13 && (tof->fakeTimeNs - lastS1S2Dtof) > 500.) {
	      lastS1S2Dtof = tof->fakeTimeNs;
	      hits++;
	      lasts = s;
	    }
	  } // for (int s = lasts; s<tofTree->GetEntries(); s++)
	} // if (globalSpillTime >= startTime && globalSpillTime <= endTime)
	grDtof->SetPoint(grDtof->GetN(), globalSpillTime, hits);
      } // for (int t = 0; t < spillTree->GetEntries(); t++)
      fdtof->Close();
      delete fdtof;
      dbFile->Close();
      delete dbFile;
    } // for (int irun = runMin; irun < runMax+1; irun++) 
  } // for (int b=0; b<blockVec.size(); b++)
  mg->Add(grDtof);
  mg->SetTitle("S1 #cap S2 dtof hits for selected 4 block utof runs; Time / s; Hits / spill");
  TCanvas *c1 = new TCanvas("c1");
  mg->Draw("AP");

  TFile *fout = new TFile(Form("%s/block4Times_out.root", saveDir), "recreate");
  fout->cd();
  mg->Write("mg");

  fout->Close();
  delete fout;

} // block4Times
