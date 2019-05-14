// block4Times.C
// Beam data about the various runs used in the analysis

// Gets timestamp from input files
TTimestamp getStampFromLine(string line) {
  stringstream ssyear(line.substr(0,4));
  stringstream ssmonth(line.substr(5,2));
  stringstream ssday(line.substr(8,2));
  stringstream sshour(line.substr(11,2));
  stringstream ssmin(line.substr(14,2));
  stringstream sssec(line.substr(17,2));
  int year = 0;
  int month = 0;
  int day = 0;
  int hour = 0;
  int min = 0;
  int sec = 0;
  int count = 0;
  ssyear >> year;
  ssmonth >> month;
  ssday >> day;
  sshour >> hour;
  ssmin >> min;
  sssec >> sec;
  // Times are in CEST so require a 2 hour offset to get to UTC
  TTimeStamp stamp(year, month, day, hour, min, sec, 0, "kTRUE", -7200);
  return stamp;
}

void datasetInfo(const char* saveDir,
		 const char* ustofDir="/zfs_home/sjones/mylinktoutof/",
		 const char* dstofDir="/scratch0/dbrailsf/temp/mylinktodtof/",
		 const char* spillDir="/scratch2/sjones/spillDB/",
		 const char* scintFile="/scratch2/sjones/DsTOFanalysis/files/Scintillator_info.txt") 
{
  TFile *fout = new TFile(Form("%s/block4Times_out.root", saveDir), "recreate");
  gROOT->SetBatch(kTRUE);
  // Get the scintillator values from the text file and put into a vector
  std::vector<int> countVec;
  std::vector<TTimeStamp> stampVec;
  string line;
  ifstream inFile;
  inFile.open(scintFile);
  while(!inFile.eof()) {
    getline(inFile, line);
    if (line.length()>20) {
      stringstream sscount(line.substr(24));
      int count = 0;
      sscount >> count;
      TTimestamp = getStampFromLine(line);
      cout<<stamp.GetSec()<<endl;
      countVec.push_back(count);
      stampVec.push_back(stamp);
    } // if (line.length()>20)
  } // while(!inFile.eof())
  // Do same thing for BZH01 (beam momentum magnet) and BZH03 (bending magnet) values
  cout<<"Scintillator file has "<<countVec.size()<<" spills recorded"<<endl;

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
  TMultiGraph *mgScint = new TMultiGraph();
  TMultiGraph *mgScintAvg = new TMultiGraph();
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
    TTimeStamp utofStartStamp(startTimeUtof);
    TTimeStamp utofEndStamp(endTimeUtof);

    TGraph *grUtof = new TGraph();
    TGraph *grScint = new TGraph();
    TGraph *grScintAvg = new TGraph();
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
      grScint->SetMarkerColor(kBlack);
      grScintAvg->SetMarkerColor(kBlack);
      grScintAvg->SetMarkerStyle(3);
    }
    else if (b==1) {
      grUtof->SetMarkerColor(kRed);
      grUtof->SetMarkerStyle(47);
      grScint->SetMarkerColor(kRed);
      grScintAvg->SetMarkerColor(kRed);
      grScintAvg->SetMarkerStyle(3);
    }
    else if (b==2) {
      grUtof->SetMarkerColor(kBlue);
      grUtof->SetMarkerStyle(20);
      grScint->SetMarkerColor(kBlue);
      grScintAvg->SetMarkerColor(kBlue);
      grScintAvg->SetMarkerStyle(3);
    }
    else if (b==3) {
      grUtof->SetMarkerColor(kCyan+1);
      grUtof->SetMarkerStyle(22);
      grScint->SetMarkerColor(kCyan+1);
      grScintAvg->SetMarkerColor(kCyan+1);
      grScintAvg->SetMarkerStyle(3);
    }
    else if (b==4) {
      grUtof->SetMarkerColor(kGreen+2);
      grUtof->SetMarkerStyle(34);
      grScint->SetMarkerColor(kGreen+2);
      grScintAvg->SetMarkerColor(kGreen+1);
      grScintAvg->SetMarkerStyle(3);
    }
    else if (b==5) {
      grUtof->SetMarkerColor(kOrange+1);
      grUtof->SetMarkerStyle(49);
      grScint->SetMarkerColor(kOrange+1);
      grScintAvg->SetMarkerColor(kOrange+1);
      grScintAvg->SetMarkerStyle(3);
    }
    else {
      grUtof->SetMarkerColor(kMagenta+1);
      grUtof->SetMarkerStyle(39);
      grScint->SetMarkerColor(kMagenta+1);
      grScintAvg->SetMarkerColor(kMagenta+1);
      grScintAvg->SetMarkerStyle(3);
    }
    
    //grUtof->SetMarkerStyle(20);
    grUtof->SetMarkerSize(2);
    mg->Add(grUtof);
    futof->Close();
    delete futof;
    fout->cd();
    grUtof->Write(Form("grUtof_%s", blockVec[b]));

    int nScint = 1;
    int nCount = 1;
    // Find appropriate scintillator values and put into graph
    for (int s=0; s<stampVec.size(); ++s) {
      if (stampVec[s] > utofEndStamp) break;
      if (stampVec[s] < utofStartStamp) continue;

      if (stampVec[s] > utofStartStamp && stampVec[s] < utofEndStamp) {
	grScint->SetPoint(grScint->GetN(), stampVec[s], countVec[s]);
	nScint++;
	nCount+=countVec[s];
	if (nScint % 20 == 0) {
	  grScintAvg->SetPoint(grScintAvg->GetN(), stampVec[s], (double)nCount/20.); 
	} // if (nScint % 20 == 0)
      } // if (stampVec[s] > utofStartStamp && stampVec[s] < utofEndStamp)
    } // for (int s=0; s<stampVec.size(); ++s)
    mg->Add(grScint);
    mgScint->Add(grScint);
    mgScintAvg->Add(grScintAvg);
    fout->cd();
    grScint->Write(Form("grScint_%s", blockVec[b]));
    grScintAvg->Write(Form("grScintAvg_%s", blockVec[b]));

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
    int nDtof = 1;
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
	nDtof++;

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

  fout->cd();
  mg->Write("mg");
  mgScint->Write("mgScint");
  mgScintAvg->Write("mgScintAvg");

  fout->Close();
  delete fout;

} // block4Times
