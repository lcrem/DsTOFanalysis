// deadtimeCorrection.C
#include "UsefulFunctions.C"

const double ustofDelay = 184.7;
const double cutOffT    = 0.29;

void deadtimeCorrection(const char* outfile)
{
  gROOT->SetBatch(kTRUE);

  TFile *fout = new TFile(outfile, "recreate");

  // Loop over samples
  for (int b=0; b<5; b++) {
    cout<<"========================="<<endl;
    cout<<b<<" blocks"<<endl;
    cout<<"========================="<<endl;

    vector<const char*> utofFiles;
    vector<double> startTimes;
    vector<double> endTimes;

    if (b == 0) {
      startTimes.push_back(start0Block);
      endTimes.push_back(end0Block);
      utofFiles.push_back(str0Block);
    }
    else if (b == 1) {
      startTimes.push_back(start1Block);
      endTimes.push_back(end1Block);
      utofFiles.push_back(str1Block);
    }
    else if (b == 2) {
      startTimes.push_back(start2Block);
      endTimes.push_back(end2Block);
      utofFiles.push_back(str2Block);
    }
    else if (b == 3) {
      startTimes.push_back(start3Block);
      endTimes.push_back(end3Block);
      utofFiles.push_back(str3Block);
    }
    else if (b == 4) {
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
      utofFiles = str4BlockVec;
    }

    TTree *hitTreeD = new TTree(Form("hitTreeD%d", b), "Spill hits");
    int dtofS1S2;
    double spillTimeD;
    hitTreeD->Branch("dtofS1S2", &dtofS1S2);
    hitTreeD->Branch("spillTime", &spillTimeD);

    TTree *hitTreeU = new TTree(Form("hitTreeU%d", b), "Spill hits");
    int utofS1S2;
    double spillTimeU;
    hitTreeU->Branch("utofS1S2", &utofS1S2);
    hitTreeU->Branch("spillTime", &spillTimeU);

    // Loop over subsamples
    for (int sub=0; sub<startTimes.size(); sub++) {
      double startTime = startTimes.at(sub);
      double endTime   = endTimes.at(sub);

      vector<double> utofTimes;
      vector<double> dtofTimes;

      // Find correct dtof files
      Int_t runMin = -1;
      Int_t runMax = -1;
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

      // Now loop through dtof runs and get the hits
      for (int irun=runMin; irun<runMax+1; irun++) {
	TFile *dbFile = new TFile(Form("/scratch0/sjones/spillDB/spillDB_run%d_run%d.root", irun, irun), "read");
	TFile *rawFile = new TFile(Form("%srun%d/DsTOFtreeRun%d_tdc1.root", dstofDir, irun, irun), "read");
	RawDsTofHeader *tof = NULL;
	TTree *tofTree = (TTree*)rawFile->Get("tofTree");
	tofTree->SetBranchAddress("tof", &tof);
	tofTree->GetEntry(0);
	double fileStart = tof->unixTime;
	TTree *spillTree = (TTree*)dbFile->Get("spillTree");
	double globalSpillTime;
	double ustofSpillTime;
	spillTree->SetBranchAddress("globalSpillTime", &globalSpillTime);
	spillTree->SetBranchAddress("ustofSpillTime", &ustofSpillTime);
	int laste = 0;
	// Loop over spills
	for (int t=0; t<spillTree->GetEntries(); t++) {
	  spillTree->GetEntry(t);
	  dtofS1S2 = 0;
	  spillTimeD = globalSpillTime;
	  if (globalSpillTime<startTime) continue;
	  if (globalSpillTime>endTime) break;
	  utofTimes.push_back(ustofSpillTime);
	  dtofTimes.push_back(globalSpillTime);
	  cout.precision(19);
	  cout<<"Spill time = "<<spillTimeD<<endl;
	  // Loop over tof entries
	  for (int e=laste; e<tofTree->GetEntries(); e++) {
	    tofTree->GetEntry(e);
	    double hitTime = tof->fakeTimeNs/1e9 + fileStart;
	    if (hitTime < globalSpillTime) continue;
	    if (hitTime > globalSpillTime + 1.) {
	      laste = e;
	      break;
	    }
	    // laste = e;
	    if (tof->channel == 13) { // Is an S1 S2 coincidence
	      dtofS1S2++;
	    } 
	  } // Loop over raw tof entries
	  fout->cd();
	  hitTreeD->Fill();
	} // Loop over spills
	dbFile->Close();
	delete dbFile;
	rawFile->Close();
	delete rawFile;
      } // Loop over runs

      // Now do a similar sort of process for the UToF

    } // Loop over subsamples
    fout->cd();
    hitTreeD->Write();
  } // Loop over samples

  fout->Close();
  delete fout;
} // deadtimeCorrection
