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

    // Loop over subsamples
    for (int sub=0; sub<startTimes.size(); sub++) {
      TTree *hitTree = new TTree(Form("hitTree%d_%d", b, sub), "Spill hits");
      int dtofS1S2, utofS1S2;
      int utofS1;
      int isGood;
      double spillTimeD, spillTimeU;
      hitTree->Branch("spillTimeD", &spillTimeD);
      hitTree->Branch("spillTimeU", &spillTimeU);
      hitTree->Branch("s1s2D", &dtofS1S2);
      hitTree->Branch("s1s2U", &utofS1S2);
      hitTree->Branch("s1U", &utofS1);
      hitTree->Branch("isGood", &isGood);

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

	// DToF files
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

	// UToF files
	TFile *utofFile = new TFile(Form("%s/%s", ustofDir, utofFiles.at(sub)), "read");
	TTree *utofTree = (TTree*)utofFile->Get("tree");
	double tS1, tTrig;
	utofTree->SetBranchAddress("tS1", &tS1);
	utofTree->SetBranchAddress("tTrig", &tTrig);
	TNamed *start = 0;
	TNamed *end   = 0;                                                                          
	utofFile->GetObject("start_of_run", start);
	const char* startchar = start->GetTitle();
	std::string startstr(startchar);
	std::string unixstart = startstr.substr(25,10);
	int utofFileStart = stoi(unixstart);
	utofTree->GetEntry(utofTree->GetEntries() - 1);
	double utofFileEnd = utofFileStart + (tS1/1e9);

	int laste = 0;
	int lastu = 0;
	// Loop over spills
	for (int t=0; t<spillTree->GetEntries(); t++) {
	  spillTree->GetEntry(t);
	  dtofS1S2 = 0;
	  spillTimeD = globalSpillTime;
	  spillTimeU = ustofSpillTime;
	  if (globalSpillTime<startTime) continue;
	  if (globalSpillTime>endTime) break;
	  utofTimes.push_back(ustofSpillTime);
	  dtofTimes.push_back(globalSpillTime);
	  cout.precision(19);
	  // cout<<"Spill time = "<<spillTimeD<<endl;
	  // Loop over tof entries
	  for (int e=laste; e<tofTree->GetEntries(); e++) {
	    tofTree->GetEntry(e);
	    double hitTime = tof->fakeTimeNs/1e9 + fileStart;
	    if (hitTime < globalSpillTime) continue;
	    if (hitTime > globalSpillTime + 1.) {
	      laste = e;
	      break;
	    }

	    if (tof->channel == 13) { // Is an S1 S2 coincidence
	      dtofS1S2++;
	    } 
	  } // Loop over raw tof entries

	  // if (ustofSpillTime < utofFileStart) continue;
	  // if (ustofSpillTime > utofFileEnd) break;
	  utofS1S2 = 0;
	  utofS1 = 0;
	  isGood = 0;
	  // Now loop through UToF file and do something similar
	  for (int e=lastu; e<utofTree->GetEntries(); e++) {
	    utofTree->GetEntry(e);
	    double hitTime = tS1/1e9 + utofFileStart;
	    if (hitTime < ustofSpillTime) continue;
	    if (hitTime > ustofSpillTime + 1.) {
	      lastu = e;
	      break;
	    }
	    utofS1++;
	    if (tTrig != 0) utofS1S2++;
	    if ((tTrig/1e9) + (double)utofFileStart > ustofSpillTime + 0.47
		&& tTrig != 0) isGood = 1;
	  } // Loop over utof tree entries

	  fout->cd(); 
	  hitTree->Fill();
	} // Loop over spills
	dbFile->Close();
	delete dbFile;
	rawFile->Close();
	delete rawFile;
	utofFile->Close();
	delete utofFile;
      } // Loop over runs

      fout->cd();
      // Now loop through the hit tree and make some TGraphs 
      TGraphErrors *gS1S2 = new TGraphErrors();
      gS1S2->SetTitle("Utof S1S2 vs. Dtof S1S2; DToF S1S2; DToF S1S2; UToF S1S2");
      TGraphErrors *gS1S2Rat = new TGraphErrors();
      gS1S2Rat->SetTitle("Utof S1S2 vs. Dtof S1S2; DToF S1S2; (UToF S1S2)/(DToF S1S2)");
      TGraphErrors *gS1U  = new TGraphErrors();
      gS1U->SetTitle("Utof S1 vs. Dtof S1S2; DToF S1S2; UToF S1");
      int p=0;
      for (int t=0; t<hitTree->GetEntries(); t++) {
	hitTree->GetEntry(t);
	if (isGood == 1 && utofS1S2>10) {
	  double eUtofS1S2 = sqrt(utofS1S2);
	  double eDtofS1S2 = sqrt(dtofS1S2);
	  double eUtofS1   = sqrt(utofS1);
	  double rS1S2  = (double)utofS1S2 / (double)dtofS1S2;
	  double erS1S2 = rS1S2 * sqrt( pow(eUtofS1S2/(double)utofS1S2, 2) + pow(eDtofS1S2/(double)dtofS1S2, 2) );
	  gS1S2->SetPoint(p, dtofS1S2, utofS1S2);
	  gS1S2->SetPointError(p, eDtofS1S2, eUtofS1S2);
	  gS1S2Rat->SetPoint(p, dtofS1S2, rS1S2);
	  gS1S2Rat->SetPointError(p, eDtofS1S2, erS1S2);
	  gS1U->SetPoint(p, dtofS1S2, utofS1);
	  gS1U->SetPointError(p, eDtofS1S2, eUtofS1);
	  p++;
	}
      }
      TF1 *f0 = new TF1(Form("f0_%d_%d", b, sub), "pol0", 200, 2200);
      TF1 *f1 = new TF1(Form("f1_%d_%d", b, sub), "pol1", 200, 2200);
      TF1 *f2 = new TF1(Form("f2_%d_%d", b, sub), "pol2", 200, 2200);
      f1->SetParameter(0, 0.5);
      f1->SetParameter(1, -0.0003);
      f2->SetParameter(0, 0.5);
      f2->SetParameter(1, -0.0003);
      f2->SetParameter(2, 1e-7);
      gS1S2->Write(Form("gS1S2_%d_%d", b, sub));
      gS1S2Rat->Write(Form("gS1S2Rat_%d_%d", b, sub));
      gS1U->Write(Form("gS1U_%d_%d", b, sub));
      gS1S2Rat->Fit(f0, "r");
      gS1S2Rat->Fit(f1, "r");
      gS1S2Rat->Fit(f2, "r");
      f0->Write();
      f1->Write();
      f2->Write();
    hitTree->Write();
    } // Loop over subsamples
  } // Loop over samples

  fout->Close();
  delete fout;
} // deadtimeCorrection
