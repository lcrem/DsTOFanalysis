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
      double lastS1S2 = 0.;
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
	  // Is in a spill && is separated by more than 400 ns
	  if ((ustofTemp - beamNs) > 0. && (ustofTemp - beamNs) < 1e9 && (ustofTemp - lastS1S2) > 500.) {
	    nS1S2dtof++;
	    lastS1S2 = ustofTemp;
	  }
	}
      } // for (int i=0; i<tof->GetEntries(); i++)
      delete tof;
      tofFile->Close();
      //      delete ustofTree;
      delete beamTree;
    } // for (int irun = runMin; irun < runMax+1; irun++)

    cout<<"Dtof S1 S2 "<<nS1S2dtof<<", Utof S1 S2 "<<nS1S2utof<<", Ratio "<<(double)nS1S2utof/(double)nS1S2dtof<<endl;
  } // for (int nBlocks = 0; nBlocks <=4; nBlocks++)
  
} // deadtime

// Same thing but do it for a given utof file with the appropriate matching dtof file
void deadtimeTimestamp(const char* utofFile,	    
		       const char* saveDir,
		       //const char* spillDB,
		       const char* dstofDir="/scratch0/dbrailsf/temp/mylinktodtof",
		       const char* ustofDir="/zfs_home/sjones/mylinktoutof") 
{
  gROOT->SetBatch(kTRUE);

  TH1D *hDeltatDtof = new TH1D(Form("hDeltaDtof_%s", utofFile), Form("S1 #cap S2 #deltat as measured in dtof: %s", utofFile), 300, 500, 1000000);
  TH1D *hDeltatUtof = new TH1D(Form("hDeltaUtof_%s", utofFile), Form("S1 #cap S2 #deltat as measured in utof: %s", utofFile), 300, 500, 1000000);
  TH1D *hUtofInSpill = new TH1D(Form("hUtofInSpill_%s", utofFile), Form("S1 #cap S2 time since start of spill: %s", utofFile), 400, 0.1, .7);

  // Ustof-dstof cable delay
  const double ustofDelay = 184.7;

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
  for (int irun=1000; irun<1400; irun++) {
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
  
  std::vector<int> nS1S2dtofVec;
  // Now go over the spill DB files for these runs and 
  // make two vectors of the appropriate matching spill times
  for (int irun = runMin; irun < runMax+1; irun++) {
    cout<<"Getting hits for dtof run "<<irun<<endl;
    // if (irun == 1058) continue;
    TFile *dbFile = new TFile(Form("~/spillDB/spillDB_run%d_run%d.root", irun, irun), "read");
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
      tempDtofTimes.push_back(globalSpillTime);
    } // for (int t = 0; t < spillTree->GetEntries(); t++)
    dbFile->Close();
    // Now open the corresponding dtof file and count the appropriate spills
    // Utof and beam signals go into both TDCs so only need to run over one
    TFile *tofFile = new TFile(Form("%s/run%d/DsTOFtreeRun%d_tdc1.root", dstofDir, irun, irun));
    RawDsTofHeader *tof = NULL;
    TTree *tofTree = (TTree*)tofFile->Get("tofTree");
    tofTree->SetBranchAddress("tof", &tof);
    tofTree->GetEntry(0);
    UInt_t firstTime = tof->unixTime;
    vector<int> nS1S2dtofTemp;
    nS1S2dtofTemp.resize(tempDtofTimes.size(), 0);
    // Loop over only the spill times in this particular dtof file
    // Saves us having to open loads of dtof files to check if the spill is there
    double lastS1S2Dtof = 0.;
    int lastdt = 0;
    for (int spill = 0; spill < tempDtofTimes.size(); spill++) {
      for (int t = lastdt; t < tofTree->GetEntries(); t++) {
	tofTree->GetEntry(t);
	if ((tof->fakeTimeNs/1e9 - tempDtofTimes[spill] + firstTime) > 1.) break;
	if ((tof->fakeTimeNs/1e9 - tempDtofTimes[spill] + firstTime) < 0.) continue;
	// Check that hit is up to 1s after the chosen spill
	if (tof->channel == 13 && 
	    (tof->fakeTimeNs/1e9 - tempDtofTimes[spill] + firstTime) < 1. &&
	    (tof->fakeTimeNs/1e9 - tempDtofTimes[spill] + firstTime) > 0. &&
	    (tof->fakeTimeNs - lastS1S2Dtof) > 500.) {
	  nS1S2dtofTemp[spill]++;
	  hDeltatDtof->Fill(tof->fakeTimeNs - lastS1S2Dtof);
	  lastS1S2Dtof = tof->fakeTimeNs;
	  lastdt = t;
	}
      } // for (int t = 0; t < tofTree->GetEntries(); t++)
    } // for (int spill = 0; spill < tempDtofTimes.size(); spill++)
    // Add this to the larger spill vector 

    for (int i = 0; i < nS1S2dtofTemp.size(); i++) {
      nS1S2dtofVec.push_back(nS1S2dtofTemp[i]);
    }

    delete tof;
    tofFile->Close();
  } // for (int irun = runMin; irun < runMax+1; irun++)

  int nSpills = utofTimes.size();
  cout<<"Counting over "<<nSpills<<" matched spills"<<endl;

  int nS1S2utof = 0;
  int nS1S2dtof = 0;
  std::vector<int> nS1S2utofVec;
  nS1S2utofVec.resize(nSpills, 0);
  //nS1S2dtofVec.resize(nSpills, 0);
  // Count number of S1 x S2 hits
  // Only do this for the spills in the spillDB
  int lastut=0;
  TFile *fout = new TFile(Form("%s/%s_out.root", saveDir, utofFile), "recreate");
  for (int spill = 0; spill < nSpills; spill++) {
    TGraph *grTimeSinceSpill = new TGraph();
    double tempUtofSpill = utofTimes[spill];
    double tempDtofSpill = dtofTimes[spill];
    for (int t=lastut; t<tree->GetEntries(); t++ ) {
      tree->GetEntry(t);
      if ((tSoSd/1e9 + startTime) > tempUtofSpill) break;

      if (tTrig !=0 && (tSoSd/1e9 + startTime) == tempUtofSpill) {
	nS1S2utof++;
	nS1S2utofVec[spill]++;
	hDeltatUtof->Fill(tTrig - lastS1S2utof);
	hUtofInSpill->Fill(tTrig/1e9 + startTime - tempUtofSpill);
	grTimeSinceSpill->SetPoint(grTimeSinceSpill->GetN(), tTrig/1e9 + startTime - tempUtofSpill, grTimeSinceSpill->GetN());
	lastS1S2utof = tTrig;
	lastut=t;
      }
    }
    cout<<"Utof spill at "<<tempUtofSpill<<", hits "<<nS1S2utofVec[spill]<<", Dtof spill at "<<tempDtofSpill<<", hits "<<nS1S2dtofVec[spill]<<endl;
    grTimeSinceSpill->SetTitle(Form("Event number versus time since spill start, spill %d; Time since spill start / s; Event no.", spill));
    grTimeSinceSpill->Write(Form("grTimeSinceSpill%d",spill));
  } // for (int spill = 0; spill < nSpills; spill++) 
  
  cout<<"Dtof S1 S2 "<<nS1S2dtof<<", Utof S1 S2 "<<nS1S2utof<<", Ratio "<<(double)nS1S2utof/(double)nS1S2dtof<<endl;

  cout<<"There are "<<nDtofSpills<<" dtof spills, "<<nUtofSpills<<" utof spills"<<endl;

  gStyle->SetOptFit(1);
  TCanvas *cdtof = new TCanvas("cdtof");
  TF1 *f1 = new TF1("f1", "exp([0]+[1]*x)", 500, 1000000);
  hDeltatDtof->Fit("f1", "R");
  hDeltatDtof->Draw("hist");
  f1->Draw("same");
  cdtof->SetGridx();
  cdtof->SetGridy();
  cdtof->Print(Form("%s/%s_deltatDtof.png", saveDir, utofFile));
  cdtof->Print(Form("%s/%s_deltatDtof.pdf", saveDir, utofFile));

  TCanvas *cutof = new TCanvas("cutof");
  TF1 *f2 = new TF1("f2", "exp([0]+[1]*x)", 40000, 1000000);
  hDeltatUtof->Fit("f2", "R");
  cout<<"Utof integral below 40us "<<f2->Integral(0, 40000)/3200<<endl;
  cout<<"Utof integral above 40us (TF1) "<<f2->Integral(40000, 1000000)/3200<<endl;
  cout<<"Utof integral above 40us (hist) "<<hDeltatUtof->Integral(12, 300)<<endl;
  hDeltatUtof->Draw("hist");
  f2->Draw("same");
  cutof->SetGridx();
  cutof->SetGridy();
  cutof->Print(Form("%s/%s_deltatUtof.png", saveDir, utofFile));
  cutof->Print(Form("%s/%s_deltatUtof.pdf", saveDir, utofFile));

  TCanvas *cComb = new TCanvas("cComb");
  cComb->SetLogy();
  hDeltatDtof->SetLineColor(kRed);
  hDeltatDtof->Draw("hist");
  hDeltatUtof->Draw("hist same");
  cComb->SetGridx();
  cComb->SetGridy();
  cComb->Print(Form("%s/%s_deltatComb.png", saveDir, utofFile));
  cComb->Print(Form("%s/%s_deltatComb.pdf", saveDir, utofFile));

  TGraph *grRatioDtof = new TGraph();
  TGraph *grUtofDtof  = new TGraph();
  TTree *s1s2Tree = new TTree("s1s2Tree", "S1S2 Hits");
  s1s2Tree->SetDirectory(0);
  double utof;
  double dtof;
  s1s2Tree->Branch("utof", &utof);
  s1s2Tree->Branch("dtof", &dtof);
  for (int n=0; n < nSpills; n++) {
    grUtofDtof->SetPoint(grRatioDtof->GetN(), nS1S2dtofVec[n], nS1S2utofVec[n]);
    grRatioDtof->SetPoint(grRatioDtof->GetN(), nS1S2dtofVec[n], (double)nS1S2utofVec[n]/(double)nS1S2dtofVec[n]);
  }

  TCanvas *cRatioDtof = new TCanvas("cRatioDtof");
  grRatioDtof->SetTitle(Form("Utof/Dtof S1 #cap S2: %s; S1 #cap S2 dtof; S1 #cap S2 utof/dtof", utofFile));
  grRatioDtof->Draw("AP*");
  cRatioDtof->Print(Form("%s/%s_RatioDtof.png", saveDir, utofFile));
  cRatioDtof->Print(Form("%s/%s_RatioDtof.pdf", saveDir, utofFile));

  TCanvas *cUtofDtof = new TCanvas("cUtofDtof");
  grUtofDtof->SetTitle(Form("Utof vs. Dtof S1 #cap S2: %s; S1 #cap S2 dtof; S1 #cap S2 utof", utofFile));
  grUtofDtof->Draw("AP*");
  cUtofDtof->Print(Form("%s/%s_UtofDtof.png", saveDir, utofFile));
  cUtofDtof->Print(Form("%s/%s_UtofDtof.pdf", saveDir, utofFile));

  TCanvas *cUtofInSpill = new TCanvas("cUtofInSpill");
  hUtofInSpill->Draw("hist");
  cUtofInSpill->Print(Form("%s/%s_UtofInSpill.png", saveDir, utofFile));
  cUtofInSpill->Print(Form("%s/%s_UtofInSpill.pdf", saveDir, utofFile));

  hDeltatDtof->Write();
  hDeltatUtof->Write();
  hUtofInSpill->Write();
  grRatioDtof->Write("grRatioDtof");
  grUtofDtof->Write("grUtofDtof");
  fout->Close();
} // deadtimeTimestamp
