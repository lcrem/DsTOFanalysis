// hitMatch.C

// clockComp
// Makes plot of ustof and dstof beam signal times, hopefully showing that they are matched up correctly
void clockComp (const char* ustofFile, const char* dstofFile) {

  gSystem->Load("libdstof.so");
  
  TFile *dstofIn = new TFile(dstofFile, "read");
  TTree *dstofTree = (TTree*)dstofIn->Get("tofCoinTree");

  RawDsTofCoincidence *tempcoin = NULL;
  dstofTree->SetBranchAddress("tofCoin", &tempcoin);
  dstofTree->GetEntry(0);
  const int dstofStart = tempcoin->unixTime[0];
  size_t lastEntry = dstofTree->GetEntries() - 1;
  dstofTree->GetEntry(lastEntry);
  const int dstofEnd = tempcoin->unixTime[0];

  std::cout<<"Dstof start, end "<<dstofStart<<", "<<dstofEnd<<std::endl;
    
  // Load ustof files
  TFile *ustofIn = new TFile(ustofFile, "read");
  TTree *ustofTree = (TTree*)ustofIn->Get("tree");

  double tSoSd;
  double tTrig;
  double tS1;

  TNamed *start = 0;
  TNamed *end   = 0;
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
  
  std::cout<<"Ustof start, end "<<ustofStart<<", "<<ustofEnd<<std::endl;
  
  ustofTree->SetBranchAddress("tSoSd", &tSoSd);
  ustofTree->SetBranchAddress("tTrig", &tTrig);
  ustofTree->SetBranchAddress("tS1", &tS1);

  // Dstof beam times
  std::vector<double> dstofBeamT;
  // ustof beam times
  std::vector<double> ustofBeamT;

  double lastDstofSpill = 0.;
  // Loop over dstof events
  for (size_t t=0; t<dstofTree->GetEntries(); t++) {
    dstofTree->GetEntry(t);
    double delayedBeamUnix = (tempcoin->lastDelayedBeamSignal / 1e9) + dstofStart;
    if (tempcoin->lastDelayedBeamSignal!=lastDstofSpill && delayedBeamUnix>ustofStart && delayedBeamUnix<ustofEnd && tempcoin->lastDelayedBeamSignal-lastDstofSpill>1e9) {
      dstofBeamT.push_back(tempcoin->lastDelayedBeamSignal - (ustofStart - dstofStart)*1e9 );
      lastDstofSpill = tempcoin->lastDelayedBeamSignal;
    }
    //else if (tempcoin->usTofSignal!=lastUstofInDstof && 
  } // dstof events

  double lastUstofSpill = 0.;
  // Ustof events
  for (size_t t=0; t<ustofTree->GetEntries(); t++) {
    ustofTree->GetEntry(t);
    double delayedBeamUnix = (tSoSd / 1e9) + ustofStart;
    if (tSoSd!=lastUstofSpill && tSoSd-lastUstofSpill>1e9 && delayedBeamUnix >= dstofStart && delayedBeamUnix <= dstofEnd) {
      ustofBeamT.push_back(tSoSd);
      lastUstofSpill = tSoSd;
      // hustofT->Fill(tSoSd);
      continue;
    }  
  } // ustof events

  const double dstofFirst = dstofBeamT[0];
  const double ustofFirst = ustofBeamT[0];
  
  // Rewrite with shifted versions
  for (int i = 0; i < dstofBeamT.size(); i++) {
    dstofBeamT[i] -= dstofFirst;
    //    hdstofT->Fill(dstofBeamT[i]);    
  }
  
  for (int i = 0; i < ustofBeamT.size(); i++) {
    ustofBeamT[i] -= ustofFirst;
    //    hustofT->Fill(ustofBeamT[i]);
  }

  TGraph *g1 = new TGraph;
  for (int i=0; i<ustofBeamT.size(); i++) {
    g1->SetPoint(g1->GetN(), ustofBeamT[i], ustofBeamT[i] - dstofBeamT[i]);
  }
  
  TCanvas *c2 = new TCanvas("c2");
  g1->SetTitle(Form("Difference between beam signal times in ustof and dstof (%d to %d); Time since the start of beam matching / ns; #Delta t / ns", dstofStart, dstofEnd));
  g1->Fit("pol1");
  g1->Draw("AP*");
  c2->Print("clock_drift.png");
  c2->Print("clock_drift.pdf");
}

// class for dstof hits
class dstofTemp {

public:
  double pmtTime[2];
  double hitTime;
  int bar;
  double xpos;

  //~dstofTemp();
  dstofTemp() {
    pmtTime[0] = 0.;
    pmtTime[1] = 0.;
    hitTime = 0.;
    bar = 0;
    xpos = 0.; // in cm
  }
  
  void setHitTime() {
    hitTime = min(pmtTime[0], pmtTime[1]) - (10. - abs(pmtTime[0] - pmtTime[1]));
    xpos = (pmtTime[0] - pmtTime[1]) * 7.;
  }  
}; // dstofTemp

// class for ustof hits
class ustofTemp {
  
public:
  int nhit;
  double tS1;
  double tTrig;
  double tToF[50];
  float xToF[50];
  float yToF[50];

  //~ustofTemp();
  ustofTemp() {
    nhit = 0;
    tS1 = 0.;
    tTrig = 0.;
    for (int i=0; i < 50; i++) {
      tToF[i] = 0.;
      xToF[i] = 0.;
      yToF[i] = 0.;
    }
  }
};

// Measure drift
double driftCalc (const char* ustofFile, const char* dstofFile) {
  gSystem->Load("libdstof.so");
  
  std::cout<<"Measuring clock drift for entire run..."<<std::endl;
    
  TFile *dstofIn = new TFile(dstofFile, "read");
  TTree *dstofTree = (TTree*)dstofIn->Get("tofCoinTree");

  RawDsTofCoincidence *tempcoin = NULL;
  dstofTree->SetBranchAddress("tofCoin", &tempcoin);
  dstofTree->GetEntry(0);
  const int dstofStart = tempcoin->unixTime[0];
  size_t lastEntry = dstofTree->GetEntries() - 1;
  dstofTree->GetEntry(lastEntry);
  const int dstofEnd = tempcoin->unixTime[0];

  // Load ustof files
  TFile *ustofIn = new TFile(ustofFile, "read");
  TTree *ustofTree = (TTree*)ustofIn->Get("tree");

  double tSoSd;

  ustofTree->SetBranchAddress("tSoSd", &tSoSd);
  
  TNamed *start = 0;
  TNamed *end   = 0;
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

  // Dstof beam times
  std::vector<double> dstofBeamT;
  // ustof beam times
  std::vector<double> ustofBeamT;

  double lastDstofSpill = 0.;
  // Loop over dstof events
  for (size_t t=0; t<dstofTree->GetEntries(); t++) {
    dstofTree->GetEntry(t);
    double delayedBeamUnix = (tempcoin->lastDelayedBeamSignal / 1e9) + dstofStart;
    if (tempcoin->lastDelayedBeamSignal!=lastDstofSpill && delayedBeamUnix>ustofStart && delayedBeamUnix<ustofEnd && tempcoin->lastDelayedBeamSignal-lastDstofSpill>1e9) {
      dstofBeamT.push_back(tempcoin->lastDelayedBeamSignal - (ustofStart - dstofStart)*1e9 );
      lastDstofSpill = tempcoin->lastDelayedBeamSignal;
    }
     
  } // dstof events

  double lastUstofSpill = 0.;
  // Ustof events
  for (size_t t=0; t<ustofTree->GetEntries(); t++) {
    ustofTree->GetEntry(t);
    double delayedBeamUnix = (tSoSd / 1e9) + ustofStart;
    if (tSoSd!=lastUstofSpill && tSoSd-lastUstofSpill>1e9 && delayedBeamUnix >= dstofStart && delayedBeamUnix <= dstofEnd) {
      ustofBeamT.push_back(tSoSd);
      lastUstofSpill = tSoSd;
      // hustofT->Fill(tSoSd);
      continue;
    }  
  } // ustof events

  const double dstofFirst = dstofBeamT[0];
  const double ustofFirst = ustofBeamT[0];
  
  // Rewrite with shifted versions
  for (int i = 0; i < dstofBeamT.size(); i++) {
    dstofBeamT[i] -= dstofFirst;  
  }
  
  for (int i = 0; i < ustofBeamT.size(); i++) {
    ustofBeamT[i] -= ustofFirst;
  }

  TGraph *g1 = new TGraph;
  for (int i=0; i<ustofBeamT.size(); i++) {
    g1->SetPoint(g1->GetN(), ustofBeamT[i], ustofBeamT[i] - dstofBeamT[i]);
  }
  TF1 *f1 = new TF1("f1","pol1");
  g1->Fit(f1);
  double intercept = f1->GetParameter(0);
  double grad = f1->GetParameter(1);
  return grad;
}

// beamUstofVec
// Returns vector of ustof beam signals
std::vector<double> beamUstofVec (TTree *ustofTree, const int ustofStart, const int dstofStart, const int dstofEnd) {
  std::cout<<"Building ustof vector of beam times"<<std::endl;
  
  double tSoSd;
  ustofTree->SetBranchAddress("tSoSd", &tSoSd);

  double lastUstofSpill = 0.;
  std::vector<double> ustofVec;
  
  for (size_t t=0; t<ustofTree->GetEntries(); t++) {
    ustofTree->GetEntry(t);
    double delayedBeamUnix = (tSoSd / 1e9) + ustofStart;
    if (tSoSd!=lastUstofSpill && tSoSd-lastUstofSpill>1e9 && delayedBeamUnix >= dstofStart && delayedBeamUnix <= dstofEnd) {
      ustofVec.push_back(tSoSd);
      lastUstofSpill = tSoSd;
    }
  }
  
  return ustofVec;
}

// beamDstofVec
// Returns vector of dstof beam signals
std::vector<double> beamDstofVec (TTree *dstofTree, const int dstofStart, const int ustofStart, const int ustofEnd) {
  std::cout<<"Building dstof vector of beam times"<<std::endl;
  
  RawDsTofCoincidence *tempcoin = NULL;
  dstofTree->SetBranchAddress("tofCoin", &tempcoin);
  
  double lastDstofSpill = 0.;
  std::vector<double> dstofVec;
  
  // Loop over dstof events
  for (size_t t=0; t<dstofTree->GetEntries(); t++) {
    dstofTree->GetEntry(t);
    double delayedBeamUnix = (tempcoin->lastDelayedBeamSignal / 1e9) + dstofStart;
    if (tempcoin->lastDelayedBeamSignal!=lastDstofSpill && delayedBeamUnix>ustofStart && delayedBeamUnix<ustofEnd && tempcoin->lastDelayedBeamSignal-lastDstofSpill>1e9) {
      dstofVec.push_back(tempcoin->lastDelayedBeamSignal);
      lastDstofSpill = tempcoin->lastDelayedBeamSignal;
    } 
  } // dstof events

  return dstofVec;
}

// dstofHitVec 
// Return vector of dstof hits in a given spill
std::vector<dstofTemp> dstofHitVec (TTree *dstofTree, const double spillTime) {
  RawDsTofCoincidence *tempcoin = NULL;
  dstofTree->SetBranchAddress("tofCoin", &tempcoin);
  std::vector<dstofTemp> dstofHitVec;
  dstofTemp dstofHit;
  for (int i=0; i<dstofTree->GetEntries(); i++) {
    dstofTree->GetEntry(i);
    // Is within spill
    if(tempcoin->fakeTimeNs[0] - spillTime < 1e9 && tempcoin->fakeTimeNs[0] > spillTime) {
      // with shift from beam signal
      // spill time = 0
      dstofHit.pmtTime[0] = tempcoin->fakeTimeNs[0] - spillTime;
      dstofHit.pmtTime[1] = tempcoin->fakeTimeNs[1] - spillTime;
      dstofHit.bar = tempcoin->bar;
      dstofHit.setHitTime();
      dstofHitVec.push_back(dstofHit);
    }
  }
  return dstofHitVec;
}

// ustofHitVec
// Return vector of ustof hits in a given spill
std::vector<ustofTemp> ustofHitVec (TTree *ustofTree, const double spillTime, const double drift) {
  double tToF[50];
  float xToF[50];
  float yToF[50];
  double tTrig;
  double tS1;
  int nhit;
  ustofTree->SetBranchAddress("tToF", tToF);
  ustofTree->SetBranchAddress("xToF", xToF);
  ustofTree->SetBranchAddress("yToF", yToF);
  ustofTree->SetBranchAddress("tTrig", &tTrig);
  ustofTree->SetBranchAddress("tS1", &tS1);
  ustofTree->SetBranchAddress("nhit", &nhit);
  std::vector<ustofTemp> ustofHitVec;
  ustofTemp ustofHit;
  for (int i=0; i<ustofTree->GetEntries(); i++) {
    ustofTree->GetEntry(i);
    // Is within spill
    if(tToF[0] > spillTime && tToF[0] - spillTime < 1e9) {
      // with shift from beam signal
      // spill time = 0
      ustofHit.nhit = nhit;
      for (int h = 0; h < nhit; h++) {
	ustofHit.tToF[h] = (tToF[h] - spillTime) / (1. + drift);
	ustofHit.xToF[h] = xToF[h];
	ustofHit.yToF[h] = yToF[h];
      }
      ustofHit.tTrig = (tTrig - spillTime) / (1. + drift);
      ustofHit.tS1 = (tS1 - spillTime) / (1. + drift);
      ustofHitVec.push_back(ustofHit);
    }
  }
  return ustofHitVec;
}

// hitMatch
void hitMatch (const char* ustofFile, const char* dstofFile1, const char* dstofFile2, const int spillNo) {

  gSystem->Load("libdstof.so");

  // Load dstof, get start end unix times
  TFile *dstofIn1 = new TFile(dstofFile1, "read");
  TTree *dstofTree1 = (TTree*)dstofIn1->Get("tofCoinTree");
  TFile *dstofIn2 = new TFile(dstofFile2, "read");
  TTree *dstofTree2 = (TTree*)dstofIn2->Get("tofCoinTree");
  
  RawDsTofCoincidence *tempcoin1 = NULL;
  dstofTree1->SetBranchAddress("tofCoin", &tempcoin1);
  dstofTree1->GetEntry(0);
  int run = tempcoin1->run;
  const int dstofStart = tempcoin1->unixTime[0];
  size_t lastEntry = dstofTree1->GetEntries() - 1;
  dstofTree1->GetEntry(lastEntry);
  const int dstofEnd = tempcoin1->unixTime[0];

  RawDsTofCoincidence *tempcoin2 = NULL;
  dstofTree2->SetBranchAddress("tofCoin", &tempcoin2);
  
  // Load ustof files
  TFile *ustofIn = new TFile(ustofFile, "read");
  TTree *ustofTree = (TTree*)ustofIn->Get("tree");

  TNamed *start = 0;
  TNamed *end   = 0;
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
  
  std::cout<<"Ustof start, end "<<ustofStart<<", "<<ustofEnd<<std::endl;
  std::cout<<"Dstof start, end "<<dstofStart<<", "<<dstofEnd<<std::endl;

  // Find the clock drift
  const double drift = driftCalc(ustofFile, dstofFile1);
  std::cout<<"Drift = "<<drift<<std::endl;

  // Now first spill should be aligned for ustof and dstof up to ns level
  // So try and match some hits
  // Get vectors of beam times
  std::vector<double> dstofBeamT = beamDstofVec(dstofTree1, dstofStart, ustofStart, ustofEnd);
  std::vector<double> ustofBeamT = beamUstofVec(ustofTree, ustofStart, dstofStart, dstofEnd);
  std::cout<<"Dstof, Ustof beam signals "<<dstofBeamT.size()<<", "<<ustofBeamT.size()<<std::endl;
  if (dstofBeamT.size() != ustofBeamT.size()) {
    std::cerr<<"Error: Different number of beam signals in this window!"<<std::endl;
  }
  else { }
  
  // Select the spill number that we want
  double myDstofSpill = dstofBeamT[spillNo - 1];
  double myUstofSpill = ustofBeamT[spillNo - 1];
  // Get all the dstof hits up to 1 second after the chosen spill (and in the correct time frame)
  std::cout<<"Getting all the dstof hits in the chosen spill..."<<std::endl;
  // TDC 1
  std::cout<<"For dstof TDC1"<<std::endl;
  std::vector<dstofTemp> dstofHit1 = dstofHitVec(dstofTree1, myDstofSpill);
  // TDC 2
  std::cout<<"For dstof TDC2"<<std::endl;
  std::vector<dstofTemp> dstofHit2 = dstofHitVec(dstofTree2, myDstofSpill);
  std::cout<<"Found "<<(dstofHit1.size() + dstofHit2.size())<<" dstof hits in spill "<<spillNo<<std::endl;

  std::cout<<"Getting all the ustof hits in the chosen spill"<<std::endl;
  std::vector<ustofTemp> ustofHit = ustofHitVec(ustofTree, myUstofSpill, drift);
  std::cout<<"Found "<<ustofHit.size()<<" ustof hits in spill "<<spillNo<<std::endl;
  
  TCanvas *c1 = new TCanvas("c1");
  TH1D *hdstofHits = new TH1D ("dstofHits", Form("Run %d: Aligned hits in ustof and dstof, spill %d; Time since spills start / ns; Hits", run, spillNo), 400, .1e9, .65e9);
  TH1D *hustofHits = new TH1D ("ustofHits", Form("Run %d: Aligned hits in ustof and dstof, spill %d; Time since spills start / ns; Hits", run, spillNo), 400, .1e9, .65e9);
  c1->SetLogy();
  for (int i=0; i<dstofHit1.size(); i++) {   
    hdstofHits->Fill(dstofHit1[i].hitTime);
  }
  for (int i=0; i<dstofHit2.size(); i++) {
    hdstofHits->Fill(dstofHit2[i].hitTime);
  }
  for (int i=0; i<ustofHit.size(); i++) {
    for (int j=0; j<ustofHit[j].nhit; j++) {
      hustofHits->Fill(ustofHit[i].tToF[j]);
    }
  }
  hdstofHits->SetLineColor(kBlue);
  hdstofHits->SetFillColor(kBlue);
  hdstofHits->SetFillStyle(3005);
  hustofHits->SetLineColor(kMagenta);
  hustofHits->SetFillColor(kMagenta);
  hustofHits->SetFillStyle(3004);
  hdstofHits->Draw("hist");
  hustofHits->Draw("hist same");
  TLegend *leg1 = new TLegend(0.7, 0.75, .85, .9);
  leg1->AddEntry(hdstofHits, "DsToF", "f");
  leg1->AddEntry(hustofHits, "UsToF", "f");
  leg1->Draw();
  c1->Print(Form("Run%dspill%d_alignTest.png", run, spillNo));
  c1->Print(Form("Run%dspill%d_alignTest.pdf", run, spillNo));
  // Let's see if we can find out how far apart these hits are
  // Way fewer ustof hits so let's try and match each of these with a dstof hit
  TGraph *g2 = new TGraph;
  TH2D *hdiff2d = new TH2D("hdiff2d", Form("Run %d, spill %d: UsToF - DsToF hit matching over spill; Time since spill / ns; #Deltat / ns", run, spillNo), 100, 0.1e9, 0.65e9, 100, -100000, 100000); 
  TH1D *hdiffOut = new TH1D("hdiffOut", Form("Run %d, spill %d: Time difference between UsToF hit and closest DsToF hit; #Deltat / ns; Events", run, spillNo), 150, -100000, 100000);
  TH1D *hdiff = new TH1D("hdiff", Form("Run %d, spill %d: Time difference between UsToF hit and closest DsToF hit; #Deltat / ns; Events", run, spillNo), 75, -20000, 20000);
  TH1D *hdiffZoom = new TH1D("hdiffZoom", Form("Run %d, spill %d: Time difference between UsToF hit and closest DsToF hit; #Deltat / ns; Events", run, spillNo), 200, -200, 200);
  TH1D *hS1diff = new TH1D("hS1diff", Form("Run %d, spill %d: Time difference between S1 hit and closest DsToF hit; #Deltat / ns; Events", run, spillNo), 200, -20, 200);

  TH2D *goodHitsUstof = new TH2D("goodHitsUstof", Form("Run %d, spill %d: Spatial distribution in S3 for hits in S1, S3, S4; x / cm; y / cm", run, spillNo), 75, -20, 180, 30, 0, 120);
  TH2D *goodHitsDstof = new TH2D("goodHitsDstof", Form("Run %d, spill %d: Distribution in S4 for hits in S1, S3, S4; x / cm; y / cm", run, spillNo), 60, 0, 140, 10, 0, 77.5);

  TH2D *goodHitsDiff = new TH2D("goodHitsDiff", Form("Run %d, spill %d: Spatial difference between S3 and S4 hits; #Delta x / cm; #Delta y / cm", run, spillNo), 60, -150, 150, 30, -100, 100);
  // Plots plots plots
  for (int us=0; us<ustofHit.size(); us++) {
    if (us % 100 == 1) {
      std::cout<<"Ustof event no. "<<us<<std::endl;
    }
    double diffLowS1 = 999999999.;
    double diffS1 = 9999999999.;
    for (int hit=0; hit<ustofHit[us].nhit; hit++) {
      double diff = 999999.;
      double diffLow = 999999999.;
      int hitLow = -1;
      int tdc = 0;
      for (int ds=0; ds<dstofHit1.size(); ds++) {
	diff = ustofHit[us].tToF[hit] - dstofHit1[ds].hitTime;
	//diffS1 = ustofHit[us].tS1 - dstofHit1[ds].hitTime;
	if (abs(diff) < abs(diffLow)) {
	  diffLow = diff;
	  diffLowS1 = ustofHit[us].tS1 - dstofHit1[ds].hitTime;
	  hitLow = ds;
	  tdc = 1;
	}
      }
      for (int ds=0; ds<dstofHit2.size(); ds++) {
	diff = ustofHit[us].tToF[hit] - dstofHit2[ds].hitTime;
	//diffS1 = ustofHit[us].tS1 - dstofHit1[ds].hitTime;
	if (abs(diff) < abs(diffLow)) {
	  diffLow = diff;
	  diffLowS1 = ustofHit[us].tS1 - dstofHit2[ds].hitTime;
	  hitLow = ds;
	  tdc = 2;
	}
      }
      g2->SetPoint(g2->GetN(), ustofHit[us].tToF[hit], diffLow);
      hdiffZoom->Fill(diffLow);
      hdiff->Fill(diffLow);
      hdiffOut->Fill(diffLow);
      hdiff2d->Fill(ustofHit[us].tToF[hit], diffLow);
      // Pick out the best hits and find the x and y positions
      if (diffLow > -100 && diffLow < 200) {
	goodHitsUstof->Fill(ustofHit[us].xToF[hit], ustofHit[us].yToF[hit]);
	if(tdc==1) {
	  goodHitsDstof->Fill((dstofHit1[hitLow].pmtTime[0] - dstofHit1[hitLow].pmtTime[1])*(7./2.)+70., (dstofHit1[hitLow].bar*7.5) - 2.5);
	  goodHitsDiff->Fill(ustofHit[us].xToF[hit] - ((dstofHit1[hitLow].pmtTime[0] - dstofHit1[hitLow].pmtTime[1])*(7. / 2.)+ 70.), ustofHit[us].yToF[hit]-((dstofHit1[hitLow].bar*7.5)-2.5));
	}
	else if(tdc==2) {
	  goodHitsDstof->Fill((dstofHit2[hitLow].pmtTime[0] - dstofHit2[hitLow].pmtTime[1]) * (7. / 2.) + 70, (dstofHit2[hitLow].bar * 7.5) - 2.5);
	  goodHitsDiff->Fill(ustofHit[us].xToF[hit] - ((dstofHit1[hitLow].pmtTime[0] - dstofHit1[hitLow].pmtTime[1])*(7. / 2.)+ 70.), ustofHit[us].yToF[hit]-((dstofHit1[hitLow].bar*7.5)-2.5));
	}
	else {
	  std::cout<<"This should not be happening!"<<std::endl;
	}
      }
    } // nhit
    hS1diff->Fill(diffLowS1);
  } // ustof entries

  TCanvas *c2_1 = new TCanvas("c2_1");
  hdiff->SetLineWidth(2);
  hdiff->Draw("hist same");
  c2_1->Print(Form("Run%dspill%d_tDiff_zoom.png", run, spillNo));
  c2_1->Print(Form("Run%dspill%d_tDiff_zoom.pdf", run, spillNo));
    
  TCanvas *c2_2 = new TCanvas("c2_2");
  hdiffOut->SetLineWidth(2);
  hdiffOut->Draw("hist");
  c2_2->Print(Form("Run%dspill%d_tDiff.png", run, spillNo));
  c2_2->Print(Form("Run%dspill%d_tDiff.pdf", run, spillNo));

  TCanvas *c2_3 = new TCanvas("c2_3");
  gStyle->SetPalette(55);
  hdiff2d->Draw("colz");
  c2_3->Print(Form("Run%dspill%d_tDiff2d.png", run, spillNo));
  c2_3->Print(Form("Run%dspill%d_tDiff2d.pdf", run, spillNo));

  g2->SetTitle(Form("Run %d, spill %d: UsToF - DsToF hit matching over spill; Time since spill / ns; #Deltat / ns", run, spillNo));
  TCanvas *c2_4 = new TCanvas("c2_4");
  g2->Draw("AP*");

  TCanvas *c2_5 = new TCanvas("c2_5");
  hS1diff->Fit("gaus");
  hS1diff->Draw("hist same");
  c2_5->Print(Form("Run%dspill%d_matchedS1.png", run, spillNo));
  c2_5->Print(Form("Run%dspill%d_matchedS1.pdf", run, spillNo));

  TCanvas *c2_6 = new TCanvas("c2_6");
  hdiffZoom->Draw("hist");
  c2_6->Print(Form("Run%dspill%d_tDiffzoomzoom.png", run, spillNo));
  c2_6->Print(Form("Run%dspill%d_tDiffzoomzoom.pdf", run, spillNo));

  TCanvas *c3_1 = new TCanvas("c3_1");
  goodHitsUstof->Draw("box");
  c3_1->Print(Form("Run%dspill%d_s3GoodHits.png", run, spillNo));
  c3_1->Print(Form("Run%dspill%d_s3GoodHits.pdf", run, spillNo));
  TCanvas *c3_2 = new TCanvas("c3_2");
  goodHitsDstof->Draw("box");
  c3_2->Print(Form("Run%dspill%d_s4GoodHits.png", run, spillNo));
  c3_2->Print(Form("Run%dspill%d_s4GoodHits.pdf", run, spillNo));
  TCanvas *c3_3 = new TCanvas("c3_3");
  goodHitsDiff->Draw("box");
  c3_3->Print(Form("Run%dspill%d_GoodHitsDiff.png", run, spillNo));
  c3_3->Print(Form("Run%dspill%d_GoodHitsDiff.pdf", run, spillNo));

}
