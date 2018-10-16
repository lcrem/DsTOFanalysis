// mergeFiles.C

// Two DsTOFtree*.root files and DsTOFcoincidence*.root files as input
void merge(const char* infile1, const char* infile2, const char* incoinfile1, const char* incoinfile2) {
  gSystem->Load("libdstof.so");
  
  //Raw file load
  TFile *in1 = new TFile(infile1, "read");
  TFile *in2 = new TFile(infile2, "read");
  
  TTree *tofTree1 = (TTree*)in1->Get("tofTree");
  TTree *tofTree2 = (TTree*)in2->Get("tofTree");

  RawDsTofHeader *temptof1 = NULL;
  RawDsTofHeader *temptof2 = NULL;
  
  tofTree1->SetBranchAddress("tof", &temptof1);
  tofTree2->SetBranchAddress("tof", &temptof2);
  // Find the offset between the two TDCs
  // Assume that the average is good enough - unsure if we can do better
  double lastBeamTDC1, lastBeamTDC2 = 0.;
  std::vector<double> beamTDC1, beamTDC2;
  double avgDiff = 0.; // Positive value indicates that TDC1 is ahead of TDC2
  // Loop over both TDCs
  for (size_t i = 0; i < tofTree1->GetEntries(); i++) {
    tofTree1->GetEntry(i);
    // 14 is delayed beam
    if (temptof1->channel == 14) {
      lastBeamTDC1 = temptof1->fakeTimeNs;
      beamTDC1.push_back(lastBeamTDC1);
    }
  } // TDC1 loop
  for (size_t i = 0; i < tofTree2->GetEntries(); i++) {
    tofTree2->GetEntry(i);
    // 14 is delayed beam
    if (temptof2->channel == 14) {
      lastBeamTDC2 = temptof2->fakeTimeNs;
      beamTDC2.push_back(lastBeamTDC2);
    }
  } // TDC2 loop
  std::cout <<"Vector sizes 1, 2: "<<beamTDC1.size()<<", "<<beamTDC2.size()<<std::endl;
  if (beamTDC1.size() != beamTDC2.size()) {
    std::cerr<<"Error: There are different numbers of beam signals in the two TDCs!"<<std::endl;
  }
  double sumDiff = 0.;
  for (size_t i = 0; i < beamTDC1.size(); i++) {
    sumDiff += (beamTDC1[i] - beamTDC2[i]);
  }
  avgDiff = sumDiff / beamTDC1.size();
  // Coincidence file load
  TFile *incoin1 = new TFile(incoinfile1, "read");
  TFile *incoin2 = new TFile(incoinfile2, "read");
  
  TTree *tofTreeCoin1 = (TTree*)incoin1->Get("tofCoinTree");
  TTree *tofTreeCoin2 = (TTree*)incoin2->Get("tofCoinTree");

  RawDsTofCoincidence *tempCoin1 = NULL;
  RawDsTofCoincidence *tempCoin2 = NULL;
  
  tofTreeCoin1->SetBranchAddress("tofCoin", &tempCoin1);
  tofTreeCoin2->SetBranchAddress("tofCoin", &tempCoin2);

  tofTreeCoin1->GetEntry(0);
  int run = tempCoin1->run;
  
  // Create new combined tree
  RawDsTofCoincidence *tempmerge = NULL;
  // New merged tree
  TTree *tofMergedCoin = new TTree("tofMergedCoin", "HPTPC Time Of Flight");
  tofMergedCoin->SetDirectory(0);
  tofMergedCoin->Branch("tofMergedCoin", &tempmerge);
  tempmerge = new RawDsTofCoincidence();
  // Fill after having applied offset in the correct order
  // Treat TDC1 as the master and shift to its time
  int index2 = -1;
  for (size_t i = 0; i < tofTreeCoin1->GetEntries(); i++) {
    tofTreeCoin1->GetEntry(i);
    for (size_t j = index2 + 1; j < tofTreeCoin2->GetEntries(); j++) {
      tofTreeCoin2->GetEntry(j);
      if (tempCoin2->fakeTimeNs[0] + avgDiff < tempCoin1->fakeTimeNs[0]) {
	tempmerge->run = tempCoin2->run;
	tempmerge->tdc = tempCoin2->tdc;
	tempmerge->bar = tempCoin2->bar;
	tempmerge->inSpill = tempCoin2->inSpill;
	tempmerge->fakeTimeNs[0] = tempCoin2->fakeTimeNs[0] + avgDiff;
 	tempmerge->fakeTimeNs[1] = tempCoin2->fakeTimeNs[1] + avgDiff;
	tempmerge->unixTime[0] = tempCoin2->unixTime[0];
 	tempmerge->unixTime[1] = tempCoin2->unixTime[1];
	tempmerge->lastRawBeamSignal = tempCoin2->lastRawBeamSignal + avgDiff;
	tempmerge->lastDelayedBeamSignal = tempCoin2->lastDelayedBeamSignal + avgDiff;
	tempmerge->usTofSignal = tempCoin2->usTofSignal + avgDiff;

	tofMergedCoin->Fill();
	index2 = j;
      }
      else {
	break;
      }
    } // TDC coin 2
    // Now put in 
    tempmerge->run = tempCoin1->run;
    tempmerge->tdc = tempCoin1->tdc;
    tempmerge->bar = tempCoin1->bar;
    tempmerge->inSpill = tempCoin1->inSpill;
    tempmerge->fakeTimeNs[0] = tempCoin1->fakeTimeNs[0];
    tempmerge->fakeTimeNs[1] = tempCoin1->fakeTimeNs[1];
    tempmerge->unixTime[0] = tempCoin1->unixTime[0];
    tempmerge->unixTime[1] = tempCoin1->unixTime[1];
    tempmerge->lastRawBeamSignal = tempCoin1->lastRawBeamSignal;
    tempmerge->lastDelayedBeamSignal = tempCoin1->lastDelayedBeamSignal;
    tempmerge->usTofSignal = tempCoin1->usTofSignal;

    tofMergedCoin->Fill();
.    if (i % 10000 == 0) {
      std::cout<<i<<" TDC1 entries done, "<<index2<<" TDC2 entries done"<<std::endl;
    }
  } // TDC coin 1
  // Put in remaining TDC2 hits
  for (size_t j = index2 + 1; j < tofTreeCoin2->GetEntries(); j++) {
    tofTreeCoin2->GetEntry(j);
    tempmerge->run = tempCoin2->run;
    tempmerge->tdc = tempCoin2->tdc;
    tempmerge->bar = tempCoin2->bar;
    tempmerge->inSpill = tempCoin2->inSpill;
    tempmerge->fakeTimeNs[0] = tempCoin2->fakeTimeNs[0] + avgDiff;
    tempmerge->fakeTimeNs[1] = tempCoin2->fakeTimeNs[1] + avgDiff;
    tempmerge->unixTime[0] = tempCoin2->unixTime[0];
    tempmerge->unixTime[1] = tempCoin2->unixTime[1];
    tempmerge->lastRawBeamSignal = tempCoin2->lastRawBeamSignal + avgDiff;
    tempmerge->lastDelayedBeamSignal = tempCoin2->lastDelayedBeamSignal + avgDiff;
    tempmerge->usTofSignal = tempCoin2->usTofSignal + avgDiff;
    
    tofMergedCoin->Fill();
  }
  TFile *fout = new TFile(Form("DsTOFcoincidenceRun%d_merged.root", run), "recreate");
  std::cout<<"Writing to file"<<std::endl;
  tofMergedCoin->Write("tofMergedCoin");
  fout->Close();
} // merge

// Compares times we have for the delayed beam signal between the two TDCs 
void compareTimes(const char* inFileTDC1, const char* inFileTDC2) {
  gSystem->Load("libdstof.so");

  TFile *in1 = new TFile(inFileTDC1, "read");
  TFile *in2 = new TFile(inFileTDC2, "read");
  
  TTree *tofTree1 = (TTree*)in1->Get("tofTree");
  TTree *tofTree2 = (TTree*)in2->Get("tofTree");

  RawDsTofHeader *temptof1 = NULL;
  RawDsTofHeader *temptof2 = NULL;
  
  tofTree1->SetBranchAddress("tof", &temptof1);
  tofTree2->SetBranchAddress("tof", &temptof2);

  tofTree1->GetEntry(0);
  int run = temptof1->run;

  double lastBeamTDC1, lastBeamTDC2 = 0.;
  
  std::vector<double> beamTDC1, beamTDC2;
  
  // Loop over both TDCs
  for (size_t i = 0; i < tofTree1->GetEntries(); i++) {
    tofTree1->GetEntry(i);
    // 14 is delayed beam
    if (temptof1->channel == 14) {
      lastBeamTDC1 = temptof1->fakeTimeNs;
      beamTDC1.push_back(lastBeamTDC1);
      //      continue;
    }
  } // TDC1 loop
  for (size_t i = 0; i < tofTree2->GetEntries(); i++) {
    tofTree2->GetEntry(i);
    // 14 is delayed beam
    if (temptof2->channel == 14) {
      lastBeamTDC2 = temptof2->fakeTimeNs;
      beamTDC2.push_back(lastBeamTDC2);
      //      continue;
    }
  } // TDC2 loop
  
  std::cout <<"Vector sizes 1, 2: "<<beamTDC1.size()<<", "<<beamTDC2.size()<<std::endl;
  // Now compare the two vectors
  TH1D *hDiffFit = new TH1D("hDiffFit", Form("Run %d: Time difference between beam signals for TDCs; fakeTimeNs [ns]; #Delta t", run), 20, 1520e3, 1520.0008);
  TH1D *hDiff = new TH1D("hDiff", Form("Run %d: Time difference between beam signals for TDCs; fakeTimeNs [ns]; #Delta t", run), 200, 1520e3, 1520.0008);
  TGraph *gDiff = new TGraph;
  gDiff->SetTitle( Form("Run %d: Time difference between beam signals for TDCs; fakeTimeNs [ns]; #Delta t", run) );
  for (size_t i = 0; i < beamTDC1.size(); i++) {
    gDiff->SetPoint( gDiff->GetN(), beamTDC1[i], (beamTDC1[i] - beamTDC2[i]) );
    hDiff->Fill(beamTDC1[i] - beamTDC2[i]);
    hDiffFit->Fill(beamTDC1[i] - beamTDC2[i]);
  }
  new TCanvas;
  gDiff->SetMarkerSize(2);
  gDiff->Draw("AP");
  new TCanvas;
  hDiffFit->Fit("gaus");
  hDiffFit->Draw("hist same");
  new TCanvas;
  hDiff->Draw("hist");
}
