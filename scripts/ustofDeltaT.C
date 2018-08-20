// ustofDeltaT.C
// Plots the time between consecutive ustof hits
// Takes a processed TDC file - doesn't matter if it's TDC 1 or 2, ustof signal
// goes into the same channel on both
// Requires libRawDsTofHeader.so be loaded
#include "RawDsTofHeader.h"

void ustofDeltaT(const char* inFile){
  
  //gSystem->Load("libRawDsTofHeader.so");
    
  TFile *tofIn = new TFile(inFile, "read");
  TTree *tofTree = (TTree*)tofIn->Get("tofTree");

  RawDsTofHeader *temptof = NULL;
  tofTree->SetBranchAddress("tof", &temptof);
  
  tofTree->GetEntry(0);
  int run = temptof->run;
  TH1D *hUstofDeltaT = new TH1D("hUstofDeltaT", "", 100, 0, 1e6); 
  
  double ustofOld = 0.;
  double ustofNew = 0.;
  for(size_t i=0; i<tofTree->GetEntries(); i++){
    tofTree->GetEntry(i);
    if(temptof->channel==13){ // ustof channel
      ustofNew = temptof->fakeTimeNs;
      double deltaT = ustofNew - ustofOld;
      hUstofDeltaT->Fill(deltaT);
      ustofOld = ustofNew;
    }
    else { }
  }
  TCanvas *c1 = new TCanvas("c1", "c1");
  hUstofDeltaT->SetTitle(Form("Run %d: Time between ustof signals; Time / ns; Number", run));
  hUstofDeltaT->Draw("hist");
  c1->Print(Form("Run%d_ustofDeltaT.png", run));
  c1->Print(Form("Run%d_ustofDeltaT.pdf", run));
}
