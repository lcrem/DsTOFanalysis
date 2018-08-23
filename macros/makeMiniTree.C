// makeMiniTree.C
// Makes tree with just ustof time and time since last beam spill

//#include "RawDsTofHeader.h"

void makeMiniTree(const char* inFile1) {

  const Double_t spillWindow = 1e9;
  
  gSystem->Load("libdstof.so");

  Double_t lastDelayedBeamSpillNs = 0.;

  TFile *tof1 = new TFile(inFile1, "read");

  TTree *tofTree1 = (TTree*)tof1->Get("tofTree");
  
  RawDsTofHeader *temptof1 = NULL;

  tofTree1->SetBranchAddress("tof", &temptof1);

  tofTree1->GetEntry(0);
  Int_t run = temptof1->run;
  int tdc = temptof1->tdc;

  RawTofMini *tempmini = NULL;
  
  // New tree structure
  TTree *tofMiniTree = new TTree("tofMiniTree", "HPTPC Time Of Flight");
  tofMiniTree->SetDirectory(0);
  tofMiniTree->Branch("tofMiniTree", &tempmini);

  tempmini = new RawTofMini();
  for(size_t i=0; i<tofTree1->GetEntries(); i++){
    tofTree1->GetEntry(i);
    if(temptof1->channel==14) {
      lastDelayedBeamSpillNs=temptof1->fakeTimeNs;
      continue;
    }
    // Only write ustof events
    if(temptof1->channel==13) {
      //RawTofMini *tofminievt = new RawTofMini(); 
      tempmini->lastDelayedBeamSpillNs = lastDelayedBeamSpillNs;
      tempmini->fakeTimeNs = temptof1->fakeTimeNs;
      tempmini->run=run;
      if (tempmini->fakeTimeNs - lastDelayedBeamSpillNs < spillWindow) {
	tempmini->inSpill = true;	
      }
      else {
	tempmini->inSpill = false;
      }
      tofMiniTree->Fill();
      //std::cout<<tempmini->fakeTimeNs - lastDelayedBeamSpillNs<<"\t"<<tempmini->inSpill<<std::endl;
    }
  }
  tof1->Close();
  delete temptof1;

  TFile *fout = new TFile(Form("DsTOFminitreeRun%d_tdc%d.root", run, tdc), "recreate");
  std::cout<<"Writing to file"<<std::endl;
  tofMiniTree->Write("tofMiniTree");
  fout->Close();
}
