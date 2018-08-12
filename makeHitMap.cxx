#include "DsTOFConventions.h"
#include "RawDsTofHeader.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2.h"
#include "TColor.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

using namespace std;

int main(int argc, char *argv[]){

  Int_t run;
  string baseDir;
  
  if((argc!=2)&&(argc!=3)){
    std::cerr << "Usage 1: " << argv[0] << " [irun] (basedir)" << std::endl;
    return 1;
  } else {
    run = atoi(argv[1]);
  if (argc==3) baseDir += argv[2];
  else  baseDir +=  whereIsMyTOFdata;
  }

  string dirname = Form("%s/run%d/", baseDir.c_str(), run);

  string filename[2];
  for (int itdc=0; itdc<2; itdc++){
    filename[itdc] = Form("%s/DsTOFtreeRun%d_tdc%d.root", dirname.c_str(), run, itdc+1);
  }
  
  if (filename[0]=="nothing" || filename[1]=="nothing"){
    cout << "Could not find parsed files in directory " << dirname << endl;
  } else {
    cout << "Path to TDC1 parsed file " << filename[0] << endl;
    cout << "Path to TDC2 parsed file " << filename[1] << endl;
  }

  string line;
  double coincidenceWindow = 20;

  Double_t lastFakeTimeNs[2][10]; // last fake time ns 0,1 are for PMT A and B, and 0-10 are the bar number
  
  int pmtSide   = 0;
  int barNumber = 0;
  Double_t deltat = 0;

  TH2D* mapHits = new TH2D("mapHits", "", 2, -0.5, 1.5, 10, 0.5, 10.5);
  TH2D* mapTimeDifference = new TH2D("mapTimeDifference", "", 100, -20, +20, 10, 0.5, 10.5);


  
  for (int itdc=0; itdc<2; itdc++){  
  
    TFile *tofFile1 = new TFile(filename[itdc].c_str(), "read");
    TTree *tofTree1 = (TTree*)tofFile1->Get("tofTree");
    tofTree1->SetName("tofTree1");
    RawDsTofHeader *tof = NULL;
    tofTree1->SetBranchAddress("tof",       &tof       );

    memset(lastFakeTimeNs, 0, sizeof(lastFakeTimeNs));
    pmtSide=0;
    barNumber=0;
    deltat=0;
    cout << "Entries TDC " << itdc+1 << " " << tofTree1->GetEntries() << endl;

    double tdcpmtmap[10], tdcbarmap[10];
    if (itdc==0){
      for (int i=0; i<10; i++){
	tdcpmtmap[i]=tdc1pmt[i];
	tdcbarmap[i]=tdc1bar[i];
      }
    }else{
      for (int i=0; i<10; i++){
	tdcpmtmap[i]=tdc2pmt[i];
	tdcbarmap[i]=tdc2bar[i];
      }
    }
    
    for (int ientry=0; ientry< tofTree1->GetEntries(); ientry++){

      tofTree1->GetEntry(ientry);

      if (tof->channel==15){
	continue;
      }
      if (tof->channel==14){
	continue;
      }
      
      pmtSide   = tdcpmtmap[tof->channel-1];
      barNumber = tdcbarmap[tof->channel-1];

      if (pmtSide==0) deltat = tof->fakeTimeNs - lastFakeTimeNs[1][barNumber-1];
      else deltat = lastFakeTimeNs[0][barNumber-1] - tof->fakeTimeNs;
    
      mapHits->Fill(pmtSide, barNumber);
      mapTimeDifference->Fill(deltat, barNumber);
      //      cout << tof->channel << " " << pmtSide << " " << barNumber << endl;
      lastFakeTimeNs[pmtSide][barNumber-1] = tof->fakeTimeNs;    
    }

    tofFile1->Close();
  }

  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("c");
  c->SetLogz();
  mapHits->GetXaxis()->SetBinLabel(1,"PMT A");
  mapHits->GetXaxis()->SetBinLabel(2,"PMT B");
  mapHits->SetTitle(Form("Run %d: hit map;;Bar number", run));
  mapHits->Draw("colz");
  c->Print(Form("%s/Run%d_hitMap.png", dirname.c_str(), run));
  c->Print(Form("%s/Run%d_hitMap.pdf", dirname.c_str(), run));

  TCanvas *c2 = new TCanvas("c2");
  mapTimeDifference->SetTitle(Form("Run %d: Coincidence Map;(time_{PMT A} - time_{PMT B})[ns];Bar number", run));
  mapTimeDifference->Draw("colz");
  c2->Print(Form("%s/Run%d_coincidenceMap.png", dirname.c_str(), run));
  c2->Print(Form("%s/Run%d_coincidenceMap.pdf", dirname.c_str(), run));

  TFile *fout = new TFile(Form("%s/Run%d_histos.root", dirname.c_str(), run), "recreate");
  mapHits->Write();
  mapTimeDifference->Write();
  fout->Close();
}


