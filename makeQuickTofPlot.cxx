#include "DsTOFConventions.h"
#include "RawDsTofHeader.h"
#include "RawDsTofCoincidence.h"

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
#include "TLegend.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

using namespace std;


int main(int argc, char *argv[]){

  
  string baseDir;
  string plotname;
  string whereToSave;
  UInt_t startTime;
  UInt_t endTime;
  
  if ( (argc!=5) && (argc!=6)){
    std::cerr << "Usage: " << argv[0] << " [start Timestamp] [end Timestamp] [plot name] [whereToSaveStuff] [optional baseDir]" << std::endl;
    return 1;
  } else {
    startTime = atoi(argv[1]);
    endTime = atoi(argv[2]);
    plotname += argv[3];
    whereToSave += argv[4];
    if (argc==6) baseDir += argv[5];
    else  baseDir +=  whereIsMyTOFdata;
  }

  if (startTime>endTime){
    cout << "Start time is greater then end time, check your timestamps! " << endl;
    return -1;
  }


  Int_t runMin=-1;
  Int_t runMax=-1;
  

  for (int irun=1000; irun<1090; irun++){
    
    TFile *fin = new TFile(Form("%s/run%d/DsTOFcoincidenceRun%d_tdc1.root", baseDir.c_str(), irun, irun), "read");
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

  }

  if (runMin==-1 || runMax==-1){
    cout << " Couldn't start or end run, check your timestamps!! " << endl;
    return -1;
  }

  cout << "Min and max runs are " << runMin << " " << runMax << endl;

  TH1D* hTof = new TH1D("hTof", "", 300, 0, 250);

  for (int itdc=0; itdc<2; itdc++){
    TChain *tofCoinChain = new TChain("tofCoinTree");

    for (int irun=runMin; irun<runMax+1; irun++){
      tofCoinChain->Add(Form("%s/run%d/DsTOFcoincidenceRun%d_tdc%d.root", baseDir.c_str(), irun, irun, itdc+1));
    }

    RawDsTofCoincidence *tofCoin = NULL;
    tofCoinChain->SetBranchAddress("tofCoin", &tofCoin);

    double dstofHitT, deltat;

    for (int ientry=0; ientry<tofCoinChain->GetEntries(); ientry++){
      tofCoinChain->GetEntry(ientry);
      if (tofCoin->unixTime[0]<startTime) continue;
      if (tofCoin->unixTime[0]>endTime) break;
      
      
      deltat = TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1]  );
      dstofHitT = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (10. - TMath::Abs(deltat) / 2 );
      hTof->Fill(dstofHitT - tofCoin->usTofSignal);
      
    }
  }
  
  TCanvas *c1 = new TCanvas("c1");
  hTof->SetTitle(Form("%s: Time of flight, DsToF - UsToF; DsToF - UsToF / ns; Events", plotname.c_str()));
  hTof->SetLineColor(kRed);
  hTof->SetFillColor(kRed);
  hTof->SetFillStyle(3005);
  hTof->Draw("hist");
  c1->Print(Form("%s/%s.png", whereToSave.c_str(), plotname.c_str()));
  c1->Print(Form("%s/%s.pdf", whereToSave.c_str(), plotname.c_str()));

  TFile *fout = new TFile (Form("%s/%s.root", whereToSave.c_str(), plotname.c_str()), "recreate");
  hTof->Write(Form("htof_%s", plotname.c_str()));

  fout->Close();


  return 0;
}


