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

  double renorm=1;
  Double_t lastBeamSpillNs=0;
  
  TFile *toftf = new TFile(filename[0].c_str(), "read");
  TTree *toftemp = (TTree*)toftf->Get("tofTree");
  toftemp->SetName("tofTree1");
  RawDsTofHeader *temptof = NULL;
  toftemp->SetBranchAddress("tof",       &temptof       );
  toftemp->GetEntry(0);
  UInt_t firstTime = temptof->unixTime;
  toftemp->GetEntry(toftemp->GetEntries()-1);
  UInt_t lastTime = temptof->unixTime;
  renorm = (lastTime-firstTime)*1./100.;
  delete temptof;
  toftf->Close();

  TH2D* mapHits = new TH2D("mapHits", "", 2, -0.5, 1.5, 10, 0.5, 10.5);
  TH2D* mapHitsTime = new TH2D("mapHitsTime", "", 100, firstTime, lastTime, 20, 1.5, 21.5);
  TH2D* mapTimeDifference = new TH2D("mapTimeDifference", "", 100, -20, +20, 10, 0.5, 10.5);
  TH2D* mapHitsBeamSpill = new TH2D("mapHitsBeamSpill", "", 100, 0, 20e9, 20, 1.5, 21.5);
  
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

    lastBeamSpillNs=0;
    for (int ientry=0; ientry< tofTree1->GetEntries(); ientry++){

      tofTree1->GetEntry(ientry);

      if (tof->channel==15){
	lastBeamSpillNs=tof->fakeTimeNs;
	continue;
      }
      if (tof->channel==14){
	continue;
      }
      
      pmtSide   = tdcpmtmap[tof->channel-1];
      barNumber = tdcbarmap[tof->channel-1];

      if (pmtSide==0) deltat = tof->fakeTimeNs - lastFakeTimeNs[1][barNumber-1];
      else deltat = lastFakeTimeNs[0][barNumber-1] - tof->fakeTimeNs;

      mapHitsTime->Fill(tof->unixTime, barNumber*2+pmtSide, 1./renorm);
      mapHitsBeamSpill->Fill(tof->fakeTimeNs-lastBeamSpillNs, barNumber*2+pmtSide, 1./renorm);
      mapHits->Fill(pmtSide, barNumber);
      mapTimeDifference->Fill(deltat, barNumber);
      //      cout << tof->channel << " " << pmtSide << " " << barNumber << endl;
      lastFakeTimeNs[pmtSide][barNumber-1] = tof->fakeTimeNs;    
    }

    tofFile1->Close();
  }

  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("c");
  c->SetRightMargin(0.18);
  c->SetLogz();
  mapHits->GetXaxis()->SetBinLabel(1,"PMT A");
  mapHits->GetXaxis()->SetBinLabel(2,"PMT B");
  mapHits->SetTitle(Form("Run %d: hit map;;Bar number;Hz", run));
  // Change z axis into Hz
  mapHits->Scale(1./(lastTime-firstTime));
  mapHits->SetMarkerSize(2.0);
  mapHits->Draw("colz text");
  c->Print(Form("%s/Run%d_hitMap.png", dirname.c_str(), run));
  c->Print(Form("%s/Run%d_hitMap.pdf", dirname.c_str(), run));

  TCanvas *c2 = new TCanvas("c2");
  c2->SetRightMargin(0.18);
  // Renormalise time difference map
  double sum, cont;
  for (int iy=1; iy<=mapTimeDifference->GetYaxis()->GetNbins(); iy++){
    sum=0;
    for (int ix=1; ix<=mapTimeDifference->GetXaxis()->GetNbins(); ix++){
      sum+= mapTimeDifference->GetBinContent(ix, iy);
    }
    for (int ix=1; ix<=mapTimeDifference->GetXaxis()->GetNbins(); ix++){
      cont = mapTimeDifference->GetBinContent(ix, iy);
      mapTimeDifference->SetBinContent(ix, iy, cont/sum);
    }

  }
  
  mapTimeDifference->SetTitle(Form("Run %d: Coincidence Map;(time_{PMT A} - time_{PMT B})[ns];Bar number;Fraction", run));
  mapTimeDifference->Draw("colz");
  c2->Print(Form("%s/Run%d_coincidenceMap.png", dirname.c_str(), run));
  c2->Print(Form("%s/Run%d_coincidenceMap.pdf", dirname.c_str(), run));

  
  TCanvas *c3 = new TCanvas("c3");
  c3->SetRightMargin(0.18);
  c3->SetLogz();
  mapHitsTime->GetXaxis()->SetNdivisions(5);
  mapHitsTime->GetXaxis()->SetTimeDisplay(1);
  mapHitsTime->GetXaxis()->SetTimeFormat("%b-%d %H:%M");
  mapHitsTime->GetXaxis()->SetTimeOffset(0,"gmt");
  for (int i=0; i<10; i++){
    mapHitsTime->GetYaxis()->SetBinLabel(2*i+1,Form("%dA", i+1));
    mapHitsTime->GetYaxis()->SetBinLabel(2*i+2,Form("%dB", i+1));
  }
  mapHitsTime->SetTitle(Form("Run %d: hit time map;Unix time UTC;PMT;Hz", run));
  mapHitsTime->Draw("colz");
  c3->Print(Form("%s/Run%d_hitTimeMap.png", dirname.c_str(), run));
  c3->Print(Form("%s/Run%d_hitTimeMap.pdf", dirname.c_str(), run));

  TCanvas *c4 = new TCanvas("c4");
  c4->SetRightMargin(0.18);
  c4->SetLogz();
  for (int i=0; i<10; i++){
    mapHitsBeamSpill->GetYaxis()->SetBinLabel(2*i+1,Form("%dA", i+1));
    mapHitsBeamSpill->GetYaxis()->SetBinLabel(2*i+2,Form("%dB", i+1));
  }
  mapHitsBeamSpill->SetTitle(Form("Run %d: time from last beam spill;Time from last beam spill [ns];PMT;Hz", run));
  mapHitsBeamSpill->Draw("colz");
  c4->Print(Form("%s/Run%d_hitBeamSpillMap.png", dirname.c_str(), run));
  c4->Print(Form("%s/Run%d_hitBeamSpillMap.pdf", dirname.c_str(), run));


  

  TFile *fout = new TFile(Form("%s/Run%d_histos.root", dirname.c_str(), run), "recreate");
  mapHits->Write();
  mapTimeDifference->Write();
  mapHitsTime->Write();
  mapHitsBeamSpill->Write();
  fout->Close();
}


