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

Double_t allBeamSpillTimes[2][2000];
Int_t countSpills[2]={0,0};

// Cable delays for the ustof and dstof
Double_t ustofDelay = 184.7;
Double_t dstofDelay = 61.6;

bool isInSpill( double timeA, int itdc);

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
  Double_t coincidenceWindow = 20;

  Double_t lastFakeTimeNs[2][10]; // last fake time ns 0,1 are for PMT A and B, and 0-10 are the bar number
  UInt_t lastUnixTime[2][10];   
  
  Int_t pmtSide   = 0;
  Int_t barNumber = 0;
  Double_t deltat = 0;

  double renorm=1;
  Double_t lastRawBeamSpillNs=0;
  Double_t lastDelayedBeamSpillNs=0;
  Double_t lastUsTofNs=0;
  Double_t tempFakeTimeNs=0;
  Bool_t inSpill=false;
  int countUsTof[2]={0,0};

  Double_t usTofNs;
  
  TTree *usTofTree[2];
  usTofTree[0]= new TTree ("usTofTree_0", "us TOF TDC 1");
  usTofTree[1]= new TTree ("usTofTree_1", "us TOF TDC 2");
  usTofTree[0]->Branch("usTofNs", &usTofNs, "usTofNs/D");
  usTofTree[1]->Branch("usTofNs", &usTofNs, "usTofNs/D");
  
  UInt_t firstTime, lastTime;
  for (int itdc=0; itdc<2; itdc++){
    TFile *toftf = new TFile(filename[itdc].c_str(), "read");
    TTree *toftemp = (TTree*)toftf->Get("tofTree");
    toftemp->SetName("tofTree1");
    RawDsTofHeader *temptof = NULL;
    toftemp->SetBranchAddress("tof",       &temptof       );
    toftemp->GetEntry(0);
    firstTime = temptof->unixTime;
    for (int i=0; i<toftemp->GetEntries(); i++){
      toftemp->GetEntry(i);
      if (temptof->channel==15){
	allBeamSpillTimes[itdc][countSpills[itdc]]=temptof->fakeTimeNs;
	countSpills[itdc]++;
      }else if (temptof->channel==13){
	usTofNs=temptof->fakeTimeNs;
	usTofTree[itdc]->Fill();
      }
    }
  
    toftemp->GetEntry(toftemp->GetEntries()-1);
    lastTime = temptof->unixTime;
    renorm = (lastTime-firstTime)*1./100.;
    delete temptof;
    toftf->Close();
  }

  usTofTree[0]->BuildIndex("usTofNs");
  usTofTree[1]->BuildIndex("usTofNs");
  cout << "Filled usTof Tree " << endl;
  
  TH2D* mapHits = new TH2D("mapHits", "", 2, -0.5, 1.5, 10, 0.5, 10.5);
  TH2D* mapHitsTime = new TH2D("mapHitsTime", "", 100, firstTime, lastTime, 20, 1.5, 21.5);
  TH2D* mapTimeDifference = new TH2D("mapTimeDifference", "", 100, -coincidenceWindow, +coincidenceWindow, 10, 0.5, 10.5);
  TH2D* mapHitsBeamSpill = new TH2D("mapHitsBeamSpill", "", 100, 0, 20e9, 20, 1.5, 21.5);
  TH2D* hitsInSpill = new TH2D("hitsInSpill", "", 100, 0, 1e9, 20, 1.5, 21.5);

  TH2D* barEff = new TH2D("barEff", "", 1, 0.5, 1.5, 10, 0.5, 10.5);
  TH1D* usTof  = new TH1D("usTof", "", 100, 0, 1e9);
  TH1D* coincidenceInSpill = new TH1D("coincidenceInSpill", "", 100, 0, 1e9);
  TH2D* coincidenceInSpillBar = new TH2D("coincidenceInSpillBar", "", 100, 0, 1e9, 10, 0.5, 10.5);

  TH1D* timeBetweenHits = new TH1D("timeBetweenHits", "", 1000, 0, 1e4);
  TH1D* tempHist[2][10];
  for (int ip=0; ip<2; ip++){
    for (int ibar=0; ibar<10; ibar++){
      tempHist[ip][ibar]=new TH1D(Form("h_%d_%d", ip, ibar), "", 1000, 0, 1e4);
    }
  }

  int ustofentry=0;
  for (int itdc=0; itdc<2; itdc++){  

    RawDsTofCoincidence *tofCoin = NULL;
    TTree *tofCoinTree = new TTree("tofTree","HPTPC Time Of Flight");
    tofCoinTree->SetDirectory(0);
    tofCoinTree->Branch("tofCoin", &tofCoin);

    
    TFile *tofFile1 = new TFile(filename[itdc].c_str(), "read");
    TTree *tofTree1 = (TTree*)tofFile1->Get("tofTree");
    tofTree1->SetName("tofTree1");
    RawDsTofHeader *tof = NULL;
    tofTree1->SetBranchAddress("tof",       &tof       );

    
    memset(lastFakeTimeNs, 0, sizeof(lastFakeTimeNs));
    memset(lastUnixTime, 0, sizeof(lastUnixTime));
    pmtSide=0;
    barNumber=0;
    deltat=0;
    lastRawBeamSpillNs=0;
    lastDelayedBeamSpillNs=0;
    lastUsTofNs=0;
    cout << "Entries TDC " << itdc+1 << " " << tofTree1->GetEntries() << endl;

    double tdcpmtmap[10], tdcbarmap[10];
    memset(tdcpmtmap, 0, sizeof(tdcpmtmap));
    memset(tdcbarmap, 0, sizeof(tdcbarmap));

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
	lastRawBeamSpillNs=tof->fakeTimeNs;
	continue;
      }

      if (tof->channel==14){
	lastDelayedBeamSpillNs=tof->fakeTimeNs;
	continue;
      }
      
      if (tof->channel==13){
	lastUsTofNs=tof->fakeTimeNs;
	usTof->Fill(lastUsTofNs-lastDelayedBeamSpillNs);
	countUsTof[itdc]++;
	continue;
      }
      
      pmtSide   = tdcpmtmap[tof->channel-1];
      barNumber = tdcbarmap[tof->channel-1];

      if (pmtSide!=0 && pmtSide!=1) continue;
      //      cout << ientry << " " << pmtSide << " " << barNumber << " " << tof->channel-1 << endl;
      //      cout << tof->channel << " " << endl;

      if (pmtSide==0) deltat = tof->fakeTimeNs - lastFakeTimeNs[1][barNumber-1];
      else deltat = lastFakeTimeNs[0][barNumber-1] - tof->fakeTimeNs;

      mapHitsTime->Fill(tof->unixTime, barNumber*2+pmtSide, 1./renorm);
      mapHitsBeamSpill->Fill(tof->fakeTimeNs-lastDelayedBeamSpillNs, barNumber*2+pmtSide, 1./renorm);
      mapHits->Fill(pmtSide, barNumber);
      mapTimeDifference->Fill(deltat, barNumber);
      //      cout << tof->channel << " " << pmtSide << " " << barNumber << endl;
      timeBetweenHits->Fill(tof->fakeTimeNs - lastFakeTimeNs[pmtSide][barNumber-1]);
      tempHist[pmtSide][barNumber-1]->Fill(tof->fakeTimeNs - lastFakeTimeNs[pmtSide][barNumber-1]);
      lastFakeTimeNs[pmtSide][barNumber-1] = tof->fakeTimeNs;
      lastUnixTime[pmtSide][barNumber-1] = tof->unixTime;

      if ( (tof->fakeTimeNs>lastDelayedBeamSpillNs) && (tof->fakeTimeNs< (lastDelayedBeamSpillNs+1e9))  ){
	inSpill=true;
	hitsInSpill->Fill(tof->fakeTimeNs-tofCoin->lastDelayedBeamSignal, barNumber*2+pmtSide);
      }else{
	inSpill=false;
      }
      
      

      tempFakeTimeNs = tof->fakeTimeNs;
      if (TMath::Abs(deltat)<coincidenceWindow){
	if(tofCoin)
	  delete tofCoin;

	tofCoin = new RawDsTofCoincidence();
	tofCoin->run=tof->run;	  
	tofCoin->tdc=(Short_t)itdc+1;
	tofCoin->bar=(Short_t)barNumber;

	if (pmtSide==0){
	  tofCoin->fakeTimeNs[0] = tof->fakeTimeNs;
	  tofCoin->fakeTimeNs[1] = lastFakeTimeNs[1][barNumber-1];
	  tofCoin->unixTime[0]   = tof->unixTime;
	  tofCoin->unixTime[1]   = lastUnixTime[1][barNumber-1];;
	} else{
	  tofCoin->fakeTimeNs[0] = lastFakeTimeNs[0][barNumber-1];
	  tofCoin->fakeTimeNs[1] = tof->fakeTimeNs;
	  tofCoin->unixTime[0]   = lastUnixTime[0][barNumber-1];;
	  tofCoin->unixTime[1]   = tof->unixTime;
	}
	
	//	tofCoin->inSpill = isInSpill(tof->fakeTimeNs, itdc);
	tofCoin->inSpill = inSpill;
	tofCoin->lastRawBeamSignal = lastRawBeamSpillNs;
	tofCoin->lastDelayedBeamSignal = lastDelayedBeamSpillNs;
	// tofCoin->usTofSignal = lastUsTofNs;
	ustofentry=usTofTree[itdc]->GetEntryNumberWithBestIndex(tof->fakeTimeNs);
	if (ustofentry==-1) {
	  //	  cout << " WE HAVE A PROBLEM, US TOF ENTRY NOT FOUND " << endl;
	  tofCoin->usTofSignal = -1;
	} else {
	  usTofTree[itdc]->GetEntry(ustofentry);
	  tofCoin->usTofSignal = usTofNs;
	}
	tofCoinTree->Fill();

	if ( tofCoin->inSpill==true ){
	  barEff->Fill(1, barNumber) ;
	  coincidenceInSpill->Fill(tof->fakeTimeNs-tofCoin->lastDelayedBeamSignal);
	  coincidenceInSpillBar->Fill(tof->fakeTimeNs-tofCoin->lastDelayedBeamSignal, barNumber);
	}
	lastFakeTimeNs[0][barNumber-1] = lastFakeTimeNs[1][barNumber-1] = 0;
      }
    }

    tofFile1->Close();

    
    cout << "Writing output " << Form("%s/DsTOFcoincidenceRun%d_tdc%d.root", dirname.c_str(), run, itdc+1) << endl;
    TFile *fout = new TFile(Form("%s/DsTOFcoincidenceRun%d_tdc%d.root", dirname.c_str(), run, itdc+1), "recreate");
    tofCoinTree->Write("tofCoinTree");
    fout->Close();
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


  TCanvas *c5 = new TCanvas("c5");
  c5->SetRightMargin(0.18);
  barEff->SetTitle(Form("Run %d: bar coincidence efficiency;;Bar;Efficiency", run));
  barEff->SetMarkerSize(2.0);
  double norm;
  cout << "usTof hits : " << countUsTof[0] << " " << countUsTof[1] << endl;
  for (int ibar=0; ibar<10; ibar++){
    if (ibar<5) norm = 1./countUsTof[1];
    else norm = 1./countUsTof[0];
    barEff->SetBinContent(1, ibar+1, barEff->GetBinContent(1, ibar+1)*norm);

  }
  barEff->Draw("colz text");
  c5->Print(Form("%s/Run%d_barEfficiency.png", dirname.c_str(), run));
  c5->Print(Form("%s/Run%d_barEfficiency.pdf", dirname.c_str(), run));

  TCanvas *c6 = new TCanvas("c6");
  c6->SetRightMargin(0.18);
  coincidenceInSpill->SetLineWidth(2);
  coincidenceInSpill->SetLineColor(kBlue);
  usTof->SetLineWidth(2);
  usTof->SetLineColor(kRed);
  coincidenceInSpill->Scale(1/renorm);
  usTof->Scale(1/renorm);

  double ymax=0;
  if (usTof->GetMaximum()>ymax) ymax = usTof->GetMaximum();
  if (coincidenceInSpill->GetMaximum()>ymax) ymax = coincidenceInSpill->GetMaximum();
  coincidenceInSpill->SetMaximum(ymax*1.1);
  TLegend *leg = new TLegend(0.6, 0.6, 0.89, 0.89);
  leg->AddEntry(usTof, "Upstream TOF", "l");
  leg->AddEntry( coincidenceInSpill, "Downstream TOF", "l");
  coincidenceInSpill->SetTitle(Form("Run %d: bar coincidences in spill;Time from last beam spill [ns];Hz", run));
  coincidenceInSpill->Draw("histo");
  cout << "Total number of coincidences in spill " << coincidenceInSpill->Integral() << endl;
  usTof->Draw("histo same");
  leg->Draw();
  c6->Print(Form("%s/Run%d_coincidenceInSpill.png", dirname.c_str(), run));
  c6->Print(Form("%s/Run%d_coincidenceInSpill.pdf", dirname.c_str(), run));

  TCanvas *c7 = new TCanvas("c7");
  c7->SetRightMargin(0.18);
  c7->SetLogz();
  coincidenceInSpillBar->SetTitle(Form("Run %d: bar coincidences in spill;Time from last beam spill [ns];Bar;Hz", run));
  coincidenceInSpillBar->Scale(1/renorm);
  coincidenceInSpillBar->Draw("colz");
  c7->Print(Form("%s/Run%d_coincidenceInSpillBar.png", dirname.c_str(), run));
  c7->Print(Form("%s/Run%d_coincidenceInSpillBar.pdf", dirname.c_str(), run));

  TCanvas *c8 = new TCanvas("c8");
  c8->SetRightMargin(0.18);
  timeBetweenHits->SetMaximum(300);
  timeBetweenHits->Draw();
  c8->Print(Form("%s/Run%d_timeBetweenHits.png", dirname.c_str(), run));
  c8->Print(Form("%s/Run%d_timeBetweenHits.pdf", dirname.c_str(), run));

  
  TCanvas *c9 = new TCanvas("c9");
  c9->SetRightMargin(0.18);
  c9->SetLogz();
  hitsInSpill->SetTitle(Form("Run %d: hits in spill;Time from last beam spill [ns];PMT;Hz", run));
  hitsInSpill->Scale(1/renorm);
  for (int i=0; i<10; i++){
    hitsInSpill->GetYaxis()->SetBinLabel(2*i+1,Form("%dA", i+1));
    hitsInSpill->GetYaxis()->SetBinLabel(2*i+2,Form("%dB", i+1));
  }
  hitsInSpill->Draw("colz");
  c9->Print(Form("%s/Run%d_hitsInSpillBar.png", dirname.c_str(), run));
  c9->Print(Form("%s/Run%d_hitsInSpillBar.pdf", dirname.c_str(), run));
  
  TFile *fout = new TFile(Form("%s/Run%d_histos.root", dirname.c_str(), run), "recreate");
  mapHits->Write();
  mapTimeDifference->Write();
  mapHitsTime->Write();
  mapHitsBeamSpill->Write();
  barEff->Write();
  coincidenceInSpill->Write();
  coincidenceInSpillBar->Write();
  hitsInSpill->Write();
  timeBetweenHits->Write();

  fout->Close();
}



bool isInSpill( double timeA, int itdc){

  bool inspill = false;

  for (int i=0; i<countSpills[itdc]-1; i++){
    if (timeA>allBeamSpillTimes[itdc][i] && timeA<allBeamSpillTimes[itdc][i+1]){

      if ((allBeamSpillTimes[itdc][i+1]-allBeamSpillTimes[itdc][i])<4e9) return true;
      else return false;

    }
    
  }
  

  return false;

}
