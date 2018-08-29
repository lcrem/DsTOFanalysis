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
  
  if((argc!=2)){
    std::cerr << "Usage : " << argv[0] << " [path to directory containing voltages.txt and currents.txt]" << std::endl;
    return 1;
  } else {
    baseDir += argv[1];
  }

   gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("c", "", 1200, 450);
  c->SetRightMargin(0.18);

  Double_t min[2] = { 1000, 100};
  Double_t max[2] = { 1300, 350};
  
  TGraph *g[2][20]; // 0: voltages, 1: currents


  string vars[2] = {"voltages", "currents"};

  for (int ivar=0; ivar<2; ivar++){

    ifstream f((baseDir+"/"+vars[ivar]+".txt").c_str());

    Double_t timestamp[2000];
    Double_t value[25][2000];
    int count = 0 ;
    string line;
    
    if (f.is_open()){
      
      while (getline (f,line)){

	stringstream ss(line);
	ss >> timestamp[count];
	
	for (int ipmt=0; ipmt<20; ipmt++){
	  ss >> value[ipmt][count];
	  // cout << value[ipmt][count] << " " ;
	}
	// cout << endl;
	count++;
	continue;
      }

    }

    TLegend *leg = new TLegend (0.85, 0.1, 0.99, 0.9);
    
    int color = 51;
    for (int ipmt = 0; ipmt<10; ipmt++){
      g[ivar][ipmt] = new TGraph (count, timestamp, value[ipmt]);
      g[ivar][ipmt]->SetLineColor(color);
      g[ivar][ipmt]->SetLineWidth(3);
      if (ipmt==0){
	g[ivar][ipmt]->SetTitle(vars[ivar].c_str());
	g[ivar][ipmt]->GetYaxis()->SetRangeUser(min[ivar], max[ivar]);
	g[ivar][ipmt]->GetXaxis()->SetNdivisions(5);
	g[ivar][ipmt]->GetXaxis()->SetTimeDisplay(1);
	g[ivar][ipmt]->GetXaxis()->SetTimeFormat("%b-%d %H:%M");
	g[ivar][ipmt]->Draw("Al");
      }else{
	g[ivar][ipmt]->Draw("l");
      }
      if (ipmt%2==0)	leg->AddEntry(g[ivar][ipmt], Form("PMT %dA", ipmt/2+1), "l");
      else leg->AddEntry(g[ivar][ipmt], Form("PMT %dB", ipmt/2+1), "l");
      color+=5;

    }
    leg->Draw();
    c->Print(Form("%s/%s_0.png", baseDir.c_str(), vars[ivar].c_str()));
    
    delete leg;
    leg = new TLegend (0.85, 0.1, 0.99, 0.9);
    color = 51;
    for (int ipmt = 10; ipmt<20; ipmt++){
      g[ivar][ipmt] = new TGraph (count, timestamp, value[ipmt]);
      g[ivar][ipmt]->SetLineColor(color);
      g[ivar][ipmt]->SetLineWidth(3);
      if (ipmt==10){
	g[ivar][ipmt]->SetTitle(vars[ivar].c_str());
	g[ivar][ipmt]->GetYaxis()->SetRangeUser(min[ivar], max[ivar]);
	g[ivar][ipmt]->GetXaxis()->SetNdivisions(5);
	g[ivar][ipmt]->GetXaxis()->SetTimeDisplay(1);
	g[ivar][ipmt]->GetXaxis()->SetTimeFormat("%b-%d %H:%M");
	g[ivar][ipmt]->Draw("Al");
      }else{
	g[ivar][ipmt]->Draw("l");
      }
      if (ipmt%2==0)	leg->AddEntry(g[ivar][ipmt], Form("PMT %dA", ipmt/2+1), "l");
      else leg->AddEntry(g[ivar][ipmt], Form("PMT %dB", ipmt/2+1), "l");
      color+=5;

    }
    leg->Draw();
    c->Print(Form("%s/%s_1.png", baseDir.c_str(), vars[ivar].c_str()));

  }

  return 0;
}
