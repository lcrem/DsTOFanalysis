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
  UInt_t startTime;
  UInt_t endTime;
  
  if ((argc!=4)|| (argc!=5)){
    std::cerr << "Usage: " << argv[0] << " [start Timestamp] [end Timestamp] [plot name] [optional baseDir]" << std::endl;
    return 1;
  } else {
    startTime = atoi(argv[1]);
    endTime = atoi(argv[2]);
    plotname += argv[3];
    if (argc==5) baseDir += argv[4];
    else  baseDir +=  whereIsMyTOFdata;
  }


  Int_t runMin=-1;
  Int_t runMax=-1;
  
  for (int irun=850; irun<1090; irun++){
    
    TFile *fin = new TFile(Form("%s/run%d/DsTOFcoincidenceRun%d_tdc1.root", baseDir.c_str(), irun, irun), "read");
    RawDsTofCoincidence *tofCoin = NULL;
    TTree *tree = (TTree*) fin->Get("tofCoinTree");
    tree->SetBranchAddress("tofCoin", &tofCoin);
    tree->GetEntry(0);
    UInt_t firstTemp = tofCoin->unixTime[0];
    tree->GetEntry(tree->GetEntries()-1);
    UInt_t lastTemp = tofCoin->unixTime[0];

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

  cout << "Min and max runs are " << runMin << " " << runMax << endl;

  return 0;
}


