// makeSpillDB.cxx

#include "DsTOFConventions.h"
#include "RawDsTofHeader.h"
#include "RawDsTofCoincidence.h"

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2.h"
#include "TColor.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TString.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

using namespace std;

const int maxRunDs = 1442; // Last dstof run
// File directories
const char* usDir = "/nfs/dbhome/jwalding/Dropbox (Royal Holloway)/HPTPC/utofdatabackup_firsthitpinnedtounixtime/";
const char* dsDir = "/nfs/dbhome/jwalding/Dropbox (Royal Holloway)/HPTPC/dtofdatabackup/www.hep.ucl.ac.uk/~lindac/dstof/hptpctof/";

// Checks if a given unix time is within a supplied ustof file
bool checkUstofFile(unsigned int time, TFile *f);


int main() {

  // The final word in beam spill databases
  TTree *spillTree = new TTree("spillTree", "HPTPC beam spills");
  spillTree->SetDirectory(0);
  // All times are in seconds unless otherwise stated
  // Run specific things are outputted for all spills in a given file
  double globalSpillTime; // Define the dtof to be the global spill time
  int unixRunStartDs;     // Timestamp of start of dstof run
  double nsTimeDs;        // Time in ns since the dstof run start
  bool stopOutDs;         // Simple test to see if there are S1,4 coincidences in the spill
  int runDs;              // Run number of this dstof spill 
  double ustofSpillTime; 
  int unixRunStartUs;     // Timestamp of start of ustof run
  double nsTimeUs;        // Time in ns since the ustof run start

  spillTree->Branch("globalSpillTime", &globalSpillTime);
  spillTree->Branch("unixRunStartDs", &unixRunStartDs);
  spillTree->Branch("nsTimeDs", &nsTimeDs);
  spillTree->Branch("runDs", &runDs);
  spillTree->Branch("ustofSpillTime", &ustofSpillTime);
  spillTree->Branch("unixRunStartUs", &unixRunStartUs);
  spillTree->Branch("nsTimeUs", &nsTimeUs);

  // Loop over all the dstof files
  for (int file = 750; file <= maxRunDs; file++) {
    cout<<"Opening run "<<file<<endl;
    const char* filename = Form("%crun%d/DsTOFcoincidenceRun%d_tdc1.root", dsDir, file, file);
    if (!gSystem->AccessPathName(filename)) {
      TFile *inFile = TFile::Open(Form("%crun%d/DsTOFcoincidenceRun%d_tdc1.root", dsDir, file, file), "read");
      TTree *coinTree = (TTree*)inFile->Get("tofCoinTree");

      delete coinTree;
      delete inFile;
    }

  } // Loop over dstof files
} // main

bool checkUstofFiles(unsigned int time, TFile *f) {
  bool isHere = false;

  return isHere;
}
