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
#include <string>
#include <fstream>
#include <sstream>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <algorithm>

using namespace std;

const int maxRunDs = 1442; // Last dstof run
// File directories
const char* usDir = "/home/sjones/mylinktoutof/";
const char* dsDir = "/scratch0/dbrailsf/temp/mylinktodtof/";

// Makes tree of all the ustof beam spill times
TTree *ustofSpillTree(int i);
// Builds vector of ustof file start/end times
vector<pair<vector<double>, TString>> ustofFileVec();
// Checks all ustof files to see if they contain a spill at a given unix time
TString checkUstofFiles(double time, vector<pair<vector<double>, TString>> fileVec);
// For a given utof file and time returns utof spill time if there is a spill within a certain time frame
// Returns -1 otherwise
double findSpill(double time, TString fileName, double window = 3.);

int main(int argc, char *argv[]) {

  unsigned int firstRun;
  unsigned int lastRun;
  if(argc!=3) {
    cerr<<"Usage: "<<argv[0]<<"[First dstof run] [Last dstof run]"<<endl;
    return 1;
  }
  else {
    firstRun = atoi(argv[1]);
    lastRun  = atoi(argv[2]);
  }

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

  vector<pair<vector<double>, TString>> usFiles = ustofFileVec();
  // Loop over all the dstof files
  for (unsigned int file = firstRun; file <= lastRun; file++) {
    const char* filename = Form("%srun%d/DsTOFcoincidenceRun%d_tdc1.root", dsDir, file, file);
    if (!gSystem->AccessPathName(filename)) {
      cout<<"Opening run "<<file<<endl;
      TFile *inFile1 = new TFile(Form("%srun%d/DsTOFcoincidenceRun%d_tdc1.root", dsDir, file, file), "read");
      TFile *inFile2 = new TFile(Form("%srun%d/DsTOFcoincidenceRun%d_tdc2.root", dsDir, file, file), "read");
      TTree *coinTree1 = (TTree*)inFile1->Get("tofCoinTree");
      TTree *coinTree2 = (TTree*)inFile2->Get("tofCoinTree");
      RawDsTofCoincidence *tempcoin1 = NULL;
      RawDsTofCoincidence *tempcoin2 = NULL;
      coinTree1->SetBranchAddress("tofCoin", &tempcoin1);
      coinTree2->SetBranchAddress("tofCoin", &tempcoin2);
      coinTree1->GetEntry(0);
      int runStart = tempcoin1->unixTime[0];
      cout<<"Start: "<<runStart<<endl;
      // Find first spill in the file
      double lastdelayed = 0.;
      double firstSpillNs   = 0.;
      double firstSpillUnix = 0.;
      TString goodFile;
      for (int t = 0; t < coinTree1->GetEntries(); t++) {
	coinTree1->GetEntry(t);
	if (tempcoin1->lastDelayedBeamSignal > 0. && tempcoin1->lastDelayedBeamSignal - tempcoin1->lastRawBeamSignal > 0.5e9 && tempcoin1->lastDelayedBeamSignal - tempcoin1->lastRawBeamSignal < 1e9) {
	  firstSpillNs = tempcoin1->lastDelayedBeamSignal;
	  firstSpillUnix = (tempcoin1->lastDelayedBeamSignal / 1e9) + runStart;
	  cout.precision(17);
	  cout<<firstSpillUnix<<endl;
	  // Now search through ustof files to find in which one this lies
	  goodFile = checkUstofFiles(firstSpillUnix, usFiles); 
	  cout<<goodFile<<endl;
	  break;
	}
      } // Loop over entries
      // Try and find matching spill within the ustof file
      if (goodFile != "nofile") {
	double spill = findSpill(firstSpillUnix, goodFile);
      }

      // Find entry with tSoSd within 3s of given spill
      delete tempcoin1;
      delete coinTree1;
      delete inFile1;
      delete tempcoin2;
      delete coinTree2;
      delete inFile2;
    }

  } // Loop over dstof files
} // main

double findSpill(double time, TString fileName, double window) {
  double spillT = -1.;
  cout<<"Looking for matching spill in ustof file "<<fileName<<endl;
  // Open file
  TString path = fileName;
  path.Prepend(usDir);
  TFile *ustofFile = new TFile(path, "read");
  TTree *tree = (TTree*)ustofFile->Get("tree");
  TNamed *start = 0;
  ustofFile->GetObject("start_of_run", start);
  const char* startchar = start->GetTitle();
  string startstr(startchar);
  string unixstart = startstr.substr(25,10);
  double ustofStart = stod(unixstart);	 
  // Branch vars
  double tSoSd;
  double tToF[50];
  int nhit;
  tree->SetBranchAddress("tSoSd", &tSoSd);
  tree->SetBranchAddress("tToF", tToF);
  tree->SetBranchAddress("nhit", &nhit);
  double lastUnix = 0.;
  int nMatch = 0;
  // Loop over entries
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
    // Find tSoSd that match within window
    double ustofUnix = ustofStart + (tSoSd / 1e9);
    if (ustofUnix >= (time - window) && ustofUnix <= (time + window) && ustofUnix != lastUnix) {
      cout<<"Spill at "<<ustofUnix<<endl;
      spillT = ustofUnix;
      nMatch++;
      lastUnix = ustofUnix;
    }
  }
  cout<<nMatch<<" matching spills"<<"\n"<<endl;
  delete tree;
  delete ustofFile;

  return spillT;
}

TString checkUstofFiles(double time, vector<pair<vector<double>, TString> > fileVec) {
  TString ustofFile = "nofile";
  for (int i=0; i<fileVec.size(); i++) {
    if (time >= fileVec[i].first[0] && time <= fileVec[i].first[1]) {
      ustofFile = fileVec[i].second;
    }
  }
  return ustofFile;
}

vector<pair<vector<double>, TString> >  ustofFileVec() {
  cout<<"Making vector of ustof start/end times"<<endl;
  const char* ext   = ".root";
  const char* pref  = "Data";
  const char* indir = "/home/sjones/mylinktoutof/";

  TString str;
  const char *entry;
  vector<pair<vector<double>, TString> > fileVec;

  char *dir  = gSystem->ExpandPathName(indir);
  void *dirp = gSystem->OpenDirectory(dir);
  TString bad("8_24_b7_800MeV_4block_ch46");

  while (entry = (char*)gSystem->GetDirEntry(dirp)) { 
    str = entry;
    if (str.EndsWith(ext) && !str.Contains(bad) && str.Contains(pref)) {
      TString strTemp = str;
      str.Prepend(indir);
      
      TFile *ustofFile = new TFile(str, "read");
      TTree *tree = (TTree*)ustofFile->Get("tree");
      TNamed *start = 0;
      TNamed *end   = 0;
      ustofFile->GetObject("start_of_run", start);
      ustofFile->GetObject("end_of_run", end);
      const char* startchar = start->GetTitle();
      const char* endchar   = end->GetTitle();
      string startstr(startchar);
      string endstr(endchar);
      string unixstart = startstr.substr(25,10);
      string unixend   = endstr.substr(23,10);
      double ustofStart = stod(unixstart);	 
      double ustofEnd   = stod(unixend);	 
      vector<double> ustofTmp =  {ustofStart, ustofEnd};
      fileVec.push_back(make_pair(ustofTmp, strTemp));
     
      delete start;
      delete end;
      delete tree;
      delete ustofFile;
    }
  }
  sort(begin(fileVec), end(fileVec));
  cout<<"Files: "<<fileVec.size()<<endl;
  for (int i=0; i<fileVec.size(); i++) {
    cout<<fileVec[i].first[0]<<" "<<fileVec[i].first[1]<<" "<<fileVec[i].second<<endl;
  }
  return fileVec;
} // ustofFileVec
