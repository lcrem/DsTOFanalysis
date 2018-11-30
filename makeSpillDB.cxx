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
// Makes of tree of all ustof file start end times
TTree *ustofFileTree();

vector<pair<vector<int>, TString>> ustofFileVec();
// Checks if a given unix time is within ustof files
string checkUstofFiles(double time, TTree *t);

string checkUstofFilesVec(double time, vector<pair<vector<int>, TString>> fileVec);

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

  //  TTree *usFiles = ustofFileTree();
  vector<pair<vector<int>, TString>> usFiles;
  // Loop over all the dstof files
  for (unsigned int file = firstRun; file <= lastRun; file++) {

    const char* filename = Form("%srun%d/DsTOFcoincidenceRun%d_tdc1.root", dsDir, file, file);
    if (!gSystem->AccessPathName(filename)) {
      cout<<"Opening run "<<file<<endl;
      TFile *inFile = new TFile(Form("%srun%d/DsTOFcoincidenceRun%d_tdc1.root", dsDir, file, file), "read");
      TTree *coinTree = (TTree*)inFile->Get("tofCoinTree");
      RawDsTofCoincidence *tempcoin = NULL;
      coinTree->SetBranchAddress("tofCoin", &tempcoin);
      coinTree->GetEntry(0);
      int runStart = tempcoin->unixTime[0];
      cout<<"Start: "<<runStart<<endl;
      // Find first spill in the file
      double lastdelayed = 0.;
      double firstSpillNs   = 0.;
      double firstSpillUnix = 0.;
      for (int t = 0; t < coinTree->GetEntries(); t++) {
	coinTree->GetEntry(t);
	if (tempcoin->lastDelayedBeamSignal > 0. && tempcoin->lastDelayedBeamSignal - tempcoin->lastRawBeamSignal > 0.5e9 && tempcoin->lastDelayedBeamSignal - tempcoin->lastRawBeamSignal < 1e9) {
	  firstSpillNs = tempcoin->lastDelayedBeamSignal;
	  firstSpillUnix = (tempcoin->lastDelayedBeamSignal / 1e9) + runStart;
	  cout.precision(17);
	  cout<<firstSpillUnix<<endl;
	  // Now search through ustof files to find in which one this lies
	  string goodFile = checkUstofFilesVec(firstSpillUnix, usFiles);
	  if (goodFile != "nofile") {
	    cout<<goodFile<<endl;
	  }
	  break;
	}
      } // Loop over entries
            
      delete coinTree;
      delete inFile;
    }

  } // Loop over dstof files
} // main

string checkUstofFilesVec(double time, vector<pair<vector<int>, TString> > fileVec) {
  string ustofFile = "nofile";
  for (int i=0; i<fileVec.size(); i++) {
    if (time >= fileVec[i].first[0] && time <= fileVec[i].first[1]) {
      ustofFile = fileVec[i].second;
    }
  }
  return ustofFile;
}

string checkUstofFiles(double time, TTree *fileTree) {
  string ustofFile = "nofile";
  string fileName;
  int start;
  int end;
  fileTree->SetBranchAddress("fileName", &fileName);
  fileTree->SetBranchAddress("start", &start);
  fileTree->SetBranchAddress("end", &end);
  // Loop over file tree
  for (int t = 0; t < fileTree->GetEntries(); t++) {
    fileTree->GetEntry(t);
    if (time >= start && time <= end) {
      ustofFile = fileName;
    }
  }

  return ustofFile;
}

TTree *ustofSpillTree() {
  const char* ext   = ".root";
  const char* pref  = "Data";
  const char* indir = "/home/sjones/mylinktoutof/";

  TString str;
  const char *entry;
  vector<pair<int, TString> > fileVec;

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
      TNamed *end = 0;
      ustofFile->GetObject("start_of_run", start);
      ustofFile->GetObject("end_of_run", end);
      const char* startchar = start->GetTitle();
      std::string startstr(startchar);
      std::string unixstart = startstr.substr(25,10);
      const int ustofStart = std::stoi(unixstart);	 
      
      fileVec.push_back(make_pair(ustofStart, strTemp));
      
      delete start;
      delete end;
      delete tree;
      delete ustofFile;
    }
  }

  // Sort files by start time
  sort(begin(fileVec), end(fileVec));
  cout<<"Files: "<<fileVec.size()<<endl;

  // New ustof beam spill tree
  TTree *ustofBeamTree = new TTree("ustofBeamTree", "HPTPC ustof beam spill");
  ustofBeamTree->SetDirectory(0);

  string fileName;
  int unixRunStart;
  double nsTime;
  double fakeUnixTime; // unix + some time in s

  ustofBeamTree->Branch("fileName", &fileName);
  ustofBeamTree->Branch("unixRunStart", &unixRunStart);
  ustofBeamTree->Branch("nsTime", &nsTime);
  ustofBeamTree->Branch("fakeUnixTime", &fakeUnixTime);

  for (int i=0; i<fileVec.size(); i++) {
    std::cout<<"Opening file "<<i<<std::endl;
    TString opendir = (fileVec[i].second).Prepend(indir);
    if(!gSystem->AccessPathName(opendir)){//.c_str())){
      TFile *inFile = new TFile(opendir,"read");//.c_str(), "read");
      TTree *inTree = (TTree*)inFile->Get("tree");
      double tSoSd;
      inTree->SetBranchAddress("tSoSd", &tSoSd);

      double lastSpill = 0.;
      for (int t = 0; t < inTree->GetEntries(); t++) {
	inTree->GetEntry(t);
	if (tSoSd > 0 && tSoSd != lastSpill && tSoSd - lastSpill > 1e9) { //&& tSoSd - lastSpill > 0.) {
	  fileName = fileVec[i].second;
	  unixRunStart = fileVec[i].first;
	  lastSpill = tSoSd;
	  nsTime = tSoSd;
	  fakeUnixTime = fileVec[i].first + (tSoSd / 1e9);
	  ustofBeamTree->Fill();
	}
      }
      
      delete inTree;
      delete inFile;
    }
    else {
      cout<<"Couldn't open file!"<<endl;
    }
  }
  return ustofBeamTree;
}

TTree *ustofFileTree() {
  cout<<"Making vector of ustof start/end times"<<endl;
  const char* ext   = ".root";
  const char* pref  = "Data";
  const char* indir = "/home/sjones/mylinktoutof/";

  TString str;
  const char *entry;
  vector<pair<vector<int>, TString> > fileVec;

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
      int ustofStart = stoi(unixstart);	 
      int ustofEnd   = stoi(unixend);	 
      vector<int> ustofTmp =  {ustofStart, ustofEnd};
      fileVec.push_back(make_pair(ustofTmp, strTemp));
     
      delete start;
      delete end;
      delete tree;
      delete ustofFile;
    }
  }
  sort(begin(fileVec), end(fileVec));
  cout<<"Files: "<<fileVec.size()<<endl;


  // Make new tree
  TTree *ustofFileTree = new TTree("ustofFileTree", "HPTPC ustof files");
  ustofFileTree->SetDirectory(0);
  string fileName;
  int start;
  int end;
  ustofFileTree->Branch("fileName", &fileName);
  ustofFileTree->Branch("start", &start);
  ustofFileTree->Branch("end", &end);
  for (int i=0; i < fileVec.size(); i++) {
    fileName = fileVec[i].second;
    start = fileVec[i].first[0];
    end = fileVec[i].first[1];
    ustofFileTree->Fill();
  }

  for (int t=0; t<ustofFileTree->GetEntries(); t++) {
    ustofFileTree->GetEntry(t);
    cout<<fileName<<" "<<start<<" "<<end<<endl;
  }

  return ustofFileTree;
} // ustofFileTree

vector<pair<vector<int>, TString> >  ustofFileVec() {
  cout<<"Making vector of ustof start/end times"<<endl;
  const char* ext   = ".root";
  const char* pref  = "Data";
  const char* indir = "/home/sjones/mylinktoutof/";

  TString str;
  const char *entry;
  vector<pair<vector<int>, TString> > fileVec;

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
      int ustofStart = stoi(unixstart);	 
      int ustofEnd   = stoi(unixend);	 
      vector<int> ustofTmp =  {ustofStart, ustofEnd};
      fileVec.push_back(make_pair(ustofTmp, strTemp));
     
      delete start;
      delete end;
      delete tree;
      delete ustofFile;
    }
  }
  sort(begin(fileVec), end(fileVec));
  cout<<"Files: "<<fileVec.size()<<endl;

  return fileVec;
} // ustofFileVec
