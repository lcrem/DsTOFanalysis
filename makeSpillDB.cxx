// makeSpillDB.cxx

#include "DsTOFConventions.h"
#include "RawDsTofHeader.h"
#include "RawDsTofCoincidence.h"

#include "TFile.h"
#include "TChain.h"
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
const double stdDrift = -2.8e-6;


// Makes tree of all the ustof beam spill times
TTree *ustofSpillTree(int i);
// Builds vector of ustof file start/end times
vector<pair<vector<double>, TString>> ustofFileVec();
// Checks all ustof files to see if they contain a spill at a given unix time
TString checkUstofFiles(double time, vector<pair<vector<double>, TString>> fileVec);
// For a given utof file and time returns utof spill time if there is a spill within a certain time frame
// Returns -1 otherwise
// First element: Utof Ns time
// Second element: Utof unix time
// Third element: Dtof Ns time
// Fourth element: Dtof unix time
vector<double> findSpill(double timeDs, double timeNsDs, TTree *dsFile1, TTree* dsFile2, TString fileName, double window = 3.);
// Checks if the hits match within a given spill
bool spillMatch(double ustofSpillT, double dstofSpillT);
// Go through spills, having found offset and attempt to match
// First pair dstof elements, second pair ustof elements. First is ns time, second is unix time
vector<pair< pair<double,double>, pair<double,double> >> spillCount(double lastSpillDs, 
								    double lastSpillUs,
								    TTree *dsFile,
								    TChain chain,
								    TString ustofFile,
								    double window = 0.3);

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

  TChain usCh("tree");

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
  for(int t=0; t<usFiles.size(); t++) {
    // Get full file path
    TString path = usFiles[t].second;
    path.Prepend(usDir);
    // Add this to the TChain
    usCh.Add(path);
  }
  cout<<"Have made TChain with "<<usCh.GetEntries()<<" entries"<<endl;


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
	  // cout<<goodFile<<endl; 
	  break;
	}
      } // Loop over entries
      // Try and find matching spill within the ustof file
      if (goodFile != "nofile") {
	vector<double> spill = findSpill(firstSpillUnix, firstSpillNs, coinTree1, coinTree2, goodFile);
	if (spill.at(2) != -1.) {
	  globalSpillTime = spill.at(1);
	  unixRunStartDs  = runStart;
	  nsTimeDs        = spill.at(0);
	  stopOutDs       = true;
	  runDs           = file;
	  unixRunStartUs  = 0.;
	  ustofSpillTime  = spill.at(3);
	  nsTimeUs        = spill.at(2);
	  spillTree->Fill();
	  // Start counting with the spill offsets
	  vector<pair< pair<double, double>, pair<double, double> >> spillList = spillCount(spill.at(1), spill.at(3), coinTree1, goodFile);
	}
	// Spill does not match
	else {

	}
      }
      
      // Try with the next spill
      else {

      }

      // Find entry with tSoSd within 3s of given spill
      delete tempcoin1;
      delete coinTree1;
      delete inFile1;
      delete tempcoin2;
      delete coinTree2;
      delete inFile2;
    }  
    cout<<" "<<endl;
  } // Loop over dstof files
} // main

double hitTime(double pmtA, double pmtB) {
  double hitT = min(pmtA, pmtB) - (10. - abs(pmtA - pmtB));
}
// countSpill implementation
vector<pair< pair<double,double>, pair<double, double> >> spillCount(const double lastSpillDs, 
								     const double lastSpillUs, 
								     TTree *dsFile, 
								     TChain chain, 
								     TString ustofFile,
								     int index,
								     double window) 
{
  vector<pair< pair<double,double>, pair<double, double> >> spills;
  // Open ustof file
  // TString path = fileName;
  // path.Prepend(usDir);
  // TFile *ustofFile = new TFile(path, "read");
  // TTree *tree = (TTree*)ustofFile->Get("tree");
  // // Get start and end times
  // TNamed *start = 0;
  // ustofFile->GetObject("start_of_run", start);
  // const char* startchar = start->GetTitle();
  // string startstr(startchar);
  // string unixstart = startstr.substr(25,10);
  // double ustofStart = stod(unixstart);	 
  // TNamed *end = 0;
  // ustofFile->GetObject("end_of_run", end);
  // const char* endchar = end->GetTitle();
  // string endstr(endchar);
  // string unixend = endstr.substr(23,10);
  // double ustofEnd = stod(unixend);

  // double tSoSd;
  // tree->SetBranchAddress("tSoSd", &tSoSd);

  double offset = lastSpillDs - lastSpillUs;
  // Go through the dstof file and find the next spill
  RawDsTofCoincidence *tempcoin = NULL;
  dsFile->SetBranchAddress("tofCoin", &tempcoin);
  // Get run start
  dsFile->GetEntry(0);
  int runStart = tempcoin->unixTime[0];

  int lastj = index;
  // Go and find next dstof signal
  double lastSpillNsDs = 0.;
  double lastSpillNsUs = 0.;
  for (int i=0; i<dsFile->GetEntries(); i++) {
    dsFile->GetEntry(i);
    double lastSpillUnixDs = runStart + (tempcoin->lastDelayedBeamSignal / 1e9);
    // If the ustof file we are going through has changed we need to stop and re-pin
    usCh->GetEntry(lastj);
    TString currentFile = usCh->GetCurrentFile()->GetName();
    if (currentFile != ustofFile) {
      break;
    }
    else if (lastSpillUnixDs > lastSpillDs && tempcoin->lastDelayedBeamSignal - tempcoin->lastRawBeamSignal > 0.5e9 && tempcoin->lastDelayedBeamSignal - tempcoin->lastRawBeamSignal < 1e9 && lastSpillNsDs != tempcoin->lastDelayedBeamSignal) {
      lastSpillNsDs = tempcoin->lastDelayedBeamSignal;
      // We have found a dstof spill signal
      for (int j=lastj; j < chain->GetEntries(); j++) {
	chain->GetEntry(j);
	double lastSpillUnixUs = ustofStart + (tSoSd / 1e9);
	// Find spill signal 
	if (lastSpillUnixUs > lastSpillUs && lastSpillUnixUs+offset >= lastSpillUnixDs - window && lastSpillUnixUs + offset <= lastSpillUnixDs + window && lastSpillNsUs != tSoSd && tSoSd - lastSpillNsUs > 1e9) {
	  lastSpillNsUs = tSoSd;
	  cout<<"Spill matched. Ustof "<<lastSpillUnixUs<<", Dstof "<<lastSpillUnixDs<<", Offset "<<(lastSpillUnixDs - lastSpillUnixUs)<<endl;
	  pair<double, double> dstofPair = make_pair(tempcoin->lastDelayedBeamSignal, lastSpillUnixDs);
	  pair<double, double> ustofPair = make_pair(tSoSd, lastSpillUnixUs);
	  pair <pair<double, double>, pair<double, double>> tempPair = make_pair(dstofPair, ustofPair);
	  spills.push_back(tempPair);

	  lastj = j; // We're always going in ascending order
	  break;     // Once we've found one we don't need to go any further
	}
      }
    }
  }
	 
  return spills;

}

// findSpill implementation
vector<double> findSpill(const double timeDs, const double timeNsDs, 
			 TTree *dsFile1, TTree *dsFile2, 
			 TChain fileName, double window) {
  vector<double> spillTimes (4);
  cout<<"Reserved ok"<<endl;
  spillTimes.at(0) = timeNsDs;
  spillTimes.at(1) = timeDs;
  spillTimes.at(2) = -1.;
  spillTimes.at(3) = -1.;
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
  vector<double> spills;
  vector<double> spillsNs;
  // Loop over entries
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
    // Find tSoSd that match within window
    double ustofUnix = ustofStart + (tSoSd / 1e9);
    if (ustofUnix >= (timeDs - window) && ustofUnix <= (timeDs + window) && ustofUnix != lastUnix) {
      spillTimes.at(3) = ustofUnix;
      nMatch++;
      lastUnix = ustofUnix;
      // Make a vector of these spills
      spills.push_back(ustofUnix);
      spillsNs.push_back(tSoSd);
    }
  }
  cout<<nMatch<<" matching spills"<<endl;
  // Get dstof hits in both TDCs
  RawDsTofCoincidence *tempcoin1 = NULL;
  RawDsTofCoincidence *tempcoin2 = NULL;
  dsFile1->SetBranchAddress("tofCoin", &tempcoin1);
  dsFile2->SetBranchAddress("tofCoin", &tempcoin2);
  vector<double> dsHits1;
  vector<double> dsHits2;
  for (int ds1=0; ds1<dsFile1->GetEntries(); ds1++) {
    dsFile1->GetEntry(ds1);
    // In spill
    if ((tempcoin1->fakeTimeNs[0] - timeNsDs) < 1e9 && tempcoin1->fakeTimeNs[0] > timeNsDs) {
      //      dsHits1.push_back(hitTime(tempcoin1->fakeTimeNs[0]-timeNsDs, tempcoin1->fakeTimeNs[1]-timeNsDs));
      dsHits1.push_back(tempcoin1->fakeTimeNs[0] - timeNsDs);
    }
  }
  for (int ds2=0; ds2<dsFile2->GetEntries(); ds2++) {
    dsFile2->GetEntry(ds2);
    // In spill
    if ((tempcoin2->fakeTimeNs[0] - timeNsDs) < 1e9 && tempcoin2->fakeTimeNs[0] > timeNsDs) {
      //      dsHits2.push_back(hitTime((tempcoin2->fakeTimeNs[0]-timeNsDs), (tempcoin2->fakeTimeNs[1]-timeNsDs)));
      dsHits2.push_back(tempcoin2->fakeTimeNs[0] - timeNsDs);
    }
  }
  // See if the hits match between ustof and dstof
  for (int s=0; s < spills.size(); s++) {
    // Get ustof hits
    vector<double> usHits;
    TGraph* gr1 = new TGraph();
    TH1D *hds = new TH1D("hds", "", 400, 0, 1e9);
    TH1D *hus = new TH1D("hus", "", 400, 0, 1e9);
    hds->SetLineColor(kMagenta);
    hds->SetFillColor(kMagenta);
    hds->SetFillStyle(3005);
    int nMatchedHits = 0;
    for (int i = 0; i<tree->GetEntries(); i++) {
      tree->GetEntry(i);
      // In spill
      if (tToF[0] > spillsNs[s] && tToF[0] - spillsNs[s] < 1e9) {
	for (int h=0; h<nhit; h++) {
	  // Use crude drift correction. Should work...
	  usHits.push_back((tToF[h] - spillsNs[s]) / (1. + stdDrift));	  
	} // hits
      }
    } // ustof tree

    // Match each ustof hit to a dstof hit
    for (int us=0; us<usHits.size(); us++) {
      double diffLow1 = 9999999999.;
      double diff1    = 9999999999.;
      double diffLow2 = 9999999999.;
      double diff2    = 9999999999.;
      double lowest   = 9999999999999.;
      int tdc = 0;
      hus->Fill(usHits[us]);
      for (int ds=0; ds<dsHits1.size(); ds++) {
	diff1 = usHits[us] - dsHits1[ds];
	hds->Fill(dsHits1[ds]);
	if (abs(diff1) < abs(diffLow1)) {
	  diffLow1 = diff1;
	  tdc = 1;
	}
      } // ds1 hits
      for (int ds=0; ds<dsHits2.size(); ds++) {
	diff2 = usHits[us] - dsHits2[ds];
	hds->Fill(dsHits2[ds]);
	if (abs(diff2) < abs(diffLow2)) {
	  diffLow2 = diff2;
	  tdc = 2;
	}
      } // ds2 hits
      if(abs(diffLow1) < abs(diffLow2)) {
	lowest = diffLow1;
      }
      else {
	lowest = diffLow2;
      }
      gr1->SetPoint(gr1->GetN(), usHits[us], lowest);
      // Put the cut on number of matched hits required to be matched at 75
      // Count the number of matched (<500ns separation) hits
      if (abs(lowest) <= 750) {
	nMatchedHits++;
      }
    }
    cout<<nMatchedHits<<" matched hits: ";
    
    if(nMatchedHits >= 70) {
      spillTimes.at(2) = spillsNs[s];
      spillTimes.at(3) = spills[s];
      cout<<"PASSED"<<endl;
    }
    else {
      cout<<"FAILED"<<endl;
    }
    /*
    new TCanvas;
    hds->Draw("hist");
    hus->Draw("hist");
    gPad->Print("test.pdf");
    gPad->Print("test.png");
    new TCanvas;
    gr1->GetYaxis()->SetRangeUser(-2000, 2000);
    gr1->Draw("AP*");
    gPad->Print("gr.png");
    */
  } // spills

  delete tree;
  delete ustofFile;

  return spillTimes;
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

  double entries = 2.;

  while (entry = (char*)gSystem->GetDirEntry(dirp)) { 
    str = entry;
    if (str.EndsWith(ext) && !str.Contains(bad) && str.Contains(pref)) {
      TString strTemp = str;
      str.Prepend(indir);

      TFile *ustofFile = new TFile(str, "read");
      TTree *tree = (TTree*)ustofFile->Get("tree");
      int subentries = tree->GetEntries();
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
      vector<double> ustofTmp =  {ustofStart, ustofEnd, entries};
      fileVec.push_back(make_pair(ustofTmp, strTemp));
      entries += subentries;
      delete start;
      delete end;
      delete tree;
      delete ustofFile;
    }
  }
  sort(begin(fileVec), end(fileVec));
  cout<<"Files: "<<fileVec.size()<<endl;
  /*
  for (int i=0; i<fileVec.size(); i++) {
    cout<<fileVec[i].first[0]<<" "<<fileVec[i].first[1]<<" "<<fileVec[i].second<<endl;
  }
  */
  return fileVec;
} // ustofFileVec
