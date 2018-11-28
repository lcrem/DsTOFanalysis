// makeSpillDB.cxx

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

int main() {

  // The final word in beam spill databases
  TTree *spillTree = new TTree("spillTree", "HPTPC beam spills");
  spillTree->SetDirectory(0);

  // All times are in seconds unless otherwise stated
  double globalSpillTime; // Define the dtof to be the global spill time
  double nsTimeDs;        // Time in ns since the dstof run start
  int runDs;              // Run number of this dstof spill 
  double ustofSpillTime; 
  double nsTimeUs;        // Time in ns since the ustof run start


}
