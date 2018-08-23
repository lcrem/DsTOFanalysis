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


string findFilename(string dirname, int itdc);
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
    filename[itdc] = findFilename(dirname, itdc+1);
  }
  
  if (filename[0]=="nothing" || filename[1]=="nothing"){
    cout << "Could not find parsed files in directory " << dirname << endl;
  } else {
    cout << "Path to TDC1 parsed file " << filename[0] << endl;
    cout << "Path to TDC2 parsed file " << filename[1] << endl;
  }

  string line;

  for (int itdc=0; itdc<2; itdc++){

    ifstream f((filename[itdc]).c_str());


    Int_t channel1;
    Int_t ticks1;
    UInt_t unixTime1;
    Double_t fakeTimeNs1;    
    Double_t lastFakeTimeNs1;
    Int_t count25k=0;
    Int_t lastClockTdc=0;
    Int_t tdcDiff;
    
    string temp;
    RawDsTofHeader *tof = NULL;

    cout << "Writing output " << Form("%s/DsTOFtreeRun%d_tdc%d.root", dirname.c_str(), run, itdc+1) << endl;
    TFile *fout = new TFile(Form("%s/DsTOFtreeRun%d_tdc%d.root", dirname.c_str(), run, itdc+1), "recreate");
    TH1D *hClock = new TH1D("hClock", "", 100, -TMath::Power(2,21), TMath::Power(2,21));
    TTree *tofTree = new TTree("tofTree","HPTPC Time Of Flight");
    tofTree->Branch("tof", &tof);
  
    if (f.is_open()){

      while (getline (f,line)){
	if (line.find("unix") != std::string::npos){
	  stringstream ss(line);
	  ss >> temp >> unixTime1;
	  continue;
	}
    
	stringstream ss(line);
	ss >> channel1 >> temp >> ticks1;

	if(channel1==0) {
	  count25k++;
	  hClock->Fill(ticks1-lastClockTdc);
	  lastClockTdc=ticks1;
	  fakeTimeNs1=count25k*40e3;
	}  else {

	  tdcDiff=ticks1-lastClockTdc;

	  if(tdcDiff<0) tdcDiff+=TMath::Power(2,21);

	  fakeTimeNs1=count25k*40e3  + tdcDiff*clockTicksNs;
	  //	  fakeTimeNs1    = ( countClock1*TMath::Power(2, 21) +ticks1)*clockTicksNs;

	  // cout << run << " " << channel1 << endl;
	  
	  if(tof)
	    delete tof;

	  tof = new RawDsTofHeader();
	  tof->run=(Short_t)run;	  
	  tof->tdc=(Short_t)itdc+1; 
	  tof->channel=(Short_t)channel1; 
	  tof->ticks=ticks1; 
	  tof->count25k=count25k; 
	  tof->unixTime=unixTime1; 
	  tof->fakeTimeNs=fakeTimeNs1; 
	  tofTree->Fill();
	}
      }

      tofTree->Write("tofTree");
      hClock->Write("hClock");
      fout->Close();

      f.close();
    }else{
      cout << "Don't know file : " << filename[itdc] << endl;
    }

    cout << "DONE TDC " << itdc+1 << endl;
  }
    
  return 0;
}


string findFilename(string dirname, int itdc){
  
  DIR *dp;
  dirent *d;
  
  
  if((dp = opendir(dirname.c_str())) == NULL)
    perror("opendir");

  while((d = readdir(dp)) != NULL){
    if(!strcmp(d->d_name,".") || !strcmp(d->d_name,".."))
      continue;
    
    string temp = d->d_name;
    if (temp.find(Form("parsed_%d_", itdc)) != std::string::npos){
      if (temp.find(".root") != std::string::npos) continue;
      else return dirname+"/"+temp;
      
    }
  }  
  return "nothing";
}
