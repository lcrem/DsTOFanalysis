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
  Int_t itdc;
  string baseDir;
  
  if((argc!=3)&&(argc!=4)){
    std::cerr << "Usage 1: " << argv[0] << " [irun] [itdc] (basedir)" << std::endl;
    return 1;
  } else {
    run = atoi(argv[1]);
    itdc = atoi(argv[2]);
    if (argc==4) baseDir += argv[3];
    else  baseDir +=  whereIsMyTOFdata;
  }

  string dirname = Form("%s/run%d/", baseDir.c_str(), run);

  string filename[2];
  filename[itdc] = findFilename(dirname, itdc+1);
  
  if (filename[0]=="nothing" || filename[1]=="nothing"){
    cout << "Could not find parsed files in directory " << dirname << endl;
  } else {
    cout << "Path to TDC parsed file " << filename[0] << endl;
  }

  string line;

  ifstream f((filename[itdc]).c_str());


  Int_t channel1;
  Int_t ticks1;
  Int_t countClock1=0;
  Int_t beamSpill1;
  Int_t unixTime1;
  Int_t usTof1;
  Double_t triggerTimeNs1;
  Double_t fakeTimeNs1;
  Int_t oldticks=0;
    
  string temp;
  RawDsTofHeader *tof = NULL;

  TFile *fout = new TFile(Form("%s/DsTOFtreeRun%d_tdc%d.root", dirname.c_str(), run, itdc+1), "recreate");
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

      if (channel1==0){
	if (ticks1<oldticks)
	  countClock1++;
	oldticks=ticks1;
      } else {

	// if (channel1==15) beamSpill1 = ticks1;
	// if (channel1==13) usTof1 = ticks1;
	// triggerTimeNs1 = ticks1*clockTicksNs;
	fakeTimeNs1    = ( countClock1*TMath::Power(2, 21) +ticks1)*clockTicksNs;
	// cout << run << " " << channel1 << endl;
	  
	if(tof)
	  delete tof;

	tof = new RawDsTofHeader();
	tof->run=(Short_t)run;	  
	tof->tdc=(Short_t)itdc+1; 
	tof->channel=(Short_t)channel1; 
	tof->ticks=ticks1; 
	tof->clockCounter=countClock1; 
	tof->unixTime=unixTime1; 
	tof->fakeTimeNs=fakeTimeNs1; 
	tofTree->Fill();
      }
    }

    cout << "Writing output " << Form("%s/DsTOFtreeRun%d_tdc%d.root", dirname.c_str(), run, itdc+1) << endl;
    tofTree->Write("tofTree");
    fout->Close();

    f.close();
  }else{
    cout << "Don't know file : " << filename[itdc] << endl;
  }

  cout << "DONE TDC " << itdc+1 << endl;
  
    
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
