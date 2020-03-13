// mcAcrylicComp.C
// Look at files with addtiional 1cm piece of acrylic
#include "UsefulFunctions.C"

const char* stem = "/scratch0/tnonnenm/FixedMC/Prototype-HPTPC-MC/wdir/moderator_data";

void mcAcrylicComp(/*const char* outfile*/)
{
  gROOT->SetBatch(kTRUE);

  //  TFile *fout = new TFile(outfile, "recreate");

  for (int b=0; b<5; b++) {
    cout<<b<<" blocks"<<endl;
    TFile *fin  = new TFile(Form("%s/NewGeom%dBlocks.root", stem, b), "read");
    TTree *tree = (TTree*)fin->Get("stepTree");
    UInt_t NSteps;
    float PostStepMom[40000];
    float PostStepPosX[40000];
    float PostStepPosY[40000];
    float PostStepPosZ[40000];
    int PreStepVolumeID[40000];
    int PostStepVolumeID[40000];
    int PDG[40000];
    tree->SetBranchAddress("NSteps", &NSteps);
    tree->SetBranchAddress("PostStepMom", PostStepMom);
    tree->SetBranchAddress("PostStepPosX", PostStepPosX);
    tree->SetBranchAddress("PostStepPosY", PostStepPosY);
    tree->SetBranchAddress("PostStepPosZ", PostStepPosZ);
    tree->SetBranchAddress("PreStepVolumeID", PreStepVolumeID);
    tree->SetBranchAddress("PostStepVolumeID", PostStepVolumeID);
    tree->SetBranchAddress("PDG", PDG);

    int nS3 = tree->GetEntries();
    int nS4 = 0;

    for (int i=0; i<tree->GetEntries(); i++) {
      tree->GetEntry(i);
      // Loop over steps
      for (int n=0; n<NSteps; n++) {
	// Condition for a proton being in S4
	if (PreStepVolumeID[n]!=7 && PostStepVolumeID[n]==7 && PostStepMom[n] > 30. 
	    && PDG[n]==2212) {
	  nS4++;
	}
      }
    } 
    fin->Close();
    delete fin;

    // Now with the extra acrylic
    TFile *finAcr  = new TFile(Form("%s/%dblock1.0cmacrylic.root", stem, b), "read");
    TTree *treeAcr = (TTree*)finAcr->Get("stepTree");
    treeAcr->SetBranchAddress("NSteps", &NSteps);
    treeAcr->SetBranchAddress("PostStepMom", PostStepMom);
    treeAcr->SetBranchAddress("PostStepPosX", PostStepPosX);
    treeAcr->SetBranchAddress("PostStepPosY", PostStepPosY);
    treeAcr->SetBranchAddress("PostStepPosZ", PostStepPosZ);
    treeAcr->SetBranchAddress("PreStepVolumeID", PreStepVolumeID);
    treeAcr->SetBranchAddress("PostStepVolumeID", PostStepVolumeID);
    treeAcr->SetBranchAddress("PDG", PDG);

    int nS4Acr = 0;

    for (int i=0; i<treeAcr->GetEntries(); i++) {
      treeAcr->GetEntry(i);
      // Loop over steps
      for (int n=0; n<NSteps; n++) {
	// Condition for a proton being in S4
	if (PreStepVolumeID[n]!=7 && PostStepVolumeID[n]==7 && PostStepMom[n] > 30. 
	    && PDG[n]==2212) {
	  nS4Acr++;
	}
      }
    } 
    finAcr->Close();
    delete finAcr;

    double errS3    = sqrt(nS3);
    double errS4    = sqrt(nS4);
    double errS4Acr = sqrt(nS4Acr);

    double r    = (double)nS4/(double)nS3;
    double rAcr = (double)nS4Acr/(double)nS3;
    double rErr    = r * sqrt( pow(errS4/(double)nS4, 2) + pow(errS3/(double)nS3, 2) );
    double rAcrErr = rAcr * sqrt( pow(errS4Acr/(double)nS4Acr, 2) + pow(errS3/(double)nS3, 2) );
    cout<<"Old ratio "<<r<<" +- "<<rErr<<endl;
    cout<<"With acrylic: "<<rAcr<<" +- "<<rAcrErr<<endl;
    double diff = (r-rAcr)/r * 100.;
    cout<<"%age difference "<<diff<<endl;
  } // Loop over blocks
} // mcAcrylicComp
