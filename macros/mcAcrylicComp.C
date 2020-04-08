// mcAcrylicComp.C
// Look at files with addtiional 1cm piece of acrylic
#include "UsefulFunctions.C"

const char* stem = "/scratch0/tnonnenm/FixedMC/Prototype-HPTPC-MC/wdir/moderator_data";

void mcAcrylicComp(const char* outfile)
{
  gROOT->SetBatch(kTRUE);

  TFile *fout = new TFile(outfile, "recreate");

  TVector3 mcS4_TL = globalToMCCoords(S4_TL);
  TVector3 mcS4_BR = globalToMCCoords(S4_BR);
  double x1 = mcS4_BR.X();
  double x2 = mcS4_TL.X();
  double y1 = mcS4_BR.Y();
  double y2 = mcS4_TL.Y();

  cout<<"x1, x2 "<<x1<<", "<<x2<<endl;
  cout<<"y1, y2 "<<y1<<", "<<y2<<endl;

  TH1D *hDiff1cm = new TH1D("hDiff1cm", "1cm acryclic; N. blocks; |Fractional difference|", 5, -0.5, 4.5);
  setHistAttr(hDiff1cm);
  hDiff1cm->SetLineColor(kBlue);
  TH1D *hDiff2cm = new TH1D("hDiff2cm", "2cm acryclic; N. blocks; |Fractional difference|", 5, -0.5, 4.5);
  setHistAttr(hDiff2cm);
  hDiff2cm->SetLineColor(kRed);
  THStack *hsDiff = new THStack("hsDiff", "Fractional change in S4 proton count due to additional acrylic; N. blocks; |Fractional difference|");

  for (int b=0; b<5; b++) {
    cout<<b<<" blocks"<<endl;
    TFile *fin  = new TFile(Form("%s/NewGeom%dBlocks.root", stem, b), "read");
    TTree *tree = (TTree*)fin->Get("stepTree");
    UInt_t NSteps;

    TH2D *hNom = new TH2D(Form("hNom%d", b),"No acrylic; X_{MC}; Y_{MC}", 40, -0.7518, 0.5496, 23, -0.29995, 0.48555);
    setHistAttr(hNom);
    TH2D *h1cm = new TH2D(Form("h1cm%d", b),"1 cm acrylic; X_{MC}; Y_{MC}", 40, -0.7518, 0.5496, 23, -0.29995, 0.48555);
    setHistAttr(h1cm);
    TH2D *h2cm = new TH2D(Form("h2cm%d", b),"2 cm acrylic; X_{MC}; Y_{MC}", 40, -0.7518, 0.5496, 23, -0.29995, 0.48555);
    setHistAttr(h2cm);

    int nSpills = 0;

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
    tree->AddFriend(Form("protonTree%dBlocks", b), "/scratch0/sjones/plots/angularDistS3_newSample/s3ProtonWeightTrees.root");
    int spill;
    tree->SetBranchAddress("spill", &spill);
    tree->GetEntry(tree->GetEntries()-1);
    nSpills = spill;

    int nS3 = tree->GetEntries();
    int nS4 = 0;

    cout<<"Looping"<<endl;
    for (int i=0; i<tree->GetEntries(); i++) {
      tree->GetEntry(i);
      // Loop over steps
      for (int n=0; n<NSteps; n++) {
	// Condition for a proton being in S4
	if (PreStepVolumeID[n]!=7 && PostStepVolumeID[n]==7 && PostStepMom[n] > 30. 
	    && PDG[n]==2212) {
	  nS4++;
	  hNom->Fill(PostStepPosX[n]/1000, PostStepPosY[n]/1000.);
	}
      }
    } 
    cout<<"Old file done"<<endl;

    // Now with the extra acrylic
    TFile *finAcr  = new TFile(Form("%s/%dblock1.0cmacrylic.root", stem, b), "read");
    TTree *treeAcr = (TTree*)finAcr->Get("stepTree");
    UInt_t NSteps1;
    float PostStepMom1[40000];
    float PostStepPosX1[40000];
    float PostStepPosY1[40000];
    float PostStepPosZ1[40000];
    float PreStepPosX1[40000];
    float PreStepPosY1[40000];
    float PreStepPosZ1[40000];
    int PreStepVolumeID1[40000];
    int PostStepVolumeID1[40000];
    int PDG1[40000];
    treeAcr->SetBranchAddress("NSteps", &NSteps1);
    treeAcr->SetBranchAddress("PostStepMom", PostStepMom1);
    treeAcr->SetBranchAddress("PostStepPosX", PostStepPosX1);
    treeAcr->SetBranchAddress("PostStepPosY", PostStepPosY1);
    treeAcr->SetBranchAddress("PostStepPosZ", PostStepPosZ1);
    treeAcr->SetBranchAddress("PreStepPosX", PreStepPosX1);
    treeAcr->SetBranchAddress("PreStepPosY", PreStepPosY1);
    treeAcr->SetBranchAddress("PreStepPosZ", PreStepPosZ1);
    treeAcr->SetBranchAddress("PreStepVolumeID", PreStepVolumeID1);
    treeAcr->SetBranchAddress("PostStepVolumeID", PostStepVolumeID1);
    treeAcr->SetBranchAddress("PDG", PDG1);

    int nS4Acr = 0;

    cout<<"Doing 1cm file"<<endl;
    for (int e=0; e<treeAcr->GetEntries(); e++) {
      treeAcr->GetEntry(e);
      //cout<<NSteps1<<endl;
      for (int n1=0; n1<NSteps1; n1++) {
	// Condition for a proton being in S4
	if (PreStepVolumeID1[n1]!=7 && PostStepVolumeID1[n1]==7 && PostStepMom1[n1] > 30. 
	    && PDG1[n1]==2212) {
	  // cout<<n1<<endl;
	  // cout<<(PostStepPosX1[n1]/1000.)<<", "<<(PostStepPosY1[n1]/1000.)<<endl;
	  h1cm->Fill((double)(PostStepPosX1[n1]/1000.), (double)(PostStepPosY1[n1]/1000.));
	  nS4Acr++;
	}
      }
    }
    cout<<"1 cm done"<<endl;

    TFile *fin2Acr  = new TFile(Form("%s/%dblock2.0cmacrylic.root", stem, b), "read");
    TTree *tree2Acr = (TTree*)fin2Acr->Get("stepTree");
    UInt_t NSteps2;
    float PostStepMom2[40000];
    float PostStepPosX2[40000];
    float PostStepPosY2[40000];
    float PostStepPosZ2[40000];
    int PreStepVolumeID2[40000];
    int PostStepVolumeID2[40000];
    int PDG2[40000];
    tree2Acr->SetBranchAddress("NSteps", &NSteps2);
    tree2Acr->SetBranchAddress("PostStepMom", PostStepMom2);
    tree2Acr->SetBranchAddress("PostStepPosX", PostStepPosX2);
    tree2Acr->SetBranchAddress("PostStepPosY", PostStepPosY2);
    tree2Acr->SetBranchAddress("PostStepPosZ", PostStepPosZ2);
    tree2Acr->SetBranchAddress("PreStepVolumeID", PreStepVolumeID2);
    tree2Acr->SetBranchAddress("PostStepVolumeID", PostStepVolumeID2);
    tree2Acr->SetBranchAddress("PDG", PDG2);

    int nS42Acr = 0;

    cout<<"Doing 2cm"<<endl;
    for (int i=0; i<tree2Acr->GetEntries(); i++) {
      tree2Acr->GetEntry(i);
      // Loop over steps
      for (int n=0; n<NSteps2; n++) {
    	// Condition for a proton being in S4
    	if (PreStepVolumeID2[n]!=7 && PostStepVolumeID2[n]==7 && PostStepMom2[n] > 30. 
    	    && PDG2[n]==2212) {
    	  nS42Acr++;
    	  h2cm->Fill(PostStepPosX2[n]/1000., PostStepPosY2[n]/1000.);
    	}
      }
    } 
    cout<<"2 cm done"<<endl;

    double errS3     = sqrt(nS3);
    double errS4     = sqrt(nS4);
    double errS4Acr  = sqrt(nS4Acr);
    double errS42Acr = sqrt(nS42Acr);

    double r    = (double)nS4/(double)nS3;
    double rAcr = (double)nS4Acr/(double)nS3;
    double r2Acr = (double)nS42Acr/(double)nS3;
    double rErr     = r * sqrt( pow(errS4/(double)nS4, 2) + pow(errS3/(double)nS3, 2) );
    double rAcrErr  = rAcr * sqrt( pow(errS4Acr/(double)nS4Acr, 2) + pow(errS3/(double)nS3, 2) );
    double r2AcrErr = r2Acr * sqrt( pow(errS42Acr/(double)nS42Acr, 2) + pow(errS3/(double)nS3, 2) );
    cout<<"Old ratio "<<r<<" +- "<<rErr<<endl;
    cout<<"With 1cm acrylic: "<<rAcr<<" +- "<<rAcrErr<<endl;
    double diff = (r-rAcr)/r * 100.;
    cout<<"%age difference "<<diff<<endl;
    cout<<"With 2cm acrylic: "<<r2Acr<<" +- "<<r2AcrErr<<endl;
    double diff2 = (r-r2Acr)/r * 100.;
    cout<<"%age difference "<<diff2<<endl;
    cout<<"================================"<<endl;

    hDiff1cm->Fill(b, abs(((double)nS4-(double)nS4Acr)/(double)nS4));
    hDiff2cm->Fill(b, abs(((double)nS4-(double)nS42Acr)/(double)nS4));

    hNom->Scale(1./nSpills);
    h1cm->Scale(1./nSpills);
    h2cm->Scale(1./nSpills);
    fout->cd();
    hNom->Write();
    h1cm->Write();
    h2cm->Write();

    fin->Close();
    delete fin;
    finAcr->Close();
    delete finAcr;
    fin2Acr->Close();
    delete fin2Acr;
  } // Loop over blocks

  fout->cd();
  hDiff1cm->Write();
  hDiff2cm->Write();
  hsDiff->Add(hDiff1cm);
  hsDiff->Add(hDiff2cm);
  hsDiff->Write();

  fout->Close();
  delete fout;
} // mcAcrylicComp
