// s1s2MCSysts.C
// Attempt to figure out systematics on MC through moving S1 and S2 positions
#include "UsefulFunctions.C"

const char* stem = "/scratch0/tnonnenm/FixedMC/Prototype-HPTPC-MC/wdir/moderator_data/";

void s1s2MCSysts(const char* outfile)
{
  gROOT->SetBatch(kTRUE);

  TFile *fout = new TFile(outfile, "recreate");

  vector<double> nS4 = {0., 0., 0., 0., 0.};
  vector<double> nS4_S1Only = {0., 0., 0., 0., 0.};
  vector<double> nS4_S2Plus5X = {0., 0., 0., 0., 0.};
  vector<double> nS4_S2Minus5X = {0., 0., 0., 0., 0.};
  vector<double> nS4_S1Plus2XY = {0., 0., 0., 0., 0.};
  vector<double> nS4_S1Minus2XY = {0., 0., 0., 0., 0.};
  vector<double> nS3 = {0., 0., 0., 0., 0.};

  for (int b=0; b<5; b++) {
    cout<<"==============================="<<endl;
    cout<<b<<" blocks"<<endl;
    cout<<"==============================="<<endl;

    TFile *fin   = new TFile(Form("%s/NewGeom%dBlocks.root", stem, b), "read");
    TTree *t = (TTree*)fin->Get("stepTree");
    UInt_t NSteps;
    float Mom[50000];
    int PreVolID[50000];
    int PostVolID[50000];
    int PDG[50000];
    t->SetBranchAddress("NSteps", &NSteps);
    t->SetBranchAddress("PostStepMom", Mom);
    t->SetBranchAddress("PreStepVolumeID", PreVolID);
    t->SetBranchAddress("PostStepVolumeID", PostVolID);
    t->SetBranchAddress("PDG", PDG);
    nS3.at(b) = t->GetEntries();
    cout<<"Looping through original MC file"<<endl;
    for(int e=0; e<t->GetEntries(); e++) {
      t->GetEntry(e);
      for (UInt_t n=0; n<NSteps; n++) {
	if (PreVolID[n]!=7 && PostVolID[n]==7 && Mom[n] > 140. && PDG[n]==2212) {
	  nS4.at(b)++;
	  break;
	}
      }
    }
    fin->Close();
    delete fin;

    TFile *finS1Only = new TFile(Form("%s/%dblockS1only.root", stem, b), "read");
    TTree *tS1Only = (TTree*)finS1Only->Get("stepTree");
    UInt_t NSteps_S1Only;
    float Mom_S1Only[50000];
    int PreVolID_S1Only[50000];
    int PostVolID_S1Only[50000];
    int PDG_S1Only[50000];
    tS1Only->SetBranchAddress("NSteps", &NSteps_S1Only);
    tS1Only->SetBranchAddress("PostStepMom", Mom_S1Only);
    tS1Only->SetBranchAddress("PreStepVolumeID", PreVolID_S1Only);
    tS1Only->SetBranchAddress("PostStepVolumeID", PostVolID_S1Only);
    tS1Only->SetBranchAddress("PDG", PDG_S1Only);
    cout<<"Looping through S1 only file"<<endl;
    for (int e=0; e<tS1Only->GetEntries(); e++) {
      tS1Only->GetEntry(e);
      for (UInt_t n=0; n<NSteps_S1Only; n++) {
	if (PreVolID_S1Only[n]!=7 && PostVolID_S1Only[n]==7 && Mom_S1Only[n] > 140. && 
	    PDG_S1Only[n]==2212) {
	  nS4_S1Only.at(b)++;
	  break;
	}
      } // Loop over steps
    }
    finS1Only->Close();
    delete finS1Only;

    TFile *finS2Plus5X = new TFile(Form("%s/%dblockS2+5X.root", stem, b), "read");
    TTree *tS2Plus5X = (TTree*)finS2Plus5X->Get("stepTree");
    UInt_t NSteps_S2Plus5X;
    float Mom_S2Plus5X[50000];
    int PreVolID_S2Plus5X[50000];
    int PostVolID_S2Plus5X[50000];
    int PDG_S2Plus5X[50000];
    tS2Plus5X->SetBranchAddress("NSteps", &NSteps_S2Plus5X);
    tS2Plus5X->SetBranchAddress("PostStepMom", Mom_S2Plus5X);
    tS2Plus5X->SetBranchAddress("PreStepVolumeID", PreVolID_S2Plus5X);
    tS2Plus5X->SetBranchAddress("PostStepVolumeID", PostVolID_S2Plus5X);
    tS2Plus5X->SetBranchAddress("PDG", PDG_S2Plus5X);
    cout<<"Looping through S2 plus 5X file"<<endl;
    for (int e=0; e<tS2Plus5X->GetEntries(); e++) {
      tS2Plus5X->GetEntry(e);
      for (UInt_t n=0; n<NSteps_S2Plus5X; n++) {
	if (PreVolID_S2Plus5X[n]!=7 && PostVolID_S2Plus5X[n]==7 && Mom_S2Plus5X[n] > 140. && 
	    PDG_S2Plus5X[n]==2212) {
	  nS4_S2Plus5X.at(b)++;
	  break;
	}
      } // Loop over steps
    }
    cout<<"S1 only file done"<<endl;
    finS2Plus5X->Close();
    delete finS2Plus5X;

    if (b != 1) {
      TFile *finS2Minus5X = new TFile(Form("%s/%dblockS2-5X.root", stem, b), "read");
      TTree *tS2Minus5X = (TTree*)finS2Minus5X->Get("stepTree");
      UInt_t NSteps_S2Minus5X;
      float Mom_S2Minus5X[50000];
      int PreVolID_S2Minus5X[50000];
      int PostVolID_S2Minus5X[50000];
      int PDG_S2Minus5X[50000];
      tS2Minus5X->SetBranchAddress("NSteps", &NSteps_S2Minus5X);
      tS2Minus5X->SetBranchAddress("PostStepMom", Mom_S2Minus5X);
      tS2Minus5X->SetBranchAddress("PreStepVolumeID", PreVolID_S2Minus5X);
      tS2Minus5X->SetBranchAddress("PostStepVolumeID", PostVolID_S2Minus5X);
      tS2Minus5X->SetBranchAddress("PDG", PDG_S2Minus5X);
      cout<<"Looping through S2 minus 5X file"<<endl;
      for (int e=0; e<tS2Minus5X->GetEntries(); e++) {
	tS2Minus5X->GetEntry(e);
	for (UInt_t n=0; n<NSteps_S2Minus5X; n++) {
	  if (PreVolID_S2Minus5X[n]!=7 && PostVolID_S2Minus5X[n]==7 && Mom_S2Minus5X[n] > 140. && 
	      PDG_S2Minus5X[n]==2212) {
	    nS4_S2Minus5X.at(b)++;
	    break;
	  }
	} // Loop over steps
      }
      cout<<"S1 only file done"<<endl;
      finS2Minus5X->Close();
      delete finS2Minus5X;
    }
    else {
      TFile *finS2Minus5X = new TFile(Form("%s/%dblockS2-5XRERUN.root", stem, b), "read");
      TTree *tS2Minus5X = (TTree*)finS2Minus5X->Get("stepTree");
      UInt_t NSteps_S2Minus5X;
      float Mom_S2Minus5X[50000];
      int PreVolID_S2Minus5X[50000];
      int PostVolID_S2Minus5X[50000];
      int PDG_S2Minus5X[50000];
      tS2Minus5X->SetBranchAddress("NSteps", &NSteps_S2Minus5X);
      tS2Minus5X->SetBranchAddress("PostStepMom", Mom_S2Minus5X);
      tS2Minus5X->SetBranchAddress("PreStepVolumeID", PreVolID_S2Minus5X);
      tS2Minus5X->SetBranchAddress("PostStepVolumeID", PostVolID_S2Minus5X);
      tS2Minus5X->SetBranchAddress("PDG", PDG_S2Minus5X);
      cout<<"Looping through S2 minus 5X file"<<endl;
      for (int e=0; e<tS2Minus5X->GetEntries(); e++) {
	tS2Minus5X->GetEntry(e);
	for (UInt_t n=0; n<NSteps_S2Minus5X; n++) {
	  if (PreVolID_S2Minus5X[n]!=7 && PostVolID_S2Minus5X[n]==7 && Mom_S2Minus5X[n] > 140. && 
	      PDG_S2Minus5X[n]==2212) {
	    nS4_S2Minus5X.at(b)++;
	    break;
	  }
	} // Loop over steps
      }
      cout<<"S1 only file done"<<endl;
      finS2Minus5X->Close();
      delete finS2Minus5X;
    }

    TFile *finS1Plus2XY = new TFile(Form("%s/%dblockS1+2XY.root", stem, b), "read");
    TTree *tS1Plus2XY = (TTree*)finS1Plus2XY->Get("stepTree");
    UInt_t NSteps_S1Plus2XY;
    float Mom_S1Plus2XY[50000];
    int PreVolID_S1Plus2XY[50000];
    int PostVolID_S1Plus2XY[50000];
    int PDG_S1Plus2XY[50000];
    tS1Plus2XY->SetBranchAddress("NSteps", &NSteps_S1Plus2XY);
    tS1Plus2XY->SetBranchAddress("PostStepMom", Mom_S1Plus2XY);
    tS1Plus2XY->SetBranchAddress("PreStepVolumeID", PreVolID_S1Plus2XY);
    tS1Plus2XY->SetBranchAddress("PostStepVolumeID", PostVolID_S1Plus2XY);
    tS1Plus2XY->SetBranchAddress("PDG", PDG_S1Plus2XY);
    cout<<"Looping through S1 plus 2XY file"<<endl;
    for (int e=0; e<tS1Plus2XY->GetEntries(); e++) {
      tS1Plus2XY->GetEntry(e);
      for (UInt_t n=0; n<NSteps_S1Plus2XY; n++) {
	if (PreVolID_S1Plus2XY[n]!=7 && PostVolID_S1Plus2XY[n]==7 && Mom_S1Plus2XY[n] > 140. && 
	    PDG_S1Plus2XY[n]==2212) {
	  nS4_S1Plus2XY.at(b)++;
	  break;
	}
      } // Loop over steps
    }
    finS1Plus2XY->Close();
    delete finS1Plus2XY;

    TFile *finS1Minus2XY = new TFile(Form("%s/%dblockS1-2XY.root", stem, b), "read");
    TTree *tS1Minus2XY = (TTree*)finS1Minus2XY->Get("stepTree");
    UInt_t NSteps_S1Minus2XY;
    float Mom_S1Minus2XY[50000];
    int PreVolID_S1Minus2XY[50000];
    int PostVolID_S1Minus2XY[50000];
    int PDG_S1Minus2XY[50000];
    tS1Minus2XY->SetBranchAddress("NSteps", &NSteps_S1Minus2XY);
    tS1Minus2XY->SetBranchAddress("PostStepMom", Mom_S1Minus2XY);
    tS1Minus2XY->SetBranchAddress("PreStepVolumeID", PreVolID_S1Minus2XY);
    tS1Minus2XY->SetBranchAddress("PostStepVolumeID", PostVolID_S1Minus2XY);
    tS1Minus2XY->SetBranchAddress("PDG", PDG_S1Minus2XY);
    cout<<"Looping through S2 minus 2XY file"<<endl;
    for (int e=0; e<tS1Minus2XY->GetEntries(); e++) {
      tS1Minus2XY->GetEntry(e);
      for (UInt_t n=0; n<NSteps_S1Minus2XY; n++) {
	if (PreVolID_S1Minus2XY[n]!=7 && PostVolID_S1Minus2XY[n]==7 && Mom_S1Minus2XY[n] > 140. && 
	    PDG_S1Minus2XY[n]==2212) {
	  nS4_S1Minus2XY.at(b)++;
	  break;
	}
      } // Loop over steps
    }
    finS1Minus2XY->Close();
    delete finS1Minus2XY;
  } // Loop over blocks

  TH1D *hS1Only    = new TH1D("hS1Only", "MC fractional changes in S4 proton count due to S1 only; Blocks; |Fractional change|", 5, -0.5, 4.5);
  TH1D *hS2Plus5X  = new TH1D("hS2Plus5X", "MC fractional changes in S4 proton count due to +5cm S2 X shift; Blocks; |Fractional change|", 5, -0.5, 4.5);
  TH1D *hS2Minus5X = new TH1D("hS2Minus5X", "MC fractional changes in S4 proton count due to -5cm S2 X shift; Blocks; |Fractional change|", 5, -0.5, 4.5);
  TH1D *hS1Plus2XY = new TH1D("hS1Plus2XY", "MC fractional changes in S4 proton count due to +2cm S1 XY shift; Blocks; |Fractional change|", 5, -0.5, 4.5);
  TH1D *hS1Minus2XY = new TH1D("hS1Minus2XY", "MC fractional changes in S4 proton count due to -2cm S1 XY shift; Blocks; |Fractional change|", 5, -0.5, 4.5);
  setHistAttr(hS1Only);
  setHistAttr(hS2Plus5X);
  setHistAttr(hS2Minus5X);
  setHistAttr(hS1Plus2XY);
  setHistAttr(hS1Minus2XY);
  hS1Only->SetLineColor(kBlack);
  hS2Plus5X->SetLineColor(kRed);
  hS2Minus5X->SetLineColor(kRed);
  hS1Plus2XY->SetLineColor(kBlue);
  hS1Minus2XY->SetLineColor(kBlue);
  hS2Minus5X->SetLineStyle(kDashed);
  hS1Minus2XY->SetLineStyle(kDashed);
  hS1Only->SetTitle("S1 only");
  hS2Plus5X->SetTitle("+5cm S2 X shift");
  hS2Minus5X->SetTitle("-5cm S2 X shift");
  hS1Plus2XY->SetTitle("+2cm S1 XY shift");
  hS1Minus2XY->SetTitle("-2cm S1 XY shift");

  TGraph *gS1Only     = new TGraph();
  TGraph *gS2Plus5X   = new TGraph();
  TGraph *gS2Minus5X  = new TGraph();
  TGraph *gS1Plus2XY  = new TGraph();
  TGraph *gS1Minus2XY = new TGraph();
  setHistAttr(gS1Only, kBlack);
  setHistAttr(gS2Plus5X, kRed);
  setHistAttr(gS2Minus5X, kRed);
  gS2Minus5X->SetMarkerStyle(21);
  setHistAttr(gS1Plus2XY, kBlue);
  setHistAttr(gS1Minus2XY, kBlue);
  gS1Minus2XY->SetMarkerStyle(21);
  gS1Only->SetTitle("S1 only");
  gS2Plus5X->SetTitle("+5cm S2 X shift");
  gS2Minus5X->SetTitle("-5cm S2 X shift");
  gS1Plus2XY->SetTitle("+2cm S1 XY shift");
  gS1Minus2XY->SetTitle("-2cm S1 XY shift");

  for(int i=0; i<nS4.size(); i++) {
    cout<<"Original S4/S3 "<<nS4.at(i)/nS3.at(i)<<endl;
    cout<<"S1 only S4/S3 "<<nS4_S1Only.at(i)/nS3.at(i)<<endl;
    cout<<"S2 plus 5X S4/S3 "<<nS4_S2Plus5X.at(i)/nS3.at(i)<<endl;
    cout<<"S2 minus 5X S4/S3 "<<nS4_S2Minus5X.at(i)/nS3.at(i)<<endl;
    cout<<"S1 plus 2XY S4/S3 "<<nS4_S1Plus2XY.at(i)/nS3.at(i)<<endl;
    cout<<"S1 minus 2XY S4/S3 "<<nS4_S1Minus2XY.at(i)/nS3.at(i)<<endl;
    cout<<"\n"<<endl;

    hS1Only->Fill(i, abs(nS4[i]-nS4_S1Only[i])/nS4[i]);
    hS2Plus5X->Fill(i, abs(nS4[i]-nS4_S2Plus5X[i])/nS4[i]);
    hS2Minus5X->Fill(i, abs(nS4[i]-nS4_S2Minus5X[i])/nS4[i]);
    hS1Plus2XY->Fill(i, abs(nS4[i]-nS4_S1Plus2XY[i])/nS4[i]);
    hS1Minus2XY->Fill(i, abs(nS4[i]-nS4_S1Minus2XY[i])/nS4[i]);

    gS1Only->SetPoint(gS1Only->GetN(), i, (nS4.at(i)-nS4_S1Only.at(i))/nS4.at(i));
    gS2Plus5X->SetPoint(gS2Plus5X->GetN(), i, (nS4.at(i)-nS4_S2Plus5X.at(i))/nS4.at(i));
    gS2Minus5X->SetPoint(gS2Minus5X->GetN(), i, (nS4.at(i)-nS4_S2Minus5X.at(i))/nS4.at(i));
    gS1Plus2XY->SetPoint(gS1Plus2XY->GetN(), i, (nS4.at(i)-nS4_S1Plus2XY.at(i))/nS4.at(i));
    gS1Minus2XY->SetPoint(gS1Minus2XY->GetN(), i, (nS4.at(i)-nS4_S1Minus2XY.at(i))/nS4.at(i));
  }

  THStack *hs = new THStack("hs", "Fractional change in S4 proton MC count with position changes; Blocks; |Fractional change|");
  hs->Add(hS1Only);
  hs->Add(hS2Plus5X);
  hs->Add(hS2Minus5X);
  hs->Add(hS1Plus2XY);
  hs->Add(hS1Minus2XY);

  TMultiGraph *mg = new TMultiGraph("mg", "Fractional change in S4 proton MC count with position changes; Blocks; Fractional change");
  mg->Add(gS1Only);
  mg->Add(gS2Plus5X);
  mg->Add(gS2Minus5X);
  mg->Add(gS1Plus2XY);
  mg->Add(gS1Minus2XY);

  fout->cd();
  hS1Only->Write();
  hS2Plus5X->Write();
  hS2Minus5X->Write();
  hS1Plus2XY->Write();
  hS1Minus2XY->Write();
  hs->Write();

  gS1Only->Write("gS1Only");
  gS2Plus5X->Write("gS2Plus5X");
  gS2Minus5X->Write("gS2Minus5X");
  gS1Plus2XY->Write("gS1Plus2XY");
  gS1Minus2XY->Write("gS1Minus2XY");
  mg->Write();

  fout->Close();
  delete fout;
} // s1s2MCSysts
