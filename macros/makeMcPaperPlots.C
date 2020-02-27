// makeMcPaperPlots.C
#include "UsefulFunctions.C"

const char* mcDir = "/scratch0/tnonnenm/FixedMC/Prototype-HPTPC-MC/wdir/moderator_data/";

void makeMcPaperPlots(const char* outfile)
{
  THStack *hsMomS4 = new THStack("hsMomS4", "Proton momentum in S4; Proton momentum / (GeV c^{-1}); Events / spill");
  THStack *hsMomTpc = new THStack("hsMomTpc", "Proton momentum in TPC active region; Proton momentum / (GeV c^{-1}); Events / spill");
  THStack *hsKeTpc = new THStack("hsKeTpc", "Proton kinetic energy in TPC active region; Proton kinetic energy / MeV; Events / spill");

  TLegend *legTpc = new TLegend(.15, .6, .5, .85);
  TLegend *legS4  = new TLegend(.15, .6, .5, .85);

  TFile *fout = new TFile(outfile, "recreate");

  for (int b=0; b<5; b++) {
    cout<<b<<" blocks"<<endl;
    TH1D *hMomS4  = new TH1D(Form("hMomS4_%d", b), Form("%d blocks", b), 66, 0.14, 0.8);
    TH1D *hMomTpc = new TH1D(Form("hMomTpc_%d", b), Form("%d blocks", b), 80, 0., 0.8);
    TH1D *hKeTpc  = new TH1D(Form("hKeTpc_%d", b), Form("%d blocks", b), 80, 0., 300.);
    setHistAttr(hMomS4);
    setHistAttr(hMomTpc);
    setHistAttr(hKeTpc);

    TFile *fin = new TFile(Form("%sNewGeom%dBlocks.root", mcDir, b), "read");
    TTree *tree = (TTree*)fin->Get("stepTree");
    UInt_t NSteps;
    float PostStepMom[40000];
    float PostStepKE[40000];
    float PostStepPosX[40000];
    float PostStepPosY[40000];
    float PostStepPosZ[40000];
    int PreStepVolumeID[40000];
    int PostStepVolumeID[40000];
    int PDG[40000];
    tree->SetBranchAddress("NSteps", &NSteps);
    tree->SetBranchAddress("PostStepMom", PostStepMom);
    tree->SetBranchAddress("PostStepKE", PostStepKE);
    tree->SetBranchAddress("PostStepPosX", PostStepPosX);
    tree->SetBranchAddress("PostStepPosY", PostStepPosY);
    tree->SetBranchAddress("PostStepPosZ", PostStepPosZ);
    tree->SetBranchAddress("PreStepVolumeID", PreStepVolumeID);
    tree->SetBranchAddress("PostStepVolumeID", PostStepVolumeID);
    tree->SetBranchAddress("PDG", PDG);

    tree->AddFriend(Form("protonTree%dBlocks", b), "/scratch0/sjones/plots/angularDistS3_newSample/s3ProtonWeightTrees.root");
    int spill;
    double weight;
    tree->SetBranchAddress("spill", &spill);
    tree->SetBranchAddress("weight", &weight);

    tree->GetEntry(tree->GetEntries()-1);
    int nSpills = spill;

    for (int e=0; e<tree->GetEntries(); e++) {
      tree->GetEntry(e);
      for (int n=0; n<NSteps; n++) {
	// Is and S4 proton
	if (PreStepVolumeID[n]!=7 && PostStepVolumeID[n]==7 && 
	    PostStepMom[n]>140 && PDG[n]==2212 && PostStepPosY[n]<360) {
	  hMomS4->Fill(PostStepMom[n]*0.001, weight);
	}
	if (PreStepVolumeID[n]!=6 && PostStepVolumeID[n]==6 && PostStepMom[n]>5 && PDG[n]==2212) {
	  hKeTpc->Fill(PostStepKE[n], weight);
	  hMomTpc->Fill(PostStepMom[n]*0.001, weight);
	}
      }
    }
    fin->Close();
    delete fin;

    hMomS4->Scale(1./nSpills);
    hMomTpc->Scale(1./nSpills);
    hKeTpc->Scale(1./nSpills);

    double es4  = 0.;
    double etpc = 0.;

    if (b==0) {
      hMomS4->SetLineColor(getColourFromBlock(0));
      hMomTpc->SetLineColor(getColourFromBlock(0));
      hKeTpc->SetLineColor(getColourFromBlock(0));
      double is4  = hMomS4->IntegralAndError(1, hMomS4->GetNbinsX(), es4);
      double itpc = hMomTpc->IntegralAndError(1, hMomTpc->GetNbinsX(), etpc);
      legS4->AddEntry(hMomS4, Form("0 blocks - %.3g #pm %.2g / spill", is4, es4), "l");
      legTpc->AddEntry(hMomTpc, Form("0 blocks - %.3g #pm %.1g / spill", itpc, etpc), "l");
    }
    else if (b==1) {
      hMomS4->SetLineColor(getColourFromBlock(1));
      hMomTpc->SetLineColor(getColourFromBlock(1));
      hKeTpc->SetLineColor(getColourFromBlock(1));
      double is4  = hMomS4->IntegralAndError(1, hMomS4->GetNbinsX(), es4);
      double itpc = hMomTpc->IntegralAndError(1, hMomTpc->GetNbinsX(), etpc);
      legS4->AddEntry(hMomS4, Form("1 block - %.4g #pm %.2g / spill", is4, es4), "l");
      legTpc->AddEntry(hMomTpc, Form("1 block - %.3g #pm %.1g / spill", itpc, etpc), "l");
    }
    else if (b==2) {
      hMomS4->SetLineColor(getColourFromBlock(2));
      hMomTpc->SetLineColor(getColourFromBlock(2));
      hKeTpc->SetLineColor(getColourFromBlock(2));
      double is4  = hMomS4->IntegralAndError(1, hMomS4->GetNbinsX(), es4);
      double itpc = hMomTpc->IntegralAndError(1, hMomTpc->GetNbinsX(), etpc);
      legS4->AddEntry(hMomS4, Form("2 blocks - %.3g #pm %.2g / spill", is4, es4), "l");
      legTpc->AddEntry(hMomTpc, Form("2 blocks - %.3g #pm %.2g / spill", itpc, etpc), "l");
    }
    else if (b==3) {
      hMomS4->SetLineColor(getColourFromBlock(3));
      hMomTpc->SetLineColor(getColourFromBlock(3));
      hKeTpc->SetLineColor(getColourFromBlock(3));
      double is4  = hMomS4->IntegralAndError(1, hMomS4->GetNbinsX(), es4);
      double itpc = hMomTpc->IntegralAndError(1, hMomTpc->GetNbinsX(), etpc);
      legS4->AddEntry(hMomS4, Form("3 blocks - %.3g #pm %.2g / spill", is4, es4), "l");
      legTpc->AddEntry(hMomTpc, Form("3 blocks - %.3g #pm %.2g / spill", itpc, etpc), "l");
    }
    else if (b==4) {
      hMomS4->SetLineColor(getColourFromBlock(4));
      hMomTpc->SetLineColor(getColourFromBlock(4));
      hKeTpc->SetLineColor(getColourFromBlock(4));
      double is4  = hMomS4->IntegralAndError(1, hMomS4->GetNbinsX(), es4);
      double itpc = hMomTpc->IntegralAndError(1, hMomTpc->GetNbinsX(), etpc);
      legS4->AddEntry(hMomS4, Form("4 blocks - %.3g #pm %.1g / spill", is4, es4), "l");
      legTpc->AddEntry(hMomTpc, Form("4 blocks - %.2g #pm %.1g / spill", itpc, etpc), "l");
    }

    hsMomS4->Add(hMomS4);
    hsMomTpc->Add(hMomTpc);
    hsKeTpc->Add(hKeTpc);

    fout->cd();
    hMomS4->Write();
    hMomTpc->Write();
    hKeTpc->Write();
  } // Loop over blocks

  fout->cd();
  hsMomS4->Write();
  hsMomTpc->Write();
  hsKeTpc->Write();
  legTpc->Write("legTpc");
  legS4->Write("legS4");

  fout->Close();
} // makeMcPaperPlots
