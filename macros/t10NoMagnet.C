// t10NoMagnet.C
#include "UsefulFunctions.C"

const TVector3 vS1 = (D1_ULB+D1_ULT+D1_URB+D1_URT)*0.25;

void t10NoMagnet(const char* outfile)
{
  const int hitMax = 50;
  // Edges of S3 in beam coordinate system
  const double s3StartX = -0.601158 + 0.026;
  const double s3EndX   = 1.08086 + 0.026;
  const double s3s1StartY = 9.07239 + 1.767;
  const double s3s1EndY   = 8.89731 + 1.767;
  // Measured from the wire chamber
  const double s3StartWcX = -0.601158 + 0.052;
  const double s3EndWcX   = 1.08086 + 0.052;
  const double s3StartWcY = 9.07239 + 2.054;
  const double s3EndWcY   = 8.89731 + 2.054;
  // z coordinates
  const double s3BarTop    = 62.; 
  const double s3BarBottom = -60.;

  // This is before the shift is applied
  const double proLow  = 53;
  const double proHi  = 80.;
  const double piLow = 35.75;
  const double piHi  = 37.75;

  TH1D *hUtof = new TH1D("hUtof", "Time of flight; t_{S3} - t_{S1} / ns; Events", 250, 25, 125);
  setHistAttr(hUtof);
  hUtof->SetLineColor(kBlack);
  TH2D *hHitsXY = new TH2D("h2dHitsXY", "Distribution of hits within S3; x / m; y / m; Hits", 100, -1., 0.6, 22, -0.6, 0.61);
  setHistAttr(hHitsXY);
  TH1D *hMom = new TH1D("hMom", "Proton momentum measured in S3; Proton momentum [GeV/c]; Events", 120, 0.3, 0.9);
  setHistAttr(hMom);
  hMom->SetLineColor(kBlack);
  TH1D *hKE = new TH1D("hKE", "Proton kinetic energy measured in S3; Proton kinetic energy / MeV; Events", 100, 40., 350.);
  setHistAttr(hKE);
  hKE->SetLineColor(kBlack);
  TH1D *hThetaS1pro = new TH1D("hThetaS1pro", "Angular distribution of proton hits in S3 (S1 trigger only); #theta / degrees; Events / spill / degree", binnum, binsTheta);
  setHistAttr(hThetaS1pro);
  TH1D *hPhiS1pro = new TH1D("hPhiS1pro", "Angular distribution of proton hits in S3 (S1 trigger only); #phi / degrees; Events / spill / degree",  22, -3.22, 3.35);
  setHistAttr(hPhiS1pro);
  TH1D *hThetaS1pi = new TH1D("hThetaS1pi", "Angular distribution of pion hits in S3 (S1 trigger only); #theta / degrees; Events / spill / degree", binnum, binsTheta);
  setHistAttr(hThetaS1pi);
  TH1D *hPhiS1pi = new TH1D("hPhiS1pi", "Angular distribution of pion hits in S3 (S1 trigger only); #phi / degrees; Events / spill / degree", 22, -3.22, 3.35);
  setHistAttr(hPhiS1pi);
  TH1D *hThetaS1Ratio = new TH1D("hThetaS1Ratio", "S1 #cap S3 angular distribution of proton/MIP ratio; #theta / degrees; Protons/MIPs", binnum, binsTheta);
  setHistAttr(hThetaS1Ratio);
  TH1D *hPhiS1Ratio = new TH1D("hPhiS1ratio", "S1 #cap S3 angular distribution of proton/MIP; #phi / degrees; Protons/MIPs", 22, -3.22, 3.35);
  setHistAttr(hPhiS1Ratio);

  // File with no magnet bend. Beam on nominal axis
  TChain *tree = new TChain("tree");
  // tree->Add(Form("%s/Data_2018_8_17_b8_800MeV_0block.root", ustofDir));
  // tree->Add(Form("%s/Data_2018_8_17_b9_800MeV_0block.root", ustofDir));
  // tree->Add(Form("%s/Data_2018_8_17_b10_800MeV_0block.root", ustofDir));
  tree->Add(Form("%s/Data_2018_9_9_b2_800MeV_0block_nobend.root", ustofDir));
  
  int nhit;
  float xToF[hitMax];
  float yToF[hitMax];
  float A1ToF[hitMax];
  float A2ToF[hitMax];
  int nBar[hitMax];
  double tToF[hitMax];
  double tTrig;
  double tS1;
  double tSoSd;

  tree->SetBranchAddress("nhit", &nhit);
  tree->SetBranchAddress("nBar", nBar);
  tree->SetBranchAddress("xToF", xToF);
  tree->SetBranchAddress("yToF", yToF);
  tree->SetBranchAddress("A1ToF", A1ToF);
  tree->SetBranchAddress("A2ToF", A2ToF);
  tree->SetBranchAddress("tS1", &tS1);
  tree->SetBranchAddress("tToF", tToF);
  tree->SetBranchAddress("tTrig", &tTrig);
  tree->SetBranchAddress("tSoSd", &tSoSd);

  for (int t=0; t<tree->GetEntries(); t++) {
    tree->GetEntry(t);
    for (int nh=0; nh<nhit; nh++) {
      double tofCalc = tToF[nh] - tS1;
      // Calculate x, y z positions relative to S1
      TVector3 utofCoords = GetUtofGlobalCoords(xToF[nh], yToF[nh]);
      double positionX = (xToF[nh]/168)*(s3EndX - s3StartX) + s3StartX;
      double positionY = (xToF[nh]/168.)*(s3s1EndY - s3s1StartY) + s3s1StartY;
      double positionZ = (yToF[nh] + s3BarBottom + 2.75) / 100.;
      double angleTheta = getThetaFromGlobal(utofCoords);
      double anglePhi   = getPhiFromGlobal(utofCoords);
      hUtof->Fill(tofCalc);
      hHitsXY->Fill(positionX, positionZ);
      if ( tofCalc > piLow && tofCalc < piHi ) {
	hThetaS1pi->Fill(angleTheta);
	hPhiS1pi->Fill(anglePhi);
      }
      else if (tofCalc > proLow && tofCalc < proHi && 
	       A1ToF[nh] > A1CutVec[nBar[nh]] && A2ToF[nh] > A2CutVec[nBar[nh]]) {
	hThetaS1pro->Fill(angleTheta);
	hPhiS1pro->Fill(anglePhi);
	hMom->Fill(momFromTime(0.938, vS1, utofCoords, tofCalc)/1000.);
	hKE->Fill(keFromTime(0.938, vS1, utofCoords, tofCalc));
	// hMom->Fill(momFromTime(0.938, 10.9, tofCalc)/1000.);
	// hKE->Fill(keFromTime(0.938, 10.9, tofCalc));
      }
    } // Loop over hits
  } // Loop over entries

  hThetaS1Ratio->Divide(hThetaS1pro, hThetaS1pi, 1., 1., "B");
  hPhiS1Ratio->Divide(hPhiS1pro, hPhiS1pi, 1., 1., "B");

  TFile *fout = new TFile(outfile, "recreate");
  fout->cd();
  hUtof->Write();
  hMom->Write();
  hKE->Write();
  hHitsXY->Write();
  hThetaS1Ratio->Write();
  hPhiS1Ratio  ->Write();
  hThetaS1pro->Write();
  hPhiS1pro  ->Write();
  hThetaS1pi ->Write();
  hPhiS1pi   ->Write();

  fout->Close();
  delete fout;
}
