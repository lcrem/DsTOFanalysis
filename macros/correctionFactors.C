// correctionFactors.C

void correctionFactors(const char* saveDir,
		       const char* ustofDir="/zfs_home/sjones/mylinktoutof")
{
  gROOT->SetBatch(kTRUE);
  // Edges of S3 in beam coordinate system
  const double s3StartX = -0.5168;
  const double s3EndX   = 0.9970;
  const double s3s1StartY = 9.0569 + 1.77;
  const double s3s1EndY   = 8.9146 + 1.77;
  // z coordinates
  const double s3BarTop    = 62.; 
  const double s3BarBottom = -60.;
  // Define the runs to be used for varying number of blocks
  const char* str0Block = "Data_2018_8_31_b2_800MeV_0block.root";
  const char* str1Block = "Data_2018_9_1_b4_800MeV_1block_bend4cm.root";
  const char* str2Block = "Data_2018_9_1_b2_800MeV_2block_bend4cm.root";
  const char* str3Block = "Data_2018_9_1_b3_800MeV_3block_bend4cm.root";
  const char* str4Block = "Data_2018_9_1_b8_800MeV_4block_bend4cm.root";

  // Integrate this hist between the angular limits of S2 then scale by this number
  // S2 angular extremities
  const double s2ThetaLow = -0.359;
  const double s2ThetaHi  = 3.957;
  const double s2PhiLow = -2.010;
  const double s2PhiHi  = 2.942;
  // S3 angular extremities
  const double s3ThetaLow = -3.3587;
  const double s3ThetaHi  =  6.2416;
  const double s3PhiLow = -3.20449;
  const double s3PhiHi  =  3.33388;
  // S4 angular extremities
  const double s4ThetaLow = 0.401;
  const double s4ThetaHi  = 6.083;
  const double s4PhiLow = -1.427;
  const double s4PhiHi  = 1.771;

  TFile *fout = new TFile(Form("%s/correctionFactorsPlots.root", saveDir), "recreate");

  for (int nBlocks = 0; nBlocks <= 4; nBlocks++) {

    TH2D *hAngleS1S2 = new TH2D(Form("hAngleS1S2%d", nBlocks), Form("Angular distribution of hits in S3 (S1 & S2 triggers), %d blocks; #theta / degrees; #phi / degrees; Events / spill", nBlocks), 100, -3.8, 6.2, 22, -3.15, 3.25);
    TH2D *hAngleS1   = new TH2D(Form("hAngleS1%d", nBlocks), Form("Angular distribution of hits in S3 (S1 trigger only), %d blocks; #theta / degrees; #phi / degrees; Events / spill", nBlocks), 100, -3.8, 6.2, 22, -3.15, 3.25);
    TH2D *hDetS1S2 = new TH2D(Form("hDetS1S2%d", nBlocks), Form("Distribution of hits in S3 (S1 & S2 triggers), %d blocks; x / cm; Bar; Events / spill", nBlocks), 100, -10., 170., 22, 2., 120.);
    TH2D *hDetS1   = new TH2D(Form("hDetS1%d", nBlocks), Form("Distribution of hits in S3 (S1 trigger only), %d blocks; x / cm ; Bar; Events / spill", nBlocks), 100, -10., 170., 22, 2., 120.);

    const char* nustof;
    if (nBlocks == 0) nustof = Form("%s/%s", ustofDir, str0Block);
    else if (nBlocks==1) nustof = Form("%s/%s", ustofDir, str1Block);
    else if (nBlocks==2) nustof = Form("%s/%s", ustofDir, str2Block);
    else if (nBlocks==3) nustof = Form("%s/%s", ustofDir, str3Block);
    else if (nBlocks==4) nustof = Form("%s/%s", ustofDir, str4Block);

    // Read in ustof file
    TFile *finustof = new TFile(nustof, "read");
    TTree *utree = (TTree*)finustof->Get("tree");
    double tTrig;
    double tS1;
    double tSoSd;
    float xToF[50];
    float yToF[50];
    int nhit;
    int nBar[50];

    utree->SetBranchAddress("tTrig", &tTrig);
    utree->SetBranchAddress("xToF", xToF);
    utree->SetBranchAddress("yToF", yToF);
    utree->SetBranchAddress("nhit", &nhit);
    utree->SetBranchAddress("nBar", nBar);

    for (int t=0; t<utree->GetEntries(); t++) {
      utree->GetEntry(t);
      // Has an S1 and S2 hit
      if (tTrig != 0 ) {
	for (int n=0; n<nhit; n++) {
	  double positionX = ((xToF[n] - 4.) / 152.)*(s3EndX - s3StartX) + s3StartX;
	  double positionY = ((xToF[n] - 4.) / 152.)*(s3s1EndY - s3s1StartY)+s3s1StartY;
	  double angleTheta = TMath::ATan2(positionX, positionY) * (180. / TMath::Pi());
	  double positionZ  = (yToF[n] + s3BarBottom + 2.75) / 100.;
	  double anglePhi   = TMath::ATan2(positionZ, positionY) * (180./TMath::Pi());
	  hAngleS1S2->Fill(angleTheta, anglePhi);
	  hDetS1S2->Fill(xToF[n], yToF[n]);
	} // for (int n=0; n<nhit; n++) 
      } // if (tTrig !=0 ) 
      else {
	for (int n=0; n<nhit; n++) {
	  double positionX = ((xToF[n] - 4.) / 152.)*(s3EndX - s3StartX) + s3StartX;
	  double positionY = ((xToF[n] - 4.) / 152.)*(s3s1EndY - s3s1StartY)+s3s1StartY;
	  double angleTheta = TMath::ATan2(positionX, positionY) * (180. / TMath::Pi());
	  double positionZ  = (yToF[n] + s3BarBottom + 2.75) / 100.;
	  double anglePhi   = TMath::ATan2(positionZ, positionY) * (180./TMath::Pi());
	  hAngleS1->Fill(angleTheta, anglePhi);
	  hDetS1->Fill(xToF[n], yToF[n]);
	} // for (int n=0; n<nhit; n++) 
      } // S1 trigger only
    } // for (int t=0; t<utree->GetEntries(); t++) 
    double s3Int = hAngleS1->Integral();
    hAngleS1->Scale(1. / s3Int);
    cout<<nBlocks<<" blocks"<<endl;
    cout<<"S2 integral = "<<hAngleS1->Integral(hAngleS1->GetXaxis()->FindBin(s2ThetaLow), hAngleS1->GetXaxis()->FindBin(s2ThetaHi), hAngleS1->GetYaxis()->FindBin(s2PhiLow),  hAngleS1->GetYaxis()->FindBin(s2PhiHi))<<endl;
    
    cout<<"S2 S4 overlap integral = "<<hAngleS1->Integral(hAngleS1->GetXaxis()->FindBin(s4ThetaLow), hAngleS1->GetXaxis()->FindBin(s2ThetaHi), hAngleS1->GetYaxis()->FindBin(s4PhiLow),  hAngleS1->GetYaxis()->FindBin(s4PhiHi))<<endl;
    
    cout<<"S4 integral = "<<hAngleS1->Integral(hAngleS1->GetXaxis()->FindBin(s4ThetaLow), hAngleS1->GetXaxis()->FindBin(s4ThetaHi), hAngleS1->GetYaxis()->FindBin(s4PhiLow),  hAngleS1->GetYaxis()->FindBin(s4PhiHi))<<endl;
    
    fout->cd();
    hAngleS1->Write();
    hAngleS1S2->Write();
    hDetS1->Write();
    hDetS1S2->Write();


  } // for (int nBlocks = 0; nBlocks <= 4; nBlocks++)

  fout->Close();
}
