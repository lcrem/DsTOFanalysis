// angularDistS3.C

// Outputs momentum in GeV/c
double momFromTime(const double mass, const double baseline, const double time)
{
  double mom = 0.;
  mom = (mass / 3e8) * baseline * ( 1. / TMath::Sqrt( pow(time*1e-9, 2) - (pow((baseline / 9e16), 2)) ) );
  return mom;
}

void angularDistS3(const char* saveDir,
		   const char* ustofDir="/zfs_home/sjones/mylinktoutof/") 
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

  // Time of flight cuts for S1 to S3
  // Particles travelling at c should cross the distance in between about 36.9ns
  // and 35.6ns
  const double tLight = 36.4; // Expected time in ns of light speed particles according to AK
  // This is before the shift is applied
  const double proLow  = -12.5;
  const double proHi   = 40.5;
  const double piLow = -29.06;
  const double piHi  = -28.06;
  // S3 amplitude cut for protons
  // Apply to A1ToF and A2ToF
  const double ACut = 0.25;

  THStack *hsThetaS1pro   = new THStack("hsThetaS1pro", "Angular distribution of S3 proton hits with S1 trigger only; #theta / degrees; Events / spill");
  THStack *hsThetaS1S2pro = new THStack("hsThetaS1S2pro", "Angular distribution of S3 proton hits with S1 & S2; #theta / degrees; Events / spill");
  THStack *hsPhiS1pro     = new THStack("hsPhiS1pro", "Angular distribution of S3 proton hits with S1 trigger only; #phi / degrees; Events / spill");
  THStack *hsPhiS1S2pro   = new THStack("hsPhiS1S2pro", "Angular distribution of S3 proton hits with S1 & S2; #phi / degrees; Events / spill");

  THStack *hsThetaS1pi   = new THStack("hsThetaS1pi", "Angular distribution of S3 pion hits with S1 trigger only; #theta / degrees; Events / spill");
  THStack *hsThetaS1S2pi = new THStack("hsThetaS1S2pi", "Angular distribution of S3 pion hits with S1 & S2; #theta / degrees; Events / spill");
  THStack *hsPhiS1pi     = new THStack("hsPhiS1pi", "Angular distribution of S3 pion hits with S1 trigger only; #phi / degrees; Events / spill");
  THStack *hsPhiS1S2pi   = new THStack("hsPhiS1S2pi", "Angular distribution of S3 pion hits with S1 & S2; #phi / degrees; Events / spill");

  THStack *hsPhiS1S2ratio   = new THStack("hsPhiS1S2ratio", "Angular distribution of S3 proton/MIP ratio (S1 & S2 trigger); #phi / degrees;  Protons/MIPs");
  THStack *hsThetaS1S2ratio = new THStack("hsThetaS1S2ratio", "Angular distribution of S3 proton/MIP ratio (S1 & S2 trigger); #theta / degrees;  Protons/MIPs");
  THStack *hsPhiS1ratio   = new THStack("hsPhiS1ratio", "Angular distribution of S3 proton/MIP ratio (S1 trigger only); #phi / degrees;  Protons/MIPs");
  THStack *hsThetaS1ratio = new THStack("hsThetaS1ratio", "Angular distribution of S3 proton/MIP ratio (S1 trigger only); #theta / degrees;  Protons/MIPs");

  THStack *hsMomS1S2 = new THStack("hsMomS1S2", "Proton momentum measured in S3 (S1 & S2 trigger); Proton momentum [GeV/c]; Events / spill");
  THStack *hsMomS1 = new THStack("hsMomS1", "Proton momentum measured in S3 (S1 trigger only); Proton momentum [GeV/c]; Events / spill");

  THStack *hsutof1dS1S2 = new THStack("hsutof1dS1S2", "Time of flight as measured in S3 (S1 & S2 trigger); Time of flight / ns; Events / spill");
  THStack *hsutof1dS1   = new THStack("hsutof1dS1", "Time of flight as measured in S3 (S1 trigger only); Time of flight / ns; Events / spill");

  TFile *fout = new TFile(Form("%s/angularDistS3plots.root", saveDir), "recreate");

  TLegend *leg = new TLegend(0.15, 0.55, 0.3, 0.85);
  TLegend *legTof = new TLegend(0.71, 0.53, 0.88, 0.85);

  TH2D *hMom2D_0blkQ = new TH2D("hMom2D_0blkQ", "Quick peak 0 block; x / cm; Bar", 100, -10, 170, 22, 0.5, 22.5);
  TH2D *hMom2D_0blkS = new TH2D("hMom2D_0blkS", "Slow peak 0 block; x / cm; Bar", 100, -10, 170, 22, 0.5, 22.5);
  TH2D *hMom2D_1blkQ = new TH2D("hMom2D_1blkQ", "Quick peak 1 block; x / cm; Bar", 100, -10, 170, 22, 0.5, 22.5);
  TH2D *hMom2D_1blkS = new TH2D("hMom2D_1blkS", "Slow peak 1 block; x / cm; Bar", 100, -10, 170, 22, 0.5, 22.5);
  TH2D *hMom2D_2blkQ = new TH2D("hMom2D_2blkQ", "Quick peak 2 block; x / cm; Bar", 100, -10, 170, 22, 0.5, 22.5);
  TH2D *hMom2D_2blkS = new TH2D("hMom2D_2blkS", "Slow peak 2 block; x / cm; Bar", 100, -10, 170, 22, 0.5, 22.5);

  for (int nBlocks = 0; nBlocks <= 4; nBlocks++) {
    int nSpills = 0;
    TH1D *hThetaS1pro   = new TH1D(Form("hThetaS1pro%d", nBlocks), Form("Angular distribution of proton hits in S3 (S1 trigger only), %d blocks; #theta / degrees; Events / spill", nBlocks), 100, -3.8, 6.2);
    hThetaS1pro->Sumw2();
    TH1D *hThetaS1S2pro = new TH1D(Form("hThetaS1S2pro%d", nBlocks), Form("Angular distribution of proton hits in S3 (S1 & S2 triggers), %d blocks; #theta / degrees; Events / spill", nBlocks), 100, -3.8, 6.2);
    hThetaS1S2pro->Sumw2();
    TH1D *hPhiS1pro   = new TH1D(Form("hPhiS1pro%d", nBlocks), Form("Angular distribution of proton hits in S3 (S1 trigger only), %d blocks; #phi / degrees; Events / spill", nBlocks), 22, -3.15, 3.25);
    hPhiS1pro->Sumw2();
    TH1D *hPhiS1S2pro = new TH1D(Form("hPhiS1S2pro%d", nBlocks), Form("Angular distribution of proton hits in S3 (S1 & S2 triggers), %d blocks; #phi / degrees; Events / spill", nBlocks), 22, -3.15, 3.25);
    hPhiS1S2pro->Sumw2();

    TH1D *hThetaS1pi   = new TH1D(Form("hThetaS1pi%d", nBlocks), Form("Angular distribution of pion hits in S3 (S1 trigger only), %d blocks; #theta / degrees; Events / spill", nBlocks), 100, -3.8, 6.2);
    hThetaS1pi->Sumw2();
    TH1D *hThetaS1S2pi = new TH1D(Form("hThetaS1S2pi%d", nBlocks), Form("Angular distribution of pion hits in S3 (S1 & S2 triggers), %d blocks; #theta / degrees; Events / spill", nBlocks), 100, -3.8, 6.2);
    hThetaS1S2pi->Sumw2();
    TH1D *hPhiS1pi   = new TH1D(Form("hPhiS1pi%d", nBlocks), Form("Angular distribution of pion hits in S3 (S1 trigger only), %d blocks; #phi / degrees; Events / spill", nBlocks), 22, -3.15, 3.25);
    hPhiS1pi->Sumw2();
    TH1D *hPhiS1S2pi = new TH1D(Form("hPhiS1S2pi%d", nBlocks), Form("Angular distribution of pion hits in S3 (S1 & S2 triggers), %d blocks; #phi / degrees; Events / spill", nBlocks), 22, -3.15, 3.25);
    hPhiS1S2pi->Sumw2();

    TH1D *hThetaS1ratio = new TH1D(Form("hThetaS1ratio%d", nBlocks), Form("Angular distribution of proton/MIP ratio in S3 (S1 trigger only), %d blocks; #phi / degrees; Protons/MIPs", nBlocks), 100, -3.8, 6.2);
    TH1D *hPhiS1ratio   = new TH1D(Form("hPhiS1ratio%d", nBlocks), Form("Angular distribution of proton/MIP in S3 (S1 trigger only), %d blocks; #phi / degrees; Protons/MIPs", nBlocks), 22, -3.15, 3.25);
    TH1D *hThetaS1S2ratio = new TH1D(Form("hThetaS1S2ratio%d", nBlocks), Form("Angular distribution of proton/MIP ratio in S3 (S1 & S2 triggers), %d blocks; #phi / degrees; Protons/MIPs", nBlocks), 100, -3.8, 6.2);
    TH1D *hPhiS1S2ratio   = new TH1D(Form("hPhiS1S2ratio%d", nBlocks), Form("Angular distribution of proton/MIP in S3 (S1 & S2 triggers), %d blocks; #phi / degrees; Protons/MIPs", nBlocks), 22, -3.15, 3.25);

    TH1D *hutof1dS1 = new TH1D(Form("hutof1dS1_%d",nBlocks), Form("Time of flight, %d blocks (S1 trigger only); S3 - S1 / ns; Events / spill", nBlocks), 250, 25, 120);
    hutof1dS1->Sumw2();
    TH1D *hutof1dS1S2 = new TH1D(Form("hutof1dS1S2_%d",nBlocks), Form("Time of flight, %d blocks (S1 & S2 trigger); S3 - S1 / ns; Events / spill", nBlocks), 250, 25, 120);
    hutof1dS1S2->Sumw2();

    TH1D *hMomS1S2 = new TH1D(Form("hMomS1S2_%d",nBlocks), Form("Proton momentum measured in S3, %d blocks; Proton momentum [GeV/c]; Events / spill", nBlocks), 100, 0.3, 0.7);
    hMomS1S2->Sumw2();
    TH1D *hMomS1 = new TH1D(Form("hMomS1_%d",nBlocks), Form("Proton momentum measured in S3, %d blocks; Proton momentum [GeV/c]; Events / spill", nBlocks), 100, 0.3, 0.7);
    hMomS1->Sumw2();

    //    TH2D *hMomZS1 = new TH2D(Form("hMomZS1Q_%d",nBlocks), Form("Proton momentum measured in S3 compared to hit position (S1 trigger only), %d blocks; Proton momentum [GeV/c]; Bar in S3",nBlocks), 100, 0.3, 0.7, 22, 0.5, 22.5);
    //    TH2D *hMomYS1 = new TH2D(Form("hMomYS1Q_%d",nBlocks), Form("Proton momentum measured in S3 compared to hit position (S1 trigger only), %d blocks; Proton momentum [GeV/c]; x / cm",nBlocks), 100, 0.3, 0.7, 100, -10., 170.);
    TH2D *hMomZS1 = new TH2D(Form("hMomZS1_%d",nBlocks), Form("Proton momentum measured in S3 compared to hit position (S1 trigger only), %d blocks; Proton momentum [GeV/c]; Bar in S3",nBlocks), 100, 0.3, 0.7, 22, 0.5, 22.5);
    TH2D *hMomYS1 = new TH2D(Form("hMomYS1_%d",nBlocks), Form("Proton momentum measured in S3 compared to hit position (S1 trigger only), %d blocks; Proton momentum [GeV/c]; x / cm",nBlocks), 100, 0.3, 0.7, 100, -10., 170.);
    TH2D *hMomZS12 = new TH2D(Form("hMomZS12_%d",nBlocks), Form("Proton momentum measured in S3 compared to hit position (S1 & S2 trigger), %d blocks; Proton momentum [GeV/c]; Bar in S3",nBlocks), 100, 0.3, 0.7, 22, 0.5, 22.5);
    TH2D *hMomYS12 = new TH2D(Form("hMomYS12_%d",nBlocks), Form("Proton momentum measured in S3 compared to hit position (S1 & S2 trigger), %d blocks; Proton momentum [GeV/c]; x / cm",nBlocks), 100, 0.3, 0.7, 100, -10., 170.);

    // Number of protons and number of MIPs
    int nP  = 0;
    int nPi = 0;
    // Define signal and background functions to be fitted
    TF1 *sPro = new TF1("sPro", "gaus", (tLight - (piLow+piHi)/2.) + proLow, (tLight - (piLow+piHi)/2.) + proHi);
    TF1 *sPi  = new TF1("sPi", "gaus", (tLight - (piLow+piHi)/2.) + piLow, (tLight - (piLow+piHi)/2.) + piHi);

    const char* nustof;
    if (nBlocks==0) nustof = Form("%sData_2018_8_31_b2_800MeV_0block.root", ustofDir);
    else if (nBlocks==1) nustof = Form("%sData_2018_9_1_b4_800MeV_1block_bend4cm.root", ustofDir);
    else if (nBlocks==2) nustof = Form("%sData_2018_9_1_b2_800MeV_2block_bend4cm.root", ustofDir);
    else if (nBlocks==3) nustof = Form("%sData_2018_9_1_b3_800MeV_3block_bend4cm.root", ustofDir);
    else if (nBlocks==4) nustof = Form("%sData_2018_9_1_b8_800MeV_4block_bend4cm.root", ustofDir);

    TFile *futof = new TFile(nustof, "read");

    double tToF[50];
    float xToF[50];
    float yToF[50];
    float A1ToF[50];
    float A2ToF[50];
    double tTrig;
    double tS1;
    double tSoSd;
    int nhit;
    int nBar[50];

    TTree *tree = (TTree*)futof->Get("tree");

    tree->SetBranchAddress("xToF", xToF);
    tree->SetBranchAddress("yToF", yToF);
    tree->SetBranchAddress("A1ToF", A1ToF);
    tree->SetBranchAddress("A2ToF", A2ToF);
    tree->SetBranchAddress("nhit", &nhit);
    tree->SetBranchAddress("tS1", &tS1);
    tree->SetBranchAddress("tToF", tToF);
    tree->SetBranchAddress("tTrig", &tTrig);
    tree->SetBranchAddress("tSoSd", &tSoSd);
    tree->SetBranchAddress("nBar", nBar);

    double lastSpill = 0.; 

    for (int t=0; t<tree->GetEntries(); t++) {
      tree->GetEntry(t);

      // Count number of spills
      if (tSoSd != lastSpill && tSoSd -lastSpill > 1e9) {
	nSpills++;
	lastSpill = tSoSd;
      } // if (tSoSd != lastSpill && tSoSd -lastSpill > 1e9) 

      // Only select those events with multiplicity of 1 - advice from AK
      if (nhit == 1) {
	double tofCalc = tToF[0] - tS1 + (tLight - (piLow + piHi)/2.);
	// Calculate x, y z positions relative to S1
	double positionX = ((xToF[0] - 4.) / 152.)*(s3EndX - s3StartX) + s3StartX;;
	double positionY = ((xToF[0] - 4.) / 152.)*(s3s1EndY - s3s1StartY) + s3s1StartY;
	double positionZ = (yToF[0] + s3BarBottom + 2.75) / 100.;
	double angleTheta = TMath::ATan(positionX / positionY) * (180./TMath::Pi());
	double anglePhi   = TMath::ATan(positionZ / positionY) * (180./TMath::Pi());
	// S1 trigger only
	if (tTrig == 0) {
	  hutof1dS1->Fill(tofCalc);
	  // Separate protons and MIPs using timing and amplitude cuts
	  // Is a MIP
	  if ( tofCalc > (tLight - (piLow+piHi)/2.) + piLow && tofCalc < (tLight - (piLow+piHi)/2.) + piHi ) {
	    nPi++;
	    hThetaS1pi->Fill(angleTheta);
	    hPhiS1pi->Fill(anglePhi);
	  } // if ( tofCalc > (tLight - (piLow+piHi)/2.) + piLow && tofCalc < (tLight - (piLow+piHi)/2.) + piHi )
	  // Is a proton
	  else if ( tofCalc > (tLight - (piLow+piHi)/2.) + proLow && tofCalc < (tLight - (piLow+piHi)/2.) + proHi && A1ToF[0] > ACut && A2ToF[0] > ACut) {
	    nP++;
	    hThetaS1pro->Fill(angleTheta);
	    hPhiS1pro->Fill(anglePhi);
	    hMomS1->Fill(momFromTime(0.938, 10.9, tofCalc));
	    hMomZS1->Fill(momFromTime(0.938, 10.9, tofCalc), nBar[0]);
	    hMomYS1->Fill(momFromTime(0.938, 10.9, tofCalc), xToF[0]);
	    if (nBlocks == 0 && momFromTime(0.938, 10.9, tofCalc) > 0.595) hMom2D_0blkQ->Fill(xToF[0], nBar[0]);
	    else if (nBlocks == 0 && momFromTime(0.938, 10.9, tofCalc) < 0.595) hMom2D_0blkS->Fill(xToF[0], nBar[0]);
	    else if (nBlocks == 1 && momFromTime(0.938, 10.9, tofCalc) > 0.570) hMom2D_1blkQ->Fill(xToF[0], nBar[0]);
	    else if (nBlocks == 1 && momFromTime(0.938, 10.9, tofCalc) < 0.570) hMom2D_1blkS->Fill(xToF[0], nBar[0]);
	    else if (nBlocks == 2 && momFromTime(0.938, 10.9, tofCalc) > 0.525) hMom2D_2blkQ->Fill(xToF[0], nBar[0]);
	    else if (nBlocks == 2 && momFromTime(0.938, 10.9, tofCalc) < 0.525) hMom2D_2blkS->Fill(xToF[0], nBar[0]);
	    
	  } // else if ( tofCalc > (tLight - (piLow+piHi)/2.) + proLow && tofCalc < (tLight - (piLow+piHi)/2.) + proHi )
	}
	// S1 & S2 trigger
	else {
	  hutof1dS1S2->Fill(tofCalc);
	  // Separate protons and MIPs using timing and amplitude cuts
	  // Is a MIP
	  if ( tofCalc > (tLight - (piLow+piHi)/2.) + piLow && tofCalc < (tLight - (piLow+piHi)/2.) + piHi ) {
	    nPi++;
	    hThetaS1S2pi->Fill(angleTheta);
	    hPhiS1S2pi->Fill(anglePhi);
	  } // if ( tofCalc > (tLight - (piLow+piHi)/2.) + piLow && tofCalc < (tLight - (piLow+piHi)/2.) + piHi )
	  // Is a proton
	  else if ( tofCalc > (tLight - (piLow+piHi)/2.) + proLow && tofCalc < (tLight - (piLow+piHi)/2.) + proHi && A1ToF[0] > ACut && A2ToF[0] > ACut) {
	    nP++;
	    hThetaS1S2pro->Fill(angleTheta);
	    hPhiS1S2pro->Fill(anglePhi);
	    hMomS1S2->Fill(momFromTime(0.938, 10.9, tofCalc));
	    hMomZS12->Fill(momFromTime(0.938, 10.9, tofCalc), nBar[0]);
	    hMomYS12->Fill(momFromTime(0.938, 10.9, tofCalc), xToF[0]);
	  } // else if ( tofCalc > (tLight - (piLow+piHi)/2.) + proLow && tofCalc < (tLight - (piLow+piHi)/2.) + proHi )
	} // S1 + S2 trigger
      } // if (nhit == 1)
    } // for (int t=0; t<tree->GetEntries(); t++)

    fout->cd();
    hThetaS1S2pro->Write();
    hThetaS1S2pi->Write();
    hPhiS1S2pro->Write();
    hPhiS1S2pi->Write();
    hThetaS1pro->Write();
    hThetaS1pi->Write();
    hPhiS1pro->Write();
    hPhiS1pi->Write();
    hMomS1S2->Write();
    hMomS1->Write();

    hutof1dS1S2->Write();
    hutof1dS1->Write();

    hThetaS1S2ratio->Divide(hThetaS1S2pro, hThetaS1S2pi, 1., 1., "B");
    hPhiS1S2ratio->Divide(hPhiS1S2pro, hPhiS1S2pi, 1., 1., "B");
    hThetaS1ratio->Divide(hThetaS1pro, hThetaS1pi, 1., 1., "B");
    hPhiS1ratio->Divide(hPhiS1pro, hPhiS1pi, 1., 1., "B");
    hThetaS1S2ratio->Write();
    hPhiS1S2ratio->Write();
    hThetaS1ratio->Write();
    hPhiS1ratio->Write();

    if (nBlocks==0) {
      hThetaS1S2pro->SetLineColor(kBlue);
      hThetaS1S2pi->SetLineColor(kBlue);
      hPhiS1S2pro->SetLineColor(kBlue);
      hPhiS1S2pi->SetLineColor(kBlue);
      hPhiS1S2ratio->SetLineColor(kBlue);
      hThetaS1S2ratio->SetLineColor(kBlue);
      hThetaS1pro->SetLineColor(kBlue);
      hThetaS1pi->SetLineColor(kBlue);
      hPhiS1pro->SetLineColor(kBlue);
      hPhiS1pi->SetLineColor(kBlue);
      hPhiS1ratio->SetLineColor(kBlue);
      hThetaS1ratio->SetLineColor(kBlue);
      hMomS1S2->SetLineColor(kBlue);
      hMomS1->SetLineColor(kBlue);
      hutof1dS1S2->SetLineColor(kBlue);
      hutof1dS1->SetLineColor(kBlue);
      leg->AddEntry(hThetaS1S2pro, "0 blocks", "l");
      legTof->AddEntry(hutof1dS1, "0 blocks", "l");

      hMom2D_0blkQ->Scale(1. / (double)nSpills);
      hMom2D_0blkS->Scale(1. / (double)nSpills);
      hMom2D_0blkQ->Write();
      hMom2D_0blkS->Write();
    }
    if (nBlocks==1) {
      hThetaS1S2pro->SetLineColor(kRed);
      hThetaS1S2pi->SetLineColor(kRed);
      hPhiS1S2pro->SetLineColor(kRed);
      hPhiS1S2pi->SetLineColor(kRed);
      hPhiS1S2ratio->SetLineColor(kRed);
      hThetaS1S2ratio->SetLineColor(kRed);

      hThetaS1pro->SetLineColor(kRed);
      hThetaS1pi->SetLineColor(kRed);
      hPhiS1pro->SetLineColor(kRed);
      hPhiS1pi->SetLineColor(kRed);
      hPhiS1ratio->SetLineColor(kRed);
      hThetaS1ratio->SetLineColor(kRed);

      hMomS1S2->SetLineColor(kRed);
      hMomS1->SetLineColor(kRed);

      hutof1dS1S2->SetLineColor(kRed);
      hutof1dS1->SetLineColor(kRed);

      leg->AddEntry(hThetaS1S2pro, "1 block", "l");
      legTof->AddEntry(hutof1dS1, "1 block", "l");

      hMom2D_1blkQ->Scale(1. / (double)nSpills);
      hMom2D_1blkS->Scale(1. / (double)nSpills);
      hMom2D_1blkQ->Write();
      hMom2D_1blkS->Write();
    }
    if (nBlocks==2) {
      hThetaS1S2pro->SetLineColor(kBlack);
      hThetaS1S2pi->SetLineColor(kBlack);
      hPhiS1S2pro->SetLineColor(kBlack);
      hPhiS1S2pi->SetLineColor(kBlack);
      hPhiS1S2ratio->SetLineColor(kBlack);
      hThetaS1S2ratio->SetLineColor(kBlack);

      hThetaS1pro->SetLineColor(kBlack);
      hThetaS1pi->SetLineColor(kBlack);
      hPhiS1pro->SetLineColor(kBlack);
      hPhiS1pi->SetLineColor(kBlack);
      hPhiS1ratio->SetLineColor(kBlack);
      hThetaS1ratio->SetLineColor(kBlack);

      hMomS1S2->SetLineColor(kBlack);
      hMomS1->SetLineColor(kBlack);

      hutof1dS1S2->SetLineColor(kBlack);
      hutof1dS1->SetLineColor(kBlack);

      leg->AddEntry(hThetaS1S2pro, "2 blocks", "l");
      legTof->AddEntry(hutof1dS1, "2 blocks", "l");

      hMom2D_2blkQ->Scale(1. / (double)nSpills);
      hMom2D_2blkS->Scale(1. / (double)nSpills);
      hMom2D_2blkQ->Write();
      hMom2D_2blkS->Write();
    }
    if (nBlocks==3) {
      hThetaS1S2pro->SetLineColor(kGreen+2);
      hThetaS1S2pi->SetLineColor(kGreen+2);
      hPhiS1S2pro->SetLineColor(kGreen+2);
      hPhiS1S2pi->SetLineColor(kGreen+2);
      hPhiS1S2ratio->SetLineColor(kGreen+2);
      hThetaS1S2ratio->SetLineColor(kGreen+2);

      hThetaS1pro->SetLineColor(kGreen+2);
      hThetaS1pi->SetLineColor(kGreen+2);
      hPhiS1pro->SetLineColor(kGreen+2);
      hPhiS1pi->SetLineColor(kGreen+2);
      hPhiS1ratio->SetLineColor(kGreen+2);
      hThetaS1ratio->SetLineColor(kGreen+2);

      hMomS1S2->SetLineColor(kGreen+2);
      hMomS1->SetLineColor(kGreen+2);

      hutof1dS1S2->SetLineColor(kGreen+2);
      hutof1dS1->SetLineColor(kGreen+2);

      leg->AddEntry(hThetaS1S2pro, "3 blocks", "l");
      legTof->AddEntry(hutof1dS1, "3 blocks", "l");
    }
    if (nBlocks==4) {
      hThetaS1S2pro->SetLineColor(kMagenta);
      hThetaS1S2pi->SetLineColor(kMagenta);
      hPhiS1S2pro->SetLineColor(kMagenta);
      hPhiS1S2pi->SetLineColor(kMagenta);
      hPhiS1S2ratio->SetLineColor(kMagenta);
      hThetaS1S2ratio->SetLineColor(kMagenta);

      hThetaS1pro->SetLineColor(kMagenta);
      hThetaS1pi->SetLineColor(kMagenta);
      hPhiS1pro->SetLineColor(kMagenta);
      hPhiS1pi->SetLineColor(kMagenta);
      hPhiS1ratio->SetLineColor(kMagenta);
      hThetaS1ratio->SetLineColor(kMagenta);

      hMomS1S2->SetLineColor(kMagenta);
      hMomS1->SetLineColor(kMagenta);

      hutof1dS1S2->SetLineColor(kMagenta);
      hutof1dS1->SetLineColor(kMagenta);

      leg->AddEntry(hThetaS1S2pro, "4 blocks", "l");
      legTof->AddEntry(hutof1dS1, "4 blocks", "l");
    }

    hPhiS1S2pro->Scale(1. / (double)nSpills);
    hPhiS1S2pi->Scale(1. / (double)nSpills);
    hThetaS1S2pro->Scale(1. / (double)nSpills);
    hThetaS1S2pi->Scale(1. / (double)nSpills);

    hPhiS1pro->Scale(1. / (double)nSpills);
    hPhiS1pi->Scale(1. / (double)nSpills);
    hThetaS1pro->Scale(1. / (double)nSpills);
    hThetaS1pi->Scale(1. / (double)nSpills);

    hMomS1S2->Scale(1. / (double)nSpills);
    hMomS1->Scale(1. / (double)nSpills);

    hutof1dS1S2->Scale(1. / (double)nSpills);
    hutof1dS1->Scale(1. / (double)nSpills);

    hMomYS1->Scale(1. / (double)nSpills);
    hMomYS12->Scale(1. / (double)nSpills);
    hMomZS1->Scale(1. / (double)nSpills);
    hMomZS12->Scale(1. / (double)nSpills);

    hMomYS1->Write();
    hMomYS12->Write();
    hMomZS1->Write();
    hMomZS12->Write();

    hsThetaS1S2pro->Add(hThetaS1S2pro);
    hsThetaS1S2pi->Add(hThetaS1S2pi);
    hsPhiS1S2pro->Add(hPhiS1S2pro);
    hsPhiS1S2pi->Add(hPhiS1S2pi);
    hsThetaS1S2ratio->Add(hThetaS1S2ratio);
    hsPhiS1S2ratio->Add(hPhiS1S2ratio);

    hsThetaS1pro->Add(hThetaS1pro);
    hsThetaS1pi->Add(hThetaS1pi);
    hsPhiS1pro->Add(hPhiS1pro);
    hsPhiS1pi->Add(hPhiS1pi);
    hsThetaS1ratio->Add(hThetaS1ratio);
    hsPhiS1ratio->Add(hPhiS1ratio);

    hsutof1dS1S2->Add(hutof1dS1S2);
    hsutof1dS1->Add(hutof1dS1);

    hsMomS1S2->Add(hMomS1S2);
    hsMomS1->Add(hMomS1);
  } // for (int nBlocks = 0; nBlocks <= 4; nBlocks++) 

  fout->cd();
  leg->Write("leg");
  hsThetaS1S2pro->Write();
  hsThetaS1S2pi->Write();
  hsPhiS1S2pro->Write();
  hsPhiS1S2pi->Write();
  hsThetaS1S2ratio->Write();
  hsPhiS1S2ratio->Write();

  hsutof1dS1S2->Write();
  hsutof1dS1->Write();

  hsThetaS1pro->Write();
  hsThetaS1pi->Write();
  hsPhiS1pro->Write();
  hsPhiS1pi->Write();
  hsThetaS1ratio->Write();
  hsPhiS1ratio->Write();

  hsMomS1S2->Write();
  hsMomS1->Write();

  TCanvas *c1_Log1 = new TCanvas("c1_Log1");
  c1_Log1->SetLogy();
  hsThetaS1S2pro->Draw("hist e nostack");
  leg->Draw();
  c1_Log1->Print(Form("%s/thetaS12proLog.png", saveDir));
  c1_Log1->Print(Form("%s/thetaS12proLog.pdf", saveDir));
  TCanvas *c1_Log2 = new TCanvas("c1_Log2");
  c1_Log2->SetLogy();
  hsThetaS1S2pi->Draw("hist e nostack");
  leg->Draw();
  c1_Log2->Print(Form("%s/thetaS12piLog.png", saveDir));
  c1_Log2->Print(Form("%s/thetaS12piLog.pdf", saveDir));
  TCanvas *c1_Log3 = new TCanvas("c1_Log3");
  c1_Log3->SetLogy();
  hsPhiS1S2pro->Draw("hist e nostack");
  leg->Draw();
  c1_Log3->Print(Form("%s/phiS12proLog.png", saveDir));
  c1_Log3->Print(Form("%s/phiS12proLog.pdf", saveDir));
  TCanvas *c1_Log4 = new TCanvas("c1_Log4");
  c1_Log4->SetLogy();
  hsPhiS1S2pi->Draw("hist e nostack");
  leg->Draw();
  c1_Log4->Print(Form("%s/phiS12piLog.png", saveDir));
  c1_Log4->Print(Form("%s/phiS12piLog.pdf", saveDir));
  TCanvas *c1_Log5 = new TCanvas("c1_Log5");
  c1_Log5->SetLogy();
  hsThetaS1S2ratio->Draw("hist e nostack");
  leg->Draw();
  c1_Log5->Print(Form("%s/thetaS12ratioLog.png", saveDir));
  c1_Log5->Print(Form("%s/thetaS12ratioLog.pdf", saveDir));
  TCanvas *c1_Log6 = new TCanvas("c1_Log6");
  c1_Log6->SetLogy();
  hsPhiS1S2ratio->Draw("hist e nostack");
  leg->Draw();
  c1_Log6->Print(Form("%s/phiS12ratioLog.png", saveDir));
  c1_Log6->Print(Form("%s/phiS12ratioLog.pdf", saveDir));

  TCanvas *c1_1 = new TCanvas("c1_1");
  hsThetaS1S2pro->Draw("hist e nostack");
  leg->Draw();
  c1_1->Print(Form("%s/thetaS12pro.png", saveDir));
  c1_1->Print(Form("%s/thetaS12pro.pdf", saveDir));
  TCanvas *c1_2 = new TCanvas("c1_2");
  hsThetaS1S2pi->Draw("hist e nostack");
  leg->Draw();
  c1_2->Print(Form("%s/thetaS12pi.png", saveDir));
  c1_2->Print(Form("%s/thetaS12pi.pdf", saveDir));
  TCanvas *c1_3 = new TCanvas("c1_3");
  hsPhiS1S2pro->Draw("hist e nostack");
  leg->Draw();
  c1_3->Print(Form("%s/phiS12pro.png", saveDir));
  c1_3->Print(Form("%s/phiS12pro.pdf", saveDir));
  TCanvas *c1_4 = new TCanvas("c1_4");
  hsPhiS1S2pi->Draw("hist e nostack");
  leg->Draw();
  c1_4->Print(Form("%s/phiS12pi.png", saveDir));
  c1_4->Print(Form("%s/phiS12pi.pdf", saveDir));
  TCanvas *c1_5 = new TCanvas("c1_5");
  hsThetaS1S2ratio->Draw("hist e nostack");
  leg->Draw();
  c1_5->Print(Form("%s/thetaS12ratio.png", saveDir));
  c1_5->Print(Form("%s/thetaS12ratio.pdf", saveDir));
  TCanvas *c1_6 = new TCanvas("c1_6");
  hsPhiS1S2ratio->Draw("hist e nostack");
  leg->Draw();
  c1_6->Print(Form("%s/phiS12ratio.png", saveDir));
  c1_6->Print(Form("%s/phiS12ratio.pdf", saveDir));

  TCanvas *c1_Logs11 = new TCanvas("c1_Logs11");
  c1_Logs11->SetLogy();
  hsThetaS1pro->Draw("hist e nostack");
  legTof->Draw();
  c1_Logs11->Print(Form("%s/thetaS1proLog.png", saveDir));
  c1_Logs11->Print(Form("%s/thetaS1proLog.pdf", saveDir));
  TCanvas *c1_Logs12 = new TCanvas("c1_Logs12");
  c1_Logs12->SetLogy();
  hsThetaS1pi->Draw("hist e nostack");
  legTof->Draw();
  c1_Logs12->Print(Form("%s/thetaS1piLog.png", saveDir));
  c1_Logs12->Print(Form("%s/thetaS1piLog.pdf", saveDir));
  TCanvas *c1_Logs13 = new TCanvas("c1_Logs13");
  c1_Logs13->SetLogy();
  hsPhiS1pro->Draw("hist e nostack");
  legTof->Draw();
  c1_Logs13->Print(Form("%s/phiS1proLog.png", saveDir));
  c1_Logs13->Print(Form("%s/phiS1proLog.pdf", saveDir));
  TCanvas *c1_Logs14 = new TCanvas("c1_Logs14");
  c1_Logs14->SetLogy();
  hsPhiS1pi->Draw("hist e nostack");
  legTof->Draw();
  c1_Logs14->Print(Form("%s/phiS1piLog.png", saveDir));
  c1_Logs14->Print(Form("%s/phiS1piLog.pdf", saveDir));
  TCanvas *c1_Logs15 = new TCanvas("c1_Logs15");
  c1_Logs15->SetLogy();
  hsThetaS1ratio->Draw("hist e nostack");
  leg->Draw();
  c1_Logs15->Print(Form("%s/thetaS1ratioLog.png", saveDir));
  c1_Logs15->Print(Form("%s/thetaS1ratioLog.pdf", saveDir));
  TCanvas *c1_Logs16 = new TCanvas("c1_Logs16");
  c1_Logs16->SetLogy();
  hsPhiS1ratio->Draw("hist e nostack");
  leg->Draw();
  c1_Logs16->Print(Form("%s/phiS1ratioLog.png", saveDir));
  c1_Logs16->Print(Form("%s/phiS1ratioLog.pdf", saveDir));

  TCanvas *c1_s11 = new TCanvas("c1_s11");
  hsThetaS1pro->Draw("hist e nostack");
  leg->Draw();
  c1_s11->Print(Form("%s/thetaS1pro.png", saveDir));
  c1_s11->Print(Form("%s/thetaS1pro.pdf", saveDir));
  TCanvas *c1_s12 = new TCanvas("c1_s12");
  hsThetaS1pi->Draw("hist e nostack");
  leg->Draw();
  c1_s12->Print(Form("%s/thetaS1pi.png", saveDir));
  c1_s12->Print(Form("%s/thetaS1pi.pdf", saveDir));
  TCanvas *c1_s13 = new TCanvas("c1_s13");
  hsPhiS1pro->Draw("hist e nostack");
  leg->Draw();
  c1_s13->Print(Form("%s/phiS1pro.png", saveDir));
  c1_s13->Print(Form("%s/phiS1pro.pdf", saveDir));
  TCanvas *c1_s14 = new TCanvas("c1_s14");
  hsPhiS1pi->Draw("hist e nostack");
  leg->Draw();
  c1_s14->Print(Form("%s/phiS1pi.png", saveDir));
  c1_s14->Print(Form("%s/phiS1pi.pdf", saveDir));
  TCanvas *c1_s15 = new TCanvas("c1_s15");
  hsThetaS1ratio->Draw("hist e nostack");
  leg->Draw();
  c1_s15->Print(Form("%s/thetaS1ratio.png", saveDir));
  c1_s15->Print(Form("%s/thetaS1ratio.pdf", saveDir));
  TCanvas *c1_s16 = new TCanvas("c1_s16");
  hsPhiS1ratio->Draw("hist e nostack");
  leg->Draw();
  c1_s16->Print(Form("%s/phiS1ratio.png", saveDir));
  c1_s16->Print(Form("%s/phiS1ratio.pdf", saveDir));

  TCanvas *cMomS1S2 = new TCanvas("cMomS1S2");
  hsMomS1S2->Draw("hist e nostack");
  hsMomS1S2->GetXaxis()->SetLabelSize(0.05);
  hsMomS1S2->GetYaxis()->SetLabelSize(0.05);
  hsMomS1S2->GetXaxis()->SetTitleSize(0.05);
  hsMomS1S2->GetYaxis()->SetTitleSize(0.05);
  leg->Draw();
  cMomS1S2->Print(Form("%s/proMomS1S2.png", saveDir));
  cMomS1S2->Print(Form("%s/proMomS1S2.pdf", saveDir));

  TCanvas *cMomS1 = new TCanvas("cMomS1");
  hsMomS1->Draw("hist e nostack");
  hsMomS1->GetXaxis()->SetLabelSize(0.05);
  hsMomS1->GetYaxis()->SetLabelSize(0.05);
  hsMomS1->GetXaxis()->SetTitleSize(0.05);
  hsMomS1->GetYaxis()->SetTitleSize(0.05);
  leg->Draw();
  cMomS1->Print(Form("%s/proMomS1.png", saveDir));
  cMomS1->Print(Form("%s/proMomS1.pdf", saveDir));

  TCanvas *cutofS1S2 = new TCanvas("cutofS1S2");
  hsutof1dS1S2->Draw("hist e nostack");
  hsutof1dS1S2->GetXaxis()->SetLabelSize(0.05);
  hsutof1dS1S2->GetYaxis()->SetLabelSize(0.05);
  hsutof1dS1S2->GetXaxis()->SetTitleSize(0.05);
  hsutof1dS1S2->GetYaxis()->SetTitleSize(0.05);
  legTof->Draw();
  cutofS1S2->Print(Form("%s/utof1dS1S2.png", saveDir));
  cutofS1S2->Print(Form("%s/utof1dS1S2.pdf", saveDir));
  TCanvas *cutofS1 = new TCanvas("cutofS1");
  hsutof1dS1->Draw("hist e nostack");
  hsutof1dS1->GetXaxis()->SetLabelSize(0.05);
  hsutof1dS1->GetYaxis()->SetLabelSize(0.05);
  hsutof1dS1->GetXaxis()->SetTitleSize(0.05);
  hsutof1dS1->GetYaxis()->SetTitleSize(0.05);
  legTof->Draw();
  cutofS1->Print(Form("%s/utof1dS1.png", saveDir));
  cutofS1->Print(Form("%s/utof1dS1.pdf", saveDir));

  TCanvas *cutofS1S2Log = new TCanvas("cutofS1S2Log");
  cutofS1S2Log->SetLogy();
  hsutof1dS1S2->Draw("hist nostack");
  hsutof1dS1S2->GetXaxis()->SetLabelSize(0.05);
  hsutof1dS1S2->GetYaxis()->SetLabelSize(0.05);
  hsutof1dS1S2->GetXaxis()->SetTitleSize(0.05);
  hsutof1dS1S2->GetYaxis()->SetTitleSize(0.05);
  cutofS1S2Log->SetLeftMargin(0.13);
  cutofS1S2Log->SetBottomMargin(0.13);
  legTof->Draw();
  cutofS1S2Log->Print(Form("%s/utof1dS1S2Log.png", saveDir));
  cutofS1S2Log->Print(Form("%s/utof1dS1S2Log.pdf", saveDir));
  TCanvas *cutofS1Log = new TCanvas("cutofS1Log");
  cutofS1Log->SetLogy();
  hsutof1dS1->Draw("hist nostack");
  hsutof1dS1->GetXaxis()->SetLabelSize(0.05);
  hsutof1dS1->GetYaxis()->SetLabelSize(0.05);
  hsutof1dS1->GetXaxis()->SetTitleSize(0.05);
  hsutof1dS1->GetYaxis()->SetTitleSize(0.05);
  cutofS1Log->SetLeftMargin(0.13);
  cutofS1Log->SetBottomMargin(0.13);
  legTof->Draw();
  cutofS1Log->Print(Form("%s/utof1dS1Log.png", saveDir));
  cutofS1Log->Print(Form("%s/utof1dS1Log.pdf", saveDir));

  fout->Close();
} // angularDistS3
