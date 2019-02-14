// angularDistS3.C

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
  const double proLow  = -8.;
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

  TFile *fout = new TFile(Form("%s/angularDistS3plots.root", saveDir), "recreate");

  TLegend *leg = new TLegend(0.15, 0.6, 0.3, 0.8);

  for (int nBlocks = 0; nBlocks <= 4; nBlocks++) {

    TH1D *hThetaS1pro   = new TH1D(Form("hThetaS1pro%d", nBlocks), Form("Angular distribution of proton hits in S3 (S1 trigger only), %d blocks; #theta / degrees; Events / spill", nBlocks), 100, -3.8, 6.2);
    TH1D *hThetaS1S2pro = new TH1D(Form("hThetaS1S2pro%d", nBlocks), Form("Angular distribution of proton hits in S3 (S1 & S2 triggers), %d blocks; #theta / degrees; Events / spill", nBlocks), 100, -3.8, 6.2);
    hThetaS1S2pro->Sumw2();
    TH1D *hPhiS1pro   = new TH1D(Form("hPhiS1pro%d", nBlocks), Form("Angular distribution of proton hits in S3 (S1 trigger only), %d blocks; #phi / degrees; Events / spill", nBlocks), 22, -3.15, 3.25);
    TH1D *hPhiS1S2pro = new TH1D(Form("hPhiS1S2pro%d", nBlocks), Form("Angular distribution of proton hits in S3 (S1 & S2 triggers), %d blocks; #phi / degrees; Events / spill", nBlocks), 22, -3.15, 3.25);
    hPhiS1S2pro->Sumw2();

    TH1D *hThetaS1pi   = new TH1D(Form("hThetaS1pi%d", nBlocks), Form("Angular distribution of pion hits in S3 (S1 trigger only), %d blocks; #theta / degrees; Events / spill", nBlocks), 100, -3.8, 6.2);
    TH1D *hThetaS1S2pi = new TH1D(Form("hThetaS1S2pi%d", nBlocks), Form("Angular distribution of pion hits in S3 (S1 & S2 triggers), %d blocks; #theta / degrees; Events / spill", nBlocks), 100, -3.8, 6.2);
    hThetaS1S2pi->Sumw2();
    TH1D *hPhiS1pi   = new TH1D(Form("hPhiS1pi%d", nBlocks), Form("Angular distribution of pion hits in S3 (S1 trigger only), %d blocks; #phi / degrees; Events / spill", nBlocks), 22, -3.15, 3.25);
    TH1D *hPhiS1S2pi = new TH1D(Form("hPhiS1S2pi%d", nBlocks), Form("Angular distribution of pion hits in S3 (S1 & S2 triggers), %d blocks; #phi / degrees; Events / spill", nBlocks), 22, -3.15, 3.25);
    hPhiS1S2pi->Sumw2();

    TH1D *hThetaS1ratio = new TH1D(Form("hThetaS1ratio%d", nBlocks), Form("Angular distribution of proton/MIP ratio in S3 (S1 trigger only), %d blocks; #phi / degrees; Protons/MIPs", nBlocks), 100, -3.8, 6.2);
    TH1D *hPhiS1ratio   = new TH1D(Form("hPhiS1ratio%d", nBlocks), Form("Angular distribution of proton/MIP in S3 (S1 trigger only), %d blocks; #phi / degrees; Protons/MIPs", nBlocks), 22, -3.15, 3.25);
    TH1D *hThetaS1S2ratio = new TH1D(Form("hThetaS1S2ratio%d", nBlocks), Form("Angular distribution of proton/MIP ratio in S3 (S1 & S2 triggers), %d blocks; #phi / degrees; Protons/MIPs", nBlocks), 100, -3.8, 6.2);
    TH1D *hPhiS1S2ratio   = new TH1D(Form("hPhiS1S2ratio%d", nBlocks), Form("Angular distribution of proton/MIP in S3 (S1 & S2 triggers), %d blocks; #phi / degrees; Protons/MIPs", nBlocks), 22, -3.15, 3.25);

    TH1D *hutof1dS1 = new TH1D(Form("hutof1dS1_%d",nBlocks), Form("Time of flight, %d blocks (S1 trigger only); S3 - S1 / ns; Events / spill", nBlocks), 260, 25, 155);
    TH1D *hutof1dS1S2 = new TH1D(Form("hutof1dS1S2_%d",nBlocks), Form("Time of flight, %d blocks (S1 & S2 trigger); S3 - S1 / ns; Events / spill", nBlocks), 260, 5, 155);

    int nP  = 0;
    int nPi = 0;

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
    int nhit;

    TTree *tree = (TTree*)futof->Get("tree");

    tree->SetBranchAddress("xToF", xToF);
    tree->SetBranchAddress("yToF", yToF);
    tree->SetBranchAddress("A1ToF", A1ToF);
    tree->SetBranchAddress("A2ToF", A2ToF);
    tree->SetBranchAddress("nhit", &nhit);
    tree->SetBranchAddress("tS1", &tS1);
    tree->SetBranchAddress("tToF", tToF);
    tree->SetBranchAddress("tTrig", &tTrig);

    for (int t=0; t<tree->GetEntries(); t++) {
      tree->GetEntry(t);
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
	  else if ( tofCalc > (tLight - (piLow+piHi)/2.) + proLow && tofCalc < (tLight - (piLow+piHi)/2.) + proHi &&
		    A1ToF[0] > ACut && A2ToF[0] > ACut) {
	    nP++;
	    hThetaS1S2pro->Fill(angleTheta);
	    hPhiS1S2pro->Fill(anglePhi);
	  } // else if ( tofCalc > (tLight - (piLow+piHi)/2.) + proLow && tofCalc < (tLight - (piLow+piHi)/2.) + proHi )
	} // S1 + S2 trigger
      } // if (nhit == 1)
    } // for (int t=0; t<tree->GetEntries(); t++)

    fout->cd();
    hThetaS1S2pro->Write();
    hThetaS1S2pi->Write();
    hPhiS1S2pro->Write();
    hPhiS1S2pi->Write();

    hThetaS1S2ratio->Divide(hThetaS1S2pro, hThetaS1S2pi, 1., 1., "B");
    hPhiS1S2ratio->Divide(hPhiS1S2pro, hPhiS1S2pi, 1., 1., "B");
    hThetaS1S2ratio->Write();
    hPhiS1S2ratio->Write();

    if (nBlocks==0) {
      hThetaS1S2pro->SetLineColor(kBlue);
      hThetaS1S2pi->SetLineColor(kBlue);
      hPhiS1S2pro->SetLineColor(kBlue);
      hPhiS1S2pi->SetLineColor(kBlue);
      hPhiS1S2ratio->SetLineColor(kBlue);
      hThetaS1S2ratio->SetLineColor(kBlue);
      leg->AddEntry(hThetaS1S2pro, "0 blocks", "l");
    }
    if (nBlocks==1) {
      hThetaS1S2pro->SetLineColor(kRed);
      hThetaS1S2pi->SetLineColor(kRed);
      hPhiS1S2pro->SetLineColor(kRed);
      hPhiS1S2pi->SetLineColor(kRed);
      hPhiS1S2ratio->SetLineColor(kRed);
      hThetaS1S2ratio->SetLineColor(kRed);
      leg->AddEntry(hThetaS1S2pro, "1 block", "l");
    }
    if (nBlocks==2) {
      hThetaS1S2pro->SetLineColor(kBlack);
      hThetaS1S2pi->SetLineColor(kBlack);
      hPhiS1S2pro->SetLineColor(kBlack);
      hPhiS1S2pi->SetLineColor(kBlack);
      hPhiS1S2ratio->SetLineColor(kBlack);
      hThetaS1S2ratio->SetLineColor(kBlack);
      leg->AddEntry(hThetaS1S2pro, "2 blocks", "l");
    }
    if (nBlocks==3) {
      hThetaS1S2pro->SetLineColor(kGreen+2);
      hThetaS1S2pi->SetLineColor(kGreen+2);
      hPhiS1S2pro->SetLineColor(kGreen+2);
      hPhiS1S2pi->SetLineColor(kGreen+2);
      hPhiS1S2ratio->SetLineColor(kGreen+2);
      hThetaS1S2ratio->SetLineColor(kGreen+2);
      leg->AddEntry(hThetaS1S2pro, "3 blocks", "l");
    }
    if (nBlocks==4) {
      hThetaS1S2pro->SetLineColor(kMagenta);
      hThetaS1S2pi->SetLineColor(kMagenta);
      hPhiS1S2pro->SetLineColor(kMagenta);
      hPhiS1S2pi->SetLineColor(kMagenta);
      hPhiS1S2ratio->SetLineColor(kMagenta);
      hThetaS1S2ratio->SetLineColor(kMagenta);
      leg->AddEntry(hThetaS1S2pro, "4 blocks", "l");
    }

    hsThetaS1S2pro->Add(hPhiS1S2pro);
    hsThetaS1S2pi->Add(hPhiS1S2pi);
    hsPhiS1S2pro->Add(hPhiS1S2pro);
    hsPhiS1S2pi->Add(hPhiS1S2pi);
    hsPhiS1S2ratio->Add(hPhiS1S2ratio);
    hsThetaS1S2ratio->Add(hThetaS1S2ratio);
  } // for (int nBlocks = 0; nBlocks <= 4; nBlocks++) 

  fout->cd();
  leg->Write();
  hsThetaS1S2pro->Write();
  hsThetaS1S2pi->Write();
  hsPhiS1S2pro->Write();
  hsPhiS1S2pi->Write();
  hsThetaS1S2ratio->Write();
  hsPhiS1S2ratio->Write();

  TCanvas *c1_1 = new TCanvas("c1_1");
  c1_1->SetLogy();
  hsThetaS1S2pro->Draw("hist e nostack");
  leg->Draw();
  c1_1->Print(Form("%s/thetaS12pro.png", saveDir));
  c1_1->Print(Form("%s/thetaS12pro.pdf", saveDir));
  TCanvas *c1_2 = new TCanvas("c1_2");
  c1_2->SetLogy();
  hsThetaS1S2pi->Draw("hist e nostack");
  leg->Draw();
  c1_2->Print(Form("%s/thetaS12pi.png", saveDir));
  c1_2->Print(Form("%s/thetaS12pi.pdf", saveDir));
  TCanvas *c1_3 = new TCanvas("c1_3");
  c1_3->SetLogy();
  hsPhiS1S2pro->Draw("hist e nostack");
  leg->Draw();
  c1_3->Print(Form("%s/phiS12pro.png", saveDir));
  c1_3->Print(Form("%s/phiS12pro.pdf", saveDir));
  TCanvas *c1_4 = new TCanvas("c1_4");
  c1_4->SetLogy();
  hsPhiS1S2pi->Draw("hist e nostack");
  leg->Draw();
  c1_4->Print(Form("%s/phiS12pi.png", saveDir));
  c1_4->Print(Form("%s/phiS12pi.pdf", saveDir));
  TCanvas *c1_5 = new TCanvas("c1_5");
  c1_5->SetLogy();
  hsThetaS1S2ratio->Draw("hist e nostack");
  leg->Draw();
  c1_5->Print(Form("%s/thetaS12ratio.png", saveDir));
  c1_5->Print(Form("%s/thetaS12ratio.pdf", saveDir));
  TCanvas *c1_6 = new TCanvas("c1_6");
  c1_6->SetLogy();
  hsPhiS1S2ratio->Draw("hist e nostack");
  leg->Draw();
  c1_6->Print(Form("%s/phiS12ratio.png", saveDir));
  c1_6->Print(Form("%s/phiS12ratio.pdf", saveDir));

  fout->Close();
} // angularDistS3
