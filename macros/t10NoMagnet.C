// t10NoMagnet.C

double momFromTime(const double mass, const double baseline, const double time)
{
  double mom = 0.;
  mom = mass * (baseline/(time*1e-9))*(1/TMath::Sqrt((pow(3e8*1e-9*time,2) - pow(baseline,2))/pow(time*1e-9,2)));
  return mom;
}

double keFromTime(const double mass, const double baseline, const double time)
{
  double mom = momFromTime(mass, baseline, time);
  double ke = TMath::Sqrt( pow(mom, 2) + pow(mass, 2) ) - mass;
  return ke;
}

void t10NoMagnet(const char* saveDir,
		 const char* utofDir="/nfs/scratch0/dbrailsf/data_backup/utof_backup_firsthitpinnedtounixtime/Data_root_v3_wo_walk_corr/")
{
  const int hitMax = 50;
  // Edges of S3 in beam coordinate system
  const double s3StartX = 0.601158;
  const double s3EndX   = -1.08086;
  const double s3s1StartY = 9.07239 + 1.767;
  const double s3s1EndY   = 8.89731 + 1.767;
  // z coordinates
  const double s3BarTop    = 62.; 
  const double s3BarBottom = -60.;
  // Time cut for hits interacting in neighbouring S3 bars
  const double twoBarCut = 0.5; // ns
  const double proLow  = 52;
  const double proHi = 80.;
  const double piLow = 35.75;
  const double piHi  = 37.75;
  double binsTheta[] = {-3.1, -3.,
			-2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, -2., 
			-1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.,
			-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 
			0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
			1.0, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 
			2.0, 2.2, 2.4, 2.6, 2.8, 
			3.0, 3.2, 3.4, 3.6, 3.8,  
			4.0, 4.25, 4.5, 4.75, 
			5.0, 5.25, 5.8,
		        6.}; 
  int binnum = sizeof(binsTheta)/sizeof(double) - 1;

  // SJ's bar by bar amplitude cut (by eye)
  // For A1 
  const std::vector<double> A1CutVec = {0.25, 0.25, 0.275, 0.2, 0.25, 0.25, 0.225, 0.25, 0.3, 
					0.3, 0.3,
					0.2, 0.2, 0.25, 0.25, 0.25, 0.25, 0.3, 0.25, 
					0.225, 0.3, 0.3};
  // For A2
  const std::vector<double> A2CutVec = {0.3, 0.275, 0.25, 0.125, 0.25, 
					0.25, 0.15, 0.225, 0.25, 0.25, 0.225,
					0.25, 0.25, 0.2, 0.2, 0.225, 
					0.225, 0.25, 0.25, 0.225, 0.3, 0.25};

  TH1D *hUtof = new TH1D("hUtof", "Time of flight; t_{S3} - t_{S1} / ns; Events", 250, 25, 125);
  hUtof->Sumw2();
  hUtof->SetLineWidth(2);
  hUtof->SetLineColor(kBlack);
  hUtof->GetXaxis()->SetLabelSize(.05);
  hUtof->GetYaxis()->SetLabelSize(.05);
  hUtof->GetXaxis()->SetTitleSize(.05);
  hUtof->GetYaxis()->SetTitleSize(.05);
  TH2D *hHitsXY = new TH2D("h2dHitsXY", "Distribution of hits within S3; x / m; y / m; Hits", 100, -1., 0.6, 22, -0.6, 0.61);
  hHitsXY->GetXaxis()->SetLabelSize(.05);
  hHitsXY->GetYaxis()->SetLabelSize(.05);
  hHitsXY->GetXaxis()->SetTitleSize(.05);
  hHitsXY->GetYaxis()->SetTitleSize(.05);
  hHitsXY->GetZaxis()->SetLabelSize(.05);
  hHitsXY->GetZaxis()->SetTitleSize(.05);
  TH1D *hMom = new TH1D("hMom", "Proton momentum measured in S3; Proton momentum [GeV/c]; Events", 70, 0.55, 0.91);
  hMom->Sumw2();
  hMom->SetLineWidth(2);
  hMom->SetLineColor(kBlack);
  hMom->GetXaxis()->SetLabelSize(.05);
  hMom->GetYaxis()->SetLabelSize(.05);
  hMom->GetXaxis()->SetTitleSize(.05);
  hMom->GetYaxis()->SetTitleSize(.05);
  TH1D *hKE = new TH1D("hKE", "Proton kinetic energy measured in S3; Proton kinetic energy / MeV; Events", 70, 150., 350.);
  hKE->Sumw2();
  hKE->SetLineWidth(2);
  hKE->SetLineColor(kBlack);
  hKE->GetXaxis()->SetLabelSize(.05);
  hKE->GetYaxis()->SetLabelSize(.05);
  hKE->GetXaxis()->SetTitleSize(.05);
  hKE->GetYaxis()->SetTitleSize(.05);
  TH1D *hThetaS1pro = new TH1D("hThetaS1pro", "Angular distribution of proton hits in S3 (S1 trigger only); #theta / degrees; Events / spill / degree", binnum, binsTheta);
  hThetaS1pro->Sumw2();
  hThetaS1pro->GetXaxis()->SetLabelSize(.05);
  hThetaS1pro->GetXaxis()->SetTitleSize(.05);
  hThetaS1pro->GetYaxis()->SetLabelSize(.05);
  hThetaS1pro->GetYaxis()->SetTitleSize(.05);
  TH1D *hPhiS1pro = new TH1D("hPhiS1pro", "Angular distribution of proton hits in S3 (S1 trigger only); #phi / degrees; Events / spill / degree",  22, -3.22, 3.35);
  hPhiS1pro->Sumw2();
  hPhiS1pro->GetXaxis()->SetLabelSize(.05);
  hPhiS1pro->GetXaxis()->SetTitleSize(.05);
  hPhiS1pro->GetYaxis()->SetLabelSize(.05);
  hPhiS1pro->GetYaxis()->SetTitleSize(.05);
  TH1D *hThetaS1pi = new TH1D("hThetaS1pi", "Angular distribution of pion hits in S3 (S1 trigger only); #theta / degrees; Events / spill / degree", binnum, binsTheta);
  hThetaS1pi->Sumw2();
  hThetaS1pi->GetXaxis()->SetLabelSize(.05);
  hThetaS1pi->GetXaxis()->SetTitleSize(.05);
  hThetaS1pi->GetYaxis()->SetLabelSize(.05);
  hThetaS1pi->GetYaxis()->SetTitleSize(.05);
  TH1D *hPhiS1pi = new TH1D("hPhiS1pi", "Angular distribution of pion hits in S3 (S1 trigger only); #phi / degrees; Events / spill / degree", 22, -3.22, 3.35);
  hPhiS1pi->Sumw2();
  hPhiS1pi->GetXaxis()->SetLabelSize(.05);
  hPhiS1pi->GetXaxis()->SetTitleSize(.05);
  hPhiS1pi->GetYaxis()->SetLabelSize(.05);
  hPhiS1pi->GetYaxis()->SetTitleSize(.05);
  TH1D *hThetaRatio = new TH1D("hThetaRatio", "S1 #cap S3 angular distribution of proton/MIP ratio; #theta / degrees; Protons/MIPs", binnum, binsTheta);
  hThetaRatio->Sumw2();
  hThetaRatio->SetLineColor(kBlack);
  hThetaRatio->SetLineWidth(2);
  hThetaRatio->GetXaxis()->SetLabelSize(.05);
  hThetaRatio->GetXaxis()->SetTitleSize(.05);
  hThetaRatio->GetYaxis()->SetLabelSize(.05);
  hThetaRatio->GetYaxis()->SetTitleSize(.05);
  TH1D *hPhiRatio = new TH1D("hPhiS1ratio", "S1 #cap S3 angular distribution of proton/MIP; #phi / degrees; Protons/MIPs", 22, -3.22, 3.35);
  hPhiRatio->Sumw2();
  hPhiRatio->SetLineColor(kBlack);
  hPhiRatio->SetLineWidth(2);
  hPhiRatio->GetXaxis()->SetLabelSize(.05);
  hPhiRatio->GetXaxis()->SetTitleSize(.05);
  hPhiRatio->GetYaxis()->SetLabelSize(.05);
  hPhiRatio->GetYaxis()->SetTitleSize(.05);

  // File with no magnet bend. Beam on nominal axis
  TFile *fin = new TFile(Form("%s/Data_2018_8_17_b8_800MeV_0block.root", utofDir), "read");
  TTree *tree = (TTree*)fin->Get("tree");

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
      // Checks if the next S3 hit is in a neighbouring bar
      // Only want to count one of them if it is
      // Check from this hit onwards to see if there are any double hit candidates
      bool isDouble = false;
      if (nh < nhit-1) {
	for (int nh2=nh+1; nh2<nhit; nh2++) {
	  // Timing cut and checks if in the neighbouring bar
	  if (abs(tToF[nh] - tToF[nh2]) < twoBarCut && 
	      (nBar[nh] == nBar[nh2]-1 || nBar[nh] == nBar[nh2]+1)) {
	    isDouble = true;
	    break;
	  }
	} // for (int nh2=nh; nh2<nhit; nh2++)
      } // if (nh < nh-1)
      if (!isDouble) {
	double tofCalc = tToF[nh] - tS1;
	// Calculate x, y z positions relative to S1
	double positionX = (xToF[nh]/168)*(s3EndX - s3StartX) + s3StartX;
	double positionY = (xToF[nh]/168.)*(s3s1EndY - s3s1StartY) + s3s1StartY;
	double positionZ = (yToF[nh] + s3BarBottom + 2.75) / 100.;
	double angleTheta = TMath::ATan(-positionX / positionY) * (180./TMath::Pi());
	double anglePhi   = TMath::ATan(positionZ / positionY) * (180./TMath::Pi());
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
	  hMom->Fill(momFromTime(0.938, 10.8, tofCalc));
	  hKE->Fill(keFromTime(0.938, 10.8, tofCalc)*1000);
	}
      } // Is not a double hit
    } // Loop over hits
  } // Loop over entries
  gStyle->SetOptStat(0);
  gStyle->SetPalette(55);

  hThetaRatio->Divide(hThetaS1pro, hThetaS1pi, 1., 1., "B");
  hPhiRatio->Divide(hPhiS1pro, hPhiS1pi, 1., 1., "B");

  TFile *fout = new TFile(Form("%s", saveDir), "recreate");
  fout->cd();
  hUtof->Write();
  hMom->Write();
  hKE->Write();
  hHitsXY->Write();
  hThetaRatio->Write();
  hPhiRatio  ->Write();
  hThetaS1pro->Write();
  hPhiS1pro  ->Write();
  hThetaS1pi ->Write();
  hPhiS1pi   ->Write();

  fout->Close();
  delete fout;
  fin->Close();
  delete fin;
}
