// angularDistS4.C
// Angular distribution of protons and pions for different moderator blocks
// Use efficiency calculation and background subtraction

// Mass squared from the time
// momentum and mass in GeV
double massFromTime(const double time, const double mom, const double base) {
  double mass = pow(mom,2) * (pow((time*1e-9)/base, 2)*pow(3e8, 2) - 1);
  return mass;
}

void angularDistS4_newSample(const char* saveDir, 
			     const char* dstofDir="/nfs/scratch0/dbrailsf/data_backup/dtof_backup/",
			     const char* ustofDir="/nfs/scratch0/dbrailsf/data_backup/utof_backup_firsthitpinnedtounixtime/Data_root_v3_wo_walk_corr/",
			     const char* spillDBDir="/scratch0/sjones/spillDB/") 
{
 
  gROOT->SetBatch(kTRUE);
  // Unix timestamps for variable block moves
  // 0.8GeV/c, 0 blocks
  // 31/08/2018
  const double start0Block = 1535713289;
  const double end0Block   = 1535716132;
  // 0.8GeV/c, 1 block
  // 01/09/2018
  const double start1Block = 1535796057;
  const double end1Block   = 1535799112;
  // 0.8GeV/c, 2 blocks
  // 01/09/2018
  const double start2Block = 1535789157;
  const double end2Block   = 1535792026;
  // 0.8GeV/c, 3 block
  // 01/09/2018
  const double start3Block = 1535792404;
  const double end3Block   = 1535795300;
  // const double end3Block   = 1535798437;
  // 0.8GeV/c, 4 block
  // 4 moderator blocks with -4cm bend
  const double start4Block = 1535836129;
  const double end4Block   = 1535879634;
  // Most runs were in this configuration so don't need to use necessarily
  //  const double start4Block = 1536537600; 
  //  const double end4Block   = 1536669600;
  // Timing cuts
  const double piLow  = 36.;
  const double piHi   = 51.;
  const double proLow = 62.;
  const double proHi  = 101.;
  // S1 -> S4 positions
  const double baselineS1S4End   = 13.9426;
  const double baselineS1S4Start = 14.0069;
  const double s4OffAxisStartX   = 0.121;
  const double s4OffAxisEndX     = 1.4224;
  // Edges of S3 in beam coordinate system
  const double s3StartX   = -0.5168;
  const double s3EndX     = 0.9970;
  const double s3s1StartY = 9.0569 + 1.77;
  const double s3s1EndY   = 8.9146 + 1.77;
  // Define the runs to be used for varying number of blocks for ustof
  const char* str0Block = "Data_2018_8_31_b2_800MeV_0block.root";
  const char* str1Block = "Data_2018_9_1_b4_800MeV_1block_bend4cm.root";
  const char* str2Block = "Data_2018_9_1_b2_800MeV_2block_bend4cm.root";
  const char* str3Block = "Data_2018_9_1_b3_800MeV_3block_bend4cm.root";
  const char* str4Block = "Data_2018_9_1_b8_800MeV_4block_bend4cm.root";
  // New datasets to combine for the 4 block data
  const char* str4Block0 = "Data_2018_8_28_b5.root";
  const char* str4Block1 = "Data_2018_8_30_b1.root";
  const char* str4Block2 = "Data_2018_8_29_b4.root";
  const char* str4Block3 = "Data_2018_8_29_b1.root";
  std::vector<const char*> str4BlockVec = {str4Block0, str4Block1, str4Block2, str4Block3};
  // z coordinates
  const double s4BarTop    = 0.43; 
  const double s4BarBottom = -0.35;
  // Shift in ns required to to pion peak at speed of light
  const double dstofShift = 44.;
  // Ustof-dstof cable delay
  const double ustofDelay = 184.7;
  // Minimum number of hits in a spill required for it to count
  const int minHits = 200;
  // S2 to S4 distance - used for m^2 plot
  const double s2s4Dist = 12.5;

  TFile *fout = new TFile(Form("%s/angularDistS4Plots.root", saveDir), "recreate");

  THStack *hsBkgSub = new THStack("hsBkgSub", "Background subtracted and efficiency corrected S4 ToF spectrum; S4 - S2 / ns; Events / spill");
  THStack *hsDtof   = new THStack("hsDtof", "S4 ToF spectrum; S4 - S2 / ns; Events / spill");

  THStack *hsProS4Vert = new THStack("hsProS4Vert", "S1 #cap S2 #cap S4 angular distribution of proton hits; #phi / degrees; Events / spill");
  THStack *hsPiS4Vert  = new THStack("hsPiS4Vert", "S1 #cap S2 #cap S4 angular distribution of MIP hits; #phi / degrees; Events / spill");
  THStack *hsProS4Horz = new THStack("hsProS4Horz", "S1 #cap S2 #cap S4 angular distribution of proton hits; #theta / degrees; Events / spill");
  THStack *hsPiS4Horz  = new THStack("hsPiS4Horz", "S1 #cap S2 #cap S4 angular distribution of MIP hits; #theta / degrees; Events / spill");
  THStack *hsRatioS4Vert = new THStack("hsRatioS4Vert", "S1 #cap S2 #cap S4 angular distribution of proton/MIP ratio; #phi / degrees; Protons/MIPs");
  THStack *hsRatioS4Horz = new THStack("hsRatioS4Horz", "S1 #cap S2 #cap S4 angular distribution of proton/MIP ratio; #theta / degrees; Protons/MIPs");

  TLegend *leg = new TLegend(0.65, 0.55, 0.85, .85);

  TLegend *legProS4Horz = new TLegend(0.55, 0.55, 0.88, .85);
  TLegend *legPiS4Horz  = new TLegend(0.55, 0.55, 0.88, .85);
  TLegend *legProS4Vert = new TLegend(0.60, 0.55, 0.88, .85);
  TLegend *legPiS4Vert  = new TLegend(0.60, 0.55, 0.88, .85);

  TLegend *legRatioVert = new TLegend(0.15, 0.5, 0.4, 0.8);

  for (int nBlocks = 0; nBlocks <= 4; nBlocks++) {
    // 1D histograms in the horizontal and vertical direction of S4
    // Horizontal
    TH1D *hProS4Horz = new TH1D(Form("hProS4Horz%d",nBlocks), Form("Horizontal angular distribution of proton hits in S4, %d blocks; #theta / degrees; Events / spill",nBlocks), 20, 0., 6.);
    hProS4Horz->Sumw2();
    TH1D *hPiS4Horz  = new TH1D(Form("hPiS4Horz%d",nBlocks), Form("Horizontal angular distribution of MIP hits in S4, %d blocks; #theta / degrees; Events / spill",nBlocks), 20, 0., 6.);
    hPiS4Horz->Sumw2();
    TH1D *hProPiRatioS4Horz  = new TH1D(Form("hProPiRatioS4Horz%d",nBlocks), Form("Horizontal angular distribution of proton/MIP ratio in S4, %d blocks; #theta / degrees; Protons/MIPs",nBlocks), 20, 0., 6.);
    // Vertical
    TH1D *hProS4Vert = new TH1D(Form("hProS4Vert%d",nBlocks), Form("Vertical angular distribution of proton hits in S4, %d blocks; #phi / degrees; Events / spill",nBlocks), 10, -1.5, 1.8);
    hProS4Vert->Sumw2();
    TH1D *hPiS4Vert  = new TH1D(Form("hPiS4Vert%d",nBlocks), Form("Vertical angular distribution of MIP hits in S4, %d blocks; #phi / degrees; Events / spill",nBlocks), 10, -1.5, 1.8);
    hPiS4Vert->Sumw2();
    TH1D *hProPiRatioS4Vert  = new TH1D(Form("hProPiRatioS4Vert%d",nBlocks), Form("Vertical angular distribution of proton/MIP ratio in S4, %d blocks; #phi / degrees; Protons/MIPs",nBlocks), 10, -1.5, 1.8);
    TH1D *hAllS4Horz = new TH1D(Form("hAllS4Horz%d",nBlocks), Form("Horizontal angular distribution of hit in S4, %d blocks; #theta / degrees; Events / spill", nBlocks), 20, 0., 6.);
    TH1D *hMSq = new TH1D(Form("hMSq%d", nBlocks), Form("Particle mass distribution, %d blocks; M^{2} [GeV^{2} / c^{2}]; Events", nBlocks), 260, -0.5, 4.5);
    hMSq->Sumw2();

    if (nBlocks != 4) {
      // Number of signal particles using just cut and count
      double nP  = 0.;
      double nPi = 0.;
      // Define signal and background functions to be fitted
      // Signals are gaussians
      TF1 *sPro = new TF1("sPro", "gaus", proLow, proHi);
      TF1 *sPi  = new TF1("sPi", "gaus", piLow, piHi);
      // Exponential background
      TF1 *fBkgExp = new TF1("fBkgExp","expo", 30, 160);
      sPro->SetLineColor(kGreen+2);
      sPi->SetLineColor(kRed);
      TF1 *fSplusBExp = new TF1("signal+bkg exp", "gaus(0)+gaus(3)+expo(6)", 30, 160);
      fSplusBExp->SetParNames("const 1", "mean 1", "sigma 1",
			      "const 2", "mean 2", "sigma 2",
			      "bkgconst", "bkgdecay");
      fSplusBExp->SetLineColor(kBlack);
      // For spill counting normalisation
      int nSpills = 0;
      int nSpillsTrue = 0;
      double lastSpill = 0.;
      // 1D ToF for background subtraction
      TH1D *hdtof1d = new TH1D(Form("hdtof1d_%d",nBlocks), Form("Time of flight, %d blocks; S4 - S1 / ns; Events", nBlocks), 260, 30, 160);
      // Cosmics hists
      TH2D *h2Cosmics = new TH2D(Form("h2Cosmics%d",nBlocks), Form("Cosmic flux, %d blocks; x / cm; Bar; Hz",nBlocks), 20, 0, 140, 10, 0.5, 10.5);
      TH1D *hCosmicsVert = new TH1D(Form("hCosmicsVert%d",nBlocks), Form("Cosmic flux, %d blocks; x / cm; Hz",nBlocks), 10, 0.5, 10.5);
      hCosmicsVert->Sumw2();
      TH1D *hCosmicsHorz = new TH1D(Form("hCosmicsHorz%d",nBlocks), Form("Cosmic flux, %d blocks; x / cm; Hz",nBlocks), 20, 0, 140);
      hCosmicsHorz->Sumw2();
      // Find the correct dstof files
      Int_t runMin=-1;
      Int_t runMax=-1;

      double startTime = 0;
      double endTime   = 0;

      if (nBlocks == 0) {
	startTime = start0Block;
	endTime   = end0Block;
      }
      else if (nBlocks == 1) {
	startTime = start1Block;
	endTime   = end1Block;
      }
      else if (nBlocks == 2) {
	startTime = start2Block;
	endTime   = end2Block;
      }
      else if (nBlocks == 3) {
	startTime = start3Block;
	endTime   = end3Block;
      }
      else if (nBlocks == 4) {
	startTime = start4Block;
	endTime   = end4Block;
      }

      for (int irun=950; irun<1400; irun++) {
	TFile *fin = new TFile(Form("%srun%d/DsTOFcoincidenceRun%d_tdc1.root", dstofDir, irun, irun), "read");
	RawDsTofCoincidence *tofCoinTemp = NULL;
	TTree *tree = (TTree*) fin->Get("tofCoinTree");
	tree->SetBranchAddress("tofCoin", &tofCoinTemp);
	tree->GetEntry(0);
	UInt_t firstTemp = tofCoinTemp->unixTime[0];
	tree->GetEntry(tree->GetEntries()-1);
	UInt_t lastTemp = tofCoinTemp->unixTime[0];
      
	fin->Close();
	delete fin;
      
	if (firstTemp>endTime){
	  break;
	}
      
	if (firstTemp<startTime && lastTemp>startTime){
	  runMin = irun;
	}
      
	if (firstTemp<endTime && lastTemp>endTime){
	  runMax = irun;
	}   
      } // for (int irun=950; irun<1400; irun++) 
    
      cout << "Min and max runs are " << runMin << " " << runMax << endl;

      // Need these to calculate bar efficiencies
      TH1D *hCoins = new TH1D(Form("hCoins_%d",nBlocks), Form("Bar coincidences + S_{1,2} coincidences, %d blocks; Bar; Events",nBlocks), 10, 0.5, 10.5);
      TH1D *hHits  = new TH1D(Form("hHits_%d",nBlocks), Form("PMT hits + S_{1,2} coincidences, %d blocks; Bar; Events",nBlocks), 10, 0.5, 10.5);
      TH1D *hEff   = new TH1D(Form("hEff_%d",nBlocks), Form("S4 bar efficiencies %d blocks; Bar; Events",nBlocks), 10, 0.5, 10.5);
      hHits->Sumw2();
      hCoins->Sumw2();

      // In this loop calculate the bar-by-bar efficiencies
      for (int itdc=0; itdc<2; itdc++) {
	// Need both the coincidence and the raw trees for this
	double tempUstof;
	double ustofNs;
	for (int irun=runMin; irun<runMax+1; irun++) {
	  // Load input files
	  TFile *tofCoinFile = new TFile(Form("%srun%d/DsTOFcoincidenceRun%d_tdc%d.root", dstofDir, irun, irun, itdc+1));
	  TFile *tofFile     = new TFile(Form("%srun%d/DsTOFtreeRun%d_tdc%d.root", dstofDir, irun, irun, itdc+1));
	  TTree *tofCoinTree = (TTree*)tofCoinFile->Get("tofCoinTree");
	  TTree *tofTree = (TTree*)tofFile->Get("tofTree");
	  RawDsTofCoincidence *tofCoin = NULL;
	  RawDsTofHeader *tof = NULL;
	  tofCoinTree->SetBranchAddress("tofCoin", &tofCoin);
	  tofTree->SetBranchAddress("tof", &tof);
	  // Create new TTree for ustof signals
	  TTree *ustofTree = new TTree("ustofTree", "ustof");
	  ustofTree->SetDirectory(0);
	  ustofTree->Branch("ustofNs", &ustofNs, "ustofNs/D");
	  tempUstof=0;
	  // Build tree of all the ustof hit times
	  for (int i=0; i<tofTree->GetEntries(); i++) {
	    tofTree->GetEntry(i);
	    if (tof->unixTime < startTime) continue;
	    if (tof->unixTime > endTime) break;

	    if (tof->channel == 13) {
	      ustofNs = tof->fakeTimeNs - ustofDelay;
	      // ustof signals shouldn't be coming closer than 500ns
	      if ( (ustofNs - tempUstof) < 500.) continue;
	      ustofTree->Fill();
	      tempUstof = ustofNs;
	    } // if (tof->channel == 13) 
	  } // for (int i=0; i<tof->GetEntries(); i++)
	  ustofTree->BuildIndex("ustofNs");
	  // Now loop over coincidence file to find number of coincidences for each bar
	  for (int t=0; t<tofCoinTree->GetEntries(); t++) {
	    tofCoinTree->GetEntry(t);
	    if (tofCoin->unixTime[0] < startTime) continue;
	    if (tofCoin->unixTime[0] > endTime) break;
	    double deltat = TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1]);
	    double dstofHitT = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (10. - TMath::Abs(deltat) / 2. );
	    double tofCalc = dstofHitT - tofCoin->usTofSignal;
	    if (tofCalc > 70. && tofCalc < 200. && tofCoin->bar != 10) {
	      hCoins->Fill(tofCoin->bar);
	    } // if (tofCalc > 70. && tofCalc < 200.)
	  } // for (int t=0; t<tofCoinTree->GetEntries(); t++) 
	  // Loop over hit file to find number of hits in coincidence 
	  for (int t=0; t<tofTree->GetEntries(); t++) {
	    tofTree->GetEntry(t);
	    if (tof->unixTime < startTime) continue;
	    if (tof->unixTime > endTime) break;
	    if (tof->channel > 10 || tof->channel == 0) continue;
	    int ustofEntry = ustofTree->GetEntryNumberWithBestIndex(tof->fakeTimeNs);
	    ustofTree->GetEntry(ustofEntry);

	    if (tof->fakeTimeNs - ustofNs < 260. && tof->fakeTimeNs - ustofNs > 130.) {
	      if (itdc == 0) { // TDC1
		if (tof->channel % 2 == 1) { // PMT B
		  hHits->Fill(((tof->channel + 1)/ -2) + 11);
		} // if (tof->channel % 2 == 1)
		else { // PMT A
		  hHits->Fill((tof->channel / -2) + 11);
		}
	      } // TDC 1
	      else { // TDC 2
		if (tof->channel % 2 == 1) { // PMT B
		  hHits->Fill(((tof->channel + 1)/ -2) + 6);
		} // if (tof->channel % 2 == 1)
		else { // PMT A
		  hHits->Fill((tof->channel / -2) + 6);
		}
	      } // TDC 2
	    } // if (tof->fakeTimeNs - ustofNs < 260. && tof->fakeTimeNs - ustofNs > 130.)
	  } // for (int t=0; t<tofTree->GetEntries(); t++)
	  delete tof;
	  delete tofCoin;
	  delete tofTree;
	  delete tofCoinTree;
	  delete ustofTree;
	  tofFile->Close();
	  tofCoinFile->Close();
	} // for (int irun=runMin; irun<runMax+1; irun++) 
      } // for (int itdc=0; itdc<2; itdc++)
      fout->cd();
      TCanvas *cHits = new TCanvas(Form("%d_cHits",nBlocks));
      hHits->Draw("hist e");
      hHits->Write();
      cHits->Print(Form("%s/%d_barHits.png",saveDir,nBlocks));
      cHits->Print(Form("%s/%d_barHits.pdf",saveDir,nBlocks));
      TCanvas *cCoins = new TCanvas(Form("%d_cCoins",nBlocks));
      hCoins->Draw("hist e");
      hCoins->Write();
      cCoins->Print(Form("%s/%d_barCoins.png",saveDir,nBlocks));
      cCoins->Print(Form("%s/%d_barCoins.pdf",saveDir,nBlocks));
      TCanvas *cEff = new TCanvas(Form("%d_cEff",nBlocks));
      hHits->Add(hCoins, -1.);
      hEff->Divide(hCoins, hHits, 1., 1., "B");
      hEff->Draw("hist e");
      hEff->Write();
      cEff->Print(Form("%s/%d_barEff.png",saveDir,nBlocks));
      cEff->Print(Form("%s/%d_barEff.pdf",saveDir,nBlocks));
      cEff->Print(Form("%s/%d_barEff.tex",saveDir,nBlocks));

      // Now loop over the coincidence files again and calculate the angular distributions
      for (int itdc=0; itdc<2; itdc++) {
	// First of all tchain the relevant files together
	TChain *tofCoinChain = new TChain("tofCoinTree");
	for (int irun=runMin; irun<runMax+1; irun++){
	  tofCoinChain->Add(Form("%srun%d/DsTOFcoincidenceRun%d_tdc%d.root", dstofDir, irun, irun, itdc+1));
	} // for (int irun=runMin; irun<runMax+1; irun++)
	RawDsTofCoincidence *tofCoin = NULL;
	tofCoinChain->SetBranchAddress("tofCoin", &tofCoin);
	int lasth = 0;
	// Use the spills recorded in the spill DB
	for (int h=0; h<tofCoinChain->GetEntries(); h++) {
	  tofCoinChain->GetEntry(h);
	  if (tofCoin->unixTime[0]<startTime) continue;
	  if (tofCoin->unixTime[0]>endTime) break;
	  if (tofCoin->lastDelayedBeamSignal != lastSpill && itdc == 0) {
	    lastSpill = tofCoin->lastDelayedBeamSignal;
	    nSpills++;
	    int nTrueHits = 0;
	    for (int sp=h; sp<tofCoinChain->GetEntries(); sp++) {
	      tofCoinChain->GetEntry(sp);
	      if (tofCoin->unixTime[0]<startTime) continue;
	      if (tofCoin->unixTime[0]>endTime) break;
	      if (tofCoin->fakeTimeNs[0] > lastSpill + 3e9) break;

	      if ((tofCoin->fakeTimeNs[0] - tofCoin->usTofSignal) < 200. &&
		  (tofCoin->fakeTimeNs[0] - tofCoin->usTofSignal) > 70. &&
		  (tofCoin->fakeTimeNs[0] - lastSpill) < 1e9 &&
		  (tofCoin->fakeTimeNs[0] - lastSpill) > 0.) {
		nTrueHits++;
		if (nTrueHits > 25) {
		  nSpillsTrue++;
		  break;
		}
	      }
	    } // for (int sp=h; sp<tofCoinChain->GetEntries(); sp++) 
	  } // if (tofCoin->lastDelayedBeamSignal != lastSpill && itdc == 0) 
	  // Need to calculate total signal hits here
	  // Weight these by the efficiency of calculated above
	  double deltat = TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1]  );
	  double dstofHitT = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (10. - TMath::Abs(deltat) / 2. );
	  double tofCalc = dstofHitT - tofCoin->usTofSignal - dstofShift;
	  if (tofCalc < 160. && tofCalc > 30. && tofCoin->bar != 10) {
	    hdtof1d->Fill(tofCalc, 1. / hEff->GetBinContent(tofCoin->bar));
	    hMSq->Fill(massFromTime(tofCalc, 0.8, s2s4Dist), 1. / hEff->GetBinContent(tofCoin->bar));
	    double positionXP = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 70.));
	    if (tofCalc < piHi & tofCalc > piLow) { 
	      nPi += (1. / /*h2CosEff->GetBinContent( h2CosEff->GetXaxis()->FindBin(positionXP), tofCoin->bar)*/hEff->GetBinContent(tofCoin->bar));
	      // Calculate position of hit in global coordinates
	      double positionX = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 65.) / 130.) * (s4OffAxisEndX - s4OffAxisStartX) + s4OffAxisStartX;
	      double positionY = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 65.) / 130.) * (baselineS1S4End - baselineS1S4Start) + baselineS1S4Start;
	      double positionZ = (tofCoin->bar * 0.075) - 0.375 - 0.01;
	      // Calculate the angles relative to the nominal beamline
	      double angleTheta = TMath::ATan(positionX / positionY) * (180./TMath::Pi());
	      double anglePhi   = TMath::ATan(positionZ / positionY) * (180./TMath::Pi());
	      hPiS4Horz->Fill(angleTheta, 1. / /*h2CosEff->GetBinContent( h2CosEff->GetXaxis()->FindBin(positionXP), tofCoin->bar)*/ hEff->GetBinContent(tofCoin->bar));
	      hPiS4Vert->Fill(anglePhi,   1. / /*h2CosEff->GetBinContent( h2CosEff->GetXaxis()->FindBin(positionXP), tofCoin->bar)*/hEff->GetBinContent(tofCoin->bar));
	    } // if (tofCalc < piHi & tofCalc > piLow)
	    else if (tofCalc < proHi & tofCalc > proLow) {
	      nP += (1. / /*h2CosEff->GetBinContent( h2CosEff->GetXaxis()->FindBin(positionXP), tofCoin->bar)*/hEff->GetBinContent(tofCoin->bar));
	      // Calculate position of hit in global coordinates
	      double positionX = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 65.) / 130.) * (s4OffAxisEndX - s4OffAxisStartX) + s4OffAxisStartX;
	      double positionY = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 65.) / 130.) * (baselineS1S4End - baselineS1S4Start) + baselineS1S4Start;
	      double positionZ = (tofCoin->bar * 0.075) - 0.375 - 0.01;
	      // Calculate the angles relative to the nominal beamline
	      double angleTheta = TMath::ATan(positionX / positionY) * (180./TMath::Pi());
	      double anglePhi   = TMath::ATan(positionZ / positionY) * (180./TMath::Pi());
	      hProS4Horz->Fill(angleTheta, 1. / /*h2CosEff->GetBinContent( h2CosEff->GetXaxis()->FindBin(positionXP), tofCoin->bar)*/hEff->GetBinContent(tofCoin->bar));
	      hProS4Vert->Fill(anglePhi,   1. / /*h2CosEff->GetBinContent( h2CosEff->GetXaxis()->FindBin(positionXP), tofCoin->bar)*/hEff->GetBinContent(tofCoin->bar));
	    } // else if (tofCalc < proHi & tofCalc > proLow) 
	  } // if (tofCalc < 160. && tofCalc > 30.) 
	} // for (int h=0; h<tofCoinChain->GetEntries(); h++)
	delete tofCoin;
	delete tofCoinChain;
      }
      fout->cd();
      TCanvas *c2d_exp = new TCanvas(Form("%d_c2d_exp",nBlocks));
      c2d_exp->SetLogy();
      hdtof1d->Fit(sPi, "R");
      hdtof1d->Fit(sPro, "R");
      hdtof1d->Fit(fBkgExp, "R");
      Double_t parExp[8];
      sPro->GetParameters(&parExp[0]);
      sPi->GetParameters(&parExp[3]);
      fBkgExp->GetParameters(&parExp[6]);
      fSplusBExp->SetParameters(parExp);
      hdtof1d->Fit(fSplusBExp, "R");
      hdtof1d->Draw("hist");
      fSplusBExp->Draw("same");
      hdtof1d->Write();
      c2d_exp->Print(Form("%s/%d_dtof1d.png",saveDir,nBlocks));
      c2d_exp->Print(Form("%s/%d_dtof1d.pdf",saveDir,nBlocks));
      c2d_exp->Print(Form("%s/%d_dtof1d.tex",saveDir,nBlocks));

      // Now we have the fit values, loop over again and subtract the background
      // For each bin, find the fraction of each particle type which are background and 
      // this fraction from the bin
      TF1 *fSub = new TF1("fSub", "exp([0]+[1]*x)", 30, 160);
      fSub->SetParameter(0, fSplusBExp->GetParameter("bkgconst"));  
      fSub->SetParameter(1, fSplusBExp->GetParameter("bkgdecay"));
      // Background subtracted tof spectrum
      TH1D *hdtof1d_sub = (TH1D*)hdtof1d->Clone(Form("hdtof1d_sub_%d",nBlocks));
      hdtof1d_sub->Add(fSub, -1.);
      hdtof1d_sub->Scale(1. / nSpillsTrue);
      // If the bin content drops below 0, set to 0
      for (int b = 0; b <=  hdtof1d_sub->GetNbinsX(); b++) {
	if (hdtof1d_sub->GetBinContent(b) < 0.) {
	  hdtof1d_sub->SetBinContent(b, 0);
	} // if (hdtof1d_sub->GetBinContent(b) < 0.)
      } // for (int b = 0; hdtof1d_sub->GetNbinsX(); b++)
      TCanvas *cdtof_sub = new TCanvas(Form("%d_cdtof_sub",nBlocks));

      // Integrate background function between proton and pion windows and then subtract
      double piBkg  = fSub->Integral(piLow,  piHi) / (0.5*nSpillsTrue);
      double proBkg = fSub->Integral(proLow, proHi) / (0.5*nSpillsTrue);
      cout<<"Pion backgrounds "<<piBkg<<", proton backgrounds "<<proBkg<<endl;

      cdtof_sub->SetLogy();
      hdtof1d_sub->Draw("hist");
      hdtof1d_sub->GetXaxis()->SetLabelSize(0.05);
      hdtof1d_sub->GetYaxis()->SetLabelSize(0.05);
      hdtof1d_sub->GetXaxis()->SetTitleSize(0.05);
      hdtof1d_sub->GetYaxis()->SetTitleSize(0.05);
      hdtof1d_sub->Write();
      cdtof_sub->Print(Form("%s/%d_dtof1d_bkgSub.png", saveDir, nBlocks));
      cdtof_sub->Print(Form("%s/%d_dtof1d_bkgSub.pdf", saveDir, nBlocks));

      cout<<"Spills "<<nSpills<<" ("<<nSpillsTrue<<" true)"<<endl;

      // Subtract background hits
      TCanvas *cProS4Horz = new TCanvas(Form("cProS4Horz%d",nBlocks));
      hProS4Horz->Scale(1. / nSpillsTrue);
      for (int b=1; b < 21; b++) {
	hProS4Horz->SetBinContent(b, hProS4Horz->GetBinContent(b) - (proBkg/20.));
	if (hProS4Horz->GetBinContent(b) < 0.) {
	  hProS4Horz->SetBinContent(b, 0);
	  hProS4Horz->SetBinError(b, 0);
	}
      }
      hProS4Horz->Draw("hist e");
      cProS4Horz->Print(Form("%s/%d_proS4Horz.png",saveDir,nBlocks));
      cProS4Horz->Print(Form("%s/%d_proS4Horz.pdf",saveDir,nBlocks));
      cProS4Horz->Print(Form("%s/%d_proS4Horz.tex",saveDir,nBlocks));
      TCanvas *cPiS4Horz  = new TCanvas(Form("cPiS4Horz%d",nBlocks));
      hPiS4Horz->Scale(1. / nSpillsTrue);
      for (int b=1; b < 21; b++) {
	hPiS4Horz->SetBinContent(b, hPiS4Horz->GetBinContent(b) - (piBkg/20.));
	if (hPiS4Horz->GetBinContent(b) < 0.) {
	  hPiS4Horz->SetBinError(b, 0);
	  hPiS4Horz->SetBinContent(b, 0);
	}
      }
      hPiS4Horz->Draw("hist e");
      cPiS4Horz->Print(Form("%s/%d_piS4Horz.png",saveDir,nBlocks));
      cPiS4Horz->Print(Form("%s/%d_piS4Horz.pdf",saveDir,nBlocks));
      cPiS4Horz->Print(Form("%s/%d_piS4Horz.tex",saveDir,nBlocks));
      TCanvas *cProS4Vert = new TCanvas(Form("cProS4Vert%d",nBlocks));
      hProS4Vert->Scale(1. / nSpillsTrue);
      for (int b=1; b < 10; b++) {
	hProS4Vert->SetBinContent(b, hProS4Vert->GetBinContent(b) - (proBkg/9.));
	if (hProS4Vert->GetBinContent(b) < 0.) {
	  hProS4Vert->SetBinError(b, 0);
	  hProS4Vert->SetBinContent(b, 0);
	}
      }
      hProS4Vert->Draw("hist e");
      cProS4Vert->Print(Form("%s/%d_proS4Vert.png",saveDir,nBlocks));
      cProS4Vert->Print(Form("%s/%d_proS4Vert.pdf",saveDir,nBlocks));
      cProS4Vert->Print(Form("%s/%d_proS4Vert.tex",saveDir,nBlocks));
      TCanvas *cPiS4Vert  = new TCanvas(Form("cPiS4Vert%d",nBlocks));
      hPiS4Vert->Scale(1. / nSpillsTrue);
      for (int b=1; b < 10; b++) {
	hPiS4Vert->SetBinContent(b, hPiS4Vert->GetBinContent(b) - (piBkg/9.));
	if (hPiS4Vert->GetBinContent(b) < 0.) {
	  hPiS4Vert->SetBinError(b, 0);
	  hPiS4Vert->SetBinContent(b, 0);
	}
      }
      hPiS4Vert->Draw("hist e");
      cPiS4Vert->Print(Form("%s/%d_piS4Vert.png",saveDir,nBlocks));
      cPiS4Vert->Print(Form("%s/%d_piS4Vert.pdf",saveDir,nBlocks));
      cPiS4Vert->Print(Form("%s/%d_piS4Vert.tex",saveDir,nBlocks));

      TCanvas *cRatioS4Vert = new TCanvas(Form("cRatioS4Vert%d",nBlocks));
      hProPiRatioS4Vert->Divide(hProS4Vert, hPiS4Vert, 1., 1., "B");
      hProPiRatioS4Vert->Draw("hist e");
      cRatioS4Vert->Print(Form("%s/%d_proPiS4Vert.png", saveDir,nBlocks));
      cRatioS4Vert->Print(Form("%s/%d_proPiS4Vert.pdf", saveDir,nBlocks));
      cRatioS4Vert->Print(Form("%s/%d_proPiS4Vert.tex", saveDir,nBlocks));
      TCanvas *cRatioS4Horz = new TCanvas(Form("cRatioS4Horz%d",nBlocks));
      hProPiRatioS4Horz->Divide(hProS4Horz, hPiS4Horz, 1., 1., "B");
      //    hProPiRatioS4Horz->SetBinContent(39, 0);
      //    hProPiRatioS4Horz->SetBinContent(20, 0);
      hProPiRatioS4Horz->Draw("hist e");
      cRatioS4Horz->Print(Form("%s/%d_proPiS4Horz.png", saveDir,nBlocks));
      cRatioS4Horz->Print(Form("%s/%d_proPiS4Horz.pdf", saveDir,nBlocks));
      cRatioS4Horz->Print(Form("%s/%d_proPiS4Horz.tex", saveDir,nBlocks));

      hdtof1d_sub->SetLineWidth(2);
      hdtof1d->SetLineWidth(2);
      hProPiRatioS4Horz->SetLineWidth(2);
      hProPiRatioS4Vert->SetLineWidth(2);
      hProS4Horz->SetLineWidth(2);
      hProS4Vert->SetLineWidth(2);
      hPiS4Horz->SetLineWidth(2);
      hPiS4Vert->SetLineWidth(2);
    
      if (nBlocks == 0) {
	hdtof1d_sub->SetLineColor(kBlack);
	hdtof1d->SetLineColor(kBlack);
	hProPiRatioS4Horz->SetLineColor(kBlack);
	hProPiRatioS4Vert->SetLineColor(kBlack);
	hProS4Horz->SetLineColor(kBlack);
	hProS4Vert->SetLineColor(kBlack);
	hPiS4Horz->SetLineColor(kBlack);
	hPiS4Vert->SetLineColor(kBlack);
	leg->AddEntry(hdtof1d, "0 blocks", "l");     
	double intProS4Horz = hProS4Horz->Integral();
	double intPiS4Horz  = hPiS4Horz->Integral();
	legProS4Horz->AddEntry(hProS4Horz, Form("0 blocks - %d per spill", (int)intProS4Horz), "l");     
	legPiS4Horz->AddEntry(hPiS4Horz, Form("0 blocks - %d per spill", (int)intPiS4Horz), "l");
	legRatioVert->AddEntry(hProPiRatioS4Vert, "0 blocks", "l");
      }
      else if (nBlocks == 1) {
	hdtof1d_sub->SetLineColor(kRed);
	hdtof1d->SetLineColor(kRed);
	hProPiRatioS4Horz->SetLineColor(kRed);
	hProPiRatioS4Vert->SetLineColor(kRed);
	hProS4Horz->SetLineColor(kRed);
	hProS4Vert->SetLineColor(kRed);
	hPiS4Horz->SetLineColor(kRed);
	hPiS4Vert->SetLineColor(kRed);
	leg->AddEntry(hdtof1d, "1 block", "l");
	double intProS4Horz = hProS4Horz->Integral();
	double intPiS4Horz  = hPiS4Horz->Integral();
	legProS4Horz->AddEntry(hProS4Horz, Form("1 block - %d per spill", (int)intProS4Horz), "l");     
	legPiS4Horz->AddEntry(hPiS4Horz, Form("1 block - %d per spill", (int)intPiS4Horz), "l");     
	legRatioVert->AddEntry(hProPiRatioS4Vert, "1 block", "l");
      }
      else if (nBlocks == 2) {
	hdtof1d_sub->SetLineColor(kBlue);
	hdtof1d->SetLineColor(kBlue);
	hProPiRatioS4Horz->SetLineColor(kBlue);
	hProPiRatioS4Vert->SetLineColor(kBlue);
	hProS4Horz->SetLineColor(kBlue);
	hProS4Vert->SetLineColor(kBlue);
	hPiS4Horz->SetLineColor(kBlue);
	hPiS4Vert->SetLineColor(kBlue);
	leg->AddEntry(hdtof1d, "2 blocks", "l");
	double intProS4Horz = hProS4Horz->Integral();
	double intPiS4Horz  = hPiS4Horz->Integral();
	legProS4Horz->AddEntry(hProS4Horz, Form("2 blocks - %d per spill", (int)intProS4Horz), "l");     
	legPiS4Horz->AddEntry(hPiS4Horz, Form("2 blocks - %d per spill", (int)intPiS4Horz), "l");     
	legRatioVert->AddEntry(hProPiRatioS4Vert, "2 blocks", "l");
      }
      else if (nBlocks == 3) {
	hdtof1d_sub->SetLineColor(kCyan+1);
	hdtof1d->SetLineColor(kCyan+1);
	hProPiRatioS4Horz->SetLineColor(kCyan+1);
	hProPiRatioS4Vert->SetLineColor(kCyan+1);
	hProS4Horz->SetLineColor(kCyan+1);
	hProS4Vert->SetLineColor(kCyan+1);
	hPiS4Horz->SetLineColor(kCyan+1);
	hPiS4Vert->SetLineColor(kCyan+1);
	leg->AddEntry(hdtof1d, "3 blocks", "l");
	double intProS4Horz = hProS4Horz->Integral();
	double intPiS4Horz  = hPiS4Horz->Integral();
	legProS4Horz->AddEntry(hProS4Horz, Form("3 blocks - %d per spill", (int)intProS4Horz), "l");     
	legPiS4Horz->AddEntry(hPiS4Horz, Form("3 blocks - %d per spill", (int)intPiS4Horz), "l");     
	legRatioVert->AddEntry(hProPiRatioS4Vert, "3 blocks", "l");
      }
 
      hdtof1d->Scale(1. / nSpillsTrue);
      hsDtof->Add(hdtof1d);
      hsBkgSub->Add(hdtof1d_sub);
      hProS4Horz->Scale(20./6.);
      hPiS4Horz->Scale(20./6.);
      hProS4Vert->Scale(10./3.3);
      hPiS4Vert->Scale(10./3.3);
      hsProS4Vert->Add(hProS4Vert);
      hsProS4Horz->Add(hProS4Horz);
      hsPiS4Vert->Add(hPiS4Vert);
      hsPiS4Horz->Add(hPiS4Horz);
      hsRatioS4Vert->Add(hProPiRatioS4Vert);
      hsRatioS4Horz->Add(hProPiRatioS4Horz);
      hProS4Vert->Write();
      hPiS4Vert->Write();
      hProS4Horz->Write();
      hPiS4Horz->Write();
      hProPiRatioS4Vert->Write();
      hProPiRatioS4Horz->Write();
      hMSq->Write();

      legRatioVert->Write("legRatioVert");

      const char* nustof;
      if (nBlocks == 0) nustof = Form("%sData_2018_8_31_b2_800MeV_0block.root", ustofDir);
      else if (nBlocks==1) nustof = Form("%sData_2018_9_1_b4_800MeV_1block_bend4cm.root", ustofDir);
      else if (nBlocks==2) nustof = Form("%sData_2018_9_1_b2_800MeV_2block_bend4cm.root", ustofDir);
      else if (nBlocks==3) nustof = Form("%sData_2018_9_1_b3_800MeV_3block_bend4cm.root", ustofDir);
 
      cout<<nustof<<endl;
      // Read in ustof file
      TFile *finustof = new TFile(nustof, "read");
      TTree *utree = (TTree*)finustof->Get("tree");
      double tTrig;
      double tS1;
      double tSoSd;
      float xToF[50];
      float yToF[50];
      int nhit;
      utree->SetBranchAddress("tTrig", &tTrig);
      //    utree->SetBranchAddress("tS1", &tS1);
      utree->SetBranchAddress("xToF", xToF);
      utree->SetBranchAddress("yToF", yToF);
      utree->SetBranchAddress("nhit", &nhit);
      //    utree->SetBranchAddress("tSoSd", &tSoSd);

      hAllS4Horz->Add(hProS4Horz);
      hAllS4Horz->Add(hPiS4Horz);

      TH1D *hXAngleS1S2 = new TH1D(Form("hXAngleS1S2%d", nBlocks), Form("Angular distribution of hits in S3 (S1 & S2 triggers), %d blocks; #theta / degrees; Events / spill", nBlocks), 100, -3.8, 6.2);
      for (int t=0; t<utree->GetEntries(); t++) {
	utree->GetEntry(t);
	// Has an S1 and S2 hit
	if (tTrig != 0 ) {
	  for (int n=0; n<nhit; n++) {
	    double positionX = ((xToF[n] - 4.) / 152.)*(s3EndX - s3StartX) + s3StartX;
	    double positionY = ((xToF[n] - 4.) / 152.)*(s3s1EndY - s3s1StartY)+s3s1StartY;
	    double angleOffAxis = TMath::ATan(positionX / positionY) * (180. / TMath::Pi());
	    hXAngleS1S2->Fill(angleOffAxis);
	  } // for (int n=0; n<nhit; n++) 
	} // if (tTrig !=0 ) 
      } // for (int t=0; t<utree->GetEntries(); t++)

      // Integrate this hist between the angular limits of S2 then scale by this number
      const double s2ThetaLow = -0.359;
      const double s2ThetaHi  = 3.957;
      const double s4ThetaLow = 0.401;
      const double s4ThetaHi  = 6.083;
      double s2Int = hXAngleS1S2->Integral(hXAngleS1S2->GetXaxis()->FindBin(s2ThetaLow), hXAngleS1S2->GetXaxis()->FindBin(s2ThetaHi));
      cout<<s2Int<<endl;
      hXAngleS1S2->Scale(1. / s2Int);
      cout<<"Integral "<<hXAngleS1S2->Integral(hXAngleS1S2->GetXaxis()->FindBin(s2ThetaLow), hXAngleS1S2->GetXaxis()->FindBin(s2ThetaHi))<<endl;
      double s4Int = hXAngleS1S2->Integral(hXAngleS1S2->GetXaxis()->FindBin(s4ThetaLow), hXAngleS1S2->GetXaxis()->FindBin(s4ThetaHi));
      cout<<"S4 correction factor "<<s4Int<<endl;

      hAllS4Horz->Scale(1. / nSpillsTrue);
    } // Non 4 block data
    // Loop over all the 4 block data
    else {
      // Number of signal particles using just cut and count
      double nP  = 0.;
      double nPi = 0.;
      // Define signal and background functions to be fitted
      // Signals are gaussians
      TF1 *sPro = new TF1("sPro", "gaus", 66, 100);
      TF1 *sPi  = new TF1("sPi", "gaus", 35, 55);
      // Exponential background
      TF1 *fBkgExp = new TF1("fBkgExp","expo", 30, 160);
      sPro->SetLineColor(kGreen+2);
      sPi->SetLineColor(kRed);
      TF1 *fSplusBExp = new TF1("signal+bkg exp", "gaus(0)+gaus(3)+expo(6)", 30, 160);
      fSplusBExp->SetParNames("const 1", "mean 1", "sigma 1",
			      "const 2", "mean 2", "sigma 2",
			      "bkgconst", "bkgdecay");
      fSplusBExp->SetLineColor(kBlack);
      // For spill counting normalisation
      int nSpills = 0;
      int nSpillsTrue = 0;
      double lastSpill = 0.;
      // 1D ToF for background subtraction
      TH1D *hdtof1d = new TH1D(Form("hdtof1d_%d",nBlocks), Form("Time of flight, %d blocks; S4 - S1 / ns; Events", nBlocks), 260, 30, 160);
      // Cosmics hists
      // TH2D *h2Cosmics = new TH2D(Form("h2Cosmics%d",nBlocks), Form("Cosmic flux, %d blocks; x / cm; Bar; Hz",nBlocks), 20, 0, 140, 10, 0.5, 10.5);
      // TH1D *hCosmicsVert = new TH1D(Form("hCosmicsVert%d",nBlocks), Form("Cosmic flux, %d blocks; x / cm; Hz",nBlocks), 10, 0.5, 10.5);
      // hCosmicsVert->Sumw2();
      // TH1D *hCosmicsHorz = new TH1D(Form("hCosmicsHorz%d",nBlocks), Form("Cosmic flux, %d blocks; x / cm; Hz",nBlocks), 20, 0, 140);
      // hCosmicsHorz->Sumw2();
      TH1D *hEffTotal = new TH1D("hEffTotal", "Eff 4", 10, 0.5, 10.5);
      for (int b4=0; b4<str4BlockVec.size(); b4++) {
	cout<<"Analysing sample "<<b4<<" of "<<str4BlockVec.size()<<endl;
	int nSpillsTrueTmp = 0;
	TH1D *hDtofTmp   = new TH1D(Form("hsDtof_%s",str4BlockVec[b4]), "S4 ToF spectrum; S4 - S2 / ns; Events / spill", 130, 30., 160.);
	TH1D *hProS4VertTmp = new TH1D(Form("hProS4Vert_%s",str4BlockVec[b4]), "S1 #cap S2 #cap S4 angular distribution of proton hits; #phi / degrees; Events / spill", 10, -1.5, 1.8);
	TH1D *hPiS4VertTmp  = new TH1D(Form("hPiS4Vert_%s",str4BlockVec[b4]), "S1 #cap S2 #cap S4 angular distribution of MIP hits; #phi / degrees; Events / spill", 10, -1.5, 1.8);
	TH1D *hProS4HorzTmp = new TH1D(Form("hProS4Horz_%s",str4BlockVec[b4]), "S1 #cap S2 #cap S4 angular distribution of proton hits; #theta / degrees; Events / spill", 20, 0., 6.);
	TH1D *hPiS4HorzTmp  = new TH1D(Form("hPiS4Horz_%s",str4BlockVec[b4]), "S1 #cap S2 #cap S4 angular distribution of MIP hits; #theta / degrees; Events / spill", 20, 0., 6.);
	TH1D *hRatioS4VertTmp = new TH1D(Form("hRatioS4Vert_%s",str4BlockVec[b4]), "S1 #cap S2 #cap S4 angular distribution of proton/MIP ratio; #phi / degrees; Protons/MIPs", 10, -1.5, 1.8);
	TH1D *hRatioS4HorzTmp = new TH1D(Form("hRatioS4Horz_%s",str4BlockVec[b4]), "S1 #cap S2 #cap S4 angular distribution of proton/MIP ratio; #theta / degrees; Protons/MIPs", 20, 0., 6.);
	// Find the correct dstof files
	Int_t runMin=-1;
	Int_t runMax=-1;
	double startTime = 0;
	double endTime   = 0;
	// Get start and end times from the relevant utof files
	TFile *futofTmp = new TFile(Form("%s/%s",ustofDir, str4BlockVec.at(b4)), "read");
	TTree *treeTmp = (TTree*)futofTmp->Get("tree");
	double tS1Tmp;
	treeTmp->SetBranchAddress("tS1", &tS1Tmp);
	TNamed *start = 0;
	TNamed *end   = 0;
	futofTmp->GetObject("start_of_run", start);
	const char* startchar = start->GetTitle();
	std::string startstr(startchar);
	std::string unixstart = startstr.substr(25,10);
	startTime = stoi(unixstart);
	treeTmp->GetEntry(treeTmp->GetEntries() - 1);
	endTime = startTime + (tS1Tmp/1e9);
	futofTmp->Close();
	delete futofTmp;

	for (int irun=950; irun<1400; irun++) {
	  TFile *fin = new TFile(Form("%srun%d/DsTOFcoincidenceRun%d_tdc1.root", dstofDir, irun, irun), "read");
	  RawDsTofCoincidence *tofCoinTemp = NULL;
	  TTree *tree = (TTree*) fin->Get("tofCoinTree");
	  tree->SetBranchAddress("tofCoin", &tofCoinTemp);
	  tree->GetEntry(0);
	  UInt_t firstTemp = tofCoinTemp->unixTime[0];
	  tree->GetEntry(tree->GetEntries()-1);
	  UInt_t lastTemp = tofCoinTemp->unixTime[0];
      
	  fin->Close();
	  delete fin;
      
	  if (firstTemp>endTime){
	    break;
	  }
      
	  if (firstTemp<startTime && lastTemp>startTime){
	    runMin = irun;
	  }
      
	  if (firstTemp<endTime && lastTemp>endTime){
	    runMax = irun;
	  }   
	} // for (int irun=950; irun<1400; irun++) 
    
	cout << "Min and max runs are " << runMin << " " << runMax << endl;

	// Need these to calculate bar efficiencies
	TH1D *hCoins = new TH1D(Form("hCoins_%d_%s",nBlocks, str4BlockVec[b4]), Form("Bar coincidences + S_{1,2} coincidences, %d blocks; Bar; Events",nBlocks), 10, 0.5, 10.5);
	TH1D *hHits  = new TH1D(Form("hHits_%d_%s",nBlocks, str4BlockVec[b4]), Form("PMT hits + S_{1,2} coincidences, %d blocks; Bar; Events",nBlocks), 10, 0.5, 10.5);
	TH1D *hEff   = new TH1D(Form("hEff_%d_%s",nBlocks, str4BlockVec[b4]), Form("S4 bar efficiencies %d blocks; Bar; Events",nBlocks), 10, 0.5, 10.5);
	hHits->Sumw2();
	hCoins->Sumw2();
	cout<<"Calculating efficiencies"<<endl;
	// In this loop calculate the bar-by-bar efficiencies
	for (int itdc=0; itdc<2; itdc++) {
	  // Need both the coincidence and the raw trees for this
	  double tempUstof;
	  double ustofNs;
	  for (int irun=runMin; irun<runMax+1; irun++) {
	    // Load input files
	    TFile *tofCoinFile = new TFile(Form("%srun%d/DsTOFcoincidenceRun%d_tdc%d.root", dstofDir, irun, irun, itdc+1));
	    TFile *tofFile     = new TFile(Form("%srun%d/DsTOFtreeRun%d_tdc%d.root", dstofDir, irun, irun, itdc+1));
	    TTree *tofCoinTree = (TTree*)tofCoinFile->Get("tofCoinTree");
	    TTree *tofTree = (TTree*)tofFile->Get("tofTree");
	    RawDsTofCoincidence *tofCoin = NULL;
	    RawDsTofHeader *tof = NULL;
	    tofCoinTree->SetBranchAddress("tofCoin", &tofCoin);
	    tofTree->SetBranchAddress("tof", &tof);
	    // Create new TTree for ustof signals
	    TTree *ustofTree = new TTree("ustofTree", "ustof");
	    ustofTree->SetDirectory(0);
	    ustofTree->Branch("ustofNs", &ustofNs, "ustofNs/D");
	    tempUstof=0;
	    // Build tree of all the ustof hit times
	    for (int i=0; i<tofTree->GetEntries(); i++) {
	      tofTree->GetEntry(i);
	      if (tof->unixTime < startTime) continue;
	      if (tof->unixTime > endTime) break;

	      if (tof->channel == 13) {
		ustofNs = tof->fakeTimeNs - ustofDelay;
		// ustof signals shouldn't be coming closer than 500ns
		if ( (ustofNs - tempUstof) < 500.) continue;
		ustofTree->Fill();
		tempUstof = ustofNs;
	      } // if (tof->channel == 13) 
	    } // for (int i=0; i<tof->GetEntries(); i++)
	    ustofTree->BuildIndex("ustofNs");
	    // Now loop over coincidence file to find number of coincidences for each bar
	    for (int t=0; t<tofCoinTree->GetEntries(); t++) {
	      tofCoinTree->GetEntry(t);
	      if (tofCoin->unixTime[0] < startTime) continue;
	      if (tofCoin->unixTime[0] > endTime) break;
	      double deltat = TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1]);
	      double dstofHitT = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (10. - TMath::Abs(deltat) / 2. );
	      double tofCalc = dstofHitT - tofCoin->usTofSignal;
	      if (tofCalc > 70. && tofCalc < 200. && tofCoin->bar != 10) {
		hCoins->Fill(tofCoin->bar);
	      } // if (tofCalc > 70. && tofCalc < 200.)
	    } // for (int t=0; t<tofCoinTree->GetEntries(); t++) 
	    // Loop over hit file to find number of hits in coincidence 
	    for (int t=0; t<tofTree->GetEntries(); t++) {
	      tofTree->GetEntry(t);
	      if (tof->unixTime < startTime) continue;
	      if (tof->unixTime > endTime) break;
	      if (tof->channel > 10 || tof->channel == 0) continue;
	      int ustofEntry = ustofTree->GetEntryNumberWithBestIndex(tof->fakeTimeNs);
	      ustofTree->GetEntry(ustofEntry);

	      if (tof->fakeTimeNs - ustofNs < 260. && tof->fakeTimeNs - ustofNs > 130.) {
		if (itdc == 0) { // TDC1
		  if (tof->channel % 2 == 1) { // PMT B
		    hHits->Fill(((tof->channel + 1)/ -2) + 11);
		  } // if (tof->channel % 2 == 1)
		  else { // PMT A
		    hHits->Fill((tof->channel / -2) + 11);
		  }
		} // TDC 1
		else { // TDC 2
		  if (tof->channel % 2 == 1) { // PMT B
		    hHits->Fill(((tof->channel + 1)/ -2) + 6);
		  } // if (tof->channel % 2 == 1)
		  else { // PMT A
		    hHits->Fill((tof->channel / -2) + 6);
		  }
		} // TDC 2
	      } // if (tof->fakeTimeNs - ustofNs < 260. && tof->fakeTimeNs - ustofNs > 130.)
	    } // for (int t=0; t<tofTree->GetEntries(); t++)
	    delete tof;
	    delete tofCoin;
	    delete tofTree;
	    delete tofCoinTree;
	    delete ustofTree;
	    tofFile->Close();
	    tofCoinFile->Close();
	  } // for (int irun=runMin; irun<runMax+1; irun++) 
	} // for (int itdc=0; itdc<2; itdc++)
	fout->cd();
	TCanvas *cHits = new TCanvas(Form("%d_cHits",nBlocks));
	hHits->Draw("hist e");
	hHits->Write();
	cHits->Print(Form("%s/%d_barHits.png",saveDir,nBlocks));
	cHits->Print(Form("%s/%d_barHits.pdf",saveDir,nBlocks));
	TCanvas *cCoins = new TCanvas(Form("%d_cCoins",nBlocks));
	hCoins->Draw("hist e");
	hCoins->Write();
	cCoins->Print(Form("%s/%d_barCoins.png",saveDir,nBlocks));
	cCoins->Print(Form("%s/%d_barCoins.pdf",saveDir,nBlocks));
	TCanvas *cEff = new TCanvas(Form("%d_cEff",nBlocks));
	hHits->Add(hCoins, -1.);
	hEff->Divide(hCoins, hHits, 1., 1., "B");
	hEff->Draw("hist e");
	hEff->Write();
	cEff->Print(Form("%s/%d_barEff.png",saveDir,nBlocks));
	cEff->Print(Form("%s/%d_barEff.pdf",saveDir,nBlocks));
	cEff->Print(Form("%s/%d_barEff.tex",saveDir,nBlocks));
	cout<<"Getting the signal hits"<<endl;
	// Now loop over the coincidence files again and calculate the angular distributions
	for (int itdc=0; itdc<2; itdc++) {
	  TChain *tofCoinChain = new TChain("tofCoinTree");
	  for (int irun=runMin; irun<runMax+1; irun++){
	    tofCoinChain->Add(Form("%srun%d/DsTOFcoincidenceRun%d_tdc%d.root", dstofDir, irun, irun, itdc+1));
	  } // for (int irun=runMin; irun<runMax+1; irun++)
	  RawDsTofCoincidence *tofCoin = NULL;
	  tofCoinChain->SetBranchAddress("tofCoin", &tofCoin);
	  for (int h=0; h<tofCoinChain->GetEntries(); h++) {
	    tofCoinChain->GetEntry(h);
	    if (tofCoin->unixTime[0]<startTime) continue;
	    if (tofCoin->unixTime[0]>endTime) break;

	    if (h % 100000 == 0) cout<<"Entry "<<h<<" of "<<tofCoinChain->GetEntries()<<endl;

	    if (tofCoin->lastDelayedBeamSignal != lastSpill && itdc == 0) {
	      lastSpill = tofCoin->lastDelayedBeamSignal;
	      nSpills++;
	      int nTrueHits = 0;
	      for (int sp=h; sp<tofCoinChain->GetEntries(); sp++) {
		tofCoinChain->GetEntry(sp);
		if (tofCoin->unixTime[0]<startTime) continue;
		if (tofCoin->unixTime[0]>endTime) break;
		if (tofCoin->fakeTimeNs[0] > lastSpill + 3e9) break;

		if ((tofCoin->fakeTimeNs[0] - tofCoin->usTofSignal) < 200. &&
		    (tofCoin->fakeTimeNs[0] - tofCoin->usTofSignal) > 70. &&
		    (tofCoin->fakeTimeNs[0] - lastSpill) < 1e9 &&
		    (tofCoin->fakeTimeNs[0] - lastSpill) > 0) {
		  nTrueHits++;
		  if (nTrueHits > 25) {
		    nSpillsTrue++;
		    nSpillsTrueTmp++;
		    break;
		  }
		}
	      } // for (int sp=h; sp<tofCoinChain->GetEntries(); sp++) 
	    } // if (tofCoin->lastDelayedBeamSignal != lastSpill && itdc == 0)
	    // Need to calculate total signal hits here
	    // Weight these by the efficiency of calculated above
	    double deltat = TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1]  );
	    double dstofHitT = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (10. - TMath::Abs(deltat) / 2. );
	    double tofCalc = dstofHitT - tofCoin->usTofSignal - dstofShift;
	    if (tofCalc < 160. && tofCalc > 30. && tofCoin->bar != 10) {
	      hdtof1d->Fill(tofCalc, 1. / hEff->GetBinContent(tofCoin->bar));
	      hMSq->Fill(massFromTime(tofCalc, 0.8, s2s4Dist), 1. / hEff->GetBinContent(tofCoin->bar));
	      double positionXP = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 70.));
	      if (tofCalc < piHi & tofCalc > piLow) { 
		nPi += (1. / /*h2CosEff->GetBinContent( h2CosEff->GetXaxis()->FindBin(positionXP), tofCoin->bar)*/hEff->GetBinContent(tofCoin->bar));
		// Calculate position of hit in global coordinates
		double positionX = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 65.) / 130.) * (s4OffAxisEndX - s4OffAxisStartX) + s4OffAxisStartX;
		double positionY = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 65.) / 130.) * (baselineS1S4End - baselineS1S4Start) + baselineS1S4Start;
		double positionZ = (tofCoin->bar * 0.075) - 0.375 - 0.01;
		// Calculate the angles relative to the nominal beamline
		double angleTheta = TMath::ATan(positionX / positionY) * (180./TMath::Pi());
		double anglePhi   = TMath::ATan(positionZ / positionY) * (180./TMath::Pi());
		hPiS4Horz->Fill(angleTheta, 1. / /*h2CosEff->GetBinContent( h2CosEff->GetXaxis()->FindBin(positionXP), tofCoin->bar)*/ hEff->GetBinContent(tofCoin->bar));
		hPiS4Vert->Fill(anglePhi,   1. / /*h2CosEff->GetBinContent( h2CosEff->GetXaxis()->FindBin(positionXP), tofCoin->bar)*/hEff->GetBinContent(tofCoin->bar));
		hPiS4HorzTmp->Fill(angleTheta, 1. / hEff->GetBinContent(tofCoin->bar));
		hPiS4VertTmp->Fill(anglePhi,   1. / hEff->GetBinContent(tofCoin->bar));
	      } // if (tofCalc < piHi & tofCalc > piLow)
	      else if (tofCalc < proHi & tofCalc > proLow) {
		nP += (1. / /*h2CosEff->GetBinContent( h2CosEff->GetXaxis()->FindBin(positionXP), tofCoin->bar)*/hEff->GetBinContent(tofCoin->bar));
		// Calculate position of hit in global coordinates
		double positionX = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 65.) / 130.) * (s4OffAxisEndX - s4OffAxisStartX) + s4OffAxisStartX;
		double positionY = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 65.) / 130.) * (baselineS1S4End - baselineS1S4Start) + baselineS1S4Start;
		double positionZ = (tofCoin->bar * 0.075) - 0.375 - 0.01;
		// Calculate the angles relative to the nominal beamline
		double angleTheta = TMath::ATan(positionX / positionY) * (180./TMath::Pi());
		double anglePhi   = TMath::ATan(positionZ / positionY) * (180./TMath::Pi());
		hProS4Horz->Fill(angleTheta, 1. / /*h2CosEff->GetBinContent( h2CosEff->GetXaxis()->FindBin(positionXP), tofCoin->bar)*/hEff->GetBinContent(tofCoin->bar));
		hProS4Vert->Fill(anglePhi,   1. / /*h2CosEff->GetBinContent( h2CosEff->GetXaxis()->FindBin(positionXP), tofCoin->bar)*/hEff->GetBinContent(tofCoin->bar));
		hProS4HorzTmp->Fill(angleTheta, 1. / hEff->GetBinContent(tofCoin->bar));
		hProS4VertTmp->Fill(anglePhi,   1. / hEff->GetBinContent(tofCoin->bar));
	      } // else if (tofCalc < proHi & tofCalc > proLow) 
	    } // if (tofCalc < 160. && tofCalc > 30.) 
	  } // for (int h=0; h<tofCoinChain->GetEntries(); h++) 
	  delete tofCoin;
	  delete tofCoinChain;
	} // for (int itdc=0; itdc<2; itdc++)
	fout->cd();
	hProS4HorzTmp->Scale(1./nSpillsTrueTmp);
	hPiS4HorzTmp->Scale(1./nSpillsTrueTmp);
	hProS4VertTmp->Scale(1./nSpillsTrueTmp);
	hPiS4VertTmp->Scale(1./nSpillsTrueTmp);
	hRatioS4HorzTmp->SetLineWidth(2);
	hRatioS4HorzTmp->Divide(hProS4HorzTmp, hPiS4HorzTmp, 1., 1., "B");
	hRatioS4VertTmp->Divide(hProS4VertTmp, hPiS4VertTmp, 1., 1., "B");
	hProS4HorzTmp->Write();
	hPiS4HorzTmp->Write();
	hProS4VertTmp->Write();
	hPiS4VertTmp->Write();
	hRatioS4HorzTmp->Write();
	hRatioS4VertTmp->Write();
	cout<<"Spills "<<nSpillsTrueTmp<<" in this sample"<<endl;
	hEff->Scale(nSpillsTrueTmp);
	hEffTotal->Add(hEff);
      } // Loop over the four blocks     
      fout->cd();
      TCanvas *c2d_exp = new TCanvas(Form("%d_c2d_exp",nBlocks));
      c2d_exp->SetLogy();
      hdtof1d->Fit(sPi, "R");
      hdtof1d->Fit(sPro, "R");
      hdtof1d->Fit(fBkgExp, "R");
      Double_t parExp[8];
      sPro->GetParameters(&parExp[0]);
      sPi->GetParameters(&parExp[3]);
      fBkgExp->GetParameters(&parExp[6]);
      fSplusBExp->SetParameters(parExp);
      hdtof1d->Fit(fSplusBExp, "R");
      hdtof1d->Draw("hist");
      fSplusBExp->Draw("same");
      hdtof1d->Write();
      c2d_exp->Print(Form("%s/%d_dtof1d.png",saveDir,nBlocks));
      c2d_exp->Print(Form("%s/%d_dtof1d.pdf",saveDir,nBlocks));
      c2d_exp->Print(Form("%s/%d_dtof1d.tex",saveDir,nBlocks));

      // Now we have the fit values, loop over again and subtract the background
      // For each bin, find the fraction of each particle type which are background and 
      // this fraction from the bin
      TF1 *fSub = new TF1("fSub", "exp([0]+[1]*x)", 30, 160);
      fSub->SetParameter(0, fSplusBExp->GetParameter("bkgconst"));  
      fSub->SetParameter(1, fSplusBExp->GetParameter("bkgdecay"));
      // Background subtracted tof spectrum
      TH1D *hdtof1d_sub = (TH1D*)hdtof1d->Clone(Form("hdtof1d_sub_%d",nBlocks));
      hdtof1d_sub->Add(fSub, -1.);
      hdtof1d_sub->Scale(1. / nSpillsTrue);
      // If the bin content drops below 0, set to 0
      for (int b = 0; b <=  hdtof1d_sub->GetNbinsX(); b++) {
	if (hdtof1d_sub->GetBinContent(b) < 0.) {
	  hdtof1d_sub->SetBinContent(b, 0);
	} // if (hdtof1d_sub->GetBinContent(b) < 0.)
      } // for (int b = 0; hdtof1d_sub->GetNbinsX(); b++)
      TCanvas *cdtof_sub = new TCanvas(Form("%d_cdtof_sub",nBlocks));

      // Integrate background function between proton and pion windows and then subtract
      double piBkg  = fSub->Integral(piLow,  piHi) / (0.5*nSpillsTrue);
      double proBkg = fSub->Integral(proLow, proHi) / (0.5*nSpillsTrue);
      cout<<"Pion backgrounds "<<piBkg<<", proton backgrounds "<<proBkg<<endl;

      cdtof_sub->SetLogy();
      hdtof1d_sub->Draw("hist");
      hdtof1d_sub->GetXaxis()->SetLabelSize(0.05);
      hdtof1d_sub->GetYaxis()->SetLabelSize(0.05);
      hdtof1d_sub->GetXaxis()->SetTitleSize(0.05);
      hdtof1d_sub->GetYaxis()->SetTitleSize(0.05);
      hdtof1d_sub->Write();
      cdtof_sub->Print(Form("%s/%d_dtof1d_bkgSub.png", saveDir, nBlocks));
      cdtof_sub->Print(Form("%s/%d_dtof1d_bkgSub.pdf", saveDir, nBlocks));

      cout<<"Spills "<<nSpills<<" ("<<nSpillsTrue<<" true)"<<endl;

      // Subtract background hits
      TCanvas *cProS4Horz = new TCanvas(Form("cProS4Horz%d",nBlocks));
      hProS4Horz->Scale(1. / nSpillsTrue);
      for (int b=1; b < 21; b++) {
	hProS4Horz->SetBinContent(b, hProS4Horz->GetBinContent(b) - (proBkg/20.));
	if (hProS4Horz->GetBinContent(b) < 0.) {
	  hProS4Horz->SetBinContent(b, 0);
	  hProS4Horz->SetBinError(b, 0);
	}
      }
      hProS4Horz->Draw("hist e");
      cProS4Horz->Print(Form("%s/%d_proS4Horz.png",saveDir,nBlocks));
      cProS4Horz->Print(Form("%s/%d_proS4Horz.pdf",saveDir,nBlocks));
      cProS4Horz->Print(Form("%s/%d_proS4Horz.tex",saveDir,nBlocks));
      TCanvas *cPiS4Horz  = new TCanvas(Form("cPiS4Horz%d",nBlocks));
      hPiS4Horz->Scale(1. / nSpillsTrue);
      for (int b=1; b < 21; b++) {
	hPiS4Horz->SetBinContent(b, hPiS4Horz->GetBinContent(b) - (piBkg/20.));
	if (hPiS4Horz->GetBinContent(b) < 0.) {
	  hPiS4Horz->SetBinError(b, 0);
	  hPiS4Horz->SetBinContent(b, 0);
	}
      }
      hPiS4Horz->Draw("hist e");
      cPiS4Horz->Print(Form("%s/%d_piS4Horz.png",saveDir,nBlocks));
      cPiS4Horz->Print(Form("%s/%d_piS4Horz.pdf",saveDir,nBlocks));
      cPiS4Horz->Print(Form("%s/%d_piS4Horz.tex",saveDir,nBlocks));
      TCanvas *cProS4Vert = new TCanvas(Form("cProS4Vert%d",nBlocks));
      hProS4Vert->Scale(1. / nSpillsTrue);
      for (int b=1; b < 10; b++) {
	hProS4Vert->SetBinContent(b, hProS4Vert->GetBinContent(b) - (proBkg/9.));
	if (hProS4Vert->GetBinContent(b) < 0.) {
	  hProS4Vert->SetBinError(b, 0);
	  hProS4Vert->SetBinContent(b, 0);
	}
      }
      hProS4Vert->Draw("hist e");
      cProS4Vert->Print(Form("%s/%d_proS4Vert.png",saveDir,nBlocks));
      cProS4Vert->Print(Form("%s/%d_proS4Vert.pdf",saveDir,nBlocks));
      cProS4Vert->Print(Form("%s/%d_proS4Vert.tex",saveDir,nBlocks));
      TCanvas *cPiS4Vert  = new TCanvas(Form("cPiS4Vert%d",nBlocks));
      hPiS4Vert->Scale(1. / nSpillsTrue);
      for (int b=1; b < 10; b++) {
	hPiS4Vert->SetBinContent(b, hPiS4Vert->GetBinContent(b) - (piBkg/9.));
	if (hPiS4Vert->GetBinContent(b) < 0.) {
	  hPiS4Vert->SetBinError(b, 0);
	  hPiS4Vert->SetBinContent(b, 0);
	}
      }
      hPiS4Vert->Draw("hist e");
      cPiS4Vert->Print(Form("%s/%d_piS4Vert.png",saveDir,nBlocks));
      cPiS4Vert->Print(Form("%s/%d_piS4Vert.pdf",saveDir,nBlocks));
      cPiS4Vert->Print(Form("%s/%d_piS4Vert.tex",saveDir,nBlocks));

      TCanvas *cRatioS4Vert = new TCanvas(Form("cRatioS4Vert%d",nBlocks));
      hProPiRatioS4Vert->Divide(hProS4Vert, hPiS4Vert, 1., 1., "B");
      hProPiRatioS4Vert->Draw("hist e");
      cRatioS4Vert->Print(Form("%s/%d_proPiS4Vert.png", saveDir,nBlocks));
      cRatioS4Vert->Print(Form("%s/%d_proPiS4Vert.pdf", saveDir,nBlocks));
      cRatioS4Vert->Print(Form("%s/%d_proPiS4Vert.tex", saveDir,nBlocks));
      TCanvas *cRatioS4Horz = new TCanvas(Form("cRatioS4Horz%d",nBlocks));
      hProPiRatioS4Horz->Divide(hProS4Horz, hPiS4Horz, 1., 1., "B");
      hProPiRatioS4Horz->Draw("hist e");
      cRatioS4Horz->Print(Form("%s/%d_proPiS4Horz.png", saveDir,nBlocks));
      cRatioS4Horz->Print(Form("%s/%d_proPiS4Horz.pdf", saveDir,nBlocks));
      cRatioS4Horz->Print(Form("%s/%d_proPiS4Horz.tex", saveDir,nBlocks));

      hdtof1d_sub->SetLineWidth(2);
      hdtof1d->SetLineWidth(2);
      hProPiRatioS4Horz->SetLineWidth(2);
      hProPiRatioS4Vert->SetLineWidth(2);
      hProS4Horz->SetLineWidth(2);
      hProS4Vert->SetLineWidth(2);
      hPiS4Horz->SetLineWidth(2);
      hPiS4Vert->SetLineWidth(2);
    
      hdtof1d_sub->SetLineColor(kOrange+1);
      hdtof1d->SetLineColor(kOrange+1);
      hProPiRatioS4Horz->SetLineColor(kOrange+1);
      hProPiRatioS4Vert->SetLineColor(kOrange+1);
      hProS4Horz->SetLineColor(kOrange+1);
      hProS4Vert->SetLineColor(kOrange+1);
      hPiS4Horz->SetLineColor(kOrange+1);
      hPiS4Vert->SetLineColor(kOrange+1);
      leg->AddEntry(hdtof1d, "4 blocks", "l");
      double intProS4Horz = hProS4Horz->Integral();
      double intPiS4Horz  = hPiS4Horz->Integral();
      legProS4Horz->AddEntry(hProS4Horz, Form("4 blocks - %d per spill", (int)intProS4Horz), "l");     
      legPiS4Horz->AddEntry(hPiS4Horz, Form("4 blocks - %d per spill", (int)intPiS4Horz), "l");     
      legRatioVert->AddEntry(hProPiRatioS4Vert, "4 blocks", "l");
    
      hdtof1d->Scale(1. / nSpillsTrue);
      hsDtof->Add(hdtof1d);
      hsBkgSub->Add(hdtof1d_sub);
      hProS4Horz->Scale(20./6.);
      hPiS4Horz->Scale(20./6.);
      hProS4Vert->Scale(10./3.3);
      hPiS4Vert->Scale(10./3.3);
      hsProS4Vert->Add(hProS4Vert);
      hsProS4Horz->Add(hProS4Horz);
      hsPiS4Vert->Add(hPiS4Vert);
      hsPiS4Horz->Add(hPiS4Horz);
      hsRatioS4Vert->Add(hProPiRatioS4Vert);
      hsRatioS4Horz->Add(hProPiRatioS4Horz);
      hProS4Vert->Write();
      hPiS4Vert->Write();
      hProS4Horz->Write();
      hPiS4Horz->Write();
      hProPiRatioS4Vert->Write();
      hProPiRatioS4Horz->Write();
      hEffTotal->SetLineColor(kOrange+1);
      hEffTotal->Scale(1. / nSpillsTrue);
      hEffTotal->Write();
      hMSq->Write();
      legRatioVert->Write("legRatioVert");

      const char* nustof;
      if (nBlocks == 0) nustof = Form("%sData_2018_8_31_b2_800MeV_0block.root", ustofDir);
      else if (nBlocks==1) nustof = Form("%sData_2018_9_1_b4_800MeV_1block_bend4cm.root", ustofDir);
      else if (nBlocks==2) nustof = Form("%sData_2018_9_1_b2_800MeV_2block_bend4cm.root", ustofDir);
      else if (nBlocks==3) nustof = Form("%sData_2018_9_1_b3_800MeV_3block_bend4cm.root", ustofDir);
      else if (nBlocks==4) nustof = Form("%sData_2018_9_1_b8_800MeV_4block_bend4cm.root", ustofDir);
    

      cout<<nustof<<endl;
      // Read in ustof file
      TFile *finustof = new TFile(nustof, "read");
      TTree *utree = (TTree*)finustof->Get("tree");
      double tTrig;
      double tS1;
      double tSoSd;
      float xToF[50];
      float yToF[50];
      int nhit;
      utree->SetBranchAddress("tTrig", &tTrig);
      //    utree->SetBranchAddress("tS1", &tS1);
      utree->SetBranchAddress("xToF", xToF);
      utree->SetBranchAddress("yToF", yToF);
      utree->SetBranchAddress("nhit", &nhit);
      //    utree->SetBranchAddress("tSoSd", &tSoSd);

      hAllS4Horz->Add(hProS4Horz);
      hAllS4Horz->Add(hPiS4Horz);

      TH1D *hXAngleS1S2 = new TH1D(Form("hXAngleS1S2%d", nBlocks), Form("Angular distribution of hits in S3 (S1 & S2 triggers), %d blocks; #theta / degrees; Events / spill", nBlocks), 100, -3.8, 6.2);
      for (int t=0; t<utree->GetEntries(); t++) {
	utree->GetEntry(t);
	// Has an S1 and S2 hit
	if (tTrig != 0 ) {
	  for (int n=0; n<nhit; n++) {
	    double positionX = ((xToF[n] - 4.) / 152.)*(s3EndX - s3StartX) + s3StartX;
	    double positionY = ((xToF[n] - 4.) / 152.)*(s3s1EndY - s3s1StartY)+s3s1StartY;
	    double angleOffAxis = TMath::ATan(positionX / positionY) * (180. / TMath::Pi());
	    hXAngleS1S2->Fill(angleOffAxis);
	  } // for (int n=0; n<nhit; n++) 
	} // if (tTrig !=0 ) 
      } // for (int t=0; t<utree->GetEntries(); t++)

      // Integrate this hist between the angular limits of S2 then scale by this number
      const double s2ThetaLow = -0.359;
      const double s2ThetaHi  = 3.957;
      const double s4ThetaLow = 0.401;
      const double s4ThetaHi  = 6.083;
      double s2Int = hXAngleS1S2->Integral(hXAngleS1S2->GetXaxis()->FindBin(s2ThetaLow), hXAngleS1S2->GetXaxis()->FindBin(s2ThetaHi));
      cout<<s2Int<<endl;
      hXAngleS1S2->Scale(1. / s2Int);
      cout<<"Integral "<<hXAngleS1S2->Integral(hXAngleS1S2->GetXaxis()->FindBin(s2ThetaLow), hXAngleS1S2->GetXaxis()->FindBin(s2ThetaHi))<<endl;
      double s4Int = hXAngleS1S2->Integral(hXAngleS1S2->GetXaxis()->FindBin(s4ThetaLow), hXAngleS1S2->GetXaxis()->FindBin(s4ThetaHi));
      cout<<"S4 correction factor "<<s4Int<<endl;

      hAllS4Horz->Scale(1. / nSpillsTrue);
    } // Loop over 4 block data
  } // for (int nBlocks = 0; nBlocks < 4; nBlocks++) 
  fout->cd();
  TCanvas *csDtof = new TCanvas("csDstof");
  csDtof->SetLogy();
  hsDtof->Draw("hist nostack");
  hsDtof->GetYaxis()->SetRangeUser(1e-2, 2e2);
  csDtof->Update();
  leg->Draw();
  hsDtof->Write();
  csDtof->Print(Form("%s/s4ToF.png", saveDir));
  csDtof->Print(Form("%s/s4ToF.pdf", saveDir));
  csDtof->Print(Form("%s/s4ToF.tex", saveDir));
  TCanvas *csBkgSub = new TCanvas("csBkgSub");
  csBkgSub->SetLogy();
  hsBkgSub->Draw("hist nostack");
  hsBkgSub->GetYaxis()->SetRangeUser(1e-2, 2e2);
  hsBkgSub->GetXaxis()->SetLabelSize(0.06);
  hsBkgSub->GetYaxis()->SetLabelSize(0.06);
  hsBkgSub->GetXaxis()->SetTitleSize(0.06);
  hsBkgSub->GetYaxis()->SetTitleSize(0.06);
  csBkgSub->SetLeftMargin(0.13);
  csBkgSub->SetBottomMargin(0.13);
  csBkgSub->SetGridx();
  csBkgSub->SetGridy();
  csBkgSub->Update();
  leg->Draw();
  hsBkgSub->Write();
  csBkgSub->Print(Form("%s/s4BkgSub.png", saveDir));
  csBkgSub->Print(Form("%s/s4BkgSub.pdf", saveDir));
  csBkgSub->Print(Form("%s/s4BkgSub.tex", saveDir));

  TCanvas *cspros4vert = new TCanvas("cspros4vert");
  hsProS4Vert->Draw("hist e nostack");
  hsProS4Vert->GetXaxis()->SetLabelSize(0.06);
  hsProS4Vert->GetYaxis()->SetLabelSize(0.06);
  hsProS4Vert->GetXaxis()->SetTitleSize(0.06);
  hsProS4Vert->GetYaxis()->SetTitleSize(0.06);
  cspros4vert->SetGridx();
  cspros4vert->SetGridy();
  leg->Draw();
  hsProS4Vert->Write();
  cspros4vert->Print(Form("%s/proS4Vert.png",saveDir));
  cspros4vert->Print(Form("%s/proS4Vert.pdf",saveDir));
  cspros4vert->Print(Form("%s/proS4Vert.tex",saveDir));
  TCanvas *cspis4vert  = new TCanvas("cspis4vert");
  hsPiS4Vert->Draw("hist e nostack");
  hsPiS4Vert->GetXaxis()->SetLabelSize(0.05);
  hsPiS4Vert->GetYaxis()->SetLabelSize(0.05);
  hsPiS4Vert->GetXaxis()->SetTitleSize(0.05);
  hsPiS4Vert->GetYaxis()->SetTitleSize(0.05);
  cspis4vert->SetGridx();
  cspis4vert->SetGridy();
  leg->Draw();
  hsPiS4Vert->Write();
  cspis4vert->Print(Form("%s/piS4Vert.png",saveDir));
  cspis4vert->Print(Form("%s/piS4Vert.pdf",saveDir));
  cspis4vert->Print(Form("%s/piS4Vert.tex",saveDir));
  TCanvas *cspros4horz = new TCanvas("cspros4horz");
  hsProS4Horz->Draw("hist e nostack");
  hsProS4Horz->GetXaxis()->SetLabelSize(0.06);
  hsProS4Horz->GetYaxis()->SetLabelSize(0.06);
  hsProS4Horz->GetXaxis()->SetTitleSize(0.06);
  hsProS4Horz->GetYaxis()->SetTitleSize(0.06);
  cspros4horz->SetLeftMargin(0.13);
  cspros4horz->SetBottomMargin(0.13);
  cspros4horz->SetGridx();
  cspros4horz->SetGridy();
  legProS4Horz->Draw();
  legProS4Horz->Write("legproS4Horz");
  hsProS4Horz->Write();
  cspros4horz->Print(Form("%s/proS4Horz.png",saveDir));
  cspros4horz->Print(Form("%s/proS4Horz.pdf",saveDir));
  cspros4horz->Print(Form("%s/proS4Horz.tex",saveDir));
  TCanvas *cspis4horz  = new TCanvas("cspis4horz");
  hsPiS4Horz->Draw("hist e nostack");
  hsPiS4Horz->GetXaxis()->SetLabelSize(0.06);
  hsPiS4Horz->GetYaxis()->SetLabelSize(0.06);
  hsPiS4Horz->GetXaxis()->SetTitleSize(0.06);
  hsPiS4Horz->GetYaxis()->SetTitleSize(0.06);
  cspis4horz->SetLeftMargin(0.13);
  cspis4horz->SetBottomMargin(0.13);
  cspis4horz->SetGridx();
  cspis4horz->SetGridy();
  legPiS4Horz->Draw();
  legPiS4Horz->Write("legpiS4Horz");
  hsPiS4Horz->Write();
  cspis4horz->Print(Form("%s/piS4Horz.png",saveDir));
  cspis4horz->Print(Form("%s/piS4Horz.pdf",saveDir));
  cspis4horz->Print(Form("%s/piS4Horz.tex",saveDir));

  TCanvas *csratios4vert = new TCanvas("csratios4vert");
  hsRatioS4Vert->Draw("hist e nostack");
  hsRatioS4Vert->GetXaxis()->SetLabelSize(0.06);
  hsRatioS4Vert->GetYaxis()->SetLabelSize(0.06);
  hsRatioS4Vert->GetXaxis()->SetTitleSize(0.06);
  hsRatioS4Vert->GetYaxis()->SetTitleSize(0.06);
  csratios4vert->SetLeftMargin(0.13);
  csratios4vert->SetBottomMargin(0.13);
  csratios4vert->SetGridx();
  csratios4vert->SetGridy();
  legRatioVert->Draw();
  leg->Write("legOrdinary");
  hsRatioS4Vert->Write();
  csratios4vert->Print(Form("%s/ratioS4Vert.png", saveDir));
  csratios4vert->Print(Form("%s/ratioS4Vert.pdf", saveDir));
  csratios4vert->Print(Form("%s/ratioS4Vert.tex", saveDir));
  TCanvas *csratios4horz = new TCanvas("csratios4horz");
  hsRatioS4Horz->Draw("hist e nostack");
  hsRatioS4Horz->GetXaxis()->SetLabelSize(0.06);
  hsRatioS4Horz->GetYaxis()->SetLabelSize(0.06);
  hsRatioS4Horz->GetXaxis()->SetTitleSize(0.06);
  hsRatioS4Horz->GetYaxis()->SetTitleSize(0.06);
  csratios4horz->SetLeftMargin(0.13);
  csratios4horz->SetBottomMargin(0.13);
  csratios4horz->SetGridx();
  csratios4horz->SetGridy();
  leg->Draw();
  hsRatioS4Horz->Write();
  csratios4horz->Print(Form("%s/ratioS4Horz.png", saveDir));
  csratios4horz->Print(Form("%s/ratioS4Horz.pdf", saveDir));
  csratios4horz->Print(Form("%s/ratioS4Horz.tex", saveDir));
} // angularDistS4
