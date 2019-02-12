// angularDistS4.C
// Angular distribution of protons and pions for different moderator blocks
// Use efficiency calculation and background subtraction
void angularDistS4(const char* saveDir, 
		   const char* dstofDir="/scratch0/dbrailsf/temp/mylinktodtof/") 
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
  const double piLow  = 40.;
  const double piHi   = 55.;
  const double proLow = 66.;
  const double proHi  = 105.;
  // S1 -> S4 positions
  const double baselineS1S4End   = 13.9426;
  const double baselineS1S4Start = 14.0069;
  const double s4OffAxisStartX   = 0.121;
  const double s4OffAxisEndX     = 1.4224;
  // z coordinates
  const double s4BarTop    = 0.43; 
  const double s4BarBottom = -0.35;
  // Shift in ns required to to pion peak at speed of light
  const double dstofShift = 40.;
  // Ustof-dstof cable delay
  const double ustofDelay = 184.7;

  TFile *fout = new TFile(Form("%s/angularDistS4Plots.root", saveDir), "recreate");

  THStack *hsBkgSub = new THStack("hsBkgSub", "Background subtracted and efficiency corrected S4 ToF spectrum; S4 - S1 / ns; Events / spill");
  THStack *hsDtof   = new THStack("hsDtof", "S4 ToF spectrum; S4 - S1 / ns; Events / spill");
  TLegend *leg = new TLegend(0.65, 0.65, 0.85, .85);

  for (int nBlocks = 0; nBlocks < 4; nBlocks++) {
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
	TTree *ustofTree = new TTree("ustofTree", "ustof");;
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
    hHits->Draw("hist");
    hHits->Write();
    cHits->Print(Form("%s/%d_barHits.png",saveDir,nBlocks));
    cHits->Print(Form("%s/%d_barHits.pdf",saveDir,nBlocks));
    TCanvas *cCoins = new TCanvas(Form("%d_cCoins",nBlocks));
    hCoins->Draw("hist");
    hCoins->Write();
    cCoins->Print(Form("%s/%d_barCoins.png",saveDir,nBlocks));
    cCoins->Print(Form("%s/%d_barCoins.pdf",saveDir,nBlocks));
    TCanvas *cEff = new TCanvas(Form("%d_cEff",nBlocks));
    hHits->Add(hCoins, -1.);
    hEff->Divide(hCoins, hHits, 1., 1., "B");
    hEff->Draw("hist");
    hEff->Write();
    cEff->Print(Form("%s/%d_barEff.png",saveDir,nBlocks));
    cEff->Print(Form("%s/%d_barEff.pdf",saveDir,nBlocks));

    // 1D histograms in the horizontal and vertical direction of S4
    // Horizontal
    TH1D *hProS4Horz = new TH1D(Form("hProS4Horz%d",nBlocks), Form("Horizontal angular distribution of proton hits in S4, %d blocks; #theta / degrees; Events / spill",nBlocks), 40, 0, 6.);
    hProS4Horz->Sumw2();
    TH1D *hPiS4Horz  = new TH1D(Form("hPiS4Horz%d",nBlocks), Form("Horizontal angular distribution of MIP hits in S4, %d blocks; #theta / degrees; Events / spill",nBlocks), 40, 0, 6);
    hPiS4Horz->Sumw2();
    TH1D *hProPiRatioS4Horz  = new TH1D(Form("hProPiRatioS4Horz%d",nBlocks), Form("Horizontal angular distribution of proton/MIP ratio in S4, %d blocks; #theta / degrees; Protons/MIPs",nBlocks), 40, 0, 6);
    // Vertical
    TH1D *hProS4Vert = new TH1D(Form("hProS4Vert%d",nBlocks), Form("Vertical angular distribution of proton hits in S4, %d blocks; #phi / degrees; Events / spill",nBlocks), 10, -1.5, 1.8);
    hProS4Vert->Sumw2();
    TH1D *hPiS4Vert  = new TH1D(Form("hPiS4Vert%d",nBlocks), Form("Vertical angular distribution of MIP hits in S4, %d blocks; #phi / degrees; Events / spill",nBlocks), 10, -1.5, 1.8);
    hPiS4Vert->Sumw2();
    TH1D *hProPiRatioS4Vert  = new TH1D(Form("hProPiRatioS4Vert%d",nBlocks), Form("Vertical angular distribution of proton/MIP ratio in S4, %d blocks; #phi / degrees; Protons/MIPs",nBlocks), 10, -1.5, 1.8);

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

	if (tofCoin->lastDelayedBeamSignal != lastSpill && itdc == 0) {
	  lastSpill = tofCoin->lastDelayedBeamSignal;
	  nSpills++;
	  for (int sp=h; sp<tofCoinChain->GetEntries(); sp++) {
	    tofCoinChain->GetEntry(sp);
	    if (tofCoin->unixTime[0]<startTime) continue;
	    if (tofCoin->unixTime[0]>endTime) break;

	    if ((tofCoin->fakeTimeNs[0] - tofCoin->usTofSignal) < 200. &&
		(tofCoin->fakeTimeNs[0] - tofCoin->usTofSignal) > 70. &&
		(tofCoin->fakeTimeNs[0] - lastSpill) < 1e9 &&
		(tofCoin->fakeTimeNs[0] - lastSpill) > 0) {
	      nSpillsTrue++;
	      break;
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
	  if (tofCalc < piHi & tofCalc > piLow) { 
	    nPi += (1. / hEff->GetBinContent(tofCoin->bar));
	    // Calculate position of hit in global coordinates
	    double positionX = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 65.) / 130.) * (s4OffAxisEndX - s4OffAxisStartX) + s4OffAxisStartX;
	    double positionY = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 65.) / 130.) * (baselineS1S4End - baselineS1S4Start) + baselineS1S4Start;
	    double positionZ = (tofCoin->bar * 0.075) - 0.375 - 0.01;
	    // Calculate the angles relative to the nominal beamline
	    double angleTheta = TMath::ATan(positionX / positionY) * (180./TMath::Pi());
	    double anglePhi   = TMath::ATan(positionZ / positionY) * (180./TMath::Pi());
	    hPiS4Horz->Fill(angleTheta, hEff->GetBinContent(tofCoin->bar));
	    hPiS4Vert->Fill(anglePhi,   hEff->GetBinContent(tofCoin->bar));
	  } // if (tofCalc < piHi & tofCalc > piLow)
	  else if (tofCalc < proHi & tofCalc > proLow) {
	    nP += (1. / hEff->GetBinContent(tofCoin->bar));
	    // Calculate position of hit in global coordinates
	    double positionX = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 65.) / 130.) * (s4OffAxisEndX - s4OffAxisStartX) + s4OffAxisStartX;
	    double positionY = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 65.) / 130.) * (baselineS1S4End - baselineS1S4Start) + baselineS1S4Start;
	    double positionZ = (tofCoin->bar * 0.075) - 0.375 - 0.01;
	    // Calculate the angles relative to the nominal beamline
	    double angleTheta = TMath::ATan(positionX / positionY) * (180./TMath::Pi());
	    double anglePhi   = TMath::ATan(positionZ / positionY) * (180./TMath::Pi());
	    hProS4Horz->Fill(angleTheta, hEff->GetBinContent(tofCoin->bar));
	    hProS4Vert->Fill(anglePhi,   hEff->GetBinContent(tofCoin->bar));
	  } // else if (tofCalc < proHi & tofCalc > proLow) 
	} // if (tofCalc < 160. && tofCalc > 30.) 
      } // for (int h=0; h<tofCoinChain->GetEntries(); h++) 
      delete tofCoin;
      delete tofCoinChain;
    } // for (int itdc=0; itdc<2; itdc++)
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

    // Now we have the fit values, loop over again and subtract the background
    // For each bin, find the fraction of each particle type which are background and 
    // this fraction from the bin
    TF1 *fSub = new TF1("fSub", "exp([0]+[1]*x)", 30, 160);
    fSub->SetParameter(0, fSplusBExp->GetParameter("bkgconst"));  
    fSub->SetParameter(1, fSplusBExp->GetParameter("bkgdecay"));
    // Background subtracted tof spectrum
    TH1D *hdtof1d_sub = (TH1D*)hdtof1d->Clone(Form("hdtof1d_sub_%d",nBlocks));
    hdtof1d_sub->Add(fSub, -1.);
    // If the bin content drops below 0, set to 0
    for (int b = 0; b <=  hdtof1d_sub->GetNbinsX(); b++) {
      if (hdtof1d_sub->GetBinContent(b) < 0.) {
	hdtof1d_sub->SetBinContent(b, 0);
      } // if (hdtof1d_sub->GetBinContent(b) < 0.)
    } // for (int b = 0; hdtof1d_sub->GetNbinsX(); b++)
    TCanvas *cdtof_sub = new TCanvas(Form("%d_cdtof_sub",nBlocks));
    cdtof_sub->SetLogy();
    hdtof1d_sub->Draw("hist");
    hdtof1d_sub->Write();
    cdtof_sub->Print(Form("%s/%d_dtof1d_bkgSub.png", saveDir, nBlocks));
    cdtof_sub->Print(Form("%s/%d_dtof1d_bkgSub.pdf", saveDir, nBlocks));

    TCanvas *cProS4Horz = new TCanvas(Form("cProS4Horz%d",nBlocks));
    hProS4Horz->Scale(1. / nSpillsTrue);
    hProS4Horz->Draw("hist");
    cProS4Horz->Print(Form("%s/%d_proS4Horz.png",nBlocks,saveDir));
    cProS4Horz->Print(Form("%s/%d_proS4Horz.pdf",nBlocks,saveDir));
    TCanvas *cPiS4Horz  = new TCanvas(Form("cPiS4Horz%d",nBlocks));
    hPiS4Horz->Scale(1. / nSpillsTrue);
    hPiS4Horz->Draw("hist");
    cPiS4Horz->Print(Form("%s/%d_piS4Horz.png",nBlocks,saveDir));
    cPiS4Horz->Print(Form("%s/%d_piS4Horz.pdf",nBlocks,saveDir));
    TCanvas *cProS4Vert = new TCanvas(Form("cProS4Vert%d",nBlocks));
    hProS4Vert->Scale(1. / nSpillsTrue);
    hProS4Vert->Draw("hist");
    cProS4Vert->Print(Form("%s/%d_proS4Vert.png",nBlocks,saveDir));
    cProS4Vert->Print(Form("%s/%d_proS4Vert.pdf",nBlocks,saveDir));
    TCanvas *cPiS4Vert  = new TCanvas(Form("cPiS4Vert%d",nBlocks));
    hPiS4Vert->Scale(1. / nSpillsTrue);
    hPiS4Vert->Draw("hist");
    cPiS4Vert->Print(Form("%s/%d_piS4Vert.png",nBlocks,saveDir));
    cPiS4Vert->Print(Form("%s/%d_piS4Vert.pdf",nBlocks,saveDir));

    if (nBlocks == 0) {
      hdtof1d_sub->SetLineColor(kBlue);
      hdtof1d->SetLineColor(kBlue);
      leg->AddEntry(hdtof1d, "0 blocks", "l");     
    }
    else if (nBlocks == 1) {
      hdtof1d_sub->SetLineColor(kRed);
      hdtof1d->SetLineColor(kRed);
      leg->AddEntry(hdtof1d, "1 block", "l");
    }
    else if (nBlocks == 2) {
      hdtof1d_sub->SetLineColor(kBlack);
      hdtof1d->SetLineColor(kBlack);
      leg->AddEntry(hdtof1d, "2 blocks", "l");
    }
    else if (nBlocks == 3) {
      hdtof1d_sub->SetLineColor(kGreen+2);
      hdtof1d->SetLineColor(kGreen+2);
      leg->AddEntry(hdtof1d, "3 blocks", "l");
    }
    else if (nBlocks ==4) {
      hdtof1d_sub->SetLineColor(kMagenta);
      hdtof1d->SetLineColor(kMagenta);
      leg->AddEntry(hdtof1d, "4 blocks", "l");
    }
    hdtof1d->Scale(1. / nSpillsTrue);
    hdtof1d_sub->Scale(1. / nSpillsTrue);
    hsDtof->Add(hdtof1d);
    hsBkgSub->Add(hdtof1d_sub);
  } // for (int nBlocks = 0; nBlocks < 4; nBlocks++) 
  TCanvas *csDtof = new TCanvas("csDstof");
  csDtof->SetLogy();
  hsDtof->Draw("hist nostack");
  hsDtof->GetYaxis()->SetRangeUser(1e-2, 2e2);
  csDtof->Update();
  leg->Draw();
  csDtof->Print(Form("%s/s4ToF.png", saveDir));
  csDtof->Print(Form("%s/s4ToF.pdf", saveDir));
  TCanvas *csBkgSub = new TCanvas("csBkgSub");
  csBkgSub->SetLogy();
  hsBkgSub->Draw("hist nostack");
  hsBkgSub->GetYaxis()->SetRangeUser(1e-2, 2e2);
  //  hsBkgSub->Draw();
  csBkgSub->Update();
  leg->Draw();
  csBkgSub->Print(Form("%s/s4BkgSub.png", saveDir));
  csBkgSub->Print(Form("%s/s4BkgSub.pdf", saveDir));
} // angularDistS4
