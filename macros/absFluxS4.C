// asbFluxS4.C
// Calculate the absolute flux of particles in S4
// with no regard for if there was an upstream hit
void absFluxS4 (const char* saveDir, 
		const char* dstofDir="/scratch0/dbrailsf/temp/mylinktodtof/",
		const char* ustofDir="/zfs_home/sjones/mylinktoutof/") 
{
  gROOT->SetBatch(kTRUE);

  // Unix timestamps for variable block moves
  // File names for the utof
  // 0.8GeV/c, 0 blocks
  const double start0Block = 1535713289; 
  const double end0Block   = 1535716132;
  // 0.8GeV/c, 1 block
  const double start1Block = 1535796057;
  const double end1Block   = 1535799112;
  // 0.8GeV/c, 2 blocks
  const double start2Block = 1535789157;
  const double end2Block   = 1535792026;
  // 0.8GeV/c, 3 block
  const double start3Block = 1535792404;
  const double end3Block   = 1535795300;
  // 0.8GeV/c, 4 block
  // 4 moderator blocks with -4cm bend
  const double start4Block = 1535836129;
  const double end4Block   = 1535879634;
  // Long sample of 4 moderator block with +5cm bend 
  //  const double start4Block = 1536537600; 
  //  const double end4Block   = 1536669600;
  // S4 positions
  // S1 -> S4 baseline length
  const double baselineS1S4End   = 13.9426;
  const double baselineS1S4Start = 14.0069;
  const double s4OffAxisStartX = 0.121;
  const double s4OffAxisEndX   = 1.4224;
  // Edges of S3 in beam coordinate system
  const double s3StartX   = -0.5168;
  const double s3EndX     = 0.9970;
  const double s3s1StartY = 9.0569 + 1.77;
  const double s3s1EndY   = 8.9146 + 1.77;
  // Shift in ns required to to pion peak at speed of light
  const double dstofShift = 40.;
  // Ustof-dstof cable delay
  const double ustofDelay = 184.7;
  // S4 solid angle coverage
  const double s4SolidAngle = 0.00557;
  // Define the runs to be used for varying number of blocks for ustof
  const char* str0Block = "Data_2018_8_31_b2_800MeV_0block.root";
  const char* str1Block = "Data_2018_9_1_b4_800MeV_1block_bend4cm.root";
  const char* str2Block = "Data_2018_9_1_b2_800MeV_2block_bend4cm.root";
  const char* str3Block = "Data_2018_9_1_b3_800MeV_3block_bend4cm.root";
  const char* str4Block = "Data_2018_9_1_b8_800MeV_4block_bend4cm.root";

  TFile *fout = new TFile(Form("%s/absFluxS4Plots.root", saveDir), "recreate");
  THStack *hs = new THStack("hsX","Absolute particle flux in S4; x / m; Events / spill");
  THStack *hsAngle = new THStack("hsAngle","Absolute particle flux in S4; #theta / degrees; Events / spill / sr");
  //  THStack *hsCosmics = new THStack("hsCosmics", "Cosmic flux in S4 across detector; x / cm; Events / s");

  TLegend *legHorz = new TLegend(0.6, 0.8, 0.85, 0.6);
  TLegend *legAngle = new TLegend(0.58, 0.88, 0.88, 0.6);
  
  for (int nBlocks = 0; nBlocks <= 4; nBlocks++) {
    cout<<"Block "<<nBlocks<<endl;
    int nSpills = 0;
    int nSpillsTrue = 0;
    double lastSpill = 0.;

    TH1D *hCoins = new TH1D(Form("hCoins%d",nBlocks), Form("Coincidences in S4 bars, %d blocks; Bar in S4; Events", nBlocks), 10, 0.5, 10.5);
    hCoins->Sumw2();
    TH1D *hHits  = new TH1D(Form("hHits%d",nBlocks), Form("Hits in S4 bars, %d blocks; Bar in S4; Events", nBlocks), 10, 0.5, 10.5);
    hHits->Sumw2();
    TH1D *hEff   = new TH1D(Form("hEff%d",nBlocks), Form("Efficiencies of S4 bars, %d blocks; Bar in S4; Efficiency", nBlocks), 10, 0.5, 10.5);
    TH1D *habsFluxX = new TH1D(Form("habsFluxX%d",nBlocks), Form("Absolute particle flux in S4, %d blocks; x / m; Events / spill", nBlocks), 20, 0., 1.4);
    habsFluxX->Sumw2();
    TH1D *habsFluxXAngle = new TH1D(Form("habsFluxXAngle%d",nBlocks), Form("Absolute particle flux in S4, %d blocks; #theta / degrees; Events / spill", nBlocks), 20, 0., 6.);
    habsFluxXAngle->Sumw2();

    TH2D *h2Cosmics = new TH2D(Form("h2Cosmics%d",nBlocks), Form("Cosmic flux, %d blocks; x / cm; Bar; Hz",nBlocks), 40, 0, 140, 10, 0.5, 10.5);
    TH1D *hCosmicsVert = new TH1D(Form("hCosmicsVert%d",nBlocks), Form("Cosmic flux, %d blocks; x / cm; Hz",nBlocks), 10, 0.5, 10.5);
    hCosmicsVert->Sumw2();
    TH1D *hCosmicsHorz = new TH1D(Form("hCosmicsHorz%d",nBlocks), Form("Cosmic flux, %d blocks; x / cm; Hz",nBlocks), 40, 0, 140);
    hCosmicsHorz->Sumw2();
    // Find the correct dstof files
    Int_t runMin=-1;
    Int_t runMax=-1;

    char* nustof;
    double startTime = 0;
    double endTime   = 0;
    if (nBlocks == 0) {
      startTime = start0Block;
      endTime   = end0Block;
      nustof = Form("%s/%s", ustofDir, str0Block);
    }
    else if (nBlocks == 1) {
      startTime = start1Block;
      endTime   = end1Block;
      nustof = Form("%s/%s", ustofDir, str1Block);
    }
    else if (nBlocks == 2) {
      startTime = start2Block;
      endTime   = end2Block;
      nustof = Form("%s/%s", ustofDir, str2Block);
    }
    else if (nBlocks == 3) {
      startTime = start3Block;
      endTime   = end3Block;
      nustof = Form("%s/%s", ustofDir, str3Block);
    }
    else if (nBlocks == 4) {
      startTime = start4Block;
      endTime   = end4Block;
      nustof = Form("%s/%s", ustofDir, str4Block);
    }

    for (int irun=950; irun<1400; irun++){
      TFile *fin = new TFile(Form("%srun%d/DsTOFcoincidenceRun%d_tdc1.root", dstofDir, irun, irun), "read");
      RawDsTofCoincidence *tofCoinTemp = NULL;
      TTree *tree = (TTree*) fin->Get("tofCoinTree");
      //      tree->SetDirectory(0);
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

    // Calculating cosmic fluxes
    // Two different cases -- 0 block data with only 3 bars amplified
    // Everything else -- all bars apart from 1, 2, and 10 amplified
    cout<<"Calculating cosmic flux"<<endl;
    int cosmicRunMin = -1;
    int cosmicRunMax = -1;
    if (nBlocks == 0) {
      cosmicRunMin = 999;
      cosmicRunMax = 1019;
    }
    else {
      cosmicRunMin = 1035;
      cosmicRunMax = 1045;
    }
    // Now loop over these files and find the cosmic fluxes
    int nSpillsCos = 0;    
    for (int itdc=0; itdc<2; itdc++) {
      nSpillsCos = 0;
      double lastSpillCos = 0.;
      //for (int irun=runMin; irun<runMax+1; irun++){
      TChain *tofCoinChain = new TChain("tofCoinTree");
      for (int run = cosmicRunMin; run < cosmicRunMin+1; run++) {
	tofCoinChain->Add(Form("%srun%d/DsTOFcoincidenceRun%d_tdc%d.root", dstofDir, run, run, itdc+1));
      }
      RawDsTofCoincidence *tofCoin = NULL;
      tofCoinChain->SetBranchAddress("tofCoin", &tofCoin);
      for (int t = 0; t < tofCoinChain->GetEntries(); t++) {
	tofCoinChain->GetEntry(t);
	// Look at the distribution of cosmics that register as hits (out of spill)
	if (!tofCoin->inSpill) {
	  double cosmicPosition = (tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 70.;
	  h2Cosmics->Fill(cosmicPosition, tofCoin->bar);
	  hCosmicsVert->Fill(tofCoin->bar);
	  hCosmicsHorz->Fill(cosmicPosition);
	} // if (!tofCoin->inSpill)
      } // for (int t = 0; t < tofCoinChain->GetEntries(); t++) 
      delete tofCoin;
      delete tofCoinChain;
      if (lastSpillCos != tofCoin->lastDelayedBeamSignal && tofCoin->lastDelayedBeamSignal - lastSpillCos > 1e9) {
	lastSpillCos = tofCoin->lastDelayedBeamSignal;
	nSpillsCos++;
      }
    } // for (int itdc=0; itdc<2; itdc++)

    // Scale by running time
    h2Cosmics->Scale( 1. / ((cosmicRunMax - cosmicRunMin)*3600. - nSpillsCos) );    
    hCosmicsVert->Scale( 1. / ((cosmicRunMax - cosmicRunMin)*3600. - nSpillsCos) );    
    hCosmicsHorz->Scale( 1. / ((cosmicRunMax - cosmicRunMin)*3600. - nSpillsCos) );    
    fout->cd();

    TCanvas *c2c = new TCanvas(Form("c2c_%d",nBlocks));
    c2c->SetRightMargin(0.13);
    gStyle->SetPalette(55);
    h2Cosmics->Draw("colz");
    h2Cosmics->Write();
    c2c->Print(Form("%s/%d_cosmics.png", saveDir, nBlocks));
    c2c->Print(Form("%s/%d_cosmics.pdf", saveDir, nBlocks));
    TCanvas *ccvert = new TCanvas(Form("ccvert_%d",nBlocks));
    hCosmicsVert->Draw("hist e");
    hCosmicsVert->Write();
    ccvert->Print(Form("%s/%d_cosmicsVert.png", saveDir, nBlocks));
    ccvert->Print(Form("%s/%d_cosmicsVert.pdf", saveDir, nBlocks));
    TCanvas *cchorz = new TCanvas(Form("cchorz_%d",nBlocks));
    hCosmicsHorz->Draw("hist e");
    hCosmicsHorz->Write();
    cchorz->Print(Form("%s/%d_cosmicsHorz.png", saveDir, nBlocks));
    cchorz->Print(Form("%s/%d_cosmicsHorz.pdf", saveDir, nBlocks));

    TH1D *hCosEff = (TH1D*)hCosmicsHorz->Clone("hCosEff");
    hCosEff->Scale(1. / hCosEff->GetBinContent(hCosEff->GetMaximumBin()));
    hCosEff->Write();

    TH1D *h2CosEff = (TH1D*)h2Cosmics->Clone("h2CosEff");
    h2CosEff->Scale(1. / h2CosEff->GetBinContent(h2CosEff->GetMaximumBin()));
    h2CosEff->Write();
    for (int i = 0; i < h2CosEff->GetNbinsX(); i++) {
      for (int j = 0; j < h2CosEff->GetNbinsY(); j++) {
	if (h2CosEff->GetBinContent(i, j) == 0) {
	  h2CosEff->SetBinContent(i, j, 1);
	}
      }
    }

    int nUtofHits = 0;
    // Loop to calculate efficiencies    
    for (int itdc=0; itdc<2; itdc++) {
      double tempUstof;
      double tempBeam;
      double ustofNs;
      double beamNs;
      for (int irun=runMin; irun<runMax+1; irun++){
	double lastBeam = 0.;
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
	//	TTree *beamTree = new TTree("beamTree", "beam");
	//	beamTree->SetDirectory(0);
	//	beamTree->Branch("beamNs", &beamNs, "beamNs/D");
	//	tempBeam=0;
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
	    // Is in spill
	    if ((ustofNs - lastBeam) < 1e9 && (ustofNs - lastBeam) > 0.) {
	      nUtofHits++;
	    }
	  } // if (tof->channel == 13) 
	  else if (tof->channel == 14) {
	    lastBeam = tof->fakeTimeNs;
	  } // else if (tof->channel == 14)
	} // for (int i=0; i<tof->GetEntries(); i++)
	ustofTree->BuildIndex("ustofNs");
	// Now loop over coincidence file to find number of coincidences for each bar
	for (int t=0; t<tofCoinTree->GetEntries(); t++) {
	  tofCoinTree->GetEntry(t);
	  if (tofCoin->unixTime[0] < startTime) continue;
	  if (tofCoin->unixTime[0] > endTime) break;
	  double deltat = TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1]  );
	  double dstofHitT = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (10. - TMath::Abs(deltat) / 2. );
	  double tofCalc = dstofHitT - tofCoin->usTofSignal;
	  if (tofCalc > 70. && tofCalc < 200. && tofCoin->bar != 10) {
	    hCoins->Fill(tofCoin->bar);
	  } // if (tofCalc > 70. && tofCalc < 200.)
	} // for (int t=0; t<tofTree->GetEntries(); t++) 
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
	tofFile->Close();
	tofCoinFile->Close();
	delete tofFile;
	delete tofCoinFile;
	delete ustofTree;
      } // for (int irun=runMin; irun<runMax+1; irun++)
    } // for (int itdc=0; itdc<2; itdc++) 

    fout->cd(0);
    TCanvas *cHits = new TCanvas(Form("cHits_%d",nBlocks));
    hHits->Draw("hist e");
    hHits->Write();
    cHits->Print(Form("%s/%d_barHits.png",saveDir,nBlocks));
    cHits->Print(Form("%s/%d_barHits.pdf",saveDir,nBlocks));
    TCanvas *cCoins = new TCanvas(Form("cCoins_%d",nBlocks));
    hCoins->Draw("hist e");
    hCoins->Write();
    cCoins->Print(Form("%s/%d_barCoins.png",saveDir,nBlocks));
    cCoins->Print(Form("%s/%d_barCoins.pdf",saveDir,nBlocks));
    TCanvas *cEff = new TCanvas(Form("cEff_%d",nBlocks));
    hHits->Add(hCoins, -1.);
    hEff->Divide(hCoins, hHits, 1., 1., "B");
    hEff->Draw("hist e");
    hEff->Write();
    cEff->Print(Form("%s/%d_barEff.png",saveDir,nBlocks));
    cEff->Print(Form("%s/%d_barEff.pdf",saveDir,nBlocks));

    for (int itdc=0; itdc<2; itdc++) {
      //for (int irun=runMin; irun<runMax+1; irun++){
      TChain *tofCoinChain = new TChain("tofCoinTree");
      for (int irun=runMin; irun<runMax+1; irun++){
	tofCoinChain->Add(Form("%srun%d/DsTOFcoincidenceRun%d_tdc%d.root", dstofDir, irun, irun, itdc+1));
      }
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
	  } // for (int sp=ientry; sp<tof1CoinChain->GetEntries(); sp++) 
	} // if (tof1Coin->lastDelayedBeamSignal != lastSpill && itdc == 0)
	// Check if this coincidence is within a spill
	if (tofCoin->inSpill && tofCoin->bar != 10) {
	  double deltat = TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1]  );
	  double dstofHitT = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (10. - TMath::Abs(deltat) / 2. );
	  // For true x, y position (relative to beam axis) interpolate between two 
	  // measured positions of the ToF (see survey data)
	  double positionX = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 65.) / 130.) * (s4OffAxisEndX - s4OffAxisStartX) + s4OffAxisStartX;
	  double positionXP = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 70.));
	  double positionY = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 65.) / 130.) * (baselineS1S4End - baselineS1S4Start) + baselineS1S4Start;
	  habsFluxX->Fill(positionX, 1. / hEff->GetBinContent(tofCoin->bar));
	  double angleOffAxis = TMath::ATan(positionX / positionY) * 180. / TMath::Pi();
	  habsFluxXAngle->Fill(angleOffAxis, 1. / /*h2CosEff->GetBinContent( h2CosEff->FindBin(positionXP, tofCoin->bar)) (*/hEff->GetBinContent(tofCoin->bar) /** hCosEff->GetBinContent( hCosEff->GetXaxis()->FindBin(positionXP)) )*/);
	} // if (tofCoin->inSpill)
      } // for (int h=0; h<tofCoinChain->GetEntries(); h++) 
      delete tofCoin;
      delete tofCoinChain;
    } // for (int itdc=0; itdc<2; itdc++)
    
    // Save all of these flux plots individually
    fout->cd();
    TCanvas *c1 = new TCanvas(Form("c1_%d",nBlocks));
    // Want to know what the flux is at theta = x = 0
    // Fit to the curves but only between x = 0.4 and x = 1.2
    // theta = 1.9 and theta = 5.2
    const double cutXLow = 0.4;
    const double cutXHi  = 1.2;
    const double cutThetaLow = 1.9;
    const double cutThetaHi  = 5.2;
    TF1 *fX = new TF1(Form("fX%d", nBlocks),"pol1", cutXLow, cutXHi);
    habsFluxX->Scale(1. / nSpillsTrue);
    habsFluxX->Draw("hist e");
    habsFluxX->Fit(fX,"R");
    fX->Draw("same");
    habsFluxX->Write();
    fX->Write();
    c1->Print(Form("%s/%d_absFluxX.png", saveDir, nBlocks));
    c1->Print(Form("%s/%d_absFluxX.pdf", saveDir, nBlocks));
    c1->Print(Form("%s/%d_absFluxX.tex", saveDir, nBlocks));

    
    TCanvas *c2 = new TCanvas(Form("c2_%d",nBlocks));
    TF1 *fTheta = new TF1(Form("fTheta%d", nBlocks),"pol1", cutThetaLow, cutThetaHi);
    habsFluxXAngle->Scale(1. / (nSpillsTrue));
    habsFluxXAngle->Draw("hist");
    habsFluxXAngle->Write();
    //    habsFluxXAngle->Fit(fTheta,"R");
    //fTheta->Draw("same");
    //habsFluxXAngle->Write();
    //    fTheta->Write();
    c2->Print(Form("%s/%d_absFluxXAngle.png", saveDir, nBlocks));
    c2->Print(Form("%s/%d_absFluxXAngle.pdf", saveDir, nBlocks));

    habsFluxX->SetLineWidth(2);
    habsFluxXAngle->SetLineWidth(2);
    // And do them all together in different colours
    if (nBlocks == 0) {
      habsFluxX->SetLineColor(kBlack);
      habsFluxXAngle->SetLineColor(kBlack);
      legHorz->AddEntry(habsFluxX, "0 blocks", "l");
      hs->Add(habsFluxX);
      double intAngle = habsFluxXAngle->Integral();
      legAngle->AddEntry(habsFluxXAngle, Form("0 blocks - %d per spill ", (int)intAngle), "l"); 
    }
    else if (nBlocks == 1) {
      habsFluxX->SetLineColor(kRed);
      habsFluxXAngle->SetLineColor(kRed);
      legHorz->AddEntry(habsFluxX, "1 block", "l");
      hs->Add(habsFluxX);
      double intAngle = habsFluxXAngle->Integral();
      legAngle->AddEntry(habsFluxXAngle, Form("1 block - %d per spill ", (int)intAngle), "l"); 
    }
    else if (nBlocks == 2) {
      habsFluxX->SetLineColor(kBlue);
      habsFluxXAngle->SetLineColor(kBlue);
      legHorz->AddEntry(habsFluxX, "2 blocks", "l");
      hs->Add(habsFluxX);
      double intAngle = habsFluxXAngle->Integral();
      legAngle->AddEntry(habsFluxXAngle, Form("2 blocks - %d per spill ", (int)intAngle), "l"); 
    }
    else if (nBlocks == 3){
      habsFluxX->SetLineColor(kCyan+1);
      habsFluxXAngle->SetLineColor(kCyan+1);
      legHorz->AddEntry(habsFluxX, "3 blocks", "l");
      hs->Add(habsFluxX);
      double intAngle = habsFluxXAngle->Integral();
      legAngle->AddEntry(habsFluxXAngle, Form("3 blocks - %d per spill ", (int)intAngle), "l"); 
    }
    else {
      habsFluxX->SetLineColor(kOrange+1);
      habsFluxXAngle->SetLineColor(kOrange+1);
      legHorz->AddEntry(habsFluxX, "4 blocks", "l");
      hs->Add(habsFluxX);
      double intAngle = habsFluxXAngle->Integral();
      legAngle->AddEntry(habsFluxXAngle, Form("4 blocks - %d per spill", (int)intAngle), "l"); 
    }
    hsAngle->Add(habsFluxXAngle);

    cout<<"Total of "<<nSpills<<" ("<<nSpillsTrue<<" true) for "<<nBlocks<<" blocks"<<endl;
    /*
    double UtofPerSpill = (double)nUtofHits / (double)nSpillsTrue;
    cout<<UtofPerSpill<<" S1 x S2 hits per spill"<<endl;
    cout<<habsFluxXAngle->Integral()<<" S4 per spill (raw)"<<endl;

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

    TH1D *hXAngleS1S2 = new TH1D(Form("hXAngleS1S2%d", nBlocks), Form("Angular distribution of hits in S3 (S1 & S2 triggers), %d blocks; #theta / degrees; Events / spill", nBlocks), 100, -3.8, 6.2);
    for (int t=0; t<utree->GetEntries(); t++) {
      utree->GetEntry(t);
      // Has an S1 and S2 hit
      if (tTrig !=0 ) {
	for (int n=0; n<nhit; n++) {
	  double positionX = ((xToF[n] - 4.) / 152.)*(s3EndX - s3StartX) + s3StartX;
	  double positionY = ((xToF[n] - 4.) / 152.)*(s3s1EndY - s3s1StartY)+s3s1StartY;
	  double angleOffAxis = TMath::ATan(positionX / positionY) * 180. / TMath::Pi();
	  hXAngleS1S2->Fill(angleOffAxis);
	} // for (int n=0; n<nhit; n++) 
      } // if (tTrig !=0 ) 
    } // for (int t=0; t<utree->GetEntries(); t++)

    // Integrate this hist between the angular limits of S2 then scale by this number
    const double s2ThetaLow = -0.359;
    const double s2ThetaHi  = 3.957;
    const double s4ThetaLow = 0.401;
    const double s4ThetaHi  = 6.083;
    double s2Int = hXAngleS1S2->Integral(hXAngleS1S2->GetBin(s2ThetaLow), hXAngleS1S2->GetBin(s2ThetaHi));
    hXAngleS1S2->Scale(1. / s2Int);
    cout<<"Integral "<<hXAngleS1S2->Integral(hXAngleS1S2->GetBin(s2ThetaLow), hXAngleS1S2->GetBin(s2ThetaHi))<<endl;
    double s4Int = hXAngleS1S2->Integral(hXAngleS1S2->GetBin(s4ThetaLow), hXAngleS1S2->GetBin(s4ThetaHi));
    cout<<"S4 correction factor "<<s4Int<<endl;
    
    finustof->Close();
    delete finustof;
    */
  } // nBlocks
  TCanvas *cs = new TCanvas("cs");
  hs->Draw("hist nostack");
  legHorz->Draw();
  hs->Write();
  cs->Print(Form("%s/absFluxX.png", saveDir));
  cs->Print(Form("%s/absFluxX.pdf", saveDir));
  TCanvas *csAngle = new TCanvas("csAngle");
  hsAngle->Draw("hist e nostack");
  hsAngle->GetXaxis()->SetLabelSize(0.04);
  hsAngle->GetYaxis()->SetLabelSize(0.04);
  hsAngle->GetXaxis()->SetTitleSize(0.04);
  hsAngle->GetYaxis()->SetTitleSize(0.04);
  legAngle->Draw();
  legAngle->Write("legAngle");
  hsAngle->Write();
  csAngle->Print(Form("%s/absFluxXAngle.png", saveDir));
  csAngle->Print(Form("%s/absFluxXAngle.pdf", saveDir));
  csAngle->Print(Form("%s/absFluxXAngle.tex", saveDir));

} // absFluxS4
