// angularDist.C
// Angular distribution of protons and pions for different moderator blocks
void angularDist (/*const int nBlocks, */ const char* saveDir, bool useEffs = false, const char* dstofDir="/scratch0/dbrailsf/temp/mylinktodtof/") {

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
  // Most runs were in this configuration so don't need to use necessarily
  const double start4Block = 1535608220;
  const double end4Block   = 1535617102;
  // Dstof directory
  //  const char* dsDir = "/scratch0/dbrailsf/temp/mylinktodtof/";
  // Timing cuts
  const double piLow  = 80.;
  const double piHi   = 95.;
  const double proLow = 106.;
  const double proHi  = 134.;

  TCanvas *c_proPiVert = new TCanvas("c_proPiVert");
  TCanvas *c_proPiHorz = new TCanvas("c_proPiHorz");

  TLegend *legHorz = new TLegend(0.6, 0.8, 0.85, 0.6);
  TLegend *legVert = new TLegend(0.2, 0.8, 0.35, 0.6);

  for (int nBlocks=0; nBlocks < 4; nBlocks++) {
    double nSpills = 0.;
    double nSpillsTrue = 0.;
    double lastSpill = 0.;

    vector< vector<double> > pmtHitTimesA;
    vector< vector<double> > pmtHitTimesB;
    pmtHitTimesA.resize(10);
    pmtHitTimesB.resize(10);

    vector<int> pmtHitsA;
    vector<int> pmtHitsB;
    pmtHitsA.resize(10, 0);
    pmtHitsB.resize(10, 0);

    int nPi = 0;
    int nP  = 0;

    TH1D *htof1dHitsA = new TH1D(Form("htof1dHitsA_%d",nBlocks), Form("Time of flight, %d blocks; PMT_{A} - UsToF / ns; Events", nBlocks), 100, 130, 230);
    htof1dHitsA->Sumw2();
    TH1D *htof1dHitsB = new TH1D(Form("htof1dHitsB_%d",nBlocks), Form("Time of flight, %d blocks; PMT_{B} - UsToF / ns; Events", nBlocks), 100, 130, 230);
    htof1dHitsB->Sumw2();
    TH1D *hHitsA = new TH1D(Form("hHitsA_%d",nBlocks), Form("PMT_A + S1,2 coincidences, %d blocks; Bar; Events", nBlocks), 10, 0.5, 10.5);
    hHitsA->Sumw2();
    TH1D *hHitsB = new TH1D(Form("hHitsB_%d",nBlocks), Form("PMT_B + S1,2 coincidences, %d blocks; Bar; Events", nBlocks), 10, 0.5, 10.5);
    hHitsB->Sumw2();
    TH1D *hHits = new TH1D(Form("hHits_%d",nBlocks), Form("PMT_A or PMT_B + S1,2 coincidences, %d blocks; Bar; Events", nBlocks), 10, 0.5, 10.5);
    hHits->Sumw2();
    TH1D *hCoins = new TH1D(Form("hCoins_%d",nBlocks), Form("PMT_A + PMT_B + S1,2 coincidences, %d blocks; Bar; Events", nBlocks), 10, 0.5, 10.5);
    hCoins->Sumw2();    

    TH1D *htof1d = new TH1D(Form("htof1d_%d",nBlocks), Form("Time of flight, %d blocks; DsToF - UsToF / ns; Events / spill", nBlocks), 100, 50, 150);
    TH2D *piHitsDstof = new TH2D("piHitsDstof", Form("%d blocks, position of S4 #pi hits; x / cm; y / cm; Events / spill", nBlocks), 20, 0, 140, 10, 0, 77.5);
    TH2D *proHitsDstof = new TH2D("proHitsDstof", Form("%d blocks, position of S4 proton hits; x / cm; y / cm; Events / spill", nBlocks), 20, 0, 140, 10, 0, 77.5);
    TH2D *proPiDstof = new TH2D("proPiDstof", Form("%d blocks, proton/pion ratio in S4; x / cm; y / cm", nBlocks), 20, 0, 140, 10, 0, 77.5);
    // Binned bar by bar only
    TH2D *piHitsDstofVert = new TH2D("piHitsDstofVert", Form("%d blocks, position of S4 #pi hits; ; Bar; Events / spill", nBlocks), 1, .5, 1.5, 10, 0.5, 10.5);
    TH2D *proHitsDstofVert = new TH2D("proHitsDstofVert", Form("%d blocks, position of S4 proton hits; ; Bar; Events / spill", nBlocks), 1, .5, 1.5, 10, 0.5, 10.5);
    TH2D *proPiDstofVert = new TH2D("proPiDstofVert", Form("%d blocks, proton/pion ratio in S4; ; Bar", nBlocks), 1, .5, 1.5, 10, 0.5, 10.5);
    // Horizontal binning only
    TH2D *piHitsDstofHorz = new TH2D("piHitsDstofHorz", Form("%d blocks, position of S4 #pi hits; x / cm", nBlocks), 20, 0, 140, 1, 0.5, 1.5);
    TH2D *proHitsDstofHorz = new TH2D("proHitsDstofHorz", Form("%d blocks, position of S4 proton hits; x / cm", nBlocks), 20, 0, 140, 1, 0.5, 1.5);
    TH2D *proPiDstofHorz = new TH2D("proPiDstofHorz", Form("%d blocks, proton/pion ratio in S4; x / cm", nBlocks), 20, 0, 140, 1, 0.5, 1.5);

    // 1D hists
    TH1D *hpiHitsDstofVert = new TH1D("hpiHitsDstofVert", Form("%d blocks, position of S4 #pi hits; Bar; Events / spill", nBlocks), 10, 0.5, 10.5);
    hpiHitsDstofVert->Sumw2();
    TH1D *hproHitsDstofVert = new TH1D("hproHitsDstofVert", Form("%d blocks, position of S4 proton hits; Bar; Events / spill", nBlocks), 10, 0.5, 10.5);
    hproHitsDstofVert->Sumw2();
    TH1D *hproPiDstofVert = new TH1D(Form("hproPiDstofVert_block%d",nBlocks), "Proton/(#pi+#mu) ratio in S4; Vertical Bar in S4", 10, 0.5, 10.5);
    // Horizontal binning only
    TH1D *hpiHitsDstofHorz = new TH1D("hpiHitsDstofHorz", Form("%d blocks, position of S4 #pi hits; x / cm", nBlocks), 20, 0, 140);
    hpiHitsDstofHorz->Sumw2();
    TH1D *hproHitsDstofHorz = new TH1D("hproHitsDstofHorz", Form("%d blocks, position of S4 proton hits; x / cm", nBlocks), 20, 0, 140);
    hproHitsDstofHorz->Sumw2();
    TH1D *hproPiDstofHorz = new TH1D(Form("hproPiDstofHorz_block%d",nBlocks), "Proton/(#pi+#mu) ratio in S4; x / cm", 20, 0, 140);
    // With efficiencies taken into account
    // 1D hists
    TH1D *hpiHitsDstofVert_eff = new TH1D(Form("hpiHitsDstofVert_eff_block%d",nBlocks), Form("%d blocks, position of S4 #pi hits (eff. adjusted); Bar; Events / spill", nBlocks), 10, 0.5, 10.5);
    hpiHitsDstofVert_eff->Sumw2();
    TH1D *hproHitsDstofVert_eff = new TH1D(Form("hproHitsDstofVert_eff_block%d",nBlocks), Form("%d blocks, position of S4 proton hits (eff. adjusted); Bar; Events / spill", nBlocks), 10, 0.5, 10.5);
    hproHitsDstofVert_eff->Sumw2();
    TH1D *hproPiDstofVert_eff = new TH1D(Form("hproPiDstofVert_eff_block%d",nBlocks), "Proton/(#pi+#mu) ratio in S4 (eff. adjusted); Vertical Bar in S4", 10, 0.5, 10.5);
    // Horizontal binning only
    TH1D *hpiHitsDstofHorz_eff = new TH1D(Form("hpiHitsDstofHorz_eff_block%d",nBlocks), Form("%d blocks, position of S4 #pi hits (eff. adjusted); x / cm", nBlocks), 20, 0, 140);
    hpiHitsDstofHorz_eff->Sumw2();
    TH1D *hproHitsDstofHorz_eff = new TH1D(Form("hproHitsDstofHorz_eff_block%d",nBlocks), Form("%d blocks, position of S4 proton hits (eff. adjusted); x / cm", nBlocks), 20, 0, 140);
    hproHitsDstofHorz_eff->Sumw2();
    TH1D *hproPiDstofHorz_eff = new TH1D(Form("hproPiDstofHorz_eff_block%d",nBlocks), "Proton/(#pi+#mu) ratio in S4 (eff. adjusted); x / cm", 20, 0, 140);

    TH1D *hbarEff = new TH1D(Form("hbarEff_block%d", nBlocks), Form("%d blocks, S4 bar efficiency; Vertical bar in S4; Efficiency", nBlocks), 10, 0.5, 10.5);

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

    for (int irun=950; irun<1200; irun++){
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
    }
    
    cout << "Min and max runs are " << runMin << " " << runMax << endl;
    
    for (int itdc=0; itdc<2; itdc++) {
      TChain *tofCoinChain = new TChain("tofCoinTree");
      TChain *tofChain = new TChain("tofTree");

      for (int irun=runMin; irun<runMax+1; irun++){
	tofCoinChain->Add(Form("%srun%d/DsTOFcoincidenceRun%d_tdc%d.root", dstofDir, irun, irun, itdc+1));
	if (useEffs) {
	  tofChain->Add(Form("%srun%d/DsTOFtreeRun%d_tdc%d.root", dstofDir, irun, irun, itdc+1));
	}
      }
      
      RawDsTofCoincidence *tofCoin = NULL;
      tofCoinChain->SetBranchAddress("tofCoin", &tofCoin);
      RawDsTofHeader *tof = NULL;
      tofChain->SetBranchAddress("tof", &tof);

      for (int h=0; h<tofCoinChain->GetEntries(); h++) {
	tofCoinChain->GetEntry(h);
	if (tofCoin->unixTime[0]<startTime) continue;
	if (tofCoin->unixTime[0]>endTime) break;
	// Need to calculate total signal hits here so we 
	double deltat = TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1]  );
	double dstofHitT = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (10. - TMath::Abs(deltat) / 2. );
	double tofCalc = dstofHitT - tofCoin->usTofSignal;
	if (tofCalc < 150. && tofCalc > 50.) {
	  hCoins->Fill(tofCoin->bar); // Signal bar hits
	}
      }

      vector<double> ustofTimes;
      for (int jentry=0; jentry<tofChain->GetEntries(); jentry++) {
	tofChain->GetEntry(jentry);
	if (tof->unixTime<startTime) continue;
	if (tof->unixTime>endTime) break;
	
	if(tof->channel ==13) {
	  ustofTimes.push_back(tof->fakeTimeNs - 184.7);
	}
      }
      cout<<"Ustof signals "<<ustofTimes.size()<<endl;
      double dstofHitT, deltat;

      if (useEffs) {
	double lastUstof = 0.;
	cout<<"Matching ustof signals to pmt hits"<<endl;
	int lastu=0;
	for (int jentry=0; jentry<tofChain->GetEntries(); jentry++) {
	  if (jentry % 10000 == 0) {
	    cout<<jentry<<" of "<<tofChain->GetEntries()<<endl;
	  }
	  tofChain->GetEntry(jentry);
	  if (tof->unixTime<startTime) continue;
	  if (tof->unixTime>endTime) break;
	  // Check if each pmt hit is in coincidence with a utof hit
	  /*
	  if (tof->channel==13) {
	    lastUstof = tof->fakeTimeNs - 184.7;
	  }
	  */
	  if (tof->channel <= 10 && tof->channel > 0) {
	    for (int u = lastu; u < ustofTimes.size(); u++) {
	      double pmtTof = tof->fakeTimeNs - ustofTimes[u];
	      // Is a coincidence
	      if (pmtTof > 130 && pmtTof < 230) {
		// TDC 1
		if (itdc == 0) {
		  if (tof->channel % 2 == 1) { // PMT B
		    pmtHitsB[((tof->channel + 1) / -2) + 10]++;
		    pmtHitTimesB[((tof->channel + 1) / -2) + 10].push_back(tof->fakeTimeNs);
		    htof1dHitsB->Fill(pmtTof);
		    hHitsB->Fill(((tof->channel + 1) / -2) + 10);
		  }
		  else { // PMT A
		    pmtHitsA[((tof->channel) / -2) + 10]++;
		    pmtHitTimesA[((tof->channel) / -2) + 10].push_back(tof->fakeTimeNs);
		    htof1dHitsA->Fill(pmtTof);
		    hHitsA->Fill(((tof->channel) / -2) + 10);
		  }
		}
		// TDC 2
		else {
		  if (tof->channel % 2 == 1) { // PMT B
		    pmtHitsB[((tof->channel + 1) / -2) + 5]++;
		    pmtHitTimesB[((tof->channel + 1) / -2) + 5].push_back(tof->fakeTimeNs);
		    htof1dHitsB->Fill(pmtTof);
		    hHitsB->Fill(((tof->channel + 1) / -2) + 5);
		  }
		  else { // PMT A
		    pmtHitsA[((tof->channel) / -2) + 5]++;
		    pmtHitTimesA[((tof->channel) / -2) + 5].push_back(tof->fakeTimeNs);
		    htof1dHitsA->Fill(pmtTof);
		    hHitsA->Fill(((tof->channel) / -2) + 5);
		  }
		}
		lastu = u;
		break;
	      }
	    }
	  }
	  
	}
      } // useEffs
    } // itdc
  
    
    int lastb = 0;
    for (int bar = 0; bar < 10; bar++) {
      for (int a = 0; a < pmtHitTimesA[bar].size(); a++) {
	for (int b = 0; b < pmtHitTimesB[bar].size(); b++) {
	  if(abs(pmtHitTimesA[bar][a] - pmtHitTimesB[bar][b]) < 20.) {
	    pmtHitsA[bar]--;
	    lastb = b;
	    break;
	  }
	}
      }
      hHits->SetBinContent(bar+1, (pmtHitsA[bar]+pmtHitsB[bar]));
    } // bars
    
    hbarEff->Divide(hCoins, hHits, 1., 1., "B");

    for (int itdc=0; itdc<2; itdc++) {
      TChain *tof1CoinChain = new TChain("tofCoinTree");
      TChain *tof1Chain = new TChain("tofTree");

      for (int irun=runMin; irun<runMax+1; irun++){
	tof1CoinChain->Add(Form("%srun%d/DsTOFcoincidenceRun%d_tdc%d.root", dstofDir, irun, irun, itdc+1));
      }    

      RawDsTofCoincidence *tof1Coin = NULL;
      tof1CoinChain->SetBranchAddress("tofCoin", &tof1Coin);
     
      for (int ientry=0; ientry<tof1CoinChain->GetEntries(); ientry++){
	tof1CoinChain->GetEntry(ientry);
	if (tof1Coin->unixTime[0]<startTime) continue;
	if (tof1Coin->unixTime[0]>endTime) break;

	if (tof1Coin->lastDelayedBeamSignal != lastSpill && itdc == 0) {
	  lastSpill = tof1Coin->lastDelayedBeamSignal;
	  nSpills++;
	  for (int sp=ientry; sp<tof1CoinChain->GetEntries(); sp++) {
	    tof1CoinChain->GetEntry(sp);
	    if (tof1Coin->unixTime[0]<startTime) continue;
	    if (tof1Coin->unixTime[0]>endTime) break;

	    if ((tof1Coin->fakeTimeNs[0] - tof1Coin->usTofSignal) < 200. &&
		(tof1Coin->fakeTimeNs[0] - tof1Coin->usTofSignal) > 70. &&
		(tof1Coin->fakeTimeNs[0] - lastSpill) < 1e9 &&
		(tof1Coin->fakeTimeNs[0] - lastSpill) > 0) {
	      nSpillsTrue++;
	      break;
	    }

	  } // for (int sp=ientry; sp<tof1CoinChain->GetEntries(); sp++) 
	} // if (tof1Coin->lastDelayedBeamSignal != lastSpill && itdc == 0)
        
	double deltat1 = TMath::Abs(tof1Coin->fakeTimeNs[0]-tof1Coin->fakeTimeNs[1]  );
	double dstofHitT1 = min(tof1Coin->fakeTimeNs[0], tof1Coin->fakeTimeNs[1]) - (10. - TMath::Abs(deltat1) / 2. );
	double tofCalc = dstofHitT1 - tof1Coin->usTofSignal;

	// Is a genuine coincidence
	if (tofCalc < 150. && tofCalc > 50.) {
	  hCoins->Fill(tof1Coin->bar); // Signal bar hits
	  htof1d->Fill(dstofHitT1 - tof1Coin->usTofSignal);
	} // if (tof < 150. && tof > 50.)
	if (tofCalc > piLow && tofCalc < piHi) {
	  piHitsDstof->Fill((-1.*tof1Coin->fakeTimeNs[0] + tof1Coin->fakeTimeNs[1])*(7./2.)+70., (tof1Coin->bar*7.5) - 2.5);
	  piHitsDstofVert->Fill(1, tof1Coin->bar);
	  piHitsDstofHorz->Fill((-1.*tof1Coin->fakeTimeNs[0] + tof1Coin->fakeTimeNs[1])*(7./2.)+70., 1);
	  nPi++;
	  if (tof1Coin->bar != 10) {
	    hpiHitsDstofVert->Fill(tof1Coin->bar);
	    hpiHitsDstofVert_eff->Fill(tof1Coin->bar, 1./hbarEff->GetBinContent(tof1Coin->bar));
	    hpiHitsDstofHorz_eff->Fill((-1.*tof1Coin->fakeTimeNs[0] + tof1Coin->fakeTimeNs[1])*(7./2.)+70., 1./hbarEff->GetBinContent(tof1Coin->bar));
	    hpiHitsDstofHorz->Fill((-1.*tof1Coin->fakeTimeNs[0] + tof1Coin->fakeTimeNs[1])*(7./2.)+70.);
	  }
	}
	else if (tofCalc > proLow && tofCalc < proHi) {
	  proHitsDstof->Fill((-1.*tof1Coin->fakeTimeNs[0] + tof1Coin->fakeTimeNs[1])*(7./2.)+70., (tof1Coin->bar*7.5) - 2.5);
	  proHitsDstofVert->Fill(1, tof1Coin->bar);
	  proHitsDstofHorz->Fill((-1.*tof1Coin->fakeTimeNs[0] + tof1Coin->fakeTimeNs[1])*(7./2.)+70., 1);
	  nP++;
	  if (tof1Coin->bar != 10) {
	    hproHitsDstofHorz->Fill((-1.*tof1Coin->fakeTimeNs[0] + tof1Coin->fakeTimeNs[1])*(7./2.)+70.);
	    hproHitsDstofHorz_eff->Fill((-1.*tof1Coin->fakeTimeNs[0] + tof1Coin->fakeTimeNs[1])*(7./2.)+70., 1./hbarEff->GetBinContent(tof1Coin->bar));
	    hproHitsDstofVert->Fill(tof1Coin->bar);
	    hproHitsDstofVert_eff->Fill(tof1Coin->bar, 1./hbarEff->GetBinContent(tof1Coin->bar));
	  } // if (tof1Coin->bar != 10)
	} // else if (tof > proLow && tof < proHi)
      } // for (int ientry=0; ientry<tof1CoinChain->GetEntries(); ientry++)
    } // for (int itdc=0; itdc<2; itdc++)


    cout<<"Data from "<<nSpills<<" spills ("<<nSpillsTrue<<" true)"<<" over "<<(endTime - startTime)<<" seconds"<<endl;

    gStyle->SetOptStat(0);
    gStyle->SetPalette(55);
    
    TCanvas *cEff = new TCanvas(Form("cEff_%d",nBlocks));
    hbarEff->Draw("hist E");
    cEff->Print(Form("%s/nBlocksPlots/Eff_blocks%d.png", saveDir, nBlocks));
    cEff->Print(Form("%s/nBlocksPlots/Eff_blocks%d.pdf", saveDir, nBlocks));
    TCanvas *cCoins = new TCanvas(Form("cCoins_%d",nBlocks));
    hCoins->Draw("hist E");
    cCoins->Print(Form("%s/nBlocksPlots/coins_blocks%d.png", saveDir, nBlocks));
    cCoins->Print(Form("%s/nBlocksPlots/coins_blocks%d.pdf", saveDir, nBlocks));
    TCanvas *cHits = new TCanvas(Form("cHits_%d",nBlocks));
    hHits->Draw("hist E");
    cHits->Print(Form("%s/nBlocksPlots/Hits_blocks%d.png", saveDir, nBlocks));
    cHits->Print(Form("%s/nBlocksPlots/Hits_blocks%d.pdf", saveDir, nBlocks));
    TCanvas *cHitsA = new TCanvas(Form("cHitsA2_%d",nBlocks));
    hHitsA->Draw("hist E");
    cHitsA->Print(Form("%s/nBlocksPlots/HitsA_blocks%d.png", saveDir, nBlocks));
    cHitsA->Print(Form("%s/nBlocksPlots/HitsA_blocks%d.pdf", saveDir, nBlocks));
    TCanvas *cHitsB = new TCanvas(Form("cHitsB2_%d",nBlocks));
    hHitsB->Draw("hist E");
    cHitsB->Print(Form("%s/nBlocksPlots/HitsB_blocks%d.png", saveDir, nBlocks));
    cHitsB->Print(Form("%s/nBlocksPlots/HitsB_blocks%d.pdf", saveDir, nBlocks));

    TCanvas *ctofHitsA = new TCanvas(Form("cHitsA_%d",nBlocks));
    htof1dHitsA->Draw("hist");
    ctofHitsA->Print(Form("%s/nBlocksPlots/tofHitsA_blocks%d.png", saveDir, nBlocks));
    ctofHitsA->Print(Form("%s/nBlocksPlots/tofHitsA_blocks%d.pdf", saveDir, nBlocks));
    TCanvas *ctofHitsB = new TCanvas(Form("cHitsB_%d",nBlocks));
    htof1dHitsB->Draw("hist");
    ctofHitsB->Print(Form("%s/nBlocksPlots/tofHitsB_blocks%d.png", saveDir, nBlocks));
    ctofHitsB->Print(Form("%s/nBlocksPlots/tofHitsB_blocks%d.pdf", saveDir, nBlocks));

    TCanvas *c1 = new TCanvas("c1");
    proPiDstof->Divide(proHitsDstof, piHitsDstof);
    proPiDstof->GetZaxis()->SetRangeUser(0, 1.7);
    c1->SetRightMargin(0.13);
    proPiDstof->Draw("colz");
    c1->Print(Form("%s/nBlocksPlots/%dblocks_propiratio.png", saveDir, nBlocks));
    c1->Print(Form("%s/nBlocksPlots/%dblocks_propiratio.pdf", saveDir, nBlocks));
    TCanvas *c2 = new TCanvas("c2");
    c2->SetRightMargin(0.13);
    proHitsDstof->Scale(1./nSpillsTrue);
    proHitsDstof->GetZaxis()->SetRangeUser(0, 3.); 
    proHitsDstof->Draw("colz");
    c2->Print(Form("%s/nBlocksPlots/%dblocks_protons.png", saveDir, nBlocks));
    c2->Print(Form("%s/nBlocksPlots/%dblocks_protons.pdf", saveDir, nBlocks));
    TCanvas *c3 = new TCanvas("c3");
    c3->SetRightMargin(0.13);
    piHitsDstof->Scale(1./nSpillsTrue);
    piHitsDstof->GetZaxis()->SetRangeUser(0, 14.);
    piHitsDstof->Draw("colz");
    c3->Print(Form("%s/nBlocksPlots/%dblocks_pions.png", saveDir, nBlocks));
    c3->Print(Form("%s/nBlocksPlots/%dblocks_pions.pdf", saveDir, nBlocks));
    TCanvas *c4 = new TCanvas("c4");
    //    htof1d->Scale(1./nSpillsTrue);
    htof1d->Draw("hist");
    c4->Print(Form("%s/nBlocksPlots/%dblocks_tof.png", saveDir, nBlocks));
    c4->Print(Form("%s/nBlocksPlots/%dblocks_tof.pdf", saveDir, nBlocks));
    /*
    TCanvas *c1_1 = new TCanvas("c1_1");
    proPiDstofVert->Divide(proHitsDstofVert, piHitsDstofVert);
    proPiDstofVert->GetZaxis()->SetRangeUser(0, 1.7);
    c1_1->SetRightMargin(0.13);
    proPiDstofVert->Draw("colz text");
    c1_1->Print(Form("%s/nBlocksPlots/%dblock_propiratio_vert.png", saveDir, nBlocks));
    c1_1->Print(Form("%s/nBlocksPlots/%dblocks_propiratio_vert.pdf", saveDir, nBlocks));
    TCanvas *c1_2 = new TCanvas("c1_2");
    c1_2->SetRightMargin(0.13);
    proHitsDstofVert->Scale(1./nSpillsTrue);
    proHitsDstofVert->GetZaxis()->SetRangeUser(0, 10.); 
    proHitsDstofVert->Draw("colz text");
    c1_2->Print(Form("%s/nBlocksPlots/%dblocks_protons_vert.png", saveDir, nBlocks));
    c1_2->Print(Form("%s/nBlocksPlots/%dblocks_protons_vert.pdf", saveDir, nBlocks));
    TCanvas *c1_3 = new TCanvas("c1_3");
    c1_3->SetRightMargin(0.13);
    piHitsDstofVert->Scale(1./nSpillsTrue);
    piHitsDstofVert->GetZaxis()->SetRangeUser(0, 70.); 
    piHitsDstofVert->Draw("colz text");
    c1_3->Print(Form("%s/nBlocksPlots/%dblocks_pions_vert.png", saveDir, nBlocks));
    c1_3->Print(Form("%s/nBlocksPlots/%dblocks_pions_vert.pdf", saveDir, nBlocks));
    
    TCanvas *c2_1 = new TCanvas("c2_1");
    proPiDstofHorz->Divide(proHitsDstofHorz, piHitsDstofHorz);
    proPiDstofHorz->GetZaxis()->SetRangeUser(0, 0.8);
    c2_1->SetRightMargin(0.13);
    proPiDstofHorz->Draw("colz");
    c2_1->Print(Form("%s/nBlocksPlots/%dblocks_propiratio_horz.png", saveDir, nBlocks));
    c2_1->Print(Form("%s/nBlocksPlots/%dblocks_propiratio_horz.pdf", saveDir, nBlocks));
    TCanvas *c2_2 = new TCanvas("c2_2");
    c2_2->SetRightMargin(0.13);
    proHitsDstofHorz->Scale(1./nSpillsTrue);
    proHitsDstofHorz->Draw("colz");
    c2_2->Print(Form("%s/nBlocksPlots/%dblocks_protons_horz.png", saveDir, nBlocks));
    c2_2->Print(Form("%s/nBlocksPlots/%dblocks_protons_horz.pdf", saveDir, nBlocks));
    TCanvas *c2_3 = new TCanvas("c2_3");
    c2_3->SetRightMargin(0.13);
    piHitsDstofHorz->Scale(1./nSpillsTrue);
    piHitsDstofHorz->Draw("colz");
    c2_3->Print(Form("%s/nBlocksPlots/%dblocks_pions_horz.png", saveDir, nBlocks));
    c2_3->Print(Form("%s/nBlocksPlots/%dblocks_pions_horz.pdf", saveDir, nBlocks));
    */
    TCanvas *c3_1 = new TCanvas("c3_1");
    hproPiDstofVert->Divide(hproHitsDstofVert, hpiHitsDstofVert, 1., 1., "B");
    hproPiDstofVert->Draw("hist E");
    c3_1->Print(Form("%s/nBlocksPlots/%dblocks_propiratio_vert_hist.png", saveDir, nBlocks));
    c3_1->Print(Form("%s/nBlocksPlots/%dblocks_propiratio_vert_hist.pdf", saveDir, nBlocks));
    TCanvas *c3_2 = new TCanvas("c3_2");
    hproPiDstofHorz->Divide(hproHitsDstofHorz, hpiHitsDstofHorz, 1., 1., "B");
    hproPiDstofHorz->Draw("hist E");
    c3_2->Print(Form("%s/nBlocksPlots/%dblocks_propiratio_horz_hist.png", saveDir, nBlocks));
    c3_2->Print(Form("%s/nBlocksPlots/%dblocks_propiratio_horz_hist.pdf", saveDir, nBlocks));
    TCanvas *c3_3 = new TCanvas("c3_3");
    hproHitsDstofHorz->Scale(1./nSpillsTrue);
    hproHitsDstofHorz->Draw("hist");
    c3_3->Print(Form("%s/nBlocksPlots/%dblocks_protons_horz_hist.png", saveDir, nBlocks));
    c3_3->Print(Form("%s/nBlocksPlots/%dblocks_protons_horz_hist.pdf", saveDir, nBlocks));
    TCanvas *c3_4 = new TCanvas("c3_4");
    hpiHitsDstofHorz->Scale(1./nSpillsTrue);
    hpiHitsDstofHorz->Draw("hist");
    c3_4->Print(Form("%s/nBlocksPlots/%dblocks_pions_horz_hist.png", saveDir, nBlocks));
    c3_4->Print(Form("%s/nBlocksPlots/%dblocks_pions_horz_hist.pdf", saveDir, nBlocks));
    TCanvas *c3_5 = new TCanvas("c3_5");
    hproHitsDstofVert->Scale(1./nSpillsTrue);
    hproHitsDstofVert->Draw("hist");
    c3_5->Print(Form("%s/nBlocksPlots/%dblocks_protons_vert_hist.png", saveDir, nBlocks));
    c3_5->Print(Form("%s/nBlocksPlots/%dblocks_protons_vert_hist.pdf", saveDir, nBlocks));
    TCanvas *c3_6 = new TCanvas("c3_6");
    hpiHitsDstofVert->Scale(1./nSpillsTrue);
    hpiHitsDstofVert->Draw("hist");
    c3_6->Print(Form("%s/nBlocksPlots/%dblocks_pions_vert_hist.png", saveDir, nBlocks));
    c3_6->Print(Form("%s/nBlocksPlots/%dblocks_pions_vert_hist.pdf", saveDir, nBlocks));

    if (useEffs) {
      TCanvas *ce_1 = new TCanvas("ce_1");
      hproPiDstofVert_eff->Divide(hproHitsDstofVert_eff, hpiHitsDstofVert_eff, 1., 1., "B");
      hproPiDstofVert_eff->Draw("hist E");
      ce_1->Print(Form("%s/nBlocksPlots/%dblocks_propiratio_vert_eff_hist.png", saveDir, nBlocks));
      ce_1->Print(Form("%s/nBlocksPlots/%dblocks_propiratio_vert_eff_hist.pdf", saveDir, nBlocks));
      TCanvas *ce_2 = new TCanvas("ce_2");
      hproPiDstofHorz_eff->Divide(hproHitsDstofHorz_eff, hpiHitsDstofHorz_eff, 1., 1., "B");
      hproPiDstofHorz_eff->Draw("hist E");
      ce_2->Print(Form("%s/nBlocksPlots/%dblocks_propiratio_horz_eff_hist.png", saveDir, nBlocks));
      ce_2->Print(Form("%s/nBlocksPlots/%dblocks_propiratio_horz_eff_hist.pdf", saveDir, nBlocks));
      TCanvas *ce_3 = new TCanvas("ce_3");
      hproHitsDstofHorz_eff->Scale(1./nSpillsTrue);
      hproHitsDstofHorz_eff->Draw("hist");
      ce_3->Print(Form("%s/nBlocksPlots/%dblocks_protons_horz_eff_hist.png", saveDir, nBlocks));
      ce_3->Print(Form("%s/nBlocksPlots/%dblocks_protons_horz_eff_hist.pdf", saveDir, nBlocks));
      TCanvas *ce_4 = new TCanvas("ce_4");
      hpiHitsDstofHorz_eff->Scale(1./nSpillsTrue);
      hpiHitsDstofHorz_eff->Draw("hist");
      ce_4->Print(Form("%s/nBlocksPlots/%dblocks_pions_horz_eff_hist.png", saveDir, nBlocks));
      ce_4->Print(Form("%s/nBlocksPlots/%dblocks_pions_horz_eff_hist.pdf", saveDir, nBlocks));
      TCanvas *ce_5 = new TCanvas("ce_5");
      hproHitsDstofVert_eff->Scale(1./nSpillsTrue);
      hproHitsDstofVert_eff->Draw("hist");
      ce_5->Print(Form("%s/nBlocksPlots/%dblocks_protons_vert_eff_hist.png", saveDir, nBlocks));
      ce_5->Print(Form("%s/nBlocksPlots/%dblocks_protons_vert_eff_hist.pdf", saveDir, nBlocks));
      TCanvas *ce_6 = new TCanvas("ce_6");
      hpiHitsDstofVert_eff->Scale(1./nSpillsTrue);
      hpiHitsDstofVert_eff->Draw("hist");
      ce_6->Print(Form("%s/nBlocksPlots/%dblocks_pions_vert_eff_hist.png", saveDir, nBlocks));
      ce_6->Print(Form("%s/nBlocksPlots/%dblocks_pions_vert_eff_hist.pdf", saveDir, nBlocks));
    }
    /*
    TCanvas *c3_2 = new TCanvas("c3_2");
    c1_2->SetRightMargin(0.13);
    proHitsDstofVert->Scale(1./nSpillsTrue);
    proHitsDstofVert->GetZaxis()->SetRangeUser(0, 10.); 
    proHitsDstofVert->Draw("colz text");
    c1_2->Print(Form("%s/nBlocksPlots/%dblocks_protons_vert.png", saveDir, nBlocks));
    c1_2->Print(Form("%s/nBlocksPlots/%dblocks_protons_vert.pdf", saveDir, nBlocks));
    TCanvas *c3_3 = new TCanvas("c3_3");
    c1_3->SetRightMargin(0.13);
    piHitsDstofVert->Scale(1./nSpillsTrue);
    piHitsDstofVert->GetZaxis()->SetRangeUser(0, 70.); 
    piHitsDstofVert->Draw("colz text");
    c1_3->Print(Form("%s/nBlocksPlots/%dblocks_pions_vert.png", saveDir, nBlocks));
    c1_3->Print(Form("%s/nBlocksPlots/%dblocks_pions_vert.pdf", saveDir, nBlocks));
    */
    c_proPiVert->cd();
    hproPiDstofVert->Divide(hproHitsDstofVert, hpiHitsDstofVert, 1., 1., "B");
    hproPiDstofVert->SetTitle("Proton/(#pi+#mu) in S4; Bar number in S4; P/(#pi+#mu)"); 
    if (nBlocks == 0) {
      hproPiDstofVert->SetLineColor(kBlue);
      legVert->AddEntry(hproPiDstofVert, "0 blocks", "l");
      hproPiDstofVert->Draw("hist E");
    }
    else if (nBlocks == 1) {
      hproPiDstofVert->SetLineColor(kRed);
      legVert->AddEntry(hproPiDstofVert, "1 block", "l");
      hproPiDstofVert->Draw("hist same E");
    }
    else if (nBlocks == 2) {
      hproPiDstofVert->SetLineColor(kBlack);
      legVert->AddEntry(hproPiDstofVert, "2 blocks", "l");
      hproPiDstofVert->Draw("hist same E");
    }
    else {
      hproPiDstofVert->SetLineColor(kGreen+2);
      legVert->AddEntry(hproPiDstofVert, "3 blocks", "l");
      hproPiDstofVert->Draw("hist same E");
    }
    
    c_proPiHorz->cd();
    hproPiDstofHorz->Divide(hproHitsDstofHorz, hpiHitsDstofHorz, 1., 1., "B");
    hproPiDstofHorz->SetTitle("Proton/(#pi+#mu) in S4; x / cm; P/(#pi+#mu)"); 
    if (nBlocks == 0) {
      hproPiDstofHorz->SetLineColor(kBlue);
      legHorz->AddEntry(hproPiDstofHorz, "0 blocks", "l");
      hproPiDstofHorz->Draw("hist E");
    }
    else if (nBlocks == 1) {
      hproPiDstofHorz->SetLineColor(kRed);
      legHorz->AddEntry(hproPiDstofHorz, "1 block", "l");
      hproPiDstofHorz->Draw("hist same E");
    }
    else if (nBlocks == 2) {
      hproPiDstofHorz->SetLineColor(kBlack);
      legHorz->AddEntry(hproPiDstofHorz, "2 blocks", "l");
      hproPiDstofHorz->Draw("hist same E");
    }
    else {
      hproPiDstofHorz->SetLineColor(kGreen+2);
      legHorz->AddEntry(hproPiDstofHorz, "3 blocks", "l");
      hproPiDstofHorz->Draw("hist same E");
    }
    
    double piInt = piHitsDstof->Integral(1, 20, 1, 10);
    cout<<nBlocks<<" blocks: Have "<<nPi<<" pis and "<<nP<<" protons in "<<nSpillsTrue<<" spills"<<endl;
  }
  c_proPiVert->cd();
  c_proPiVert->SetTitle("Proton/(#pi+#mu) in S4"); 
  legVert->Draw();
  if (useEffs) {
    c_proPiVert->Print(Form("%s/nBlocksPlots/proPiVert_eff_hist.png", saveDir));
    c_proPiVert->Print(Form("%s/nBlocksPlots/proPiVert_eff_hist.pdf", saveDir));
  }
  else {
    c_proPiVert->Print(Form("%s/nBlocksPlots/proPiVert_hist.png", saveDir));
    c_proPiVert->Print(Form("%s/nBlocksPlots/proPiVert_hist.pdf", saveDir));
  }

  c_proPiHorz->cd();
  c_proPiHorz->SetTitle("Proton/(#pi+#mu) in S4"); 
  if (useEffs) {
    legHorz->Draw();
    c_proPiHorz->Print(Form("%s/nBlocksPlots/proPiHorz_eff_hist.png", saveDir));
    c_proPiHorz->Print(Form("%s/nBlocksPlots/proPiHorz_eff_hist.pdf", saveDir));
  }
  else {
    c_proPiHorz->Print(Form("%s/nBlocksPlots/proPiHorz_hist.png", saveDir));
    c_proPiHorz->Print(Form("%s/nBlocksPlots/proPiHorz_hist.pdf", saveDir));
  }
}

// Divide up Tof spectra by bar to see if there is a difference in background
void tofBar(const char* saveDir, const char* dstofDir="/scratch0/dbrailsf/temp/mylinktodtof/") {
  gROOT->SetBatch(kTRUE);

  // Unix timestamps for variable block moves
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
  const double end3Block   = 1535798437;
  // 0.8GeV/c, 4 block
  // Most runs were in this configuration so don't need to use necessarily
  const double start4Block = 1535608220;
  const double end4Block   = 1535617102;

  for (int nBlocks=0; nBlocks < 4; nBlocks++) {
    double nSpills = 0.;
    double nSpillsTrue = 0.;
    double lastSpill = 0.;

    TH2D *h2dtof = new TH2D(Form("h2dtof_%d",nBlocks), Form("S4 - S1 by vertical S4 bar: %d blocks; S4 - S1 / ns; Vertical bar in S4; Events / spill",nBlocks), 100, 50, 150, 10, 0.5, 10.5);
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

    for (int irun=950; irun<1200; irun++){
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
    }
  
    cout << "Min and max runs are " << runMin << " " << runMax << endl;

    for (int itdc=0; itdc<2; itdc++) {
      TChain *tofCoinChain = new TChain("tofCoinTree");

      for (int irun=runMin; irun<runMax+1; irun++){
	tofCoinChain->Add(Form("%srun%d/DsTOFcoincidenceRun%d_tdc%d.root", dstofDir, irun, irun, itdc+1));
      }

      RawDsTofCoincidence *tofCoin = NULL;
      tofCoinChain->SetBranchAddress("tofCoin", &tofCoin);

      double dstofHitT, deltat;

      for (int ientry=0; ientry<tofCoinChain->GetEntries(); ientry++){
	tofCoinChain->GetEntry(ientry);
	if (tofCoin->unixTime[0]<startTime) continue;
	if (tofCoin->unixTime[0]>endTime) break;
      
	if (tofCoin->lastDelayedBeamSignal != lastSpill && itdc==0) {
	  lastSpill = tofCoin->lastDelayedBeamSignal;
	  nSpills++;
	  for (int sp=ientry; sp < tofCoinChain->GetEntries(); sp++) {
	    tofCoinChain->GetEntry(sp);
	    if ((tofCoin->fakeTimeNs[0] - tofCoin->usTofSignal) < 200. && 
		(tofCoin->fakeTimeNs[0] - tofCoin->usTofSignal) > 70. &&
		(tofCoin->fakeTimeNs[0] - lastSpill) < 1e9 &&
		(tofCoin->fakeTimeNs[0] - lastSpill) > 0) {
	      nSpillsTrue++;
	      break;
	    }
	  
	    if (tofCoin->unixTime[0]>endTime) break;
	  }
	}      
      
	deltat = TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1]  );
	dstofHitT = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (10.-TMath::Abs(deltat) / 2 );
	if (tofCoin->inSpill) {
	  double tof = dstofHitT - tofCoin->usTofSignal;
	  h2dtof->Fill(tof, tofCoin->bar);  
	}
      }
    }
    gStyle->SetPalette(55);
    gStyle->SetOptStat(0);
    h2dtof->Scale(1./nSpillsTrue);
    TCanvas *c1 = new TCanvas(Form("c1_%d", nBlocks));
    c1->SetLogz();
    c1->SetRightMargin(0.13);
    h2dtof->Draw("colz");
    c1->Print(Form("%s/nBlocksPlots/2dtofVert_blocks%d.png", saveDir, nBlocks));
    c1->Print(Form("%s/nBlocksPlots/2dtofVert_blocks%d.pdf", saveDir, nBlocks));
  } // nBlocks
} // tofBar
