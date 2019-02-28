// tof1dComp.C
void tof1dComp (const char* saveDir, 
		const char* dstofDir="/scratch0/dbrailsf/temp/mylinktodtof/",
		const char* ustofDir="/home/sjones/mylinktoutof/") 
{
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
  const double end3Block   = 1535795300;
  // 0.8GeV/c, 4 block
  // Most runs were in this configuration so don't need to use necessarily
  const double start4Block = 1535836129;
  const double end4Block   = 1535879634;
  // Time of flight cuts for S1 to S3
  // Particles travelling at c should cross the distance in between about 36.9ns
  // and 35.6ns
  // This is before the shift is applied
  const double proLow  = -15.;
  const double proHi   = 42.5;
  const double piLow = -32.5;
  const double piHi  = -22.5;

  TCanvas *cd = new TCanvas("cd");
  cd->SetLogy();
  TCanvas *cu = new TCanvas("cu");
  cu->SetLogy();
  TCanvas *cuS12 = new TCanvas("cuS12");
  cuS12->SetLogy();
  THStack *hsuS1 = new THStack("hsuS1", "Time of flight for various moderator thicknesses (S1 trigger only); S3 - S1 / ns; Events / spill");
  THStack *hsuS12 = new THStack("hsuS12", "Time of flight for various moderator thicknesses (S1 & S2 trigger); S3 - S1 / ns; Events / spill");
  THStack *hsd = new THStack("hsd", "Time of flight for various moderator thicknesses; S4 - S1 / ns; Events / spill");

  TLegend *legd = new TLegend(0.75, 0.6, 0.87, 0.8);
  TLegend *legu = new TLegend(0.75, 0.6, 0.87, 0.8);
  // Find the appropriate dtof files
  for (int nBlocks=0; nBlocks <= 3; nBlocks++) {

    // THStack for the utof comparison with and without the S2 trigger
    THStack *hsutofComp = new THStack(Form("hsutofComp%d",nBlocks), Form("Time of flight as measured in S3 with and without S2 trigger, %d blocks; S3 - S1 / ns; Events / spill",nBlocks));

    int nSpills = 0;
    int nSpillsTrue = 0;
    double lastSpill = 0.;
    
    TH1D *hdtof1d = new TH1D(Form("hdtof1d_%d",nBlocks), Form("Time of flight, %d blocks; S4 - S1 / ns; Events / spill", nBlocks), 260, 70, 200);
    TH1D *hutof1dS1 = new TH1D(Form("hutof1dS1_%d",nBlocks), Form("Time of flight, %d blocks (S1 trigger only); S3 - S1 / ns; Events / spill", nBlocks), 260, -40, 90);
    TH1D *hutof1dS12 = new TH1D(Form("hutof1dS12_%d",nBlocks), Form("Time of flight, %d blocks (S1 & S2 trigger); S3 - S1 / ns; Events / spill", nBlocks), 260, -40, 90);

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

    TH2D *hdtof2d = new TH2D(Form("hdtof2d_%d",nBlocks), Form("Time of flight across run length, %d blocks; S4 - S1 / ns; Unix time / s",nBlocks), 260, 70, 200, 400, startTime, endTime);

    for (int irun=950; irun<1400; irun++){
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
      delete tofCoinTemp;
    }
    
    cout << "Min and max runs are " << runMin << " " << runMax << endl;

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

	  } // for (int sp=ientry; sp<tof1CoinChain->GetEntries(); sp++) 
	} // if (tof1Coin->lastDelayedBeamSignal != lastSpill && itdc == 0)
	// Need to calculate total signal hits here so we 
	double deltat = TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1]  );
	double dstofHitT = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (10. - TMath::Abs(deltat) / 2. );
	double tofCalc = dstofHitT - tofCoin->usTofSignal;
	if (tofCalc < 200. && tofCalc > 70.) {
	  hdtof1d->Fill(tofCalc);
	  hdtof2d->Fill(tofCalc, tofCoin->unixTime[0]);
	} // if (tofCalc < 200. && tofCalc > 70.) 
      } // for (int h=0; h<tofCoinChain->GetEntries(); h++) 
      delete tofCoin;
    } // for (int itdc=0; itdc<2; itdc++)

    hdtof1d->SetLineWidth(2);
    if (nBlocks==0) {
      hdtof1d->SetLineColor(kBlack);
      legd->AddEntry(hdtof1d, "0 blocks",  "l");
    }
    else if (nBlocks==1) {
      hdtof1d->SetLineColor(kRed);
      legd->AddEntry(hdtof1d, "1 block",  "l");
    }
    else if (nBlocks==2) {
      hdtof1d->SetLineColor(kBlue);
      legd->AddEntry(hdtof1d, "2 blocks",  "l");
    }
    else if (nBlocks==3) {
      hdtof1d->SetLineColor(kCyan+1);
      legd->AddEntry(hdtof1d, "3 blocks",  "l");
    }
    else if (nBlocks==4) {
      hdtof1d->SetLineColor(kOrange+1);
      legd->AddEntry(hdtof1d, "4 blocks",  "l");
    }

    hdtof1d->Scale(1./ (double)nSpillsTrue);
    hsd->Add(hdtof1d);

    TCanvas *c2d = new TCanvas("c2d");
    hdtof2d->Draw();
    c2d->Print(Form("%s/%d_dtof2d.png",saveDir,nBlocks));
    c2d->Print(Form("%s/%d_dtof2d.pdf",saveDir,nBlocks));

    // Now do the same thing for the utof
    // Input appropriate files by hand
    char *nustof;
    if (nBlocks==0) nustof = Form("%sData_2018_8_31_b2_800MeV_0block.root", ustofDir);
    else if (nBlocks==1) nustof = Form("%sData_2018_9_1_b4_800MeV_1block_bend4cm.root", ustofDir);
    else if (nBlocks==2) nustof = Form("%sData_2018_9_1_b2_800MeV_2block_bend4cm.root", ustofDir);
    else if (nBlocks==3) nustof = Form("%sData_2018_9_1_b3_800MeV_3block_bend4cm.root", ustofDir);
    else if (nBlocks==4) nustof = Form("%sData_2018_9_1_b8_800MeV_4block_bend4cm.root", ustofDir);

    TFile *futof = new TFile(nustof, "read");

    double tToF[50];
    double tTrig;
    double tS1;
    int nhit;

    TTree *tree = (TTree*)futof->Get("tree");

    tree->SetBranchAddress("nhit", &nhit);
    tree->SetBranchAddress("tS1", &tS1);
    tree->SetBranchAddress("tToF", tToF);
    tree->SetBranchAddress("tTrig", &tTrig);

    for (int t=0; t<tree->GetEntries(); t++) {
      tree->GetEntry(t);
      // Only select those events with multiplicity of 1 - advice from AK
      if (nhit == 1) {
	double tofCalc = tToF[0] - tS1;
	// S1 trigger only
	if (tTrig == 0) {
	  hutof1dS1->Fill(tofCalc);
	}
	// S1 & S2 trigger
	else {
	  hutof1dS12->Fill(tofCalc);
	}
      }
    } // for (int t=0; t<tree->GetEntries(); t++)

    hutof1dS1->SetLineWidth(2);
    hutof1dS12->SetLineWidth(2);
    if (nBlocks==0) {
      hutof1dS1->SetLineColor(kBlack);
      hutof1dS12->SetLineColor(kBlack);
      legu->AddEntry(hutof1dS1, "0 blocks",  "l");
    }
    else if (nBlocks==1) {
      hutof1dS1->SetLineColor(kRed);
      hutof1dS12->SetLineColor(kRed);
      legu->AddEntry(hutof1dS1, "1 block",  "l");
    }
    else if (nBlocks==2) {
      hutof1dS1->SetLineColor(kBlue);
      hutof1dS12->SetLineColor(kBlue);
      legu->AddEntry(hutof1dS1, "2 blocks",  "l");
    }
    else if (nBlocks==3) {
      hutof1dS1->SetLineColor(kCyan+1);
      hutof1dS12->SetLineColor(kCyan+1);
      legu->AddEntry(hutof1dS1, "3 blocks",  "l");
    }
    else if (nBlocks==4) {
      hutof1dS1->SetLineColor(kOrange+1);
      hutof1dS12->SetLineColor(kOrange+1);
      legu->AddEntry(hutof1dS1, "4 blocks",  "l");
    }

    hutof1dS1->Scale(1./ (double)nSpillsTrue);
    hutof1dS12->Scale(1./ (double)nSpillsTrue);

    TF1 *sPro = new TF1("sPro", "gaus", proLow, proHi);
    TF1 *sPi  = new TF1("sPi", "gaus", piLow, piHi);
    TCanvas *ctmp = new TCanvas(Form("ctmp%d",nBlocks));
    ctmp->SetLogy();
    hutof1dS1->Draw("hist");
    hutof1dS1->Fit(sPro, "R");
    hutof1dS1->Fit(sPi, "R");
    sPro->Draw("same");
    sPi->Draw("same");
    ctmp->Print(Form("%s/%d_utofS1.png",saveDir,nBlocks));
    ctmp->Print(Form("%s/%d_utofS1.pdf",saveDir,nBlocks));

    hsuS1->Add(hutof1dS1);
    hsuS12->Add(hutof1dS12);

    hutof1dS1->SetLineStyle(7);
    hsutofComp->Add(hutof1dS1);
    hsutofComp->Add(hutof1dS12);

    TCanvas *cutofComp = new TCanvas(Form("cutofComp%d",nBlocks));
    cutofComp->SetLogy();
    hsutofComp->Draw("hist nostack");
    cutofComp->Print(Form("%s/%d_utofComp.png", saveDir, nBlocks));
    cutofComp->Print(Form("%s/%d_utofComp.pdf", saveDir, nBlocks));

  } // nBlocks

  cd->cd();
  hsd->Draw("hist nostack");
  legd->Draw();
  cd->Print(Form("%s/dstofAllBlockLog.png",saveDir));
  cd->Print(Form("%s/dstofAllBlockLog.pdf",saveDir));
  cd->Print(Form("%s/dstofAllBlockLog.tex",saveDir));

  cu->cd();
  hsuS1->Draw("hist nostack");
  legu->Draw();
  cu->Print(Form("%s/ustofAllBlockLogS1.png", saveDir));
  cu->Print(Form("%s/ustofAllBlockLogS1.pdf", saveDir));
  cu->Print(Form("%s/ustofAllBlockLogS1.tex", saveDir));

  cuS12->cd();
  hsuS12->Draw("hist nostack");
  legu->Draw();
  cuS12->Print(Form("%s/ustofAllBlockLogS12.png", saveDir));
  cuS12->Print(Form("%s/ustofAllBlockLogS12.pdf", saveDir));
  cuS12->Print(Form("%s/ustofAllBlockLogS12.tex", saveDir));

}
