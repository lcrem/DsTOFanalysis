// tof1dComp.C
void tof1dComp (const char* saveDir, 
		const char* dstofDir="/scratch0/dbrailsf/temp/mylinktodtof/",
		const char* ustofDir="/home/sjones/mylinktoutof/") 
{
  //  gROOT->SetBatch(kTRUE);

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
  const double start4Block = 1535608220;
  const double end4Block   = 1535617102;

  TCanvas *cd = new TCanvas("cd");
  cd->SetLogy();
  TCanvas *cu = new TCanvas("cu");
  cu->SetLogy();
  THStack *hsu = new THStack("hsu", "Time of flight for various moderator thicknesses; S3 - S1 / ns; Events / spill");
  THStack *hsd = new THStack("hsd", "Time of flight for various moderator thicknesses; S4 - S1 / ns; Events / spill");

  TLegend *legd = new TLegend(0.75, 0.6, 0.87, 0.8);
  TLegend *legu = new TLegend(0.75, 0.6, 0.87, 0.8);
  // Find the appropriate dtof files
  for (int nBlocks=0; nBlocks < 4; nBlocks++) {

    int nSpills = 0;
    int nSpillsTrue = 0;
    double lastSpill = 0.;
    
    TH1D *hdtof1d = new TH1D(Form("hdtof1d_%d",nBlocks), Form("Time of flight, %d blocks; S4 - S1 / ns; Events / spill", nBlocks), 160, 70, 150);
    TH1D *hutof1d = new TH1D(Form("hutof1d_%d",nBlocks), Form("Time of flight, %d blocks; S3 - S1 / ns; Events / spill", nBlocks), 160, -40, 40);

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

    TH2D *hdtof2d = new TH2D(Form("hdtof2d_%d",nBlocks), Form("Time of flight across run length, %d blocks; S4 - S1 / ns; Unix time / s",nBlocks), 160, 70, 150, 400, startTime, endTime);

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
	if (tofCalc < 150. && tofCalc > 70.) {
	  hdtof1d->Fill(tofCalc);
	  hdtof2d->Fill(tofCalc, tofCoin->unixTime[0]);
	} // if (tofCalc < 150. && tofCalc > 50.) 
      } // for (int h=0; h<tofCoinChain->GetEntries(); h++) 
    } // for (int itdc=0; itdc<2; itdc++)

    if (nBlocks==0) {
      hdtof1d->SetLineColor(kRed);
      legd->AddEntry(hdtof1d, "0 blocks",  "l");
    }
    else if (nBlocks==1) {
      hdtof1d->SetLineColor(kBlue);
      legd->AddEntry(hdtof1d, "1 block",  "l");
    }
    else if (nBlocks==2) {
      hdtof1d->SetLineColor(kBlack);
      legd->AddEntry(hdtof1d, "2 blocks",  "l");
    }
    else if (nBlocks==3) {
      hdtof1d->SetLineColor(kGreen+2);
      legd->AddEntry(hdtof1d, "3 blocks",  "l");
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

    TFile *futof = new TFile(nustof, "read");

    double tToF[50];
    double tS1;
    int nhit;

    TTree *tree = (TTree*)futof->Get("tree");

    tree->SetBranchAddress("nhit", &nhit);
    tree->SetBranchAddress("tS1", &tS1);
    tree->SetBranchAddress("tToF", tToF);

    for (int t=0; t<tree->GetEntries(); t++) {
      tree->GetEntry(t);
      for (int n=0; n<nhit; n++) {
	double tofCalc = tToF[n] - tS1;
	hutof1d->Fill(tofCalc);
      }
    } // for (int t=0; t<tree->GetEntries(); t++)

    if (nBlocks==0) {
      hutof1d->SetLineColor(kRed);
      legu->AddEntry(hutof1d, "0 blocks",  "l");
    }
    else if (nBlocks==1) {
      hutof1d->SetLineColor(kBlue);
      legu->AddEntry(hutof1d, "1 block",  "l");
    }
    else if (nBlocks==2) {
      hutof1d->SetLineColor(kBlack);
      legu->AddEntry(hutof1d, "2 blocks",  "l");
    }
    else if (nBlocks==3) {
      hutof1d->SetLineColor(kGreen+2);
      legu->AddEntry(hutof1d, "3 blocks",  "l");
    }

    hutof1d->Scale(1./ (double)nSpillsTrue);

    hsu->Add(hutof1d);

  } // nBlocks

  cd->cd();
  hsd->Draw("hist nostack");
  legd->Draw();
  cd->Print(Form("%s/dstofAllBlockLog.png",saveDir));
  cd->Print(Form("%s/dstofAllBlockLog.pdf",saveDir));

  cu->cd();
  hsu->Draw("hist nostack");
  legu->Draw();
  cu->Print(Form("%s/ustofAllBlockLog.png", saveDir));
  cu->Print(Form("%s/ustofAllBlockLog.pdf", saveDir));
}
