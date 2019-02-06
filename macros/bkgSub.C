// bkgSub.C
// Performs background subtraction on the dtof spectrum

void bkgSub (const char* saveDir, 
	     const int firstBlock,
	     const int lastBlock,
	     const char* dstofDir="/scratch0/dbrailsf/temp/mylinktodtof/") 
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
  const double start4Block = 1536537600; 
  const double end4Block   = 1536669600;


  // Find the appropriate dtof files
  for (int nBlocks=firstBlock; nBlocks <= lastBlock; nBlocks++) {

    // Define signal and background functions to be fitted
    // Signals are gaussians
    TF1 *sPro = new TF1("sPro", "gaus", 106, 140);
    TF1 *sPi  = new TF1("sPi", "gaus", 75, 95);
    // Background is flat
    TF1 *fBkg    = new TF1("fBkg", "[0]", 70, 180);
    TF1 *fBkgLin = new TF1("fBkgLin","[0] + [1]*x", 70, 200);
    TF1 *fBkgExp = new TF1("fBkgExp","expo", 70, 200);
    fBkg->SetLineColor(kBlue);
    fBkgLin->SetLineColor(kBlue);
    sPro->SetLineColor(kGreen+2);
    sPi->SetLineColor(kRed);

    TF1 *fSplusB    = new TF1("signal+bkg", "gaus(0)+gaus(3)+[6]", 70, 180);
    TF1 *fSplusBLin = new TF1("signal+bkg lin", "gaus(0)+gaus(3)+[6]+[7]*x", 70, 200);
    TF1 *fSplusBExp = new TF1("signal+bkg lin", "gaus(0)+gaus(3)+expo(6)", 70, 200);
    fSplusB->SetParNames("const 1", "mean 1", "sigma 1",
			 "const 2", "mean 2", "sigma 2",
			 "bkg");
    fSplusBLin->SetParNames("const 1", "mean 1", "sigma 1",
			    "const 2", "mean 2", "sigma 2",
			    "bkgconst", "bkgslope");
    fSplusBExp->SetParNames("const 1", "mean 1", "sigma 1",
			    "const 2", "mean 2", "sigma 2",
			    "bkgconst", "bkgdecay");
    fSplusB->SetLineColor(kBlack);
    fSplusBLin->SetLineColor(kBlack);
    fSplusBExp->SetLineColor(kBlack);

    int nSpills = 0;
    int nSpillsTrue = 0;
    double lastSpill = 0.;
    
    TH1D *hdtof1d = new TH1D(Form("hdtof1d_%d",nBlocks), Form("Time of flight, %d blocks; S4 - S1 / ns; Events / spill", nBlocks), 260, 70, 200);

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
	} // if (tofCalc < 180. && tofCalc > 50.) 
      } // for (int h=0; h<tofCoinChain->GetEntries(); h++) 
    } // for (int itdc=0; itdc<2; itdc++)

    if (nBlocks==0) {
      hdtof1d->SetLineColor(kRed);
    }
    else if (nBlocks==1) {
      hdtof1d->SetLineColor(kBlue);
    }
    else if (nBlocks==2) {
      hdtof1d->SetLineColor(kBlack);
    }
    else if (nBlocks==3) {
      hdtof1d->SetLineColor(kGreen+2);
    }
    else if (nBlocks==4) {
      hdtof1d->SetLineColor(kMagenta);
    }

    hdtof1d->Scale(1./ (double)nSpillsTrue);

    TCanvas *c2d = new TCanvas(Form("%d_c2d",nBlocks));
    c2d->SetLogy();
    hdtof1d->Fit(sPi, "R");
    hdtof1d->Fit(sPro, "R");
    hdtof1d->Fit(fBkg, "R");

    Double_t par[7];
    sPro->GetParameters(&par[0]);
    sPi->GetParameters(&par[3]);
    fBkg->GetParameters(&par[6]);
    fSplusB->SetParameters(par);
    hdtof1d->Fit(fSplusB, "R");

    hdtof1d->Draw("hist");
    sPro->Draw("same");
    sPi->Draw("same");
    fBkg->Draw("same");
    fSplusB->Draw("same");

    c2d->Print(Form("%s/%d_dtofFitted.png",saveDir,nBlocks));
    c2d->Print(Form("%s/%d_dtofFitted.pdf",saveDir,nBlocks));

    TCanvas *c2d_lin = new TCanvas(Form("%d_c2d_lin",nBlocks));
    c2d_lin->SetLogy();
    hdtof1d->Fit(fBkgLin, "R");
    Double_t parLin[8];
    sPro->GetParameters(&parLin[0]);
    sPi->GetParameters(&parLin[3]);
    fBkgLin->GetParameters(&parLin[6]);
    fSplusBLin->SetParameters(parLin);
    hdtof1d->Fit(fSplusBLin, "R");

    hdtof1d->Draw("hist");
    sPro->Draw("same");
    sPi->Draw("same");
    fBkgLin->Draw("same");
    fSplusBLin->Draw("same");

    c2d_lin->Print(Form("%s/%d_dtofFittedLinear.png",saveDir,nBlocks));
    c2d_lin->Print(Form("%s/%d_dtofFittedLinear.pdf",saveDir,nBlocks));

    TCanvas *c2d_exp = new TCanvas(Form("%d_c2d_exp",nBlocks));
    c2d_exp->SetLogy();
    hdtof1d->Fit(fBkgExp, "R");
    Double_t parExp[8];
    sPro->GetParameters(&parExp[0]);
    sPi->GetParameters(&parExp[3]);
    fBkgExp->GetParameters(&parExp[6]);
    fSplusBExp->SetParameters(parExp);
    hdtof1d->Fit(fSplusBExp, "R");

    hdtof1d->Draw("hist");
    //    sPro->Draw("same");
    //    sPi->Draw("same");
    //    fBkgLin->Draw("same");
    fSplusBExp->Draw("same");

    c2d_exp->Print(Form("%s/%d_dtofFittedExp.png",saveDir,nBlocks));
    c2d_exp->Print(Form("%s/%d_dtofFittedExp.pdf",saveDir,nBlocks));

  } // nBlocks

} // bkgSub


