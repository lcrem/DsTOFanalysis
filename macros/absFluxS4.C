// asbFluxS4.C
// Calculate the absolute flux of particles in S4
// with no regard for if there was an upstream hit
void absFluxS4 (const char* saveDir, 
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
  // 4 moderator blocks with -4cm bend
  const double start4Block = 1535836129;
  const double end4Block   = 1535849634;
  // Long sample of 4 moderator block with +5cm bend 
  //  const double start4Block = 1536537600; 
  //  const double end4Block   = 1536669600;
  // S4 positions
  // S1 -> S4 baseline length
  const double baselineS1S4End   = 13.9426;
  const double baselineS1S4Start = 14.0069;
  const double s4OffAxisStartX = 0.121;
  const double s4OffAxisEndX   = 1.4224;
  // Shift in ns required to to pion peak at speed of light
  const double dstofShift = 40.;

  TFile *fout = new TFile(Form("%s/absFluxS4Plots.root", saveDir), "recreate");
  THStack *hs = new THStack("hs","Absolute particle flux in S4; x / m; Events / spill");
  THStack *hsXPrime = new THStack("hsXPrime","Absolute particle flux in S4; x / m; Events / spill");
  THStack *hsAngle = new THStack("hs","Absolute particle flux in S4; #theta / degrees; Events / spill");
  TLegend *legHorz = new TLegend(0.6, 0.8, 0.85, 0.6);
  
  for (int nBlocks = 0; nBlocks <= 4; nBlocks++) {
    cout<<"Block "<<nBlocks<<endl;
    int nSpills = 0;
    int nSpillsTrue = 0;
    double lastSpill = 0.;

    TH1D *habsFluxX = new TH1D(Form("habsFluxX%d",nBlocks), Form("Absolute particle flux in S4, %d blocks; x / m; Events / spill", nBlocks), 20, 0., 1.4);
    //    TH1D *habsFluxXPrime = new TH1D(Form("habsFluxXPrime%d",nBlocks), Form("Absolute particle flux in S4, %d blocks; x' / m; Events / spill", nBlocks), 20, 0., 1.4);
    TH1D *habsFluxXAngle = new TH1D(Form("habsFluxXAngle%d",nBlocks), Form("Absolute particle flux in S4, %d blocks; #theta / degrees; Events / spill", nBlocks), 40, 0, 6.);

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
      delete tofCoinTemp;
    } // for (int irun=950; irun<1400; irun++)
    
    cout << "Min and max runs are " << runMin << " " << runMax << endl;
    

    for (int itdc=0; itdc<2; itdc++) {
      for (int irun=runMin; irun<runMax+1; irun++){
	TChain *tofCoinChain = new TChain("tofCoinTree");
	for (int irun=runMin; irun<runMax+1; irun++){
	  tofCoinChain->Add(Form("%srun%d/DsTOFcoincidenceRun%d_tdc%d.root", dstofDir, irun, irun, itdc+1));
	  
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
	    if (tofCoin->inSpill) {
	      double deltat = TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1]  );
	      double dstofHitT = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (10. - TMath::Abs(deltat) / 2. );
	      // For true x, y position (relative to beam axis) interpolate between two 
	      // measured positions of the ToF (see survey data)
	      double positionX = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 65.) / 130.) * (s4OffAxisEndX - s4OffAxisStartX) + s4OffAxisStartX;
	      double positionY = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 65.) / 130.) * (baselineS1S4End - baselineS1S4Start) + baselineS1S4Start;
	      habsFluxX->Fill(positionX);
	      //	      habsFluxXPrime->Fill(((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 65.) / 130.);
	      double angleOffAxis = TMath::ATan(positionX / positionY) * 180. / TMath::Pi();
	      habsFluxXAngle->Fill(angleOffAxis);
	    } // if (tofCoin->inSpill)
	  } // for (int h=0; h<tofCoinChain->GetEntries(); h++) 
	} // for (int irun=runMin; irun<runMax+1; irun++)
      } // for (int irun=runMin; irun<runMax+1; irun++)
    } // for (int itdc=0; itdc<2; itdc++)
    
    // Save all of these flux plots individually
    fout->cd();
    TCanvas *c1 = new TCanvas("c1");
    habsFluxX->Scale(1. / nSpillsTrue);
    habsFluxX->Draw("hist");
    habsFluxX->Write();
    c1->Print(Form("%s/%d_absFluxX.png", saveDir, nBlocks));
    c1->Print(Form("%s/%d_absFluxX.pdf", saveDir, nBlocks));

    TCanvas *c2 = new TCanvas("c2");
    habsFluxXAngle->Scale(1. / nSpillsTrue);
    habsFluxXAngle->Draw("hist");
    habsFluxXAngle->Write();
    c2->Print(Form("%s/%d_absFluxXAngle.png", saveDir, nBlocks));
    c2->Print(Form("%s/%d_absFluxXAngle.pdf", saveDir, nBlocks));
    // And do them all together in different colours
    if (nBlocks == 0) {
      habsFluxX->SetLineColor(kBlue);
      habsFluxXAngle->SetLineColor(kBlue);
      legHorz->AddEntry(habsFluxX, "0 blocks", "l");
      hs->Add(habsFluxX);
    }
    else if (nBlocks == 1) {
      habsFluxX->SetLineColor(kRed);
      habsFluxXAngle->SetLineColor(kRed);
      legHorz->AddEntry(habsFluxX, "1 block", "l");
      hs->Add(habsFluxX);
    }
    else if (nBlocks == 2) {
      habsFluxX->SetLineColor(kBlack);
      habsFluxXAngle->SetLineColor(kBlack);
      legHorz->AddEntry(habsFluxX, "2 blocks", "l");
      hs->Add(habsFluxX);
    }
    else if (nBlocks == 3){
      habsFluxX->SetLineColor(kGreen+2);
      habsFluxXAngle->SetLineColor(kGreen+2);
      legHorz->AddEntry(habsFluxX, "3 blocks", "l");
      hs->Add(habsFluxX);
    }
    else {
      habsFluxX->SetLineColor(kMagenta);
      habsFluxXAngle->SetLineColor(kMagenta);
      legHorz->AddEntry(habsFluxX, "4 blocks", "l");
      hs->Add(habsFluxX);
    }
    hsAngle->Add(habsFluxXAngle);
    hsXPrime->Add(habsFluxXPrime);
} // nBlocks
  TCanvas *cs = new TCanvas("cs");
  hs->Draw("hist nostack");
  legHorz->Draw();
  hs->Write();
  cs->Print(Form("%s/absFluxX.png", saveDir));
  cs->Print(Form("%s/absFluxX.pdf", saveDir));
  TCanvas *csAngle = new TCanvas("csAngle");
  hsAngle->Draw("hist nostack");
  legHorz->Draw();
  hsAngle->Write();
  csAngle->Print(Form("%s/absFluxXAngle.png", saveDir));
  csAngle->Print(Form("%s/absFluxXAngle.pdf", saveDir));
  /*
  TCanvas *cXPrime = new TCanvas("cXPrime");
  hsXPrime->Draw("hist nostack");
  cXPrime->Print(Form("%s/absFluxXPrime.png", saveDir));
  cXPrime->Print(Form("%s/absFluxXPrime.pdf", saveDir));
  */
} // absFluxS4
