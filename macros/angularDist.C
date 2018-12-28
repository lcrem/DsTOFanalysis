// angularDist.C
// Angular distribution of protons and pions for different moderator blocks
void angularDist (const int nBlocks, const char* dstofDir="/scratch0/dbrailsf/temp/mylinktodtof/") {

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
  // Dstof directory
  //  const char* dsDir = "/scratch0/dbrailsf/temp/mylinktodtof/";
  // Timing cuts
  const double piLow  = 80.;
  const double piHi   = 95.;
  const double proLow = 106.;
  const double proHi  = 134.;

  double nSpills = 0.;
  double nSpillsTrue = 0.;
  double lastSpill = 0.;

  TH1D *htof1d = new TH1D("htof1d", Form("Time of flight, %d blocks; DsToF - UsToF / ns; Events / spill", nBlocks), 100, 50, 150);
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
      dstofHitT = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (10. - TMath::Abs(deltat) / 2 );
      double tof = dstofHitT - tofCoin->usTofSignal;
      htof1d->Fill(dstofHitT - tofCoin->usTofSignal);
      if (tof > piLow && tof < piHi) {
	piHitsDstof->Fill((-1.*tofCoin->fakeTimeNs[0] + tofCoin->fakeTimeNs[1])*(7./2.)+70., (tofCoin->bar*7.5) - 2.5);
	piHitsDstofVert->Fill(1, tofCoin->bar);
	piHitsDstofHorz->Fill((-1.*tofCoin->fakeTimeNs[0] + tofCoin->fakeTimeNs[1])*(7./2.)+70., 1);
      }
      else if (tof > proLow && tof < proHi) {
	proHitsDstof->Fill((-1.*tofCoin->fakeTimeNs[0] + tofCoin->fakeTimeNs[1])*(7./2.)+70., (tofCoin->bar*7.5) - 2.5);
	proHitsDstofVert->Fill(1, tofCoin->bar);
	proHitsDstofHorz->Fill((-1.*tofCoin->fakeTimeNs[0] + tofCoin->fakeTimeNs[1])*(7./2.)+70., 1);
      }
    }
  }

  cout<<"Data from "<<nSpills<<" spills ("<<nSpillsTrue<<" true)"<<" over "<<(endTime - startTime)<<" seconds"<<endl;

  gStyle->SetOptStat(0);
  gStyle->SetPalette(55);
  TCanvas *c1 = new TCanvas("c1");
  proPiDstof->Divide(proHitsDstof, piHitsDstof);
  proPiDstof->GetZaxis()->SetRangeUser(0, 1.7);
  c1->SetRightMargin(0.13);
  proPiDstof->Draw("colz");
  c1->Print(Form("../nBlocksPlots/%dblocks_propiratio.png", nBlocks));
  c1->Print(Form("../nBlocksPlots/%dblocks_propiratio.pdf", nBlocks));
  TCanvas *c2 = new TCanvas("c2");
  c2->SetRightMargin(0.13);
  proHitsDstof->Scale(1./nSpillsTrue);
  proHitsDstof->GetZaxis()->SetRangeUser(0, 3.); 
  proHitsDstof->Draw("colz");
  c2->Print(Form("../nBlocksPlots/%dblocks_protons.png", nBlocks));
  c2->Print(Form("../nBlocksPlots/%dblocks_protons.pdf", nBlocks));
  TCanvas *c3 = new TCanvas("c3");
  c3->SetRightMargin(0.13);
  piHitsDstof->Scale(1./nSpillsTrue);
  piHitsDstof->GetZaxis()->SetRangeUser(0, 14.);
  piHitsDstof->Draw("colz");
  c3->Print(Form("../nBlocksPlots/%dblocks_pions.png", nBlocks));
  c3->Print(Form("../nBlocksPlots/%dblocks_pions.pdf", nBlocks));
  TCanvas *c4 = new TCanvas("c4");
  htof1d->Scale(1./nSpillsTrue);
  htof1d->Draw("hist");
  c4->Print(Form("../nBlocksPlots/%dblocks_tof.png", nBlocks));
  c4->Print(Form("../nBlocksPlots/%dblocks_tof.pdf", nBlocks));

  TCanvas *c1_1 = new TCanvas("c1_1");
  proPiDstofVert->Divide(proHitsDstofVert, piHitsDstofVert);
  proPiDstofVert->GetZaxis()->SetRangeUser(0, 1.7);
  c1_1->SetRightMargin(0.13);
  proPiDstofVert->Draw("colz text");
  c1_1->Print(Form("../nBlocksPlots/%dblocks_propiratio_vert.png", nBlocks));
  c1_1->Print(Form("../nBlocksPlots/%dblocks_propiratio_vert.pdf", nBlocks));
  TCanvas *c1_2 = new TCanvas("c1_2");
  c1_2->SetRightMargin(0.13);
  proHitsDstofVert->Scale(1./nSpillsTrue);
  proHitsDstofVert->GetZaxis()->SetRangeUser(0, 10.); 
  proHitsDstofVert->Draw("colz text");
  c1_2->Print(Form("../nBlocksPlots/%dblocks_protons_vert.png", nBlocks));
  c1_2->Print(Form("../nBlocksPlots/%dblocks_protons_vert.pdf", nBlocks));
  TCanvas *c1_3 = new TCanvas("c1_3");
  c1_3->SetRightMargin(0.13);
  piHitsDstofVert->Scale(1./nSpillsTrue);
  piHitsDstofVert->GetZaxis()->SetRangeUser(0, 70.); 
  piHitsDstofVert->Draw("colz text");
  c1_3->Print(Form("../nBlocksPlots/%dblocks_pions_vert.png", nBlocks));
  c1_3->Print(Form("../nBlocksPlots/%dblocks_pions_vert.pdf", nBlocks));

  TCanvas *c2_1 = new TCanvas("c2_1");
  proPiDstofHorz->Divide(proHitsDstofHorz, piHitsDstofHorz);
  proPiDstofHorz->GetZaxis()->SetRangeUser(0, 0.8);
  c2_1->SetRightMargin(0.13);
  proPiDstofHorz->Draw("colz");
  c2_1->Print(Form("../nBlocksPlots/%dblocks_propiratio_horz.png", nBlocks));
  c2_1->Print(Form("../nBlocksPlots/%dblocks_propiratio_horz.pdf", nBlocks));
  TCanvas *c2_2 = new TCanvas("c2_2");
  c2_2->SetRightMargin(0.13);
  proHitsDstofHorz->Scale(1./nSpillsTrue);
  proHitsDstofHorz->Draw("colz");
  c2_2->Print(Form("../nBlocksPlots/%dblocks_protons_horz.png", nBlocks));
  c2_2->Print(Form("../nBlocksPlots/%dblocks_protons_horz.pdf", nBlocks));
  TCanvas *c2_3 = new TCanvas("c2_3");
  c2_3->SetRightMargin(0.13);
  piHitsDstofHorz->Scale(1./nSpillsTrue);
  piHitsDstofHorz->Draw("colz");
  c2_3->Print(Form("../nBlocksPlots/%dblocks_pions_horz.png", nBlocks));
  c2_3->Print(Form("../nBlocksPlots/%dblocks_pions_horz.pdf", nBlocks));

  double piInt = piHitsDstof->Integral(1, 20, 1, 10);
  cout<<"Have "<<piInt<<" pis per spill"<<endl;
}
