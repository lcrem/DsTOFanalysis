// proPiFracUnix.C
// Same as proPiFrac but takes unix timestamps as inputs

void proPiFracUnix(const int unixStart, const int unixEnd, 
		   const char* saveDir,
		   const char* dstofDir="/scratch0/dbrailsf/temp/mylinktodtof/")
{

  THStack *hsProPi = new THStack("hsProPi", "Proton and MIP profiles across spill as measured in S4; Time since beam signal / s; Events / spill");

  // Timing cuts
  const double piLow  = 80.;
  const double piHi   = 95.;
  const double proLow = 106.;
  const double proHi  = 134.;

  double nSpills = 0.;
  double nSpillsTrue = 0.;
  double lastSpill = 0.;

  TH1D *hproton = new TH1D("hproton", "Number of protons; Time since spill start / ns; Events / spill", 40, 0, 1);
  TH1D *hpion = new TH1D("hpion", "Number of MIPs; Time since spill start / ns; Events / spill", 40, 0, 1);
  hproton->Sumw2();
  hpion->Sumw2();
  TH1D *hratio = new TH1D("hratio", "Proton/MIP ratio across spill; Time since beam signal / s; Proton/MIP", 40, 0, 1);


  int startTime = unixStart;
  int endTime = unixEnd;

  // Find the correct dstof files
  Int_t runMin=-1;
  Int_t runMax=-1;

  for (int irun=1300; irun<1400; irun++){
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

      if (tof > piLow && tof < piHi && tofCoin->bar!=10 && tofCoin->inSpill) {
	hpion->Fill((dstofHitT - tofCoin->lastDelayedBeamSignal)/1e9);
      }
      else if (tof > proLow && tof < proHi && tofCoin->bar!=10 && tofCoin->inSpill) {
	hproton->Fill((dstofHitT - tofCoin->lastDelayedBeamSignal)/1e9);
      }
    }
  }
  gStyle->SetOptStat(0);

  hpion->Scale(1./(double)nSpillsTrue);
  hproton->Scale(1./(double)nSpillsTrue);

  TCanvas *c1 = new TCanvas("c1");
  hratio->Divide(hproton, hpion, 1., 1., "B");
  hratio->Draw("e");
  hratio->GetXaxis()->SetLabelSize(0.05);
  hratio->GetYaxis()->SetLabelSize(0.05);
  hratio->GetXaxis()->SetTitleSize(0.05);
  hratio->GetYaxis()->SetTitleSize(0.05);
  c1->Print(Form("%s/ratio.png", saveDir));
  c1->Print(Form("%s/ratio.pdf", saveDir));
  TCanvas *c2 = new TCanvas("c2");
  TLegend *leg = new TLegend(0.6, 0.6, 0.85, 0.85);
  leg->AddEntry(hpion, "MIPs", "f");
  leg->AddEntry(hproton, "Protons", "f");
  hpion->SetLineColor(kRed);
  hpion->SetFillColor(kRed);
  hproton->SetLineColor(kBlue);
  hproton->SetFillColor(kBlue);
  hpion->SetFillStyle(3004);
  hproton->SetFillStyle(3005);
  hsProPi->Add(hpion);
  hsProPi->Add(hproton);
  hsProPi->Draw("hist e nostack");
  leg->Draw();
  hsProPi->GetXaxis()->SetLabelSize(0.05);
  hsProPi->GetYaxis()->SetLabelSize(0.05);
  hsProPi->GetXaxis()->SetTitleSize(0.05);
  hsProPi->GetYaxis()->SetTitleSize(0.05);
  c2->Print(Form("%s/proPi.png", saveDir));
  c2->Print(Form("%s/proPi.pdf", saveDir));

} // proPiFracUnix
