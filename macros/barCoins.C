// barCoins.C
// Big matrix of coincidences between various bars

void barCoins(const char* saveDir, 
	      const char* dstofDir="/nfs/scratch0/dbrailsf/data_backup/dtof_backup/",
	      const char* ustofDir="/nfs/scratch0/dbrailsf/data_backup/utof_backup_firsthitpinnedtounixtime/Data_root_v3_wo_walk_corr/")
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
  std::vector<const char*> str4BlockVec = {/*str4Block0*/ str4Block1, str4Block2, str4Block3};

  const double startLongTime = 1535997500;
  const double endLongTime   = 1536083900;

  const double coinCut = 1.; // Cut for determining if a coincidence has occurred
  const double deadtimeCut = 185.;

  TFile *fout = new TFile(saveDir, "recreate");

  for (int sample=0; sample<6; sample++) {
    cout<<"Sample "<<sample+1<<endl;
    int totalTime = 0;
    string name;
    vector<double> startTimes;
    vector<double> endTimes;
    bool isMultiple = false;
    startTimes.clear();
    endTimes.clear();

    if (sample==0) {
      name="0";
      startTimes.push_back(start0Block);
      endTimes.push_back(end0Block);
    }
    else if (sample==1) {
      name="1";
      startTimes.push_back(start1Block);
      endTimes.push_back(end1Block);
    }
    else if (sample==2) {
      name="2";
      startTimes.push_back(start2Block);
      endTimes.push_back(end2Block);
    }
    else if (sample==3) {
      name="3";
      startTimes.push_back(start3Block);
      endTimes.push_back(end3Block);
    }
    else if (sample==4) {
      name="4";
      isMultiple = true;
      for (int b4=0; b4<str4BlockVec.size(); b4++) {
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
	int startTime = stoi(unixstart);
	treeTmp->GetEntry(treeTmp->GetEntries() - 1);
	int endTime = startTime + (tS1Tmp/1e9);
	futofTmp->Close();
	delete futofTmp;
	startTimes.push_back(startTime);
	endTimes.push_back(endTime);
      }
    }
    else if (sample==5) {
      name="long";
      startTimes.push_back(startLongTime);
      endTimes.push_back(endLongTime);
    }

    TH1D *hCosmicRate = new TH1D(Form("hCosmicRate_%s", name.c_str()), "Cosmic rate in S4; Bar; Rate / s^{-1}", 10, 0.5, 10.5);
    hCosmicRate->GetXaxis()->SetTitleSize(.05);
    hCosmicRate->GetXaxis()->SetLabelSize(.05);
    hCosmicRate->GetYaxis()->SetTitleSize(.05);
    hCosmicRate->GetYaxis()->SetLabelSize(.05);
    hCosmicRate->GetZaxis()->SetTitleSize(.05);
    hCosmicRate->GetZaxis()->SetLabelSize(.05);
    hCosmicRate->Sumw2();
    TH2D *h2CosmicRate = new TH2D(Form("h2CosmicRate_%s", name.c_str()), "Cosmic rate in S4; Bar position / cm; Bar; Rate / s^{-1}", 20, 0., 140., 10, 0.5, 10.5);
    h2CosmicRate->GetXaxis()->SetTitleSize(.05);
    h2CosmicRate->GetXaxis()->SetLabelSize(.05);
    h2CosmicRate->GetYaxis()->SetTitleSize(.05);
    h2CosmicRate->GetYaxis()->SetLabelSize(.05);
    h2CosmicRate->GetZaxis()->SetTitleSize(.05);
    h2CosmicRate->GetZaxis()->SetLabelSize(.05);
    h2CosmicRate->Sumw2();
    TH2D *h2matrix = new TH2D(Form("matrix_%s", name.c_str()), "Bar coincidences; Bar 1; Bar 2; Rate / s^{-1}", 10, 0.5, 10.5, 10, 0.5, 10.5);
    h2matrix->GetXaxis()->SetTitleSize(.05);
    h2matrix->GetXaxis()->SetLabelSize(.05);
    h2matrix->GetYaxis()->SetTitleSize(.05);
    h2matrix->GetYaxis()->SetLabelSize(.05);
    h2matrix->GetZaxis()->SetTitleSize(.05);
    h2matrix->GetZaxis()->SetLabelSize(.05);

    vector<TH1D*> cosmicRateVec;
    vector<double> cosmicsLeftVec;
    vector<double> cosmicsRightVec;
    vector<double> asymmetryVec;
    cosmicsLeftVec.resize(10, 0.);
    cosmicsRightVec.resize(10, 0.);
    asymmetryVec.resize(10, 0.);
    for (int b=0; b<10; b++) {
      TH1D *hCosmicRateBar = new TH1D(Form("hCosmicRateBar_%s_%d", name.c_str(), b+1), Form("Bar %d; Horizontal position / cm; Rate / s^{-1}", b+1), 20., 0, 140.);
      hCosmicRateBar->GetXaxis()->SetTitleSize(.05);
      hCosmicRateBar->GetXaxis()->SetLabelSize(.05);
      hCosmicRateBar->GetYaxis()->SetTitleSize(.05);
      hCosmicRateBar->GetYaxis()->SetLabelSize(.05);
      hCosmicRateBar->SetLineColor(52 + b*3);
      hCosmicRateBar->SetLineWidth(2);
      hCosmicRateBar->Sumw2();
      cosmicRateVec.push_back(hCosmicRateBar);
    }

    // Loop over samples
    for (int n=0; n<startTimes.size(); n++) {
      cout<<"Subsample "<<n+1<<" of "<<startTimes.size()<<endl;
      int startTime = startTimes.at(n);
      int endTime   = endTimes.at(n);
      int nSpills = 0;
      double lastSpill = 0.;
      vector< vector<double> > spillsVec;
      vector<double> spills1;
      vector<double> spills2;
      spillsVec.push_back(spills1);
      spillsVec.push_back(spills2);
      vector< vector<short> > runVec;
      vector<short> run1;
      vector<short> run2;
      runVec.push_back(run1);
      runVec.push_back(run2);
      vector<double> deadtimeVec;
      deadtimeVec.resize(10, 0.);

      // Find dtof runs
      int runMin = -1;
      int runMax = -1;
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

      // First match hits between bars in the same TDC
      for (int itdc = 0; itdc < 2; itdc++) {
	TChain *tofCoinChain = new TChain("tofCoinTree");
	for (int irun=runMin; irun<runMax+1; irun++) {
	  // Load input files
	  tofCoinChain->Add(Form("%srun%d/DsTOFcoincidenceRun%d_tdc%d.root", dstofDir, irun, irun, itdc+1));
	} // Loop over runs
	RawDsTofCoincidence *tofCoin = NULL;
	tofCoinChain->SetBranchAddress("tofCoin", &tofCoin);
	for (int h=0; h<tofCoinChain->GetEntries(); h++) {
	  tofCoinChain->GetEntry(h);
	  if (h % 500000 == 0) cout<<"Entry "<<h<<" of "<<tofCoinChain->GetEntries()<<endl;

	  if (tofCoin->unixTime[0]<startTime) continue;
	  if (tofCoin->unixTime[0]>endTime) break;

	  if (tofCoin->lastDelayedBeamSignal != lastSpill) {
	    lastSpill = tofCoin->lastDelayedBeamSignal;
	    spillsVec[itdc].push_back(tofCoin->lastDelayedBeamSignal);
	    runVec[itdc].push_back(tofCoin->run);
	  }

	  if (tofCoin->inSpill) continue;

	  int bar1 = tofCoin->bar;
	  double dstofHitT1 = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (10.-TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1])/2.);
	  if (!tofCoin->inSpill && (dstofHitT1 - deadtimeVec.at(bar1-1)) > deadtimeCut) {
	    double xpos = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 70.));
	    cosmicRateVec.at(bar1-1)->Fill(xpos);
	    (xpos > 70.) ? cosmicsRightVec.at(bar1-1)++ : cosmicsLeftVec.at(bar1-1)++;
	    h2CosmicRate->Fill(xpos, bar1);
	    hCosmicRate->Fill(bar1);
	    deadtimeVec.at(bar1-1) = dstofHitT1;
	    // Skip forward in this file to see if there any coincidences
	    for (int h2=h; h2<tofCoinChain->GetEntries(); h2++) {
	      tofCoinChain->GetEntry(h2);
	      double dstofHitT2 = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (10.-TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1])/2.);
	      if (dstofHitT2-dstofHitT1 > 20.) break;
	      int bar2 = tofCoin->bar;
	      if (!tofCoin->inSpill && (dstofHitT2 - dstofHitT1) < coinCut && 
		  dstofHitT2 - dstofHitT1 > 0. && bar1 != bar2 && 
		  dstofHitT2 - deadtimeVec.at(bar2-1) > deadtimeCut) {
		deadtimeVec.at(bar2-1) = dstofHitT2;
		h2matrix->Fill(bar1, bar2);
	      } // Is not in a spill
	    } // Loop over TChain for a second time
	  } // Is not in a spill
	} // Loop over entries

	deadtimeVec.clear();
	deadtimeVec.resize(10, 0.);
	lastSpill = 0.;

	delete tofCoin;
	delete tofCoinChain;
      } // Loop over TDCs

      bool spillsMatch = false;
      if (spillsVec[0].size() == spillsVec[1].size()) {
	cout<<spillsVec[0].size()<<" spills"<<endl;
	spillsMatch = true;
      }
      else {
	cout<<"Different number of spills between TDCs"<<endl;
      }

      // Now need to match between TDCs. Spill signals have a constant offset for each run 
      // Calculate this offset then apply it
      vector<double> offsetVec;
      for (int irun = runMin; irun < runMax+1; irun++) {
	cout<<"Run "<<irun<<endl;
	double avg = 0.;
	int nGoodSpills = 0;
	cout.precision(11);
	for (int i=0; i<spillsVec[0].size(); i++) {
	  if ((spillsVec[0].at(i) - spillsVec[1].at(i)) != 0 && runVec[0].at(i) == irun) {
	    avg += spillsVec[0].at(i) - spillsVec[1].at(i);
	    nGoodSpills++;
	  }
	}
	avg /= nGoodSpills;
	cout<<"Average offset = "<<avg<<endl;
	offsetVec.push_back(avg);
	for (int i=0; i<spillsVec[0].size(); i++) {
	  if (abs(spillsVec[0].at(i) - spillsVec[1].at(i) - avg) > 1. &&
	      spillsVec[0].at(i) - spillsVec[1].at(i) != 0 &&
	      runVec[0].at(i) == irun) { 
	    cout<<"Spill is far from average!!!"<<endl; 
	    cout<<"Offset: "<<(spillsVec[0].at(i) - spillsVec[1].at(i))<<" Diff: "<<(spillsVec[0].at(i) - spillsVec[1].at(i) - avg)<<endl;
	  }
	}
	// Now use calculated average offset to do the matching between TDCs
	// We are only matching in one direction so need to do this twice
	for (int itdc=0; itdc < 2; itdc++) {
	  cout<<"Cross TDC hits, TDC "<<itdc+1<<endl;
	  int lasth2 = 0;
	  int tdc1 = (itdc == 0) ? 1 : 2;
	  int tdc2 = (itdc == 0) ? 2 : 1;
	  TFile *tofFile1 = new TFile(Form("%srun%d/DsTOFcoincidenceRun%d_tdc%d.root", dstofDir, irun, irun, tdc1), "read");
	  TFile *tofFile2 = new TFile(Form("%srun%d/DsTOFcoincidenceRun%d_tdc%d.root", dstofDir, irun, irun, tdc2), "read");
	  TTree *tofTree1 = (TTree*)tofFile1->Get("tofCoinTree");
	  TTree *tofTree2 = (TTree*)tofFile2->Get("tofCoinTree");
	  RawDsTofCoincidence *tof1 = NULL;
	  RawDsTofCoincidence *tof2 = NULL;
	  tofTree1->SetBranchAddress("tofCoin", &tof1);
	  tofTree2->SetBranchAddress("tofCoin", &tof2);
	  for (int h=0; h<tofTree1->GetEntries(); h++) {
	    if (h % 100000 == 0) cout<<"Entry "<<h<<" of "<<tofTree1->GetEntries()<<endl;
	    tofTree1->GetEntry(h);
	    if (tof1->unixTime[0]<startTime) continue;
	    if (tof1->unixTime[0]>endTime) break;
	    if (tof1->inSpill) continue;
	    int bar1 = tof1->bar;
	    double dstofHitT1 = min(tof1->fakeTimeNs[0], tof1->fakeTimeNs[1]) - (10.-TMath::Abs(tof1->fakeTimeNs[0]-tof1->fakeTimeNs[1])/2.);
	    if (!tof1->inSpill && (dstofHitT1 - deadtimeVec.at(bar1-1)) > deadtimeCut) {
	      deadtimeVec.at(bar1-1) = dstofHitT1;
	      for (int h2=lasth2; h2<tofTree2->GetEntries(); h2++) {
		tofTree2->GetEntry(h2);
		double dstofHitT2 = min(tof2->fakeTimeNs[0], tof2->fakeTimeNs[1]) - (10.-TMath::Abs(tof2->fakeTimeNs[0]-tof2->fakeTimeNs[1])/2.);
		dstofHitT2 = (itdc==0) ? dstofHitT2 + avg : dstofHitT2 - avg;

		if (dstofHitT2 < dstofHitT1) {
		  lasth2 = h2;
		  continue;
		}

		if (dstofHitT2 - dstofHitT1 > 20.) break;
		int bar2 = tof2->bar;
		if (!tof2->inSpill && (dstofHitT2 - deadtimeVec.at(bar2-1)) > deadtimeCut &&
		    (dstofHitT2 - dstofHitT1) < coinCut && (dstofHitT2 - dstofHitT1) > 0.) {
		  deadtimeVec.at(bar2-1) = dstofHitT2;
		  h2matrix->Fill(bar1, bar2);
		  lasth2 = h2;
		}
	      } // Loop over second TDC entries
	    }
	  } // Loop over first TDC entries

	  delete tof1;
	  delete tof2;
	  tofFile1->Close();
	  tofFile2->Close();
	  delete tofFile1;
	  delete tofFile2;
	  deadtimeVec.clear();
	  deadtimeVec.resize(10, 0.);
	} // Loop over TDCs
      } // Loop over runs 
      totalTime += endTime - startTime - spillsVec[0].size();

    } // Loop over the sub samples within a sample
    fout->cd();
    THStack *hsBarRates = new THStack(Form("hsBarRates_%s", name.c_str()), Form("%s block cosmic rates; Horizontal position / cm; Rate / s^{-1}", name.c_str()));
    TGraph *grAsymm = new TGraph();
    grAsymm->SetTitle(Form("Cosmic ray asymmetry, %s block; Bar; Asymmetry", name.c_str()));
    int cosmicsLeft = 0;
    int cosmicsRight = 0;
    for (int i=0; i<cosmicRateVec.size(); i++) {
      cosmicRateVec.at(i)->Scale(1. / totalTime);
      cosmicRateVec.at(i)->Write();
      grAsymm->SetPoint(i, i+1, (cosmicsRightVec.at(i)-cosmicsLeftVec.at(i))/(cosmicsRightVec.at(i)+cosmicsLeftVec.at(i)));
      cosmicsLeft += cosmicsLeftVec.at(i);
      cosmicsRight += cosmicsRightVec.at(i);
      hsBarRates->Add(cosmicRateVec.at(i));
    }
    cout<<"Asymmetry: "<<((double)cosmicsRight-(double)cosmicsLeft)/((double)cosmicsRight+(double)cosmicsLeft)<<endl;

    hCosmicRate->Scale(1. / totalTime);
    h2CosmicRate->Scale(1. / totalTime);
    h2matrix->Scale(1. / totalTime);
    hCosmicRate->Write();
    h2CosmicRate->Write();
    h2matrix->Write();
    hsBarRates->Write();
    grAsymm->Write(Form("grAsymm_%s", name.c_str()));
  }
  TLine *lh = new TLine(0.5, 5.5, 10.5, 5.5);
  TLine *lv = new TLine(5.5, 0.5, 5.5, 10.5);
  lh->SetLineWidth(2);
  lh->SetLineStyle(2);
  lv->SetLineWidth(2);
  lv->SetLineStyle(2);
  lh->Write("lh");
  lv->Write("lv");

  fout->Close();
  delete fout;
} // barCoins
