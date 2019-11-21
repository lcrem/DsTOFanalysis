// mcToData.C
// Take the Monte Carlo and imply the inverse of the corrections + see what we get out

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
const char* str0Block = "Data_2018_9_1_b4_800MeV_1block_bend4cm.root";
const char* str1Block = "Data_2018_9_1_b2_800MeV_2block_bend4cm.root";
const char* str2Block = "Data_2018_9_1_b3_800MeV_3block_bend4cm.root";
// New datasets to combine for the 4 block data
const char* str4Block0 = "Data_2018_8_28_b5.root";
const char* str4Block1 = "Data_2018_8_30_b1.root";
const char* str4Block2 = "Data_2018_8_29_b4.root";
const char* str4Block3 = "Data_2018_8_29_b1.root";
std::vector<const char*> str4BlockVec = {str4Block0, str4Block1, str4Block2};
const double s4BarTime = 18.2;
const double s4OffAxisStartX   = 0.121;
const double s4OffAxisEndX     = 1.4224;
const double barOverallEff = 0.75;
const vector<double> dataS3  = {1983., 1656., 1325., 899., 136.3};
const vector<double> mcS3    = {56590, 67273, 60652, 35781, 104318};
// For smearing histograms
const std::vector<double> timeRes  = {1, 1.2, 1.4, 1.6, 1.8, 2.}; // in ns
const std::vector<double> effWidth = {0.2, 0.22, 0.24, 0.26}; // one sigma values for efficiency curves

void setHistAttr(TH1D *h) 
{
  h->SetLineWidth(2);
  h->GetXaxis()->SetTitleSize(.05);
  h->GetYaxis()->SetTitleSize(.05);
  h->GetXaxis()->SetLabelSize(.05);
  h->GetYaxis()->SetLabelSize(.05);
  h->Sumw2();
  h->SetOption("hist");
}

void setHistAttr(TH2D *h2) 
{
  h2->SetOption("colz");
  h2->GetXaxis()->SetTitleSize(.05);
  h2->GetYaxis()->SetTitleSize(.05);
  h2->GetZaxis()->SetTitleSize(.05);
  h2->GetXaxis()->SetLabelSize(.05);
  h2->GetYaxis()->SetLabelSize(.05);
  h2->GetZaxis()->SetLabelSize(.05);
  h2->Sumw2();
}

void mcToData(const char* saveDir, 
	      const char* s4Plots="/scratch0/sjones/plots/angularDistS4_newSample/withTree/includesUnweightedHists.root", 
	      const char* mcDir="/scratch0/tnonnenm/FixedMC/Prototype-HPTPC-MC/wdir/moderator_data/",
	      const char* ustofDir="/nfs/scratch0/dbrailsf/data_backup/utof_backup_firsthitpinnedtounixtime/Data_root_v3_wo_walk_corr/",
	      const char* dstofDir="/nfs/scratch0/dbrailsf/data_backup/dtof_backup/",
	      const char* smearPlots="/scratch0/sjones/plots/barEffSmearing/gauss4.root")
{
  TFile *fout = new TFile(saveDir, "recreate");
  TFile *fS4 = new TFile(s4Plots, "read"); // S4 plots
  TFile *fSmear = new TFile(smearPlots, "read"); // Smearing histograms

  THStack *hsRatio = new THStack("hsRatio", "Data / MC for multiple samples; x / m; Data / MC");

  for (int sample=0; sample<5; sample++) {
    cout<<"Sample "<<sample+1<<endl;
    string name;
    vector<double> startTimes;
    vector<double> endTimes;
    bool isMultiple = false;
    startTimes.clear();
    endTimes.clear();

    double totalTime = 0.;

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
      } // Loop over 4 block sub-samples
    }

    TH2D *h2CosmicRate = new TH2D(Form("h2CosmicRate%s", name.c_str()), "Cosmic rate in S4; x / m; y / m; Rate / s^{-1}", 20, -1., 0.4, 10, -0.35, 0.45);
    setHistAttr(h2CosmicRate);
    TH2D *h2CosmicEff = new TH2D(Form("h2CosmicEff%s", name.c_str()), "Cosmic efficiency in S4; x / m; y / m; Efficiency", 20, -1., 0.4, 10, -0.35, 0.45);
    setHistAttr(h2CosmicEff);
    TH2D *h2MCTruth = new TH2D(Form("h2MCTruth%s", name.c_str()), "Proton Monte Carlo truth; x / m; y / m; S4 / S3", 20., -1., 0.4, 10., -.35, 0.45);
    setHistAttr(h2MCTruth);
    TH2D *h2MCWeighted = new TH2D(Form("h2MCWeighted%s", name.c_str()), "Proton Monte Carlo weighted; x / m; y / m; S4 / S3", 20., -1., 0.4, 10., -.35, 0.45);
    setHistAttr(h2MCWeighted);
    TH2D *h2MCWeightedSmooth = new TH2D(Form("h2MCWeightedSmooth%s", name.c_str()), "Proton Monte Carlo weighted; x / m; y / m; S4 / S3", 20., -1., 0.4, 10., -.35, 0.45);
    setHistAttr(h2MCWeightedSmooth);
    TH1D *hMCTruth = new TH1D(Form("hMCTruth%s", name.c_str()), "Proton Monte Carlo truth; x / m; S4 / S3", 20., -1., 0.4);
    setHistAttr(hMCTruth);
    TH1D *hMCWeighted = new TH1D(Form("hMCWeighted%s", name.c_str()), "Proton Monte Carlo weighted; x / m; S4 / S3", 20., -1., 0.4);
    setHistAttr(hMCWeighted);
    TH1D *hMCWeightedSmooth = new TH1D(Form("hMCWeightedSmooth%s", name.c_str()), "Proton Monte Carlo weighted; x / m; S4 / S3", 20., -1., 0.4);
    setHistAttr(hMCWeightedSmooth);

    for (int n=0; n<startTimes.size(); n++) {
      cout<<"Subsample "<<n+1<<" of "<<startTimes.size()<<endl;
      int startTime = startTimes.at(n);
      int endTime   = endTimes.at(n);
      int nSpills = 0;
      double lastSpill = 0.;

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

	  if (tofCoin->lastDelayedBeamSignal != lastSpill && itdc == 0) {
	    lastSpill = tofCoin->lastDelayedBeamSignal;
	    nSpills++;
	  }

	  if (tofCoin->inSpill || abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1]) > s4BarTime) continue;

	  double dstofHitT1 = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (s4BarTime/2.-TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1])/2.);
	  // Is a cosmic
	  if (!tofCoin->inSpill && abs(tofCoin->fakeTimeNs[0] - tofCoin->fakeTimeNs[1]) < s4BarTime) {
	    // Apply offsets to get these in same frame as MC
	    double mcX = -((((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(0.7/s4BarTime) + .65) / 1.3) * (s4OffAxisEndX - s4OffAxisStartX) + s4OffAxisStartX) + 0.491;
	    double mcY = (tofCoin->bar * 0.075) - 0.385 + 0.0114;
	    h2CosmicRate->Fill(mcX, mcY);
	    h2CosmicEff->Fill(mcX, mcY);
	  } // Is a cosmic
	} // Loop over entries
      } // Loop over TDCs
      totalTime += endTime - startTime - nSpills;
    } // Loop over start times

    // Normalise to highest bar segment
    h2CosmicEff->Scale(1. / h2CosmicEff->GetBinContent(h2CosmicEff->GetMaximumBin()));
    h2CosmicRate->Scale(1. / totalTime);

    TH2D *h2CosmicEffSmooth = (TH2D*)h2CosmicEff->Clone(Form("h2CosmicEffSmooth%d", sample));
    h2CosmicEffSmooth->Smooth(1, "k3a");
    // Do MC
    TFile *mcFile   = new TFile(Form("%s/%dBlock.root", mcDir, sample), "read");
    TTree *stepTree = (TTree*)mcFile->Get("stepTree");
    UInt_t NSteps;
    float PostStepMom[40000];
    float PostStepPosX[40000];
    float PostStepPosY[40000];
    int PreStepVolumeID[40000];
    int PostStepVolumeID[40000];
    stepTree->SetBranchAddress("NSteps", &NSteps);
    stepTree->SetBranchAddress("PostStepMom", PostStepMom);
    stepTree->SetBranchAddress("PostStepPosX", PostStepPosX);
    stepTree->SetBranchAddress("PostStepPosY", PostStepPosY);
    stepTree->SetBranchAddress("PreStepVolumeID", PreStepVolumeID);
    stepTree->SetBranchAddress("PostStepVolumeID", PostStepVolumeID);
    // Loop over entries
    for (int i=0; i<stepTree->GetEntries(); i++) {
      stepTree->GetEntry(i);
      // Loop over steps
      for (int n=0; n<NSteps; n++) {
	// Condition for a proton being in S4
	if (PreStepVolumeID[n]!=7 && PostStepVolumeID[n]==7 && PostStepMom[n] > 140.) {
	  h2MCTruth->Fill(PostStepPosX[n]/1000., PostStepPosY[n]/1000.);	  
	  hMCTruth->Fill(PostStepPosX[n]/1000.);
	  h2MCWeighted->Fill(PostStepPosX[n]/1000., PostStepPosY[n]/1000., barOverallEff*h2CosmicEff->GetBinContent(h2CosmicEff->GetXaxis()->FindBin(PostStepPosX[n]/1000.), h2CosmicEff->GetYaxis()->FindBin(PostStepPosY[n]/1000.)));
	  hMCWeighted->Fill(PostStepPosX[n]/1000., barOverallEff*h2CosmicEff->GetBinContent(h2CosmicEff->GetXaxis()->FindBin(PostStepPosX[n]/1000.), h2CosmicEff->GetYaxis()->FindBin(PostStepPosY[n]/1000.)));
	}
      }
    } // Loop over tree entries  

    // Now scale by the appropriate S3 number (with bin width) to give comparison to data
    h2MCTruth->Scale(h2MCTruth->GetNbinsX()*h2MCTruth->GetNbinsY()/mcS3.at(sample));
    h2MCWeighted->Scale(h2MCWeighted->GetNbinsX()*h2MCWeighted->GetNbinsY()/mcS3.at(sample));
    hMCTruth->Scale(hMCTruth->GetNbinsX()/mcS3.at(sample));
    hMCWeighted->Scale(hMCWeighted->GetNbinsX()/mcS3.at(sample));

    fout->cd();
    h2MCTruth->Write();
    h2MCWeighted->Write();
    hMCTruth->Write();
    hMCWeighted->Write();
    h2CosmicEff->Write();
    h2CosmicEffSmooth->Write();
    h2CosmicRate->Write();

    // Get S4 proton tree and use this to compare
    TTree *tree = (TTree*)fS4->Get(Form("protonTree%d", sample));
    TH1D *hData = new TH1D(Form("hData%d", sample), Form("Data, %d blocks; x / m; S4 / S3", sample), 20, -1., 0.4);
    double mcX, mom;
    int spill;
    tree->SetBranchAddress("mcX", &mcX);
    tree->SetBranchAddress("spill", &spill);
    tree->SetBranchAddress("mom", &mom);
    setHistAttr(hData);
    for (int t=0; t<tree->GetEntries(); t++) {
      tree->GetEntry(t);
      if (mom > 140.) hData->Fill(mcX);
    } // Loop over entries
    // Divide by number of spills and S3 values
    tree->GetEntry(tree->GetEntries()-1);
    hData->Scale(hData->GetNbinsX() / (spill * dataS3.at(sample)));
    fout->cd();
    hData->Write();
    // Makes the stacks
    hMCWeighted->SetLineColor(kBlack);
    hData->SetLineColor(kRed);
    THStack *hsDataMC = new THStack(Form("hsDataMC%d", sample), Form("Data - MC comp, %d; x / m; S4 / S3 ", sample));
    hsDataMC->Add(hData);
    hsDataMC->Add(hMCWeighted);
    hsDataMC->Write();

    THStack *hsRatioSample = new THStack(Form("hsRatio%d", sample), Form("Data / MC %d blocks; x / m; Data / MC", sample));
    for (int t=0; t<timeRes.size(); t++) {
      THStack *hsRatioSampleRes = new THStack(Form("hsRatio%dRes%d", sample, (int)(timeRes.at(t)*10.)), Form("Data / MC: %d blocks, %.2g ns; x / m; Data / MC", sample, timeRes.at(t)));
      for (int eff=0; eff < effWidth.size(); eff++) {
	TH1D *hSmear = (TH1D*)fSmear->Get(Form("hRatioEff%dRes%d", (int)(effWidth.at(eff)*100.), (int)(timeRes.at(t)*10.)));
	THStack *hsDataMCSmear = new THStack(Form("hsDataMCSmear%dRes%dEff%d", sample, (int)(timeRes.at(t)*10.), (int)(effWidth.at(eff)*100.)), Form("Data - MC comp with %d cm eff and %.2g ns smearing, %d blocks; x / m; S4 / S3 ", (int)(effWidth.at(eff)*100), timeRes.at(t), sample));
	TH1D *hMCSmear = (TH1D*)hMCWeighted->Clone(Form("hMCSmear%dRes%dEff%d", sample, (int)(timeRes.at(t)*10.), (int)(effWidth.at(eff))));

	// Weight each bin in the weighted MC histogram by 1 / the relevant bin in the smearing hist
	for (int b=0; b<=hSmear->GetNbinsX(); b++) {
	  hMCSmear->SetBinContent(b, hMCSmear->GetBinContent(b) / hSmear->GetBinContent(b));
	}
	fout->cd();
	hsDataMCSmear->Add(hData);
	hsDataMCSmear->Add(hMCSmear);
	hsDataMCSmear->Write();

	TH1D *hRatio = new TH1D(Form("hRatio%dRes%dEff%d", sample, (int)(timeRes.at(t)*10.), (int)(effWidth.at(eff)*100.)), Form("Res = %.2g ns, eff width %.2g m; x position / m; Data / MC", timeRes.at(t), effWidth.at(eff)), 20., -1., 0.4); 
	hRatio->Divide(hData, hMCSmear, 1., 1., "B");
	hRatio->SetLineColor(52 + (eff+t) * 4);
	hRatio->SetLineWidth(2);
	hsRatioSample->Add(hRatio);
	hRatio->SetLineColor(52 + eff * 6);
	hsRatioSampleRes->Add(hRatio);
	hRatio->Write();
      }
      hsRatioSampleRes->Write();
    }
    hsRatioSample->Write();

    TH1D *hRatio = new TH1D(Form("hRatio%s", name.c_str()), Form("Data / MC ratio, %s blocks; x position / m; Data / MC", name.c_str()), 20., -1., 0.4); 
    setHistAttr(hRatio);
    hRatio->Divide(hData, hMCWeighted, 1., 1., "B");
    hRatio->SetLineColor(52 + sample * 4);
    hRatio->Write();
    hsRatio->Add(hRatio);
    hsRatioSample->Write();
  } // Loop over samples
  fout->cd();
  hsRatio->Write();

  fSmear->Close();
  fS4->Close();
  fout->Close();
  delete fout;
  delete fS4;
  delete fSmear;
}

