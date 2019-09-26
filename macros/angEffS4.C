// angEffS4.C
// angular efficiency using cosmics for S4

void angEffS4(const char* saveDir, 
	      const char* dstofDir="/nfs/scratch0/dbrailsf/data_backup/dtof_backup/",
	      const char* ustofDir="/nfs/scratch0/dbrailsf/data_backup/utof_backup_firsthitpinnedtounixtime/Data_root_v3_wo_walk_corr/",
	      const char* spillDBDir="/scratch0/sjones/spillDB/") {
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
  // 0.8GeV/c, 4 block
  // 4 moderator blocks with -4cm bend
  const double start4Block = 1535836129;
  const double end4Block   = 1535879634;

  // New datasets to combine for the 4 block data
  const char* str4Block0 = "Data_2018_8_28_b5.root";
  const char* str4Block1 = "Data_2018_8_30_b1.root";
  const char* str4Block2 = "Data_2018_8_29_b4.root";
  const char* str4Block3 = "Data_2018_8_29_b1.root";
  std::vector<const char*> str4BlockVec = {str4Block0, str4Block1, str4Block2, str4Block3};

  std::vector<double> startTimeVec = {start0Block, start1Block, start2Block, 
				      start3Block, start4Block}; 
  std::vector<double> endTimeVec   = {end0Block, end1Block, end2Block, 
				      end3Block, end4Block}; 

  TFile *fo = new TFile(Form("%s/angEffS4.root", saveDir), "recreate");

  // Get the start and end times for the remaining 4 block data sets
  for (int b4=0; b4<str4BlockVec.size(); b4++) {
    // Find the correct dstof files
    double startTime = 0;
    double endTime   = 0;
    // Get start and end times from the relevant utof files
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
    startTime = stoi(unixstart);
    treeTmp->GetEntry(treeTmp->GetEntries() - 1);
    endTime = startTime + (tS1Tmp/1e9);
    futofTmp->Close();
    delete futofTmp;

    startTimeVec.push_back(startTime);
    endTimeVec.push_back(endTime);
  } // Loop over the remaining four block data

  for (int file=0; file<startTimeVec.size(); file++) {
    Int_t runMin=-1;
    Int_t runMax=-1;
    double startTime = startTimeVec.at(file);
    double endTime   = endTimeVec.at(file);

    TH2D *h2Cosmics = new TH2D(Form("h2Cosmics_%d", file), "Cosmics in 2D S4; x / cm; Bar; Rate / Hz", 20, 0., 140., 10, 0.5, 10.5);
    TH1D *hCosmics  = new TH1D(Form("hCosmics_%d", file), "Cosmics in 1D S4; x / cm; Rate / Hz", 20, 0., 140.);
    std::vector<TH1D*> eff1dVec;
    for (int b=0; b<10; b++) {
      TH1D *hEff1d = new TH1D(Form("hEff1d_%d_bar%d", file, b+1), Form("Bar %d efficiency; x / cm; Eff", b+1), 20, 0., 140.);
      eff1dVec.push_back(hEff1d);
    }

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
    int nSpills = 0;
    double lastDelayed = 0.;
    // Now loop over the coincidence files again and calculate the angular distributions
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

	if (itdc==0 && tofCoin->lastDelayedBeamSignal != lastDelayed) {
	  nSpills++;
	  lastDelayed = tofCoin->lastDelayedBeamSignal;
	}
	// Only want hits out of spill - i.e. cosmics
	if (!tofCoin->inSpill) {
	  double positionX = (tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(7./2.) + 70.;
	  h2Cosmics->Fill(positionX, tofCoin->bar);
	  eff1dVec.at(tofCoin->bar-1)->Fill(positionX);
	  hCosmics->Fill(positionX);
	} // if (!tofCoin->InSpill)
      } // for (int h=0; h<tofCoinChain->GetEntries(); h++)
    } // Loop over both TDCs

    fo->cd();
    // Scale by number of seconds (recall spills are 1s long)
    h2Cosmics->Scale(1. / (endTime - startTime - nSpills) );
    hCosmics->Scale(1. / (endTime - startTime - nSpills) );
    h2Cosmics->Write();
    hCosmics->Write();
    for (int bar=0; bar<eff1dVec.size(); bar++) {
      eff1dVec.at(bar)->Scale(1./eff1dVec.at(bar)->GetBinContent(eff1dVec.at(bar)->GetMaximumBin()));
      eff1dVec.at(bar)->Write();
    }
  } // Loop over all the different datasets

  fo->Close();
  delete fo;
} // angEffS4
