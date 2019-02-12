// angularDistS4.C
// Angular distribution of protons and pions for different moderator blocks
// Use efficiency calculation and background subtraction
void angularDistS4(const char* saveDir, 
		   const char* dstofDir="/scratch0/dbrailsf/temp/mylinktodtof/") 
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
  // Most runs were in this configuration so don't need to use necessarily
  //  const double start4Block = 1536537600; 
  //  const double end4Block   = 1536669600;
  // Timing cuts
  const double piLow  = 40.;
  const double piHi   = 55.;
  const double proLow = 66.;
  const double proHi  = 94.;
  // S1 -> S4 positions
  const double baselineS1S4End   = 13.9426;
  const double baselineS1S4Start = 14.0069;
  const double s4OffAxisStartX   = 0.121;
  const double s4OffAxisEndX     = 1.4224;
  // Shift in ns required to to pion peak at speed of light
  const double dstofShift = 40.;
  // Ustof-dstof cable delay
  const double ustofDelay = 184.7;

  TFile *fout = new TFile(Form("%s/angularDistS4Plots.root", saveDir), "recreate");

  for (int nBlocks = 0; nBlocks < 4; nBlocks++) {
    // Number of signal particles using just cut and count
    double nP  = 0.;
    double nPi = 0.;
    // Define signal and background functions to be fitted
    // Signals are gaussians
    TF1 *sPro = new TF1("sPro", "gaus", 66, 100);
    TF1 *sPi  = new TF1("sPi", "gaus", 35, 55);
    // Exponential background
    TF1 *fBkgExp = new TF1("fBkgExp","expo", 30, 160);
    sPro->SetLineColor(kGreen+2);
    sPi->SetLineColor(kRed);
    TF1 *fSplusBExp = new TF1("signal+bkg exp", "gaus(0)+gaus(3)+expo(6)", 30, 160);
    fSplusBExp->SetParNames("const 1", "mean 1", "sigma 1",
			    "const 2", "mean 2", "sigma 2",
			    "bkgconst", "bkgdecay");
    fSplusBExp->SetLineColor(kBlack);
    // For spill counting normalisation
    int nSpills = 0;
    int nSpillsTrue = 0;
    double lastSpill = 0.;
    // 1D ToF for background subtraction
    TH1D *hdtof1d = new TH1D(Form("hdtof1d_%d",nBlocks), Form("Time of flight, %d blocks; S4 - S1 / ns; Events", nBlocks), 260, 30, 160);

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

    // Need these to calculate bar efficiencies
    TH1D *hCoins = new TH1D(Form("hCoins_%d",nBlocks), Form("Bar coincidences + S_{1,2} coincidences, %d blocks; Bar; Events",nBlocks), 10, 0.5, 10.5);
    TH1D *hHits  = new TH1D(Form("hHits_%d",nBlocks), Form("PMT hits + S_{1,2} coincidences, %d blocks; Bar; Events",nBlocks), 10, 0.5, 10.5);
    TH1D *hEff   = new TH1D(Form("hEff_%d",nBlocks), Form("S4 bar efficiencies %d blocks; Bar; Events",nBlocks), 10, 0.5, 10.5);
    hHits->Sumw2();
    hCoins->Sumw2();

    // In this loop calculate the bar-by-bar efficiencies
    for (int itdc=0; itdc<2; itdc++) {
      // Need both the coincidence and the raw trees for this
      double tempUstof;
      double ustofNs;
      for (int irun=runMin; irun<runMax+1; irun++) {
	// Load input files
	TFile *tofCoinFile = new TFile(Form("%srun%d/DsTOFcoincidenceRun%d_tdc%d.root", dstofDir, irun, irun, itdc+1));
	TFile *tofFile     = new TFile(Form("%srun%d/DsTOFtreeRun%d_tdc%d.root", dstofDir, irun, irun, itdc+1));
	TTree *tofCoinTree = (TTree*)tofCoinFile->Get("tofCoinTree");
	TTree *tofTree = (TTree*)tofFile->Get("tofTree");
	RawDsTofCoincidence *tofCoin = NULL;
	RawDsTofHeader *tof = NULL;
	tofCoinTree->SetBranchAddress("tofCoin", &tofCoin);
	tofTree->SetBranchAddress("tof", &tof);
	// Create new TTree for ustof signals
	TTree *ustofTree = new TTree("ustofTree", "ustof");;
	ustofTree->SetDirectory(0);
	ustofTree->Branch("ustofNs", &ustofNs, "ustofNs/D");
	tempUstof=0;
      } // for (int irun=runMin; irun<runMax+1; irun++) 
    } // for (int itdc=0; itdc<2; itdc++)

  } // for (int nBlocks = 0; nBlocks < 4; nBlocks++) 
} // angularDistS4
