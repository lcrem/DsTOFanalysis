// angularDistS4.C
// Angular distribution of protons and pions for different moderator blocks
// Use efficiency calculation and background subtraction
#include "UsefulFunctions.C"

// Error on the weight going into the histogram
double weightErr(const double posEff, const double posEffErr, 
		 const double barEff, const double barEffErr)
{
  double err = 0.;
  err = TMath::Sqrt( pow(posEffErr/(posEff*posEff*barEff), 2) + pow(barEffErr/(barEff*barEff*posEff), 2) + pow(1. / (posEff*barEff), 2) );
  return err;
}

double weightErr(const double weight, const double werr)
{
  double err = 0.;
  err = TMath::Sqrt( pow(1. / weight, 2) + pow(werr/(weight*weight), 2) );
  return err;
}

double weightErrBar(const double barEff, const double barEffErr)
{
  double err = 0.;
  err = TMath::Sqrt(pow(barEffErr / (barEff * barEff), 2) + pow(1. / barEff, 2));
  return err;
}

double ratioErr(const double num, const double numErr, const double denom, const double denomErr)
{
  double err = (num / denom) * TMath::Sqrt(pow(numErr/num, 2) + pow(denomErr/denom, 2));
  return err;
}

// Main macro
void angularDistS4_newSample(const char* saveDir, 
			     const char* dstofDir="/nfs/scratch0/dbrailsf/data_backup/dtof_backup/",
			     const char* ustofDir="/nfs/scratch0/dbrailsf/data_backup/utof_backup_firsthitpinnedtounixtime/Data_root_v3_wo_walk_corr/",
			     const char* spillDBDir="/scratch0/sjones/spillDB/", 
			     const char* smearHists="/scratch0/sjones/plots/fitGaussToCosmics/newSigmoids.root") 
{ 
  gROOT->SetBatch(kTRUE);

  // S1 -> S4 positions
  const double baselineS1S4End   = 13.9426;
  const double baselineS1S4Start = 14.0069;
  const double s4OffAxisStartX   = 0.121;
  const double s4OffAxisEndX     = 1.4224;
  // Edges of S3 in beam coordinate system
  const double s3StartX   = -0.5168;
  const double s3EndX     = 0.9970;
  const double s3s1StartY = 9.0569 + 1.77;
  const double s3s1EndY   = 8.9146 + 1.77;
  // z coordinates
  const double s4BarTop    = 0.43; 
  const double s4BarBottom = -0.35;
  // Shift in ns required to to pion peak at speed of light
  const double dstofShift = 44.;
  // Ustof-dstof cable delay
  const double ustofDelay = 184.7;
  // Minimum number of hits in a spill required for it to count
  const int minHits = 200;
  // S2 to S4 distance - used for m^2 plot
  const double s2s4Dist = 12.65;
  // S1 to S4 distance
  const double s1s4Dist = 14.0698;
  // Bar efficiency factor in addition to cosmic data 
  const double barOverallEff = .8;
  // For calculating ratios
  const vector<double> dataS3    = {1983., 1656., 1325., 899., 136.3};
  const vector<double> dataS3Err = {9., 6., 5., 6., 0.5};
  // For MC plots
  const double xUp = 0.44; // m
  const double xDown = -0.97; // m 

  TFile *fout = new TFile(saveDir, "recreate");

  // THStack *hsBkgSub = new THStack("hsBkgSub", "S2 to S4 time of flight spectrum; t_{S4} - t_{S2} / ns; Events / spill");
  THStack *hsDtof   = new THStack("hsDtof", "S4 ToF spectrum; t_{S4} - t_{S2} / ns; Events / spill");

  THStack *hsProS4Vert = new THStack("hsProS4Vert", "S1 #cap S2 #cap S4 angular distribution of proton hits; #phi / degrees; Events / spill / degree");
  THStack *hsPiS4Vert  = new THStack("hsPiS4Vert", "S1 #cap S2 #cap S4 angular distribution of MIP hits; #phi / degrees; Events / spill / degree");
  THStack *hsProS4Horz = new THStack("hsProS4Horz", "S1 #cap S2 #cap S4 angular distribution of proton hits; #theta / degrees; Events / spill / degree");
  THStack *hsProS4HorzSmear = new THStack("hsProS4HorzSmear", "S1 #cap S2 #cap S4 angular distribution of proton hits; #theta / degrees; Events / spill / degree");
  THStack *hsPiS4Horz  = new THStack("hsPiS4Horz", "S1 #cap S2 #cap S4 angular distribution of MIP hits; #theta / degrees; Events / spill / degree");
  THStack *hsRatioS4Vert = new THStack("hsRatioS4Vert", "S1 #cap S2 #cap S4 angular distribution of proton/MIP ratio; #phi / degrees; Protons/MIPs");
  THStack *hsRatioS4Horz = new THStack("hsRatioS4Horz", "S1 #cap S2 #cap S4 angular distribution of proton/MIP ratio; #theta / degrees; Protons/MIPs");
  // No weights
  THStack *hsProS4VertUnwgt = new THStack("hsProS4VertUnwgt", "S1 #cap S2 #cap S4 angular distribution of proton hits; #phi / degrees; Events");
  THStack *hsPiS4VertUnwgt  = new THStack("hsPiS4VertUnwgt", "S1 #cap S2 #cap S4 angular distribution of MIP hits; #phi / degrees; Events");
  THStack *hsProS4HorzUnwgt = new THStack("hsProS4HorzUnwgt", "S1 #cap S2 #cap S4 angular distribution of proton hits; #theta / degrees; Events");
  THStack *hsPiS4HorzUnwgt  = new THStack("hsPiS4HorzUnwgt", "S1 #cap S2 #cap S4 angular distribution of MIP hits; #theta / degrees; Events");
  THStack *hsRatioS4VertUnwgt = new THStack("hsRatioS4VertUnwgt", "S1 #cap S2 #cap S4 angular distribution of proton/MIP ratio; #phi / degrees; Protons/MIPs");
  THStack *hsRatioS4HorzUnwgt = new THStack("hsRatioS4HorzUnwgt", "S1 #cap S2 #cap S4 angular distribution of proton/MIP ratio; #theta / degrees; Protons/MIPs");

  THStack *hsMom = new THStack("hsMom", "Proton momentum as measured in S4; Momentum [MeV/c]; Events / spill");
  THStack *hsMomS1S4 = new THStack("hsMomS1S4", "Proton momentum as measured in S4 (assuming time is S1-S4 time); Momentum [MeV/c]; Events / spill");

  THStack *hsEff = new THStack("hsEff", "S4 bar efficiencies; Bar number; Efficiency");
  THStack *hsEffRatio = new THStack("hsEffRatio", "S4 bar efficiency comparison; Bar number; Beam / cosmic");

  TLegend *leg = new TLegend(0.65, 0.55, 0.85, .85);

  TLegend *legProS4Horz = new TLegend(0.55, 0.55, 0.88, .85);
  TLegend *legPiS4Horz  = new TLegend(0.55, 0.55, 0.88, .85);
  TLegend *legProS4Vert = new TLegend(0.60, 0.55, 0.88, .85);
  TLegend *legPiS4Vert  = new TLegend(0.60, 0.55, 0.88, .85);

  TLegend *legRatioVert = new TLegend(0.15, 0.5, 0.4, 0.8);

  TFile *fEffHists = new TFile(smearHists, "read");
  /*
  double binsCosmics[] = {0., 7., 14., 21., 28., 35., 42., 49.,
			  52., 55., 58., 
			  61., 64., 67., 
			  70., 73., 76., 79., 
			  82., 85., 88., 
			  91., 98., 105., 112., 119., 126., 133., 140.};
  */
  double binsCosmics[] = {0., 7., 14., 21., 28., 35., 42., 49.,
			  51., 53., 55., 57., 59.,
			  61., 63., 65., 67., 69.,  
			  71., 73., 75., 77., 79., 
			  81., 83., 85., 87., 89., 
			  91., 98., 105., 112., 119., 126., 133., 140.};
  int binnum = sizeof(binsCosmics)/sizeof(double) - 1;

  for (int nBlocks = 0; nBlocks <= 3; nBlocks++) {
    cout<<"=========================================="<<endl;
    cout<<nBlocks<<" blocks"<<endl;
    cout<<"=========================================="<<endl;
    vector<double> startTimes;
    vector<double> endTimes;
    startTimes.clear();
    endTimes.clear();
    fout->cd();
    TTree *protonTree = new TTree(Form("protonTree%d", nBlocks), Form("protonTree%d", nBlocks));
    double weight;
    double tof;
    double mom;
    double theta;
    double phi;
    double mcx, mcy, mcz;
    int spill;
    protonTree->Branch("weight", &weight);
    protonTree->Branch("tof", &tof);
    protonTree->Branch("mom", &mom);
    protonTree->Branch("theta", &theta);
    protonTree->Branch("phi", &phi);
    protonTree->Branch("mcX", &mcx); 
    protonTree->Branch("mcY", &mcy); 
    protonTree->Branch("mcZ", &mcz); 
    protonTree->Branch("spill", &spill);
    int doubleHits = 0;

    THStack* hsEffComp = new THStack(Form("hsEffComp%d", nBlocks), Form("Efficiency comparison, %d blocks; Bar; Efficiency", nBlocks));
    THStack* hsEffComp2dNorm = new THStack(Form("hsEffComp2dNorm%d", nBlocks), Form("Efficiency comparison, %d blocks; Bar; Efficiency", nBlocks));
    // 1D histograms in the horizontal and vertical direction of S4
    // Horizontal
    TH1D *hProS4Horz = new TH1D(Form("hProS4Horz%d",nBlocks), Form("Horizontal angular distribution of proton hits in S4, %d blocks; #theta / degrees; Events / spill",nBlocks), 20, 0., 6.);
    setHistAttr(hProS4Horz);
    vector<double> proS4HorzErr;
    proS4HorzErr.resize(hProS4Horz->GetNbinsX()+2, 0);
    TH1D *hProS4HorzSmear = new TH1D(Form("hProS4HorzSmear%d",nBlocks), Form("Horizontal angular distribution of proton hits in S4, %d blocks; #theta / degrees; Events / spill",nBlocks), 20, 0., 6.);
    setHistAttr(hProS4HorzSmear);
    TH1D *hPiS4Horz  = new TH1D(Form("hPiS4Horz%d",nBlocks), Form("Horizontal angular distribution of MIP hits in S4, %d blocks; #theta / degrees; Events / spill",nBlocks), 20, 0., 6.);
    setHistAttr(hPiS4Horz);
    vector<double> piS4HorzErr;
    piS4HorzErr.resize(hPiS4Horz->GetNbinsX()+2, 0);
    TH1D *hAllS4Horz = new TH1D(Form("hAllS4Horz%d",nBlocks), Form("Horizontal angular distribution of hit in S4, %d blocks; #theta / degrees; Events / spill", nBlocks), 20, 0., 6.);
    setHistAttr(hAllS4Horz);
    TH1D *hProPiRatioS4Horz  = new TH1D(Form("hProPiRatioS4Horz%d",nBlocks), Form("Horizontal angular distribution of proton/MIP ratio in S4, %d blocks; #theta / degrees; Protons/MIPs",nBlocks), 20, 0., 6.);
    setHistAttr(hProPiRatioS4Horz);
    // Vertical
    TH1D *hProS4Vert = new TH1D(Form("hProS4Vert%d",nBlocks), Form("Vertical angular distribution of proton hits in S4, %d blocks; #phi / degrees; Events / spill",nBlocks), 10, -1.5, 1.8);
    setHistAttr(hProS4Vert);
    vector<double> proS4VertErr;
    proS4VertErr.resize(hProS4Vert->GetNbinsX()+2, 0);
    TH1D *hPiS4Vert  = new TH1D(Form("hPiS4Vert%d",nBlocks), Form("Vertical angular distribution of MIP hits in S4, %d blocks; #phi / degrees; Events / spill",nBlocks), 10, -1.5, 1.8);
    setHistAttr(hPiS4Vert);
    vector<double> piS4VertErr;
    piS4VertErr.resize(hPiS4Vert->GetNbinsX()+2, 0);
    TH1D *hProPiRatioS4Vert  = new TH1D(Form("hProPiRatioS4Vert%d",nBlocks), Form("Vertical angular distribution of proton/MIP ratio in S4, %d blocks; #phi / degrees; Protons/MIPs",nBlocks), 10, -1.5, 1.8);
    setHistAttr(hProPiRatioS4Vert);
    TH1D *hMSq = new TH1D(Form("hMSq%d", nBlocks), Form("Particle mass distribution, %d blocks; M^{2} [GeV^{2} / c^{2}]; Events / spills", nBlocks), 109, -0.3, 4.5);
    setHistAttr(hMSq);
    vector<double> mSqErr;
    mSqErr.resize(hMSq->GetNbinsX()+2, 0);
    TH2D *h2DoubleFirst = new TH2D(Form("h2DoubleFirst%d", nBlocks), Form("%d blocks: S4 double hits, 1st bar; Bar position / cm; Bar; Events / spill", nBlocks), 20, 0., 140., 10, 0.5, 10.5);
    setHistAttr(h2DoubleFirst);
    TH2D *h2DoubleSecond = new TH2D(Form("h2DoubleSecond%d", nBlocks), Form("%d blocks: S4 double hits, 2nd bar; Bar position / cm; Bar; Events / spill", nBlocks), 20, 0., 140., 10, 0.5, 10.5);
    setHistAttr(h2DoubleSecond);
    /*
    TH1D *hX = new TH1D(Form("hX%d", nBlocks), Form("%d blocks, data; x / m; Events", nBlocks), 50, -0.2, -1.7);
    setHistAttr(hX);
    TH1D *hZ = new TH1D(Form("hZ%d", nBlocks), Form("%d blocks, data; z / m; Events", nBlocks), 50, 12., 12.5);
    setHistAttr(hZ);
    */
    // No weighting
    TH1D *hProS4HorzUnwgt = new TH1D(Form("hProS4HorzUnwgt%d",nBlocks), Form("Horizontal angular distribution of proton hits in S4, %d blocks; #theta / degrees; Events",nBlocks), 20, 0., 6.);
    setHistAttr(hProS4HorzUnwgt);
    TH1D *hPiS4HorzUnwgt  = new TH1D(Form("hPiS4HorzUnwgt%d",nBlocks), Form("Horizontal angular distribution of MIP hits in S4, %d blocks; #theta / degrees; Events",nBlocks), 20, 0., 6.);
    setHistAttr(hPiS4HorzUnwgt);
    TH1D *hProS4VertUnwgt = new TH1D(Form("hProS4VertUnwgt%d",nBlocks), Form("Vertical angular distribution of proton hits in S4, %d blocks; #phi / degrees; Events",nBlocks), 10, -1.5, 1.8);
    setHistAttr(hProS4VertUnwgt);
    TH1D *hPiS4VertUnwgt  = new TH1D(Form("hPiS4VertUnwgt%d",nBlocks), Form("Vertical angular distribution of MIP hits in S4, %d blocks; #phi / degrees; Events",nBlocks), 10, -1.5, 1.8);
    setHistAttr(hPiS4VertUnwgt);
    TH1D *hProPiRatioS4HorzUnwgt  = new TH1D(Form("hProPiRatioS4HorzUnwgt%d",nBlocks), Form("Horizontal angular distribution of proton/MIP ratio in S4, %d blocks; #theta / degrees; Protons/MIPs",nBlocks), 20, 0., 6.);
    TH1D *hProPiRatioS4VertUnwgt  = new TH1D(Form("hProPiRatioS4VertUnwgt%d",nBlocks), Form("Vertical angular distribution of proton/MIP ratio in S4, %d blocks; #phi / degrees; Protons/MIPs",nBlocks), 10, -1.5, 1.8);
    TH2D *h2dAngProS1Unwgt = new TH2D(Form("h2dAngProS1Unwgt%d", nBlocks), Form("S4 angular distribution of proton hits, %d blocks; #theta / degrees; #phi / degrees; Events", nBlocks), 20, 0., 6., 10, -1.5, 1.8);
    setHistAttr(h2dAngProS1Unwgt);
    TH2D *h2dAngPiS1Unwgt = new TH2D(Form("h2dAngPiS1Unwgt%d", nBlocks), Form("S4 angular distribution of MIP hits, %d blocks; #theta / degrees; #phi / degrees; Events", nBlocks), 20, 0., 6., 10, -1.5, 1.8);
    setHistAttr(h2dAngPiS1Unwgt);
    TH2D *h2dAngRatioS1Unwgt = new TH2D(Form("h2dAngRatioS1Unwgt%d", nBlocks), Form("S4 angular distribution of proton/MIP ratio, %d blocks; #theta / degrees; #phi / degrees; Ratio", nBlocks), 20, 0., 6., 10, -1.5, 1.8);
    setHistAttr(h2dAngRatioS1Unwgt);

    // Proton momentum
    TH1D *hMom = new TH1D(Form("hMom%d", nBlocks), Form("%d blocks; Momentum [MeV / c]; Events / spill", nBlocks), 75, 100, 850);
    hMom->SetLineWidth(2);
    vector<double> momErr;
    momErr.resize(hMom->GetNbinsX()+2, 0);
    TH1D *hMomS1S4 = new TH1D(Form("hMomS1S4%d", nBlocks), Form("%d blocks; Momentum [MeV / c]; Events / spill", nBlocks), 75, 100, 850);
    hMomS1S4->SetLineWidth(2);

    TH1D *hSmearWeight = new TH1D(Form("hSmearWeight%d", nBlocks), Form("Smear weights, %d blocks", nBlocks), 20, xDown, xUp);
    setHistAttr(hSmearWeight);

    TH2D *h2dAngProS1 = new TH2D(Form("h2dAngProS1%d", nBlocks), Form("S4 angular distribution of proton hits, %d blocks; #theta / degrees; #phi / degrees; Events / spill", nBlocks), 20, 0., 6., 10, -1.5, 1.8);
    TH2D *h2dAngPiS1 = new TH2D(Form("h2dAngPiS1%d", nBlocks), Form("S4 angular distribution of MIP hits, %d blocks; #theta / degrees; #phi / degrees; Events / spill", nBlocks), 20, 0., 6., 10, -1.5, 1.8);
    TH2D *h2dAngRatioS1 = new TH2D(Form("h2dAngRatioS1%d", nBlocks), Form("S4 angular distribution of proton/MIP ratio, %d blocks; #theta / degrees; #phi / degrees; Events / spill", nBlocks), 20, 0., 6., 10, -1.5, 1.8);

    TH2D *h2dAngProWC = new TH2D(Form("h2dAngProWC%d", nBlocks), Form("S4 angular distribution of proton hits, %d blocks; #theta / degrees; #phi / degrees; Events / spill", nBlocks), 20, 0., 6., 10, -1.5, 1.8);
    TH2D *h2dAngPiWC = new TH2D(Form("h2dAngPiWC%d", nBlocks), Form("S4 angular distribution of MIP hits, %d blocks; #theta / degrees; #phi / degrees; Events / spill", nBlocks), 20, 0., 6., 10, -1.5, 1.8);
    TH2D *h2dAngRatioWC = new TH2D(Form("h2dAngRatioWC%d", nBlocks), Form("S4 angular distribution of proton/MIP ratio, %d blocks; #theta / degrees; #phi / degrees; Events / spill", nBlocks), 20, 0., 6., 10, -1.5, 1.8);
    TH2D *h2dProMCComp = new TH2D(Form("h2dProMCComp%d", nBlocks), Form("Proton distribution, %d blocks; x / m; y / m; Events / spill", nBlocks), 20, -0.9805, 0.4183, 10, -0.3406, 0.4434);
    setHistAttr(h2dProMCComp);
    TH2D *h2dProMCCompCut = new TH2D(Form("h2dProMCCompCut%d", nBlocks), Form("Proton distribution, %d blocks; x / m; y / m; Events / spill", nBlocks), 20, -0.8805, 0.3183, 9, -0.3406, 0.3684);
    setHistAttr(h2dProMCCompCut);
 
    // 1D ToF 
    TH1D *hdtof1d = new TH1D(Form("hdtof1d_%d",nBlocks), Form("Time of flight, measured in S4, %d blocks; t_{S4} - t_{S2} / ns; Events / spill", nBlocks), 182, 30, proCutHiS4);
    setHistAttr(hdtof1d);
    hdtof1d->SetLineWidth(2);
    vector<double> dtof1dErr;
    dtof1dErr.resize(hdtof1d->GetNbinsX()+2, 0);
    // Total efficiency
    TH1D *hEffTotal = new TH1D(Form("hEffTotal_%d", nBlocks), Form("With beam data, %d blocks; Bar; Events",nBlocks), 10, 0.5, 10.5);

    // Number of signal particles using just cut and count
    double nP  = 0.;
    double nPi = 0.;

    vector<TH1D*> smearHistVec;
    for (int bar=1; bar<10; bar++) {
      TH1D *hSmear = (TH1D*)fEffHists->Get(Form("hTrueDataRatio%dbar%d", nBlocks, bar));
      smearHistVec.push_back(hSmear);
    }
    // Define signal and background functions to be fitted
    // Signals are gaussians
    TF1 *sPro = new TF1(Form("sPro%d", nBlocks), "gaus", proFitLowS4.at(nBlocks), proFitHiS4.at(nBlocks));
    TF1 *sPi  = new TF1(Form("sPi%d", nBlocks), "gaus", piLowS4, piHiS4);
    // Exponential background
    // TF1 *fBkgExp = new TF1(Form("fBkgExp%d", nBlocks),"expo", 30, proCutHiS4);
    TF1 *fBkg = new TF1(Form("fBkg%d", nBlocks),"pol0", 30, proCutHiS4);
    sPro->SetLineColor(kGreen+2);
    sPi->SetLineColor(kRed);
    /*
      TF1 *fSplusBExp = new TF1(Form("signal_plus_bkg_exp_%d", nBlocks), "gaus(0)+gaus(3)+expo(6)", 30, proCutHiS4);
      fSplusBExp->SetParNames("const 1", "mean 1", "sigma 1",
      "const 2", "mean 2", "sigma 2",
      "bkgconst", "bkgdecay");
      fSplusBExp->SetLineColor(kBlack);
    */
    TF1 *fSplusB = new TF1(Form("signal_plus_bkg_%d", nBlocks), "gaus(0)+gaus(3)+pol0(6)", 30, proCutHiS4);
    fSplusB->SetParNames("piConst", "piMean", "piSigma",
			 "proConst", "proMean", "proSigma",
			 "bkg");
    fSplusB->SetLineColor(kRed);
    // For spill counting normalisation
    int nSpills = 0;
    int nSpillsTrue = 0;
    double lastSpill = 0.;

    // Find the correct dstof files
    Int_t runMin=-1;
    Int_t runMax=-1;
    double startTime = 0;
    double endTime   = 0;

    if (nBlocks == 0) {
      startTimes.push_back(start0Block);
      endTimes.push_back(end0Block);
    }
    else if (nBlocks == 1) {
      startTimes.push_back(start1Block);
      endTimes.push_back(end1Block);
    }
    else if (nBlocks == 2) {
      startTimes.push_back(start2Block);
      endTimes.push_back(end2Block);
    }
    else if (nBlocks == 3) {
      startTimes.push_back(start3Block);
      endTimes.push_back(end3Block);
    }
    else if (nBlocks == 4) {
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

    // Loop over subsamples
    for (int sub=0; sub<startTimes.size(); sub++) {
      startTime = startTimes.at(sub);
      endTime   = endTimes.at(sub);
      // Cosmic hists -- need to be remade for each subsample
      TH2D *h2Cosmics = new TH2D(Form("h2Cosmics%d",nBlocks), Form("Cosmic flux, %d blocks; x / cm; Bar; Rate / s^{-1}", nBlocks), 20, 0, 140, 10, 0.5, 10.5);
      setHistAttr(h2Cosmics);
      TH2D *h2CosmicsEff = new TH2D(Form("h2CosmicsEff%d", nBlocks), Form("Relative efficiency, %d blocks; x / cm; Bar; Eff", nBlocks), 20, 0, 140, 10, 0.5, 10.5);
      setHistAttr(h2CosmicsEff);
      TH2D *h2CosmicsEffBins = new TH2D(Form("h2CosmicsEffBins%d", nBlocks), Form("Relative efficiency, %d blocks; x / cm; Bar; Eff", nBlocks), binnum, binsCosmics, 10, 0.5, 10.5);
      setHistAttr(h2CosmicsEffBins);
      TH2D *h2CosmicsEffMC = new TH2D(Form("h2CosmicsEffMC%d", nBlocks), Form("Relative efficiency, %d blocks; x / cm; y / cm; Eff", nBlocks), 20, -0.95, 0.41, 10, -0.3361, 0.4139);
      setHistAttr(h2CosmicsEffMC);
      TH1D *hCosmicsVertEff = new TH1D(Form("hCosmicsVertEff%d",nBlocks), Form("With cosmics, %d blocks; x / cm; Eff",nBlocks), 10, 0.5, 10.5);
      setHistAttr(hCosmicsVertEff);
      TH1D *hCosmicsVertEff2dNorm = new TH1D(Form("hCosmicsVertEff2dNorm%d",nBlocks), Form("With cosmics, %d blocks; x / cm; Eff",nBlocks), 10, 0.5, 10.5);
      setHistAttr(hCosmicsVertEff2dNorm);
      TH1D *hCosmicsHorz = new TH1D(Form("hCosmicsHorz%d",nBlocks), Form("Cosmic flux, %d blocks; x / cm; Rate / Hz",nBlocks), 20, 0, 140);
      setHistAttr(hCosmicsHorz);
      // Bar by bar efficiencies calculated by cosmics
      std::vector<TH1D*> eff1dVec;
      for (int b=0; b<10; b++) {
	TH1D *hEff1d = new TH1D(Form("hEff1d_block%d_bar%d_sub%d", nBlocks, b+1, sub), Form("%d blocks, bar %d: Efficiency; x / cm; Eff", nBlocks, b+1), 20, 0., 140.);
	setHistAttr(hEff1d);
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

      // Need these to calculate bar efficiencies
      TH1D *hCoins = new TH1D(Form("hCoins_%d_%d", nBlocks, sub), Form("Bar coincidences + S_{1,2} coincidences, %d blocks; Bar; Events",nBlocks), 10, 0.5, 10.5);
      TH1D *hHits  = new TH1D(Form("hHits_%d_%d", nBlocks, sub), Form("PMT hits + S_{1,2} coincidences, %d blocks; Bar; Events",nBlocks), 10, 0.5, 10.5);
      TH1D *hEff   = new TH1D(Form("hEff_%d_%d", nBlocks, sub), Form("With beam data, %d blocks; Bar; Events",nBlocks), 10, 0.5, 10.5);
      hHits->Sumw2();
      hCoins->Sumw2();

      // In this loop calculate the bar-by-bar efficiencies and the angular efficiencies
      double lastSpillSignal = 0.;
      int nSpillsCosmics = 0;
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
	  TTree *ustofTree = new TTree("ustofTree", "ustof");
	  ustofTree->SetDirectory(0);
	  ustofTree->Branch("ustofNs", &ustofNs, "ustofNs/D");
	  tempUstof=0;
	  // Build tree of all the ustof hit times
	  for (int i=0; i<tofTree->GetEntries(); i++) {
	    tofTree->GetEntry(i);
	    if (tof->unixTime < startTime) continue;
	    if (tof->unixTime > endTime) break;

	    if (tof->channel == 13) {
	      ustofNs = tof->fakeTimeNs - ustofDelay;
	      // ustof signals shouldn't be coming closer than 500ns
	      if ( (ustofNs - tempUstof) < 500.) continue;
	      ustofTree->Fill();
	      tempUstof = ustofNs;
	    } // if (tof->channel == 13) 
	  } // for (int i=0; i<tof->GetEntries(); i++)
	  ustofTree->BuildIndex("ustofNs");
	  // Now loop over coincidence file to find number of coincidences for each bar
	  for (int t=0; t<tofCoinTree->GetEntries(); t++) {
	    tofCoinTree->GetEntry(t);
	    if (tofCoin->unixTime[0] < startTime) continue;
	    if (tofCoin->unixTime[0] > endTime) break;
	    if (abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1]) > s4BarTime) continue;
	    double deltat = TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1]);
	    double dstofHitT = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (s4BarTime/2. - TMath::Abs(deltat) / 2. );
	    double tofCalc = dstofHitT - tofCoin->usTofSignal;
	    // Count number of spills for time normalisation
	    if (tofCoin->lastDelayedBeamSignal != lastSpillSignal && itdc==1) {
	      nSpillsCosmics++;
	      lastSpillSignal = tofCoin->lastDelayedBeamSignal;
	    }

	    if (tofCalc > 70. && tofCalc < 200./* && tofCoin->bar != 10*/) {
	      hCoins->Fill(tofCoin->bar);
	    } // if (tofCalc > 70. && tofCalc < 200.)
	      // If the hits are not in a spill then consider them cosmics 
	      // and use them for angular efficiency
	    if (!tofCoin->inSpill) {
	      if (abs(tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0]) < s4BarTime) {
		double positionX = localDtofPosition(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]);
		double positionZ = (tofCoin->bar * 0.075) - 0.375 - 0.01;
		double mcX = -positionX/100. + 0.491;
		double mcY = positionZ + 0.0114;
		eff1dVec.at(tofCoin->bar-1)->Fill(positionX);
		h2Cosmics->Fill(positionX, tofCoin->bar);
		h2CosmicsEff->Fill(positionX, tofCoin->bar);
		h2CosmicsEffBins->Fill(positionX, tofCoin->bar);
		h2CosmicsEffMC->Fill(mcX, mcY);
		hCosmicsVertEff->Fill(tofCoin->bar);
		hCosmicsVertEff2dNorm->Fill(tofCoin->bar);
	      }
	    } // Is not in a spill
	  } // for (int t=0; t<tofCoinTree->GetEntries(); t++) 
	  // Loop over hit file to find number of hits in coincidence 
	  for (int t=0; t<tofTree->GetEntries(); t++) {
	    tofTree->GetEntry(t);
	    if (tof->unixTime < startTime) continue;
	    if (tof->unixTime > endTime) break;
	    if (tof->channel > 10 || tof->channel == 0) continue;
	    int ustofEntry = ustofTree->GetEntryNumberWithBestIndex(tof->fakeTimeNs);
	    ustofTree->GetEntry(ustofEntry);

	    if (tof->fakeTimeNs - ustofNs < 260. && tof->fakeTimeNs - ustofNs > 130.) {
	      if (itdc == 0) { // TDC1
		if (tof->channel % 2 == 1) { // PMT B
		  hHits->Fill(((tof->channel + 1)/ -2) + 11);
		} // if (tof->channel % 2 == 1)
		else { // PMT A
		  hHits->Fill((tof->channel / -2) + 11);
		}
	      } // TDC 1
	      else { // TDC 2
		if (tof->channel % 2 == 1) { // PMT B
		  hHits->Fill(((tof->channel + 1)/ -2) + 6);
		} // if (tof->channel % 2 == 1)
		else { // PMT A
		  hHits->Fill((tof->channel / -2) + 6);
		}
	      } // TDC 2
	    } // if (tof->fakeTimeNs - ustofNs < 260. && tof->fakeTimeNs - ustofNs > 130.)
	  } // for (int t=0; t<tofTree->GetEntries(); t++)
	  delete tof;
	  delete tofCoin;
	  delete tofTree;
	  delete tofCoinTree;
	  delete ustofTree;
	  tofFile->Close();
	  tofCoinFile->Close();
	} // for (int irun=runMin; irun<runMax+1; irun++) 
      } // for (int itdc=0; itdc<2; itdc++)
      cout<<"Created efficiency hists"<<endl;
      fout->cd(); 
      // for (int i=0; i<h2Cosmics->GetNbinsX()+1; i++) {
      // 	for (int j=0; j<h2Cosmics->GetNbinsY()+1; j++) {
      // 	  h2CosmicsEff->SetBinContent(i, j, h2Cosmics->GetBinContent(i, j)/h2Cosmics->GetBinContent(h2Cosmics->GetMaximumBin()));
      // 	}
      // }
      // Bar-by-bar angular efficiencies
      for (int bar=0; bar<eff1dVec.size(); bar++) {
	eff1dVec.at(bar)->Scale(1./eff1dVec.at(bar)->GetBinContent(eff1dVec.at(bar)->GetMaximumBin()));
	eff1dVec.at(bar)->Write();
      }
      // Scale by largest bin in the 2d histogram (with appropriate area normalisation)
      hCosmicsVertEff->Scale(1./hCosmicsVertEff->GetBinContent(hCosmicsVertEff->GetMaximumBin()));
      h2CosmicsEff->Scale(1./((endTime - startTime - nSpillsCosmics)*maxCosmics /*h2CosmicsEff->GetBinContent(h2CosmicsEff->GetMaximumBin())*/));
      h2CosmicsEffBins->Scale(1., "width");
      h2CosmicsEffBins->Scale(1./h2CosmicsEffBins->GetBinContent(h2CosmicsEffBins->GetMaximumBin()));
      h2CosmicsEffMC->Scale(1./h2CosmicsEffMC->GetBinContent(h2CosmicsEffMC->GetMaximumBin()));
      hCosmicsVertEff2dNorm->Scale(1./(h2Cosmics->GetBinContent(h2Cosmics->GetMaximumBin())*h2Cosmics->GetNbinsX()));
      h2Cosmics->Scale(1. / (endTime - startTime - nSpillsCosmics));
      hHits->Write();
      hCoins->Write();
      hHits->Add(hCoins, -1.);
      hEff->Divide(hCoins, hHits, 1., 1., "B");
      h2CosmicsEff->Write();
      h2CosmicsEffBins->Write();
      h2CosmicsEffMC->Write();
      h2Cosmics->Write();
      // Now loop over the coincidence files again and calculate the angular distributions
      cout<<"Getting signal hits"<<endl;
      for (int itdc=0; itdc<2; itdc++) {
	// First of all tchain the relevant files together
	TChain *tofCoinChain = new TChain("tofCoinTree");
	for (int irun=runMin; irun<runMax+1; irun++){
	  tofCoinChain->Add(Form("%srun%d/DsTOFcoincidenceRun%d_tdc%d.root", dstofDir, irun, irun, itdc+1));
	} // for (int irun=runMin; irun<runMax+1; irun++)
	cout<<"Got TChain successfully"<<endl;
	RawDsTofCoincidence *tofCoin = NULL;
	tofCoinChain->SetBranchAddress("tofCoin", &tofCoin);
	cout<<"Branch address set successfully"<<endl;
	int lasth = 0;
	// Use the spills recorded in the spill DB
	for (int h=0; h<tofCoinChain->GetEntries(); h++) {
	  tofCoinChain->GetEntry(h);
	  if (tofCoin->unixTime[0]<startTime) continue;
	  if (tofCoin->unixTime[0]>endTime) break;
	  if (tofCoin->lastDelayedBeamSignal != lastSpill && itdc == 0) {
	    lastSpill = tofCoin->lastDelayedBeamSignal;
	    nSpills++;
	    int nTrueHits = 0;
	    for (int sp=h; sp<tofCoinChain->GetEntries(); sp++) {
	      tofCoinChain->GetEntry(sp);
	      if (tofCoin->unixTime[0]<startTime) continue;
	      if (tofCoin->unixTime[0]>endTime) break;
	      if (tofCoin->fakeTimeNs[0] > lastSpill + 3e9) break;

	      if ((tofCoin->fakeTimeNs[0] - tofCoin->usTofSignal) < 200. &&
		  (tofCoin->fakeTimeNs[0] - tofCoin->usTofSignal) > 70. &&
		  (tofCoin->fakeTimeNs[0] - lastSpill) < 1e9 &&
		  (tofCoin->fakeTimeNs[0] - lastSpill) > 0.) {
		nTrueHits++;
		if (nTrueHits > 25) {
		  nSpillsTrue++;
		  break;
		}
	      }
	    } // for (int sp=h; sp<tofCoinChain->GetEntries(); sp++) 
	  } // if (tofCoin->lastDelayedBeamSignal != lastSpill && itdc == 0) 
	  tofCoinChain->GetEntry(h);
	  // Need to calculate total signal hits here
	  // Weight these by the efficiency of calculated above
	  double deltat = TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1]);
	  double dstofHitT = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (s4BarTime/2. - TMath::Abs(deltat) / 2.);
	  double tofCalc = dstofHitT - tofCoin->usTofSignal - dstofShift;
	  if (tofCalc < proCutHiS4 && tofCalc > 30. && tofCoin->bar != 10 && deltat < s4BarTime) {
	    double positionXP = localDtofPosition(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]);
	    TVector3 globalCoords = GetDtofGlobalCoords(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1], tofCoin->bar);
	    double mcXForSmear = globalToMCCoords(globalCoords).X();
	    double smearWeight = smearHistVec.at(tofCoin->bar-1)->GetBinContent(smearHistVec.at(tofCoin->bar-1)->GetXaxis()->FindBin(mcXForSmear));
	    hSmearWeight->Fill(mcXForSmear, smearWeight);

	    if (smearWeight > 40) smearWeight = 40.;
	    double w = 1. / (h2CosmicsEff->GetBinContent(h2CosmicsEff->GetXaxis()->FindBin(positionXP), tofCoin->bar)*barOverallEff); 
	    // Need to count if there is a double hit associated with this one
	    int bar1 = tofCoin->bar;
	    for (int dub = h-1; dub<=h+1; dub+=2) {
	      tofCoinChain->GetEntry(dub);
	      double dstofHitT2 = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (s4BarTime/2. - TMath::Abs(deltat) / 2.);
	      int bar2 = tofCoin->bar;
	      if (dstofHitT2 - dstofHitT < 2. && dstofHitT2 > dstofHitT && abs(bar1-bar2)==1) {
		h2DoubleFirst->Fill(positionXP, bar1);
		double positionXP2 = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(70./s4BarTime) + 70.));
		h2DoubleSecond->Fill(positionXP2, bar2);
		doubleHits++;
		// If there is a double hit we should halve the weight for this event
		w /= 2;
		break;
	      }
	    } // Loop to find double hits
	    double errSq = pow(weightErr(h2CosmicsEff->GetBinContent(h2CosmicsEff->GetXaxis()->FindBin(positionXP), tofCoin->bar)*barOverallEff/**(1./smearHistVec.at(tofCoin->bar-1)->GetBinContent(smearHistVec.at(tofCoin->bar-1)->GetXaxis()->FindBin(positionXP)))*/, h2CosmicsEff->GetBinError(h2CosmicsEff->GetXaxis()->FindBin(positionXP), tofCoin->bar)), 2);
	    double errSqBar = pow(weightErrBar(hEff->GetBinContent(tofCoin->bar), hEff->GetBinError(tofCoin->bar)), 2);
	    // Calculate position of hit in global coordinates
	    TVector3 dtofCoords = GetDtofGlobalCoords(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1], tofCoin->bar);
	    double positionX = dtofCoords.X();
	    double positionY = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(70./s4BarTime) + 65.) / 130.) * (baselineS1S4End - baselineS1S4Start) + baselineS1S4Start;
	    double positionZ = dtofCoords.Y();
	    // Calculate the angles relative to the nominal beamline
	    double angleTheta = TMath::ATan(positionX / positionY) * (180./TMath::Pi());
	    double anglePhi   = TMath::ATan(positionZ / positionY) * (180./TMath::Pi());
	    // Apply offsets to get these in same frame as MC
	    TVector3 mcCoords = globalToMCCoords(dtofCoords);
	    double mcX = mcCoords.X();
	    double mcY = mcCoords.Y();
	    double mcZ = mcCoords.Z();
	    hdtof1d->Fill(tofCalc, w);
	    hMSq->Fill(massFromTime(tofCalc, 0.8, s2s4Dist), w);
	    dtof1dErr.at(hdtof1d->GetXaxis()->FindBin(positionXP)) += errSq;
	    mSqErr.at(hMSq->GetXaxis()->FindBin(massFromTime(tofCalc, 0.8, s2s4Dist))) += errSq;

	    if (tofCalc < piHiS4 & tofCalc > piLowS4) { 
	      nPi += w;
	      hPiS4Horz->Fill(angleTheta, w);
	      hPiS4Vert->Fill(anglePhi, w);
	      h2dAngPiS1->Fill(angleTheta, anglePhi, w);

	      hPiS4HorzUnwgt->Fill(angleTheta);
	      hPiS4VertUnwgt->Fill(anglePhi);
	      h2dAngPiS1Unwgt->Fill(angleTheta, anglePhi);

	      piS4HorzErr.at(hPiS4Horz->GetXaxis()->FindBin(angleTheta)) += errSq;
	      piS4VertErr.at(hPiS4Vert->GetXaxis()->FindBin(anglePhi))   += errSq;
	    } // if (tofCalc < piHiS4 & tofCalc > piLowS4)
	    else if (tofCalc < proCutHiS4 & tofCalc > proCutLowS4) {
	      tof = tofCalc;
	      mom = momFromTime(0.938, s2s4Dist, tofCalc);
	      weight = w;
	      theta = angleTheta;
	      phi = anglePhi;
	      mcx = mcX;
	      mcy = mcY;
	      mcz = mcZ;
	      spill = nSpillsTrue;
	      protonTree->Fill();

	      nP += w;
	      hProS4Horz->Fill(angleTheta, w);
	      hProS4HorzSmear->Fill(angleTheta, w*smearWeight);
	      hProS4Vert->Fill(anglePhi, w);
	      hProS4HorzUnwgt->Fill(angleTheta);
	      hProS4VertUnwgt->Fill(anglePhi);
	      h2dAngProS1->Fill(angleTheta, anglePhi, w);
	      h2dAngProS1Unwgt->Fill(angleTheta, anglePhi);
	      hMom->Fill(momFromTime(0.938, s2s4Dist, tofCalc), w);
	      hMomS1S4->Fill(momFromTime(0.938, s1s4Dist, tofCalc+4.7), w);
	      h2dProMCComp->Fill(mcX, mcY, w);
	      if (positionXP > 10. && positionXP < 130.) {
		h2dProMCCompCut->Fill(mcX, mcY, w);
	      }
	      proS4HorzErr.at(hProS4Horz->GetXaxis()->FindBin(angleTheta)) += errSq;
	      proS4VertErr.at(hProS4Vert->GetXaxis()->FindBin(anglePhi))   += errSq;
	      momErr.at(hMom->GetXaxis()->FindBin(momFromTime(0.938, s2s4Dist, tofCalc))) += errSq;
	    } // else if (tofCalc < proHi & tofCalc > proLow) 
	  } // if (tofCalc < 160. && tofCalc > 30.) 
	} // for (int h=0; h<tofCoinChain->GetEntries(); h++)
	delete tofCoin;
	delete tofCoinChain;
      } // Loop over TDCs
      cout<<"Finished getting signal hits"<<endl;
      cout<<"Found "<<doubleHits<<" double proton hits"<<endl;;
      TH1D *hEffRatio = new TH1D(Form("hEffRatio%d", nBlocks), Form("Beam data efficiency / Cosmic efficiency, %d blocks; Bar, Beam / Cosmic", nBlocks), 10, 0.5, 10.5);
      setHistAttr(hEffRatio);
      hEffRatio->Divide(hEff, hCosmicsVertEff, 1., 1., "B");
      hEffRatio->Write();
      TH1D *hEffRatio2dNorm = new TH1D(Form("hEffRatio2dNorm%d", nBlocks), Form("Beam data efficiency / Cosmic efficiency, %d blocks; Bar, Beam / Cosmic", nBlocks), 10, 0.5, 10.5);
      setHistAttr(hEffRatio2dNorm);
      hEffRatio2dNorm->Divide(hEff, hCosmicsVertEff2dNorm, 1., 1., "B");
      hEffRatio2dNorm->Write();
      cout<<"Finished subsample "<<sub<<endl;
    } // Loop over subsamples

    for (int bin=0; bin < hdtof1d->GetNbinsX()+1; bin++) {
      hdtof1d->SetBinError(bin, TMath::Sqrt(dtof1dErr.at(bin)));
    }
    for (int bin=0; bin < hMSq->GetNbinsX()+1; bin++) {
      hMSq->SetBinError(bin, TMath::Sqrt(mSqErr.at(bin)));
    }
    for (int bin=0; bin < hProS4Horz->GetNbinsX()+1; bin++) {
      hProS4Horz->SetBinError(bin, TMath::Sqrt(proS4HorzErr.at(bin)));
      hPiS4Horz->SetBinError(bin, TMath::Sqrt(piS4HorzErr.at(bin)));
    }
    for (int bin=0; bin < hProS4Vert->GetNbinsX()+1; bin++) {
      hProS4Vert->SetBinError(bin, TMath::Sqrt(proS4VertErr.at(bin)));
      hPiS4Vert->SetBinError(bin, TMath::Sqrt(piS4VertErr.at(bin)));
    }
    for (int bin=0; bin<hMom->GetNbinsX()+1; bin++) {
      hMom->SetBinError(bin, TMath::Sqrt(momErr.at(bin)));
    }

    hdtof1d->Scale(1. / nSpillsTrue);
    fout->cd();
    hdtof1d->Fit(sPi, "R");
    hdtof1d->Fit(sPro, "R");
    hdtof1d->Fit(fBkg, "R");
    Double_t par[7];
    sPi->GetParameters(&par[0]);
    sPro->GetParameters(&par[3]);
    fBkg->GetParameters(&par[6]);
    fSplusB->SetParameters(par);
    hdtof1d->Fit(fSplusB, "R");
    hdtof1d->Write();
    sPro->Write();
    sPi->Write();
    fBkg->Write();
    fSplusB->Write();

    hProS4Horz->Scale(1. / nSpillsTrue);
    hProS4HorzSmear->Scale(1. / nSpillsTrue);
    hPiS4Horz->Scale(1. / nSpillsTrue);
    hProS4Vert->Scale(1. / nSpillsTrue);
    hPiS4Vert->Scale(1. / nSpillsTrue);
    hSmearWeight->Scale(1. / nSpillsTrue);
    hProS4Horz->Scale(1., "width");
    hProS4HorzSmear->Scale(1., "width");
    hPiS4Horz->Scale(1., "width");
    hProS4Vert->Scale(1., "width");
    hPiS4Vert->Scale(1., "width");
    hSmearWeight->Scale(1., "width");
    /*
      hProS4HorzUnwgt->Scale(1. / nSpillsTrue);
      hPiS4HorzUnwgt->Scale(1. / nSpillsTrue);
      hProS4VertUnwgt->Scale(1. / nSpillsTrue);
      hPiS4VertUnwgt->Scale(1. / nSpillsTrue);
      hProS4HorzUnwgt->Scale(1., "width");
      hPiS4HorzUnwgt->Scale(1., "width");
      hProS4VertUnwgt->Scale(1., "width");
      hPiS4VertUnwgt->Scale(1., "width");
    */

    // Now we have the fit values, loop over again and subtract the background
    // For each bin, find the fraction of each particle type which are background and 
    // this fraction from the bin
    TF1 *fSub = new TF1("fSub", "pol0", 30, proCutHiS4);
    fSub->SetParameter(0, fSplusB->GetParameter("bkg"));  
    // Background subtracted tof spectrum
    /*
      TH1D *hdtof1d_sub = (TH1D*)hdtof1d->Clone(Form("hdtof1d_sub_%d",nBlocks));
      hdtof1d_sub->Add(fSub, -1.);
      // If the bin content drops below 0, set to 0
      for (int b = 0; b <=  hdtof1d_sub->GetNbinsX(); b++) {
      if (hdtof1d_sub->GetBinContent(b) < 0.) {
      hdtof1d_sub->SetBinContent(b, 0);
      hdtof1d_sub->SetBinError(b, 0);
      } // if (hdtof1d_sub->GetBinContent(b) < 0.)
      } // for (int b = 0; hdtof1d_sub->GetNbinsX(); b++)
    */
    // Integrate background function between proton and pion windows and then subtract
    double bkgPerBin = fSplusB->GetParameter("bkg");
    double bkgPerNs  = bkgPerBin / hdtof1d->GetBinWidth(5);
    cout<<"Bkg per bin, per ns "<<bkgPerBin<<", "<<bkgPerNs<<endl;
    double piBkg  = fSub->Integral(piLowS4,  piHiS4) / hdtof1d->GetBinWidth(5);
    double proBkg = fSub->Integral(proCutLowS4, proCutHiS4) / hdtof1d->GetBinWidth(5);
    cout<<"Pion background: "<<piBkg<<" per spill. Proton background: "<<proBkg<<" per spill"<<endl;
    // Need to account for bin width
    double piBkgHorz  = piBkg / hPiS4Horz->GetBinWidth(3);
    double proBkgHorz = proBkg / hProS4Horz->GetBinWidth(3);
    double piBkgVert  = piBkg / hPiS4Vert->GetBinWidth(3);
    double proBkgVert = proBkg / hProS4Vert->GetBinWidth(3);
    cout<<"Pion background vert: "<<piBkgVert<<" / spill / degree. Proton background vert "<<proBkgVert<<" / spill / degree"<<endl;
    cout<<"Pion background horz: "<<piBkgHorz<<" / spill / degree. Proton background horz "<<proBkgHorz<<" / spill / degree"<<endl;
      
    cout<<"Spills "<<nSpills<<" ("<<nSpillsTrue<<" true)"<<endl;

    // Subtract background hits
    for (int b=1; b < hProS4Horz->GetNbinsX()+1; b++) {
      hProS4Horz->SetBinContent(b, hProS4Horz->GetBinContent(b) * (1 - proBkgHorz/nP));
      if (hProS4Horz->GetBinContent(b) < 0.) {
    	hProS4Horz->SetBinContent(b, 0);
    	hProS4Horz->SetBinError(b, 0);
      }
    }

    for (int b=1; b < hPiS4Horz->GetNbinsX()+1; b++) {
      hPiS4Horz->SetBinContent(b, hPiS4Horz->GetBinContent(b) * (1 - piBkgHorz/nPi));
      if (hPiS4Horz->GetBinContent(b) < 0.) {
    	hPiS4Horz->SetBinError(b, 0);
    	hPiS4Horz->SetBinContent(b, 0);
      }
    }

    for (int b=1; b < 10; b++) {
      hProS4Vert->SetBinContent(b, hProS4Vert->GetBinContent(b) * (1 - proBkgVert/nP));
      if (hProS4Vert->GetBinContent(b) < 0.) {
    	hProS4Vert->SetBinError(b, 0);
    	hProS4Vert->SetBinContent(b, 0);
      }
    }

    for (int b=1; b < 10; b++) {
      hPiS4Vert->SetBinContent(b, hPiS4Vert->GetBinContent(b) * (1 - piBkgVert/nPi));
      if (hPiS4Vert->GetBinContent(b) < 0.) {
    	hPiS4Vert->SetBinError(b, 0);
    	hPiS4Vert->SetBinContent(b, 0);
      }
    }

    hProPiRatioS4Vert->Divide(hProS4Vert, hPiS4Vert, 1., 1., "B");
    hProPiRatioS4Horz->Divide(hProS4Horz, hPiS4Horz, 1., 1., "B");
    hProPiRatioS4VertUnwgt->Divide(hProS4VertUnwgt, hPiS4VertUnwgt, 1., 1., "B");
    hProPiRatioS4HorzUnwgt->Divide(hProS4HorzUnwgt, hPiS4HorzUnwgt, 1., 1., "B");

    // hdtof1d_sub->SetLineWidth(2);
    hdtof1d->SetLineWidth(2);
    hProPiRatioS4Horz->SetLineWidth(2);
    hProPiRatioS4Vert->SetLineWidth(2);
    hProS4Horz->SetLineWidth(2);
    hProS4Vert->SetLineWidth(2);
    hPiS4Horz->SetLineWidth(2);
    hPiS4Vert->SetLineWidth(2);
    hProPiRatioS4HorzUnwgt->SetLineWidth(2);
    hProPiRatioS4VertUnwgt->SetLineWidth(2);
    hProS4HorzUnwgt->SetLineWidth(2); 
    hProS4VertUnwgt->SetLineWidth(2);
    hPiS4HorzUnwgt->SetLineWidth(2);
    hPiS4VertUnwgt->SetLineWidth(2);

    double eProS4Horz = 0;
    double ePiS4Horz  = 0;
    double intProS4Horz = hProS4Horz->IntegralAndError(1, hProS4Horz->GetNbinsX(), eProS4Horz, "width");
    double intPiS4Horz  = hPiS4Horz->IntegralAndError(1, hProS4Horz->GetNbinsX(), ePiS4Horz, "width");

    if (nBlocks == 0) {
      // hdtof1d_sub->SetLineColor(kBlack);
      hdtof1d->SetLineColor(kBlack);
      hProPiRatioS4Horz->SetLineColor(kBlack);
      hProPiRatioS4Vert->SetLineColor(kBlack);
      hProS4Horz->SetLineColor(kBlack);
      hProS4HorzSmear->SetLineColor(kBlack);
      hProS4Vert->SetLineColor(kBlack);
      hPiS4Horz->SetLineColor(kBlack);
      hPiS4Vert->SetLineColor(kBlack);
      hProPiRatioS4HorzUnwgt->SetLineColor(kBlack);
      hProPiRatioS4VertUnwgt->SetLineColor(kBlack);
      hProS4HorzUnwgt->SetLineColor(kBlack);
      hProS4VertUnwgt->SetLineColor(kBlack);
      hPiS4HorzUnwgt->SetLineColor(kBlack);
      hPiS4VertUnwgt->SetLineColor(kBlack);
      hEffTotal->SetLineColor(kBlack);
      hMSq->SetLineColor(kBlack);
      hMom->SetLineColor(kBlack);
      hMomS1S4->SetLineColor(kBlack);
      // hCosmicsVertEff->SetLineColor(kBlack);
      // hCosmicsVertEff2dNorm->SetLineColor(kBlack);
      // hEffRatio->SetLineColor(kBlack);
      // hEffRatio2dNorm->SetLineColor(kBlack);
      leg->AddEntry(hdtof1d, "0 blocks", "l");     
      legProS4Horz->AddEntry(hProS4Horz, Form("0 blocks - %.3g #pm %.1g per spill", intProS4Horz, eProS4Horz), "l");     
      legPiS4Horz->AddEntry(hPiS4Horz, Form("0 blocks - %d #pm %d per spill", (int)intPiS4Horz, (int)ePiS4Horz), "l");
      legRatioVert->AddEntry(hProPiRatioS4Vert, "0 blocks", "l");
    }
    else if (nBlocks == 1) {
      // hdtof1d_sub->SetLineColor(kRed);
      hdtof1d->SetLineColor(kRed);
      hProPiRatioS4Horz->SetLineColor(kRed);
      hProPiRatioS4Vert->SetLineColor(kRed);
      hProS4Horz->SetLineColor(kRed);
      hProS4HorzSmear->SetLineColor(kRed);
      hProS4Vert->SetLineColor(kRed);
      hPiS4Horz->SetLineColor(kRed);
      hPiS4Vert->SetLineColor(kRed);
      hProPiRatioS4HorzUnwgt->SetLineColor(kRed);
      hProPiRatioS4VertUnwgt->SetLineColor(kRed);
      hProS4HorzUnwgt->SetLineColor(kRed);
      hProS4VertUnwgt->SetLineColor(kRed);
      hPiS4HorzUnwgt->SetLineColor(kRed);
      hPiS4VertUnwgt->SetLineColor(kRed);
      hEffTotal->SetLineColor(kRed);
      hMSq->SetLineColor(kRed);
      hMom->SetLineColor(kRed);
      hMomS1S4->SetLineColor(kRed);
      // hCosmicsVertEff->SetLineColor(kRed);
      // hCosmicsVertEff2dNorm->SetLineColor(kRed);
      // hEffRatio->SetLineColor(kRed);
      // hEffRatio2dNorm->SetLineColor(kRed);
      leg->AddEntry(hdtof1d, "1 block", "l");
      legProS4Horz->AddEntry(hProS4Horz, Form("1 block - %.3g #pm %.1g per spill", intProS4Horz, eProS4Horz), "l");     
      legPiS4Horz->AddEntry(hPiS4Horz, Form("1 block - %d #pm %d per spill", (int)intPiS4Horz, (int)ePiS4Horz), "l");     
      legRatioVert->AddEntry(hProPiRatioS4Vert, "1 block", "l");
    }
    else if (nBlocks == 2) {
      // hdtof1d_sub->SetLineColor(kBlue);
      hdtof1d->SetLineColor(kBlue);
      hProPiRatioS4Horz->SetLineColor(kBlue);
      hProPiRatioS4Vert->SetLineColor(kBlue);
      hProS4Horz->SetLineColor(kBlue);
      hProS4HorzSmear->SetLineColor(kBlue);
      hProS4Vert->SetLineColor(kBlue);
      hPiS4Horz->SetLineColor(kBlue);
      hPiS4Vert->SetLineColor(kBlue);
      hProPiRatioS4HorzUnwgt->SetLineColor(kBlue);
      hProPiRatioS4VertUnwgt->SetLineColor(kBlue);
      hProS4HorzUnwgt->SetLineColor(kBlue);
      hProS4VertUnwgt->SetLineColor(kBlue);
      hPiS4HorzUnwgt->SetLineColor(kBlue);
      hPiS4VertUnwgt->SetLineColor(kBlue);
      hEffTotal->SetLineColor(kBlue);
      hMSq->SetLineColor(kBlue);
      hMom->SetLineColor(kBlue);
      hMomS1S4->SetLineColor(kBlue);
      // hCosmicsVertEff->SetLineColor(kBlue);
      // hEffRatio->SetLineColor(kBlue);
      // hEffRatio2dNorm->SetLineColor(kBlue);
      leg->AddEntry(hdtof1d, "2 blocks", "l");
      legProS4Horz->AddEntry(hProS4Horz, Form("2 blocks - %.3g #pm %.1g per spill", intProS4Horz, eProS4Horz), "l");     
      legPiS4Horz->AddEntry(hPiS4Horz, Form("2 blocks - %d #pm %d per spill", (int)intPiS4Horz, (int)ePiS4Horz), "l");
      legRatioVert->AddEntry(hProPiRatioS4Vert, "2 blocks", "l");
    }
    else if (nBlocks == 3) {
      // hdtof1d_sub->SetLineColor(kCyan+1);
      hdtof1d->SetLineColor(kCyan+1);
      hProPiRatioS4Horz->SetLineColor(kCyan+1);
      hProPiRatioS4Vert->SetLineColor(kCyan+1);
      hProS4Horz->SetLineColor(kCyan+1);
      hProS4HorzSmear->SetLineColor(kCyan+1);
      hProS4Vert->SetLineColor(kCyan+1);
      hPiS4Horz->SetLineColor(kCyan+1);
      hPiS4Vert->SetLineColor(kCyan+1);
      hProPiRatioS4HorzUnwgt->SetLineColor(kCyan+1);
      hProPiRatioS4VertUnwgt->SetLineColor(kCyan+1);
      hProS4HorzUnwgt->SetLineColor(kCyan+1);
      hProS4VertUnwgt->SetLineColor(kCyan+1);
      hPiS4HorzUnwgt->SetLineColor(kCyan+1);
      hPiS4VertUnwgt->SetLineColor(kCyan+1);
      hEffTotal->SetLineColor(kCyan+1);
      hMSq->SetLineColor(kCyan+1);
      hMom->SetLineColor(kCyan+1);
      hMomS1S4->SetLineColor(kCyan+1);
      // hCosmicsVertEff->SetLineColor(kCyan+1);
      // hCosmicsVertEff2dNorm->SetLineColor(kCyan+1);
      // hEffRatio->SetLineColor(kCyan+1);
      // hEffRatio2dNorm->SetLineColor(kCyan+1);
      leg->AddEntry(hdtof1d, "3 blocks", "l");
      legProS4Horz->AddEntry(hProS4Horz, Form("3 blocks - %.3g #pm %.g1 per spill", intProS4Horz, eProS4Horz), "l");     
      legPiS4Horz->AddEntry(hPiS4Horz, Form("3 blocks - %d #pm %d per spill", (int)intPiS4Horz, (int)ePiS4Horz), "l");
      legRatioVert->AddEntry(hProPiRatioS4Vert, "3 blocks", "l");
    }
    else if (nBlocks == 4) {
      hdtof1d->SetLineColor(kOrange+1);
      hProPiRatioS4Horz->SetLineColor(kOrange+1);
      hProPiRatioS4Vert->SetLineColor(kOrange+1);
      hProS4Horz->SetLineColor(kOrange+1);
      hProS4HorzSmear->SetLineColor(kOrange+1);
      hProS4Vert->SetLineColor(kOrange+1);
      hPiS4Horz->SetLineColor(kOrange+1);
      hPiS4Vert->SetLineColor(kOrange+1);
      hProPiRatioS4HorzUnwgt->SetLineColor(kOrange+1);
      hProPiRatioS4VertUnwgt->SetLineColor(kOrange+1);
      hProS4HorzUnwgt->SetLineColor(kOrange+1);
      hProS4VertUnwgt->SetLineColor(kOrange+1);
      hPiS4HorzUnwgt->SetLineColor(kOrange+1);
      hPiS4VertUnwgt->SetLineColor(kOrange+1);
      hEffTotal->SetLineColor(kOrange+1);
      hMSq->SetLineColor(kOrange+1);
      hMom->SetLineColor(kOrange+1);
      hMomS1S4->SetLineColor(kOrange+1);
      // hCosmicsVertEff->SetLineColor(kOrange+1);
      // hCosmicsVertEff2dNorm->SetLineColor(kOrange+1);
      // hEffRatio->SetLineColor(kOrange+1);
      // hEffRatio2dNorm->SetLineColor(kOrange+1);
      leg->AddEntry(hdtof1d, "4 blocks", "l");
      legProS4Horz->AddEntry(hProS4Horz, Form("4 blocks - %.3g #pm %.g1 per spill", intProS4Horz, eProS4Horz), "l");     
      legPiS4Horz->AddEntry(hPiS4Horz, Form("4 blocks - %d #pm %d per spill", (int)intPiS4Horz, (int)ePiS4Horz), "l");
      legRatioVert->AddEntry(hProPiRatioS4Vert, "4 blocks", "l");
    }

    h2dAngProS1->Scale(1. / nSpillsTrue);
    h2dAngPiS1->Scale(1. / nSpillsTrue);
    /*
      h2dAngProS1Unwgt->Scale(1. / nSpillsTrue);
      h2dAngPiS1Unwgt->Scale(1. / nSpillsTrue);
    */
    hMSq->Scale(1. / nSpillsTrue);
    hMom->Scale(1. / nSpillsTrue);
    hMomS1S4->Scale(1. / nSpillsTrue);
    h2dProMCComp->Scale(1. / nSpillsTrue);
    h2dProMCCompCut->Scale(1. / nSpillsTrue);
    h2DoubleFirst->Scale(1. / nSpillsTrue);
    h2DoubleSecond->Scale(1. / nSpillsTrue);

    // Print out ratio with error
    double rerr = ratioErr(intProS4Horz, eProS4Horz, dataS3.at(nBlocks), dataS3Err.at(nBlocks));
    cout<<"Ratio = "<<(intProS4Horz/dataS3.at(nBlocks))<<" +- "<<rerr<<endl;

    // hsEffRatio->Add(hEffRatio2dNorm);
    hsDtof->Add(hdtof1d);
    // hsBkgSub->Add(hdtof1d_sub);
    hsProS4Vert->Add(hProS4Vert);
    hsProS4Horz->Add(hProS4Horz);
    hsProS4HorzSmear->Add(hProS4HorzSmear);
    hsPiS4Vert->Add(hPiS4Vert);
    hsPiS4Horz->Add(hPiS4Horz);
    hsRatioS4Vert->Add(hProPiRatioS4Vert);
    hsRatioS4Horz->Add(hProPiRatioS4Horz);
    hsProS4VertUnwgt->Add(hProS4VertUnwgt);
    hsProS4HorzUnwgt->Add(hProS4HorzUnwgt);
    hsPiS4VertUnwgt->Add(hPiS4VertUnwgt);
    hsPiS4HorzUnwgt->Add(hPiS4HorzUnwgt);
    hsRatioS4VertUnwgt->Add(hProPiRatioS4VertUnwgt);
    hsRatioS4HorzUnwgt->Add(hProPiRatioS4HorzUnwgt);
    hsMom->Add(hMom);
    hsMomS1S4->Add(hMomS1S4);
    hsEffComp->Add(hEffTotal);
    // hsEffComp->Add(hCosmicsVertEff);
    hsEffComp2dNorm->Add(hEffTotal);
    // hsEffComp2dNorm->Add(hCosmicsVertEff2dNorm);
    hsEff->Add(hEffTotal);
    hProS4Vert->Write();
    hPiS4Vert->Write();
    hProS4Horz->Write();
    hProS4HorzSmear->Write();
    hPiS4Horz->Write();
    hProPiRatioS4Vert->Write();
    hProPiRatioS4Horz->Write();
    hProS4VertUnwgt->Write();
    hPiS4VertUnwgt->Write();
    hProS4HorzUnwgt->Write();
    hPiS4HorzUnwgt->Write();
    hProPiRatioS4VertUnwgt->Write();
    hProPiRatioS4HorzUnwgt->Write();
    hMSq->Write();
    hMom->Write();
    hMomS1S4->Write();
    h2dProMCComp->Write();
    h2dProMCCompCut->Write();
    // hCosmicsVertEff->Write();
    // hCosmicsVertEff2dNorm->Write();
    // hEffRatio2dNorm->Write();
    // hEffRatio->Write();
    hEffTotal->Write();
    hsEffComp->Write();
    hsEffComp2dNorm->Write();
    protonTree->Write();
    h2DoubleFirst->Write();
    h2DoubleSecond->Write();
    hSmearWeight->Write();

    h2dAngRatioS1->Divide(h2dAngProS1, h2dAngPiS1, 1., 1.);
    h2dAngProS1->Write();
    h2dAngPiS1->Write();
    h2dAngRatioS1->Write();
    h2dAngRatioS1Unwgt->Divide(h2dAngProS1Unwgt, h2dAngPiS1Unwgt, 1., 1.);
    h2dAngProS1Unwgt->Write();
    h2dAngPiS1Unwgt->Write();
    h2dAngRatioS1Unwgt->Write();

    legRatioVert->Write("legRatioVert");

    /*
    const char* nustof;
    if (nBlocks == 0) nustof = Form("%sData_2018_8_31_b2_800MeV_0block.root", ustofDir);
    else if (nBlocks==1) nustof = Form("%sData_2018_9_1_b4_800MeV_1block_bend4cm.root", ustofDir);
    else if (nBlocks==2) nustof = Form("%sData_2018_9_1_b2_800MeV_2block_bend4cm.root", ustofDir);
    else if (nBlocks==3) nustof = Form("%sData_2018_9_1_b3_800MeV_3block_bend4cm.root", ustofDir);
 
    // Read in ustof file
    TFile *finustof = new TFile(nustof, "read");
    TTree *utree = (TTree*)finustof->Get("tree");
    double tTrig;
    double tS1;
    double tSoSd;
    float xToF[50];
    float yToF[50];
    int nhit;
    utree->SetBranchAddress("tTrig", &tTrig);
    //    utree->SetBranchAddress("tS1", &tS1);
    utree->SetBranchAddress("xToF", xToF);
    utree->SetBranchAddress("yToF", yToF);
    utree->SetBranchAddress("nhit", &nhit);
    //    utree->SetBranchAddress("tSoSd", &tSoSd);

    hAllS4Horz->Add(hProS4Horz);
    hAllS4Horz->Add(hPiS4Horz);

    TH1D *hXAngleS1S2 = new TH1D(Form("hXAngleS1S2%d", nBlocks), Form("Angular distribution of hits in S3 (S1 & S2 triggers), %d blocks; #theta / degrees; Events / spill", nBlocks), 100, -3.8, 6.2);
    for (int t=0; t<utree->GetEntries(); t++) {
      utree->GetEntry(t);
      // Has an S1 and S2 hit
      if (tTrig != 0 ) {
	for (int n=0; n<nhit; n++) {
	  double positionX = ((xToF[n] - 4.) / 152.)*(s3EndX - s3StartX) + s3StartX;
	  double positionY = ((xToF[n] - 4.) / 152.)*(s3s1EndY - s3s1StartY)+s3s1StartY;
	  double angleOffAxis = TMath::ATan(positionX / positionY) * (180. / TMath::Pi());
	  hXAngleS1S2->Fill(angleOffAxis);
	} // for (int n=0; n<nhit; n++) 
      } // if (tTrig !=0 ) 
    } // for (int t=0; t<utree->GetEntries(); t++)

      // Integrate this hist between the angular limits of S2 then scale by this number
    const double s2ThetaLow = -0.359;
    const double s2ThetaHi  = 3.957;
    const double s4ThetaLow = 0.401;
    const double s4ThetaHi  = 6.083;
    double s2Int = hXAngleS1S2->Integral(hXAngleS1S2->GetXaxis()->FindBin(s2ThetaLow), hXAngleS1S2->GetXaxis()->FindBin(s2ThetaHi));
      
    hXAngleS1S2->Scale(1. / s2Int);    
    double s4Int = hXAngleS1S2->Integral(hXAngleS1S2->GetXaxis()->FindBin(s4ThetaLow), hXAngleS1S2->GetXaxis()->FindBin(s4ThetaHi));
    cout<<"S4 correction factor "<<s4Int<<endl;
    */
  } // for (int nBlocks = 0; nBlocks < 4; nBlocks++) 

  fout->cd();
  hsEff->Write();
  hsDtof->Write();
  // hsBkgSub->Write();
  hsProS4Vert->Write();
  hsPiS4Vert->Write();
  hsProS4VertUnwgt->Write();
  hsPiS4VertUnwgt->Write();
  legProS4Horz->Write("legProS4Horz");
  hsProS4Horz->Write();
  hsProS4HorzSmear->Write();
  hsProS4HorzUnwgt->Write();
  legPiS4Horz->Write("legPiS4Horz");
  hsPiS4Horz->Write();
  hsPiS4HorzUnwgt->Write();
  hsMom->Write();
  leg->Write("legOrdinary");
  hsRatioS4Vert->Write();
  hsRatioS4Horz->Write();
  hsRatioS4VertUnwgt->Write();
  hsRatioS4HorzUnwgt->Write();

  TArrow *arPi  = new TArrow(0.0195, 8e2, 0.0195, 1e2, 0.03, "|>");
  TArrow *arPro = new TArrow(0.8804, 8e2, 0.8804, 3e1, 0.03, "|>");
  TText *tPi  = new TText(.15, 2.5e2, "pion");
  TText *tPro = new TText(1., 2.5e2, "proton"); 
  arPi->Write("piArrow");
  tPi->Write("piText");
  arPro->Write("protonArrow");
  tPro->Write("protonText");

  fout->Close();
  delete fout;  
} // angularDistS4
