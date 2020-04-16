// angularDistS3.C
// Does the analysis for the S3 flux data
#include "UsefulFunctions.C"
// Outputs momentum in GeV/c

double dtVar(const double slope, const double slopeErr, const double hits, 
	     const double constant, const double constantErr)
{
  double weight = 1. / (slope * hits + constant);
  double var = pow(weight, 2) + pow(weight, 4) * (pow(slope,2)*hits + pow(hits*slopeErr,2) + pow(constantErr,2));
  return var;
}

void angularDistS3_newSample(const char* outfile,
		   const char* ustofDir="/nfs/scratch0/dbrailsf/data_backup/utof_backup_firsthitpinnedtounixtime/Data_root_v3_wo_walk_corr/",
		   const char* dstofDir="/nfs/scratch0/dbrailsf/data_backup/dtof_backup/",
		   const char* spillDir="/scratch0/sjones/spillDB/") 
{
  gROOT->SetBatch(kTRUE);
  // Edges of S3 in beam coordinate system
  const double s3StartX = -0.601158 + 0.026;
  const double s3EndX   = 1.08086 + 0.026;
  const double s3s1StartY = 9.07239 + 1.767;
  const double s3s1EndY   = 8.89731 + 1.767;
  // Measured from the wire chamber
  const double s3StartWcX = -0.601158 + 0.052;
  const double s3EndWcX   = 1.08086 + 0.052;
  const double s3StartWcY = 9.07239 + 2.054;
  const double s3EndWcY   = 8.89731 + 2.054;
  // z coordinates
  const double s3BarTop    = 62.; 
  const double s3BarBottom = -60.;

  // Deadtime corrections
  // Just use a constant ratio for the 0 block case
  const double block0Slope    = 0.;//-0.0003738;
  const double block0SlopeErr = 0;//0.00006863;
  const double block0Const    = 0.0913247;//0.2332; 
  const double block0ConstErr = 0.00102858;//0.02436;
  const double block1Slope    = -0.0002569;
  const double block1SlopeErr = 0.00001487;
  const double block1Const    =  0.3965;
  const double block1ConstErr = 0.0146;
  const double block2Slope    = -0.0001739;
  const double block2SlopeErr = 6.511e-6;
  const double block2Const    =  0.419;
  const double block2ConstErr = 0.009931;
  const double block3Slope    = -0.0001498;
  const double block3SlopeErr = 1.038e-5;
  const double block3Const    = 0.4343;
  const double block3ConstErr = 0.01836;
  // 4 block data
  const double block4Slope1    = -0.00030855237;
  const double block4Slope1Err = 9.1525e-06;
  const double block4Const1    = 0.639835;
  const double block4Const1Err = 0.012901;
  const double block4Slope2    = -0.00036724024;
  const double block4Slope2Err = 5.9539e-6;
  const double block4Const2    = 0.81216136;
  const double block4Const2Err = 0.0087677;
  const double block4Slope3    = -0.00031197486;
  const double block4Slope3Err = 1.020427e-5;
  const double block4Const3    = 0.71958558;
  const double block4Const3Err = 0.0151212;

  std::vector<double> block4SlopeVec = {block4Slope1, block4Slope2, block4Slope3};
  std::vector<double> block4SlopeErrVec = {block4Slope1Err, block4Slope2Err, block4Slope3Err};
  std::vector<double> block4ConstVec = {block4Const1, block4Const2, block4Const3};
  std::vector<double> block4ConstErrVec = {block4Const1Err, block4Const2Err, block4Const3Err};
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
  // Time of flight cuts for S1 to S3
  // Particles travelling at c should cross the distance in between about 36.9ns
  // and 35.6ns
  const double tLight = 35.8; // Expected time in ns of light speed particles according to AK
  // This is before the shift is applied
  const double proLow  = 53;
  const double proHiOther = 115.;
  const double proHi0     = 80.;
  const double piLow = 35.75;
  const double piHi  = 37.75;
  // S3 amplitude cut for protons
  // Apply to A1ToF and A2ToF
  // AK's standard cut
  const double ACut = 0.25;
  // SJ's bar by bar amplitude cut (by eye)
  // For A1 
  const std::vector<double> A1CutVec = {0.25, 0.25, 0.275, 0.2, 0.25, 0.25, 0.225, 0.25, 0.3, 
					0.3, 0.3,
					0.2, 0.2, 0.25, 0.25, 0.25, 0.25, 0.3, 0.25, 
					0.225, 0.3, 0.3};
  // For A2
  const std::vector<double> A2CutVec = {0.3, 0.275, 0.25, 0.125, 0.25, 
					0.25, 0.15, 0.225, 0.25, 0.25, 0.225,
					0.25, 0.25, 0.2, 0.2, 0.225, 
					0.225, 0.25, 0.25, 0.225, 0.3, 0.25};
  // Time cut for hits interacting in neighbouring S3 bars
  const double twoBarCut = 0.5; // ns

  // Fit boxes for protons
  const vector<double> proFitLow1 = {54., 56., 60., 65., 80.};
  const vector<double> proFitHi1  = {57., 60., 66., 76., 110.};
  const vector<double> proFitLow2 = {57., 60., 65., 76., 80.};
  const vector<double> proFitHi2  = {63., 70., 80., 103., 110.};

  THStack *hsThetaS1pro   = new THStack("hsThetaS1pro", "S1 #cap S3 angular distribution of proton hits; #theta / degrees; Events / spill / degree");
  THStack *hsThetaS1proNoS2 = new THStack("hsThetaS1proNoS2", "S1 #cap S3 angular distribution of proton hits; #theta / degrees; Events / spill / degree");
  THStack *hsThetaS1S2pro = new THStack("hsThetaS1S2pro", "S1 #cap S2 #cap S3 angular distribution of proton hits; #theta / degrees; Events / spill / degree");
  THStack *hsPhiS1pro     = new THStack("hsPhiS1pro", "S1 #cap S3 angular distribution of proton hits; #phi / degrees; Events / spill / degree");
  THStack *hsPhiS1S2pro   = new THStack("hsPhiS1S2pro", "S1 #cap S2 #cap S3 angular distribution of proton hits; #phi / degrees; Events / spill / degree");

  THStack *hsThetaS1pi   = new THStack("hsThetaS1pi", "S1 #cap S3 angular distribution of MIP hits; #theta / degrees; Events / spill / degree");
  THStack *hsThetaS1piNoS2 = new THStack("hsThetaS1piNoS2", "S1 #cap S3 angular distribution of MIP hits; #theta / degrees; Events / spill / degree");
  THStack *hsThetaS1S2pi = new THStack("hsThetaS1S2pi", "S1 #cap S2 #cap S3 angular distribution of MIP hits; #theta / degrees; Events / spill / degree");
  THStack *hsPhiS1pi     = new THStack("hsPhiS1pi", "S1 #cap S3 angular distribution of MIP hits; #phi / degrees; Events / spill / degree");
  THStack *hsPhiS1S2pi   = new THStack("hsPhiS1S2pi", "S1 #cap S2 #cap S3 angular distribution of MIP hits; #phi / degrees; Events / spill / degree");

  THStack *hsPhiS1S2ratio   = new THStack("hsPhiS1S2ratio", "S1 #cap S2 #cap S3 angular distribution of proton/MIP ratio; #phi / degrees;  Protons/MIPs");
  THStack *hsThetaS1S2ratio = new THStack("hsThetaS1S2ratio", "S1 #cap S2 #cap S3 angular distribution of proton/MIP ratio; #theta / degrees;  Protons/MIPs");
  THStack *hsPhiS1ratio   = new THStack("hsPhiS1ratio", "S1 #cap S3 angular distribution of proton/MIP ratio; #phi / degrees;  Protons/MIPs");
  THStack *hsThetaS1ratio = new THStack("hsThetaS1ratio", "S1 #cap S3 angular distribution of proton/MIP ratio; #theta / degrees;  Protons/MIPs");
  THStack *hsThetaS1ratioNoS2 = new THStack("hsThetaS1ratioNoS2", "S1 #cap S3 angular distribution of proton/MIP ratio; #theta / degrees;  Protons/MIPs");

  THStack *hsMomS1S2 = new THStack("hsMomS1S2", "Proton momentum measured in S3 (S1 & S2 trigger); Proton momentum [GeV/c]; Events / spill");
  THStack *hsMomS1 = new THStack("hsMomS1", "Proton momentum measured in S3; Proton momentum [GeV/c]; Events / spill");
  THStack *hsMomTpc = new THStack("hsMomTpc", "Proton momentum measured in S3 for particles passing through TPC; Proton momentum [GeV/c]; Events / spill");
  THStack *hsKE = new THStack("hsKE", "Proton kinetic energy measured for protons crossing the TPC; Proton kinetic energy / MeV; Events / spill");
  THStack *hsutof1dS1S2 = new THStack("hsutof1dS1S2", "Time of flight as measured in S3 (S1 & S2 trigger); Time of flight / ns; Events / spill");
  THStack *hsutof1dS1   = new THStack("hsutof1dS1", "Time of flight as measured in S3; Time of flight / ns; Events / spill");
  THStack *hsutof1dS1NoS2   = new THStack("hsutof1dS1NoS2", "Time of flight as measured in S3 (no S2 trigger); Time of flight / ns; Events / spill");

  TFile *fout = new TFile(outfile, "recreate");

  TLegend *leg = new TLegend(0.15, 0.55, 0.3, 0.85);
  TLegend *legTheta = new TLegend(0.23, 0.5, 0.38, 0.85);
  TLegend *legTof = new TLegend(0.71, 0.53, 0.88, 0.85);
  TLegend *legPhiRatio = new TLegend(0.38, 0.56, 0.55, 0.88);

  // Legends with numbers of particles
  TLegend *legThetaS1pro = new TLegend(0.6, 0.6, 0.9, 0.9);
  TLegend *legThetaS1pi  = new TLegend(0.6, 0.6, 0.9, 0.9);
  TLegend *legKE  = new TLegend(0.6, 0.6, 0.9, 0.9);
  TLegend *legTpc = new TLegend(.6, .6, .9, .9);

  TH2D *hMom2D_0blkQ = new TH2D("hMom2D_0blkQ", "Quick peak 0 block; x / cm; Bar", 100, -10, 170, 22, 0.5, 22.5);
  TH2D *hMom2D_0blkS = new TH2D("hMom2D_0blkS", "Slow peak 0 block; x / cm; Bar", 100, -10, 170, 22, 0.5, 22.5);
  TH2D *hMom2D_1blkQ = new TH2D("hMom2D_1blkQ", "Quick peak 1 block; x / cm; Bar", 100, -10, 170, 22, 0.5, 22.5);
  TH2D *hMom2D_1blkS = new TH2D("hMom2D_1blkS", "Slow peak 1 block; x / cm; Bar", 100, -10, 170, 22, 0.5, 22.5);
  TH2D *hMom2D_2blkQ = new TH2D("hMom2D_2blkQ", "Quick peak 2 block; x / cm; Bar", 100, -10, 170, 22, 0.5, 22.5);
  TH2D *hMom2D_2blkS = new TH2D("hMom2D_2blkS", "Slow peak 2 block; x / cm; Bar", 100, -10, 170, 22, 0.5, 22.5);

  double proHi = 0.;
  for (int nBlocks = 0; nBlocks < 5; nBlocks++) {
    cout<<"=================================="<<endl;
    cout<<nBlocks<<" blocks"<<endl;
    cout<<"=================================="<<endl;

    // Tree for the S3 proton weights
    TTree *protonTree = new TTree(Form("protonTree%dBlocks", nBlocks), Form("protonTree%dBlocks", nBlocks));
    double tof, mom;
    double x, y, z;
    double weight, error;
    int spill;
    int isS2;
    protonTree->Branch("tof", &tof);
    protonTree->Branch("mom", &mom);
    protonTree->Branch("x", &x);
    protonTree->Branch("y", &y);
    protonTree->Branch("z", &z);
    protonTree->Branch("weight", &weight);
    protonTree->Branch("error", &error);
    protonTree->Branch("spill", &spill);
    protonTree->Branch("isS2", &isS2);

    int nSpills = 0;

    if (nBlocks == 0) proHi = proHi0;
    else proHi = proHiOther;

    TH1D *hKE = new TH1D(Form("hKE_%d", nBlocks), Form("Proton kinetic energy measured for protons crossing the TPC, %d blocks; Proton kinetic energy / MeV; Events", nBlocks), 100, 40., 350.);
    setHistAttr(hKE);
    vector<double> keErr;
    keErr.resize(hKE->GetNbinsX()+2, 0);

    TH1D *hThetaS1pro   = new TH1D(Form("hThetaS1pro%d", nBlocks), Form("Angular distribution of proton hits in S3 (S1 trigger only), %d blocks; #theta / degrees; Events / spill / degree", nBlocks), binnum, binsTheta);
    setHistAttr(hThetaS1pro);
    vector<double> thetaS1proErr;
    thetaS1proErr.resize(hThetaS1pro->GetNbinsX()+2, 0);

    TH1D *hThetaS1proNoS2 = new TH1D(Form("hThetaS1proNoS2_%d", nBlocks), Form("Angular distribution of proton hits in S3 (S1 trigger only), %d blocks; #theta / degrees; Events / spill / degree", nBlocks), binnum, binsTheta);
    hThetaS1proNoS2->Sumw2();
    TH1D *hThetaS1S2pro = new TH1D(Form("hThetaS1S2pro%d", nBlocks), Form("Angular distribution of proton hits in S3 (S1 & S2 triggers), %d blocks; #theta / degrees; Events / spill", nBlocks), binnum, binsTheta);
    hThetaS1S2pro->Sumw2();
    TH1D *hPhiS1pro   = new TH1D(Form("hPhiS1pro%d", nBlocks), Form("Angular distribution of proton hits in S3 (S1 trigger only), %d blocks; #phi / degrees; Events / spill / degree", nBlocks), 22, -3.22, 3.35);
    setHistAttr(hPhiS1pro);
    vector<double> phiS1proErr;
    phiS1proErr.resize(hPhiS1pro->GetNbinsX()+2, 0);

    TH1D *hPhiS1S2pro = new TH1D(Form("hPhiS1S2pro%d", nBlocks), Form("Angular distribution of proton hits in S3 (S1 & S2 triggers), %d blocks; #phi / degrees; Events / spill", nBlocks), 22, -3.22, 3.35);
    hPhiS1S2pro->Sumw2();

    TH1D *hThetaS1pi   = new TH1D(Form("hThetaS1pi%d", nBlocks), Form("Angular distribution of pion hits in S3 (S1 trigger only), %d blocks; #theta / degrees; Events / spill / degree", nBlocks), binnum, binsTheta);
    setHistAttr(hThetaS1pi);
    vector<double> thetaS1piErr;
    thetaS1piErr.resize(hThetaS1pi->GetNbinsX()+2, 0);

    TH1D *hThetaS1piNoS2 = new TH1D(Form("hThetaS1piNoS2_%d", nBlocks), Form("Angular distribution of pion hits in S3 (S1 trigger only), %d blocks; #theta / degrees; Events / spill / degree", nBlocks), binnum, binsTheta);
    hThetaS1piNoS2->Sumw2();
    TH1D *hThetaS1S2pi = new TH1D(Form("hThetaS1S2pi%d", nBlocks), Form("Angular distribution of pion hits in S3 (S1 & S2 triggers), %d blocks; #theta / degrees; Events / spill", nBlocks), binnum, binsTheta);
    hThetaS1S2pi->Sumw2();
    TH1D *hPhiS1pi   = new TH1D(Form("hPhiS1pi%d", nBlocks), Form("Angular distribution of pion hits in S3 (S1 trigger only), %d blocks; #phi / degrees; Events / spill / degree", nBlocks), 22, -3.22, 3.35);
    setHistAttr(hPhiS1pi);
    vector<double> phiS1piErr;
    phiS1piErr.resize(hPhiS1pi->GetNbinsX()+2, 0);

    TH1D *hPhiS1S2pi = new TH1D(Form("hPhiS1S2pi%d", nBlocks), Form("Angular distribution of pion hits in S3 (S1 & S2 triggers), %d blocks; #phi / degrees; Events / spill", nBlocks), 22, -3.22, 3.35);
    hPhiS1S2pi->Sumw2();

    TH1D *hThetaS1ratio = new TH1D(Form("hThetaS1ratio%d", nBlocks), Form("S1 #cap S3 angular distribution of proton/MIP ratio, %d blocks; #phi / degrees; Protons/MIPs", nBlocks), binnum, binsTheta);
    setHistAttr(hThetaS1ratio);
    TH1D *hThetaS1ratioNoS2 = new TH1D(Form("hThetaS1ratioNoS2_%d", nBlocks), Form("S1 #cap S3 angular distribution of proton/MIP ratio, %d blocks; #phi / degrees; Protons/MIPs", nBlocks), binnum, binsTheta);
    TH1D *hPhiS1ratio   = new TH1D(Form("hPhiS1ratio%d", nBlocks), Form("S1 #cap S3 angular distribution of proton/MIP, %d blocks; #phi / degrees; Protons/MIPs", nBlocks), 22, -3.22, 3.35);
    setHistAttr(hPhiS1ratio);
    TH1D *hThetaS1S2ratio = new TH1D(Form("hThetaS1S2ratio%d", nBlocks), Form("S1 #cap S2 #cap S3 angular distribution of proton/MIP ratio, %d blocks; #phi / degrees; Protons/MIPs", nBlocks), binnum, binsTheta);
    setHistAttr(hThetaS1S2ratio);
    TH1D *hPhiS1S2ratio   = new TH1D(Form("hPhiS1S2ratio%d", nBlocks), Form("S1 #cap S2 #cap S3 angular distribution of proton/MIP ratio, %d blocks; #phi / degrees; Protons/MIPs", nBlocks), 22, -3.22, 3.35);
    setHistAttr(hPhiS1S2ratio);
    TH1D *hutof1dS1 = new TH1D(Form("hutof1dS1_%d",nBlocks), Form("Time of flight, %d blocks (S1 trigger only); S3 - S1 / ns; Events / spill", nBlocks), 250, 25, 125);
    setHistAttr(hutof1dS1);
    vector<double> utof1dS1Err;
    utof1dS1Err.resize(hutof1dS1->GetNbinsX()+2, 0);

    TH1D *hutof1dS1NoS2 = new TH1D(Form("hutof1dS1NoS2_%d",nBlocks), Form("Time of flight, %d blocks (S1 trigger only); S3 - S1 / ns; Events / spill", nBlocks), 250, 25, 125);
    hutof1dS1NoS2->Sumw2();
    TH1D *hutof1dS1S2 = new TH1D(Form("hutof1dS1S2_%d",nBlocks), Form("Time of flight, %d blocks (S1 & S2 trigger); S3 - S1 / ns; Events / spill", nBlocks), 250, 25, 125);
    hutof1dS1S2->Sumw2();

    TH1D *hMomS1S2 = new TH1D(Form("hMomS1S2_%d",nBlocks), Form("Proton momentum measured in S3, %d blocks; Proton momentum [GeV/c]; Events / spill", nBlocks), 120, 0.3, 0.9);
    setHistAttr(hMomS1S2);
    TH1D *hMomS1 = new TH1D(Form("hMomS1_%d",nBlocks), Form("Proton momentum measured in S3, %d blocks; Proton momentum [GeV/c]; Events / spill", nBlocks), 120, 0.3, 0.9);
    setHistAttr(hMomS1);
    vector<double> momS1Err;
    momS1Err.resize(hMomS1->GetNbinsX()+2, 0);

    TH1D *hMomTpc = new TH1D(Form("hMomTpc_%d",nBlocks), Form("Proton momentum measured in S3 passing through TPC, %d blocks; Proton momentum [GeV/c]; Events / spill", nBlocks), 120, 0.3, 0.9);
    setHistAttr(hMomTpc);
    vector<double> momTpcErr;
    momTpcErr.resize(hMomTpc->GetNbinsX()+2, 0);

    // XY distributions in S3
    TH2D *hprotonXY = new TH2D(Form("hprotonXY%d", nBlocks), Form("S3 spatial distribution of proton hits, %d blocks; x / cm; y / cm; Events / spill", nBlocks), 105, 0., 168., 22, 0., 120.);
    TH2D *hpionXY   = new TH2D(Form("hpionXY%d", nBlocks), Form("S3 spatial distribution of MIP hits, %d blocks; x / cm; y / cm; Events / spill", nBlocks), 100, 0., 168., 22, 0., 120.);
    TH2D *hAllXY    = new TH2D(Form("hAllXY%d", nBlocks), Form("S3 spatial distribution of all hits, %d blocks; #theta / degrees; #phi / degrees; Events / spill", nBlocks), 100, -3.8, 6.2, 22, -3.22, 3.35);

    // 2D angular distributions
    // S1 origin
    TH2D *h2dAngProS1 = new TH2D(Form("h2dAngProS1_%d", nBlocks), Form("S3 angular distribution of protons hits, %d blocks; #theta / degrees; #phi / degrees; Events / spill", nBlocks), 105, -3.05, 5.95, 22, -3.15, 3.34);
    TH2D *h2dAngPiS1 = new TH2D(Form("h2dAngPiS1_%d", nBlocks), Form("S3 angular distribution of MIP hits, %d blocks; #theta / degrees; #phi / degrees; Events / spill", nBlocks), 105, -3.05, 5.95, 22, -3.15, 3.34); 
    TH2D *h2dAngRatioS1 = new TH2D(Form("h2dAngRatioS1_%d", nBlocks), Form("S3 angular distribution of proton/MIP ratio, %d blocks; #theta / degrees; #phi / degrees; Protons/MIPs", nBlocks), 105, -3.05, 5.95, 22, -3.15, 3.34);
    // Beam position monitor as the origin
    TH2D *h2dAngProWC = new TH2D(Form("h2dAngProWC_%d", nBlocks), Form("S3 angular distribution of protons hits, %d blocks; #theta / degrees; #phi / degrees; Events / spill", nBlocks), 105, -2.75, 5.95, 22, -3.07, 3.29);
    TH2D *h2dAngPiWC = new TH2D(Form("h2dAngPiWC_%d", nBlocks), Form("S3 angular distribution of MIP hits, %d blocks; #theta / degrees; #phi / degrees; Events / spill", nBlocks), 105, -2.75, 5.95, 22, -3.07, 3.29);
    TH2D *h2dAngRatioWC = new TH2D(Form("h2dAngRatioWC_%d", nBlocks), Form("S3 angular distribution of proton/MIP ratio, %d blocks; #theta / degrees; #phi / degrees; Protons/MIPs", nBlocks), 105, -2.75, 5.95, 22, -3.07, 3.29);
    // 2D time of flight vs. angle
    TH2D *h2dTofPhiS1 = new TH2D(Form("h2dTofPhiS1_%d", nBlocks), Form("S3 time of flight vs. off-axis angle, %d blocks; t_{S3} - t_{S1} / ns; #phi / degrees; Events / spill", nBlocks), 250, 25, 125, 22, -3.15, 3.34);
    TH2D *h2dTofThetaS1 = new TH2D(Form("h2dTofThetaS1_%d", nBlocks), Form("S3 time of flight vs. off-axis angle, %d blocks; t_{S3} - t_{S1} / ns; #theta / degrees; Events / spill", nBlocks), 250, 25, 125, 105, -3.05, 5.95);
    TH2D *h2dTofPhiWC = new TH2D(Form("h2dTofPhiWC_%d", nBlocks), Form("S3 time of flight vs. off-axis angle, %d blocks; t_{S3} - t_{S1} / ns; #phi / degrees; Events / spill", nBlocks), 250, 25, 125, 22, -3.15, 3.34);
    TH2D *h2dTofThetaWC = new TH2D(Form("h2dTofThetaWC_%d", nBlocks), Form("S3 time of flight vs. off-axis angle, %d blocks; t_{S3} - t_{S1} / ns; #theta / degrees; Events / spill", nBlocks), 250, 25, 125, 105, -3.05, 5.95);
    setHistAttr(h2dTofPhiWC);
    setHistAttr(h2dTofThetaWC);
    setHistAttr(h2dTofPhiS1);
    setHistAttr(h2dTofThetaS1);
    setHistAttr(h2dAngProS1);
    setHistAttr(h2dAngPiS1);
    setHistAttr(h2dAngRatioS1);
    setHistAttr(h2dAngProWC);
    setHistAttr(h2dAngPiWC);
    setHistAttr(h2dAngRatioWC);

    // Define signal and background functions to be fitted
    TF1 *sPro1 = new TF1(Form("sPro1_%d", nBlocks), "gaus", proFitLow1[nBlocks],proFitHi1[nBlocks]);
    TF1 *sPro2 = new TF1(Form("sPro2_%d", nBlocks), "gaus", proFitLow2[nBlocks],proFitHi2[nBlocks]);
    TF1 *sPi  = new TF1(Form("sPi_%d", nBlocks), "gaus", piLow-1, piHi+1);
    TF1 *fBkgExp  = new TF1(Form("fBkgExp_%d", nBlocks), "expo", piLow, proHi);
    TF1 *fBkgFlat = new TF1(Form("fBkgFlat_%d", nBlocks), "pol0", piLow, proHi);
    sPro1->SetLineColor(kGreen+2);
    sPi->SetLineColor(kRed);
    TF1 *fSplusBExp = new TF1(Form("signal_plus_bkg_exp_%d", nBlocks), "gaus(0)+gaus(3)+gaus(6)+expo(9)", piLow, proHi);
    TF1 *fSplusBFlat = new TF1(Form("signal_plus_bkg_flat_%d", nBlocks), "gaus(0)+gaus(3)+gaus(6)+pol0(9)", piLow-1., proHi);
    fSplusBExp->SetParNames("const 1", "mean 1", "sigma 1",
			    "const 2", "mean 2", "sigma 2",
			    "const 3", "mean 3", "sigma 3",
			    "bkgconst", "bkgdecay");
    fSplusBFlat->SetParNames("const 1", "mean 1", "sigma 1",
			    "const 2", "mean 2", "sigma 2",
			    "const 3", "mean 3", "sigma 3",
			    "bkgflat");
    fSplusBExp->SetLineColor(kBlack);
      // Number of protons and number of MIPs
      int nP  = 0;
      int nPi = 0;

      // Find the correct dstof files
      Int_t runMin=-1;
      Int_t runMax=-1;

      vector<double> startTimes;
      vector<double> endTimes;

      vector<const char*> nustof;
      vector<double> slope;
      vector<double> slopeErr;
      vector<double> constant;
      vector<double> constantErr;
      if (nBlocks==0) {
	nustof.push_back(str0Block);
	startTimes.push_back(start0Block);
	endTimes.push_back(end0Block);
	slope.push_back(block0Slope);
	constant.push_back(block0Const);
	slopeErr.push_back(block0SlopeErr);
	constantErr.push_back(block0ConstErr);
      }
      else if (nBlocks==1) {
	nustof.push_back(str1Block);
	startTimes.push_back(start1Block);
	endTimes.push_back(end1Block);
	slope.push_back(block1Slope);
	constant.push_back(block1Const);
	slopeErr.push_back(block1SlopeErr);
	constantErr.push_back(block1ConstErr);
      }
      else if (nBlocks==2) {
	nustof.push_back(str2Block);
	startTimes.push_back(start2Block);
	endTimes.push_back(end2Block);
	slope.push_back(block2Slope);
	constant.push_back(block2Const);
	slopeErr.push_back(block2SlopeErr);
	constantErr.push_back(block2ConstErr);
      }
      else if (nBlocks==3) {
	nustof.push_back(str3Block);
	startTimes.push_back(start3Block);
	endTimes.push_back(end3Block);
	slope.push_back(block3Slope);
	constant.push_back(block3Const);
	slopeErr.push_back(block3SlopeErr);
	constantErr.push_back(block3ConstErr);
      }
      else if (nBlocks==4) {
	nustof = str4BlockVec;
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
	slope = block4SlopeVec;
	constant = block4ConstVec;
	slopeErr = block4SlopeErrVec;
	constantErr = block4ConstErrVec;
      }

      for (int sub=0; sub<startTimes.size(); sub++) {

	// Find dtof runs
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
      
	  if (firstTemp>endTimes.at(sub)){
	    break;
	  }
      
	  if (firstTemp<startTimes.at(sub) && lastTemp>startTimes.at(sub)){
	    runMin = irun;
	  }
      
	  if (firstTemp<endTimes.at(sub) && lastTemp>endTimes.at(sub)){
	    runMax = irun;
	  }   
	} // for (int irun=950; irun<1400; irun++) 
    
	cout << "Min and max dtof runs are " << runMin << " " << runMax << endl;

	std::vector<double> dtofTimes;
	std::vector<int> dtofS1S2Hits;
	std::vector<double> utofTimes;
	std::vector<int> utofS1S2Hits;
 
	// Open the appropriate spill DB files and get the spill times
	for (int irun = runMin; irun < runMax+1; irun++) {
	  TFile *dbFile = new TFile(Form("%s/spillDB_run%d_run%d.root", spillDir, irun, irun), "read");
	  TTree *spillTree = (TTree*)dbFile->Get("spillTree");
	  double globalSpillTime;
	  double ustofSpillTime;
	  spillTree->SetBranchAddress("globalSpillTime", &globalSpillTime);
	  spillTree->SetBranchAddress("ustofSpillTime", &ustofSpillTime);
	  for (int t = 0; t < spillTree->GetEntries(); t++) {
	    spillTree->GetEntry(t);
	    if (globalSpillTime >= startTimes.at(sub) && globalSpillTime <= endTimes.at(sub)) {
	      dtofTimes.push_back(globalSpillTime);
	      utofTimes.push_back(ustofSpillTime);
	    }
	  } // for (int t = 0; t < spillTree->GetEntries(); t++)

	  dbFile->Close();
	  delete dbFile;
	} // for (int irun = runMin; irun < runMax+1; irun++) 

	dtofS1S2Hits.resize(dtofTimes.size(), 0);
	utofS1S2Hits.resize(dtofTimes.size(), 0);

	cout<<"Finding number of dtof hits in each spill"<<endl;
	int lastt = 0;
	int lastrun = 0;
	for (int s=0; s<dtofTimes.size(); s++) {
	  if (s % 100 == 0) {
	    cout<<"Spill "<<s<<" of "<<dtofTimes.size()<<endl;
	  }
	  // Loop over the all the files
	  for (int irun = runMin; irun < runMax+1; irun++) {
	
	    TFile *dtofFile = new TFile(Form("%srun%d/DsTOFtreeRun%d_tdc1.root", dstofDir, irun, irun), "read");
	    RawDsTofHeader *tof = NULL;
	    TTree *tofTree = (TTree*)dtofFile->Get("tofTree");
	    tofTree->SetBranchAddress("tof", &tof);
	    tofTree->GetEntry(0);
	    double firstTime = tof->unixTime;
	    tofTree->GetEntry(tofTree->GetEntries()-1);
	    double lastTime = tof->unixTime;
	    double lastS1S2Dtof = 0.;
	    // Spill is in this file
	    if (firstTime <= dtofTimes[s] && lastTime >= dtofTimes[s]) {
	      if (irun != lastrun) {
		lastt = 0;
		lastrun = irun;
	      }
	      // Loop over all entries and count the number of S1 S2 hits within the spill
	      for (int t=lastt; t<tofTree->GetEntries(); t++) {
		tofTree->GetEntry(t);
		if ((tof->fakeTimeNs/1e9)+firstTime < dtofTimes[s]) continue;
		if ((tof->fakeTimeNs/1e9)+firstTime > dtofTimes[s]+1.) break;
		// Is within a spill
		if ((tof->fakeTimeNs/1e9)+firstTime >= dtofTimes[s] &&
		    (tof->fakeTimeNs/1e9)+firstTime <= dtofTimes[s]+1. && 
		    tof->channel == 13 && (tof->fakeTimeNs - lastS1S2Dtof) > 500.) {
		  dtofS1S2Hits[s]++;
		  lastt = t;
		  lastS1S2Dtof = tof->fakeTimeNs;
		} // Is within the spill
	      } // for (int t=0; t<tofTree->GetEntries(); t++)
	    } // if (firstTime <= dtofTimes[s] && lastTime >=  dtofTimes[s])

	    delete tof;
	    dtofFile->Close();
	    delete dtofFile;
	  } // for (int irun = runMin; irun < runMax+1; irun++) 
	} // for (int s=0; s<dtofTimes.size(); s++)

	TFile *futof = new TFile(Form("%s/%s", ustofDir, nustof.at(sub)), "read");

	double tToF[50];
	float xToF[50];
	float yToF[50];
	float A1ToF[50];
	float A2ToF[50];
	double tTrig;
	double tS1;
	double tSoSd;
	int nhit;
	int nBar[50];

	TTree *tree = (TTree*)futof->Get("tree");

	tree->SetBranchAddress("xToF", xToF);
	tree->SetBranchAddress("yToF", yToF);
	tree->SetBranchAddress("A1ToF", A1ToF);
	tree->SetBranchAddress("A2ToF", A2ToF);
	tree->SetBranchAddress("nhit", &nhit);
	tree->SetBranchAddress("tS1", &tS1);
	tree->SetBranchAddress("tToF", tToF);
	tree->SetBranchAddress("tTrig", &tTrig);
	tree->SetBranchAddress("tSoSd", &tSoSd);
	tree->SetBranchAddress("nBar", nBar);

	double lastSpill = 0.; 

	TNamed *start = 0;
	TNamed *end   = 0;
	futof->GetObject("start_of_run", start);
	futof->GetObject("end_of_run", end);
    
	const char* startchar = start->GetTitle();
	std::string startstr(startchar);
	std::string unixstart = startstr.substr(25,10);
	int startTimeUtof = stoi(unixstart);

	int lastut = 0;
	// Loop over the spills and perform the adjustment for each spill
	for (int s = 0; s < utofTimes.size(); s++) {
	  if (s % 100 == 0) cout<<"Getting hits from spill "<<s<<" of "<<utofTimes.size()<<endl;
	  double deadtimeWeight = dtofS1S2Hits[s] * slope.at(sub) + constant.at(sub);
	  double deadtimeVar = dtVar(slope.at(sub), slopeErr.at(sub), dtofS1S2Hits[s], constant.at(sub), constantErr.at(sub));
	  // Initial data quality loop
	  bool isGood = false;
	  for (int t=lastut; t<tree->GetEntries(); t++) {
	    tree->GetEntry(t);
	    if ((tS1/1e9) + startTimeUtof < utofTimes[s]) continue;
	    if ((tS1/1e9) + startTimeUtof > utofTimes[s] + 1.) break;

	    if ((tTrig/1e9) + startTimeUtof > utofTimes[s] + 0.47) {
	      isGood = true;
	      nSpills++;
	      break;
	    }
	  } // Initial data quality loop
	  // Only do this spill if it's good
	  if (isGood) {
	    for (int t=lastut; t<tree->GetEntries(); t++) {
	      tree->GetEntry(t);
	      if ((tS1/1e9) + startTimeUtof < utofTimes[s]) continue;
	      if ((tS1/1e9) + startTimeUtof > utofTimes[s] + 1.) break;
	    
	      for (int nh=0; nh<nhit; nh++) {
		// Checks if the next S3 hit is in a neighbouring bar
		// Only want to count one of them if it is
		// Check from this hit onwards to see if there are any double hit candidates
		bool isDouble = false;
		if (nh < nhit-1) {
		  for (int nh2=nh+1; nh2<nhit; nh2++) {
		    // Timing cut and checks if in the neighbouring bar
		    if (abs(tToF[nh] - tToF[nh2]) < twoBarCut && 
			(nBar[nh] == nBar[nh2]-1 || nBar[nh] == nBar[nh2]+1)) {
		      isDouble = true;
		      break;
		    }
		  } // for (int nh2=nh; nh2<nhit; nh2++)
		} // if (nh < nh-1)
		// If it's not a double hit then go ahead
		if (!isDouble) {
		  double tofCalc = tToF[nh] - tS1;
		  // Calculate x, y z positions relative to S1
		  TVector3 utofCoords = GetUtofGlobalCoords(xToF[nh], yToF[nh]);
		  double positionX = (xToF[nh]/168)*(s3EndX - s3StartX) + s3StartX;
		  double positionY = (xToF[nh]/168.)*(s3s1EndY - s3s1StartY) + s3s1StartY;
		  double positionWcX = (xToF[nh]/168)*(s3EndWcX - s3StartWcX) + s3StartWcX;
		  double positionWcY = (xToF[nh]/168.)*(s3EndWcY - s3StartWcY) + s3StartWcY; 
		  double positionZ = (yToF[nh] + s3BarBottom + 2.75) / 100.;
		  double angleTheta = getThetaFromGlobal(utofCoords);
		  double anglePhi   = getPhiFromGlobal(utofCoords);
		  double angleWcTheta = TMath::ATan(positionWcX / positionWcY) * (180./TMath::Pi());
		  double angleWcPhi   = TMath::ATan(positionZ / positionWcY) * (180./TMath::Pi());
		  double travelledDist = TMath::Sqrt(pow(positionX,2)+pow(positionY,2)+pow(positionZ,2));
		  // All triggers
		  hAllXY->Fill(angleTheta, anglePhi, 1./deadtimeWeight);
		  hutof1dS1->Fill(tofCalc, 1./deadtimeWeight);
		  utof1dS1Err.at(hutof1dS1->GetXaxis()->FindBin(tofCalc)) += deadtimeVar;
		  h2dTofThetaS1->Fill(tofCalc, angleTheta, 1./deadtimeWeight);
		  h2dTofPhiS1->Fill(tofCalc, anglePhi, 1./deadtimeWeight);
		  h2dTofThetaWC->Fill(tofCalc, angleWcTheta, 1./deadtimeWeight);
		  h2dTofPhiWC->Fill(tofCalc, angleWcPhi, 1./deadtimeWeight);
		  if (tTrig == 0) hutof1dS1NoS2->Fill(tofCalc, 1./deadtimeWeight);
		  // Separate protons and MIPs using timing and amplitude cuts
		  // Is a MIP 
		  if ( tofCalc > piLow && tofCalc < piHi ) {
		    nPi++;
		    hThetaS1pi->Fill(angleTheta, 1./deadtimeWeight);
		    hPhiS1pi->Fill(anglePhi, 1./deadtimeWeight);
		    phiS1piErr.at(hPhiS1pi->GetXaxis()->FindBin(anglePhi)) += deadtimeVar;
		    thetaS1piErr.at(hThetaS1pi->GetXaxis()->FindBin(angleTheta)) += deadtimeVar;
		    hpionXY->Fill(xToF[nh], yToF[nh], 1./deadtimeWeight);
		    h2dAngPiS1->Fill(angleTheta, anglePhi, 1./deadtimeWeight);
		    h2dAngPiWC->Fill(angleWcTheta, angleWcPhi, 1./deadtimeWeight);
		    if (tTrig==0) hThetaS1piNoS2->Fill(angleTheta, 1./deadtimeWeight);
		    lastut = t;
		  } // if ( tofCalc > (tLight - (piLow+piHi)/2.) + piLow && tofCalc < (tLight - (piLow+piHi)/2.) + piHi )
		  // Is a proton
		  else if ( tofCalc > proLow && tofCalc < proHi && A1ToF[nh] > A1CutVec[nBar[nh]] && A2ToF[nh] > A2CutVec[nBar[nh]]) {
		    // Variables for the proton tree
		    tof = tofCalc;
		    mom = momFromTime(0.938, 10.9, tofCalc);
		    x = utofCoords.X() + 0.491;
		    y = utofCoords.Y() + 0.0114;
		    z = utofCoords.Z() - 10.829;
		    isS2 = 0;
		    weight = 1. / deadtimeWeight;
		    error = TMath::Sqrt(deadtimeVar);
		    if (tTrig != 0) isS2 = 1;
		    spill = nSpills;
		    protonTree->Fill();

		    nP++;
		    hThetaS1pro->Fill(angleTheta, 1./deadtimeWeight);
		    hPhiS1pro->Fill(anglePhi, 1./deadtimeWeight);
		    phiS1proErr.at(hPhiS1pro->GetXaxis()->FindBin(anglePhi)) += deadtimeVar;
		    thetaS1proErr.at(hThetaS1pro->GetXaxis()->FindBin(angleTheta)) += deadtimeVar;
		    if (tTrig==0) hThetaS1proNoS2->Fill(angleTheta, 1./deadtimeWeight);
		    hprotonXY->Fill(xToF[nh], yToF[nh], 1./deadtimeWeight);
		    h2dAngProS1->Fill(angleTheta, anglePhi, 1./deadtimeWeight);
		    h2dAngProWC->Fill(angleWcTheta, angleWcPhi, 1./deadtimeWeight);
		    // Remove deuteron peak in 0 block data
		    if (nBlocks != 0) {
		      hMomS1->Fill(momFromTime(0.938, 10.9, tofCalc)/1000., 1./deadtimeWeight);
		      momS1Err.at(hMomS1->GetXaxis()->FindBin(momFromTime(0.938, 10.9, tofCalc)/1000.)) += deadtimeVar;
		      // Only protons passing through TPC active area
		      if (angleTheta < tpcThetaHigh && angleTheta > tpcThetaLow &&
			  anglePhi > tpcPhiLow && anglePhi < tpcPhiHigh) {
			hMomTpc->Fill(momFromTime(0.938, 10.9, tofCalc)/1000., 1./deadtimeWeight);
			hKE->Fill(keFromTime(0.938, 10.9, tofCalc)/1000., 1./deadtimeWeight);
			momTpcErr.at(hMomTpc->GetXaxis()->FindBin(momFromTime(0.938, 10.9, tofCalc)/1000.)) += deadtimeVar;
			keErr.at(hKE->GetXaxis()->FindBin(keFromTime(0.938, 10.8, tofCalc)/1000.)) += deadtimeVar;
		      }
		    }
		    else {
		      double mom = momFromTime(0.938, 10.9, tofCalc)/1000.;
		      if (mom > 0.45) {
			hMomS1->Fill(mom, 1./deadtimeWeight);
			momS1Err.at(hMomS1->GetXaxis()->FindBin(momFromTime(0.938, 10.9, tofCalc)/1000.)) += deadtimeVar;
			// Only protons passing through TPC active area
			if (angleTheta > tpcThetaLow && angleTheta < tpcThetaHigh &&
			    anglePhi > tpcPhiLow && anglePhi < tpcPhiHigh) {
			  hMomTpc->Fill(momFromTime(0.938, 10.8, tofCalc)/1000., 1./deadtimeWeight);
			  hKE->Fill(keFromTime(0.938, 10.8, tofCalc)*1000., 1./deadtimeWeight);
			  momTpcErr.at(hMomTpc->GetXaxis()->FindBin(momFromTime(0.938, 10.8, tofCalc)/1000.)) += deadtimeVar;
			  keErr.at(hKE->GetXaxis()->FindBin(keFromTime(0.938, 10.8, tofCalc)*1000.)) += deadtimeVar;
			}
		      }
		    }
		    lastut = t;
		    if (nBlocks == 0 && momFromTime(0.938, 10.9, tofCalc)/1000. > 0.595) hMom2D_0blkQ->Fill(xToF[nh], nBar[nh]);
		    else if (nBlocks == 0 && momFromTime(0.938, 10.9, tofCalc)/1000. < 0.595) hMom2D_0blkS->Fill(xToF[nh], nBar[nh]);
		    else if (nBlocks == 1 && momFromTime(0.938, 10.9, tofCalc)/1000. > 0.570) hMom2D_1blkQ->Fill(xToF[nh], nBar[nh]);
		    else if (nBlocks == 1 && momFromTime(0.938, 10.9, tofCalc)/1000. < 0.570) hMom2D_1blkS->Fill(xToF[nh], nBar[nh]);
		    else if (nBlocks == 2 && momFromTime(0.938, 10.9, tofCalc)/1000. > 0.525) hMom2D_2blkQ->Fill(xToF[0], nBar[nh]);
		    else if (nBlocks == 2 && momFromTime(0.938, 10.9, tofCalc)/1000. < 0.525) hMom2D_2blkS->Fill(xToF[nh], nBar[nh]);
	    
		  } // else if ( tofCalc > (tLight - (piLow+piHi)/2.) + proLow && tofCalc < (tLight - (piLow+piHi)/2.) + proHi )
		  //	}
		  // S1 & S2 trigger only
		  if (tTrig !=0) {
		    hutof1dS1S2->Fill(tofCalc, 1./deadtimeWeight);
		    // Separate protons and MIPs using timing and amplitude cuts
		    // Is a MIP
		    if ( tofCalc > piLow && tofCalc < piHi ) {
		      nPi++;
		      hThetaS1S2pi->Fill(angleTheta, 1./deadtimeWeight);
		      hPhiS1S2pi->Fill(anglePhi, 1./deadtimeWeight);
		    } // if ( tofCalc > (tLight - (piLow+piHi)/2.) + piLow && tofCalc < (tLight - (piLow+piHi)/2.) + piHi )
		    // Is a proton
		    else if ( tofCalc > proLow && tofCalc < proHi && A1ToF[nh] > A1CutVec[nBar[nh]] && A2ToF[nh] > A2CutVec[nBar[nh]]) {
		      nP++;
		      hThetaS1S2pro->Fill(angleTheta, 1./deadtimeWeight);
		      hPhiS1S2pro->Fill(anglePhi, 1./deadtimeWeight);
		      hMomS1S2->Fill(momFromTime(0.938, 10.9, tofCalc)/1000., 1./deadtimeWeight);
		      // hMomZS12->Fill(momFromTime(0.938, 10.9, tofCalc), nBar[nh]);
		      // hMomYS12->Fill(momFromTime(0.938, 10.9, tofCalc), xToF[nh]);
		    } // else if ( tofCalc > (tLight - (piLow+piHi)/2.) + proLow && tofCalc < (tLight - (piLow+piHi)/2.) + proHi )
		  } // S1 + S2 trigger
		}
	      } // Loop over nhits
	    } // for (int t=0; t<tree->GetEntries(); t++)
	  } // if (isGood)
	}

	cout<<"Utof spills "<<nSpills<<" vs. "<<utofTimes.size()<<" in spill DB"<<endl;
	fout->cd();
      } // Loop over subsamples

      // Sort the errors
      for (int bin = 0; bin < hThetaS1pi->GetNbinsX()+1; bin++) {
	hThetaS1pi->SetBinError(bin, TMath::Sqrt(thetaS1piErr.at(bin)));
	hThetaS1pro->SetBinError(bin, TMath::Sqrt(thetaS1proErr.at(bin)));
      }
      for (int bin = 0; bin < hPhiS1pi->GetNbinsX()+1; bin++) {
	hPhiS1pi->SetBinError(bin, TMath::Sqrt(phiS1piErr.at(bin)));
	hPhiS1pro->SetBinError(bin, TMath::Sqrt(phiS1proErr.at(bin)));
      }
      for (int bin = 0; bin < hutof1dS1->GetNbinsX()+1; bin++) {
	hutof1dS1->SetBinError(bin, TMath::Sqrt(utof1dS1Err.at(bin)));
      }
      for (int bin = 0; bin < hMomS1->GetNbinsX()+1; bin++) {
	hMomS1->SetBinError(bin, TMath::Sqrt(momS1Err.at(bin)));
      }
      for (int bin = 0; bin < hMomTpc->GetNbinsX()+1; bin++) {
	hMomTpc->SetBinError(bin, TMath::Sqrt(momTpcErr.at(bin)));
      }
      for (int bin = 0; bin < hKE->GetNbinsX()+1; bin++) {
	hKE->SetBinError(bin, TMath::Sqrt(keErr.at(bin)));
      }

      hThetaS1pi->Sumw2();
      hThetaS1pro->Sumw2();
      hPhiS1pi->Sumw2();
      hPhiS1pro->Sumw2();
      hutof1dS1->Sumw2();
      hMomS1->Sumw2();
      hMomTpc->Sumw2();
      hKE->Sumw2();

      hThetaS1S2ratio->Divide(hThetaS1S2pro, hThetaS1S2pi, 1., 1., "B");
      hPhiS1S2ratio->Divide(hPhiS1S2pro, hPhiS1S2pi, 1., 1., "B");
      hThetaS1ratio->Divide(hThetaS1pro, hThetaS1pi, 1., 1., "B");
      hPhiS1ratio->Divide(hPhiS1pro, hPhiS1pi, 1., 1., "B");
      hThetaS1ratioNoS2->Divide(hThetaS1proNoS2, hThetaS1piNoS2, 1., 1., "B");

      hPhiS1S2pro->Scale(1. / (double)nSpills);
      hPhiS1S2pi->Scale(1. / (double)nSpills);
      hThetaS1S2pro->Scale(1. / (double)nSpills);
      hThetaS1S2pi->Scale(1. / (double)nSpills);
      hAllXY->Scale(1. / (double)nSpills);
      hPhiS1pro->Scale(1. / (double)nSpills);
      hPhiS1pi->Scale(1. / (double)nSpills);
      hPhiS1pro->Scale(1., "width");
      hPhiS1pi->Scale(1., "width");
      hPhiS1S2pro->Scale(1., "width");
      hPhiS1S2pi->Scale(1., "width");
      hThetaS1pro->Scale(1. / (double)nSpills);
      hThetaS1pi->Scale(1. / (double)nSpills);
      hThetaS1proNoS2->Scale(1. / (double)nSpills);
      hThetaS1piNoS2->Scale(1. / (double)nSpills);

      hThetaS1pi->Scale(1, "width");
      hThetaS1pro->Scale(1, "width");
      hThetaS1piNoS2->Scale(1, "width");
      hThetaS1proNoS2->Scale(1, "width");
      hThetaS1S2pi->Scale(1, "width");
      hThetaS1S2pro->Scale(1, "width");
      /*
	for (int i=1; i < hThetaS1pi->GetNbinsX(); i++) {
	double binWidth = (hThetaS1pi->GetXaxis()->GetBinUpEdge(i) - hThetaS1pi->GetXaxis()->GetBinLowEdge(i));
	hThetaS1pi->SetBinContent(i, hThetaS1pi->GetBinContent(i) / binWidth);
	hThetaS1pro->SetBinContent(i, hThetaS1pro->GetBinContent(i) / binWidth);
	hThetaS1piNoS2->SetBinContent(i, hThetaS1piNoS2->GetBinContent(i) / binWidth);
	hThetaS1proNoS2->SetBinContent(i, hThetaS1proNoS2->GetBinContent(i) / binWidth);
	hThetaS1S2pi->SetBinContent(i, hThetaS1S2pi->GetBinContent(i) / binWidth);
	hThetaS1S2pro->SetBinContent(i, hThetaS1S2pro->GetBinContent(i) / binWidth);
	}
      */

      hMomS1S2->Scale(1. / (double)nSpills);
      hMomS1->Scale(1. / (double)nSpills);
      hMomTpc->Scale(1. / (double)nSpills);
      hKE->Scale(1. / (double)nSpills);
      hutof1dS1S2->Scale(1. / (double)nSpills);
      hutof1dS1->Scale(1. / (double)nSpills);
      hutof1dS1NoS2->Scale(1. / (double)nSpills);

      double ePro = 0.;
      double ePi  = 0.;
      double eKE  = 0.;
      double eTpc = 0.;
      double iPro = hThetaS1pro->IntegralAndError(1, hThetaS1pro->GetNbinsX(), ePro, "width");
      double iPi  = hThetaS1pi->IntegralAndError(1, hThetaS1pi->GetNbinsX(), ePi, "width");
      double iKE  = hKE->IntegralAndError(1, hKE->GetNbinsX(), eKE);
      double iTpc = hMomTpc->IntegralAndError(1, hMomTpc->GetNbinsX(), eTpc);
      cout<<iPi<<" +- "<<ePi<<endl;
      cout<<iPro<<" +- "<<ePro<<endl;
      cout<<iKE<<" +- "<<eKE<<endl;

      if (nBlocks==0) {
	hThetaS1S2pro->SetLineColor(kBlack);
	hThetaS1S2pi->SetLineColor(kBlack);
	hPhiS1S2pro->SetLineColor(kBlack);
	hPhiS1S2pi->SetLineColor(kBlack);
	hPhiS1S2ratio->SetLineColor(kBlack);
	hThetaS1S2ratio->SetLineColor(kBlack);
	hThetaS1pro->SetLineColor(kBlack);
	hThetaS1pi->SetLineColor(kBlack);
	hThetaS1proNoS2->SetLineColor(kBlack);
	hThetaS1piNoS2->SetLineColor(kBlack);
	hPhiS1pro->SetLineColor(kBlack);
	hPhiS1pi->SetLineColor(kBlack);
	hPhiS1ratio->SetLineColor(kBlack);
	hThetaS1ratio->SetLineColor(kBlack);
	hThetaS1ratioNoS2->SetLineColor(kBlack);
	hMomS1S2->SetLineColor(kBlack);
	hMomS1->SetLineColor(kBlack);
	hMomTpc->SetLineColor(kBlack);
	hKE->SetLineColor(kBlack);
	hutof1dS1S2->SetLineColor(kBlack);
	hutof1dS1->SetLineColor(kBlack);
	hutof1dS1NoS2->SetLineColor(kBlack);
      
	legThetaS1pro->AddEntry(hThetaS1pro, 
				Form("0 blocks - %d #pm %d per spill", (int)iPro, (int)ePro), "le"); 
	legThetaS1pi->AddEntry(hThetaS1pi, 
			       Form("0 blocks - %d #pm %d per spill", (int)iPi, (int)ePi), "le"); 
	legKE->AddEntry(hKE, Form("0 blocks - %.3g #pm %.2g per spill", iKE, eKE), "le"); 
	legTpc->AddEntry(hMomTpc, Form("0 blocks - (%.3g #pm %.2g) / spill", iTpc, eTpc), "le");
	leg->AddEntry(hThetaS1S2pro, "0 blocks", "le");
	legTheta->AddEntry(hThetaS1ratio, "0 blocks", "le");
	legTof->AddEntry(hutof1dS1, "0 blocks", "le");
	legPhiRatio->AddEntry(hPhiS1ratio, "0 blocks", "le");

	hMom2D_0blkQ->Scale(1. / (double)nSpills);
	hMom2D_0blkS->Scale(1. / (double)nSpills);
	hMom2D_0blkQ->Write();
	hMom2D_0blkS->Write();
      }
      if (nBlocks==1) {
	hThetaS1S2pro->SetLineColor(kRed);
	hThetaS1S2pi->SetLineColor(kRed);
	hPhiS1S2pro->SetLineColor(kRed);
	hPhiS1S2pi->SetLineColor(kRed);
	hPhiS1S2ratio->SetLineColor(kRed);
	hThetaS1S2ratio->SetLineColor(kRed);
	hThetaS1pro->SetLineColor(kRed);
	hThetaS1pi->SetLineColor(kRed);
	hThetaS1proNoS2->SetLineColor(kRed);
	hThetaS1piNoS2->SetLineColor(kRed);
	hPhiS1pro->SetLineColor(kRed);
	hPhiS1pi->SetLineColor(kRed);
	hPhiS1ratio->SetLineColor(kRed);
	hThetaS1ratio->SetLineColor(kRed);
	hThetaS1ratioNoS2->SetLineColor(kRed);
	hMomS1S2->SetLineColor(kRed);
	hMomS1->SetLineColor(kRed);
	hKE->SetLineColor(kRed);
	hMomTpc->SetLineColor(kRed);
	hutof1dS1S2->SetLineColor(kRed);
	hutof1dS1->SetLineColor(kRed);
	hutof1dS1NoS2->SetLineColor(kRed);
      
	legThetaS1pro->AddEntry(hThetaS1pro, 
				Form("1 block - %d #pm %d per spill", (int)iPro, (int)ePro), "le"); 
	legThetaS1pi->AddEntry(hThetaS1pi, 
			       Form("1 block - %d #pm %d per spill", (int)iPi, (int)ePi), "le"); 
	legKE->AddEntry(hKE, Form("1 block - %.3g #pm %.2g per spill", iKE, eKE), "le"); 
	legTpc->AddEntry(hMomTpc, Form("1 block - (%.3g #pm %.2g) / spill", iTpc, eTpc), "le");
	leg->AddEntry(hThetaS1S2pro, "1 block", "le");
	legTheta->AddEntry(hThetaS1ratio, "1 block", "le");
	legTof->AddEntry(hutof1dS1, "1 block", "le");
	legPhiRatio->AddEntry(hPhiS1ratio, "1 block", "le");
	hMom2D_1blkQ->Scale(1. / (double)nSpills);
	hMom2D_1blkS->Scale(1. / (double)nSpills);
	hMom2D_1blkQ->Write();
	hMom2D_1blkS->Write();
      }
      if (nBlocks==2) {
	hThetaS1S2pro->SetLineColor(kBlue);
	hThetaS1S2pi->SetLineColor(kBlue);
	hPhiS1S2pro->SetLineColor(kBlue);
	hPhiS1S2pi->SetLineColor(kBlue);
	hPhiS1S2ratio->SetLineColor(kBlue);
	hThetaS1S2ratio->SetLineColor(kBlue);
	hThetaS1pro->SetLineColor(kBlue);
	hThetaS1pi->SetLineColor(kBlue);
	hThetaS1piNoS2->SetLineColor(kBlue);
	hThetaS1proNoS2->SetLineColor(kBlue);
	hPhiS1pro->SetLineColor(kBlue);
	hPhiS1pi->SetLineColor(kBlue);
	hPhiS1ratio->SetLineColor(kBlue);
	hThetaS1ratio->SetLineColor(kBlue);
	hThetaS1ratioNoS2->SetLineColor(kBlue);
	hMomS1S2->SetLineColor(kBlue);
	hMomS1->SetLineColor(kBlue);
	hMomTpc->SetLineColor(kBlue);
	hKE->SetLineColor(kBlue);
	hutof1dS1S2->SetLineColor(kBlue);
	hutof1dS1->SetLineColor(kBlue);
	hutof1dS1NoS2->SetLineColor(kBlue);
	legThetaS1pro->AddEntry(hThetaS1pro, 
				Form("2 blocks - %d #pm %d per spill", (int)iPro, (int)ePro), "le"); 
	legThetaS1pi->AddEntry(hThetaS1pi, 
			       Form("2 blocks - %d #pm %d per spill", (int)iPi, (int)ePi), "le"); 
	legKE->AddEntry(hKE, Form("2 blocks - %.3g #pm %.2g per spill", iKE, eKE), "le"); 
	legTpc->AddEntry(hMomTpc, Form("2 blocks - (%.3g #pm %.2g) / spill", iTpc, eTpc), "le");
	leg->AddEntry(hThetaS1S2pro, "2 blocks", "le");
	legTheta->AddEntry(hThetaS1ratio, "2 blocks", "le");
	legTof->AddEntry(hutof1dS1, "2 blocks", "le");
	legPhiRatio->AddEntry(hPhiS1ratio, "2 blocks", "le");

	hMom2D_2blkQ->Scale(1. / (double)nSpills);
	hMom2D_2blkS->Scale(1. / (double)nSpills);
	hMom2D_2blkQ->Write();
	hMom2D_2blkS->Write();
      }
      if (nBlocks==3) {
	hThetaS1S2pro->SetLineColor(kCyan+1);
	hThetaS1S2pi->SetLineColor(kCyan+1);
	hPhiS1S2pro->SetLineColor(kCyan+1);
	hPhiS1S2pi->SetLineColor(kCyan+1);
	hPhiS1S2ratio->SetLineColor(kCyan+1);
	hThetaS1S2ratio->SetLineColor(kCyan+1);
	hThetaS1pro->SetLineColor(kCyan+1);
	hThetaS1pi->SetLineColor(kCyan+1);
	hThetaS1proNoS2->SetLineColor(kCyan+1);
	hThetaS1piNoS2->SetLineColor(kCyan+1);
	hPhiS1pro->SetLineColor(kCyan+1);
	hPhiS1pi->SetLineColor(kCyan+1);
	hPhiS1ratio->SetLineColor(kCyan+1);
	hThetaS1ratio->SetLineColor(kCyan+1);
	hThetaS1ratioNoS2->SetLineColor(kCyan+1);
	hMomS1S2->SetLineColor(kCyan+1);
	hMomS1->SetLineColor(kCyan+1);
	hKE->SetLineColor(kCyan+1);
	hMomTpc->SetLineColor(kCyan+1);
	hutof1dS1S2->SetLineColor(kCyan+1);
	hutof1dS1->SetLineColor(kCyan+1);
	hutof1dS1NoS2->SetLineColor(kCyan+1);
	legThetaS1pro->AddEntry(hThetaS1pro, 
				Form("3 blocks - %d #pm %d per spill", (int)iPro, (int)ePro), "le"); 
	legThetaS1pi->AddEntry(hThetaS1pi, 
			       Form("3 blocks - %d #pm %d per spill", (int)iPi, (int)ePi), "le"); 
	legKE->AddEntry(hKE, Form("3 blocks - %.3g #pm %.2g per spill", iKE, eKE), "le");
	legTpc->AddEntry(hMomTpc, Form("3 blocks - (%.3g #pm %.2g) / spill", iTpc, eTpc), "le");
	leg->AddEntry(hThetaS1S2pro, "3 blocks", "le");
	legTheta->AddEntry(hThetaS1ratio, "3 blocks", "le");
	legTof->AddEntry(hutof1dS1, "3 blocks", "le");
	legPhiRatio->AddEntry(hPhiS1ratio, "3 blocks", "le");

	hMom2D_2blkQ->Scale(1. / (double)nSpills);
	hMom2D_2blkS->Scale(1. / (double)nSpills);
	hMom2D_2blkQ->Write();
	hMom2D_2blkS->Write();
      }
      if (nBlocks==4) {
	hThetaS1S2pro->SetLineColor(kOrange+1);
	hThetaS1S2pi->SetLineColor(kOrange+1);
	hPhiS1S2pro->SetLineColor(kOrange+1);
	hPhiS1S2pi->SetLineColor(kOrange+1);
	hPhiS1S2ratio->SetLineColor(kOrange+1);
	hThetaS1S2ratio->SetLineColor(kOrange+1);
	hThetaS1pro->SetLineColor(kOrange+1);
	hThetaS1pi->SetLineColor(kOrange+1);
	hThetaS1proNoS2->SetLineColor(kOrange+1);
	hThetaS1piNoS2->SetLineColor(kOrange+1);
	hPhiS1pro->SetLineColor(kOrange+1);
	hPhiS1pi->SetLineColor(kOrange+1);
	hPhiS1ratio->SetLineColor(kOrange+1);
	hThetaS1ratio->SetLineColor(kOrange+1);
	hThetaS1ratioNoS2->SetLineColor(kOrange+1);
	hMomS1S2->SetLineColor(kOrange+1);
	hMomS1->SetLineColor(kOrange+1);
	hMomTpc->SetLineColor(kOrange+1);
	hKE->SetLineColor(kOrange+1);
	hutof1dS1S2->SetLineColor(kOrange+1);
	hutof1dS1->SetLineColor(kOrange+1);
	legThetaS1pro->AddEntry(hThetaS1pro, 
				Form("4 blocks - %d #pm %d per spill", (int)iPro, (int)ePro), "le"); 
	legThetaS1pi->AddEntry(hThetaS1pi, 
			       Form("4 blocks - %d #pm %d per spill", (int)iPi, (int)ePi), "le"); 
	legKE->AddEntry(hKE, Form("4 blocks - %.2g #pm %.1g per spill", iKE, eKE), "le"); 
	legTpc->AddEntry(hMomTpc, Form("4 blocks - (%.3g #pm %.2g) / spill", iTpc, eTpc), "le");
	leg->AddEntry(hThetaS1S2pro, "4 blocks", "le");
	legTheta->AddEntry(hThetaS1ratio, "4 blocks", "le");
	legTof->AddEntry(hutof1dS1, "4 blocks", "le");
	legPhiRatio->AddEntry(hPhiS1ratio, "4 blocks", "le");

	hMom2D_2blkQ->Scale(1. / (double)nSpills);
	hMom2D_2blkS->Scale(1. / (double)nSpills);
	hMom2D_2blkQ->Write();
	hMom2D_2blkS->Write();
      }

      hThetaS1S2ratio->Write();
      hPhiS1S2ratio->Write();
      hThetaS1ratio->Write();
      hThetaS1ratioNoS2->Write();
      hPhiS1ratio->Write();

      hutof1dS1->Fit(sPi, "R");
      hutof1dS1->Fit(sPro1, "R");
      hutof1dS1->Fit(sPro2, "R");
      hutof1dS1->Fit(fBkgExp, "R");
      hutof1dS1->Fit(fBkgFlat, "R");
      Double_t parExp[11];
      sPi->GetParameters(&parExp[0]);
      sPro1->GetParameters(&parExp[3]);
      sPro2->GetParameters(&parExp[6]);
      fBkgExp->GetParameters(&parExp[9]);
      fSplusBExp->SetParameters(parExp);
      hutof1dS1->Fit(fSplusBExp, "R");
      hutof1dS1->Draw("hist");
      fSplusBExp->Draw("same");
      fSplusBExp->Write();    

      Double_t parFlat[10];
      sPi->GetParameters(&parFlat[0]);
      sPro1->GetParameters(&parFlat[3]);
      sPro2->GetParameters(&parFlat[6]);
      fBkgFlat->GetParameters(&parFlat[9]);
      fSplusBFlat->SetParameters(parFlat);
      hutof1dS1->Fit(fSplusBFlat, "R");
      fSplusBFlat->Write();

      hprotonXY->Scale(1. / (double)nSpills);
      hpionXY->Scale(1. / (double)nSpills);
      hprotonXY->Write();
      hpionXY->Write();
      hAllXY->Write();

      hThetaS1S2pro->Write();
      hThetaS1S2pi->Write();
      hPhiS1S2pro->Write();
      hPhiS1S2pi->Write();
      hThetaS1pro->Write();
      hThetaS1pi->Write();
      hThetaS1proNoS2->Write();
      hThetaS1piNoS2->Write();
      hPhiS1pro->Write();
      hPhiS1pi->Write();

      hMomS1S2->Write();
      hMomS1->Write();
      hMomTpc->Write();
      hKE->Write();

      hutof1dS1S2->Write();
      hutof1dS1->Write();
      hutof1dS1NoS2->Write();

      hsThetaS1S2pro->Add(hThetaS1S2pro);
      hsThetaS1S2pi->Add(hThetaS1S2pi);
      hsPhiS1S2pro->Add(hPhiS1S2pro);
      hsPhiS1S2pi->Add(hPhiS1S2pi);
      hsThetaS1S2ratio->Add(hThetaS1S2ratio);
      hsPhiS1S2ratio->Add(hPhiS1S2ratio);

      hsThetaS1pro->Add(hThetaS1pro);
      hsThetaS1pi->Add(hThetaS1pi);
      hsThetaS1proNoS2->Add(hThetaS1proNoS2);
      hsThetaS1piNoS2->Add(hThetaS1piNoS2);
      hsPhiS1pro->Add(hPhiS1pro);
      hsPhiS1pi->Add(hPhiS1pi);
      hsThetaS1ratio->Add(hThetaS1ratio);
      hsThetaS1ratioNoS2->Add(hThetaS1ratioNoS2);
      hsPhiS1ratio->Add(hPhiS1ratio);

      hsutof1dS1S2->Add(hutof1dS1S2);
      hsutof1dS1->Add(hutof1dS1);
      hsutof1dS1NoS2->Add(hutof1dS1NoS2);

      hsMomS1S2->Add(hMomS1S2);
      hsMomS1->Add(hMomS1);
      hsMomTpc->Add(hMomTpc);
      hsKE->Add(hKE);

      h2dAngPiS1->Scale(1. / (double)nSpills);
      h2dAngProS1->Scale(1. / (double)nSpills);
      h2dAngPiWC->Scale(1. / (double)nSpills);
      h2dAngProWC->Scale(1. / (double)nSpills);
      h2dAngRatioWC->Divide(h2dAngProWC, h2dAngPiWC, 1., 1.);
      h2dAngRatioS1->Divide(h2dAngProS1, h2dAngPiS1, 1., 1.);
      h2dAngPiS1->Write();
      h2dAngProS1->Write();
      h2dAngPiWC->Write();
      h2dAngProWC->Write();
      h2dAngRatioWC->Write();
      h2dAngRatioS1->Write();

      h2dTofThetaS1->Scale(1. / (double)nSpills); 
      h2dTofThetaWC->Scale(1. / (double)nSpills);
      h2dTofPhiS1->Scale(1. / (double)nSpills); 
      h2dTofPhiWC->Scale(1. / (double)nSpills);
      h2dTofThetaS1->Write();
      h2dTofThetaWC->Write();
      h2dTofPhiS1->Write();
      h2dTofPhiWC->Write();
            
      

      fout->cd();
      protonTree->Write();
  } // for (int nBlocks = 0; nBlocks <= 4; nBlocks++) 
  fout->cd();
  leg->Write("leg");

  hsThetaS1S2pro->Write();
  hsThetaS1S2pi->Write();
  hsPhiS1S2pro->Write();
  hsPhiS1S2pi->Write();
  hsThetaS1S2ratio->Write();
  hsPhiS1S2ratio->Write();

  hsutof1dS1S2->Write();
  hsutof1dS1->Write();
  hsutof1dS1NoS2->Write();

  hsThetaS1pro->Write();
  hsThetaS1pi->Write();
  hsThetaS1proNoS2->Write();
  hsThetaS1piNoS2->Write();
  hsPhiS1pro->Write();
  hsPhiS1pi->Write();
  hsThetaS1ratio->Write();
  hsThetaS1ratioNoS2->Write();
  hsPhiS1ratio->Write();

  hsMomS1S2->Write();
  hsMomS1->Write();
  hsMomTpc->Write();
  hsKE->Write();

  legThetaS1pro->Write("legThetaS1pro");
  legThetaS1pi->Write("legThetaS1pi");
  legKE->Write("legKE");
  legTof->Write("legTof");
  legTpc->Write("legTpc");

  fout->Close();
  delete fout;
} // angularDistS3
