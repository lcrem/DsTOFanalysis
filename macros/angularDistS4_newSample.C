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

double weightErr(const double weight, const double barerr, const double n)
{
  double err = 0.;
  if (n!=0) {
    err = sqrt(pow(weight * sqrt(1. + (1. / n)), 2) + barerr*barerr);
  }
  //double err = TMath::Sqrt( weight*weight + pow(werr*(weight*weight), 2) );
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
  // For calculating ratios
  const vector<double> dataS3    = {1983., 1656., 1325., 899., 136.3};
  const vector<double> dataS3Err = {9., 6., 5., 6., 0.5};
  // For MC plots
  const double xUp = 0.44; // m
  const double xDown = -0.97; // m 

  TFile *fout = new TFile(saveDir, "recreate");

  // THStack *hsBkgSub = new THStack("hsBkgSub", "S2 to S4 time of flight spectrum; t_{S4} - t_{S2} / ns; Events / spill");
  THStack *hsDtof     = new THStack("hsDtof", "S4 ToF spectrum; t_{S4} - t_{S2} / ns; Events / spill");
  THStack *hsDtofSub  = new THStack("hsDtofSub", "S4 ToF spectrum; t_{S4} - t_{S2} / ns; Events / spill / ns");
  THStack *hsDtofBins = new THStack("hsDtofBins", "S4 ToF spectrum; t_{S4} - t_{S2} / ns; Events / spill / ns");

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

  vector<double> edgesDtofVec = getDtofEdges();
  double arr[edgesDtofVec.size()-1];
  for (int i=0; i<edgesDtofVec.size(); i++) {
    arr[i] = edgesDtofVec.at(i);
    cout<<arr[i]<<endl;
  }

  for (int nBlocks = 0; nBlocks <= 4; nBlocks++) {
    cout<<"=========================================="<<endl;
    cout<<nBlocks<<" blocks"<<endl;
    cout<<"=========================================="<<endl;
    vector<double> startTimes;
    vector<double> endTimes;
    startTimes.clear();
    endTimes.clear();
    fout->cd();
    TTree *protonTree = new TTree(Form("protonTree%d", nBlocks), Form("protonTree%d", nBlocks));
    double weight, error;
    double tof;
    double mom;
    double theta;
    double phi;
    double mcx, mcy, mcz;
    int spill;
    protonTree->Branch("weight", &weight);
    protonTree->Branch("error", &error);
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
    TH1D *hProS4Horz = new TH1D(Form("hProS4Horz%d",nBlocks), Form("Horizontal angular distribution of proton hits in S4, %d blocks; #theta / degrees; Events / spill",nBlocks), nBinsS4Horz, binsS4HorzLow, binsS4HorzHigh);
    setHistAttr(hProS4Horz);
    TH1D *hProS4HorzErr = new TH1D(Form("hProS4HorzErr%d",nBlocks), Form("Horizontal angular distribution of proton hits in S4, %d blocks; #theta / degrees; Events / spill",nBlocks), nBinsS4Horz, binsS4HorzLow, binsS4HorzHigh);
    setHistAttr(hProS4HorzErr);
    TH1D *hProS4HorzSmear = new TH1D(Form("hProS4HorzSmear%d",nBlocks), Form("Horizontal angular distribution of proton hits in S4, %d blocks; #theta / degrees; Events / spill",nBlocks), nBinsS4Horz, binsS4HorzLow, binsS4HorzHigh);
    setHistAttr(hProS4HorzSmear);
    TH1D *hPiS4Horz  = new TH1D(Form("hPiS4Horz%d",nBlocks), Form("Horizontal angular distribution of MIP hits in S4, %d blocks; #theta / degrees; Events / spill",nBlocks), nBinsS4Horz, binsS4HorzLow, binsS4HorzHigh);
    setHistAttr(hPiS4Horz);
TH1D *hPiS4HorzErr = new TH1D(Form("hPiS4HorzErr%d",nBlocks), Form("Horizontal angular distribution of MIP hits in S4, %d blocks; #theta / degrees; Events / spill",nBlocks), nBinsS4Horz, binsS4HorzLow, binsS4HorzHigh);
    setHistAttr(hPiS4HorzErr);
    TH1D *hAllS4Horz = new TH1D(Form("hAllS4Horz%d",nBlocks), Form("Horizontal angular distribution of hit in S4, %d blocks; #theta / degrees; Events / spill", nBlocks), nBinsS4Horz, binsS4HorzLow, binsS4HorzHigh);
    setHistAttr(hAllS4Horz);
    TH1D *hProPiRatioS4Horz  = new TH1D(Form("hProPiRatioS4Horz%d",nBlocks), Form("Horizontal angular distribution of proton/MIP ratio in S4, %d blocks; #theta / degrees; Protons/MIPs",nBlocks), nBinsS4Horz, binsS4HorzLow, binsS4HorzHigh);
    setHistAttr(hProPiRatioS4Horz);
    // Vertical
    TH1D *hProS4Vert = new TH1D(Form("hProS4Vert%d",nBlocks), Form("Vertical angular distribution of proton hits in S4, %d blocks; #phi / degrees; Events / spill",nBlocks), nBinsS4Vert, binsS4VertLow, binsS4VertHigh);
    setHistAttr(hProS4Vert);
    TH1D *hProS4VertErr = new TH1D(Form("hProS4VertErr%d",nBlocks), Form("Vertical angular distribution of proton hits in S4, %d blocks; #phi / degrees; Events / spill",nBlocks), nBinsS4Vert, binsS4VertLow, binsS4VertHigh);
    setHistAttr(hProS4VertErr);
    TH1D *hPiS4Vert  = new TH1D(Form("hPiS4Vert%d",nBlocks), Form("Vertical angular distribution of MIP hits in S4, %d blocks; #phi / degrees; Events / spill",nBlocks), nBinsS4Vert, binsS4VertLow, binsS4VertHigh);
    setHistAttr(hPiS4Vert);
    TH1D *hPiS4VertErr = new TH1D(Form("hPiS4VertErr%d",nBlocks), Form("Vertical angular distribution of MIP hits in S4, %d blocks; #phi / degrees; Events / spill",nBlocks), nBinsS4Vert, binsS4VertLow, binsS4VertHigh);
    setHistAttr(hPiS4VertErr);
    TH1D *hProPiRatioS4Vert  = new TH1D(Form("hProPiRatioS4Vert%d",nBlocks), Form("Vertical angular distribution of proton/MIP ratio in S4, %d blocks; #phi / degrees; Protons/MIPs",nBlocks), nBinsS4Vert, binsS4VertLow, binsS4VertHigh);
    setHistAttr(hProPiRatioS4Vert);
    TH1D *hMSq = new TH1D(Form("hMSq%d", nBlocks), Form("Particle mass distribution, %d blocks; M^{2} [GeV^{2} / c^{2}]; Events / spills", nBlocks), 109, -0.3, 4.5);
    setHistAttr(hMSq);
    TH1D *hMSqErr = new TH1D(Form("hMSqErr%d", nBlocks), Form("Particle mass distribution, %d blocks; M^{2} [GeV^{2} / c^{2}]; Events / spills", nBlocks), 109, -0.3, 4.5);
    setHistAttr(hMSqErr);
    TH2D *h2DoubleFirst = new TH2D(Form("h2DoubleFirst%d", nBlocks), Form("%d blocks: S4 double hits, 1st bar; Bar position / cm; Bar; Events / spill", nBlocks), nBinsS4Horz, 0., 140., 10, 0.5, 10.5);
    setHistAttr(h2DoubleFirst);
    TH2D *h2DoubleSecond = new TH2D(Form("h2DoubleSecond%d", nBlocks), Form("%d blocks: S4 double hits, 2nd bar; Bar position / cm; Bar; Events / spill", nBlocks), nBinsS4Horz, 0., 140., 10, 0.5, 10.5);
    setHistAttr(h2DoubleSecond);

    // No weighting
    TH1D *hProS4HorzUnwgt = new TH1D(Form("hProS4HorzUnwgt%d",nBlocks), Form("Horizontal angular distribution of proton hits in S4, %d blocks; #theta / degrees; Events",nBlocks), nBinsS4Horz, binsS4HorzLow, binsS4HorzHigh);
    setHistAttr(hProS4HorzUnwgt);
    TH1D *hPiS4HorzUnwgt  = new TH1D(Form("hPiS4HorzUnwgt%d",nBlocks), Form("Horizontal angular distribution of MIP hits in S4, %d blocks; #theta / degrees; Events",nBlocks), nBinsS4Horz, binsS4HorzLow, binsS4HorzHigh);
    setHistAttr(hPiS4HorzUnwgt);
    TH1D *hProS4VertUnwgt = new TH1D(Form("hProS4VertUnwgt%d",nBlocks), Form("Vertical angular distribution of proton hits in S4, %d blocks; #phi / degrees; Events",nBlocks), nBinsS4Vert, binsS4VertLow, binsS4VertHigh);
    setHistAttr(hProS4VertUnwgt);
    TH1D *hPiS4VertUnwgt  = new TH1D(Form("hPiS4VertUnwgt%d",nBlocks), Form("Vertical angular distribution of MIP hits in S4, %d blocks; #phi / degrees; Events",nBlocks), nBinsS4Vert, binsS4VertLow, binsS4VertHigh);
    setHistAttr(hPiS4VertUnwgt);
    TH1D *hProPiRatioS4HorzUnwgt  = new TH1D(Form("hProPiRatioS4HorzUnwgt%d",nBlocks), Form("Horizontal angular distribution of proton/MIP ratio in S4, %d blocks; #theta / degrees; Protons/MIPs",nBlocks), nBinsS4Horz, binsS4HorzLow, binsS4HorzHigh);
    TH1D *hProPiRatioS4VertUnwgt  = new TH1D(Form("hProPiRatioS4VertUnwgt%d",nBlocks), Form("Vertical angular distribution of proton/MIP ratio in S4, %d blocks; #phi / degrees; Protons/MIPs",nBlocks), nBinsS4Vert, binsS4VertLow, binsS4VertHigh);
    TH2D *h2dAngProS1Unwgt = new TH2D(Form("h2dAngProS1Unwgt%d", nBlocks), Form("S4 angular distribution of proton hits, %d blocks; #theta / degrees; #phi / degrees; Events", nBlocks), nBinsS4Horz, binsS4HorzLow, binsS4HorzHigh, nBinsS4Vert, binsS4VertLow, binsS4VertHigh);
    setHistAttr(h2dAngProS1Unwgt);
    TH2D *h2dAngPiS1Unwgt = new TH2D(Form("h2dAngPiS1Unwgt%d", nBlocks), Form("S4 angular distribution of MIP hits, %d blocks; #theta / degrees; #phi / degrees; Events", nBlocks), nBinsS4Horz, binsS4HorzLow, binsS4HorzHigh, nBinsS4Vert, binsS4VertLow, binsS4VertHigh);
    setHistAttr(h2dAngPiS1Unwgt);
    TH2D *h2dAngRatioS1Unwgt = new TH2D(Form("h2dAngRatioS1Unwgt%d", nBlocks), Form("S4 angular distribution of proton/MIP ratio, %d blocks; #theta / degrees; #phi / degrees; Ratio", nBlocks), nBinsS4Horz, binsS4HorzLow, binsS4HorzHigh, nBinsS4Vert, binsS4VertLow, binsS4VertHigh);
    setHistAttr(h2dAngRatioS1Unwgt);

    // Proton momentum
    TH1D *hMom = new TH1D(Form("hMom%d", nBlocks), Form("%d blocks; Momentum [MeV / c]; Events / spill", nBlocks), 75, 100, 850);
    setHistAttr(hMom);
    TH1D *hMomErr = new TH1D(Form("hMomErr%d", nBlocks), Form("%d blocks; Momentum [MeV / c]; Events / spill", nBlocks), 75, 100, 850);
    setHistAttr(hMomErr);
    TH1D *hMomS1S4 = new TH1D(Form("hMomS1S4%d", nBlocks), Form("%d blocks; Momentum [MeV / c]; Events / spill", nBlocks), 75, 100, 850);
    setHistAttr(hMomS1S4);

    TH1D *hSmearWeight = new TH1D(Form("hSmearWeight%d", nBlocks), Form("Smear weights, %d blocks", nBlocks), nBinsS4Horz, xDown, xUp);
    setHistAttr(hSmearWeight);

    TH2D *h2dAngProS1 = new TH2D(Form("h2dAngProS1%d", nBlocks), Form("S4 angular distribution of proton hits, %d blocks; #theta / degrees; #phi / degrees; Events / spill", nBlocks), nBinsS4Horz, binsS4HorzLow, binsS4HorzHigh, nBinsS4Vert, binsS4VertLow, binsS4VertHigh);
    TH2D *h2dAngPiS1 = new TH2D(Form("h2dAngPiS1%d", nBlocks), Form("S4 angular distribution of MIP hits, %d blocks; #theta / degrees; #phi / degrees; Events / spill", nBlocks), nBinsS4Horz, binsS4HorzLow, binsS4HorzHigh, nBinsS4Vert, binsS4VertLow, binsS4VertHigh);
    TH2D *h2dAngRatioS1 = new TH2D(Form("h2dAngRatioS1%d", nBlocks), Form("S4 angular distribution of proton/MIP ratio, %d blocks; #theta / degrees; #phi / degrees; Events / spill", nBlocks), nBinsS4Horz, binsS4HorzLow, binsS4HorzHigh, nBinsS4Vert, binsS4VertLow, binsS4VertHigh);

    TH2D *h2dAngProWC = new TH2D(Form("h2dAngProWC%d", nBlocks), Form("S4 angular distribution of proton hits, %d blocks; #theta / degrees; #phi / degrees; Events / spill", nBlocks), nBinsS4Horz, binsS4HorzLow, binsS4HorzHigh, nBinsS4Vert, binsS4VertLow, binsS4VertHigh);
    TH2D *h2dAngPiWC = new TH2D(Form("h2dAngPiWC%d", nBlocks), Form("S4 angular distribution of MIP hits, %d blocks; #theta / degrees; #phi / degrees; Events / spill", nBlocks), nBinsS4Horz, binsS4HorzLow, binsS4HorzHigh, nBinsS4Vert, binsS4VertLow, binsS4VertHigh);
    TH2D *h2dAngRatioWC = new TH2D(Form("h2dAngRatioWC%d", nBlocks), Form("S4 angular distribution of proton/MIP ratio, %d blocks; #theta / degrees; #phi / degrees; Events / spill", nBlocks), nBinsS4Horz, binsS4HorzLow, binsS4HorzHigh, nBinsS4Vert, binsS4VertLow, binsS4VertHigh);
    TH2D *h2dProMCComp = new TH2D(Form("h2dProMCComp%d", nBlocks), Form("Proton distribution, %d blocks; x / m; y / m; Events / spill", nBlocks), nBinsS4Horz, -0.9805, 0.4183, 10, -0.3406, 0.4434);
    setHistAttr(h2dProMCComp);
    TH2D *h2dProMCCompCut = new TH2D(Form("h2dProMCCompCut%d", nBlocks), Form("Proton distribution, %d blocks; x / m; y / m; Events / spill", nBlocks), nBinsS4Horz, -0.8805, 0.3183, 9, -0.3406, 0.3684);
    setHistAttr(h2dProMCCompCut);
 
    // 1D ToF 
    TH1D *hdtof1d    = new TH1D(Form("hdtof1d%d",nBlocks), Form("Time of flight, measured in S4, %d blocks; t_{S4} - t_{S2} / ns; Events / spill", nBlocks), nBinsDtof, binsDtofLow, binsDtofHigh);
    setHistAttr(hdtof1d);
    TH1D *hdtof1dErr = new TH1D(Form("hdtof1dErr%d",nBlocks), Form("Time of flight, measured in S4, %d blocks; t_{S4} - t_{S2} / ns; Events / spill", nBlocks), nBinsDtof, binsDtofLow, binsDtofHigh);
    setHistAttr(hdtof1dErr);
    TH1D *hdtof1dBins = new TH1D(Form("hdtof1dBins%d",nBlocks), Form("Time of flight, measured in S4, %d blocks; t_{S4} - t_{S2} / ns; Events / spill / ns", nBlocks), edgesDtofVec.size()-1, arr);
    setHistAttr(hdtof1dBins);
    TH1D *hdtof1dBinsErr = new TH1D(Form("hdtof1dBinsErr%d",nBlocks), Form("Time of flight, measured in S4, %d blocks; t_{S4} - t_{S2} / ns; Events / spill", nBlocks), edgesDtofVec.size()-1, arr);
    setHistAttr(hdtof1dBinsErr);
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
      TH2D *h2Cosmics = new TH2D(Form("h2Cosmics%d",nBlocks), Form("Cosmic flux, %d blocks; x / cm; Bar; Rate / s^{-1}", nBlocks), nBinsS4Horz, 0, 140, 10, 0.5, 10.5);
      setHistAttr(h2Cosmics);
      TH2D *h2CosmicsEff = new TH2D(Form("h2CosmicsEff%d", nBlocks), Form("Relative efficiency, %d blocks; x / cm; Bar; Eff", nBlocks), nBinsS4Horz, 0, 140, 10, 0.5, 10.5);
      setHistAttr(h2CosmicsEff);
      TH2D *h2CosmicsErr = new TH2D(Form("h2CosmicsErr%d", nBlocks), Form("Relative efficiency, %d blocks; x / cm; Bar; Err", nBlocks), nBinsS4Horz, 0, 140, 10, 0.5, 10.5);
      setHistAttr(h2CosmicsErr);
      TH2D *h2CosmicsEffBins = new TH2D(Form("h2CosmicsEffBins%d", nBlocks), Form("Relative efficiency, %d blocks; x / cm; Bar; Eff", nBlocks), binnum, binsCosmics, 10, 0.5, 10.5);
      setHistAttr(h2CosmicsEffBins);
      TH2D *h2CosmicsEffMC = new TH2D(Form("h2CosmicsEffMC%d", nBlocks), Form("Relative efficiency, %d blocks; x / cm; y / cm; Eff", nBlocks), nBinsS4Horz, -0.95, 0.41, 10, -0.3361, 0.4139);
      setHistAttr(h2CosmicsEffMC);
      TH1D *hCosmicsVertEff = new TH1D(Form("hCosmicsVertEff%d",nBlocks), Form("With cosmics, %d blocks; x / cm; Eff",nBlocks), 10, 0.5, 10.5);
      setHistAttr(hCosmicsVertEff);
      TH1D *hCosmicsVertEff2dNorm = new TH1D(Form("hCosmicsVertEff2dNorm%d",nBlocks), Form("With cosmics, %d blocks; x / cm; Eff",nBlocks), 10, 0.5, 10.5);
      setHistAttr(hCosmicsVertEff2dNorm);
      TH1D *hCosmicsHorz = new TH1D(Form("hCosmicsHorz%d",nBlocks), Form("Cosmic flux, %d blocks; x / cm; Rate / Hz",nBlocks), nBinsS4Horz, 0, 140);
      setHistAttr(hCosmicsHorz);

      // Bar by bar efficiencies calculated by cosmics
      std::vector<TH1D*> eff1dVec;
      for (int b=0; b<10; b++) {
	TH1D *hEff1d = new TH1D(Form("hEff1d_block%d_bar%d_sub%d", nBlocks, b+1, sub), Form("%d blocks, bar %d: Efficiency; x / cm; Eff", nBlocks, b+1), nBinsS4Horz, 0., 140.);
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

      vector<double> barTimeVec;
      barTimeVec.resize(10, 0.);

      for (int itdc=0; itdc<2; itdc++) {
	// Need both the coincidence and the raw trees for this
	double tempUstof;
	double ustofNs;
	unsigned int lastRun = 0;
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
	    if (dstofHitT - barTimeVec.at(tofCoin->bar-1) < s4DeadtimeCut && tofCoin->run==lastRun) {
	      continue;
	    }
	    else if (tofCoin->run != lastRun) {
	      barTimeVec.clear();
	      barTimeVec.resize(10, 0.);
	      lastRun = tofCoin->run;
	    }
	    barTimeVec.at(tofCoin->bar-1) = dstofHitT;

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
		h2CosmicsErr->Fill(positionX, tofCoin->bar);
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

      for (int ii=0; ii<h2CosmicsEff->GetNbinsX()+1; ii++) {
	for (int jj=0; jj<h2CosmicsEff->GetNbinsY()+1; jj++) {
	  h2CosmicsErr->SetBinContent(ii, jj, h2CosmicsEff->GetBinError(ii, jj)); 
	}
      }

      hHits->Write();
      hCoins->Write();
      hHits->Add(hCoins, -1.);
      hEff->Divide(hCoins, hHits, 1., 1., "B");
      h2CosmicsEff->Write();
      h2CosmicsErr->Write();
      h2CosmicsEffBins->Write();
      h2CosmicsEffMC->Write();
      h2Cosmics->Write();
      // Now loop over the coincidence files again and calculate the angular distributions
      cout<<"Getting signal hits"<<endl;
      for (int itdc=0; itdc<2; itdc++) {
	unsigned int lastRun = 0;
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
	  double dstofHitT = dtofHitTime(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]);

	  // Make bar go dead after hit
	  if (dstofHitT - barTimeVec.at(tofCoin->bar-1) < s4DeadtimeCut && tofCoin->run==lastRun) {
	    continue;
	  }
	  else if (tofCoin->run != lastRun) {
	    barTimeVec.clear();
	    barTimeVec.resize(10, 0.);
	    lastRun = tofCoin->run;
	  }
	  barTimeVec.at(tofCoin->bar-1) = dstofHitT;

	  double tofCalc = dstofHitT - tofCoin->usTofSignal - dstofShift;
	  if (tofCalc < proCutHiS4 && tofCalc > 30. && tofCoin->bar != 10 && deltat < s4BarTime) {
	    double positionXP = localDtofPosition(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]);
	    TVector3 globalCoords = GetDtofGlobalCoords(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1], tofCoin->bar);
	    double mcXForSmear = globalToMCCoords(globalCoords).X();

	    double w = 1. / (h2CosmicsEff->GetBinContent(h2CosmicsEff->GetXaxis()->FindBin(positionXP), tofCoin->bar)*barOverallEff); 
	    // Need to count if there is a double hit associated with this one
	    int bar1 = tofCoin->bar;
	    for (int dub = h-1; dub<=h+1; dub+=2) {
	      tofCoinChain->GetEntry(dub);
	      double dstofHitT2 = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (s4BarTime/2. - TMath::Abs(deltat) / 2.);
	      int bar2 = tofCoin->bar;
	      if (dstofHitT2 - dstofHitT < 2. && dstofHitT2 > dstofHitT && abs(bar1-bar2) % 2 == 1) {
		h2DoubleFirst->Fill(positionXP, bar1);
		double positionXP2 = (((tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0])*(70./s4BarTime) + 70.));
		h2DoubleSecond->Fill(positionXP2, bar2);
		doubleHits++;
		// If there is a double hit we should halve the weight for this event
		w /= 2;
		break;
	      }
	    } // Loop to find double hits

	    double errSq    = pow(weightErr(w, barOverallEffErr, h2Cosmics->GetBinContent(h2Cosmics->GetXaxis()->FindBin(positionXP), tofCoin->bar)), 2);
	    double errSqBar = pow(weightErrBar(hEff->GetBinContent(tofCoin->bar), hEff->GetBinError(tofCoin->bar)), 2);
	    // Calculate position of hit in global coordinates
	    TVector3 dtofCoords = GetDtofGlobalCoords(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1], tofCoin->bar);
	    // Calculate the angles relative to the nominal beamline
	    double angleTheta = getThetaFromGlobal(dtofCoords);
	    double anglePhi   = getPhiFromGlobal(dtofCoords);
	    // Apply offsets to get these in same frame as MC
	    TVector3 mcCoords = globalToMCCoords(dtofCoords);
	    double mcX = mcCoords.X();
	    double mcY = mcCoords.Y();
	    double mcZ = mcCoords.Z();
	    hdtof1d->Fill(tofCalc, w);
	    hdtof1dErr->Fill(tofCalc, errSq);
	    hdtof1dBins->Fill(tofCalc, w);
	    hdtof1dBinsErr->Fill(tofCalc, errSq);
	    hMSq->Fill(massFromTime(tofCalc, 0.8, s2s4Dist), w);
	    hMSqErr->Fill(massFromTime(tofCalc, 0.8, s2s4Dist), errSq);

	    if (tofCalc < piHiS4 & tofCalc > piLowS4) { 
	      nPi += w;
	      hPiS4Horz->Fill(angleTheta, w);
	      hPiS4HorzErr->Fill(angleTheta, errSq);
	      hPiS4Vert->Fill(anglePhi, w);
	      hPiS4VertErr->Fill(anglePhi, errSq);
	      h2dAngPiS1->Fill(angleTheta, anglePhi, w);

	      hPiS4HorzUnwgt->Fill(angleTheta);
	      hPiS4VertUnwgt->Fill(anglePhi);
	      h2dAngPiS1Unwgt->Fill(angleTheta, anglePhi);
	    } // if (tofCalc < piHiS4 & tofCalc > piLowS4)
	    else if (tofCalc < proCutHiS4 & tofCalc > proCutLowS4) {
	      tof = tofCalc;
	      mom = momFromTime(0.938, s2s4Dist, tofCalc);
	      weight = w;
	      error = sqrt(errSq);
	      theta = angleTheta;
	      phi = anglePhi;
	      mcx = mcX;
	      mcy = mcY;
	      mcz = mcZ;
	      spill = nSpillsTrue;
	      protonTree->Fill();

	      nP += w;
	      hProS4Horz->Fill(angleTheta, w);
	      hProS4HorzErr->Fill(angleTheta, errSq);
	      hProS4Vert->Fill(anglePhi, w);
	      hProS4VertErr->Fill(anglePhi, errSq);
	      hProS4HorzUnwgt->Fill(angleTheta);
	      hProS4VertUnwgt->Fill(anglePhi);
	      h2dAngProS1->Fill(angleTheta, anglePhi, w);
	      h2dAngProS1Unwgt->Fill(angleTheta, anglePhi);
	      hMom->Fill(momFromTime(0.938, s2s4Dist, tofCalc), w);
	      hMomErr->Fill(momFromTime(0.938, s2s4Dist, tofCalc), errSq);
	      hMomS1S4->Fill(momFromTime(0.938, s1s4Dist, tofCalc+4.7), w);
	      h2dProMCComp->Fill(mcX, mcY, w);
	      if (positionXP > 10. && positionXP < 130.) {
		h2dProMCCompCut->Fill(mcX, mcY, w);
	      }
	      
	    } // else if (tofCalc < proHi & tofCalc > proLow) 
	  } // if (tofCalc < 160. && tofCalc > 30.) 
	} // for (int h=0; h<tofCoinChain->GetEntries(); h++)
	delete tofCoin;
	delete tofCoinChain;
      } // Loop over TDCs
      cout<<"Finished getting signal hits"<<endl;
      cout<<"Found "<<doubleHits<<" double proton hits"<<endl;;
      cout<<"Finished subsample "<<sub<<endl;

    } // Loop over subsamples

    for (int bin=0; bin < hdtof1d->GetNbinsX()+1; bin++) {
      hdtof1d->SetBinError(bin, sqrt(hdtof1dErr->GetBinContent(bin)));
    }
    for (int bin=0; bin < hdtof1dBins->GetNbinsX()+1; bin++) {
      hdtof1dBins->SetBinError(bin, sqrt(hdtof1dBinsErr->GetBinContent(bin)));
    }
    for (int bin=0; bin < hMSq->GetNbinsX()+1; bin++) {
      hMSq->SetBinError(bin, sqrt(hMSqErr->GetBinContent(bin)));
    }
    for (int bin=0; bin < hProS4Horz->GetNbinsX()+1; bin++) {
      hProS4Horz->SetBinError(bin, sqrt(hProS4HorzErr->GetBinContent(bin)));
      hPiS4Horz->SetBinError(bin, sqrt(hPiS4HorzErr->GetBinContent(bin)));
    }
    for (int bin=0; bin < hProS4Vert->GetNbinsX()+1; bin++) {
      hProS4Vert->SetBinError(bin, sqrt(hProS4VertErr->GetBinContent(bin)));
      hPiS4Vert->SetBinError(bin, sqrt(hPiS4VertErr->GetBinContent(bin)));
    }
    for (int bin=0; bin<hMom->GetNbinsX()+1; bin++) {
      hMom->SetBinError(bin, sqrt(hMomErr->GetBinContent(bin)));
    }

    hdtof1d->Scale(1. / nSpillsTrue);
    hdtof1dBins->Scale(1. / nSpillsTrue);
    hdtof1dBins->Scale(1., "width");
    fout->cd();
    TF1 *fSplusB4 = new TF1(Form("signal_plus_bkg_%d", nBlocks), "gaus(0)+pol0(3)", 30, proCutHiS4);
    if (nBlocks != 4) {
      hdtof1dBins->Fit(sPi, "R");
      hdtof1dBins->Fit(sPro, "R");
      hdtof1dBins->Fit(fBkg, "R");
      Double_t par[7];
      sPi->GetParameters(&par[0]);
      sPro->GetParameters(&par[3]);
      fBkg->GetParameters(&par[6]);
      fSplusB->SetParameters(par);
      hdtof1dBins->Fit(fSplusB, "R");
      sPro->Write();
      sPi->Write();
      fBkg->Write();
      fSplusB->Write();
    }
    else {
      hdtof1dBins->Fit(sPi, "R");
      hdtof1dBins->Fit(fBkg, "R");
      Double_t par[4];
      sPi->GetParameters(&par[0]);
      fBkg->GetParameters(&par[3]);
      fSplusB4->SetParNames("piConst", "piMean", "piSigma", "bkg");
      fSplusB4->SetLineColor(kRed);
      fSplusB4->SetParameters(par);
      hdtof1dBins->Fit(fSplusB4, "R");
      sPi->Write();
      fBkg->Write();
      fSplusB4->Write();
    }

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

    hProS4VertErr->Scale(1. / nSpillsTrue);
    hProS4HorzErr->Scale(1. / nSpillsTrue);
    hPiS4VertErr->Scale(1. / nSpillsTrue);
    hPiS4HorzErr->Scale(1. / nSpillsTrue);
    hMom->Scale(1. / nSpillsTrue);

    hProS4VertErr->Write();
    hProS4HorzErr->Write();
    hPiS4VertErr->Write();
    hPiS4HorzErr->Write();
    hMom->Write();
    hdtof1d->Write();
    hdtof1dBins->Write();

    // Now we have the fit values, loop over again and subtract the background
    // For each bin, find the fraction of each particle type which are background and 
    // this fraction from the bin
    TF1 *fSub = new TF1("fSub", "pol0", 30, proCutHiS4);
    (nBlocks !=4) ? fSub->SetParameter(0, fSplusB->GetParameter("bkg")) : fSub->SetParameter(0, fSplusB4->GetParameter("bkg"));  

    // Background subtracted tof spectrum
    TH1D *hdtof1dBins_sub = (TH1D*)hdtof1dBins->Clone(Form("hdtof1dBins_sub%d",nBlocks));
    hdtof1dBins_sub->Add(fSub, -1.);
    // If the bin content drops below 0, set to 0
    for (int b = 0; b <=  hdtof1dBins_sub->GetNbinsX(); b++) {
      if (hdtof1dBins_sub->GetBinContent(b) < 0.) {
	hdtof1dBins_sub->SetBinContent(b, 0);
	hdtof1dBins_sub->SetBinError(b, 0);
      } // if (hdtof1d_sub->GetBinContent(b) < 0.)
    } // for (int b = 0; hdtof1d_sub->GetNbinsX(); b++)
    
    // Integrate background function between proton and pion windows and then subtract
    double bkgPerNs = 0.;
    double bkgErr = 0.;
    if (nBlocks==4) {
      bkgErr   = fSplusB4->GetParError(3);
      bkgPerNs = fSplusB4->GetParameter("bkg");
    }
    else {
      bkgErr   = fSplusB->GetParError(6);
      bkgPerNs = fSplusB->GetParameter("bkg");
    }
    
    nP  /= (double)nSpillsTrue;
    nPi /= (double)nSpillsTrue;

    cout<<"Bkg / spill / ns "<<bkgPerNs<<endl;
    double piBkg    = (piHiS4 - piLowS4)*bkgPerNs;
    double piBkgErr = (piHiS4 - piLowS4)*bkgErr;
    double proBkg    = (proCutHiS4 - proCutLowS4)*bkgPerNs;
    double proBkgErr = (proCutHiS4 - proCutLowS4)*bkgErr;
    cout<<"Pion background: "<<piBkg<<" +- "<<piBkgErr<<" / spill. Proton background: "<<proBkg<<" +- "<<proBkgErr<<" / spill, compared to "<<nP<<" total"<<endl;
    cout<<"Spills "<<nSpills<<" ("<<nSpillsTrue<<" true)"<<endl;

    // Subtract background hits
    // Do this by just multiplying by background/integral
    for (int b=0; b <= hProS4Horz->GetNbinsX(); b++) {
      hProS4Horz->SetBinContent(b, hProS4Horz->GetBinContent(b) * (1 - proBkg/nP));
      if (hProS4Horz->GetBinContent(b) < 0.) {
    	hProS4Horz->SetBinContent(b, 0);
    	hProS4Horz->SetBinError(b, 0);
      }
    }

    for (int b=0; b <= hPiS4Horz->GetNbinsX(); b++) {
      hPiS4Horz->SetBinContent(b, hPiS4Horz->GetBinContent(b) * (1 - piBkg/nPi));
      if (hPiS4Horz->GetBinContent(b) < 0.) {
    	hPiS4Horz->SetBinError(b, 0);
    	hPiS4Horz->SetBinContent(b, 0);
      }
    }

    for (int b=0; b <= hProS4Vert->GetNbinsX(); b++) {
      hProS4Vert->SetBinContent(b, hProS4Vert->GetBinContent(b) * (1 - proBkg/nP));
      if (hProS4Vert->GetBinContent(b) < 0.) {
    	hProS4Vert->SetBinError(b, 0);
    	hProS4Vert->SetBinContent(b, 0);
      }
    }

    for (int b=0; b<=hPiS4Vert->GetNbinsX(); b++) {
      hPiS4Vert->SetBinContent(b, hPiS4Vert->GetBinContent(b) * (1 - piBkg/nPi));
      if (hPiS4Vert->GetBinContent(b) < 0.) {
    	hPiS4Vert->SetBinError(b, 0);
    	hPiS4Vert->SetBinContent(b, 0);
      }
    }

    for (int b=0; b<=hMSq->GetNbinsX(); b++) {
      hMSq->SetBinContent(b, hMSq->GetBinContent(b) * (1 - piBkg/nPi));
      if (hMSq->GetBinContent(b) < 0.) {
    	hMSq->SetBinError(b, 0);
    	hMSq->SetBinContent(b, 0);
      }
    }

    hProPiRatioS4Vert->Divide(hProS4Vert, hPiS4Vert, 1., 1., "B");
    hProPiRatioS4Horz->Divide(hProS4Horz, hPiS4Horz, 1., 1., "B");
    hProPiRatioS4VertUnwgt->Divide(hProS4VertUnwgt, hPiS4VertUnwgt, 1., 1., "B");
    hProPiRatioS4HorzUnwgt->Divide(hProS4HorzUnwgt, hPiS4HorzUnwgt, 1., 1., "B");

    double eProS4Horz = 0;
    double ePiS4Horz  = 0;
    double intProS4Horz = hProS4Horz->IntegralAndError(1, hProS4Horz->GetNbinsX(), eProS4Horz, "width");
    double intPiS4Horz  = hPiS4Horz->IntegralAndError(1, hPiS4Horz->GetNbinsX(), ePiS4Horz, "width");
    eProS4Horz = sqrt(eProS4Horz*eProS4Horz + proBkgErr*proBkgErr);
    ePiS4Horz  = sqrt(ePiS4Horz*ePiS4Horz + piBkgErr*piBkgErr);

    if (nBlocks == 0) {
      hdtof1d->SetLineColor(kBlack);
      hdtof1dBins->SetLineColor(kBlack);
      hdtof1dBins_sub->SetLineColor(kBlack);
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
      leg->AddEntry(hdtof1d, "0 blocks", "l");     
      legProS4Horz->AddEntry(hProS4Horz, Form("0 blocks - %.3g #pm %.1g per spill", intProS4Horz, eProS4Horz), "l");     
      legPiS4Horz->AddEntry(hPiS4Horz, Form("0 blocks - %d #pm %d per spill", (int)intPiS4Horz, (int)ePiS4Horz), "l");
      legRatioVert->AddEntry(hProPiRatioS4Vert, "0 blocks", "l");
    }
    else if (nBlocks == 1) {      
      hdtof1d->SetLineColor(kRed);
      hdtof1dBins->SetLineColor(kRed);
      hdtof1dBins_sub->SetLineColor(kRed);
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
      leg->AddEntry(hdtof1d, "1 block", "l");
      legProS4Horz->AddEntry(hProS4Horz, Form("1 block - %.3g #pm %.1g per spill", intProS4Horz, eProS4Horz), "l");     
      legPiS4Horz->AddEntry(hPiS4Horz, Form("1 block - %d #pm %d per spill", (int)intPiS4Horz, (int)ePiS4Horz), "l");     
      legRatioVert->AddEntry(hProPiRatioS4Vert, "1 block", "l");
    }
    else if (nBlocks == 2) {
      hdtof1d->SetLineColor(kBlue);
      hdtof1dBins->SetLineColor(kBlue);
      hdtof1dBins_sub->SetLineColor(kBlue);
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
      leg->AddEntry(hdtof1d, "2 blocks", "l");
      legProS4Horz->AddEntry(hProS4Horz, Form("2 blocks - %.3g #pm %.1g per spill", intProS4Horz, eProS4Horz), "l");     
      legPiS4Horz->AddEntry(hPiS4Horz, Form("2 blocks - %d #pm %d per spill", (int)intPiS4Horz, (int)ePiS4Horz), "l");
      legRatioVert->AddEntry(hProPiRatioS4Vert, "2 blocks", "l");
    }
    else if (nBlocks == 3) {
      hdtof1d->SetLineColor(kCyan+1);
      hdtof1dBins->SetLineColor(kCyan+1);
      hdtof1dBins_sub->SetLineColor(kCyan+1);
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
      leg->AddEntry(hdtof1d, "3 blocks", "l");
      legProS4Horz->AddEntry(hProS4Horz, Form("3 blocks - %.3g #pm %.1g per spill", intProS4Horz, eProS4Horz), "l");     
      legPiS4Horz->AddEntry(hPiS4Horz, Form("3 blocks - %d #pm %d per spill", (int)intPiS4Horz, (int)ePiS4Horz), "l");
      legRatioVert->AddEntry(hProPiRatioS4Vert, "3 blocks", "l");
    }
    else if (nBlocks == 4) {
      hdtof1d->SetLineColor(kOrange+1);
      hdtof1dBins->SetLineColor(kOrange+1);
      hdtof1dBins_sub->SetLineColor(kOrange+1);
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
      leg->AddEntry(hdtof1d, "4 blocks", "l");
      legProS4Horz->AddEntry(hProS4Horz, Form("4 blocks - %.3g #pm %.1g per spill", intProS4Horz, eProS4Horz), "l");     
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
    cout<<"S4 protons per spill = "<<intProS4Horz<<" +- "<<eProS4Horz<<endl;
    cout<<"Ratio = "<<(intProS4Horz/dataS3.at(nBlocks))<<" +- "<<rerr<<endl;

    hsDtof->Add(hdtof1d);
    hsDtofSub->Add(hdtof1dBins_sub);
    hsDtofBins->Add(hdtof1dBins);
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
    hsEffComp2dNorm->Add(hEffTotal);
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
  } // for (int nBlocks = 0; nBlocks < 4; nBlocks++) 

  fout->cd();
  hsDtofSub->Write();
  hsEff->Write();
  hsDtof->Write();
  hsDtofBins->Write();
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
