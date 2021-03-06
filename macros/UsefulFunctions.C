// S1
const TVector3 D1_ULB(0.0617, -0.0086, -1.7672);
const TVector3 D1_ULT(0.0615, 0.0043, -1.7681);
const TVector3 D1_URT(-0.0102, 0.0041, -1.7655);
const TVector3 D1_URB(-0.0103, -0.0087, -1.7658);
// S2
const TVector3 D2_UTL(0.0344, 0.0706, -0.3472);
const TVector3 D2_UTR(-0.0725, 0.0679, -0.3498);
const TVector3 D2_UBL(0.0344, -0.0494, -0.3472);
const TVector3 D2_UBR(-0.0725, -0.0521, -0.3498);
const TVector3 D2_C(-0.0133, -0.2928, -0.3503); // Not the centre of S2
// S3  
const TVector3 D3_UTL(0.5215, 0.6244, 9.0650);
const TVector3 D3_UTR(-0.9928, 0.6220, 8.9245);
const TVector3 D3_UBL(0.5120, -0.5993, 9.0488);
const TVector3 D3_UBR(-1.0012, -0.6012, 8.9047);
// S3 calculated edges
const TVector3 S3_TL(0.6012, 0.6244, 9.0724);
const TVector3 S3_TR(-1.0725, 0.6220, 8.9171);
const TVector3 S3_BL(0.5917, -0.5993, 9.0562);
const TVector3 S3_BR(-1.0809, -0.6012, 8.8973);
const TVector3 S3_C = (S3_TL + S3_TR + S3_BL + S3_BR)*0.25;

// S4 frame
const TVector3 D5_UTL(-0.1227, 0.7751, 12.248);
const TVector3 D5_UBL(-0.1194, -0.8390, 12.2258);
const TVector3 D5_UBR(-1.4233, -0.8456, 12.1679);
const TVector3 D5_UTR(-1.4214, 0.7774, 12.1773);
// S4 bar-by-bar
// Bar 10
const TVector3 D5_U1L(-0.2796, 0.4325, 12.1776);
const TVector3 D5_U1R(-1.1621, 0.4312, 12.1326);
// Bar 9
const TVector3 D5_D1L(-0.4040, 0.3483, 12.3797);
const TVector3 D5_D1R(-1.1282, 0.3510, 12.3385);
// Bar 8
const TVector3 D5_U2L(-0.5110, 0.2803, 12.1542);
const TVector3 D5_U2R(-1.1521, 0.2788, 12.1090);
// Bar 7
const TVector3 D5_D2L(-0.4097, 0.2015, 12.4030);
const TVector3 D5_D2R(-1.2043, 0.2027, 12.3657);
// Bar 6
const TVector3 D5_U3L(-0.4854, 0.1300, 12.1861);
const TVector3 D5_U3R(-1.1489, 0.1308, 12.1562);
// Bar 5
const TVector3 D5_D3L(-0.4573, 0.0527, 12.4040);
const TVector3 D5_D3R(-1.1475, 0.0536, 12.3721);
// Bar 4
const TVector3 D5_U4L(-0.4625, -0.0196, 12.1580);
const TVector3 D5_U4R(-1.0993, -0.0176, 12.1265);
// Bar 3
const TVector3 D5_D4L(-0.3039, -0.1032, 12.3952);
const TVector3 D5_D4R(-1.0865, -0.0953, 12.3544);
// Bar 2
const TVector3 D5_U5L(-0.4575, -0.1724, 12.1594);
const TVector3 D5_U5R(-1.0901, -0.1719, 12.1220);
// Bar 1
const TVector3 D5_D5L(-0.3797, -0.2530, 12.4000);
const TVector3 D5_D5R(-1.1521, -0.2506, 12.3576);
const std::vector<TVector3> S4BarsL = {D5_D5L, D5_U5L, D5_D4L, D5_U4L, D5_D3L, D5_U3L, 
				       D5_D2L, D5_U2L, D5_D1L, D5_U1L};
const std::vector<TVector3> S4BarsR = {D5_D5R, D5_U5R, D5_D4R, D5_U4R, D5_D3R, D5_U3R, 
				       D5_D2R, D5_U2R, D5_D1R, D5_U1R};
// S4 calculated
const TVector3 S4_TL(-0.121, 0.4325, 12.2369);
const TVector3 S4_TR(-1.4224, 0.4325, 12.1726);
const TVector3 S4_BL(-0.121, -0.3530, 12.2369);
const TVector3 S4_BR(-1.4224, -0.3530, 12.1726);
// TPC frame points
const TVector3 D4_UL(0.2009, -0.8964, 10.3024);
const TVector3 D4_DL(0.0787, -0.8963, 11.4982);
const TVector3 D4_UR(-1.8011, -0.8797, 9.9464);
const TVector3 D4_DR(-1.9547, -0.8827, 11.4366);
const TVector3 D4_LC(0.0818, -0.0504, 10.8961);
const TVector3 D4_RC(-1.4230, -0.0557, 10.7402);
// Calculated TPC active region points
const TVector3 TPCActiveBUL(-0.2778, -0.5614, 10.3048);
const TVector3 TPCActiveBDL(-0.3896, -0.5614, 11.3991);
const TVector3 TPCActiveTUL(-0.2778, 0.5386, 10.3048);
const TVector3 TPCActiveTDL(-0.3896, 0.5386, 11.3991);
const TVector3 TPCActiveBUR(-0.7290, -0.5614, 10.2587);
const TVector3 TPCActiveBDR(-0.8408, -0.5614, 11.3530);
const TVector3 TPCActiveTUR(-0.7290, 0.5386, 10.2587);
const TVector3 TPCActiveTDR(-0.8408, 0.5386, 11.3530);
const TVector3 TPCActiveC(-0.5593, -0.0114, 10.82885);

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
std::vector<const char*> str4BlockVec = {/*str4Block0, */str4Block1, str4Block2, str4Block3};
const char* dstofDir = "/nfs/scratch0/dbrailsf/data_backup/dtof_backup/";
const char* ustofDir = "/nfs/scratch0/dbrailsf/data_backup/utof_backup_firsthitpinnedtounixtime/Data_root_v3_wo_walk_corr/";
// Time for signal to travel down and S4
const double s4BarTime = 18.2;
// Maximum rate of cosmics measured for some arbitrary area -- used for normalisation
const double maxCosmics = 1.627; // Hz

// S3 amplitude cut for protons
// Apply to A1ToF and A2ToF
// AK's standard cut
const double ACut = 0.25;
// SJ's bar by bar amplitude cut (by eye)
// For A1 
const std::vector<double> A1CutVec = {0.25, 0.25, 0.25, 0.2, 0.25, 
				      0.25, 0.225, 0.25, 0.275, 0.275, 
				      0.3, 0.2, 0.2, 0.225, 0.25, 
				      0.225, 0.225, 0.3, 0.25, 0.225,
				      0.3, 0.25};
// For A2
const std::vector<double> A2CutVec = {0.3, 0.275, 0.25, 0.1, 0.25, 
				      0.225, 0.14, 0.225, 0.225, 0.225, 
				      0.21, 0.25, 0.225, 0.19, 0.19, 
				      0.225, 0.21, 0.225, 0.25, 0.21, 
				      0.275, 0.25};

// S4 fitting regions and timing cuts
// Timing cuts
const double piLowS4  = 36.;
const double piHiS4   = 51.;
// Cuts for selecting the protons
const double proCutLowS4 = 62.;
const double proCutHiS4  = 125.;
// Fit regions for protons
const vector<double> proFitLowS4 = {62., 62., 69., 75., 75.};
const vector<double> proFitHiS4  = {86., 94., 100., 105., 285.};
// Deadtime cut for S4
// Enforces a deadtime on S4 bars, necessary to avoid some weird effects
const double s4DeadtimeCut = 200.; // ns
const double barOverallEff    = .8;
const double barOverallEffErr = .1;

// Beam paper binnings and limits
const int nBinsS4Horz = 20;
const int nBinsS4HorzCosmics = 2;
const double binsCosmics[] = {0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 
			      77, 84, 91, 98, 105, 112, 119, 126, 133, 140};
const int binnumCosmics = sizeof(binsCosmics)/sizeof(double) - 1;
const double binsS4HorzLow  = -6.;
const double binsS4HorzHigh = 0.;
const int nBinsS4Vert = 9;
const double binsS4VertLow  = -1.5;
const double binsS4VertHigh = 1.42;
const int nBinsDtof = 182;
const double binsDtofLow  = 30.;
const double binsDtofHigh = proCutHiS4;
const double dtofFixedWidth = 1.4;

// For S3 paper
// double binsTheta[] = {-3.1, -3.,
// 		      -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, -2., 
// 		      -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.,
// 		      -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 
// 		      0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
// 		      1.0, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 
// 		      2.0, 2.2, 2.4, 2.6, 2.8, 
// 		      3.0, 3.2, 3.4, 3.6, 3.8,  
// 		      4.0, 4.25, 4.5, 4.75, 
// 		      5.0, 5.25, 5.8,
// 		      6.}; 
double binsTheta[] = {-6., -5.8, -5.25, -5.,
		      -4.75, -4.5, -4.25, -4.0, 
		      -3.8, -3.6, -3.4, -3.2, -3.0, 
		      -2.8, -2.6, -2.4, -2.2, -2.0, 
		      -1.875, -1.75, -1.625, -1.5, -1.375, -1.25, -1.125, -1.0, 
		      -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.,
		      0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 
		      1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 
		      2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3., 
		      3.1};
int binnum = sizeof(binsTheta)/sizeof(double) - 1;

const vector<double> s4s3MC = {0.0281, 0.0680, 0.0861, 0.0582,  0.0149};

// Deadtime corrections
// These are the ones used for the S1S2 hits
// Just use a constant ratio for the 0 block case
const double block0Slope    = 0.;/*-0.00037602;*/
const double block0SlopeErr = 0;/*5.04972e-5;*/
const double block0Const    = 0.0913247;/*0.2394;*/
const double block0ConstErr = 0.00102858;/*0.01989655;*/
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
const double block4Slope0    = -0.00054221873;
const double block4Slope0Err = 1.371236-05;
const double block4Const0    = 1.0005352;
const double block4Const0Err = 0.019797654;
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
std::vector<double> block4SlopeVec    = {/*block4Slope0, */block4Slope1, block4Slope2, block4Slope3};
std::vector<double> block4SlopeErrVec = {/*block4Slope0Err, */block4Slope1Err, block4Slope2Err, block4Slope3Err};
std::vector<double> block4ConstVec    = {/*block4Const0, */block4Const1, block4Const2, block4Const3};
std::vector<double> block4ConstErrVec = {/*block4Const0Err, */block4Const1Err, block4Const2Err, block4Const3Err};

// Deadtime corrections for S1 hits
// 0 block
const double block0SlopeS1     = 7.783e-5;
const double block0SlopeErrS1  = 1.419e-5;
const double block0ConstS1     = 0.0001268;
const double block0ConstErrS1  = 0.005448;
// 1 block
const double block1SlopeS1     = 5.949e-6;
const double block1SlopeErrS1  = 5.013e-6;
const double block1ConstS1     = 0.08174;
const double block1ConstErrS1  = 0.005121;
// 2 block
const double block2SlopeS1     = -1.887e-6;
const double block2SlopeErrS1  = 3.69e-6;
const double block2ConstS1     = 0.1529;
const double block2ConstErrS1  = 0.00577;
// 3 block
const double block3SlopeS1     = 6.578e-6;
const double block3SlopeErrS1  = 7.04e-6;
const double block3ConstS1     = 0.185;
const double block3ConstErrS1  = 0.01287;
// 4 block
const double block4SlopeS10    = -6.4762466e-6;
const double block4SlopeErrS10 = 6.1196825e-6;
const double block4ConstS10    = 0.22989897;
const double block4ConstErrS10 = 0.0087381546;
const double block4SlopeS11    = 1.195e-5;
const double block4SlopeErrS11 = 6.553e-6;
const double block4ConstS11    = 0.2036;
const double block4ConstErrS11 = 0.009115;
const double block4SlopeS12    = 1.31e-5; 
const double block4SlopeErrS12 = 3.117e-6;
const double block4ConstS12    = 0.2061;
const double block4ConstErrS12 = 0.004548;
const double block4SlopeS13    = 1.363e-5;
const double block4SlopeErrS13 = 5.827e-6;
const double block4ConstS13    = 0.2053;
const double block4ConstErrS13 = 0.008546;
std::vector<double> block4SlopeS1Vec    = {/*block4SlopeS10, */block4SlopeS11, 
					   block4SlopeS12, block4SlopeS13};
std::vector<double> block4SlopeErrS1Vec = {/*block4SlopeErrS10, */block4SlopeErrS11, 
					   block4SlopeErrS12, block4SlopeErrS13};
std::vector<double> block4ConstS1Vec    = {/*block4ConstS10, */block4ConstS11, 
					   block4ConstS12, block4ConstS13};
std::vector<double> block4ConstErrS1Vec = {/*block4ConstErrS10, */block4ConstErrS11, 
					   block4ConstErrS12, block4ConstErrS13};

/// Functions
vector<double> getDtofEdges() 
{
  std::vector<double> vec;
  for (int i=0; i<61; i++) {
    double edge = i * (dtofFixedWidth) + binsDtofLow;
    vec.push_back(edge);
  }
  for (int i=62; i<95; i++) {
    if (i % 2 == 0) {
      double edge = i * (dtofFixedWidth) + binsDtofLow;
      vec.push_back(edge);
    }
  }
  for (int i=95; i<183; i++) {
    if ((i+1) % 4==0) {
      double edge = i * (dtofFixedWidth) + binsDtofLow;
      vec.push_back(edge);
    }
  }
  return vec;
}

double dtofEdges[] = {30.0, 31.4, 32.8, 34.2, 35.6, 37.0, 38.4, 39.8, 41.2, 42.6,
		      44.0, 45.4, 46.8, 48.2, 49.6, 51.0, 52.4, 53.8, 55.2, 56.6,
		      58.0, 59.4, 60.8, 62.2, 63.6, 65.0, 66.4, 67.8, 69.2, 70.6,
		      72.0, 73.4, 74.8, 76.2, 77.6, 79.0, 80.4, 81.8, 83.2, 84.6,
		      86.0, 87.4, 88.8, 90.2, 91.6, 93.0, 94.4, 95.8, 97.2, 98.6,
		      100.0, 102.8, 105.6, 108.4, 111.2, 
		      114.0, 119.6, 125.2};//, 130.8, 136.4, 142, 147.6, 153.2, 158.8, 164.4, 
		      // 170.0, 175.6, 181.2, 186.8, 192.4, 198, 203.6, 209.2, 214.8, 220.4,
		      // 226.0, 231.6, 237.2, 242.8, 248.4, 254, 259.6, 265.2, 270.8, 276.4,
		      // 282, 285};
int binsDtofVar = sizeof(dtofEdges)/sizeof(double) - 1;

// Histogram styles
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
  h2->GetXaxis()->SetTitleSize(.05);
  h2->GetYaxis()->SetTitleSize(.05);
  h2->GetZaxis()->SetTitleSize(.05);
  h2->GetXaxis()->SetLabelSize(.05);
  h2->GetYaxis()->SetLabelSize(.05);
  h2->GetZaxis()->SetLabelSize(.05);
  h2->Sumw2();
  h2->SetOption("colz");
}

void setHistAttr(TGraph *g, const int col)
{
  g->SetLineColor(col);
  g->SetMarkerColor(col);
  g->SetMarkerStyle(8);
  g->SetLineWidth(2);
}

// Downstream TOF global coordinates
TVector3 GetDtofGlobalCoords(const double fakeTime0, const double fakeTime1, const int bar) 
{
  double posX = (((fakeTime1-fakeTime0)*(70./s4BarTime)+65.)/130.) * (S4_TR.X() - S4_TL.X()) + S4_TL.X();
  double posY = (S4BarsL.at(bar-1).Y() + S4BarsR.at(bar-1).Y())/2. - 0.05;
  double posZ = 0.;
  if (bar % 2 == 0) {
    posZ = (((fakeTime1-fakeTime0)*(70./s4BarTime)+65.)/130.) * (12.12 - 12.16) + 12.16;
  }
  else {
    posZ = (((fakeTime1-fakeTime0)*(70./s4BarTime)+65.)/130.) * (12.34 - 12.4) + 12.4;
  }
  TVector3 vec(posX, posY, posZ);
  return vec;
}

// Upstream TOF global coordinates
TVector3 GetUtofGlobalCoords(const double xtof, const double ytof) {
  double posX = (xtof / 168.) * (S3_TR.X() - S3_TL.X()) + S3_TL.X();
  double posY = (ytof / 100.) + (S3_BR.Y() + S3_BL.Y())/2.;
  double posZ = (xtof / 168.) * (S3_TR.Z() - S3_TL.Z()) + S3_TL.Z();
  TVector3 coords(posX, posY, posZ);  
  return coords;
}

// Returns the local position along the bar of a hit
double localDtofPosition(const double fakeTime0, const double fakeTime1) 
{
  double pos = ((fakeTime1-fakeTime0)*(70./s4BarTime)+70.);
  return pos;
}

// Angles for beam flux paper from survey coordinates
double getThetaFromGlobal(const TVector3 vec)
{
  const TVector3 vS1 = (D1_ULB+D1_ULT+D1_URB+D1_URT)*0.25;
  TVector3 v = vec;
  v = v - vS1;
  double theta = TMath::ATan(v.X()/v.Z())*(180./TMath::Pi());
  return theta;
}

double getPhiFromGlobal(const TVector3 vec)
{
  const TVector3 vS1 = (D1_ULB+D1_ULT+D1_URB+D1_URT)*0.25;
  TVector3 v = vec;
  v = v - vS1;
  double phi = TMath::ATan(v.Y()/v.Z())*(180./TMath::Pi());
  return phi;
}

// Converts global coordinates to the coordinate system in the MC
TVector3 globalToMCCoords(TVector3 v)
{
  TVector3 vec;
  vec.SetX(v.X() + 0.6706);
  vec.SetY(v.Y() + 0.05305);
  vec.SetZ(v.Z() - 10.8182);
  return vec;
}

// Calculating actual time of hit in S4
double dtofHitTime(const double t0, const double t1) {
  double t = min(t0, t1) - (s4BarTime/2. - TMath::Abs(t0 - t1) / 2.);
  return t;
}

// Converts MC coordinates to the global coordinate system from the survey
TVector3 MCToGlobalCoords(TVector3 v)
{
  TVector3 vec;
  vec.SetX(v.X() - 0.6706);
  vec.SetY(v.Y() - 0.05305);
  vec.SetZ(v.Z() + 10.8182);
  return vec;
}

// Get bar number from MC coordinates
int dtofBarFromMC(const TVector3 vec) 
{
  int bar=0;
  TVector3 global = MCToGlobalCoords(vec);

  for (int b=0; b<10; b++) {
    if ( global.Y() > S4BarsL.at(b).Y() - 0.1 && global.Y() < S4BarsL.at(b).Y() && b==0) {
      bar = 1;
      break;
    }
    else if ( global.Y() > S4BarsL.at(b).Y() - 0.075 && global.Y() < S4BarsL.at(b).Y() && b!=0) {
      bar = b+1;
      break;
    }
    else if ( global.Y() > S4BarsL.at(9).Y() ) {
      bar = 11;
      break;
    }
  }
  return bar;
}

// Get bar position from MC coordinates
double localDtofFromMC(const TVector3 vec)
{
  TVector3 global = MCToGlobalCoords(vec);
  double t1Minust0 = ((global.X() - S4_TL.X())*130./(S4_TR.X()-S4_TL.X()) - 65.);
  double pos = t1Minust0 + 70.;
  return pos;
}

// Mass squared from the time
// momentum and mass in GeV
// Time is in ns
double massFromTime(const double time, const double mom, const double base) {
  double mass = pow(mom,2) * (pow((time*1e-9)/base, 2)*pow(3e8, 2) - 1);
  return mass;
}

// Time in ns
// Outputs momentum in MeV/c
double momFromTime(const double mass, const double baseline, const double time)
{
  double mom = mass * 1000. * (baseline/(time*1e-9))*(1/TMath::Sqrt((pow(3e8*1e-9*time,2) - pow(baseline,2))/pow(time*1e-9,2)));
  return mom;
}
// The same but take a start and end point
double momFromTime(const double mass, const TVector3 start, const TVector3 end, const double time)
{
  const TVector3 dist = end - start;
  const double baseline = dist.Mag();
  double mom = momFromTime(mass, baseline, time);
  return mom;
}

// Time in ns
// Output K.E. in MeV
double keFromTime(const double mass, const double baseline, const double time)
{
  double mom = momFromTime(mass, baseline, time);
  double ke = 1000. * (TMath::Sqrt( pow(mom/1000., 2) + pow(mass, 2) ) - mass);
  return ke;
}
// The same but take a start and end point
double keFromTime(const double mass, const TVector3 start, const TVector3 end, const double time)
{
  double mom = momFromTime(mass, start, end, time);
  double ke = 1000. * (TMath::Sqrt( pow(mom/1000., 2) + pow(mass, 2) ) - mass);
  return ke;
}

// Get line colour from block colour
int getColourFromBlock(const int block)
{
  int col = 0;
  if (block==0) col=1;
  else if (block==1) col=632;
  else if (block==2) col=600;
  else if (block==3) col=433;
  else if (block==4) col=801;

  return col;
}

void makeEffHist(TH2D* h2) 
{
  h2->Scale(1., "width");
  h2->Scale(1. / h2->GetBinContent(h2->GetMaximumBin()));
}

const double tpcThetaLow  = getThetaFromGlobal(TPCActiveBDR);
const double tpcThetaHigh = getThetaFromGlobal(TPCActiveBUL);
const double tpcPhiHigh  = getPhiFromGlobal(TPCActiveTUL);
const double tpcPhiLow = getPhiFromGlobal(TPCActiveBUL);

