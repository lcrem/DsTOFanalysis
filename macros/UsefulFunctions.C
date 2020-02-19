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
std::vector<const char*> str4BlockVec = {str4Block1, str4Block2, str4Block3};
// Time for signal to travel down and S4
const double s4BarTime = 18.2;
// Maximum rate of cosmics measured for some arbitrary area -- used for normalisation
const double maxCosmics = 1.627; // Hz

// S4 fitting regions and timing cuts
// Timing cuts
const double piLowS4  = 36.;
const double piHiS4   = 51.;
// Cuts for selecting the protons
const double proCutLowS4 = 62.;
const double proCutHiS4  = 285.;
// Fit regions for protons
const vector<double> proFitLowS4 = {62., 62., 69., 75., 75.};
const vector<double> proFitHiS4  = {86., 94., 100., 105., 285.};
// Deadtime cut for S4
// Enforces a deadtime on S4 bars, necessary to avoid some weird effects
const double s4DeadtimeCut = 200.; // ns

/// Functions
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
  h2->SetOption("colz");
  h2->GetXaxis()->SetTitleSize(.05);
  h2->GetYaxis()->SetTitleSize(.05);
  h2->GetZaxis()->SetTitleSize(.05);
  h2->GetXaxis()->SetLabelSize(.05);
  h2->GetYaxis()->SetLabelSize(.05);
  h2->GetZaxis()->SetLabelSize(.05);
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

// Returns the local position along the bar of a hit
double localDtofPosition(const double fakeTime0, const double fakeTime1) 
{
  double pos = ((fakeTime1-fakeTime0)*(70./s4BarTime)+70.);
  return pos;
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

// Mass squared from the time
// momentum and mass in GeV
// Time is in ns
double massFromTime(const double time, const double mom, const double base) {
  double mass = pow(mom,2) * (pow((time*1e-9)/base, 2)*pow(3e8, 2) - 1);
  return mass;
}

// Time in ns
// Outputs momentum in GeV/c
double momFromTime(const double mass, const double baseline, const double time)
{
  double mom = 0.;
  mom = mass * 1000. * (baseline/(time*1e-9))*(1/TMath::Sqrt((pow(3e8*1e-9*time,2) - pow(baseline,2))/pow(time*1e-9,2)));
  return mom;
}
// Time in ns
// Output K.E. in GeV
double keFromTime(const double mass, const double baseline, const double time)
{
  double mom = momFromTime(mass, baseline, time);
  double ke = TMath::Sqrt( pow(mom, 2) + pow(mass, 2) ) - mass;
  return ke;
}
