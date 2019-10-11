// toyS4Cosmics.C
// Macro to estimate the number of coincidences we should see in each bar

void toyS4Cosmics(const char* saveDir, const int nAttempts = 10000000)
{
  const double cosmicStartY = 10.;
  const double cosmicLowX = -12.;
  const double cosmicHiX  = 12.;
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

  const double yOffset = 0.3530;
  double avgX = 0.;
  for (int i=0; i<S4BarsL.size(); i++) {
    avgX += S4BarsL.at(i).Z();
  }
  avgX /= S4BarsL.size();
  // Because we are doing this in 2D for now, we can imagine that we are in a plane 
  // which is directly end on to the bars
  std::vector<double> s4BarTops;
  std::vector<double> s4BarBottoms;
  std::vector<double> s4BarXs;
  for (int i=0; i<S4BarsL.size(); i++) {
    s4BarTops.push_back(S4BarsL.at(i).Y() + yOffset);
    s4BarBottoms.push_back(S4BarsL.at(i).Y() - 0.1 + yOffset);
    s4BarXs.push_back(S4BarsL.at(i).Z() - avgX);
  }
  TRandom3 *rand = new TRandom3(6789);

  TFile *fout = new TFile(saveDir, "recreate");
  // 2D approximation - x is such that the bars are 1d in x, y is vertical
  // Cosmics are considered to be distributed uniformly in x and all start at same height
  TH1D *hCosDist = new TH1D("hCosDist", "Cosmic Distribution; #theta / #pi", 100., -0.5, 0.5);
  TH1D *hCosXDist = new TH1D("hCosXDist", "Cosmic Distribution; X / m", 100., cosmicLowX, cosmicHiX);
  TH2D *h2CosXAngDist = new TH2D("h2CosXAngDist", "Cosmic Distribution; X / m", 100., cosmicLowX, cosmicHiX, 100., -0.5, 0.5);
  TH1D *hCosEndX = new TH1D("hCosEndX", "Cosmic end point; X / m; Events", 100., cosmicLowX, cosmicHiX);
  for (int n = 0; n < nAttempts; n++) {
    double flatpi = -TMath::Pi()/2. + rand->Rndm() * TMath::Pi();
    double flat = rand->Rndm();
    if (flat < (TMath::Cos(flatpi) * TMath::Cos(flatpi))) {
      hCosDist->Fill(flatpi/TMath::Pi());
      double xStart = cosmicLowX + (cosmicHiX - cosmicLowX) * rand->Rndm();
      hCosXDist->Fill(xStart);
      h2CosXAngDist->Fill(xStart, flatpi/TMath::Pi());
      hCosEndX->Fill(xStart + cosmicStartY * TMath::Tan(flatpi));
    }
  }
  //  hCosDist->Draw("hist");
  //hCosXDist->Draw("hist");

  fout->cd();
  h2CosXAngDist->Write();
  hCosXDist->Write();
  hCosDist->Write();
  hCosEndX->Write();
  fout->Close();
  delete fout;

} // toyS4Cosmics
