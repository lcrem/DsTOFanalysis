// toyS4Cosmics.C
// Macro to estimate the number of coincidences we should see in each bar

const double startRadius  = 7.;

TVector3 intersectionPnt(const TVector3 startPnt, const TVector3 dir, 
			 const TVector3 barVec, const TVector3 barPnt)
{
  double d = ((barPnt - startPnt).Dot(barVec))/dir.Dot(barVec);
  TVector3 intersection = d * dir + startPnt;
  return intersection;
}

bool crossesBar(const TVector3 startPnt, const TVector3 dir, const TVector3 barVec,  
		const TVector3 bar1, const TVector3 bar2, const TVector3 bar3)
{
  TVector3 intPnt = intersectionPnt(startPnt, dir, barVec, bar1);
  bool cross = false;
  if ((bar2 * (bar1-bar2)) <= (intPnt * (bar1-bar2)) && 
      (bar1 * (bar1-bar2)) >= (intPnt * (bar1-bar2)) &&
      (bar2 * (bar3-bar2)) <= (intPnt * (bar3-bar2)) && 
      (bar3 * (bar3-bar2)) >= (intPnt * (bar3-bar2)))
    cross = true;
  return cross;
}

void toyS4Cosmics(const char* saveDir, const int nAttempts = 1e8)
{
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
  const std::vector<TVector3> s4BarsL = {D5_D5L, D5_U5L, D5_D4L, D5_U4L, D5_D3L, D5_U3L, 
					 D5_D2L, D5_U2L, D5_D1L, D5_U1L};
  const std::vector<TVector3> s4BarsR = {D5_D5R, D5_U5R, D5_D4R, D5_U4R, D5_D3R, D5_U3R, 
					 D5_D2R, D5_U2R, D5_D1R, D5_U1R};

  const double yOffset = 0.3530;
  double avgX = 0.;
  double avgZ = 0.;
  for (int i=0; i<s4BarsL.size(); i++) {
    avgX += s4BarsL.at(i).X();
    avgX += s4BarsR.at(i).X();
    avgZ += s4BarsL.at(i).Z();
    avgZ += s4BarsR.at(i).Z();
  }
  avgX /= (s4BarsL.size()+s4BarsR.size());
  avgZ /= (s4BarsL.size()+s4BarsR.size());
  cout<<"Avg x, z "<<avgX<<", "<<avgZ<<endl;
  vector<double> s4BarXs;
  vector<double> s4BarTops;
  vector<double> s4BarBottoms;
  TVector3 shiftVector(avgX, 0., avgZ);

  vector<TVector3> s4LU1, s4RU1, s4LD1, s4RD1, s4LU2, s4RU2, s4LD2, s4RD2;
  vector<TVector3> s4BarVec1, s4BarVec2;
  for (int i=0; i<s4BarsL.size(); i++) {
    s4BarTops.push_back(s4BarsL.at(i).Y() + yOffset);
    s4BarBottoms.push_back(s4BarsL.at(i).Y() - 0.1 + yOffset);
    s4BarXs.push_back(s4BarsL.at(i).Z() - avgZ);
    TVector3 LU, RU, LD, RD; 
    LU = s4BarsL.at(i);
    LU.SetY(LU.Y()+yOffset);
    RU = s4BarsR.at(i);
    RU.SetY(RU.Y()+yOffset);
    LD = s4BarsL.at(i);
    LD.SetY(LD.Y()+yOffset-0.1);
    RD = s4BarsR.at(i);
    RD.SetY(RD.Y()+yOffset-0.1);
    s4LU1.push_back(LU);
    s4RU1.push_back(RU);
    s4LD1.push_back(LD);
    s4RD1.push_back(RD);
    s4BarVec1.push_back(((LU - RU).Cross(RD - RU)).Unit());
    LU.SetZ(LU.Z()-0.01);
    RU.SetZ(RU.Z()-0.01);
    LD.SetZ(LD.Z()-0.01);
    RD.SetZ(RD.Z()-0.01);
    s4LU2.push_back(LU);
    s4RU2.push_back(RU);
    s4LD2.push_back(LD);
    s4RD2.push_back(RD);
    s4BarVec2.push_back(((LU - RU).Cross(RD - RU)).Unit());
  }
  vector< vector<TVector3> > s4Pnts1, s4Pnts2;
  s4Pnts1.push_back(s4LU1);
  s4Pnts1.push_back(s4RU1);
  s4Pnts1.push_back(s4LD1);
  s4Pnts1.push_back(s4RD1);
  s4Pnts2.push_back(s4LU2);
  s4Pnts2.push_back(s4RU2);
  s4Pnts2.push_back(s4LD2);
  s4Pnts2.push_back(s4RD2);

  TRandom3 *rand = new TRandom3(6789);

  TFile *fout = new TFile(saveDir, "recreate");
  // 2D approximation - x is such that the bars are 1d in x, y is vertical
  // Cosmics are considered to be distributed uniformly in spatial phi and theta
  // However, need cos squared theta for direction
  TH1D *hCosDist = new TH1D("hCosDist", "Cosmic Distribution; #theta / #pi", 100., -0.5, 0.5);
  hCosDist->Sumw2();
  hCosDist->GetXaxis()->SetTitleSize(.05);
  hCosDist->GetYaxis()->SetTitleSize(.05);
  hCosDist->GetXaxis()->SetLabelSize(.05);
  hCosDist->GetYaxis()->SetLabelSize(.05);
  TH1D *hCosPhiDist = new TH1D("hCosPhiDist", "Cosmic Distribution; #phi / #pi", 100., 0, 2.);
  hCosPhiDist->Sumw2();
  hCosPhiDist->GetXaxis()->SetTitleSize(.05);
  hCosPhiDist->GetYaxis()->SetTitleSize(.05);
  hCosPhiDist->GetXaxis()->SetLabelSize(.05);
  hCosPhiDist->GetYaxis()->SetLabelSize(.05);
  TH2D *h2CosPhiThetaDist = new TH2D("h2CosPhiThetaDist", "Cosmic Distribution; #phi / #pi; #theta / #pi", 100., 0., 2., 100., -0.5, 0.5);
  h2CosPhiThetaDist->Sumw2();
  h2CosPhiThetaDist->GetXaxis()->SetTitleSize(.05);
  h2CosPhiThetaDist->GetYaxis()->SetTitleSize(.05);
  h2CosPhiThetaDist->GetXaxis()->SetLabelSize(.05);
  h2CosPhiThetaDist->GetYaxis()->SetLabelSize(.05);
  TH2D *h2PhiThetaSpatialDist = new TH2D("h2PhiThetaSpatialDist", "Cosmic Distribution; #phi / #pi; #theta / #pi", 100., -1., 0., 100., -1., 0.);
  h2PhiThetaSpatialDist->Sumw2();
  h2PhiThetaSpatialDist->GetXaxis()->SetTitleSize(.05);
  h2PhiThetaSpatialDist->GetYaxis()->SetTitleSize(.05);
  h2PhiThetaSpatialDist->GetXaxis()->SetLabelSize(.05);
  h2PhiThetaSpatialDist->GetYaxis()->SetLabelSize(.05);
  TH1D *hBarHits = new TH1D("hBarHits", "Cosmic hits in each bar; Bar; Hits", 10, 0.5, 10.5);
  hBarHits->Sumw2();
  hBarHits->GetXaxis()->SetTitleSize(.05);
  hBarHits->GetYaxis()->SetTitleSize(.05);
  hBarHits->GetXaxis()->SetLabelSize(.05);
  hBarHits->GetYaxis()->SetLabelSize(.05);
  TH2D *h2BarCoins = new TH2D("h2BarCoins", "Bar coincidences; First bar; Second bar; Fraction", 10, 0.5, 10.5, 10., 0.5, 10.5);
  h2BarCoins->Sumw2();
  h2BarCoins->GetXaxis()->SetTitleSize(.05);
  h2BarCoins->GetYaxis()->SetTitleSize(.05);
  h2BarCoins->GetXaxis()->SetLabelSize(.05);
  h2BarCoins->GetYaxis()->SetLabelSize(.05);
  h2BarCoins->GetZaxis()->SetTitleSize(.05);
  h2BarCoins->GetZaxis()->SetLabelSize(.05);
  TH2D *h2EndPos = new TH2D("h2EndPos", "Cosmics end positions; X / m; Z / m; Number", 100, -8.+shiftVector.X(), 8.+shiftVector.X(), 100, -8.+shiftVector.Z(), 8.+shiftVector.Z());
  h2EndPos->Sumw2();
  h2EndPos->GetXaxis()->SetTitleSize(.05);
  h2EndPos->GetYaxis()->SetTitleSize(.05);
  h2EndPos->GetXaxis()->SetLabelSize(.05);
  h2EndPos->GetYaxis()->SetLabelSize(.05);
  h2EndPos->GetZaxis()->SetTitleSize(.05);
  h2EndPos->GetZaxis()->SetLabelSize(.05);
  TH3D *h3StartPos = new TH3D("h3StartPos", "Cosmic ray start positions; X / m; Y / m; Z / m", 50, -6.+shiftVector.X(), 6.+shiftVector.X(), 50, -6., 6., 50, -6.+shiftVector.Z(), 6.+shiftVector.Z());
  h3StartPos->Sumw2();
  h3StartPos->GetXaxis()->SetTitleSize(.05);
  h3StartPos->GetYaxis()->SetTitleSize(.05);
  h3StartPos->GetXaxis()->SetLabelSize(.05);
  h3StartPos->GetYaxis()->SetLabelSize(.05);
  h3StartPos->GetZaxis()->SetTitleSize(.05);
  h3StartPos->GetZaxis()->SetLabelSize(.05);
  TH3D *h3StartDir = new TH3D("h3StartDir", "Cosmic ray start directions; X / m; Y / m; Z / m", 50, -1, 1, 50, -1, 1, 50, -1, 1);
  h3StartDir->Sumw2();
  h3StartDir->GetXaxis()->SetTitleSize(.05);
  h3StartDir->GetYaxis()->SetTitleSize(.05);
  h3StartDir->GetXaxis()->SetLabelSize(.05);
  h3StartDir->GetYaxis()->SetLabelSize(.05);
  h3StartDir->GetZaxis()->SetTitleSize(.05);
  h3StartDir->GetZaxis()->SetLabelSize(.05);
  TH2D *h2StartXY = new TH2D("h2StartXY", "Cosmic ray start positions; X / m; Y / m", 50, -6.+shiftVector.X(), 6.+shiftVector.X(), 50, -6., 6.);
  h2StartXY->Sumw2();
  h2StartXY->GetXaxis()->SetTitleSize(.05);
  h2StartXY->GetYaxis()->SetTitleSize(.05);
  h2StartXY->GetXaxis()->SetLabelSize(.05);
  h2StartXY->GetYaxis()->SetLabelSize(.05);
  h2StartXY->GetZaxis()->SetTitleSize(.05);
  h2StartXY->GetZaxis()->SetLabelSize(.05);
  TH2D *h2StartXZ = new TH2D("h2StartXZ", "Cosmic ray start positions; X / m; Z / m", 50, -6+shiftVector.X(), 6.+shiftVector.X(), 25, -6.+shiftVector.Z(), 6+shiftVector.Z());
  h2StartXZ->Sumw2();
  h2StartXZ->GetXaxis()->SetTitleSize(.05);
  h2StartXZ->GetYaxis()->SetTitleSize(.05);
  h2StartXZ->GetXaxis()->SetLabelSize(.05);
  h2StartXZ->GetYaxis()->SetLabelSize(.05);
  h2StartXZ->GetZaxis()->SetTitleSize(.05);
  h2StartXZ->GetZaxis()->SetLabelSize(.05);

  vector<TH3D*> h3BarVec;
  for (int b=0; b<10; b++) {
    TH3D *h3Bar = new TH3D(Form("h3Bar%d", b+1), Form("Bar %d; X / m; Z / m; Y / m", b+1), 50, -1.+shiftVector.X(), 1.+shiftVector.X(), 50, -2, 2., 50, -.5+shiftVector.Z(), .5+shiftVector.Z());
    h3BarVec.push_back(h3Bar);
  }
  TH3D *h3BarTot = new TH3D("h3BarTot", "Total; X / m; Z / m; Y / m", 50, -1.+shiftVector.X(), 1.+shiftVector.X(), 50, -2, 2., 50, -.5+shiftVector.Z(), .5+shiftVector.Z());

  TVector3 floorVec(0., 1., 0.);
  for (int n = 0; n < nAttempts; n++) {
    double dirTheta = -TMath::Pi()/2. + rand->Rndm() * TMath::Pi();
    double flat = rand->Rndm();
    if (flat < (TMath::Cos(dirTheta) * TMath::Cos(dirTheta))) {
      hCosDist->Fill(dirTheta/TMath::Pi());
      double dirPhi = TMath::Pi() * rand->Rndm();
      // Generate uniform points on sphere
      double x = -startRadius + 2. * startRadius * rand->Rndm();
      double z = -startRadius + 2. * startRadius * rand->Rndm();
      double r = TMath::Sqrt(pow(x, 2)+pow(z, 2));
      if (r < startRadius) {
	double y = TMath::Sqrt(pow(startRadius, 2) - pow(x, 2) - pow(z, 2));
	TVector3 v(x, y, z);
	hCosPhiDist->Fill(dirPhi / TMath::Pi());
	h2PhiThetaSpatialDist->Fill(v.Phi()/TMath::Pi(), v.Theta()/TMath::Pi());
	h2CosPhiThetaDist->Fill(dirPhi/TMath::Pi(), dirTheta/TMath::Pi());
	v += shiftVector;
	cout<<"x, y, z "<<x<<", "<<y<<", "<<z<<endl;
	h3StartPos->Fill(v.X(), v.Y(), v.Z());
	h2StartXY->Fill(v.X(), v.Y());
	h2StartXZ->Fill(v.X(), v.Z());
	TVector3 dir;
	dir.SetMagThetaPhi(1., dirPhi, dirTheta);
	h3StartDir->Fill(-dir.X(), -dir.Y(), -dir.Z());
	TVector3 f(0., 0., 0.);
	TVector3 floorInt = intersectionPnt(v, -dir, floorVec, f);
	h2EndPos->Fill(floorInt.X(), floorInt.Z());
	// Now test to see which bar (if any) it passes through
	for (int b=0; b < s4Pnts1[0].size(); b++) {	
	  TVector3 int11 = intersectionPnt(v, -dir, s4BarVec1.at(b), s4Pnts1.at(0).at(b));
	  TVector3 int12 = intersectionPnt(v, -dir, s4BarVec2.at(b), s4Pnts2.at(0).at(b));
	  if ((crossesBar(v, -dir, s4BarVec1.at(b), s4Pnts1.at(1).at(b), 
			  s4Pnts1.at(0).at(b), s4Pnts1.at(2).at(b)) && int11.Y() <= v.Y()) ||
	      (crossesBar(v, -dir, s4BarVec2.at(b), s4Pnts2.at(1).at(b), 
			  s4Pnts2.at(0).at(b), s4Pnts2.at(2).at(b)) && int12.Y() <= v.Y())) {
	    hBarHits->Fill(b+1);
	    for (int b2=b; b2 < s4Pnts1[0].size(); b2++) {
	      TVector3 int21 = intersectionPnt(v, -dir, s4BarVec1.at(b2), s4Pnts1.at(0).at(b2));
	      TVector3 int22 = intersectionPnt(v, -dir, s4BarVec2.at(b2), s4Pnts2.at(0).at(b2));
	      if (((crossesBar(v, -dir, s4BarVec1.at(b2), s4Pnts1.at(1).at(b2), 
			       s4Pnts1.at(0).at(b2), s4Pnts1.at(2).at(b2)) && int21.Y() <= v.Y()) ||
		   (crossesBar(v, -dir, s4BarVec2.at(b2), s4Pnts2.at(1).at(b2), 
			       s4Pnts2.at(0).at(b2), s4Pnts2.at(2).at(b2)) && int22.Y() <= v.Y())) && b!=b2) {
		if (int11.Y() > int21.Y()) h2BarCoins->Fill(b+1, b2+1);
		else h2BarCoins->Fill(b2+1, b+1);
		h3BarTot->Fill(int11.X(), int11.Z(), int11.Y());
		h3BarTot->Fill(int21.X(), int21.Z(), int21.Y());
		h3BarVec.at(b)->Fill(int11.X(), int11.Z(), int11.Y());
		h3BarVec.at(b2)->Fill(int21.X(), int21.Z(), int21.Y());
	      }
	    }
	  }
	} // Loop over bars
      }
    }
  }
  // h2BarCoins->Scale(1. / h2BarCoins->Integral());

  // Some graphs to help show bar positions
  TMultiGraph *mg = new TMultiGraph("mgBars", "S4 bar positions in MC; X / m; Y / m");
  for (int i=0; i<s4BarTops.size(); i++) {
    TGraph *grBar = new TGraph();
    grBar->SetTitle(Form("Bar %d", i+1));
    grBar->SetLineWidth(2);
    grBar->SetLineColor(52 + i * 3);
    grBar->SetPoint(grBar->GetN(), s4BarXs.at(i), s4BarTops.at(i));
    grBar->SetPoint(grBar->GetN(), s4BarXs.at(i), s4BarTops.at(i)-0.1);
    grBar->SetPoint(grBar->GetN(), s4BarXs.at(i)-0.01, s4BarTops.at(i)-0.1);
    grBar->SetPoint(grBar->GetN(), s4BarXs.at(i)-0.01, s4BarTops.at(i));
    grBar->SetPoint(grBar->GetN(), s4BarXs.at(i), s4BarTops.at(i));
    mg->Add(grBar);
  }

  fout->cd();
  h3BarTot->Write();
  for(int b=0; b<10; b++) {
    h3BarVec.at(b)->Write();
  }
  hBarHits->Write();
  h2CosPhiThetaDist->Write();
  hCosPhiDist->Write();
  hCosDist->Write();
  h2PhiThetaSpatialDist->Write();
  h2EndPos->Write();
  h3StartPos->Write();
  h3StartDir->Write();
  h2StartXY->Write();
  h2StartXZ->Write();
  h2BarCoins->Write();
  mg->Write();
  fout->Close();
  delete fout;

} // toyS4Cosmics
