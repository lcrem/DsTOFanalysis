// overlapCalcLines.C

double atanErr(const double opp, const double adj, const double inerr) {
  double err = TMath::Sqrt(pow(adj/(pow(adj, 2)+pow(opp, 2)), 2) * pow(inerr, 2) +
			   pow(opp/(pow(adj, 2)+pow(opp, 2)), 2) * pow(inerr, 2));
  return err;
}

void overlapCalcLines(const char* outDir) {

  gSystem->Load("libPhysics.so");

  TFile *fout = new TFile(Form("%s/overlapCalcLines.root", outDir), "recreate");
  
  // Wire chamber points
  TVector3 WC_ULB(-2.0563, 0.1632, -0.0801);
  TVector3 WC_ULT(-2.0554, 0.1634,  0.0636);
  TVector3 WC_UTL(-2.0547, 0.1287,  0.1013);
  TVector3 WC_UBL(-2.0557, 0.1159, -0.1230);
  TVector3 WC_UTR(-2.0523, -0.0111, 0.1020);
  TVector3 WC_URT(-2.0519, -0.0583, 0.0655);
  TVector3 WC_URB(-2.0526, -0.0599, -0.0829);
  TVector3 WC_UBR(-2.0535, -0.0242, -0.1232);
  // Calculated centre for wire chamber
  TVector3 WC_C = (WC_ULB + WC_ULT + WC_UTL + WC_UBL + WC_UTR + WC_URT + WC_URB + WC_UBR) * 0.125;
  // S1 points
  // TVector3 vs1Centre(-1.765, 0.036, -0.002);
  TVector3 vs1_ULB(-1.7672, 0.0617, -0.0086);
  TVector3 vs1_ULT(-1.7681, 0.0615, 0.0043);
  TVector3 vs1_UTL(-1.7689, 0.0491, 0.0337);
  TVector3 vs1_UBL(-1.7654, 0.0370, -0.0383);
  TVector3 vs1_UTR(-1.7669, 0.0116, 0.0334);
  TVector3 vs1_URT(-1.7655, -0.0102, 0.0041);
  TVector3 vs1_URB(-1.7658, -0.0103, -0.0087);
  TVector3 vs1_UBR(-1.7646, 0.0037, -0.0382);
  // Outline box of S1
  TVector3 vs1TopLeft(-1.7672, 0.0615, 0.0337);
  TVector3 vs1TopRight(-1.7669, -0.0103, 0.0337);
  TVector3 vs1BottomLeft(-1.7672, 0.0617, -0.0383);
  TVector3 vs1BottomRight(-1.7672, -0.0103, -0.0383);
  // Calculated centre of S1
  TVector3 vs1_C = (vs1_ULB+vs1_ULT+vs1_UTL+vs1_UBL+vs1_UTR+vs1_URT+vs1_URB+vs1_UBR) * 0.125;
  
  TVector3 vs2TopLeft(-0.3472, 0.0344, 0.0706);
  TVector3 vs2_activeBL(-0.3472, 0.0344, -0.0496);
  TVector3 vs2TopRight(-0.3498, -0.0725, 0.0679);
  TVector3 vs2_activeBR(-0.3498, -0.0725, -0.0521);
  TVector3 vs2Bottom(-0.3503, -0.0133, -0.2928);
  TVector3 vs2ActiveCentre = (vs2TopLeft+vs2TopRight+vs2_activeBL+vs2_activeBR)*0.25;
  cout<<vs2ActiveCentre.X()<<", "<<vs2ActiveCentre.Y()<<", "<<vs2ActiveCentre.Z()<<endl;
  TVector3 D3_UTR(8.9245, -0.9928, 0.6220);
  TVector3 D3_UBR(8.9047, -1.0012, -0.6012);
  TVector3 D3_UBL(9.0488, 0.5120, -0.5993);
  TVector3 D3_UTL(9.0650, 0.5215, 0.6244);
  double l1 = TMath::Sqrt(pow(D3_UTR.X()-D3_UTL.X(), 2) + pow(D3_UTR.Y()-D3_UTL.Y(), 2));
  double l2 = TMath::Sqrt(pow(D3_UBR.X()-D3_UBL.X(), 2) + pow(D3_UBR.Y()-D3_UBL.Y(), 2));
  double grad = (D3_UTR.Y()-D3_UTL.Y())/(D3_UTR.X()-D3_UTL.X());
  double xExt = 0.08/TMath::Sqrt(1+pow(grad,2));
  double yExt = xExt * grad;
  cout<<"x, y extensions "<<xExt<<", "<<yExt<<endl;
  TVector3 s3Ext(xExt, yExt, 0.);
  TVector3 vs3TopLeft = D3_UTL + s3Ext;
  TVector3 vs3TopRight = D3_UTR - s3Ext;
  TVector3 vs3BottomLeft = D3_UBL + s3Ext;
  TVector3 vs3BottomRight = D3_UBR - s3Ext;

  cout<<"S3 TL X, Y, Z "<<vs3TopLeft.X()<<", "<<vs3TopLeft.Y()<<", "<<vs3TopLeft.Z()<<endl;
  cout<<"S3 TR X, Y, Z "<<vs3TopRight.X()<<", "<<vs3TopRight.Y()<<", "<<vs3TopRight.Z()<<endl;
  cout<<"S3 BL X, Y, Z "<<vs3BottomLeft.X()<<", "<<vs3BottomLeft.Y()<<", "<<vs3BottomLeft.Z()<<endl;
  cout<<"S3 BR X, Y, Z "<<vs3BottomRight.X()<<", "<<vs3BottomRight.Y()<<", "<<vs3BottomRight.Z()<<endl;

  TVector3 vs4_UTL(12.2480, -0.1227, 0.7751);
  TVector3 vs4_UBL(12.2258, -0.1194, -0.8390);
  TVector3 vs4_UBR(12.1679, -1.4233, -0.8456);
  TVector3 vs4_UTR(12.1773, -1.4214, 0.7774);
  TVector3 vs4_U1L(12.1776, -0.2796, 0.4325);
  TVector3 vs4_U1R(12.1326, -1.1621, 0.4312);
  TVector3 vs4_U2L(12.1542, -0.5110, 0.2803);
  TVector3 vs4_U2R(12.1090, -1.1521, 0.2788);
  TVector3 vs4_U5L(12.1594, -0.4575, -0.1724);
  TVector3 vs4_U5R(12.1220, -1.0901, -0.1719);
  TVector3 vs4_D1L(12.3797, -0.4040, 0.3483);
  TVector3 vs4_D1R(12.3385, -1.1282, 0.3510);
  TVector3 vs4_D5L(12.4000, -0.3797, -0.2530);
  TVector3 vs4_D5R(12.3576, -1.1521, -0.2506);
  // Vectors for active area of S4
  TVector3 vs4_activeTL(12.28, -0.0727, 0.432);
  TVector3 vs4_activeTR(12.28, -1.4715, 0.432);
  TVector3 vs4_activeBL(12.28, -0.0727, -0.352);
  TVector3 vs4_activeBR(12.28, -1.4715, -0.352);
  // Vectors for TPC drift volume
  TVector3 TPCActiveBottomUL(10.3048, -0.2778, -0.5614);
  TVector3 TPCActiveBottomDL(11.3991, -0.3896, -0.5614);
  TVector3 TPCActiveTopUL(10.3048, -0.2778, 0.5386);
  TVector3 TPCActiveTopDL(11.3991, -0.3896, 0.5386);
  TVector3 TPCActiveBottomUR(10.2587, -0.7290, -0.5614);
  TVector3 TPCActiveBottomDR(11.3530, -0.8408, -0.5614);
  TVector3 TPCActiveTopUR(10.2587, -0.7290, 0.5386);
  TVector3 TPCActiveTopDR(11.3530, -0.8408, 0.5386);
  TVector3 TPCActiveC(10.82885, -0.5593, -0.0114);
  // Vectors for vessel used in simulation
  TVector3 vesselC(10.83543, -0.49096, -0.0114);
  TVector3 vesselTopLeft(10.89599, 0.13863, 0.6886);
  TVector3 vesselTopRight(10.77487, -1.12055, 0.6886);
  TVector3 vesselBottomLeft(10.89599, 0.13863, -0.7114);
  TVector3 vesselBottomRight(10.77487, -1.12055, -0.7114);

  std::vector<TVector3> vesselVec = {vesselTopLeft, vesselTopRight, 
				     vesselBottomRight, vesselBottomLeft};

  std::vector<TVector3> tpcUsVec = {TPCActiveTopUL, TPCActiveTopUR,
				    TPCActiveBottomUR, TPCActiveBottomUL};

  std::vector<TVector3> tpcDsVec = {TPCActiveTopDL, TPCActiveTopDR,
				    TPCActiveBottomDR, TPCActiveBottomDL};
  
  std::vector<TVector3> s1Vec = {vs1TopLeft, vs1TopRight, vs1BottomRight, vs1BottomLeft}; 
 
  std::vector<TVector3> s2Vec = {vs2TopLeft, vs2_activeBL, vs2_activeBR, vs2TopRight};
  
  std::vector<TVector3> s3Vec = {vs3TopLeft, vs3TopRight, vs3BottomRight, vs3BottomLeft};
  
  std::vector<TVector3> s4Vec;
  s4Vec.push_back(vs4_UTL);
  s4Vec.push_back(vs4_UBL);
  s4Vec.push_back(vs4_UTR);
  s4Vec.push_back(vs4_UBR);

  std::vector<TVector3> s4ActiveVec = {vs4_activeTL, vs4_activeTR, vs4_activeBR, vs4_activeBL,};

  std::cout<<"\n";
  std::cout<<"WC as the origin ("<<WC_C.X()<<", "<<WC_C.Y()<<", "<<WC_C.Z()<<")"<<std::endl;
  TMultiGraph *mg_noProj = new TMultiGraph("mg_noProj", "Positions of objects in beamline; -x / m; y / m");
  TGraph *grs1_noProj = new TGraph();
  grs1_noProj->SetLineColor(kBlack);
  grs1_noProj->SetLineWidth(2);
  for (int i=0; i<s1Vec.size(); i++) {
    grs1_noProj->SetPoint(grs1_noProj->GetN(), s1Vec.at(i).Y()*-1, s1Vec.at(i).Z());
  }
  grs1_noProj->SetPoint(grs1_noProj->GetN(), s1Vec.at(0).Y()*-1, s1Vec.at(0).Z());
  mg_noProj->Add(grs1_noProj);
  
  TGraph *grs2_noProj = new TGraph();
  grs2_noProj->SetLineColor(kOrange);
  grs2_noProj->SetLineWidth(2);
  for (int i=0; i<s2Vec.size(); i++) {
    grs2_noProj->SetPoint(grs2_noProj->GetN(), s2Vec.at(i).Y()*-1, s2Vec.at(i).Z());
  }
  grs2_noProj->SetPoint(grs2_noProj->GetN(), s2Vec.at(0).Y()*-1, s2Vec.at(0).Z());
  mg_noProj->Add(grs2_noProj);

  TGraph *grs3_noProj = new TGraph();
  grs3_noProj->SetLineColor(kGreen+2);
  grs3_noProj->SetLineWidth(2);
  for (int i=0; i<s3Vec.size(); i++) {
    grs3_noProj->SetPoint(grs3_noProj->GetN(), s3Vec.at(i).Y()*-1, s3Vec.at(i).Z());
  }
  grs3_noProj->SetPoint(grs3_noProj->GetN(), s3Vec.at(0).Y()*-1, s3Vec.at(0).Z());
  mg_noProj->Add(grs3_noProj);

  TGraph *grs4_noProj = new TGraph();
  grs4_noProj->SetLineColor(kRed+2);
  grs4_noProj->SetLineWidth(2);
  for (int i=0; i<s4ActiveVec.size(); i++) {
    grs4_noProj->SetPoint(grs4_noProj->GetN(), s4ActiveVec.at(i).Y()*-1, s4ActiveVec.at(i).Z());
  }
  grs4_noProj->SetPoint(grs4_noProj->GetN(), s4ActiveVec.at(0).Y()*-1, s4ActiveVec.at(0).Z());
  mg_noProj->Add(grs4_noProj);

  TGraph *grTpcUs_noProj = new TGraph();
  grTpcUs_noProj->SetLineColor(kBlue);
  grTpcUs_noProj->SetLineWidth(2);
  for (int i=0; i<tpcUsVec.size(); i++) {
    grTpcUs_noProj->SetPoint(grTpcUs_noProj->GetN(), tpcUsVec.at(i).Y()*-1, tpcUsVec.at(i).Z());
  }
  grTpcUs_noProj->SetPoint(grTpcUs_noProj->GetN(), tpcUsVec.at(0).Y()*-1, tpcUsVec.at(0).Z());
  mg_noProj->Add(grTpcUs_noProj);

  TGraph *grTpcDs_noProj = new TGraph();
  grTpcDs_noProj->SetLineColor(kBlue+2);
  grTpcDs_noProj->SetLineWidth(2);
  for (int i=0; i<tpcDsVec.size(); i++) {
    grTpcDs_noProj->SetPoint(grTpcDs_noProj->GetN(), tpcDsVec.at(i).Y()*-1, tpcDsVec.at(i).Z());
  }
  grTpcDs_noProj->SetPoint(grTpcDs_noProj->GetN(), tpcDsVec.at(0).Y()*-1, tpcDsVec.at(0).Z());
  mg_noProj->Add(grTpcDs_noProj);

  // Angular positions with origin at the wire chamber
  TMultiGraph *mg_Ang = new TMultiGraph("mg_Ang", "Positions of objects in beamline (wire chamber origin); #theta / degrees; #phi / degrees");
  TGraph *grs1_Ang = new TGraph();
  grs1_Ang->SetTitle("S1");
  grs1_Ang->SetLineColor(kBlack);
  grs1_Ang->SetLineWidth(2);
  for (int i=0; i<s1Vec.size(); i++) {
    grs1_Ang->SetPoint(grs1_Ang->GetN(),
		       TMath::ATan2((s1Vec.at(i).Y()-WC_C.Y())*-1, s1Vec.at(i).X()-WC_C.X())*180./TMath::Pi(),
		       TMath::ATan2(s1Vec.at(i).Z()-WC_C.Z(), s1Vec.at(i).X()-WC_C.X())*180./TMath::Pi());
    std::cout<<"S1 theta, phi "<<TMath::ATan2((s1Vec.at(i).Y()-WC_C.Y())*-1, s1Vec.at(i).X()-WC_C.X())*180./TMath::Pi()<<" +/- "<<atanErr((s1Vec.at(i).Y()-WC_C.Y())*-1, s1Vec.at(i).X()-WC_C.X(), 0.00071)*180./TMath::Pi()<<", "<<TMath::ATan2(s1Vec.at(i).Z()-WC_C.Z(), s1Vec.at(i).X()-WC_C.X())*180./TMath::Pi()<<" +/- "<<atanErr((s1Vec.at(i).Z()-WC_C.Z())*-1, s1Vec.at(i).X()-WC_C.X(), 0.00071)*180./TMath::Pi()<<std::endl;
  }
  grs1_Ang->SetPoint(grs1_Ang->GetN(),
		     TMath::ATan2((s1Vec.at(0).Y()-WC_C.Y())*-1, s1Vec.at(0).X()-WC_C.X())*180./TMath::Pi(),
		     TMath::ATan2(s1Vec.at(0).Z()-WC_C.Z(), s1Vec.at(0).X()-WC_C.X())*180./TMath::Pi());
  mg_Ang->Add(grs1_Ang);
  std::cout<<"\n";
  
  TGraph *grs2_Ang = new TGraph();
  grs2_Ang->SetTitle("S2");
  grs2_Ang->SetLineColor(kOrange);
  grs2_Ang->SetLineWidth(2);
  for (int i=0; i<s2Vec.size(); i++) {
    grs2_Ang->SetPoint(grs2_Ang->GetN(),
		       TMath::ATan2((s2Vec.at(i).Y()-WC_C.Y())*-1, s2Vec.at(i).X()-WC_C.X())*180./TMath::Pi(),
		       TMath::ATan2(s2Vec.at(i).Z()-WC_C.Z(), s2Vec.at(i).X()-WC_C.X())*180./TMath::Pi());
    std::cout<<"S2 theta, phi "<<TMath::ATan2((s2Vec.at(i).Y()-WC_C.Y())*-1, s2Vec.at(i).X()-WC_C.X())*180./TMath::Pi()<<" +/- "<<atanErr((s2Vec.at(i).Y()-WC_C.Y())*-1, s2Vec.at(i).X()-WC_C.X(), 0.00071)*180./TMath::Pi()<<", "<<TMath::ATan2(s2Vec.at(i).Z()-WC_C.Z(), s2Vec.at(i).X()-WC_C.X())*180./TMath::Pi()<<" +/- "<<atanErr((s2Vec.at(i).Z()-WC_C.Z())*-1, s2Vec.at(i).X()-WC_C.X(), 0.00071)*180./TMath::Pi()<<std::endl;
  }
  grs2_Ang->SetPoint(grs2_Ang->GetN(),
		     TMath::ATan2((s2Vec.at(0).Y()-WC_C.Y())*-1, s2Vec.at(0).X()-WC_C.X())*180./TMath::Pi(),
		     TMath::ATan2(s2Vec.at(0).Z()-WC_C.Z(), s2Vec.at(0).X()-WC_C.X())*180./TMath::Pi());
  mg_Ang->Add(grs2_Ang);
  std::cout<<"\n";

  TGraph *grs3_Ang = new TGraph();
  grs3_Ang->SetTitle("S3");
  grs3_Ang->SetLineColor(kCyan+2);
  grs3_Ang->SetLineWidth(2);
  for (int i=0; i<s3Vec.size(); i++) {
    grs3_Ang->SetPoint(grs3_Ang->GetN(),
		       TMath::ATan2((s3Vec.at(i).Y()-WC_C.Y())*-1, s3Vec.at(i).X()-WC_C.X())*180./TMath::Pi(),
		       TMath::ATan2(s3Vec.at(i).Z()-WC_C.Z(), s3Vec.at(i).X()-WC_C.X())*180./TMath::Pi());
    std::cout<<"S3 theta, phi "<<TMath::ATan2((s3Vec.at(i).Y()-WC_C.Y())*-1, s3Vec.at(i).X()-WC_C.X())*180./TMath::Pi()<<" +/- "<<atanErr((s3Vec.at(i).Y()-WC_C.Y())*-1, s3Vec.at(i).X()-WC_C.X(), 0.00071)*180./TMath::Pi()<<", "<<TMath::ATan2(s3Vec.at(i).Z()-WC_C.Z(), s3Vec.at(i).X()-WC_C.X())*180./TMath::Pi()<<" +/- "<<atanErr((s3Vec.at(i).Z()-WC_C.Z())*-1, s3Vec.at(i).X()-WC_C.X(), 0.00071)*180./TMath::Pi()<<std::endl;
  }
  grs3_Ang->SetPoint(grs3_Ang->GetN(),
		     TMath::ATan2((s3Vec.at(0).Y()-WC_C.Y())*-1, s3Vec.at(0).X()-WC_C.X())*180./TMath::Pi(),
		     TMath::ATan2(s3Vec.at(0).Z()-WC_C.Z(), s3Vec.at(0).X()-WC_C.X())*180./TMath::Pi());
  mg_Ang->Add(grs3_Ang);
  std::cout<<"\n";
  
  TGraph *grs4_Ang = new TGraph();
  grs4_Ang->SetTitle("S4");
  grs4_Ang->SetLineColor(kRed+2);
  grs4_Ang->SetLineWidth(2);
  for (int i=0; i<s4ActiveVec.size(); i++) {
    grs4_Ang->SetPoint(grs4_Ang->GetN(),
		       TMath::ATan2((s4ActiveVec.at(i).Y()-WC_C.Y())*-1, s4ActiveVec.at(i).X()-WC_C.X())*180./TMath::Pi(),
		       TMath::ATan2(s4ActiveVec.at(i).Z()-WC_C.Z(), s4ActiveVec.at(i).X()-WC_C.X())*180./TMath::Pi());
    std::cout<<"S4 theta, phi "<<TMath::ATan2((s4ActiveVec.at(i).Y()-WC_C.Y())*-1, s4ActiveVec.at(i).X()-WC_C.X())*180./TMath::Pi()<<" +/- "<<atanErr((s4ActiveVec.at(i).Y()-WC_C.Y())*-1, s4ActiveVec.at(i).X()-WC_C.X(), 0.00071)*180./TMath::Pi()<<", "<<TMath::ATan2(s4ActiveVec.at(i).Z()-WC_C.Z(), s4ActiveVec.at(i).X()-WC_C.X())*180./TMath::Pi()<<" +/- "<<atanErr((s4ActiveVec.at(i).Z()-WC_C.Z())*-1, s4ActiveVec.at(i).X()-WC_C.X(), 0.00071)*180./TMath::Pi()<<std::endl;
  }
  grs4_Ang->SetPoint(grs4_Ang->GetN(), TMath::ATan2((s4ActiveVec.at(0).Y()-WC_C.Y())*-1, s4ActiveVec.at(0).X()-WC_C.X())*180./TMath::Pi(), TMath::ATan2(s4ActiveVec.at(0).Z()-WC_C.Z(), s4ActiveVec.at(0).X()-WC_C.X())*180./TMath::Pi());
  mg_Ang->Add(grs4_Ang);
  std::cout<<"\n";

  TGraph *grtpcUs_Ang = new TGraph();
  grtpcUs_Ang->SetTitle("TPC US");
  grtpcUs_Ang->SetLineColor(kMagenta+1);
  grtpcUs_Ang->SetLineStyle(7);
  grtpcUs_Ang->SetLineWidth(2);
  for (int i=0; i<tpcUsVec.size(); i++) {
    grtpcUs_Ang->SetPoint(grtpcUs_Ang->GetN(), TMath::ATan2((tpcUsVec.at(i).Y()-WC_C.Y())*-1, tpcUsVec.at(i).X()-WC_C.X())*180./TMath::Pi(), TMath::ATan2(tpcUsVec.at(i).Z()-WC_C.Z(), tpcUsVec.at(i).X()-WC_C.X())*180./TMath::Pi());
    std::cout<<"TPC US theta, phi "<<TMath::ATan2((tpcUsVec.at(i).Y()-WC_C.Y())*-1, tpcUsVec.at(i).X()-WC_C.X())*180./TMath::Pi()<<" +/- "<<atanErr((tpcUsVec.at(i).Y()-WC_C.Y())*-1, tpcUsVec.at(i).X()-WC_C.X(), 0.00071)*180./TMath::Pi()<<", "<<TMath::ATan2(tpcUsVec.at(i).Z()-WC_C.Z(), tpcUsVec.at(i).X()-WC_C.X())*180./TMath::Pi()<<" +/- "<<atanErr((tpcUsVec.at(i).Z()-WC_C.Z())*-1, tpcUsVec.at(i).X()-WC_C.X(), 0.00071)*180./TMath::Pi()<<std::endl;
  }
  grtpcUs_Ang->SetPoint(grtpcUs_Ang->GetN(), TMath::ATan2((tpcUsVec.at(0).Y()-WC_C.Y())*-1, tpcUsVec.at(0).X()-WC_C.X())*180./TMath::Pi(), TMath::ATan2(tpcUsVec.at(0).Z()-WC_C.Z(), tpcUsVec.at(0).X()-WC_C.X())*180./TMath::Pi());
  mg_Ang->Add(grtpcUs_Ang);
  std::cout<<"\n";
  
  TGraph *grtpcDs_Ang = new TGraph();
  grtpcDs_Ang->SetTitle("TPC DS");
  grtpcDs_Ang->SetLineColor(kMagenta+1);
  grtpcDs_Ang->SetLineWidth(2);
  for (int i=0; i<tpcDsVec.size(); i++) {
    grtpcDs_Ang->SetPoint(grtpcDs_Ang->GetN(), TMath::ATan2((tpcDsVec.at(i).Y()-WC_C.Y())*-1, tpcDsVec.at(i).X()-WC_C.X())*180./TMath::Pi(), TMath::ATan2(tpcDsVec.at(i).Z()-WC_C.Z(), tpcDsVec.at(i).X()-WC_C.X())*180./TMath::Pi());
    std::cout<<"TPC DS theta, phi "<<TMath::ATan2((tpcDsVec.at(i).Y()-WC_C.Y())*-1, tpcDsVec.at(i).X()-WC_C.X())*180./TMath::Pi()<<" +/- "<<atanErr((tpcDsVec.at(i).Y()-WC_C.Y())*-1, tpcDsVec.at(i).X()-WC_C.X(), 0.00071)*180./TMath::Pi()<<", "<<TMath::ATan2(tpcDsVec.at(i).Z()-WC_C.Z(), tpcDsVec.at(i).X()-WC_C.X())*180./TMath::Pi()<<" +/- "<<atanErr((tpcDsVec.at(i).Z()-WC_C.Z())*-1, tpcDsVec.at(i).X()-WC_C.X(), 0.00071)*180./TMath::Pi()<<std::endl;
  }
  grtpcDs_Ang->SetPoint(grtpcDs_Ang->GetN(), TMath::ATan2((tpcDsVec.at(0).Y()-WC_C.Y())*-1, tpcDsVec.at(0).X()-WC_C.X())*180./TMath::Pi(), TMath::ATan2(tpcDsVec.at(0).Z()-WC_C.Z(), tpcDsVec.at(0).X()-WC_C.X())*180./TMath::Pi());
  mg_Ang->Add(grtpcDs_Ang);  

  TGraph *grvessel_Ang = new TGraph();
  grvessel_Ang->SetTitle("Vessel");
  grvessel_Ang->SetLineColor(kGreen+2);
  grvessel_Ang->SetLineWidth(2);
  for (int i=0; i<vesselVec.size(); i++) {
    grvessel_Ang->SetPoint(grvessel_Ang->GetN(), TMath::ATan2((vesselVec.at(i).Y()-WC_C.Y())*-1, vesselVec.at(i).X()-WC_C.X())*180./TMath::Pi(), TMath::ATan2(vesselVec.at(i).Z()-WC_C.Z(), vesselVec.at(i).X()-WC_C.X())*180./TMath::Pi());
  }
  grvessel_Ang->SetPoint(grvessel_Ang->GetN(), TMath::ATan2((vesselVec.at(0).Y()-WC_C.Y())*-1, vesselVec.at(0).X()-WC_C.X())*180./TMath::Pi(), TMath::ATan2(vesselVec.at(0).Z()-WC_C.Z(), vesselVec.at(0).X()-WC_C.X())*180./TMath::Pi());
  mg_Ang->Add(grvessel_Ang); 

  std::cout<<"\n";
  // Now do the same thing but using S1 as the origin
  std::cout<<"Now with S1 as the origin ("<<vs1_C.X()<<", "<<vs1_C.Y()<<", "<<vs1_C.Z()<<")"<<std::endl;
  TMultiGraph *mg_AngS1 = new TMultiGraph("mg_AngS1", "Positions of objects in beamline (S1 origin); #theta / degrees; #phi / degrees");
  TGraph *grs2_AngS1 = new TGraph();
  grs2_AngS1->SetTitle("S2");
  grs2_AngS1->SetLineColor(kOrange);
  grs2_AngS1->SetLineWidth(2);
  for (int i=0; i<s2Vec.size(); i++) {
    grs2_AngS1->SetPoint(grs2_AngS1->GetN(),
			 TMath::ATan2((s2Vec.at(i).Y()-vs1_C.Y())*-1, s2Vec.at(i).X()-vs1_C.X())*180./TMath::Pi(),
			 TMath::ATan2(s2Vec.at(i).Z()-vs1_C.Z(), s2Vec.at(i).X()-vs1_C.X())*180./TMath::Pi());
    std::cout<<"S2 theta, phi "<<TMath::ATan2((s2Vec.at(i).Y()-vs1_C.Y())*-1, s2Vec.at(i).X()-vs1_C.X())*180./TMath::Pi()<<" +/- "<<atanErr((s2Vec.at(i).Y()-vs1_C.Y())*-1, s2Vec.at(i).X()-vs1_C.X(), 0.00071)*180./TMath::Pi()<<", "<<TMath::ATan2(s2Vec.at(i).Z()-vs1_C.Z(), s2Vec.at(i).X()-vs1_C.X())*180./TMath::Pi()<<" +/- "<<atanErr((s2Vec.at(i).Z()-vs1_C.Z())*-1, s2Vec.at(i).X()-vs1_C.X(), 0.00071)*180./TMath::Pi()<<std::endl;
  }
  grs2_AngS1->SetPoint(grs2_AngS1->GetN(),
		       TMath::ATan2((s2Vec.at(0).Y()-vs1_C.Y())*-1, s2Vec.at(0).X()-vs1_C.X())*180./TMath::Pi(),
		       TMath::ATan2(s2Vec.at(0).Z()-vs1_C.Z(), s2Vec.at(0).X()-vs1_C.X())*180./TMath::Pi());
  mg_AngS1->Add(grs2_AngS1);
  std::cout<<"\n";

  TGraph *grs3_AngS1 = new TGraph();
  grs3_AngS1->SetTitle("S3");
  grs3_AngS1->SetLineColor(kCyan+2);
  grs3_AngS1->SetLineWidth(2);
  for (int i=0; i<s3Vec.size(); i++) {
    grs3_AngS1->SetPoint(grs3_AngS1->GetN(),
			 TMath::ATan2((s3Vec.at(i).Y()-vs1_C.Y())*-1, s3Vec.at(i).X()-vs1_C.X())*180./TMath::Pi(),
			 TMath::ATan2(s3Vec.at(i).Z()-vs1_C.Z(), s3Vec.at(i).X()-vs1_C.X())*180./TMath::Pi());
    std::cout<<"S3 theta, phi "<<TMath::ATan2((s3Vec.at(i).Y()-vs1_C.Y())*-1, s3Vec.at(i).X()-vs1_C.X())*180./TMath::Pi()<<" +/- "<<atanErr((s3Vec.at(i).Y()-vs1_C.Y())*-1, s3Vec.at(i).X()-vs1_C.X(), 0.00071)*180./TMath::Pi()<<", "<<TMath::ATan2(s3Vec.at(i).Z()-vs1_C.Z(), s3Vec.at(i).X()-vs1_C.X())*180./TMath::Pi()<<" +/- "<<atanErr((s3Vec.at(i).Z()-vs1_C.Z())*-1, s3Vec.at(i).X()-vs1_C.X(), 0.00071)*180./TMath::Pi()<<std::endl;
  }
  grs3_AngS1->SetPoint(grs3_AngS1->GetN(),
		       TMath::ATan2((s3Vec.at(0).Y()-vs1_C.Y())*-1, s3Vec.at(0).X()-vs1_C.X())*180./TMath::Pi(),
		       TMath::ATan2(s3Vec.at(0).Z()-vs1_C.Z(), s3Vec.at(0).X()-vs1_C.X())*180./TMath::Pi());
  mg_AngS1->Add(grs3_AngS1);
  std::cout<<"\n";
  
  TGraph *grs4_AngS1 = new TGraph();
  grs4_AngS1->SetTitle("S4");
  grs4_AngS1->SetLineColor(kRed+2);
  grs4_AngS1->SetLineWidth(2);
  for (int i=0; i<s4ActiveVec.size(); i++) {
    grs4_AngS1->SetPoint(grs4_AngS1->GetN(),
			 TMath::ATan2((s4ActiveVec.at(i).Y()-vs1_C.Y())*-1, s4ActiveVec.at(i).X()-vs1_C.X())*180./TMath::Pi(),
			 TMath::ATan2(s4ActiveVec.at(i).Z()-vs1_C.Z(), s4ActiveVec.at(i).X()-vs1_C.X())*180./TMath::Pi());
    std::cout<<"S4 theta, phi "<<TMath::ATan2((s4ActiveVec.at(i).Y()-vs1_C.Y())*-1, s4ActiveVec.at(i).X()-vs1_C.X())*180./TMath::Pi()<<" +/- "<<atanErr((s4ActiveVec.at(i).Y()-vs1_C.Y())*-1, s4ActiveVec.at(i).X()-vs1_C.X(), 0.00071)*180./TMath::Pi()<<", "<<TMath::ATan2(s4ActiveVec.at(i).Z()-vs1_C.Z(), s4ActiveVec.at(i).X()-vs1_C.X())*180./TMath::Pi()<<" +/- "<<atanErr((s4ActiveVec.at(i).Z()-vs1_C.Z())*-1, s4ActiveVec.at(i).X()-vs1_C.X(), 0.00071)*180./TMath::Pi()<<std::endl;
  }
  grs4_AngS1->SetPoint(grs4_AngS1->GetN(), TMath::ATan2((s4ActiveVec.at(0).Y()-vs1_C.Y())*-1, s4ActiveVec.at(0).X()-vs1_C.X())*180./TMath::Pi(), TMath::ATan2(s4ActiveVec.at(0).Z()-vs1_C.Z(), s4ActiveVec.at(0).X()-vs1_C.X())*180./TMath::Pi());
  mg_AngS1->Add(grs4_AngS1);
  std::cout<<"\n";

  TGraph *grtpcUs_AngS1 = new TGraph();
  grtpcUs_AngS1->SetTitle("TPC US");
  grtpcUs_AngS1->SetLineColor(kMagenta+1);
  grtpcUs_AngS1->SetLineStyle(7);
  grtpcUs_AngS1->SetLineWidth(2);
  for (int i=0; i<tpcUsVec.size(); i++) {
    grtpcUs_AngS1->SetPoint(grtpcUs_AngS1->GetN(), TMath::ATan2((tpcUsVec.at(i).Y()-vs1_C.Y())*-1, tpcUsVec.at(i).X()-vs1_C.X())*180./TMath::Pi(), TMath::ATan2(tpcUsVec.at(i).Z()-vs1_C.Z(), tpcUsVec.at(i).X()-vs1_C.X())*180./TMath::Pi());
    std::cout<<"TPC US theta, phi "<<TMath::ATan2((tpcUsVec.at(i).Y()-vs1_C.Y())*-1, tpcUsVec.at(i).X()-vs1_C.X())*180./TMath::Pi()<<" +/- "<<atanErr((tpcUsVec.at(i).Y()-vs1_C.Y())*-1, tpcUsVec.at(i).X()-vs1_C.X(), 0.00071)*180./TMath::Pi()<<", "<<TMath::ATan2(tpcUsVec.at(i).Z()-vs1_C.Z(), tpcUsVec.at(i).X()-vs1_C.X())*180./TMath::Pi()<<" +/- "<<atanErr((tpcUsVec.at(i).Z()-vs1_C.Z())*-1, tpcUsVec.at(i).X()-vs1_C.X(), 0.00071)*180./TMath::Pi()<<std::endl;
  }
  grtpcUs_AngS1->SetPoint(grtpcUs_AngS1->GetN(), TMath::ATan2((tpcUsVec.at(0).Y()-vs1_C.Y())*-1, tpcUsVec.at(0).X()-vs1_C.X())*180./TMath::Pi(), TMath::ATan2(tpcUsVec.at(0).Z()-vs1_C.Z(), tpcUsVec.at(0).X()-vs1_C.X())*180./TMath::Pi());
  mg_AngS1->Add(grtpcUs_AngS1);
  std::cout<<"\n";
  
  TGraph *grtpcDs_AngS1 = new TGraph();
  grtpcDs_AngS1->SetTitle("TPC DS");
  grtpcDs_AngS1->SetLineColor(kMagenta+1);
  grtpcDs_AngS1->SetLineWidth(2);
  for (int i=0; i<tpcDsVec.size(); i++) {
    grtpcDs_AngS1->SetPoint(grtpcDs_AngS1->GetN(), TMath::ATan2((tpcDsVec.at(i).Y()-vs1_C.Y())*-1, tpcDsVec.at(i).X()-vs1_C.X())*180./TMath::Pi(), TMath::ATan2(tpcDsVec.at(i).Z()-vs1_C.Z(), tpcDsVec.at(i).X()-vs1_C.X())*180./TMath::Pi());
    std::cout<<"TPC DS theta, phi "<<TMath::ATan2((tpcDsVec.at(i).Y()-vs1_C.Y())*-1, tpcDsVec.at(i).X()-vs1_C.X())*180./TMath::Pi()<<" +/- "<<atanErr((tpcDsVec.at(i).Y()-vs1_C.Y())*-1, tpcDsVec.at(i).X()-vs1_C.X(), 0.00071)*180./TMath::Pi()<<", "<<TMath::ATan2(tpcDsVec.at(i).Z()-vs1_C.Z(), tpcDsVec.at(i).X()-vs1_C.X())*180./TMath::Pi()<<" +/- "<<atanErr((tpcDsVec.at(i).Z()-vs1_C.Z())*-1, tpcDsVec.at(i).X()-vs1_C.X(), 0.00071)*180./TMath::Pi()<<std::endl;
  }
  grtpcDs_AngS1->SetPoint(grtpcDs_AngS1->GetN(), TMath::ATan2((tpcDsVec.at(0).Y()-vs1_C.Y())*-1, tpcDsVec.at(0).X()-vs1_C.X())*180./TMath::Pi(), TMath::ATan2(tpcDsVec.at(0).Z()-vs1_C.Z(), tpcDsVec.at(0).X()-vs1_C.X())*180./TMath::Pi());
  mg_AngS1->Add(grtpcDs_AngS1);  

  TGraph *grvessel_AngS1 = new TGraph();
  grvessel_AngS1->SetTitle("Vessel");
  grvessel_AngS1->SetLineColor(kGreen+2);
  grvessel_AngS1->SetLineWidth(2);
  for (int i=0; i<vesselVec.size(); i++) {
    grvessel_AngS1->SetPoint(grvessel_AngS1->GetN(), TMath::ATan2((vesselVec.at(i).Y()-vs1_C.Y())*-1, vesselVec.at(i).X()-vs1_C.X())*180./TMath::Pi(), TMath::ATan2(vesselVec.at(i).Z()-vs1_C.Z(), vesselVec.at(i).X()-vs1_C.X())*180./TMath::Pi());
  }
  grvessel_AngS1->SetPoint(grvessel_AngS1->GetN(), TMath::ATan2((vesselVec.at(0).Y()-vs1_C.Y())*-1, vesselVec.at(0).X()-vs1_C.X())*180./TMath::Pi(), TMath::ATan2(vesselVec.at(0).Z()-vs1_C.Z(), vesselVec.at(0).X()-vs1_C.X())*180./TMath::Pi());
  mg_AngS1->Add(grvessel_AngS1); 
  
  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->SetGridx();
  c1->SetGridy();
  mg_noProj->Draw("al");

  TCanvas *c2 = new TCanvas("c2", "c2");
  mg_Ang->Draw("al");
  mg_Ang->GetXaxis()->SetTitleSize(0.05);
  mg_Ang->GetXaxis()->SetLabelSize(0.05);
  mg_Ang->GetYaxis()->SetTitleSize(0.05);
  mg_Ang->GetYaxis()->SetLabelSize(0.05);
  c2->SetLeftMargin(0.12);
  c2->SetBottomMargin(0.12);
  c2->Update();

  TCanvas *c3 = new TCanvas("c3", "c3");
  mg_AngS1->Draw("al");
  mg_AngS1->GetXaxis()->SetTitleSize(0.05);
  mg_AngS1->GetXaxis()->SetLabelSize(0.05);
  mg_AngS1->GetYaxis()->SetTitleSize(0.05);
  mg_AngS1->GetYaxis()->SetLabelSize(0.05);
  c3->SetLeftMargin(0.12);
  c3->SetBottomMargin(0.12);
  c3->Update();
  
  fout->cd();
  mg_noProj->Write();
  mg_Ang->Write();
  mg_AngS1->Write();
  
  fout->Close();
  delete fout;
}
