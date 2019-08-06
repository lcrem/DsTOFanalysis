// overlapCalcLines.C

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
  TVector3 vs1Centre(-1.765, 0.036, -0.002);
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
  
  TVector3 vs2TopLeft(-0.3472, 0.0344, 0.0706);
  TVector3 vs2_activeBL(-0.3472, 0.0344, -0.0496);
  TVector3 vs2TopRight(-0.3498, -0.0725, 0.0679);
  TVector3 vs2_activeBR(-0.3498, -0.0725, -0.0521);
  TVector3 vs2Bottom(-0.3503, -0.0133, -0.2928);

  TVector3 vs3TopLeft(9.0750, 0.6015, 0.6244);
  TVector3 vs3TopRight(8.9145, -1.0728, 0.6220);
  TVector3 vs3BottomLeft(9.0488, 0.5920, -0.5993);
  TVector3 vs3BottomRight(8.9147, -1.0812, -0.6012);

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
  }
  grs1_Ang->SetPoint(grs1_Ang->GetN(),
		     TMath::ATan2((s1Vec.at(0).Y()-WC_C.Y())*-1, s1Vec.at(0).X()-WC_C.X())*180./TMath::Pi(),
		     TMath::ATan2(s1Vec.at(0).Z()-WC_C.Z(), s1Vec.at(0).X()-WC_C.X())*180./TMath::Pi());
  mg_Ang->Add(grs1_Ang);
  
  TGraph *grs2_Ang = new TGraph();
  grs2_Ang->SetTitle("S2");
  grs2_Ang->SetLineColor(kOrange);
  grs2_Ang->SetLineWidth(2);
  for (int i=0; i<s2Vec.size(); i++) {
    grs2_Ang->SetPoint(grs2_Ang->GetN(),
		       TMath::ATan2((s2Vec.at(i).Y()-WC_C.Y())*-1, s2Vec.at(i).X()-WC_C.X())*180./TMath::Pi(),
		       TMath::ATan2(s2Vec.at(i).Z()-WC_C.Z(), s2Vec.at(i).X()-WC_C.X())*180./TMath::Pi());
  }
  grs2_Ang->SetPoint(grs2_Ang->GetN(),
		     TMath::ATan2((s2Vec.at(0).Y()-WC_C.Y())*-1, s2Vec.at(0).X()-WC_C.X())*180./TMath::Pi(),
		     TMath::ATan2(s2Vec.at(0).Z()-WC_C.Z(), s2Vec.at(0).X()-WC_C.X())*180./TMath::Pi());
  mg_Ang->Add(grs2_Ang);

  TGraph *grs3_Ang = new TGraph();
  grs3_Ang->SetTitle("S3");
  grs3_Ang->SetLineColor(kCyan+2);
  grs3_Ang->SetLineWidth(2);
  for (int i=0; i<s3Vec.size(); i++) {
    grs3_Ang->SetPoint(grs3_Ang->GetN(),
		       TMath::ATan2((s3Vec.at(i).Y()-WC_C.Y())*-1, s3Vec.at(i).X()-WC_C.X())*180./TMath::Pi(),
		       TMath::ATan2(s3Vec.at(i).Z()-WC_C.Z(), s3Vec.at(i).X()-WC_C.X())*180./TMath::Pi());
  }
  grs3_Ang->SetPoint(grs3_Ang->GetN(),
		     TMath::ATan2((s3Vec.at(0).Y()-WC_C.Y())*-1, s3Vec.at(0).X()-WC_C.X())*180./TMath::Pi(),
		     TMath::ATan2(s3Vec.at(0).Z()-WC_C.Z(), s3Vec.at(0).X()-WC_C.X())*180./TMath::Pi());
  mg_Ang->Add(grs3_Ang);

  TGraph *grs4_Ang = new TGraph();
  grs4_Ang->SetTitle("S4");
  grs4_Ang->SetLineColor(kRed+2);
  grs4_Ang->SetLineWidth(2);
  for (int i=0; i<s4ActiveVec.size(); i++) {
    grs4_Ang->SetPoint(grs4_Ang->GetN(),
		       TMath::ATan2((s4ActiveVec.at(i).Y()-WC_C.Y())*-1, s4ActiveVec.at(i).X()-WC_C.X())*180./TMath::Pi(),
		       TMath::ATan2(s4ActiveVec.at(i).Z()-WC_C.Z(), s4ActiveVec.at(i).X()-WC_C.X())*180./TMath::Pi());
  }
  grs4_Ang->SetPoint(grs4_Ang->GetN(),
		     TMath::ATan2((s4ActiveVec.at(0).Y()-WC_C.Y())*-1, s4ActiveVec.at(0).X()-WC_C.X())*180./TMath::Pi(),
		     TMath::ATan2(s4ActiveVec.at(0).Z()-WC_C.Z(), s4ActiveVec.at(0).X()-WC_C.X())*180./TMath::Pi());
  mg_Ang->Add(grs4_Ang);

  TGraph *grtpcUs_Ang = new TGraph();
  grtpcUs_Ang->SetTitle("TPC US");
  grtpcUs_Ang->SetLineColor(kMagenta+1);
  grtpcUs_Ang->SetLineStyle(7);
  grtpcUs_Ang->SetLineWidth(2);
  for (int i=0; i<tpcUsVec.size(); i++) {
    grtpcUs_Ang->SetPoint(grtpcUs_Ang->GetN(),
			  TMath::ATan2((tpcUsVec.at(i).Y()-WC_C.Y())*-1, tpcUsVec.at(i).X()-WC_C.X())*180./TMath::Pi(),
			  TMath::ATan2(tpcUsVec.at(i).Z()-WC_C.Z(), tpcUsVec.at(i).X()-WC_C.X())*180./TMath::Pi());
  }
  grtpcUs_Ang->SetPoint(grtpcUs_Ang->GetN(),
			TMath::ATan2((tpcUsVec.at(0).Y()-WC_C.Y())*-1, tpcUsVec.at(0).X()-WC_C.X())*180./TMath::Pi(),
			TMath::ATan2(tpcUsVec.at(0).Z()-WC_C.Z(), tpcUsVec.at(0).X()-WC_C.X())*180./TMath::Pi());
  mg_Ang->Add(grtpcUs_Ang);

  TGraph *grtpcDs_Ang = new TGraph();
  grtpcDs_Ang->SetTitle("TPC DS");
  grtpcDs_Ang->SetLineColor(kMagenta+1);
  grtpcDs_Ang->SetLineWidth(2);
  for (int i=0; i<tpcDsVec.size(); i++) {
    grtpcDs_Ang->SetPoint(grtpcDs_Ang->GetN(),
			  TMath::ATan2((tpcDsVec.at(i).Y()-WC_C.Y())*-1, tpcDsVec.at(i).X()-WC_C.X())*180./TMath::Pi(),
			  TMath::ATan2(tpcDsVec.at(i).Z()-WC_C.Z(), tpcDsVec.at(i).X()-WC_C.X())*180./TMath::Pi());
  }
  grtpcDs_Ang->SetPoint(grtpcDs_Ang->GetN(),
			TMath::ATan2((tpcDsVec.at(0).Y()-WC_C.Y())*-1, tpcDsVec.at(0).X()-WC_C.X())*180./TMath::Pi(),
			TMath::ATan2(tpcDsVec.at(0).Z()-WC_C.Z(), tpcDsVec.at(0).X()-WC_C.X())*180./TMath::Pi());
  mg_Ang->Add(grtpcDs_Ang);  
  
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
  
  fout->cd();
  mg_noProj->Write();
  mg_Ang->Write();
  
  fout->Close();
  delete fout;
}
