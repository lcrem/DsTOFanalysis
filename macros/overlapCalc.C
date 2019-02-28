// overlapCalc.C

void overlapCalc(const char* saveDir)
{
  gSystem->Load("libPhysics.so");
  
  TVector3 vs1Centre(-1.765, 0.036, -0.002);
  TVector3 vs1_ULB(-1.7672, 0.0617, -0.0086);
  TVector3 vs1_ULT(-1.7681, 0.0615, 0.0043);
  TVector3 vs1_UTL(-1.7689, 0.0491, 0.0337);
  TVector3 vs1_UBL(-1.7654, 0.0370, -0.0383);
  TVector3 vs1_UTR(-1.7669, 0.0116, 0.0334);
  TVector3 vs1_URT(-1.7655, -0.0102, 0.0041);
  TVector3 vs1_URB(-1.7658, -0.0103, -0.0087);
  TVector3 vs1_UBR(-1.7646, 0.0037, -0.0382);
  
  TVector3 vs2TopLeft(-0.3472, 0.0344, 0.0706);
  TVector3 vs2_activeBL(-0.3472, 0.0344, -0.0496);
  TVector3 vs2TopRight(-0.3498, -0.0725, 0.0679);
  TVector3 vs2_activeBR(-0.3498, -0.0725, -0.0521);
  TVector3 vs2Bottom(-0.3503, -0.0133, -0.2928);

  TVector3 vs3TopLeft(9.0650, 0.5215, 0.6244);
  TVector3 vs3TopRight(8.9245, -0.9928, 0.6220);
  TVector3 vs3BottomLeft(9.0488, 0.5120, -0.5993);
  TVector3 vs3BottomRight(8.9047, -1.0012, -0.6012);

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
  // Vectprs for active area of S4
  TVector3 vs4_activeTL(12.28, -0.0727, 0.432);
  TVector3 vs4_activeTR(12.28, -1.4715, 0.432);
  TVector3 vs4_activeBL(12.28, -0.0727, -0.352);
  TVector3 vs4_activeBR(12.28, -1.4715, -0.352);
  
  std::vector<TVector3> s1Vec;
  s1Vec.push_back(vs1_ULB);
  s1Vec.push_back(vs1_ULT);
  s1Vec.push_back(vs1_UTL);
  s1Vec.push_back(vs1_UBL);
  s1Vec.push_back(vs1_UTR);
  s1Vec.push_back(vs1_URT);
  s1Vec.push_back(vs1_URB);
  s1Vec.push_back(vs1_UBR);
  
  std::vector<TVector3> s2Vec;
  s2Vec.push_back(vs2TopLeft);
  s2Vec.push_back(vs2_activeBL);
  s2Vec.push_back(vs2_activeBR);
  s2Vec.push_back(vs2TopRight);
  s2Vec.push_back(vs2Bottom);
  
  std::vector<TVector3> s3Vec;
  s3Vec.push_back(vs3TopLeft);
  s3Vec.push_back(vs3TopRight);
  s3Vec.push_back(vs3BottomLeft);
  s3Vec.push_back(vs3BottomRight);
  
  std::vector<TVector3> s4Vec;
  s4Vec.push_back(vs4_UTL);
  s4Vec.push_back(vs4_UBL);
  s4Vec.push_back(vs4_UTR);
  s4Vec.push_back(vs4_UBR);
  /*
  s4Vec.push_back(vs4_U1R);
  s4Vec.push_back(vs4_U1L);
  s4Vec.push_back(vs4_U2R);
  s4Vec.push_back(vs4_U2L);
  s4Vec.push_back(vs4_U5R);
  s4Vec.push_back(vs4_U5L);
  s4Vec.push_back(vs4_D1L);
  s4Vec.push_back(vs4_D1R);
  s4Vec.push_back(vs4_D5L);
  s4Vec.push_back(vs4_D5R);
  */
  std::vector<TVector3> s4ActiveVec;
  s4ActiveVec.push_back(vs4_activeTL);
  s4ActiveVec.push_back(vs4_activeTR);
  s4ActiveVec.push_back(vs4_activeBL);
  s4ActiveVec.push_back(vs4_activeBR);
  
  TMultiGraph *mg_noProj = new TMultiGraph();
  TGraph *grs1_noProj = new TGraph();
  grs1_noProj->SetMarkerColor(kBlack);
  grs1_noProj->SetMarkerSize(2);
  grs1_noProj->SetMarkerStyle(49);
  TGraph *grs2_noProj = new TGraph();
  grs2_noProj->SetMarkerColor(kOrange);
  grs2_noProj->SetMarkerSize(2);
  grs2_noProj->SetMarkerStyle(49);
  TGraph *grs3_noProj = new TGraph();
  grs3_noProj->SetMarkerColor(kGreen);
  grs3_noProj->SetMarkerSize(2);
  grs3_noProj->SetMarkerStyle(49);
  TGraph *grs4_noProj = new TGraph();
  grs4_noProj->SetMarkerColor(kRed);
  grs4_noProj->SetMarkerSize(2);
  grs4_noProj->SetMarkerStyle(49);
  TGraph *grs4Act_noProj = new TGraph();
  grs4Act_noProj->SetMarkerColor(kRed+2);
  grs4Act_noProj->SetMarkerSize(2);
  grs4Act_noProj->SetMarkerStyle(49);

  double vs1X = 0.;
  double vs1Y = 0.;
  double vs1Z = 0.;
  
  for (int i=0; i<s1Vec.size(); i++) {
    grs1_noProj->SetPoint(grs1_noProj->GetN(), s1Vec[i].Y()*-1, s1Vec[i].Z());
    vs1X += s1Vec[i].X();
    vs1Y += s1Vec[i].Y();
    vs1Z += s1Vec[i].Z();
  }
  vs1X *= (1. / (double)s1Vec.size());
  vs1Y *= (1. / (double)s1Vec.size());
  vs1Z *= (1. / (double)s1Vec.size());
  cout<<vs1X<<" "<<vs1Y<<" "<<vs1Z<<endl;
  TVector3 vs1_avg(vs1X, vs1Y, vs1Z);

  // Want to project these objects down the beamline
  // Have S2 on the drawing plane, i.e. at the origin
  TVector3 vs2Centre(0., 0., 0.);
  vs2Centre = (vs2TopLeft + vs2TopRight + vs2Bottom) * (1. / 3.);
  cout<<"S2 centre "<<(vs2Centre-vs1_avg).X()<<" "<<(vs2Centre-vs1_avg).Y()*-1<<" "<<(vs2Centre-vs1_avg).Z()<<endl;

  //  s2Vec.push_back(vs2Centre);
  // Move S1 to the origin and rotate the axis in xy plane so it points through the centre of S2
  for (int i=0; i<s2Vec.size(); i++) {
    grs2_noProj->SetPoint(grs2_noProj->GetN(), (s2Vec[i].Y()-vs1_avg.Y())*-1, s2Vec[i].Z()-vs1_avg.Z());
    s2Vec[i] -= vs1_avg;
    s2Vec[i].SetY(s2Vec[i].Y()*-1.);
  }
  grs2_noProj->SetPoint(grs2_noProj->GetN(), (vs2Centre.Y()-vs1_avg.Y())*-1., vs2Centre.Z());
  vs2Centre -= vs1_avg;
  vs2Centre.SetY(vs2Centre.Y()*-1.);
  for (int i=0; i<s3Vec.size(); i++) {
    grs3_noProj->SetPoint(grs3_noProj->GetN(), (s3Vec[i].Y()-vs1_avg.Y())*-1, s3Vec[i].Z()-vs1_avg.Z());
    s3Vec[i] -= vs1_avg;
    s3Vec[i].SetY(s3Vec[i].Y()*-1.);
  }
  for (int i=0; i<s4Vec.size(); i++) {
    grs4_noProj->SetPoint(grs4_noProj->GetN(), (s4Vec[i].Y()-vs1_avg.Y())*-1, s4Vec[i].Z()-vs1_avg.Z());
    s4Vec[i] -= vs1_avg;
    s4Vec[i].SetY(s4Vec[i].Y()*-1.);
  }
  for (int i=0; i<s4ActiveVec.size(); i++) {
    grs4Act_noProj->SetPoint(grs4Act_noProj->GetN(), (s4ActiveVec[i].Y()-vs1_avg.Y())*-1, s4ActiveVec[i].Z()-vs1_avg.Z());
    s4ActiveVec[i] -= vs1_avg;
    s4ActiveVec[i].SetY(s4ActiveVec[i].Y()*-1.);
  }
  mg_noProj->Add(grs1_noProj);
  mg_noProj->Add(grs2_noProj);
  mg_noProj->Add(grs3_noProj);
  mg_noProj->Add(grs4_noProj);
  mg_noProj->Add(grs4Act_noProj);
  mg_noProj->SetTitle("Nominal 'Beam's-eye view' of T10 area; -y / m; z / m");
  
  TCanvas *cnoProj = new TCanvas("cnoProj");
  TLegend *legnoProj = new TLegend(0.1, 0.5, 0.27, 0.77);
  legnoProj->AddEntry(grs1_noProj, "S1", "p");
  legnoProj->AddEntry(grs2_noProj, "S2", "p");
  legnoProj->AddEntry(grs3_noProj, "S3", "p");
  legnoProj->AddEntry(grs4_noProj, "S4", "p");
  legnoProj->AddEntry(grs4Act_noProj, "S4 Active", "p");
  mg_noProj->Draw("AP");
  legnoProj->Draw();
  cnoProj->Print(Form("%s/beamlineOrig.png", saveDir));
  cnoProj->Print(Form("%s/beamlineOrig.pdf", saveDir));   

  // To find angle to rotate
  // Vector pointing down beamline
  TVector3 vXY(1., 0., 0.);
  TVector3 vs2Angle(vs2Centre.X(), vs2Centre.Y(), 0.);
  double s2Angle = vXY.Angle(vs2Angle);
  cout<<"Angle to S2 from S1: "<<s2Angle<<" "<<s2Angle*(180./TMath::Pi())<<endl;

  TMultiGraph *mg_rot = new TMultiGraph();
  TGraph *grs1_rot = new TGraph();
  grs1_rot->SetMarkerColor(kBlack);
  grs1_rot->SetMarkerSize(2);
  grs1_rot->SetMarkerStyle(49);
  TGraph *grs2_rot = new TGraph();
  grs2_rot->SetMarkerColor(kOrange);
  grs2_rot->SetMarkerSize(2);
  grs2_rot->SetMarkerStyle(49);
  TGraph *grs3_rot = new TGraph();
  grs3_rot->SetMarkerColor(kGreen);
  grs3_rot->SetMarkerSize(2);
  grs3_rot->SetMarkerStyle(49);
  TGraph *grs4_rot = new TGraph();
  grs4_rot->SetMarkerColor(kRed);
  grs4_rot->SetMarkerSize(2);
  grs4_rot->SetMarkerStyle(49);
  TGraph *grs4Active_rot = new TGraph();
  grs4Active_rot->SetMarkerColor(kRed+2);
  grs4Active_rot->SetMarkerSize(2);
  grs4Active_rot->SetMarkerStyle(49);
  // Rotate about z axis by this angle
  for (int i=0; i<s1Vec.size(); i++) {
    s1Vec[i].RotateZ(-s2Angle);
  }
  for (int i=0; i<s2Vec.size(); i++) {
    s2Vec[i].RotateZ(-s2Angle);
    grs2_rot->SetPoint(grs2_rot->GetN(), s2Vec[i].Y(), s2Vec[i].Z());
  }
  vs2Centre.RotateZ(-s2Angle);
  grs2_rot->SetPoint(grs2_rot->GetN(), vs2Centre.Y(), vs2Centre.Z());
  for(int i=0; i<s3Vec.size(); i++) {
    s3Vec[i].RotateZ(-s2Angle);
    grs3_rot->SetPoint(grs3_rot->GetN(), s3Vec[i].Y(), s3Vec[i].Z());
  }
  for (int i=0; i<s4Vec.size(); i++) {
    s4Vec[i].RotateZ(-s2Angle);
    grs4_rot->SetPoint(grs4_rot->GetN(), s4Vec[i].Y(), s4Vec[i].Z());
  }
  for (int i=0; i<s4ActiveVec.size(); i++) {
    s4ActiveVec[i].RotateZ(-s2Angle);
    grs4Active_rot->SetPoint(grs4Active_rot->GetN(), s4ActiveVec[i].Y(), s4ActiveVec[i].Z());
  }
  mg_rot->Add(grs2_rot);
  mg_rot->Add(grs3_rot);
  mg_rot->Add(grs4_rot);
  mg_rot->Add(grs4Active_rot);
  mg_rot->SetTitle("T10 points rotated to S1-S2 axis; -y / m; z / m");
  TLegend *legRot = new TLegend(0.1, 0.5, 0.29, 0.75);
  legRot->AddEntry(grs2_noProj, "S2", "p");
  legRot->AddEntry(grs3_noProj, "S3", "p");
  legRot->AddEntry(grs4_noProj, "S4", "p");
  legRot->AddEntry(grs4Act_noProj, "S4 Active", "p");
  TCanvas *crot = new TCanvas("crot");
  mg_rot->Draw("AP");
  legRot->Draw();
  crot->Print(Form("%s/beamlineRotated.png",saveDir));
  crot->Print(Form("%s/beamlineRotated.pdf",saveDir));
  // Now project down the z axis with the paddle in the drawing plane
  // Origin still at S1
  TMultiGraph *mg_proj = new TMultiGraph();
  TGraph *grs1_proj = new TGraph();
  grs1_proj->SetMarkerColor(kBlack);
  grs1_proj->SetMarkerSize(2);
  grs1_proj->SetMarkerStyle(49);
  TGraph *grs2_proj = new TGraph();
  grs2_proj->SetMarkerColor(kOrange);
  grs2_proj->SetMarkerSize(2);
  grs2_proj->SetMarkerStyle(49);
  TGraph *grs3_proj = new TGraph();
  grs3_proj->SetMarkerColor(kGreen);
  grs3_proj->SetMarkerSize(2);
  grs3_proj->SetMarkerStyle(49);
  TGraph *grs4_proj = new TGraph();
  grs4_proj->SetMarkerColor(kRed);
  grs4_proj->SetMarkerSize(2);
  grs4_proj->SetMarkerStyle(49);
  TGraph *grs4Act_proj = new TGraph();
  grs4Act_proj->SetMarkerColor(kRed+2);
  grs4Act_proj->SetMarkerSize(2);
  grs4Act_proj->SetMarkerStyle(49);
  
  const double d = vs2Centre.X();
  for (int i=0; i<s2Vec.size(); i++) {
    grs2_proj->SetPoint(grs2_proj->GetN(), (s2Vec[i].Y()*d)/s2Vec[i].X(), (s2Vec[i].Z()*d)/s2Vec[i].X());
  }
  grs2_proj->SetPoint(grs2_proj->GetN(), (vs2Centre.Y()*d)/vs2Centre.X(), (vs2Centre.Z()*d)/vs2Centre.X());
  for (int i=0; i<s3Vec.size(); i++) {
    grs3_proj->SetPoint(grs3_proj->GetN(), (s3Vec[i].Y()*d)/s3Vec[i].X(), (s3Vec[i].Z()*d)/s3Vec[i].X());
  }
  for (int i=0; i<s4Vec.size(); i++) {
    grs4_proj->SetPoint(grs4_proj->GetN(), (s4Vec[i].Y()*d)/s4Vec[i].X(), (s4Vec[i].Z()*d)/s4Vec[i].X());
  }
  for (int i=0; i<s4ActiveVec.size(); i++) {
    grs4Act_proj->SetPoint(grs4Act_proj->GetN(), (s4ActiveVec[i].Y()*d)/s4ActiveVec[i].X(), (s4ActiveVec[i].Z()*d)/s4ActiveVec[i].X());
  }
  mg_proj->Add(grs2_proj);
  mg_proj->Add(grs3_proj);
  mg_proj->Add(grs4_proj);
  mg_proj->Add(grs4Act_proj);
  TCanvas *cproj = new TCanvas("cproj");
  mg_proj->SetTitle("T10 points projected along S1-S2 axis (S1 origin); -y / m; z / m");
  mg_proj->Draw("AP");
  TLegend *legProj = new TLegend(0.1, 0.15, 0.30, 0.35);
  legProj->AddEntry(grs2_noProj, "S2", "p");
  legProj->AddEntry(grs3_noProj, "S3", "p");
  legProj->AddEntry(grs4_noProj, "S4", "p");
  legProj->AddEntry(grs4Act_noProj, "S4 Active", "p");
  legProj->Draw();
  cproj->Print(Form("%s/beamlineProjected.png",saveDir));
  cproj->Print(Form("%s/beamlineProjected.pdf",saveDir));

  // Project but without draw distance to S2
  TMultiGraph *mg_nproh = new TMultiGraph();
  TGraph *grs1_nproh = new TGraph();
  grs1_nproh->SetMarkerColor(kBlack);
  grs1_nproh->SetMarkerSize(2);
  grs1_nproh->SetMarkerStyle(49);
  TGraph *grs2_nproh = new TGraph();
  grs2_nproh->SetMarkerColor(kOrange);
  grs2_nproh->SetMarkerSize(2);
  grs2_nproh->SetMarkerStyle(49);
  TGraph *grs3_nproh = new TGraph();
  grs3_nproh->SetMarkerColor(kGreen);
  grs3_nproh->SetMarkerSize(2);
  grs3_nproh->SetMarkerStyle(49);
  TGraph *grs4_nproh = new TGraph();
  grs4_nproh->SetMarkerColor(kRed);
  grs4_nproh->SetMarkerSize(2);
  grs4_nproh->SetMarkerStyle(49);
  TGraph *grs4Act_nproh = new TGraph();
  grs4Act_nproh->SetMarkerColor(kRed+2);
  grs4Act_nproh->SetMarkerSize(2);
  grs4Act_nproh->SetMarkerStyle(49);
  
  for (int i=0; i<s2Vec.size(); i++) {
    grs2_nproh->SetPoint(grs2_nproh->GetN(), s2Vec[i].Y()/s2Vec[i].X(), s2Vec[i].Z()/s2Vec[i].X());
  }
  grs2_nproh->SetPoint(grs2_nproh->GetN(), vs2Centre.Y()/vs2Centre.X(), vs2Centre.Z()/vs2Centre.X());
  for (int i=0; i<s3Vec.size(); i++) {
    grs3_nproh->SetPoint(grs3_nproh->GetN(), s3Vec[i].Y()/s3Vec[i].X(), s3Vec[i].Z()/s3Vec[i].X());
  }
  for (int i=0; i<s4Vec.size(); i++) {
    grs4_nproh->SetPoint(grs4_nproh->GetN(), s4Vec[i].Y()/s4Vec[i].X(), s4Vec[i].Z()/s4Vec[i].X());
  }
  for (int i=0; i<s4ActiveVec.size(); i++) {
    grs4Act_nproh->SetPoint(grs4Act_nproh->GetN(), s4ActiveVec[i].Y()/s4ActiveVec[i].X(), s4ActiveVec[i].Z()/s4ActiveVec[i].X());
  }
  mg_nproh->Add(grs2_nproh);
  mg_nproh->Add(grs3_nproh);
  mg_nproh->Add(grs4_nproh);
  mg_nproh->Add(grs4Act_nproh);
  TCanvas *cnproh = new TCanvas("cnproh");
  mg_nproh->SetTitle("T10 points projected along S1-S2 axis (S1 origin); -y / m; z / m");
  mg_nproh->Draw("AP");
  legProj->Draw();
  cnproh->Print(Form("%s/beamlineProjectedNew.png",saveDir));
  cnproh->Print(Form("%s/beamlineProjectedNew.pdf",saveDir));

  // Change the axes so we're doing this in terms of angle
  TMultiGraph *mg_ang = new TMultiGraph();
  TGraph *grs1_ang = new TGraph();
  grs1_ang->SetMarkerColor(kBlack);
  grs1_ang->SetMarkerSize(2);
  grs1_ang->SetMarkerStyle(49);
  TGraph *grs2_ang = new TGraph();
  grs2_ang->SetMarkerColor(kOrange);
  grs2_ang->SetMarkerSize(2);
  grs2_ang->SetMarkerStyle(49);
  TGraph *grs3_ang = new TGraph();
  grs3_ang->SetMarkerColor(kGreen);
  grs3_ang->SetMarkerSize(2);
  grs3_ang->SetMarkerStyle(49);
  TGraph *grs4_ang = new TGraph();
  grs4_ang->SetMarkerColor(kRed);
  grs4_ang->SetMarkerSize(2);
  grs4_ang->SetMarkerStyle(49);
  TGraph *grs4Act_ang = new TGraph();
  grs4Act_ang->SetMarkerColor(kRed+2);
  grs4Act_ang->SetMarkerSize(2);
  grs4Act_ang->SetMarkerStyle(49);
  
  for (int i=0; i<s2Vec.size(); i++) {
    grs2_ang->SetPoint(grs2_ang->GetN(), (TMath::ATan(s2Vec[i].Y()/s2Vec[i].X())+s2Angle)*(180./TMath::Pi()), TMath::ATan(s2Vec[i].Z()/s2Vec[i].X())*(180./TMath::Pi()));
  }
  grs2_ang->SetPoint(grs2_ang->GetN(), (TMath::ATan(vs2Centre.Y()/vs2Centre.X())+s2Angle)*(180./TMath::Pi()), TMath::ATan(vs2Centre.Z()/vs2Centre.X())*(180./TMath::Pi()));
  for (int i=0; i<s3Vec.size(); i++) {
    grs3_ang->SetPoint(grs3_ang->GetN(), (TMath::ATan(s3Vec[i].Y()/s3Vec[i].X())+s2Angle)*(180./TMath::Pi()), TMath::ATan(s3Vec[i].Z()/s3Vec[i].X())*(180./TMath::Pi()));
  }
  for (int i=0; i<s4Vec.size(); i++) {
    grs4_ang->SetPoint(grs4_ang->GetN(), (TMath::ATan(s4Vec[i].Y()/s4Vec[i].X())+s2Angle)*(180./TMath::Pi()), TMath::ATan(s4Vec[i].Z()/s4Vec[i].X())*(180./TMath::Pi()));
  }
  for (int i=0; i<s4ActiveVec.size(); i++) {
    grs4Act_ang->SetPoint(grs4Act_ang->GetN(), (TMath::ATan(s4ActiveVec[i].Y()/s4ActiveVec[i].X())+s2Angle)*(180./TMath::Pi()), TMath::ATan(s4ActiveVec[i].Z()/s4ActiveVec[i].X())*(180./TMath::Pi()) );
  }
  mg_ang->Add(grs2_ang);
  mg_ang->Add(grs3_ang);
  mg_ang->Add(grs4_ang);
  mg_ang->Add(grs4Act_ang);
  TCanvas *cang = new TCanvas("cang");
  mg_ang->SetTitle("Angular distribution of T10 (S1 origin); #theta / degrees; #phi / degrees");
  mg_ang->Draw("AP");
  legProj->Draw();
  cang->Print(Form("%s/beamlineAng.png",saveDir));
  cang->Print(Form("%s/beamlineAng.pdf",saveDir));

  cout<<"S4 angular positions"<<endl;
  for (int i=0; i<s4ActiveVec.size(); i++) {
    cout<<"Theta: "<<(TMath::ATan(s4ActiveVec[i].Y()/s4ActiveVec[i].X())+s2Angle)*(180./TMath::Pi())<<", Phi: "<<TMath::ATan(s4ActiveVec[i].Z()/s4ActiveVec[i].X())*(180./TMath::Pi())<<endl;
  }
  cout<<"S3 angular positions"<<endl;
  for (int i=0; i<s3Vec.size(); i++) {
    cout<<"Theta: "<<(TMath::ATan(s3Vec[i].Y()/s3Vec[i].X())+s2Angle)*(180./TMath::Pi())<<", Phi: "<<TMath::ATan(s3Vec[i].Z()/s3Vec[i].X())*(180./TMath::Pi())<<endl;
  }
  cout<<"S2 angular positions"<<endl;
  for (int i=0; i<s2Vec.size(); i++) {
    cout<<"Theta: "<<(TMath::ATan(s2Vec[i].Y()/s2Vec[i].X())+s2Angle)*(180./TMath::Pi())<<", Phi: "<<TMath::ATan(s2Vec[i].Z()/s2Vec[i].X())*(180./TMath::Pi())<<endl;
  }

  std::vector<TVector3> vectorVec;
  vectorVec.push_back(vs2TopLeft);
  vectorVec.push_back(vs2TopRight);
  vectorVec.push_back(vs2Bottom);
  vectorVec.push_back(vs3TopLeft);
  vectorVec.push_back(vs3TopRight);
  vectorVec.push_back(vs3BottomLeft);
  vectorVec.push_back(vs3BottomRight);
  vectorVec.push_back(vs4_UTL);
  vectorVec.push_back(vs4_UBL);
  vectorVec.push_back(vs4_UTR);
  vectorVec.push_back(vs4_UBR);
  vectorVec.push_back(vs4_U1R);
  vectorVec.push_back(vs4_U1L);
  vectorVec.push_back(vs4_D5L);
  vectorVec.push_back(vs4_D5R);

} // overlapCalc
