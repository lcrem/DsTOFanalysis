// overlapCalc.C

void overlapCalc()
{
  gSystem->Load("libPhysics.so");
  
  TVector3 vs1Centre(-1.765, 0.036, -0.002);
  
  TVector3 vs2TopLeft(-0.3472, 0.0344, 0.0706);
  TVector3 vs2TopRight(-0.3498, -0.0725, 0.0679);
  TVector3 vs2Bottom(-0.3503, -0.0133, -0.2928);

  TVector3 vs3TopLeft(9.0650, 0.5215, 0.6244);
  TVector3 vs3TopRight(8.9245, -0.9928, 0.6220);
  TVector3 vs3BottomLeft(9.0488, 0.5120, -0.5993);
  TVector3 vs3BottomRight(8.9047, -1.0012, -0.6012);

  /*
  TVector3 vs4TopLeft(
  TVector3 vs4TopRight(
  TVector3 vs4BottomLeft(
  TVector3 vs4BottomRight(
  */

  std::vector<TVector3> vectorVec;
  vectorVec.push_back(vs2TopLeft);
  vectorVec.push_back(vs2TopRight);
  vectorVec.push_back(vs2Bottom);
  vectorVec.push_back(vs3TopLeft);
  vectorVec.push_back(vs3TopRight);
  vectorVec.push_back(vs3BottomLeft);
  vectorVec.push_back(vs3BottomRight);
  // Want to project these objects down the beamline
  // Have S2 on the drawing plane, ie at the origin
  TVector3 vs2Centre(0., 0., 0.);
  vs2Centre = (vs2TopLeft + vs2TopRight + vs2Bottom) * (1. / 3.);
  cout<<vs2Centre.X()<<" "<<vs2Centre.Y()<<" "<<vs2Centre.Z()<<endl;

  // We want S1 to be at (-(x distance from S1 to S2), 0, 0)
  // Vector to correct for this
  TVector3 correctVec(vs2TopLeft.X()*-1. , vs1Centre.Y()*-1., vs1Centre.Z()*-1.);

  cout<<"Now S1 centre is at "<<(vs1Centre.X()+correctVec.X())<<" "<<(vs1Centre.Y()+correctVec.Y())<<" "<<(vs1Centre.Z()+correctVec.Z())<<endl;
  // Now apply this correction to all the other vectors
  for (int i=0; i<vectorVec.size(); i++) {
    vectorVec[i] += correctVec;
  }

  TMultiGraph *mg = new TMultiGraph();
  TGraph *grs2 = new TGraph();
  grs2->SetMarkerColor(kOrange);
  grs2->SetMarkerSize(2);
  TGraph *grs3 = new TGraph();
  grs3->SetMarkerColor(kGreen);
  grs3->SetMarkerSize(2);
  
  const double d = vs1Centre.X() * -1.;
  // y' = y*d/(x+d)
  grs2->SetPoint(grs2->GetN(), (vs2TopLeft.Y()*d)/(vs2TopLeft.X()+d),  (vs2TopLeft.Z()*d)/(vs2TopLeft.X()+d));
  grs2->SetPoint(grs2->GetN(), (vs2TopRight.Y()*d)/(vs2TopRight.X()+d), (vs2TopRight.Z()*d)/(vs2TopRight.X()+d));
  grs2->SetPoint(grs2->GetN(), (vs2Bottom.Y()*d)/(vs2Bottom.X()+d),   (vs2Bottom.Z()*d)/(vs2Bottom.X()+d));

  grs3->SetPoint(grs3->GetN(), (vs3TopLeft.Y()*d)/(vs3TopLeft.X()+d),  (vs3TopLeft.Z()*d)/(vs3TopLeft.X()+d));
  grs3->SetPoint(grs3->GetN(), (vs3TopRight.Y()*d)/(vs3TopRight.X()+d),  (vs3TopRight.Z()*d)/(vs3TopRight.X()+d));
  grs3->SetPoint(grs3->GetN(), (vs3BottomLeft.Y()*d)/(vs3BottomLeft.X()+d),  (vs3BottomLeft.Z()*d)/(vs3BottomLeft.X()+d));
  grs3->SetPoint(grs3->GetN(), (vs3BottomRight.Y()*d)/(vs3BottomRight.X()+d),  (vs3BottomRight.Z()*d)/(vs3BottomRight.X()+d));
  
  mg->Add(grs2);
  mg->Add(grs3);
  mg->Draw("AP*");
} // overlapCalc
