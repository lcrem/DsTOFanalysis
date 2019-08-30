// distancesAndErrors.C
// Calculate all of these upstream, downstream and off-axis angles
// Calculate some kind of error band on these measurements
// Argument is output txt file so these can be easily put in a LaTeX doc
void distancesAndErrors(const char *saveFile) {
  
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
  //TVector3 vs1_C = (vs1_ULB+vs1_ULT+vs1_UTL+vs1_UBL+vs1_UTR+vs1_URT+vs1_URB+vs1_UBR) * 0.125;
  TVector3 vs1_C = (vs1TopLeft + vs1TopRight + vs1BottomLeft + vs1BottomRight)*0.25; 
    
  TVector3 vs2TopLeft(-0.3472, 0.0344, 0.0706);
  TVector3 vs2_activeBL(-0.3472, 0.0344, -0.0496);
  TVector3 vs2TopRight(-0.3498, -0.0725, 0.0679);
  TVector3 vs2_activeBR(-0.3498, -0.0725, -0.0521);
  TVector3 vs2Bottom(-0.3503, -0.0133, -0.2928);
  TVector3 vs2_C = (vs2TopLeft+vs2TopRight+vs2_activeBL+vs2_activeBR)*0.25;

  TVector3 D3_UTR(8.9245, -0.9928, 0.6220);
  TVector3 D3_UBR(8.9047, -1.0012, -0.6012);
  TVector3 D3_UBL(9.0488, 0.5120, -0.5993);
  TVector3 D3_UTL(9.0650, 0.5215, 0.6244);
  double l1 = TMath::Sqrt(pow(D3_UTR.X()-D3_UTL.X(), 2) + pow(D3_UTR.Y()-D3_UTL.Y(), 2));
  double l2 = TMath::Sqrt(pow(D3_UBR.X()-D3_UBL.X(), 2) + pow(D3_UBR.Y()-D3_UBL.Y(), 2));
  double grad = (D3_UTR.Y()-D3_UTL.Y())/(D3_UTR.X()-D3_UTL.X());
  double xExt = 0.08/TMath::Sqrt(1+pow(grad,2));
  double yExt = xExt * grad;

  TVector3 s3Ext(xExt, yExt, 0.);
  TVector3 vs3TopLeft = D3_UTL + s3Ext;
  TVector3 vs3TopRight = D3_UTR - s3Ext;
  TVector3 vs3BottomLeft = D3_UBL + s3Ext;
  TVector3 vs3BottomRight = D3_UBR - s3Ext;

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
  // Vectors for camera plane
  // 0.67 m to right of cathode
  TVector3 cameraTop(10.8059, -1.4549, 0.5386);
  TVector3 cameraBottom(10.8059, -1.4549, -0.5614);
  std::vector<TVector3> cameraVec = {cameraTop, cameraBottom};
  
  std::vector<TVector3> vesselVec = {vesselTopLeft, vesselTopRight, 
				     vesselBottomRight, vesselBottomLeft};
  std::vector<TVector3> tpcUsVec = {TPCActiveTopUL, TPCActiveTopUR,
				    TPCActiveBottomUR, TPCActiveBottomUL};
  std::vector<TVector3> tpcDsVec = {TPCActiveTopDL, TPCActiveTopDR,
-				    TPCActiveBottomDR, TPCActiveBottomDL};
  std::vector<TVector3> tpcVec = {TPCActiveTopUL, TPCActiveTopUR,
				  TPCActiveBottomUR, TPCActiveBottomUL,
				  TPCActiveTopDL, TPCActiveTopDR,
				  TPCActiveBottomDR, TPCActiveBottomDL};
  //std::vector<TVector3> s1Vec = {vs1_ULB,vs1_ULT,vs1_UTL,vs1_UBL,vs1_UTR,vs1_URT,vs1_URB,vs1_UBR};
  std::vector<TVector3> s1Vec = {vs1TopLeft, vs1TopRight, vs1BottomLeft, vs1BottomRight};
  std::vector<TVector3> s2Vec = {vs2TopLeft, vs2_activeBL, vs2_activeBR, vs2TopRight};  
  std::vector<TVector3> s3Vec = {vs3TopLeft, vs3TopRight, vs3BottomRight, vs3BottomLeft};
  std::vector<TVector3> s4Vec = {vs4_UTL, vs4_UBL, vs4_UTR, vs4_UBR};
  std::vector<TVector3> s4ActiveVec = {vs4_activeTL, vs4_activeTR, vs4_activeBR, vs4_activeBL,};

  std::ofstream outfile;
  outfile.open(saveFile);
  // Do distances first
  TVector3 vs3_C = (vs3TopLeft + vs3TopRight + vs3BottomLeft + vs3BottomRight) * 0.25;
  TVector3 vs4_C = (vs4_activeTL + vs4_activeTR + vs4_activeBR + vs4_activeBL) * 0.25;
  TVector3 tpc_C = (TPCActiveTopUL+TPCActiveTopUR+TPCActiveBottomUR+TPCActiveBottomUL+
		    TPCActiveTopDL+TPCActiveTopDR+TPCActiveBottomDR+TPCActiveBottomDL)*0.125;
  double wcs1Dist = (WC_C - vs1_C).Mag();
  double s1s2Dist = (vs2_C - vs1_C).Mag();
  double s1s3Dist = (vs3_C - vs1_C).Mag();
  double s1s4Dist = (vs4_C - vs1_C).Mag();
  double s2s4Dist = (vs4_C - vs2_C).Mag();
  double s3tpcDist = (tpc_C - vs3_C).Mag();
  double s4tpcDist = (vs4_C - tpc_C).Mag();
  // Now calculate errors
  // For each of the measured points calculate the distance to every other measured point
  double highests1s2 = 0;
  double lowests1s2 = 999999;
  for (int i = 0; i < s1Vec.size(); i++) {
    for (int j = 0; j < s2Vec.size(); j++) {
      double dist = (s2Vec.at(j) - s1Vec.at(i)).Mag();
      if (dist > highests1s2) highests1s2 = dist;
      if (dist < lowests1s2) lowests1s2 = dist;
    }
  }

  cout<<"Point to point: "<<s1s2Dist<<", Lowest S1S2 "<<lowests1s2-s1s2Dist<<", Highest S1S2 "<<highests1s2-s1s2Dist<<endl;
  
  outfile<<"\\hline"<<endl;
  outfile<<"Distance between various objects in beamline\\\ "<<endl;
  outfile<<"\\hline"<<endl;
  outfile<<"$\\mathit{S1}-\\mathit{S2}$ & $\\mathit{S1}-\\mathit{S3}$ & $\\mathit{S3}$ - TPC drift volume & TPC drift volume - $\\mathit{S4}$ & $\\mathit{S2}-\\mathit{S4}$  \\\ "<<endl;
  outfile<<"$("<<s1s2Dist<<"\\pm)~\\text{m}$ & "<<"$("<<s1s3Dist<<"\\pm)~\\text{m}$ & "<<"$("<<s3tpcDist<<"\\pm)~\\text{m}$ & "<<"$("<<s4tpcDist<<"\\pm)~\\text{m}$ & "<<"$("<<s2s4Dist<<"\\pm)~\\text{m}$"<<endl;
  outfile.close();
}

double atanErr(const double opp, const double adj, const double inerr) {
  double err = TMath::Sqrt(pow(adj/(pow(adj, 2)+pow(opp, 2)), 2) * pow(inerr, 2) +
			   pow(opp/(pow(adj, 2)+pow(opp, 2)), 2) * pow(inerr, 2));
  return err;
}
