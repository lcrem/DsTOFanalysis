// mcGeomComp.C
// Geometry check for the MC. Want to make sure that everything is where I expect it to be
#include "UsefulFunctions.C"

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

const double phimin = -4.;
const double phimax = 4.;
const double thetamin = -7.;
const double thetamax = 3.5;
const int phibins = 26;
const int thetabins = 70;

TVector3 TPCActiveBUL(-0.2778, -0.5614, 10.3048);
TVector3 TPCActiveBDL(-0.3896, -0.5614, 11.3991);
TVector3 TPCActiveTUL(-0.2778, 0.5386, 10.3048);
TVector3 TPCActiveTDL(-0.3896, 0.5386, 11.3991);
TVector3 TPCActiveBUR(-0.7290, -0.5614, 10.2587);
TVector3 TPCActiveBDR(-0.8408, -0.5614, 11.3530);
TVector3 TPCActiveTUR(-0.7290, 0.5386, 10.2587);
TVector3 TPCActiveTDR(-0.8408, 0.5386, 11.3530);
TVector3 TPCActiveC(-0.5593, -0.0114, 10.82885);

void mcGeomComp(const char* outfile)
{
  gROOT->SetBatch(kTRUE);

  const char* mcDir = "/scratch0/tnonnenm/FixedMC/Prototype-HPTPC-MC/wdir/moderator_data/";

  const TVector3 vS1 = (D1_ULB+D1_ULT+D1_URB+D1_URT)*0.25;

  TFile *fout = new TFile(outfile, "recreate");

  // Loop over blocks
  for (int nb=0; nb<5; nb++) {
    cout<<"=========================================="<<endl;
    cout<<nb<<" blocks"<<endl;
    cout<<"=========================================="<<endl;

    TFile *fin = new TFile(Form("%s/NewGeom%dBlocks.root", mcDir, nb), "read");
    TTree *stepTree = (TTree*)fin->Get("stepTree");
    UInt_t NSteps;
    float PostStepMom[40000];
    float PostStepPosX[40000];
    float PostStepPosY[40000];
    float PostStepPosZ[40000];
    float PreStepPosX[40000];
    float PreStepPosY[40000];
    float PreStepPosZ[40000];
    int PreStepVolumeID[40000];
    int PostStepVolumeID[40000];
    int PDG[40000];
    stepTree->SetBranchAddress("NSteps", &NSteps);
    stepTree->SetBranchAddress("PostStepMom", PostStepMom);
    stepTree->SetBranchAddress("PostStepPosX", PostStepPosX);
    stepTree->SetBranchAddress("PostStepPosY", PostStepPosY);
    stepTree->SetBranchAddress("PostStepPosZ", PostStepPosZ);
    stepTree->SetBranchAddress("PreStepPosX", PreStepPosX);
    stepTree->SetBranchAddress("PreStepPosY", PreStepPosY);
    stepTree->SetBranchAddress("PreStepPosZ", PreStepPosZ);
    stepTree->SetBranchAddress("PreStepVolumeID", PreStepVolumeID);
    stepTree->SetBranchAddress("PostStepVolumeID", PostStepVolumeID);
    stepTree->SetBranchAddress("PDG", PDG);

    // Histograms
    // S3
    TH2D *hS3 = new TH2D(Form("hS3_%d", nb), Form("MC hits in S3, %d blocks; #theta / degrees; #phi / degrees; Events / spill", nb), thetabins, thetamin, thetamax, phibins, phimin, phimax);
    setHistAttr(hS3);
    TH2D *hS3withS2 = new TH2D(Form("hS3withS2_%d", nb), Form("MC hits in S3, %d blocks (S2 trigger); #theta / degrees; #phi / degrees; Events / spill", nb), thetabins, thetamin, thetamax, phibins, phimin, phimax);
    setHistAttr(hS3withS2);
    TH2D *hS3noS2 = new TH2D(Form("hS3noS2_%d", nb), Form("MC hits in S3, %d blocks (no S2 trigger); #theta / degrees; #phi / degrees; Events / spill", nb), thetabins, thetamin, thetamax, phibins, phimin, phimax);
    setHistAttr(hS3noS2);
    // S4
    TH2D *hS4 = new TH2D(Form("hS4_%d", nb), Form("MC hits in S4, %d blocks; #theta / degrees; #phi / degrees; Events / spill", nb), thetabins, thetamin, thetamax, phibins, phimin, phimax);
    setHistAttr(hS4);
    TH2D *hS4withS2 = new TH2D(Form("hS4withS2_%d", nb), Form("MC hits in S4, %d blocks (S2 trigger); #theta / degrees; #phi / degrees; Events / spill", nb), thetabins, thetamin, thetamax, phibins, phimin, phimax);
    setHistAttr(hS4withS2);
    TH2D *hS4noS2 = new TH2D(Form("hS4noS2_%d", nb), Form("MC hits in S4, %d blocks (no S2 trigger); #theta / degrees; #phi / degrees; Events / spill", nb), thetabins, thetamin, thetamax, phibins, phimin, phimax);
    setHistAttr(hS4noS2);
    // Global XY positions
    TH2D *hS4XY = new TH2D(Form("hS4XY_%d", nb), Form("MC hits in S4, %d blocks; X / m; Y / m", nb), 100, 0, -2., 100, -0.6, 0.6);
    setHistAttr(hS4XY);
    TH2D *hS4XZ = new TH2D(Form("hS4XZ_%d", nb), Form("MC hits in S4, %d blocks; X / m; Z / m", nb), 100, 0, -2., 100, 11., 12.5);
    setHistAttr(hS4XZ);
    TH2D *hS4YZ = new TH2D(Form("hS4YZ_%d", nb), Form("MC hits in S4, %d blocks; Z / m; Y / m", nb), 100, 0.6, -0.6, 100, 11., 12.5);
    setHistAttr(hS4YZ);
    // Vessel doors
    TH2D *hDoors = new TH2D(Form("hDoors_%d", nb), Form("MC hits in vessel doors, %d blocks; #theta / degrees; #phi / degrees; Events / spill", nb), thetabins, thetamin, thetamax, phibins, phimin, phimax);
    setHistAttr(hDoors);
    TH2D *hDoorsWithS2 = new TH2D(Form("hDoorsWithS2_%d", nb), Form("MC hits in vessel doors, %d blocks (S2 trigger); #theta / degrees; #phi / degrees; Events / spill", nb), thetabins, thetamin, thetamax, phibins, phimin, phimax);
    setHistAttr(hDoorsWithS2);
    TH2D *hDoorsNoS2 = new TH2D(Form("hDoorsNoS2_%d", nb), Form("MC hits in vessel doors, %d blocks (no S2 trigger); #theta / degrees; #phi / degrees; Events / spill", nb), thetabins, thetamin, thetamax, phibins, phimin, phimax);
    setHistAttr(hDoorsNoS2);
    TH2D *hDoor1XY = new TH2D(Form("hDoor1XY_%d", nb), Form("MC hits in vessel door 1, %d blocks; X / m; Y / m; Events", nb), 100, -1.5, 0.4, 100, -0.8, 0.8);
    setHistAttr(hDoor1XY);
    TH2D *hDoor1XZ = new TH2D(Form("hDoor1XZ_%d", nb), Form("MC hits in vessel door 1, %d blocks; X / m; Z / m; Events", nb), 100, -1.5, 0.4, 100, 10., 11.8);
    setHistAttr(hDoor1XZ);
    TH2D *hDoor1YZ = new TH2D(Form("hDoor1YZ_%d", nb), Form("MC hits in vessel door 1, %d blocks; Y / m; Z / m; Events", nb), 100, -0.8, 0.8, 100, 10., 11.8);
    setHistAttr(hDoor1YZ);
    TH2D *hDoor2XY = new TH2D(Form("hDoor2XY_%d", nb), Form("MC hits in vessel door 2, %d blocks; X / m; Y / m; Events", nb), 100, -1.5, 0.4, 100, -0.8, 0.8);
    setHistAttr(hDoor2XY);
    TH2D *hDoor2XZ = new TH2D(Form("hDoor2XZ_%d", nb), Form("MC hits in vessel door 2, %d blocks; X / m; Z / m; Events", nb), 100, -1.5, 0.4, 100, 10., 11.8);
    setHistAttr(hDoor2XZ);
    TH2D *hDoor2YZ = new TH2D(Form("hDoor2YZ_%d", nb), Form("MC hits in vessel door 2, %d blocks; Y / m; Z / m; Events", nb), 100, -0.8, 0.8, 100, 10., 11.8);
    setHistAttr(hDoor2XZ);
    // TPC active volume
    TH2D *hTpc = new TH2D(Form("hTpc_%d", nb), Form("MC hits in TPC active areas, %d blocks; #theta / degrees; #phi / degrees; Events / spill", nb), thetabins, thetamin, thetamax, phibins, phimin, phimax);
    setHistAttr(hTpc);

    TH2D *hSteelXY = new TH2D(Form("hSteelXY_%d", nb), Form("MC hits in the steel, %d blocks; X / m; Y / m", nb), 100., -1.5, 0.4, 100, -0.8, 0.8);
    setHistAttr(hSteelXY);
    TH2D *hSteelXZ = new TH2D(Form("hSteelXZ_%d", nb), Form("MC hits in the steel, %d blocks; X / m; Z / m", nb), 100., -1.5, 0.4, 100, 10.,11.8);
    setHistAttr(hSteelXZ);
    TH2D *hSteelYZ = new TH2D(Form("hSteelYZ_%d", nb), Form("MC hits in the steel, %d blocks; Y / m; Z / m", nb), 100, -0.8, 0.8, 100, 10., 11.8);
    setHistAttr(hSteelYZ);

    stepTree->AddFriend(Form("protonTree%dBlocks", nb), "/scratch0/sjones/plots/angularDistS3_newSample/s3ProtonWeightTrees.root"); // Tree for weights and spills added as friend    
    double weight, error;
    int spill;
    int isS2;
    stepTree->SetBranchAddress("spill", &spill);
    stepTree->SetBranchAddress("isS2",  &isS2);
    stepTree->SetBranchAddress("weight", &weight);
    stepTree->SetBranchAddress("error",  &error);

    for (int ii=0; ii<stepTree->GetEntries(); ii++) {
      stepTree->GetEntry(ii);
      TVector3 Hit(PreStepPosX[0]/1000., PreStepPosY[0]/1000., PreStepPosZ[0]/1000.);
      TVector3 hitPos = MCToGlobalCoords(Hit);
      hitPos = hitPos - vS1;

      double theta = TMath::ATan(hitPos.X()/hitPos.Z())*(180./TMath::Pi());
      double phi   = TMath::ATan(hitPos.Y()/hitPos.Z())*(180./TMath::Pi());

      hS3->Fill(theta, phi, weight);
      if (isS2) hS3withS2->Fill(theta, phi, weight);
      else hS3noS2->Fill(theta, phi, weight);

      // Loop over steps to get the S4 hits
      for (int n=0; n<NSteps; n++) {
	if (PreStepVolumeID[n]!=7 && PostStepVolumeID[n]==7 && PDG[n]==2212 && PostStepMom[n] > 30.) {
	  TVector3 mcHitS4(PostStepPosX[n]/1000., PostStepPosY[n]/1000., PostStepPosZ[n]/1000.);
	  TVector3 globalHitS4 = MCToGlobalCoords(mcHitS4);
	  // std::cout<<"x, y, z "<<globalHitS4.X()<<", "<<globalHitS4.Y()<<", "<<globalHitS4.Z()<<std::endl;
	  hS4->Fill(getThetaFromGlobal(globalHitS4), getPhiFromGlobal(globalHitS4), weight);
	  hS4XY->Fill(globalHitS4.X(), globalHitS4.Y(), weight);
	  hS4XZ->Fill(globalHitS4.X(), globalHitS4.Z(), weight);
	  hS4YZ->Fill(globalHitS4.Z(), globalHitS4.Y(), weight);
	  if (isS2) hS4withS2->Fill(getThetaFromGlobal(globalHitS4), getPhiFromGlobal(globalHitS4), weight);
	  else hS4noS2->Fill(getThetaFromGlobal(globalHitS4), getPhiFromGlobal(globalHitS4), weight);
	} // Is in S4
	else if (((PreStepVolumeID[n]!=8 && PostStepVolumeID[n]==8) ||
		  (PreStepVolumeID[n]!=10 && PostStepVolumeID[n]==10)) &&
		 PDG[n]==2212 && PostStepMom[n] > 30.) {
	  TVector3 mcHit(PostStepPosX[n]/1000., PostStepPosY[n]/1000., PostStepPosZ[n]/1000.);
	  TVector3 globalHit = MCToGlobalCoords(mcHit);
	  hDoors->Fill(getThetaFromGlobal(globalHit), getPhiFromGlobal(globalHit), weight);

	  if (isS2) hDoorsWithS2->Fill(getThetaFromGlobal(globalHit), getPhiFromGlobal(globalHit), weight);
	  else hDoorsNoS2->Fill(getThetaFromGlobal(globalHit), getPhiFromGlobal(globalHit), weight);

	  // if (PreStepVolumeID[n]!=8 && PostStepVolumeID[n]==8) {
	    hDoor1XY->Fill(globalHit.X(), globalHit.Y(), weight);
	    hDoor1XZ->Fill(globalHit.X(), globalHit.Z(), weight);
	    hDoor1YZ->Fill(globalHit.Y(), globalHit.Z(), weight);
	  // }
	  // else {
	  //   hDoor2XY->Fill(globalHit.X(), globalHit.Y(), weight);
	  //   hDoor2XZ->Fill(globalHit.X(), globalHit.Z(), weight);
	  //   hDoor2YZ->Fill(globalHit.Y(), globalHit.Z(), weight); 
	  // }
	} // Is in vessel doors
	else if (PreStepVolumeID[n]!=6 && PostStepVolumeID[n]==6 &&
		 PDG[n]==2212 && PostStepMom[n] > 30.) {
	  TVector3 mcHit(PostStepPosX[n]/1000., PostStepPosY[n]/1000., PostStepPosZ[n]/1000.);
	  TVector3 globalHit = MCToGlobalCoords(mcHit);
	  hTpc->Fill(getThetaFromGlobal(globalHit), getPhiFromGlobal(globalHit), weight);
	} // Is in TPC active area

	if (((PreStepVolumeID[n]!=8 && PostStepVolumeID[n]==8) ||
	     (PreStepVolumeID[n]!=10 && PostStepVolumeID[n]==10) ||
	     (PreStepVolumeID[n]!=3 && PostStepVolumeID[n]==3)) &&
	    PDG[n]==2212 && PostStepMom[n] > 30.) {
	  TVector3 mcHit(PostStepPosX[n]/1000., PostStepPosY[n]/1000., PostStepPosZ[n]/1000.);
	  TVector3 globalHit = MCToGlobalCoords(mcHit);
	  hSteelXY->Fill(globalHit.X(), globalHit.Y(), weight);
	  hSteelXZ->Fill(globalHit.X(), globalHit.Z(), weight);
	  hSteelYZ->Fill(globalHit.Y(), globalHit.Z(), weight);
	}
      } // Loop over number of steps
    } // Loop over entries

    stepTree->GetEntry(stepTree->GetEntries()-1);
    int lastSpill = spill;
    hS3->Scale(1./lastSpill);
    hS3withS2->Scale(1./lastSpill);
    hS3noS2->Scale(1./lastSpill);

    hS4->Scale(1./lastSpill);
    hS4withS2->Scale(1./lastSpill);
    hS4noS2->Scale(1./lastSpill);

    hDoors->Scale(1./lastSpill);
    hDoorsWithS2->Scale(1./lastSpill);
    hDoorsNoS2->Scale(1./lastSpill);

    hDoor1XY->Scale(1./lastSpill);
    hDoor1YZ->Scale(1./lastSpill);
    hDoor1XZ->Scale(1./lastSpill);

    hDoor2XY->Scale(1./lastSpill);
    hDoor2YZ->Scale(1./lastSpill);
    hDoor2XZ->Scale(1./lastSpill);

    hSteelXY->Scale(1./lastSpill);
    hSteelYZ->Scale(1./lastSpill);
    hSteelXZ->Scale(1./lastSpill);

    hS4XY->Scale(1./lastSpill);
    hS4YZ->Scale(1./lastSpill);
    hS4XZ->Scale(1./lastSpill);

    hTpc->Scale(1./lastSpill);

    fout->cd();

    hS3->Write();
    hS3withS2->Write();
    hS3noS2->Write();

    hS4->Write();
    hS4withS2->Write();
    hS4noS2->Write();

    hDoors->Write();
    hDoorsWithS2->Write();
    hDoorsNoS2->Write();

    hDoor1XY->Write();
    hDoor1XZ->Write();
    hDoor1YZ->Write();

    hDoor2XY->Write();
    hDoor2XZ->Write();
    hDoor2YZ->Write();

    hSteelXY->Write();
    hSteelXZ->Write();
    hSteelYZ->Write();

    hTpc->Write();

    hS4XY->Write();
    hS4YZ->Write();
    hS4XZ->Write();

    fin->Close();
    delete fin;
  } //  Loop over blocks

  // Lines for drawing detector components
  TGraph *lS2 = new TGraph();
  lS2->SetLineColor(kOrange+1);
  lS2->SetLineWidth(3);
  lS2->SetPoint(lS2->GetN(), getThetaFromGlobal(D2_UTL), getPhiFromGlobal(D2_UTL));
  lS2->SetPoint(lS2->GetN(), getThetaFromGlobal(D2_UTR), getPhiFromGlobal(D2_UTR));
  lS2->SetPoint(lS2->GetN(), getThetaFromGlobal(D2_UBR), getPhiFromGlobal(D2_UBR));
  lS2->SetPoint(lS2->GetN(), getThetaFromGlobal(D2_UBL), getPhiFromGlobal(D2_UBL));
  lS2->SetPoint(lS2->GetN(), getThetaFromGlobal(D2_UTL), getPhiFromGlobal(D2_UTL));

  TGraph *lS3 = new TGraph();
  lS3->SetLineColor(kCyan+2);
  lS3->SetLineWidth(3);
  lS3->SetPoint(lS3->GetN(), getThetaFromGlobal(S3_TL), getPhiFromGlobal(S3_TL));
  lS3->SetPoint(lS3->GetN(), getThetaFromGlobal(S3_TR), getPhiFromGlobal(S3_TR));
  lS3->SetPoint(lS3->GetN(), getThetaFromGlobal(S3_BR), getPhiFromGlobal(S3_BR));
  lS3->SetPoint(lS3->GetN(), getThetaFromGlobal(S3_BL), getPhiFromGlobal(S3_BL));
  lS3->SetPoint(lS3->GetN(), getThetaFromGlobal(S3_TL), getPhiFromGlobal(S3_TL));

  TGraph *lS4 = new TGraph();
  lS4->SetLineColor(kRed+2);
  lS4->SetLineWidth(3);
  lS4->SetPoint(lS4->GetN(), getThetaFromGlobal(S4_TL), getPhiFromGlobal(S4_TL));
  lS4->SetPoint(lS4->GetN(), getThetaFromGlobal(S4_TR), getPhiFromGlobal(S4_TR));
  lS4->SetPoint(lS4->GetN(), getThetaFromGlobal(S4_BR), getPhiFromGlobal(S4_BR));
  lS4->SetPoint(lS4->GetN(), getThetaFromGlobal(S4_BL), getPhiFromGlobal(S4_BL));
  lS4->SetPoint(lS4->GetN(), getThetaFromGlobal(S4_TL), getPhiFromGlobal(S4_TL));

  TGraph *lS4XY = new TGraph();
  lS4XY->SetLineColor(kRed+2);
  lS4XY->SetLineWidth(3);
  lS4XY->SetPoint(lS4XY->GetN(), S4_TL.X(), S4_TL.Y());
  lS4XY->SetPoint(lS4XY->GetN(), S4_TR.X(), S4_TR.Y());
  lS4XY->SetPoint(lS4XY->GetN(), S4_BR.X(), S4_BR.Y());
  lS4XY->SetPoint(lS4XY->GetN(), S4_BL.X(), S4_BL.Y());
  lS4XY->SetPoint(lS4XY->GetN(), S4_TL.X(), S4_TL.Y());
  lS4XY->Write("lS4XY");

  TGraph *lS4XZ = new TGraph();
  lS4XZ->SetLineColor(kRed+2);
  lS4XZ->SetLineWidth(3);
  lS4XZ->SetPoint(lS4XZ->GetN(), S4_TL.X(), S4_TL.Z());
  lS4XZ->SetPoint(lS4XZ->GetN(), S4_TR.X(), S4_TR.Z());
  lS4XZ->SetPoint(lS4XZ->GetN(), S4_BR.X(), S4_BR.Z());
  lS4XZ->SetPoint(lS4XZ->GetN(), S4_BL.X(), S4_BL.Z());
  lS4XZ->SetPoint(lS4XZ->GetN(), S4_TL.X(), S4_TL.Z());
  lS4XZ->Write("lS4XZ");

  TGraph *lS4YZ = new TGraph();
  lS4YZ->SetLineColor(kRed+2);
  lS4YZ->SetLineWidth(3);
  lS4YZ->SetPoint(lS4YZ->GetN(), S4_TL.Z(), S4_TL.Y());
  lS4YZ->SetPoint(lS4YZ->GetN(), S4_TR.Z(), S4_TR.Y());
  lS4YZ->SetPoint(lS4YZ->GetN(), S4_BR.Z(), S4_BR.Y());
  lS4YZ->SetPoint(lS4YZ->GetN(), S4_BL.Z(), S4_BL.Y());
  lS4YZ->SetPoint(lS4YZ->GetN(), S4_TL.Z(), S4_TL.Y());
  lS4YZ->Write("lS4YZ");

  TMarker *m1 = new TMarker(getThetaFromGlobal(D4_LC), getPhiFromGlobal(D4_LC), 8);
  m1->SetMarkerColor(kGreen+2);
  m1->SetMarkerSize(2);
  m1->Write("m1");
  TMarker *m2 = new TMarker(getThetaFromGlobal(D4_RC), getPhiFromGlobal(D4_RC), 8);
  m2->SetMarkerColor(kGreen+2);
  m2->SetMarkerSize(2);
  m2->Write("m2");

  TGraph *lTpcUs = new TGraph();
  lTpcUs->SetLineColor(kGreen+2);
  lTpcUs->SetLineStyle(7);
  lTpcUs->SetLineWidth(3);
  lTpcUs->SetPoint(lTpcUs->GetN(), getThetaFromGlobal(TPCActiveBUL), getPhiFromGlobal(TPCActiveBUL));
  lTpcUs->SetPoint(lTpcUs->GetN(), getThetaFromGlobal(TPCActiveBUR), getPhiFromGlobal(TPCActiveBUR));
  lTpcUs->SetPoint(lTpcUs->GetN(), getThetaFromGlobal(TPCActiveTUR), getPhiFromGlobal(TPCActiveTUR));
  lTpcUs->SetPoint(lTpcUs->GetN(), getThetaFromGlobal(TPCActiveTUL), getPhiFromGlobal(TPCActiveTUL));
  lTpcUs->SetPoint(lTpcUs->GetN(), getThetaFromGlobal(TPCActiveBUL), getPhiFromGlobal(TPCActiveBUL));

  TGraph *lTpcDs = new TGraph();
  lTpcDs->SetLineColor(kGreen+2);
  lTpcDs->SetLineWidth(3);
  lTpcDs->SetPoint(lTpcDs->GetN(), getThetaFromGlobal(TPCActiveBDL), getPhiFromGlobal(TPCActiveBDL));
  lTpcDs->SetPoint(lTpcDs->GetN(), getThetaFromGlobal(TPCActiveBDR), getPhiFromGlobal(TPCActiveBDR));
  lTpcDs->SetPoint(lTpcDs->GetN(), getThetaFromGlobal(TPCActiveTDR), getPhiFromGlobal(TPCActiveTDR));
  lTpcDs->SetPoint(lTpcDs->GetN(), getThetaFromGlobal(TPCActiveTDL), getPhiFromGlobal(TPCActiveTDL));
  lTpcDs->SetPoint(lTpcDs->GetN(), getThetaFromGlobal(TPCActiveBDL), getPhiFromGlobal(TPCActiveBDL));

  TMarker *m1GlobalXY = new TMarker(D4_LC.X(), D4_LC.Y(), 8);
  m1GlobalXY->SetMarkerSize(2);
  m1GlobalXY->SetMarkerColor(kGreen+2);
  TMarker *m1GlobalXZ = new TMarker(D4_LC.X(), D4_LC.Z(), 8);
  m1GlobalXZ->SetMarkerSize(2);
  m1GlobalXZ->SetMarkerColor(kGreen+2);
  TMarker *m1GlobalYZ = new TMarker(D4_LC.Y(), D4_LC.Z(), 8);
  m1GlobalYZ->SetMarkerSize(2);
  m1GlobalYZ->SetMarkerColor(kGreen+2);

  TMarker *m2GlobalXY = new TMarker(D4_RC.X(), D4_RC.Y(), 8);
  m2GlobalXY->SetMarkerSize(2);
  m2GlobalXY->SetMarkerColor(kGreen+2);
  TMarker *m2GlobalXZ = new TMarker(D4_RC.X(), D4_RC.Z(), 8);
  m2GlobalXZ->SetMarkerSize(2);
  m2GlobalXZ->SetMarkerColor(kGreen+2);
  TMarker *m2GlobalYZ = new TMarker(D4_RC.Y(), D4_RC.Z(), 8);
  m2GlobalYZ->SetMarkerSize(2);
  m2GlobalYZ->SetMarkerColor(kGreen+2);

  m1GlobalXY->Write("m1GlobalXY");
  m1GlobalXZ->Write("m1GlobalXZ");
  m1GlobalYZ->Write("m1GlobalYZ");
  m2GlobalXY->Write("m2GlobalXY");
  m2GlobalXZ->Write("m2GlobalXZ");
  m2GlobalYZ->Write("m2GlobalYZ");

  lS2->Write("lS2");
  lS3->Write("lS3");
  lS4->Write("lS4");
  lTpcUs->Write("lTpcUs");
  lTpcDs->Write("lTpcDs");

  fout->Close();
  delete fout;
} // mcGeomComp
