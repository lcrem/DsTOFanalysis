// absFluxS3.C
// Similar to absFluxS4.C but for S3
// Look at difference between only using S1 and using S1+S2

void absFluxS3(const char* saveDir,
	       const char* ustofDir = "/zfs_home/sjones/mylinktoutof/")
{
  gROOT->SetBatch(kTRUE);
 TFile *fout = new TFile(Form("%s/absFluxS3Plots.root",saveDir), "recreate");

  // Edges of S3 in beam coordinate system
  const double s3StartX = -0.5168;
  const double s3EndX   = 0.9970;
  const double s3s1StartY = 9.0569 + 1.77;
  const double s3s1EndY   = 8.9146 + 1.77;

  // Define the runs to be used for varying number of blocks
  const char* str0Block = "Data_2018_8_31_b2_800MeV_0block.root";
  const char* str1Block = "Data_2018_9_1_b4_800MeV_1block_bend4cm.root";
  const char* str2Block = "Data_2018_9_1_b2_800MeV_2block_bend4cm.root";
  const char* str3Block = "Data_2018_9_1_b3_800MeV_3block_bend4cm.root";
  const char* str4Block = "Data_2018_9_1_b8_800MeV_4block_bend4cm.root";

  THStack *hsXAngleS1 = new THStack("hsXAngleS1", "Angular distribution of S3 hits with S1 trigger only; #theta / degrees; Events / spill");
  THStack *hsXAngleS1S2 = new THStack("hsXAngleS1S2", "Angular distribution of S3 hits with S1 & S2; #theta / degrees; Events / spill");
  THStack *hsXS1 = new THStack("hsXS1", "x distribution of S3 hits with S1 trigger only; x / m; Events / spill");
  THStack *hsXS1S2 = new THStack("hsXS1S2", "x distribution of S3 hits with S1 & S2; x / m; Events / spill");
  TLegend *leg = new TLegend(0.65, 0.6, 0.88, 0.8);

  for (int nBlocks = 0; nBlocks <=4; nBlocks++) {
    int nSpills = 0;
    char* nustof;
    if (nBlocks == 0) nustof = Form("%s/%s", ustofDir, str0Block);
    else if (nBlocks == 1) nustof = Form("%s/%s", ustofDir, str1Block);
    else if (nBlocks == 2) nustof = Form("%s/%s", ustofDir, str2Block);
    else if (nBlocks == 3) nustof = Form("%s/%s", ustofDir, str3Block);
    else if (nBlocks == 4) nustof = Form("%s/%s", ustofDir, str4Block);
    
    // Histograms to be filled 
    TH1D *hXS1 = new TH1D(Form("hXS1%d", nBlocks), Form("x coordinate distribution of hits in S3 (S1 trigger only), %d blocks; x / m; Events / spill", nBlocks), 100, -.7, 1.15);
    TH1D *hXS1S2 = new TH1D(Form("hXS1S2%d", nBlocks), Form("x coordinate distribution of hits in S3 (S1 & S2 triggers), %d blocks; x / m; Events / spill", nBlocks), 100, -.7, 1.15);
    TH1D *hXAngleS1 = new TH1D(Form("hXAngleS1%d", nBlocks), Form("Angular distribution of hits in S3 (S1 trigger only), %d blocks; #theta / degrees; Events / spill", nBlocks), 100, -3.8, 6.2);
    TH1D *hXAngleS1S2 = new TH1D(Form("hXAngleS1S2%d", nBlocks), Form("Angular distribution of hits in S3 (S1 & S2 triggers), %d blocks; #theta / degrees; Events / spill", nBlocks), 100, -3.8, 6.2);
    // Open appropriate ustof file
    TFile *fin = new TFile(nustof, "read");
    TTree *tree = (TTree*)fin->Get("tree");
    double tTrig;
    double tS1;
    double tSoSd;
    float xToF[50];
    float yToF[50];
    int nhit;
    tree->SetBranchAddress("tTrig", &tTrig);
    tree->SetBranchAddress("tS1", &tS1);
    tree->SetBranchAddress("xToF", xToF);
    tree->SetBranchAddress("yToF", yToF);
    tree->SetBranchAddress("nhit", &nhit);
    tree->SetBranchAddress("tSoSd", &tSoSd);

    double lastSpill = 0.;

    for (int t=0; t < tree->GetEntries(); t++) {
      tree->GetEntry(t);
      // Assume the width of S3 is 180cm -- need to check this
      // Distance between the surveyed end points is 152cm
      // Assume S3 juts out 14cm each side of the surveyed points
      // Select events where hit multiplicity is 1? -- suggestion from A. Korzenev
      for (int n=0; n < nhit; n++) {
	// There is an S1 and S2 hit associated with this event
	if (tTrig != 0) {
	  double positionX = ((xToF[n] - 4.) / 152.)*(s3EndX - s3StartX) + s3StartX;
	  double positionY =((yToF[n] - 4.) / 152.)*(s3s1EndY - s3s1StartY) + s3s1StartY;
	  hXS1S2->Fill(positionX);
	  double angleOffAxis = TMath::ATan(positionX / positionY) * 180. / TMath::Pi();
	  hXAngleS1S2->Fill(angleOffAxis);
	} // if (tTrig != 0)
	// There is an S1 hit only associated with this event
	else {
	  double positionX = ((xToF[n] - 4.) / 152.)*(s3EndX - s3StartX) + s3StartX;	 
	  double positionY = ((xToF[n] - 4.) / 152.)*(s3s1EndY - s3s1StartY) + s3s1StartY;
	  hXS1->Fill(positionX);
	  double angleOffAxis = TMath::ATan(positionX / positionY) * 180. / TMath::Pi();
	  hXAngleS1->Fill(angleOffAxis);
	}
      } // for (int n=0; n < nhit; n++)

      // Count number of spills
      if (tSoSd != lastSpill && tSoSd -lastSpill > 1e9) {
	nSpills++;
	lastSpill = tSoSd;
      } // if (tSoSd != lastSpill && tSoSd -lastSpill > 1e9) 

    } // for (int t=0; t < tree->GetEntries(); t++)

    fin->Close();
    delete fin;
    fout->cd();

    TCanvas *cXS1S2 = new TCanvas(Form("cXS1S2_%d", nBlocks));
    hXS1S2->Scale(1. / nSpills);
    hXS1S2->Draw("hist");
    hXS1S2->Write();
    cXS1S2->Print(Form("%s/%d_absFluxXS123.png", saveDir, nBlocks));
    cXS1S2->Print(Form("%s/%d_absFluxXS123.pdf", saveDir, nBlocks));
    TCanvas *cXAngleS1S2 = new TCanvas(Form("cXAngleS1S2_%d", nBlocks));
    hXAngleS1S2->Scale(1. / nSpills);
    hXAngleS1S2->Draw("hist");
    hXAngleS1S2->Write();
    cXAngleS1S2->Print(Form("%s/%d_absFluxXAngleS123.png", saveDir, nBlocks));
    cXAngleS1S2->Print(Form("%s/%d_absFluxXAngleS123.pdf", saveDir, nBlocks));

    TCanvas *cXS1 = new TCanvas(Form("cXS1_%d", nBlocks));
    hXS1->Scale(1. / nSpills);
    hXS1->Draw("hist");
    hXS1->Write();
    cXS1->Print(Form("%s/%d_absFluxXS13.png", saveDir, nBlocks));
    cXS1->Print(Form("%s/%d_absFluxXS13.pdf", saveDir, nBlocks));
    TCanvas *cXAngleS1 = new TCanvas(Form("cXAngleS1_%d", nBlocks));
    hXAngleS1->Scale(1. / nSpills);
    hXAngleS1->Draw("hist");
    hXAngleS1->Write();
    cXAngleS1->Print(Form("%s/%d_absFluxXAngleS13.png", saveDir, nBlocks));
    cXAngleS1->Print(Form("%s/%d_absFluxXAngleS13.pdf", saveDir, nBlocks));

    // And do them all together in different colours
    if (nBlocks == 0) {
      hXAngleS1S2->SetLineColor(kBlue);
      hXAngleS1->SetLineColor(kBlue);
      hXS1S2->SetLineColor(kBlue);
      hXS1->SetLineColor(kBlue);
      leg->AddEntry(hXAngleS1S2, "0 blocks", "l");
    }
    else if (nBlocks == 1) {
      hXAngleS1S2->SetLineColor(kRed);
      hXAngleS1->SetLineColor(kRed);
      hXS1S2->SetLineColor(kRed);
      hXS1->SetLineColor(kRed);
      leg->AddEntry(hXAngleS1S2, "1 block", "l");
    }
    else if (nBlocks == 2) {
      hXAngleS1S2->SetLineColor(kBlack);
      hXAngleS1->SetLineColor(kBlack);
      hXS1S2->SetLineColor(kBlack);
      hXS1->SetLineColor(kBlack);
      leg->AddEntry(hXAngleS1S2, "2 blocks", "l");
    }
    else if (nBlocks == 3){
      hXAngleS1S2->SetLineColor(kGreen+2);
      hXAngleS1->SetLineColor(kGreen+2);
      hXS1S2->SetLineColor(kGreen+2);
      hXS1->SetLineColor(kGreen+2);
      leg->AddEntry(hXAngleS1S2, "3 blocks", "l");
    }
    else {
      hXAngleS1S2->SetLineColor(kMagenta);
      hXAngleS1->SetLineColor(kMagenta);
      hXS1S2->SetLineColor(kMagenta);
      hXS1->SetLineColor(kMagenta);
      leg->AddEntry(hXAngleS1S2, "4 blocks", "l");
    }
    hsXAngleS1S2->Add(hXAngleS1S2);
    hsXAngleS1->Add(hXAngleS1);
    hsXS1S2->Add(hXS1S2);
    hsXS1->Add(hXS1);

  } // for (int nBlocks = 0; nBlocks <=4; nBlocks++)
  fout->cd();

  TCanvas *cXAngleS1S2 = new TCanvas("cXAngleS1S2");
  hsXAngleS1S2->Draw("hist nostack");
  leg->Draw();
  hsXAngleS1S2->Write();
  leg->Write();
  cXAngleS1S2->Print(Form("%s/absFluxXAngleS123.png", saveDir));
  cXAngleS1S2->Print(Form("%s/absFluxXAngleS123.pdf", saveDir));
  TCanvas *cXAngleS1 = new TCanvas("cXAngleS1");
  hsXAngleS1->Draw("hist nostack");
  leg->Draw();
  hsXAngleS1->Write();
  cXAngleS1->Print(Form("%s/absFluxXAngleS13.png", saveDir));
  cXAngleS1->Print(Form("%s/absFluxXAngleS13.pdf", saveDir));
  TCanvas *cXS1S2 = new TCanvas("cXS1S2");
  hsXS1S2->Draw("hist nostack");
  leg->Draw();
  hsXS1S2->Write();
  cXS1S2->Print(Form("%s/absFluxXS123.png", saveDir));
  cXS1S2->Print(Form("%s/absFluxXS123.pdf", saveDir));
  TCanvas *cXS1 = new TCanvas("cXS1");
  hsXS1->Draw("hist nostack");
  leg->Draw();
  hsXS1->Write();
  cXS1->Print(Form("%s/absFluxXS13.png", saveDir));
  cXS1->Print(Form("%s/absFluxXS13.pdf", saveDir));

  fout->Close();
}
