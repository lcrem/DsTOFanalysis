// s3CutOpt.C
// Determine optimal amplitude cut for each of the bars

void s3CutOpt(const char* saveDir,
	      const int block,
	      const char* ustofDir)
{
  gROOT->SetBatch(kTRUE);


  const char* nustof;
  if (block==0) nustof = Form("%sData_2018_8_31_b2_800MeV_0block.root", ustofDir);
  else if (block==1) nustof = Form("%sData_2018_9_1_b4_800MeV_1block_bend4cm.root", ustofDir);
  else if (block==2) nustof = Form("%sData_2018_9_1_b2_800MeV_2block_bend4cm.root", ustofDir);
  else if (block==3) nustof = Form("%sData_2018_9_1_b3_800MeV_3block_bend4cm.root", ustofDir);
  else if (block==4) nustof = Form("%sData_2018_9_1_b8_800MeV_4block_bend4cm.root", ustofDir);

  // Open file
  TFile *futof = new TFile(nustof, "read");

  double tToF[50];
  float xToF[50];
  float yToF[50];
  float A1ToF[50];
  float A2ToF[50];
  double tTrig;
  double tS1;
  double tSoSd;
  int nhit;
  int nBar[50];

  TTree *tree = (TTree*)futof->Get("tree");

  tree->SetBranchAddress("xToF", xToF);
  tree->SetBranchAddress("yToF", yToF);
  tree->SetBranchAddress("A1ToF", A1ToF);
  tree->SetBranchAddress("A2ToF", A2ToF);
  tree->SetBranchAddress("nhit", &nhit);
  tree->SetBranchAddress("tS1", &tS1);
  tree->SetBranchAddress("tToF", tToF);
  tree->SetBranchAddress("tTrig", &tTrig);
  tree->SetBranchAddress("tSoSd", &tSoSd);
  tree->SetBranchAddress("nBar", nBar);

  // Create vector of 2D hists
  std::vector<TH2D*> histVec;
  for (int i=0; i < 22; i++) {
    TH2D *hAmpTime = new TH2D(Form("hAmpTime%d", i), Form("Bar %d; #delta t / ns; Amplitude / V", i), 100, -40, 60, 100, 0, 0.5);
    histVec.push_back(hAmpTime);
  } // for (int i=0; i < 22; i++) 

  for (int t=0; t<tree->GetEntries(); t++) {
    tree->GetEntry(t);
    if (nhit == 0) {
      histVec[nBar]->Fill((tToF[0] - tS1), A1ToF[0]);
    }
  }
  TFile* fout = new TFile(Form("%s/s3CutOpt.root",saveDir), "recreate");

}
