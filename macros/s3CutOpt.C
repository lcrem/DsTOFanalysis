// s3CutOpt.C
// Determine optimal amplitude cut for each of the bars

void s3CutOpt(const char* saveDir,
	      const int block,
	      const char* ustofDir="/nfs/scratch0/dbrailsf/data_backup/utof_backup_firsthitpinnedtounixtime/Data_root_v3_wo_walk_corr/")
{
  gROOT->SetBatch(kTRUE);


  const char* nustof;
  if (block==0) nustof = Form("%sData_2018_8_31_b2_800MeV_0block.root", ustofDir);
  else if (block==1) nustof = Form("%sData_2018_9_1_b4_800MeV_1block_bend4cm.root", ustofDir);
  else if (block==2) nustof = Form("%sData_2018_9_1_b2_800MeV_2block_bend4cm.root", ustofDir);
  else if (block==3) nustof = Form("%sData_2018_9_1_b3_800MeV_3block_bend4cm.root", ustofDir);
  else if (block==4) nustof = Form("%sData_2018_8_29_b4.root", ustofDir);

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
  std::vector<TH2D*> histA1Vec;
  std::vector<TH2D*> histA2Vec;
  for (int i=0; i < 22; i++) {
    TH2D *hA1Time = new TH2D(Form("hA1Time%d", i), Form("Bar %d; #delta t / ns; A1 / V", i), 100, 20, 120, 100, 0, 0.45);
    TH2D *hA2Time = new TH2D(Form("hA2Time%d", i), Form("Bar %d; #delta t / ns; A2 / V", i), 100, 20, 120, 100, 0, 0.45);
    histA1Vec.push_back(hA1Time);
    histA2Vec.push_back(hA2Time);
  } // for (int i=0; i < 22; i++) 
  // Fill appropriate bar histogram
  for (int t=0; t<tree->GetEntries(); t++) {
    tree->GetEntry(t);
    if (nhit == 1) {
      histA1Vec.at(nBar[0])->Fill((tToF[0] - tS1), A1ToF[0]);
      histA2Vec.at(nBar[0])->Fill((tToF[0] - tS1), A2ToF[0]);
    }
  } // Loop over entries in tree

  TFile* fout = new TFile(Form("%s/s3CutOpt%dblock.root",saveDir,block), "recreate");
  cout<<"Writing to "<<Form("%s/s3CutOpt%dblock.root",saveDir,block)<<endl;
  // Write the histograms to file
  for (int i=0; i<22; i++) {
    histA1Vec[i]->Write();
    histA2Vec[i]->Write();
  }
  fout->Close();
}
