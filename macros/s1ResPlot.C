// s1ResPlot.C
// Tests/remakes a plot of the S1 timing resolution
#include "UsefulFunctions.C"

void s1ResPlot(const char* outFile, 
	       const char* utofDir="/scratch0/dbrailsf/data_backup/utof_backup_firsthitpinnedtounixtime/Data_root_v3_wo_walk_corr/")
{
  TFile *fout = new TFile(outFile, "recreate");

  TChain *tree = new TChain("tree");
  // tree->Add(Form("%s/Data_2018_9_12_b1.root", utofDir));
  // tree->Add(Form("%s/Data_2018_9_13_b1.root", utofDir));
  // tree->Add(Form("%s/Data_2018_9_13_bc2.root", utofDir));
  tree->Add(Form("%s/Data_2018_9_14_b1.root", utofDir));
  tree->Add(Form("%s/Data_2018_9_15_b1.root", utofDir));
  tree->Add(Form("%s/Data_2018_9_15_b2.root", utofDir));
  tree->Add(Form("%s/Data_2018_9_16_b1.root", utofDir));
  tree->Add(Form("%s/Data_2018_9_16_b2.root", utofDir));

  double tS1_0, tS1_1, tS1_2, tS1_3;
  tree->SetBranchAddress("tS1_0", &tS1_0);
  tree->SetBranchAddress("tS1_1", &tS1_1);
  tree->SetBranchAddress("tS1_2", &tS1_2);
  tree->SetBranchAddress("tS1_3", &tS1_3);

  TH1D *hDiff = new TH1D("hDiff", "Measurement of the difference in S1 PMT trigger times; #Deltat for S1 PMTs / ns; Events", 45, -0.2, 0.2);
  setHistAttr(hDiff);
  hDiff->SetLineColor(kBlack);

  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
    hDiff->Fill(((tS1_1+tS1_2)-(tS1_0+tS1_3))/4.);
  }

  fout->cd();
  hDiff->Write();

  // fin->Close();
  // delete fin;
  fout->Close();
  delete fout;
}
