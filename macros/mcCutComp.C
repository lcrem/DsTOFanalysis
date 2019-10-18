// mcCutComp.C

const int nVertSteps = 10;
const int nMomSteps  = 10;
const int momMin = 0.15;
const int momMin = 0.3;
vector<double> mcS3   = {56590, 67273, 60652, 35781, 104318};
vector<double> dataS3 = {1983, 1656, 1325, 899, 136.3};
vector<double> momCuts = {0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325};

void mcCutComp(const char* outFile, const char* dataFile="/scratch0/sjones/plots/angularDistS4_newSample/withTree/angularDistS4Plots.root", const char* mcDir="/scratch0/tnonnenm/FixedMC/Prototype-HPTPC-MC/wdir/moderator_data/")
{

  TFile *fout = new TFile(outFile, "read");
  TFile *dataFile = new TFile(dataFile, "read");
  // Loop over blocks
  for (int b=0; b<4; b++) {
    TGraph *gData = new TGraph();
    TGraph *gMC   = new TGraph();
    gData->SetName(Form("gData%d", b));
    gMC->SetName(Form("gMC%d", b));
    vector<double> nProData, nProMC;
    nProData.resize(momCuts.size(), 0.);
    nProMC.resize(momCuts.size(), 0.);
    cout<<"Doing "<<b<<" block sample"<<endl;
    TFile *mcFile   = new TFile(Form("%s/%dBlock.root", mcDir, b), "read");
    TTree *dataTree = (TTree*)dataFile->Get(Form("protonTree%d",b));
    double mom;
    double mcX, mcY, mcZ;
    double weight;
    dataTree->SetBranchAddress("mom", &mom);
    dataTree->SetBranchAddress("weight", &weight);
    dataTree->SetBranchAddress("mcX", &mcX);
    dataTree->SetBranchAddress("mcY", &mcY);
    dataTree->SetBranchAddress("mcZ", &mcZ);
    dataTree->SetBranchAddress("spill", &spill);
    // Loop over proton tree
    dataTree->GetEntry(dataTree->GetEntries()-1);
    int nSpills = spill;
    for (int i=0; i<dataTree->GetEntries(); i++) {
      dataTree->GetEntry(i);
      for (int c=0; c<momCuts.size(); c++) {
	if (mom > momCuts.at(c)) nProData.at(c) += weight;
      }      
    } // Loop over proton tree
    // Per spill numbers
    for (int c=0; c<momCuts.size(); c++) {
      nProData.at(c) /= nSpills;
      gData->SetPoint(gData->GetN(), momCuts.at(c), nProData.at(c)/dataS3.at(b)); 
    }

    fout->cd();
    gData->Write();
    gMC->Write();

    mcFile->Close();
    delete mcFile;
  } // Loop over blocks

  fout->Close();
  dataFile->Close();
  delete fout;
  delete dataFile;
} // mcCutComp
