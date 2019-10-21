// mcCutComp.C

const int nVertSteps = 10;
const int nMomSteps  = 10;
const double momMin = 0.1;
const double momMax = 0.325;
const vector<double> mcS3    = {56590, 67273, 60652, 35781, 104318};
const vector<double> dataS3  = {1983, 1656, 1325, 899, 136.3};
const vector<double> momCuts = {100, 125, 150, 175, 200, 225, 250, 275, 300, 325};
const vector<double> sliceEdges = {-1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, 
				   -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4};

void mcCutComp(const char* outFile, const char* dataFileDir="/scratch0/sjones/plots/angularDistS4_newSample/withTree/angularDistS4Plots.root", const char* mcDir="/scratch0/tnonnenm/FixedMC/Prototype-HPTPC-MC/wdir/moderator_data/")
{

  TFile *fout = new TFile(outFile, "recreate");
  TFile *dataFile = new TFile(dataFileDir, "read");
  // Loop over blocks
  for (int b=0; b<5; b++) {
    cout<<"Doing "<<b<<" block sample"<<endl;
    TMultiGraph *mgDataMCComp = new TMultiGraph(Form("mgDataMCComp%d",b), Form("Data - MC comparison for various minimum momenta, %d blocks; Momentum cut [GeV/c]; S4 / S3; ", b));
    THStack *hsDataMCCompVert = new THStack(Form("hsDataMCCompVert%d", b), Form("Data - MC comparison across S4, %d blocks; x position / m; Ratio", b));
    TGraph *gData = new TGraph();
    TGraph *gMC   = new TGraph();
    gData->SetName(Form("gData%d", b));
    gMC->SetName(Form("gMC%d", b));
    vector<double> nProData, nProMC;
    nProData.resize(momCuts.size(), 0.);
    nProMC.resize(momCuts.size(), 0.);
    TH1D *hDataVert = new TH1D(Form("hDataVert%d",b), Form("S4/S3 ratio in S4 data, %d blocks; x position / m; S4 / S3", b), 14, -1., 0.4);
    TH1D *hMCVert = new TH1D(Form("hMCVert%d",b), Form("S4/S3 ratio in MC, %d blocks; x position / m; S4 / S3", b), 14, -1., 0.4);
    hMCVert->SetLineWidth(2);
    hDataVert->SetLineWidth(2);
    hMCVert->SetLineColor(kRed);
    hDataVert->SetLineColor(kBlack);
    // Get proton data
    TTree *dataTree = (TTree*)dataFile->Get(Form("protonTree%d",b));
    double mom;
    double mcX, mcY, mcZ;
    double weight;
    int spill;
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
      hDataVert->Fill(mcX, weight);
      for (int c=0; c<momCuts.size(); c++) {
	if (mom > momCuts.at(c)) nProData.at(c) += weight;
      }      
    } // Loop over proton tree
    // Per spill numbers and divide by S3 number
    for (int c=0; c<momCuts.size(); c++) {
      cout<<nProData.at(c)<<endl;
      nProData.at(c) /= (double)nSpills;
      gData->SetPoint(gData->GetN(), momCuts.at(c), nProData.at(c)/dataS3.at(b)); 
    }

    // Now do MC
    TFile *mcFile   = new TFile(Form("%s/%dBlock.root", mcDir, b), "read");
    TTree *stepTree = (TTree*)mcFile->Get("stepTree");
    UInt_t NSteps;
    float PostStepMom[40000];
    float PostStepPosX[40000];
    int PreStepVolumeID[40000];
    int PostStepVolumeID[40000];
    stepTree->SetBranchAddress("NSteps", &NSteps);
    stepTree->SetBranchAddress("PostStepMom", PostStepMom);
    stepTree->SetBranchAddress("PostStepPosX", PostStepPosX);
    stepTree->SetBranchAddress("PreStepVolumeID", PreStepVolumeID);
    stepTree->SetBranchAddress("PostStepVolumeID", PostStepVolumeID);
    // Loop over entries
    for (int i=0; i<stepTree->GetEntries(); i++) {
      stepTree->GetEntry(i);
      // Loop over steps
      for (int n=0; n<NSteps; n++) {
	// Condition for a proton being in S4
	if (PreStepVolumeID[n]!=7 && PostStepVolumeID[n]==7) {
	  hMCVert->Fill(PostStepPosX[n]/1000.);
	  // If it passes a momentum cut add it
	  for (int c=0; c<momCuts.size(); c++) {
	    if (PostStepMom[n] > momCuts.at(c)) nProMC.at(c)++;
	  }
	}
      }
    } // Loop over entries

    for (int c=0; c<momCuts.size(); c++) {
      gMC->SetPoint(gMC->GetN(), momCuts.at(c), nProMC.at(c)/mcS3.at(b)); 
    }

    // Scale hists so the ratios match the values we usually see
    hDataVert->Scale(hDataVert->GetNbinsX() / (dataS3.at(b) * nSpills));
    hMCVert->Scale(hMCVert->GetNbinsX() / mcS3.at(b));

    gMC->SetLineColor(kRed);
    gMC->SetMarkerColor(kRed);
    gMC->SetLineWidth(2);
    gData->SetLineColor(kBlack);
    gData->SetMarkerColor(kBlack);
    gData->SetLineWidth(2);
    mgDataMCComp->Add(gMC);
    mgDataMCComp->Add(gData);

    hsDataMCCompVert->Add(hMCVert);
    hsDataMCCompVert->Add(hDataVert);

    TLegend *leg = new TLegend(0.2, 0.4, 0.65, 0.65);
    leg->AddEntry(gMC, "MC", "l");
    leg->AddEntry(gData, "Data", "l");

    fout->cd();
    leg->Write("leg");
    gData->Write();
    gMC->Write();
    mgDataMCComp->Write();
    hsDataMCCompVert->Write();
    hDataVert->Write();
    hMCVert->Write();

    mcFile->Close();
    delete mcFile;
  } // Loop over blocks

  fout->Close();
  dataFile->Close();
  delete fout;
  delete dataFile;
} // mcCutComp
