// intProtons.C
// Calculate total number of protons per spill across the data with no electron target

void intProtons(const char* saveDir,
		const char* ustofDir="/zfs_home/sjones/mylinktoutof/",
		const char* dstofDir="/scratch0/dbrailsf/temp/mylinktodtof/") 
{
  gROOT->SetBatch(kTRUE);

  // Define the runs to be used for varying number of blocks
  const char* str0Block = "Data_2018_8_31_b2_800MeV_0block.root";
  const char* str1Block = "Data_2018_9_1_b4_800MeV_1block_bend4cm.root";
  const char* str2Block = "Data_2018_9_1_b2_800MeV_2block_bend4cm.root";
  const char* str3Block = "Data_2018_9_1_b3_800MeV_3block_bend4cm.root";
  const char* str4Block = "Data_2018_9_1_b8_800MeV_4block_bend4cm.root";

  // Unix timestamp for start of electron target data
  const double eTargetT = 1536145800.;

  // Time of flight cuts for S1 to S3
  // Particles travelling at c should cross the distance in between about 36.9ns
  // and 35.6ns
  const double tLight = 35.8; // Expected time in ns of light speed particles according to AK
  // This is before the shift is applied
  const double proLow  = -12.5;
  const double proHi   = 60.;
  const double piLow = -29.06;
  const double piHi  = -28.06;
  // S3 amplitude cut for protons
  // Apply to A1ToF and A2ToF
  // AK's standard cut
  const double ACut = 0.25;
  // SJ's bar by bar amplitude cut (by eye)
  // For A1 
  const std::vector<double> A1CutVec = {0.25, 0.25, 0.275, 0.2, 0.25, 0.25, 0.25, 0.275, 0.325, 
					0.3, 0.3, 0.2, 0.2, 0.25, 0.225, 0.25, 0.25, 0.3, 0.3, 
					0.25, 0.3, 0.3};
  // For A2
  const std::vector<double> A2CutVec = {0.3, 0.275, 0.25, 0.125, 0.3, 0.3, 0.15, 0.225, 0.25, 0.25,
					0.225, 0.225, 0.2, 0.225, 0.225, 0.225, 0.2, 0.275, 0.25, 
					0.225, 0.3, 0.3};

  // Have each entry as a new spill
  TFile *fout = new TFile(Form("%s/intProtons.root", saveDir), "recreate");
  TTree *outTree = new TTree("outTree", "Particles in spill");
  int nP;
  int nPi;
  double spillTime;
  TString file;
  outTree->Branch("nP", &nP);
  outTree->Branch("nPi", &nPi);
  outTree->Branch("spillTime", &spillTime);
  outTree->Branch("file", &file);

  const char* ext   = ".root";
  const char* pref  = "Data";
  
  TString str;
  const char *entry;

  char *dir  = gSystem->ExpandPathName(ustofDir);
  void *dirp = gSystem->OpenDirectory(dir);

  vector<pair<double, TString> > fileVec;

  cout.precision(13);

  // Go through directory and open each utof file
  while (entry = (char*)gSystem->GetDirEntry(dirp)) {
    str = entry;
    if (str.EndsWith(ext) && str.Contains(pref) && !str.Contains("3block") && !str.Contains("2block") && !str.Contains("1block") && !str.Contains("0block")) {
      cout<<"Checking utof file "<<str<<endl;
      TString strTemp = str;
      str.Prepend(ustofDir);
      TFile *ustofFile = new TFile(str, "read");
      TTree *tree = (TTree*)ustofFile->Get("tree");
      TNamed *start = 0;
      TNamed *end   = 0;
      ustofFile->GetObject("start_of_run", start);
      ustofFile->GetObject("end_of_run", end);
      const char* startchar = start->GetTitle();
      const char* endchar   = end->GetTitle();
      string startstr(startchar);
      string endstr(endchar);
      string unixstart = startstr.substr(25,10);
      string unixend   = endstr.substr(23,10);
      double ustofStart = stod(unixstart);	 
      double ustofEnd   = stod(unixend);
      fileVec.push_back(make_pair(ustofStart, strTemp));
            
    } // if (str.EndsWith(ext) && str.Contains(pref) && str.Contains("4block"))
  } // while (entry = (char*)gSystem->GetDirEntry(dirp))
  sort(begin(fileVec), end(fileVec));
  cout<<"Files: "<<fileVec.size()<<endl;

  for (int i=0; i<fileVec.size(); i++){ 
    cout<<fileVec[i].second<<", "<<fileVec[i].first<<endl;
  }

  TGraph *grNonE_nP    = new TGraph();
  TGraph *grNonE_nPi   = new TGraph();
  TGraph *grNonE_ratio = new TGraph();
  TGraph *grAll_nP    = new TGraph();
  TGraph *grAll_nPi   = new TGraph();
  TGraph *grAll_ratio = new TGraph();

  // Now go through files and count the protons
  for (int i=0; i<fileVec.size(); i++) {
    TString path = fileVec[i].second;
    TString pathNoSuff = path.ReplaceAll(".root", "");
    TH1D *hFileTof = new TH1D(Form("hFileTof_%s", pathNoSuff.Data()), Form("%s time of flight", pathNoSuff.Data()), 300, -100, 200);
    path.Append(".root");
    path.Prepend(ustofDir);
    TFile *ustofFile = new TFile(path, "read");
    TTree *tree = (TTree*)ustofFile->Get("tree");
    // Go through file and count how many protons and pions in each spill
    double tSoSd;
    int nhit;
    int nBar[50];
    double tToF[50];
    float A1ToF[50];
    float A2ToF[50];
    double tS1;
    double tTrig;
    tree->SetBranchAddress("tSoSd", &tSoSd);
    tree->SetBranchAddress("tS1", &tS1);
    tree->SetBranchAddress("tTrig", &tTrig);
    tree->SetBranchAddress("tToF", tToF);
    tree->SetBranchAddress("nhit", &nhit);
    tree->SetBranchAddress("nBar", nBar);
    tree->SetBranchAddress("A1ToF", A1ToF);
    tree->SetBranchAddress("A2ToF", A2ToF);

    spillTime = fileVec[i].first;
    fout->cd();
    double lastSpill=0.;
    for (int t = 0; t < tree->GetEntries(); t++) {
      tree->GetEntry(t);
      file = pathNoSuff;
      if ((tSoSd - lastSpill) > 2e9) {
	cout<<"Spill at "<<spillTime<<", nP, nPi, "<<nP<<", "<<nPi<<endl;
	if (spillTime < eTargetT && spillTime >1534.8e6) {
	  grNonE_nP->SetPoint(grNonE_nP->GetN(), spillTime, nP);
	  if (nPi < 5000) {
	    grNonE_nPi->SetPoint(grNonE_nPi->GetN(), spillTime, nPi);
	  }
	  if (nPi!=0) {
	    grNonE_ratio->SetPoint(grNonE_ratio->GetN(), spillTime, (double)nP/(double)nPi);	  
	  } // if (nPi!=0)
	} // if (spillTime < eTargetT) 
	if (spillTime > 1534.8e6) {
	  grAll_nP->SetPoint(grAll_nP->GetN(), spillTime, nP);
	  if (nPi < 5000) {
	    grAll_nPi->SetPoint(grAll_nPi->GetN(), spillTime, nPi);
	  }
	  if (nPi!=0) {
	    grAll_ratio->SetPoint(grAll_ratio->GetN(),spillTime,(double)nP/(double)nPi);
	  } // if (nPi!=0)
	}
	outTree->Fill();
	lastSpill = tSoSd;
	spillTime = fileVec[i].first + (tSoSd/1e9);
	nPi = 0;
	nP  = 0;
      } // if ((tSoSd - lastSpill) < 2e9)

      for (int n = 0; n < nhit; n++) {
	double tofCalc = tToF[n] - tS1; //+ (tLight - (piLow + piHi)/2.);
	hFileTof->Fill(tofCalc);
	if (tofCalc > /*(tLight - (piLow+piHi)/2.) +*/ piLow && 
	    tofCalc < /*(tLight - (piLow+piHi)/2.) +*/ piHi &&
	    (tToF[n] - tSoSd) < 1e9) 
	  {
	    nPi++;
	  } // Is a pion
	else if (tofCalc > /*(tLight - (piLow + piHi)/2.) +*/ proLow && 
		 tofCalc < /*(tLight - (piLow + piHi)/2.) +*/ proHi && 
		 A1ToF[n] > A1CutVec[nBar[n]] && A2ToF[n] > A2CutVec[nBar[n]] &&
		 (tToF[n] - tSoSd) < 1e9) 
	  {
	    nP++;
	  } // Is a proton	    
      } // for (int n = 0; n < hit; n++)
    } // for (int t = 0; t < tree->GetEntries(); t++)

    hFileTof->Write();
    TCanvas *c1 = new TCanvas(Form("c1_%s", pathNoSuff.Data()));
    hFileTof->Draw("hist");
    c1->Print(Form("%s/%sTof.png",saveDir,pathNoSuff.Data())); 
    c1->Print(Form("%s/%sTof.pdf",saveDir,pathNoSuff.Data()));
    ustofFile->Close();
  }
  fout->cd();
  grNonE_nP->Write("grNonE_nP");
  grNonE_nPi->Write("grNonE_nPi");
  grNonE_ratio->Write("grNonE_ratio");
  grAll_nP->Write("grAll_nP");
  grAll_nPi->Write("grAll_nPi");
  grAll_ratio->Write("grAll_ratio");
  outTree->Write();
  cout<<outTree->GetEntries()<<" recorded"<<endl;
  fout->Close();

} // intProtons
