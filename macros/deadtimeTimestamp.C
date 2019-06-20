// Same thing but do it for a given utof file with the appropriate matching dtof file
void deadtimeTimestamp(const char* utofFile,	    
		       //	       const char* saveDir,
		       //const char* spillDB,
		       const char* dstofDir="/nfs/scratch0/dbrailsf/data_backup/dtof_backup/";
		       const char* ustofDir="/nfs/scratch0/dbrailsf/data_backup/utof_backup_firsthitpinnedtounixtime/Data_root_v3_wo_walk_corr/") 
{
  gROOT->SetBatch(kTRUE);

  gSystem->Exec(Form("mkdir /scratch0/sjones/plots/deadtime/%s/", utofFile));

  TFile *fout = new TFile(Form("/scratch0/sjones/plots/deadtime/%s/%s_out.root", utofFile, utofFile), "recreate");

  TH1D *hDeltatDtof = new TH1D(Form("hDeltaDtof_%s", utofFile), Form("S1 #cap S2 #deltat as measured in dtof: %s", utofFile), 300, 500, 1000000);
  TH1D *hDeltatUtof = new TH1D(Form("hDeltaUtof_%s", utofFile), Form("S1 #cap S2 #deltat as measured in utof: %s", utofFile), 300, 500, 1000000);
  TH1D *hDeltatDtofEarly = new TH1D(Form("hDeltatDtofEarly_%s", utofFile), Form("S1 #cap S2 #deltat measured in dtof (first 0.29s of spill): %s", utofFile), 300, 500, 1000000);
  TH1D *hDeltatUtofEarly = new TH1D(Form("hDeltatUtofEarly_%s", utofFile), Form("S1 #cap S2 #deltat measured in utof (first 0.29s of spill): %s", utofFile), 300, 500, 1000000);
  TH1D *hUtofInSpill = new TH1D(Form("hUtofInSpill_%s", utofFile), Form("S1 #cap S2 time since start of spill: %s", utofFile), 400, 0.1, .7);

  // Ustof-dstof cable delay
  const double ustofDelay = 184.7;
  // Apparent cut off time from spill profiles
  const double cutOffT = 0.29;

  double startTime = 0;
  double endTime   = 0;

  // Go through utof files and count the number of hits
  TFile *futof = new TFile(Form("%s/%s.root", ustofDir, utofFile), "read");

  double tToF[50];
  double tTrig;
  double tS1;
  double tSoSd;
  int nhit;
  int nBar[50];

  TTree *tree = (TTree*)futof->Get("tree");

  tree->SetBranchAddress("nhit", &nhit);
  tree->SetBranchAddress("tS1", &tS1);
  tree->SetBranchAddress("tToF", tToF);
  tree->SetBranchAddress("tTrig", &tTrig);
  tree->SetBranchAddress("tSoSd", &tSoSd);
  tree->SetBranchAddress("nBar", nBar);

  double lastS1S2utof = 0.;

  double lastUtofSpill = 0.;
  int nUtofSpills = 0;
  int nDtofSpills = 0;

  std::vector<double> utofTimes;
  std::vector<double> dtofTimes;

  // Get the start and end time of the file
  TNamed *start = 0;
  TNamed *end   = 0;
  futof->GetObject("start_of_run", start);
  futof->GetObject("end_of_run", end);

  const char* startchar = start->GetTitle();
  std::string startstr(startchar);
  std::string unixstart = startstr.substr(25,10);
  startTime = stoi(unixstart);
  
  const char* endchar = end->GetTitle();
  std::string endstr(endchar);
  std::string unixend = endstr.substr(23,10);
  endTime = stoi(unixend);

  cout.precision(13);
  cout<<"Utof file start, end "<<startTime<<", "<<endTime<<endl;

  Int_t runMin=-1;
  Int_t runMax=-1;

  // Find the correct dtof files
  for (int irun=900; irun<1400; irun++) {
    TFile *fin = new TFile(Form("%s/run%d/DsTOFcoincidenceRun%d_tdc1.root", dstofDir, irun, irun), "read");
    RawDsTofCoincidence *tofCoinTemp = NULL;
    TTree *tree = (TTree*) fin->Get("tofCoinTree");
    tree->SetBranchAddress("tofCoin", &tofCoinTemp);
    tree->GetEntry(0);
    UInt_t firstTemp = tofCoinTemp->unixTime[0];
    tree->GetEntry(tree->GetEntries()-1);
    UInt_t lastTemp = tofCoinTemp->unixTime[0];
      
    fin->Close();
    delete fin;
      
    if (firstTemp>endTime){
      break;
    }
      
    if (firstTemp<startTime && lastTemp>startTime){
      runMin = irun;
    }
      
    if (firstTemp<endTime && lastTemp>endTime){
      runMax = irun;
    }   
  } // for (int irun=950; irun<1400; irun++) 
    
  cout << "Min and max runs are " << runMin << " " << runMax << endl;
  
  std::vector<int> nS1S2dtofVec;
  int nS1S2dtof      = 0;
  int nS1S2dtofEarly = 0;
  // Now go over the spill DB files for these runs and 
  // make two vectors of the appropriate matching spill times
  int spillCount = 0;

  TTree *dtofSpillTree = new TTree("dtofSpillTree", "Dtof spills");
  TTree *utofSpillTree = new TTree("utofSpillTree", "Utof spills");
  int dtofHits;
  int dtofHitsEarly;
  dtofSpillTree->Branch("dtofHits", &dtofHits);
  dtofSpillTree->Branch("dtofHitsEarly", &dtofHitsEarly);
  int utofHits;
  int utofHitsEarly;
  utofSpillTree->Branch("utofHits", &utofHits);
  utofSpillTree->Branch("utofHitsEarly", &utofHitsEarly);

  for (int irun = runMin; irun < runMax+1; irun++) {
    cout<<"Getting hits for dtof run "<<irun<<endl;

    TFile *dbFile = new TFile(Form("/scratch0/sjones/spillDB/spillDB_run%d_run%d.root", irun, irun), "read");
    std::vector<double> tempDtofTimes;
    TTree *spillTree = (TTree*)dbFile->Get("spillTree");
    double globalSpillTime;
    double ustofSpillTime;
    spillTree->SetBranchAddress("globalSpillTime", &globalSpillTime);
    spillTree->SetBranchAddress("ustofSpillTime", &ustofSpillTime);
    // Loop over these entries
    // Only add them to the vectors if they are between
    // correct unix start and end times
    for (int t = 0; t < spillTree->GetEntries(); t++) {
      spillTree->GetEntry(t);

      if (globalSpillTime < startTime) continue;
      if (globalSpillTime > endTime) break;

      utofTimes.push_back(ustofSpillTime);
      dtofTimes.push_back(globalSpillTime);
      tempDtofTimes.push_back(globalSpillTime);
    } // for (int t = 0; t < spillTree->GetEntries(); t++)
    dbFile->Close();
    // Now open the corresponding dtof file and count the appropriate spills
    // Utof and beam signals go into both TDCs so only need to run over one
    TFile *tofFile = new TFile(Form("%s/run%d/DsTOFtreeRun%d_tdc1.root", dstofDir, irun, irun));
    RawDsTofHeader *tof = NULL;
    TTree *tofTree = (TTree*)tofFile->Get("tofTree");
    tofTree->SetBranchAddress("tof", &tof);
    tofTree->GetEntry(0);
    UInt_t firstTime = tof->unixTime;
    vector<int> nS1S2dtofTemp;
    nS1S2dtofTemp.resize(tempDtofTimes.size(), 0);
    // Loop over only the spill times in this particular dtof file
    // Saves us having to open loads of dtof files to check if the spill is there
    double lastS1S2Dtof = 0.;
    double lastS1S2DtofEarly = 0.;
    int lastdt = 0;

    for (int spill = 0; spill < tempDtofTimes.size(); spill++) {
      dtofHits = 0;
      dtofHitsEarly = 0;
      TGraph *grTimeSinceSpillDtof = new TGraph();
      for (int t = lastdt; t < tofTree->GetEntries(); t++) {
	tofTree->GetEntry(t);
	if ((tof->fakeTimeNs/1e9 - tempDtofTimes[spill] + firstTime) > 1.) break;
	if ((tof->fakeTimeNs/1e9 - tempDtofTimes[spill] + firstTime) < 0.) continue;
	// Check that hit is up to 1s after the chosen spill
	if (tof->channel == 13 && 
	    (tof->fakeTimeNs/1e9 - tempDtofTimes[spill] + firstTime) < 1. &&
	    (tof->fakeTimeNs/1e9 - tempDtofTimes[spill] + firstTime) > 0. &&
	    (tof->fakeTimeNs - lastS1S2Dtof) > 500.) {
	  nS1S2dtofTemp[spill]++;
	  nS1S2dtof++;
	  dtofHits++;
	  hDeltatDtof->Fill(tof->fakeTimeNs - lastS1S2Dtof);
	  lastS1S2Dtof = tof->fakeTimeNs;
	  grTimeSinceSpillDtof->SetPoint(grTimeSinceSpillDtof->GetN(), tof->fakeTimeNs/1e9 + firstTime - tempDtofTimes[spill], grTimeSinceSpillDtof->GetN());
	  lastdt = t;
	  if ((tof->fakeTimeNs/1e9 - tempDtofTimes[spill] + firstTime) < .29) {
	    nS1S2dtofEarly++;
	    dtofHitsEarly++;
	    hDeltatDtofEarly->Fill(tof->fakeTimeNs - lastS1S2DtofEarly);
	    lastS1S2DtofEarly = tof->fakeTimeNs;
	  }
	} // Is an S1S2 coincidence
      } // for (int t = 0; t < tofTree->GetEntries(); t++)
      fout->cd();
      grTimeSinceSpillDtof->SetTitle(Form("Event number versus time since spill start in dtof, spill %d; Time since spill start / s; Event no.", spillCount));
      grTimeSinceSpillDtof->Write(Form("grTimeSinceSpillDtof%d", spillCount));
      spillCount++;
      dtofSpillTree->Fill();
    } // for (int spill = 0; spill < tempDtofTimes.size(); spill++)
    // Add this to the larger spill vector 

    for (int i = 0; i < nS1S2dtofTemp.size(); i++) {
      nS1S2dtofVec.push_back(nS1S2dtofTemp[i]);
    }

    delete tof;
    tofFile->Close();
  } // for (int irun = runMin; irun < runMax+1; irun++)

  int nSpills = utofTimes.size();
  cout<<"Counting over "<<nSpills<<" matched spills"<<endl;

  int nS1S2utof = 0;
  int nS1S2utofEarly = 0;
  std::vector<int> nS1S2utofVec;
  nS1S2utofVec.resize(nSpills, 0);
  //nS1S2dtofVec.resize(nSpills, 0);
  // Count number of S1 x S2 hits
  // Only do this for the spills in the spillDB
  int lastut=0;
  double lastS1S2utofEarly = 0.;
  for (int spill = 0; spill < nSpills; spill++) {
    TGraph *grTimeSinceSpill = new TGraph();
    double tempUtofSpill = utofTimes[spill];
    double tempDtofSpill = dtofTimes[spill];
    utofHits = 0;
    utofHitsEarly = 0;
    for (int t=lastut; t<tree->GetEntries(); t++ ) {
      tree->GetEntry(t);
      if ((tSoSd/1e9 + startTime) > tempUtofSpill) break;

      if (tTrig !=0 && (tSoSd/1e9 + startTime) == tempUtofSpill) {
	nS1S2utof++;
	nS1S2utofVec[spill]++;
	utofHits++;
	hDeltatUtof->Fill(tTrig - lastS1S2utof);
	hUtofInSpill->Fill(tTrig/1e9 + startTime - tempUtofSpill);
	grTimeSinceSpill->SetPoint(grTimeSinceSpill->GetN(), tTrig/1e9 + startTime - tempUtofSpill, grTimeSinceSpill->GetN());
	lastS1S2utof = tTrig;
	lastut=t;
	if ((tTrig - tSoSd) < 0.29e9) {
	  nS1S2utofEarly++;
	  utofHitsEarly++;
	  hDeltatUtofEarly->Fill(tTrig - lastS1S2utofEarly);
	  lastS1S2utofEarly = tTrig;
	} // if ((tTrig - tSoSd) < 0.29e9)
      }
    }
    fout->cd();
    //cout<<"Utof spill at "<<tempUtofSpill<<", hits "<<nS1S2utofVec[spill]<<", Dtof spill at "<<tempDtofSpill<<", hits "<<nS1S2dtofVec[spill]<<endl;
    grTimeSinceSpill->SetTitle(Form("Event number versus time since spill start, spill %d; Time since spill start / s; Event no.", spill));
    grTimeSinceSpill->Write(Form("grTimeSinceSpill%d",spill));
    utofSpillTree->Fill();
  } // for (int spill = 0; spill < nSpills; spill++) 
  
  cout<<"Dtof S1 S2 "<<nS1S2dtof<<", Utof S1 S2 "<<nS1S2utof<<", Ratio "<<(double)nS1S2utof/(double)nS1S2dtof<<endl;
  cout<<"Early hits, Dtof S1 S2 "<<nS1S2dtofEarly<<", Utof S1 S2 "<<nS1S2utofEarly<<", Ratio "<<(double)nS1S2utofEarly/(double)nS1S2dtofEarly<<endl;

  // cout<<"There are "<<nDtofSpills<<" dtof spills, "<<nUtofSpills<<" utof spills"<<endl;
  fout->cd();
  gStyle->SetOptFit(1);
  TCanvas *cdtof = new TCanvas("cdtof");
  hDeltatDtof->Write();
  TF1 *f1 = new TF1("f1", "exp([0]+[1]*x)", 500, 1000000);
  hDeltatDtof->Fit("f1", "R");
  hDeltatDtof->Draw("hist");
  f1->Draw("same");
  cdtof->SetGridx();
  cdtof->SetGridy();
  cdtof->Print(Form("/scratch0/sjones/plots/deadtime/%s/%s_deltatDtof.png", utofFile, utofFile));
  cdtof->Print(Form("/scratch0/sjones/plots/deadtime/%s/%s_deltatDtof.pdf", utofFile, utofFile));

  TCanvas *cutof = new TCanvas("cutof");
  hDeltatUtof->Write();
  TF1 *f2 = new TF1("f2", "exp([0]+[1]*x)", 40000, 1000000);
  hDeltatUtof->Fit("f2", "R");
  cout<<"Utof integral below 40us "<<f2->Integral(0, 40000)/3200<<endl;
  cout<<"Utof integral above 40us (TF1) "<<f2->Integral(40000, 1000000)/3200<<endl;
  cout<<"Utof integral above 40us (hist) "<<hDeltatUtof->Integral(12, 300)<<endl;
  hDeltatUtof->Draw("hist");
  f2->Draw("same");
  cutof->SetGridx();
  cutof->SetGridy();
  cutof->Print(Form("/scratch0/sjones/plots/deadtime/%s/%s_deltatUtof.png", utofFile, utofFile));
  cutof->Print(Form("/scratch0/sjones/plots/deadtime/%s/%s_deltatUtof.pdf", utofFile, utofFile));

  TCanvas *cdtofE = new TCanvas("cdtofE");
  hDeltatDtofEarly->Write();
  TF1 *f3 = new TF1("f3", "exp([0]+[1]*x)", 500, 1000000);
  hDeltatDtofEarly->Fit("f3", "R");
  hDeltatDtofEarly->Draw("hist");
  f3->Draw("same");
  cdtofE->SetGridx();
  cdtofE->SetGridy();
  cdtofE->Print(Form("/scratch0/sjones/plots/deadtime/%s/%s_deltatDtofEarly.png", utofFile, utofFile));
  cdtofE->Print(Form("/scratch0/sjones/plots/deadtime/%s/%s_deltatDtofEarly.pdf", utofFile, utofFile));

  TCanvas *cutofE = new TCanvas("cutofE");
  hDeltatUtofEarly->Write();
  TF1 *f4 = new TF1("f4", "exp([0]+[1]*x)", 40000, 1000000);
  hDeltatUtofEarly->Fit("f4", "R");
  hDeltatUtofEarly->Draw("hist");
  f4->Draw("same");
  cutofE->SetGridx();
  cutofE->SetGridy();
  cutofE->Print(Form("/scratch0/sjones/plots/deadtime/%s/%s_deltatUtofEarly.png", utofFile, utofFile));
  cutofE->Print(Form("/scratch0/sjones/plots/deadtime/%s/%s_deltatUtofEarly.pdf", utofFile, utofFile));

  TCanvas *cComb = new TCanvas("cComb");
  cComb->SetLogy();
  hDeltatDtof->SetLineColor(kRed);
  hDeltatDtof->Draw("hist");
  hDeltatUtof->Draw("hist same");
  cComb->SetGridx();
  cComb->SetGridy();
  cComb->Print(Form("/scratch0/sjones/plots/deadtime/%s/%s_deltatComb.png", utofFile, utofFile));
  cComb->Print(Form("/scratch0/sjones/plots/deadtime/%s/%s_deltatComb.pdf", utofFile, utofFile));

  TGraph *grRatioDtof = new TGraph();
  TGraph *grRatioDtofFit = new TGraph();
  TGraph *grUtofDtof  = new TGraph();
  /*
  TTree *s1s2Tree = new TTree("s1s2Tree", "S1S2 Hits");
  s1s2Tree->SetDirectory(0);
  double utof;
  double dtof;
  s1s2Tree->Branch("utof", &utof);
  s1s2Tree->Branch("dtof", &dtof);
  */
  for (int n=0; n < nSpills; n++) {
    grUtofDtof->SetPoint(grRatioDtof->GetN(), nS1S2dtofVec[n], nS1S2utofVec[n]);
    if (nS1S2dtofVec[n] != 0) {
      grRatioDtof->SetPoint(grRatioDtof->GetN(), nS1S2dtofVec[n], (double)nS1S2utofVec[n]/(double)nS1S2dtofVec[n]);
      // Is in the region we want to fit
      if ((((double)nS1S2utofVec[n]/(double)nS1S2dtofVec[n]) > 0.06) &&
	  (nS1S2dtofVec[n] > 300)) {
	grRatioDtofFit->SetPoint(grRatioDtofFit->GetN(), nS1S2dtofVec[n], (double)nS1S2utofVec[n]/(double)nS1S2dtofVec[n]);
      } // if (((double)nS1S2utofVec[n]/(double)nS1S2dtofVec[n]) > 0.05)
    } // if (nS1S2dtofVec[n] != 0) 
  } // for (int n=0; n < nSpills; n++)

  TCanvas *cRatioDtof = new TCanvas("cRatioDtof");
  grRatioDtof->SetTitle(Form("Utof/Dtof S1 #cap S2: %s; S1 #cap S2 dtof; S1 #cap S2 utof/dtof", utofFile));
  grRatioDtof->Draw("AP*");
  cRatioDtof->Print(Form("/scratch0/sjones/plots/deadtime/%s/%s_RatioDtof.png", utofFile, utofFile));
  cRatioDtof->Print(Form("/scratch0/sjones/plots/deadtime/%s/%s_RatioDtof.pdf", utofFile, utofFile));

  TF1 *fFit = new TF1("fFit", "[0]+[1]*x", 300., 3000.);
  grRatioDtofFit->Fit("fFit", "R");

  TCanvas *cUtofDtof = new TCanvas("cUtofDtof");
  grUtofDtof->SetTitle(Form("Utof vs. Dtof S1 #cap S2: %s; S1 #cap S2 dtof; S1 #cap S2 utof", utofFile));
  grUtofDtof->Draw("AP*");
  cUtofDtof->Print(Form("/scratch0/sjones/plots/deadtime/%s/%s_UtofDtof.png", utofFile, utofFile));
  cUtofDtof->Print(Form("/scratch0/sjones/plots/deadtime/%s/%s_UtofDtof.pdf", utofFile, utofFile));

  TCanvas *cUtofInSpill = new TCanvas("cUtofInSpill");
  hUtofInSpill->Draw("hist");
  cUtofInSpill->Print(Form("/scratch0/sjones/plots/deadtime/%s/%s_UtofInSpill.png", utofFile, utofFile));
  cUtofInSpill->Print(Form("/scratch0/sjones/plots/deadtime/%s/%s_UtofInSpill.pdf", utofFile, utofFile));

  fout->cd();  
  TVectorD rms(1);
  rms[0] = grRatioDtofFit->GetRMS(2);
  rms.Write("rms");
  TVectorD slope(1);
  slope[0] = fFit->GetParameter(1);
  slope.Write("slope");
  TVectorD intercept(1);
  intercept[0] = fFit->GetParameter(0);
  intercept.Write("intercept");
  TObject writeS1S2dtof;
  writeS1S2dtof.SetUniqueID(nS1S2dtof);
  writeS1S2dtof.Write("nS1S2dtof");
  TObject writeS1S2utof;
  writeS1S2utof.SetUniqueID(nS1S2utof);  
  writeS1S2utof.Write("nS1S2utof");
  
  start->Write("start_of_run");
  end->Write("end_of_run");
  utofSpillTree->Write();
  dtofSpillTree->Write();
  hDeltatDtof->Write();
  hDeltatUtof->Write();
  hUtofInSpill->Write();
  grRatioDtof->Write("grRatioDtof");
  grRatioDtofFit->Write("grRatioDtofFit");
  grUtofDtof->Write("grUtofDtof");
  fout->Close();
  delete fout;
} // deadtimeTimestamp
