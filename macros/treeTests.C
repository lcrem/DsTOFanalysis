// treeTests.C

// spillDiff
// Calculate minimum time difference between closest spills in utof and dtof
void spillDiff(const char* ustofDB, const char* dstofDB) {
  TFile *usIn   = new TFile(ustofDB, "read");
  TFile *dsIn   = new TFile(dstofDB, "read");
  TTree *usTree = (TTree*)usIn->Get("utofBeamTree");
  TTree *dsTree = (TTree*)dsIn->Get("dtofBeamTree");
  int unixRunStartUs;
  double nsTimeUs;
  double fakeUnixTimeUs;
  int unixRunStartDs;
  double nsTimeDs;
  double fakeUnixTimeDs;
  usTree->SetBranchAddress("unixRunStart", &unixRunStartUs);
  usTree->SetBranchAddress("nsTime", &nsTimeUs);
  usTree->SetBranchAddress("fakeUnixTime", &fakeUnixTimeUs);
  dsTree->SetBranchAddress("unixRunStart", &unixRunStartDs);
  dsTree->SetBranchAddress("nsTime", &nsTimeDs);
  dsTree->SetBranchAddress("fakeUnixTime", &fakeUnixTimeDs);

  TH1D *htDiffds = new TH1D("htDiffds", "Time between successive spills; Time / s", 300, 0, 50);
  TH1D *htDiffus = new TH1D("htDiffus", "Time between successive spills; Time / s", 300, 0, 50);
  TH1D *hspillSepDs  = new TH1D("hspillSepDs", "#Deltat between dstof spill & closest ustof spill; Time / s", 300, -50, 50);
  TH2D *h2spillSepDs = new TH2D("h2spillSepDs", "#Deltat between dstof spill & closest ustof spill; Run time / s; #Deltat / s", 250, 1534.85e6, 1537.195e6, 250, -50, 50);
  TH1D *hspillSepUs  = new TH1D("hspillSepUs", "#Deltat between ustof spill & closest dstof spill; Time / s", 300, -50, 50);
  TH2D *h2spillSepUs = new TH2D("h2spillSepUs", "#Deltat between ustof spill & closest dstof spill; Run time / s; #Deltat / s", 250, 1534.85e6, 1537.195e6, 250, -50, 50);
  dsTree->GetEntry(0);
  usTree->GetEntry(0);
  double lastSpillus = fakeUnixTimeUs;
  double lastSpillds = fakeUnixTimeDs;
  TH1D *hunmatchedDs = new TH1D("hunmatchedDs", "Run time of unmatched spills DsToF spills; Run time / s", 250, 1534.85e6, 1537.195e6);
  TH1D *hunmatchedUs = new TH1D("hunmatchedUs", "Run time of unmatched spills UsToF spills; Run time / s", 250, 1534.85e6, 1537.195e6);
  TGraph *grDs = new TGraph();
  TGraph *grUs = new TGraph();
  
  for(int t=1; t<dsTree->GetEntries(); t++) {
    double lowestDiff = 99999999999.;
    if(t%1000 == 0) cout<<"Dstof spill "<<t<<endl;
    dsTree->GetEntry(t);
    grDs->SetPoint(grDs->GetN(), fakeUnixTimeDs, grDs->GetN());
    htDiffds->Fill(fakeUnixTimeDs-lastSpillds);
    lastSpillds = fakeUnixTimeDs;
    for(int u=0; u<usTree->GetEntries(); u++) {
      usTree->GetEntry(u);
      double diff = fakeUnixTimeDs - fakeUnixTimeUs;
      if (abs(diff) < abs(lowestDiff)) {
	lowestDiff = diff;
      }
    }
    hspillSepDs->Fill(lowestDiff);
    h2spillSepDs->Fill(fakeUnixTimeDs, lowestDiff);
    // Spills are unmatched
    if (abs(lowestDiff) > 2) {
      hunmatchedDs->Fill(fakeUnixTimeDs);
    }
  }
  
  for(int t=1; t<usTree->GetEntries(); t++) {
    double lowestDiff2 = 99999999999.;
    if(t%1000 == 0) cout<<"Ustof spill "<<t<<endl;
    usTree->GetEntry(t);    
    grUs->SetPoint(grUs->GetN(), fakeUnixTimeUs, grUs->GetN());
    htDiffus->Fill(fakeUnixTimeUs-lastSpillus);
    lastSpillus = fakeUnixTimeUs;
    for(int u=0; u<dsTree->GetEntries(); u++) {
      dsTree->GetEntry(u);
      double diff2 = fakeUnixTimeUs - fakeUnixTimeDs;
      if (abs(diff2) < abs(lowestDiff2)) {
	lowestDiff2 = diff2;
      }
    }
    hspillSepUs->Fill(lowestDiff2);
    h2spillSepUs->Fill(fakeUnixTimeUs, lowestDiff2);
    // Spills are unmatched
    if (abs(lowestDiff2) > 2) {
      hunmatchedUs->Fill(fakeUnixTimeUs);
    }
  }
  
  TCanvas *c1 = new TCanvas("c1");
  htDiffds->SetLineColor(kRed);
  htDiffds->SetFillColor(kRed);
  htDiffds->SetFillStyle(3004);
  htDiffus->SetLineColor(kBlue);
  htDiffus->SetFillColor(kBlue);
  htDiffus->SetFillStyle(3005);
  htDiffds->Draw("hist");
  htDiffus->Draw("hist same");
  TLegend * leg1 = new TLegend(0.6, 0.6, 0.8, 0.8);
  leg1->AddEntry(htDiffds, "DsToF", "f");
  leg1->AddEntry(htDiffus, "UsToF", "f");
  leg1->Draw();
  c1->Print("spillSepBoth.png");
  c1->Print("spillSepBoth.pdf");

  TCanvas *c2 = new TCanvas("c2");
  hspillSepDs->Draw("hist");
  c2->Print("deltaT_Ds.png");
  c2->Print("deltaT_Ds.pdf");
  
  gStyle->SetPalette(55);
  TCanvas *c3 = new TCanvas("c3");
  h2spillSepDs->Draw("colz");
  c3->Print("2dspillSep_Ds.png");
  c3->Print("2dspillSep_Ds.pdf");

  TCanvas *c4 = new TCanvas("c4");
  hunmatchedDs->Draw("hist");
  c4->Print("unmatchedSpillTimes_Ds.png");
  c4->Print("unmatchedSpillTimes_Ds.pdf");

  TCanvas *c2_1 = new TCanvas("c2_1");
  hspillSepUs->Draw("hist");
  c2_1->Print("deltaT_Us.png");
  c2_1->Print("deltaT_Us.pdf");
  
  gStyle->SetPalette(55);
  TCanvas *c3_1 = new TCanvas("c3_1");
  h2spillSepUs->Draw("colz");
  c3_1->Print("2dspillSep_Us.png");
  c3_1->Print("2dspillSep_Us.pdf");

  TCanvas *c4_1 = new TCanvas("c4_1");
  hunmatchedUs->Draw("hist");
  c4_1->Print("unmatchedSpillTimes_Us.png");
  c4_1->Print("unmatchedSpillTimes_Us.pdf");
  
  TCanvas *c5 = new TCanvas("c5");
  TMultiGraph *mg = new TMultiGraph();
  grUs->SetMarkerStyle(6);
  grDs->SetMarkerStyle(6);
  grUs->SetMarkerColor(kBlue);
  grDs->SetMarkerColor(kRed);
  mg->SetTitle("Beam spill number in ToFs over time; Time / s; Spill number");
  mg->Add(grUs);
  mg->Add(grDs);
  mg->Draw("ap");
  TLegend *leg2 = new TLegend(.2, .6, .4, .8);
  leg2->AddEntry(grUs, "UsToF", "p");
  leg2->AddEntry(grDs, "DsToF", "p");
  leg2->Draw();
  c5->Print("spillNumberComp.png");
  c5->Print("spillNumberComp.pdf");
} // spillDiff

// dstofComp
void dstofComp(const char* us, const char* dstofDB) {
  // Ustof file
  TFile *usIn = new TFile(us, "read");
  TTree *ustree = (TTree*)usIn->Get("tree");
  double tSoSd;
  ustree->SetBranchAddress("tSoSd", &tSoSd);
  TNamed *start = 0;
  TNamed *end   = 0;
  usIn->GetObject("start_of_run", start);
  usIn->GetObject("end_of_run", end);
  
  const char* startchar = start->GetTitle();
  std::string startstr(startchar);
  std::string unixstart = startstr.substr(25,10);
  const int ustofStart = stoi(unixstart);
  
  const char* endchar = end->GetTitle();
  std::string endstr(endchar);
  std::string unixend = endstr.substr(23,10);
  const int ustofEnd = stoi(unixend);

  cout<<"Ustof start, end "<<ustofStart<<", "<<ustofEnd<<endl;
  // Dstof file
  TFile *dstofIn = new TFile(dstofDB, "read");
  TTree *dstofTree = (TTree*)dstofIn->Get("beamTree");
  double fakeUnixNs;
  dstofTree->SetBranchAddress("fakeUnixNs", &fakeUnixNs);

  TH1D *hdstof = new TH1D("hdstof","", 1200, 1536490758, 1536491758);
  TH1D *hustof = new TH1D("hustof","", 1200, 1536490758, 1536491758);
  
  vector<double> dstofvec;
  double dstofSpill;
  for(int t=0; t<dstofTree->GetEntries(); t++) {
    dstofTree->GetEntry(t);
    if(fakeUnixNs>=ustofStart && fakeUnixNs<=ustofEnd) {
      dstofSpill = fakeUnixNs;
      dstofvec.push_back(dstofSpill);
      hdstof->Fill(dstofSpill);
    }
  }
  vector<double>ustofvec;
  double secSpill=0.;
  double lastSpill=0.;
  for (int t=0; t<ustree->GetEntries(); t++) {
    ustree->GetEntry(t);
    if(tSoSd != lastSpill && tSoSd - lastSpill > 1e9) {
      lastSpill = tSoSd;
      double secSpill = (lastSpill / 1e9) + ustofStart;
      ustofvec.push_back(secSpill);
      hustof->Fill(secSpill);
    }
  }
  cout<<"Ustof, dstof "<<ustofvec.size()<<", "<<dstofvec.size()<<endl;

 TCanvas *c1 = new TCanvas("c1");
 hustof->SetLineColor(kRed);
 hustof->Draw("hist");
 hdstof->Draw("hist same");
}
// utofDouble
void utofDouble() {
  TFile *utofFile = new TFile("utofBeamSpills.root", "read");
  TTree *tree = (TTree*)utofFile->Get("utofBeamTree");

  string fileName;
  double fakeUnixNs;
  int unixRunStart;
  double nsTime;
  
  tree->SetBranchAddress("fileName", &fileName);
  tree->SetBranchAddress("fakeUnixNs", &fakeUnixNs);
  tree->SetBranchAddress("unixRunStart", &unixRunStart);
  tree->SetBranchAddress("nsTime", &nsTime);

  int lastSpill = 0;
  for (int t=0; t<tree->GetEntries(); t++) {
    tree->GetEntry(t);
    if(fakeUnixNs > 1535.95e6 && fakeUnixNs < 1536.025e6 && unixRunStart != lastSpill) {
      cout<<unixRunStart<<", "<<fakeUnixNs<<", Entry: "<<t<<endl;
      lastSpill = unixRunStart;
      
    }
  }
}
// signalSep
void signalSep() {
  TFile *utofFile = new TFile("beamSpills.root", "read");
  TTree *tree = (TTree*)utofFile->Get("beamTree");

  double fakeUnixNs;
  int unixRunStart;
  double nsTime;
  
  tree->SetBranchAddress("fakeUnixNs", &fakeUnixNs);
  tree->SetBranchAddress("unixRunStart", &unixRunStart);
  tree->SetBranchAddress("nsTime", &nsTime);

  TCanvas *c1 = new TCanvas("c1");
  TH1D *hdiff = new TH1D("hdiff", "Time separation between successive beam signals in DsToF; Time / s; Signals", 150, 0, 40);


  double lastSpill = 0.;
  for(int t=0; t<tree->GetEntries(); t++) {
    tree->GetEntry(t);
    hdiff->Fill(fakeUnixNs - lastSpill);
    lastSpill = fakeUnixNs;
  }
  hdiff->Draw("hist");
  c1->Print("spillSep.png");
  c1->Print("spillSep.pdf");
}

void overlapComp(const char* us1, const char* us2, const char* dstofDB) {
  // Dstof files
  TFile *dstofIn = new TFile(dstofDB, "read");
  TTree *dstofTree = (TTree*)dstofIn->Get("beamTree");
  int run;
  int unixRunStart;
  double nsTime;
  double fakeUnixNs;
  dstofTree->SetBranchAddress("run", &run);
  dstofTree->SetBranchAddress("unixRunStart", &unixRunStart);
  dstofTree->SetBranchAddress("nsTime", &nsTime);
  dstofTree->SetBranchAddress("fakeUnixNs", &fakeUnixNs);
  
  // Ustof files
  TFile *in1 = new TFile(us1, "read");
  TFile *in2 = new TFile(us2, "read");

  TTree *tree1 = (TTree*)in1->Get("tree");
  TTree *tree2 = (TTree*)in2->Get("tree");

  double tSoSd1;
  double tSoSd2;
  
  //  TFile *spillfile = new TFile("utofBeamSpill.root", "read");
  //  TTree *beamTree = (TTree*)spillfile->Get("utofBeamTree");

  tree1->SetBranchAddress("tSoSd", &tSoSd1);
  tree2->SetBranchAddress("tSoSd", &tSoSd2);

  TNamed *start1 = 0;
  TNamed *end1   = 0;
  in1->GetObject("start_of_run", start1);
  in1->GetObject("end_of_run", end1);

  TNamed *start2 = 0;
  TNamed *end2   = 0;
  in2->GetObject("start_of_run", start2);
  in2->GetObject("end_of_run", end2);

  const char* startchar1 = start1->GetTitle();
  std::string startstr1(startchar1);
  std::string unixstart1 = startstr1.substr(25,10);
  const int ustofStart1 = stoi(unixstart1);
  
  const char* endchar1 = end1->GetTitle();
  std::string endstr1(endchar1);
  std::string unixend1 = endstr1.substr(23,10);
  const int ustofEnd1 = stoi(unixend1);

  const char* startchar2 = start2->GetTitle();
  std::string startstr2(startchar2);
  std::string unixstart2 = startstr2.substr(25,20);
  const int ustofStart2 = stoi(unixstart2);
  
  const char* endchar2 = end2->GetTitle();
  std::string endstr2(endchar2);
  std::string unixend2 = endstr2.substr(23,20);
  const int ustofEnd2 = stoi(unixend2);

  vector<double> spills1;
  vector<double> spills2;

  cout<<"File 1 start, end "<<ustofStart1<<", "<<ustofEnd1<<endl;
  cout<<"File 2 start, end "<<ustofStart2<<", "<<ustofEnd2<<endl;
  // Make two vectors of beam spills in the windown
  double lastSpill1 = 0.; // ins
  for (int t=0; t<tree1->GetEntries(); t++) {
    tree1->GetEntry(t);
    if(tSoSd1 != lastSpill1 && tSoSd1 - lastSpill1 > 1e9 && ustofStart1 + (tSoSd1 / 1e9) >= ustofStart2 && ustofStart1 +(tSoSd1 / 1e9) <= ustofEnd2) {
      lastSpill1 = tSoSd1;
      double secSpill1 = (lastSpill1 / 1e9) + ustofStart1;
      spills1.push_back(secSpill1);
    }
  }

  double lastSpill2 = 0.; // in s
  for (int t=0; t<tree2->GetEntries(); t++) {
    tree2->GetEntry(t);
    if(tSoSd2 != lastSpill2 && tSoSd2 - lastSpill2 > 1e9 && ustofStart2 + (tSoSd2 / 1e9) >= ustofStart1 && ustofStart2 +(tSoSd2 / 1e9) <= ustofEnd1) {
      lastSpill2 = tSoSd2;
      double secSpill2 = (lastSpill2 / 1e9) + ustofStart2;
      spills2.push_back(secSpill2);
    }
  }

  vector<double> spillsDstof;
  // Dstof beam tree
  for (int t=0; t<dstofTree->GetEntries(); t++) {
    dstofTree->GetEntry(t);
    if (fakeUnixNs>=ustofStart1 && fakeUnixNs>=ustofStart2 && fakeUnixNs<=ustofEnd1 && fakeUnixNs<=ustofEnd2) {
      double dstofSpill = fakeUnixNs;
      spillsDstof.push_back(dstofSpill);
    }
  }
  
  cout<<"File 1 "<<spills1.size()<<", File 2 "<<spills2.size()<<", Dstof "<<spillsDstof.size()<<endl;
  TH1D *hfile1 = new TH1D("file1", "file 1 beam signals; Time", 250, 1536821499, 1536821699);
  TH1D *hfile2 = new TH1D("file2", "file 2 beam signals; Time", 250, 1536821499, 1536821699);
  TH1D *hdstof = new TH1D("hdstof", "dstof beam signals; Time", 250, 1536821499, 1536821699);
  TH1D *hdiff1 = new TH1D("hdiff1", "Time between signals; Time / s", 40, 0, 40);
  TH1D *hdiff2 = new TH1D("hdiff2", "file 2 time between signals; Time / s", 40, 0, 40);
  TH1D *hdiffdstof = new TH1D("hdiffdstof", "dstof time between signals; Time / s", 40, 0, 40);
  
  const double firstspill1 = spills1[0];
  const double firstspill2 = spills2[0];
  const double firstSpillDstof = spillsDstof[0];
  cout.precision(17);
  cout<<"First spill1, spill2, dstof "<<firstspill1<<", "<<firstspill2<<", "<<firstSpillDstof<<endl;

  double recSpill1 = firstspill1;
  double recSpill2 = firstspill2;
  
  for(int i=0; i<spills1.size(); i++) {
    hfile1->Fill(spills1[i]);// - firstspill1);
    hdiff1->Fill(spills1[i] - recSpill1);
    recSpill1 = spills1[i];
  }
  for(int i=0; i<spills2.size(); i++) {
    hfile2->Fill(spills2[i]); //- firstspill2);
    hdiff2->Fill(spills2[i] - recSpill2);
    recSpill2 = spills2[i];
  }
  double recSpillDstof = firstSpillDstof;
  for(int i=0; i<spillsDstof.size(); i++) {
    hdstof->Fill(spillsDstof[i]);
    hdiffdstof->Fill(spillsDstof[i] - recSpillDstof);
    recSpillDstof = spillsDstof[i];
  }
  
  hfile1->SetLineColor(kRed);
  hfile1->SetFillColor(kRed);
  hfile1->SetFillStyle(3005);
  hfile2->SetLineColor(kBlue);
  hfile2->SetFillColor(kBlue);
  hfile2->SetFillStyle(3004);
  hdstof->SetLineColor(kBlack);
  hdstof->SetFillColor(kBlack);
  hdstof->SetFillStyle(3006);
  TCanvas *c1 = new TCanvas("c1");
  hfile1->SetTitle("Time of beam signals in various files/systems; Time / s");
  hfile1->Draw("hist");
  hfile2->Draw("hist same");
  hdstof->Draw("hist same");
  TLegend *leg1 = new TLegend(0, 0.83, 0.2, 1);
  leg1->AddEntry(hfile1, "Data_2018_9_12_b1.root", "f");
  leg1->AddEntry(hfile2, "Data_2018_9_13_b1.root", "f");
  leg1->AddEntry(hdstof, "Dstof beam database", "f");
  leg1->Draw();
  c1->Print("spillTimeComp.png");
  c1->Print("spillTimeComp.pdf");
  
  hdiff1->SetLineColor(kRed);
  hdiff1->SetFillColor(kRed);
  hdiff1->SetFillStyle(3005);
  hdiff2->SetLineColor(kBlue);
  hdiff2->SetFillColor(kBlue);
  hdiff2->SetFillStyle(3004);
  hdiffdstof->SetLineColor(kBlack);
  hdiffdstof->SetFillColor(kBlack);
  hdiffdstof->SetFillStyle(3006);
  TCanvas *c2 = new TCanvas("c2");
  hdiff1->Draw("hist");
  hdiff2->Draw("hist same");
  hdiffdstof->Draw("hist same");
  leg1->Draw();
  c1->Print("spillTimeDiff.png");
  c1->Print("spillTimeDiff.pdf");
}


