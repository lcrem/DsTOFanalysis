// overRun
void overRun(const char* saveDir,
	     const char* deadTimeDir="/zfs_home/sjones/plots/deadtime/outFiles/", 
	     const char* protonDir="/zfs_home/sjones/plots/intProtons/")
{
  gROOT->SetBatch(kTRUE);

  // Unix timestamp for start of electron target data
  const double eTargetT = 1536145800.;

  TFile *fout = new TFile(Form("%s/overRun_out.root", saveDir), "recreate");
  
  // Go through deadtime directory and get the files
  const char* ext  = ".root";
  const char* pref = "Data";
  TString str;
  const char *entry;
  char *dir  = gSystem->ExpandPathName(deadTimeDir);
  void *dirp = gSystem->OpenDirectory(dir);

  //  TGraphErrors *dtGr = new TGraphErrors();

  vector<double> rmsVec;
  vector<double> interceptVec;
  vector<double> slopeVec;
  vector<double> dtVec;
  vector<double> timeVec;
  vector<double> zeroVec;

  while (entry = (char*)gSystem->GetDirEntry(dirp)) {
    str = entry;
    if (str.EndsWith(ext) && str.Contains(pref) && !str.Contains("3block") && !str.Contains("2block") && !str.Contains("1block") && !str.Contains("0block")) {
      cout<<str<<endl;
      str.Prepend(deadTimeDir);
      TFile *dtFile = new TFile(str, "read");
      TNamed *start = 0;
      TNamed *end   = 0;
      dtFile->GetObject("start_of_run", start);
      dtFile->GetObject("end_of_run", end);
      const char* startchar = start->GetTitle();
      const char* endchar   = end->GetTitle();
      string startstr(startchar);
      string endstr(endchar);
      string unixstart = startstr.substr(25,10);
      string unixend   = endstr.substr(23,10);
      double dtStart = stod(unixstart);	 
      double dtEnd   = stod(unixend);
      cout<<dtStart<<" "<<dtEnd<<endl;
      TObject* nS1S2dtof = 0;
      TObject* nS1S2utof = 0;
      dtFile->GetObject("nS1S2dtof", nS1S2dtof);
      dtFile->GetObject("nS1S2utof", nS1S2utof);

      double ratio = (double)nS1S2utof->GetUniqueID()/(double)nS1S2dtof->GetUniqueID();
      cout<<nS1S2dtof->GetUniqueID()<<" "<<nS1S2utof->GetUniqueID()<<" "<<ratio<<endl;

      TVectorD *rms = (TVectorD*)dtFile->Get("rms");
      TVectorD *intercept = (TVectorD*)dtFile->Get("intercept");
      TVectorD *slope = (TVectorD*)dtFile->Get("slope");
      rmsVec.push_back((*rms)[0]);
      interceptVec.push_back((*intercept)[0]);
      slopeVec.push_back((*slope)[0]);
      dtVec.push_back(ratio);
      zeroVec.push_back(0.);
      timeVec.push_back((dtStart+dtEnd)/2.);

      //dtGr->SetPoint(dtGr->GetN(), (dtStart+dtEnd)/2., ratio);
      cout<<(*rms)[0]<<" "<<(*intercept)[0]<<endl;

      delete nS1S2dtof;
      delete nS1S2utof;
      dtFile->Close();
      delete dtFile;
    } // Is correct file
  } // while (entry = (char*)gSystem->GetDirEntry(dirp))

  TGraphErrors *dtGr = new TGraphErrors(rmsVec.size(), &timeVec[0], &dtVec[0], &zeroVec[0], &rmsVec[0]);
  TGraph *interceptGr = new TGraph(rmsVec.size(), &timeVec[0], &interceptVec[0]);
  TGraph *slopeGr     = new TGraph(rmsVec.size(), &timeVec[0], &slopeVec[0]);

  fout->cd();
  dtGr->SetTitle("S1 #cap S2 (utof / dtof) for each utof run; Time; (S1 #cap S2 utof)/(S1 #cap S2 dtof)");
  //  dtGr->GetXaxis()->SetRangeUser(1534.8e6, 1536.2e6);
  dtGr->Write("dtGr");
  // Open protons per spill file
  TFile *filePro = new TFile(Form("%s/intProtons.root", protonDir));
  TGraph *proGr = new TGraph();
  TGraph *piGr  = new TGraph();
  TGraph *ratioGr = new TGraph();
  // Go through all the non-electron target data and every 100 spills 
  // print the average protons per spill to a TGraph
  TTree *outTree = (TTree*)filePro->Get("outTree");
  int nP;
  int nPi;
  double spillTime;
  TString file;
  outTree->SetBranchAddress("nP", &nP);
  outTree->SetBranchAddress("nPi", &nPi);
  outTree->SetBranchAddress("spillTime", &spillTime);
  outTree->SetBranchAddress("file", &file);

  int sumP  = 0;
  int sumPi = 0;
  int nSpills = 0;
  double sumTime = 0.;
  for (int t = 0; t < outTree->GetEntries(); t++) {
    outTree->GetEntry(t);
    if (file != "Data_2018_8_30_b2" && nPi < 5000 && spillTime < eTargetT && spillTime >1534.8e6 && nP < 150) {
      sumP  += nP;
      sumPi += nPi;
      sumTime += spillTime;
      nSpills++;
      if (nSpills == 50) {
	double proPi = (double)sumP/(double)sumPi;
	proGr->SetPoint(proGr->GetN(), sumTime/(double)nSpills, (double)sumP/(double)nSpills);
	piGr->SetPoint(piGr->GetN(), sumTime/(double)nSpills, (double)sumPi/(double)nSpills);
	ratioGr->SetPoint(ratioGr->GetN(), sumTime/(double)nSpills, proPi);
	sumP  = 0;
	sumPi = 0;
	nSpills = 0;
	sumTime = 0;
      } // if (nSpills == 100) 
    } // if (file != "Data_2018_8_30_b2")
  } // for (int t = 0; t < outTree->GetEntries(); t++)
  proGr->SetTitle("50 spill average of number of protons per spill for 4 block data; Time; Protons / spill"); 
  piGr->SetTitle("50 spill average of number of MIPs per spill for 4 block data; Time; MIPs / spill"); 
  ratioGr->SetTitle("50 spill average of number of proton/MIP ratio for 4 block data; Time; Protons / MIPs"); 
  // Now get the S1S2 hits for the same time period
  TTree *s1s2OutTree = (TTree*)filePro->Get("s1s2OutTree");
  int nS1S2;
  double dtofSpillTime;
  s1s2OutTree->SetBranchAddress("nS1S2", &nS1S2);
  s1s2OutTree->SetBranchAddress("dtofSpillTime", &dtofSpillTime);
  int sumS1S2 = 0;
  nSpills = 0;
  TGraph *piS1S2Gr = new TGraph();
  for (int i=0; i<s1s2OutTree->GetEntries(); i++) {
    s1s2OutTree->GetEntry(i);
    sumS1S2 += nS1S2;
    sumTime += dtofSpillTime;
    nSpills++;
    if (nSpills % 50 == 0) {
      if ((sumTime/50.) < 1537e6) {
	piS1S2Gr->SetPoint(piS1S2Gr->GetN(), sumTime/50., (double)sumS1S2/50.);
      }
      sumTime = 0.;
      sumS1S2 = 0;
    } // if (nSpills == 50) 
  } // for (int i=0; i<s1s2OutTree->GetEntries(); i++)

  fout->cd();
  proGr->GetXaxis()->SetRangeUser(1534.8e6, 1536.2e6);
  piGr->GetXaxis()->SetRangeUser(1534.8e6, 1536.2e6);
  piS1S2Gr->GetXaxis()->SetRangeUser(1534.8e6, 1536.2e6);
  ratioGr->GetXaxis()->SetRangeUser(1534.8e6, 1536.2e6);
  proGr->Write("proGr");
  piGr->Write("piGr");
  ratioGr->Write("ratioGr");
  slopeGr->Write("slopeGr");
  interceptGr->Write("interceptGr");
  piS1S2Gr->Write("piS1S2Gr");

  fout->Close();
  filePro->Close();
} // overRun
