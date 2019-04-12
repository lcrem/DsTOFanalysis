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

  TGraph *dtGr = new TGraph();

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

      TObject* nS1S2dtof = 0;
      TObject* nS1S2utof = 0;
      dtFile->GetObject("nS1S2dtof", nS1S2dtof);
      dtFile->GetObject("nS1S2utof", nS1S2utof);

      double ratio = (double)nS1S2utof->GetUniqueID()/(double)nS1S2dtof->GetUniqueID();
      cout<<nS1S2dtof->GetUniqueID()<<" "<<nS1S2utof->GetUniqueID()<<" "<<ratio<<endl;

      dtGr->SetPoint(dtGr->GetN(), (dtStart+dtEnd)/2., ratio);

      delete nS1S2dtof;
      delete nS1S2utof;
      dtFile->Close();
    } // Is correct file
  } // while (entry = (char*)gSystem->GetDirEntry(dirp))

  fout->cd();
  dtGr->SetTitle("S1 #cap S2 (utof / dtof) for each utof run; Time; (S1 #cap S2 utof)/(S1 #cap S2 dtof)");
  dtGr->GetXaxis()->SetRangeUser(1534.8e6, 1536.2e6);
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

  fout->cd();
  proGr->GetXaxis()->SetRangeUser(1534.8e6, 1536.2e6);
  piGr->GetXaxis()->SetRangeUser(1534.8e6, 1536.2e6);
  ratioGr->GetXaxis()->SetRangeUser(1534.8e6, 1536.2e6);
  proGr->Write("proGr");
  piGr->Write("piGr");
  ratioGr->Write("ratioGr");

  fout->Close();
  filePro->Close();
} // overRun
