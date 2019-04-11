// overRun
void overRun(const char* saveDir,
	     const char* deadTimeDir="/zfs_home/sjones/plots/deadtime/outFiles/", 
	     const char* protonDir="/zfs_home/sjones/plots/intProtons/")
{
  gROOT->SetBatch(kTRUE);

  TFile *filePro = new TFile(Form("%s/intProtons.root", protonDir));
  
  // Go through deadtime directory and get the files
  const char* ext  = ".root";
  const char* pref = "Data";
  TString str;
  const char *entry;
  char *dir  = gSystem->ExpandPathName(deadTimeDir);
  void *dirp = gSystem->OpenDirectory(dir);

  TGraph *dtGraph = new TGraph();

  while (entry = (char*)gSystem->GetDirEntry(dirp)) {
    str = entry;
    cout<<str<<endl;
    if (str.EndsWith(ext) && str.Contains(pref)) {
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

      cout<<nS1S2dtof<<" "<<nS1S2utof<<endl;

      dtFile->Close();
    }
  } // while (entry = (char*)gSystem->GetDirEntry(dirp))

  filePro->Close();
} // overRun
