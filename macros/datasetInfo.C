// block4Times.C
// Beam data about the various runs used in the analysis

// Gets timestamp from input files
TTimeStamp getStampFromLine(string line) {
  stringstream ssyear(line.substr(0,4));
  stringstream ssmonth(line.substr(5,2));
  stringstream ssday(line.substr(8,2));
  stringstream sshour(line.substr(11,2));
  stringstream ssmin(line.substr(14,2));
  stringstream sssec(line.substr(17,2));
  int year = 0;
  int month = 0;
  int day = 0;
  int hour = 0;
  int min = 0;
  int sec = 0;
  int count = 0;
  ssyear >> year;
  ssmonth >> month;
  ssday >> day;
  sshour >> hour;
  ssmin >> min;
  sssec >> sec;
  // Times are in CEST so require a 2 hour offset to get to UTC
  TTimeStamp stamp(year, month, day, hour, min, sec, 0, "kTRUE", -7200);
  return stamp;
}

void datasetInfo(const char* saveDir,
		 const char* ustofDir="/zfs_home/sjones/mylinktoutof/",
		 const char* dstofDir="/scratch0/dbrailsf/temp/mylinktodtof/",
		 const char* spillDir="/scratch2/sjones/spillDB/",
		 const char* scintFile="/scratch2/sjones/DsTOFanalysis/files/Scintillator_info.txt",
		 const char* bzh01File="/scratch2/sjones/DsTOFanalysis/files/BZH01Values.txt",
		 const char* bzh03File="/scratch2/sjones/DsTOFanalysis/files/BZH03Values.txt") 
{
  const int cutTime = 1535.8698e6;
  TFile *fout = new TFile(Form("%s/datasetInfo_out.root", saveDir), "recreate");
  gROOT->SetBatch(kTRUE);
  // Get the scintillator values from the text file and put into a vector
  TGraph *grScintAll = new TGraph();
  std::vector<int> countVec;
  std::vector<TTimeStamp> stampVec;
  string line;
  ifstream inFile;
  inFile.open(scintFile);
  while(!inFile.eof()) {
    getline(inFile, line);
    if (line.length()>20) {
      stringstream sscount(line.substr(24));
      int count = 0;
      sscount >> count;
      TTimeStamp stamp = getStampFromLine(line);
      countVec.push_back(count);
      stampVec.push_back(stamp);
      grScintAll->SetPoint(grScintAll->GetN(), stamp, count);
    } // if (line.length()>20)
  } // while(!inFile.eof())
  inFile.close();
  cout<<"Scintillator file has "<<countVec.size()<<" spills recorded"<<endl;

  // Do same thing for BZH01 (beam momentum magnet) and BZH03 (bending magnet) values
  TGraph *grbzh03 = new TGraph();
  std::vector<int> bzh03Vec;
  std::vector<TTimeStamp> bzh03StampVec;
  ifstream bzh03Stream;
  bzh03Stream.open(bzh03File);
  int nLines = 0;
  while(!bzh03Stream.eof()) {
    getline(bzh03Stream, line);
    if (nLines!=0 && line.length()>20) {
      stringstream ssbzh03(line.substr(24));
      int bzh03 = 0;
      ssbzh03 >> bzh03;
      TTimeStamp stamp = getStampFromLine(line);
      bzh03Vec.push_back(bzh03);
      bzh03StampVec.push_back(stamp);
      if (bzh03 != 0.) {
	grbzh03->SetPoint(grbzh03->GetN(), stamp, bzh03);
      } 
    } 
    nLines++;
  }
  bzh03Stream.close();
  // BZH01
  TGraph *grbzh01 = new TGraph();
  std::vector<int> bzh01Vec;
  std::vector<TTimeStamp> bzh01StampVec;
  ifstream bzh01Stream;
  bzh01Stream.open(bzh01File);
  nLines = 0;
  while(!bzh01Stream.eof()) {
    getline(bzh01Stream, line);
    if (nLines!=0 && line.length()>20) {
      stringstream ssbzh01(line.substr(24));
      double bzh01 = 0;
      ssbzh01 >> bzh01;
      TTimeStamp stamp = getStampFromLine(line);
      bzh01Vec.push_back(bzh01);
      bzh01StampVec.push_back(stamp);
      //      cout<<"bzh01 "<<stamp<<endl;
      if (bzh01 != 0.) {
	grbzh01->SetPoint(grbzh01->GetN(), stamp, bzh01);
      }
    } 
    nLines++;
  }
  bzh01Stream.close();
  cout<<"BZH01 file has "<<bzh01Vec.size()<<" spills recorded"<<endl;

  fout->cd();
  grbzh01->Write("grbzh01All");
  grbzh03->Write("grbzh03All");

  // For each set of data get the start and end times and put in an individual graph
  const char* str4Block1 = "Data_2018_9_1_b8_800MeV_4block_bend4cm.root";
  const char* str4Block2 = "Data_2018_9_3_b2_800MeV_4block_bend4cm.root";
  const char* str4Block3 = "Data_2018_9_1_b5_800MeV_4block.root";
  const char* str4Block4 = "Data_2018_9_4_b1_800MeV_4block.root";
  const char* str4Block5 = "Data_2018_9_4_b2_800MeV_4block.root";
  const char* str4Block6 = "Data_2018_9_3_b4_800MeV_4block_1mUSB.root";
  const char* str4Block7 = "Data_2018_8_30_b1.root";
  const char* str4Block8 = "Data_2018_8_29_b4.root";
  //  const char* str4Block9 = "Data_2018_8_29_b2.root";
  const char* str4Block10 = "Data_2018_8_29_b1.root";
  const char* str4Block11 = "Data_2018_8_28_b6.root";
  const char* str4Block12 = "Data_2018_8_28_b5.root";
  //  const char* str4Block13 = "Data_2018_8_27_bc2_800MeV_4block.root";
  //  const char* str4Block14 = "Data_2018_8_27_b1_800MeV_4block.root";
  const char* str4Block15 = "Data_2018_8_26_b3_800MeV_4block.root";
  const char* str4Block16 = "Data_2018_8_26_b1_800MeV_4block.root";
  const char* str4Block17 = "Data_2018_8_25_b3_800MeV_4block.root";
  const char* str4Block18 = "Data_2018_8_25_b2_800MeV_4block.root";
  const char* str0Block = "Data_2018_8_31_b2_800MeV_0block.root";
  const char* str1Block = "Data_2018_9_1_b4_800MeV_1block_bend4cm.root";
  const char* str2Block = "Data_2018_9_1_b2_800MeV_2block_bend4cm.root";
  const char* str3Block = "Data_2018_9_1_b3_800MeV_3block_bend4cm.root";
  std::vector<const char*> blockVec={str4Block1, str4Block2, str4Block3, str4Block4, 
				     str4Block5, str4Block6, str4Block7, str4Block8,
				     /*str4Block9,*/ 
				     str4Block10, str4Block11, str4Block12,
				     /*str4Block13, str4Block14,*/
				     str4Block15, str4Block16,
				     str4Block17, str4Block18,
				     str0Block, str1Block, str2Block, str3Block};
  TMultiGraph *mg = new TMultiGraph();
  TMultiGraph *mgbzh01 = new TMultiGraph();
  TMultiGraph *mgbzh03 = new TMultiGraph();
  TMultiGraph *mgbzh01Avg = new TMultiGraph();
  TMultiGraph *mgbzh03Avg = new TMultiGraph();
  TMultiGraph *mgMagnet = new TMultiGraph();
  TMultiGraph *mgMagnetAvg = new TMultiGraph();
  TMultiGraph *mgScint = new TMultiGraph();
  TMultiGraph *mgScintAvg = new TMultiGraph();
  TGraph *grDtof = new TGraph();
  for (int b=0; b<blockVec.size(); b++) {
    TFile *futof = new TFile(Form("%s/%s",ustofDir,blockVec[b]), "read");
    TTree *tree = (TTree*)futof->Get("tree");
    double tS1;
    tree->SetBranchAddress("tS1", &tS1);
    TNamed *start = 0;
    TNamed *end   = 0;
    futof->GetObject("start_of_run", start);
    futof->GetObject("end_of_run", end);
    const char* startchar = start->GetTitle();
    std::string startstr(startchar);
    std::string unixstart = startstr.substr(25,10);
    int startTimeUtof = stoi(unixstart);
    tree->GetEntry(tree->GetEntries()-1);
    int endTimeUtof = (tS1/1e9) + startTimeUtof;
    TTimeStamp utofStartStamp(startTimeUtof);
    TTimeStamp utofEndStamp(endTimeUtof);

    TGraph *grUtof = new TGraph();
    TGraph *grScint = new TGraph();
    TGraph *grbzh01Run = new TGraph();
    TGraph *grbzh03Run = new TGraph();
    TGraph *grbzh01RunAvg = new TGraph();
    TGraph *grbzh03RunAvg = new TGraph();
    TGraph *grScintAvg = new TGraph();
    if (b < 15) {
      grUtof->SetPoint(grUtof->GetN(), startTimeUtof, 3500);
      grUtof->SetPoint(grUtof->GetN(), endTimeUtof, 3500);
    }
    else {
      grUtof->SetPoint(grUtof->GetN(), startTimeUtof, 3000);
      grUtof->SetPoint(grUtof->GetN(), endTimeUtof, 3000);
    }

    if (b==0) {
      grUtof->SetMarkerColor(kBlack);
      grUtof->SetMarkerStyle(21);
      grScint->SetMarkerColor(kBlack);
      grScintAvg->SetMarkerColor(kBlack);
      grScintAvg->SetMarkerStyle(3);
      grbzh01Run->SetMarkerColor(kBlack);
      grbzh01Run->SetMarkerStyle(4);
      grbzh03Run->SetMarkerColor(kBlack);
      grbzh03Run->SetMarkerStyle(5);
      grbzh01RunAvg->SetMarkerColor(kBlack);
      grbzh01RunAvg->SetMarkerStyle(4);
      grbzh03RunAvg->SetMarkerColor(kBlack);
      grbzh03RunAvg->SetMarkerStyle(5);
    }
    else if (b==1) {
      grUtof->SetMarkerColor(kRed);
      grUtof->SetMarkerStyle(47);
      grScint->SetMarkerColor(kRed);
      grScintAvg->SetMarkerColor(kRed);
      grScintAvg->SetMarkerStyle(3);
      grbzh01Run->SetMarkerColor(kRed);
      grbzh01Run->SetMarkerStyle(4);
      grbzh03Run->SetMarkerColor(kRed);
      grbzh03Run->SetMarkerStyle(5);
      grbzh01RunAvg->SetMarkerColor(kRed);
      grbzh01RunAvg->SetMarkerStyle(4);
      grbzh03RunAvg->SetMarkerColor(kRed);
      grbzh03RunAvg->SetMarkerStyle(5);
    }
    else if (b==2) {
      grUtof->SetMarkerColor(kBlue);
      grUtof->SetMarkerStyle(20);
      grScint->SetMarkerColor(kBlue);
      grScintAvg->SetMarkerColor(kBlue);
      grScintAvg->SetMarkerStyle(3);
      grbzh01Run->SetMarkerColor(kBlue);
      grbzh01Run->SetMarkerStyle(4);
      grbzh03Run->SetMarkerColor(kBlue);
      grbzh03Run->SetMarkerStyle(5);
      grbzh01RunAvg->SetMarkerColor(kBlue);
      grbzh01RunAvg->SetMarkerStyle(4);
      grbzh03RunAvg->SetMarkerColor(kBlue);
      grbzh03RunAvg->SetMarkerStyle(5);
    }
    else if (b==3) {
      grUtof->SetMarkerColor(kCyan+1);
      grUtof->SetMarkerStyle(22);
      grScint->SetMarkerColor(kCyan+1);
      grScintAvg->SetMarkerColor(kCyan+1);
      grScintAvg->SetMarkerStyle(3);
      grbzh01Run->SetMarkerColor(kCyan+1);
      grbzh01Run->SetMarkerStyle(4);
      grbzh03Run->SetMarkerColor(kCyan+1);
      grbzh03Run->SetMarkerStyle(5);
      grbzh01RunAvg->SetMarkerColor(kCyan+1);
      grbzh01RunAvg->SetMarkerStyle(4);
      grbzh03RunAvg->SetMarkerColor(kCyan+1);
      grbzh03RunAvg->SetMarkerStyle(5);
    }
    else if (b==4) {
      grUtof->SetMarkerColor(kGreen+2);
      grUtof->SetMarkerStyle(34);
      grScint->SetMarkerColor(kGreen+2);
      grScintAvg->SetMarkerColor(kGreen+2);
      grScintAvg->SetMarkerStyle(3);
      grbzh01Run->SetMarkerColor(kGreen+2);
      grbzh01Run->SetMarkerStyle(4);
      grbzh03Run->SetMarkerColor(kGreen+2);
      grbzh03Run->SetMarkerStyle(5);
      grbzh01RunAvg->SetMarkerColor(kGreen+2);
      grbzh01RunAvg->SetMarkerStyle(4);
      grbzh03RunAvg->SetMarkerColor(kGreen+2);
      grbzh03RunAvg->SetMarkerStyle(5);
    }
    else if (b==5) {
      grUtof->SetMarkerColor(kOrange+1);
      grUtof->SetMarkerStyle(49);
      grScint->SetMarkerColor(kOrange+1);
      grScintAvg->SetMarkerColor(kOrange+1);
      grScintAvg->SetMarkerStyle(3);
      grbzh01Run->SetMarkerColor(kOrange+1);
      grbzh01Run->SetMarkerStyle(4);
      grbzh03Run->SetMarkerColor(kOrange+1);
      grbzh03Run->SetMarkerStyle(5);
      grbzh01RunAvg->SetMarkerColor(kOrange+1);
      grbzh01RunAvg->SetMarkerStyle(4);
      grbzh03RunAvg->SetMarkerColor(kOrange+1);
      grbzh03RunAvg->SetMarkerStyle(5);
    }
    else if (b==6) {
      grUtof->SetMarkerColor(kSpring);
      grUtof->SetMarkerStyle(23);
      grScint->SetMarkerColor(kSpring);
      grScintAvg->SetMarkerColor(kSpring);
      grScintAvg->SetMarkerStyle(3);
      grbzh01Run->SetMarkerColor(kSpring);
      grbzh01Run->SetMarkerStyle(4);
      grbzh03Run->SetMarkerColor(kSpring);
      grbzh03Run->SetMarkerStyle(5);
      grbzh01RunAvg->SetMarkerColor(kSpring);
      grbzh01RunAvg->SetMarkerStyle(4);
      grbzh03RunAvg->SetMarkerColor(kSpring);
      grbzh03RunAvg->SetMarkerStyle(5);
    }
    else if (b==7) {
      grUtof->SetMarkerColor(kBlack);
      grUtof->SetMarkerStyle(21);
      grScint->SetMarkerColor(kBlack);
      grScintAvg->SetMarkerColor(kBlack);
      grScintAvg->SetMarkerStyle(3);
      grbzh01Run->SetMarkerColor(kBlack);
      grbzh01Run->SetMarkerStyle(4);
      grbzh03Run->SetMarkerColor(kBlack);
      grbzh03Run->SetMarkerStyle(5);
      grbzh01RunAvg->SetMarkerColor(kBlack);
      grbzh01RunAvg->SetMarkerStyle(4);
      grbzh03RunAvg->SetMarkerColor(kBlack);
      grbzh03RunAvg->SetMarkerStyle(5);
    }
    else if (b==8) {
      grUtof->SetMarkerColor(kRed);
      grUtof->SetMarkerStyle(47);
      grScint->SetMarkerColor(kRed);
      grScintAvg->SetMarkerColor(kRed);
      grScintAvg->SetMarkerStyle(3);
      grbzh01Run->SetMarkerColor(kRed);
      grbzh01Run->SetMarkerStyle(4);
      grbzh03Run->SetMarkerColor(kRed);
      grbzh03Run->SetMarkerStyle(5);
      grbzh01RunAvg->SetMarkerColor(kRed);
      grbzh01RunAvg->SetMarkerStyle(4);
      grbzh03RunAvg->SetMarkerColor(kRed);
      grbzh03RunAvg->SetMarkerStyle(5);
    }
    else if (b==9) {
      grUtof->SetMarkerColor(kBlue);
      grUtof->SetMarkerStyle(20);
      grScint->SetMarkerColor(kBlue);
      grScintAvg->SetMarkerColor(kBlue);
      grScintAvg->SetMarkerStyle(3);
      grbzh01Run->SetMarkerColor(kBlue);
      grbzh01Run->SetMarkerStyle(4);
      grbzh03Run->SetMarkerColor(kBlue);
      grbzh03Run->SetMarkerStyle(5);
      grbzh01RunAvg->SetMarkerColor(kBlue);
      grbzh01RunAvg->SetMarkerStyle(4);
      grbzh03RunAvg->SetMarkerColor(kBlue);
      grbzh03RunAvg->SetMarkerStyle(5);
    }
    else if (b==10) {
      grUtof->SetMarkerColor(kCyan+1);
      grUtof->SetMarkerStyle(22);
      grScint->SetMarkerColor(kCyan+1);
      grScintAvg->SetMarkerColor(kCyan+1);
      grScintAvg->SetMarkerStyle(3);
      grbzh01Run->SetMarkerColor(kCyan+1);
      grbzh01Run->SetMarkerStyle(4);
      grbzh03Run->SetMarkerColor(kCyan+1);
      grbzh03Run->SetMarkerStyle(5);
      grbzh01RunAvg->SetMarkerColor(kCyan+1);
      grbzh01RunAvg->SetMarkerStyle(4);
      grbzh03RunAvg->SetMarkerColor(kCyan+1);
      grbzh03RunAvg->SetMarkerStyle(5);
    }
    else if (b==11) {
      grUtof->SetMarkerColor(kGreen+2);
      grUtof->SetMarkerStyle(34);
      grScint->SetMarkerColor(kGreen+2);
      grScintAvg->SetMarkerColor(kGreen+2);
      grScintAvg->SetMarkerStyle(3);
      grbzh01Run->SetMarkerColor(kGreen+2);
      grbzh01Run->SetMarkerStyle(4);
      grbzh03Run->SetMarkerColor(kGreen+2);
      grbzh03Run->SetMarkerStyle(5);
      grbzh01RunAvg->SetMarkerColor(kGreen+2);
      grbzh01RunAvg->SetMarkerStyle(4);
      grbzh03RunAvg->SetMarkerColor(kGreen+2);
      grbzh03RunAvg->SetMarkerStyle(5);
    }
    else if (b==12) {
      grUtof->SetMarkerColor(kOrange+1);
      grUtof->SetMarkerStyle(49);
      grScint->SetMarkerColor(kOrange+1);
      grScintAvg->SetMarkerColor(kOrange+1);
      grScintAvg->SetMarkerStyle(3);
      grbzh01Run->SetMarkerColor(kOrange+1);
      grbzh01Run->SetMarkerStyle(4);
      grbzh03Run->SetMarkerColor(kOrange+1);
      grbzh03Run->SetMarkerStyle(5);
      grbzh01RunAvg->SetMarkerColor(kOrange+1);
      grbzh01RunAvg->SetMarkerStyle(4);
      grbzh03RunAvg->SetMarkerColor(kOrange+1);
      grbzh03RunAvg->SetMarkerStyle(5);
    }
    else if (b==13) {
      grUtof->SetMarkerColor(kSpring);
      grUtof->SetMarkerStyle(23);
      grScint->SetMarkerColor(kSpring);
      grScintAvg->SetMarkerColor(kSpring);
      grScintAvg->SetMarkerStyle(3);
      grbzh01Run->SetMarkerColor(kSpring);
      grbzh01Run->SetMarkerStyle(4);
      grbzh03Run->SetMarkerColor(kSpring);
      grbzh03Run->SetMarkerStyle(5);
      grbzh01RunAvg->SetMarkerColor(kSpring);
      grbzh01RunAvg->SetMarkerStyle(4);
      grbzh03RunAvg->SetMarkerColor(kSpring);
      grbzh03RunAvg->SetMarkerStyle(5);
    }
    else if (b==14) {
      grUtof->SetMarkerColor(kBlack);
      grUtof->SetMarkerStyle(21);
      grScint->SetMarkerColor(kBlack);
      grScintAvg->SetMarkerColor(kBlack);
      grScintAvg->SetMarkerStyle(3);
      grbzh01Run->SetMarkerColor(kBlack);
      grbzh01Run->SetMarkerStyle(4);
      grbzh03Run->SetMarkerColor(kBlack);
      grbzh03Run->SetMarkerStyle(5);
      grbzh01RunAvg->SetMarkerColor(kBlack);
      grbzh01RunAvg->SetMarkerStyle(4);
      grbzh03RunAvg->SetMarkerColor(kBlack);
      grbzh03RunAvg->SetMarkerStyle(5);
    }
    /*
    else if (b==15) {
      grUtof->SetMarkerColor(kRed);
      grUtof->SetMarkerStyle(47);
      grScint->SetMarkerColor(kRed);
      grScintAvg->SetMarkerColor(kRed);
      grScintAvg->SetMarkerStyle(3);
      grbzh01Run->SetMarkerColor(kRed);
      grbzh01Run->SetMarkerStyle(4);
      grbzh03Run->SetMarkerColor(kRed);
      grbzh03Run->SetMarkerStyle(5);
      grbzh01RunAvg->SetMarkerColor(kRed);
      grbzh01RunAvg->SetMarkerStyle(4);
      grbzh03RunAvg->SetMarkerColor(kRed);
      grbzh03RunAvg->SetMarkerStyle(5);
    }
    else if (b==16) {
      grUtof->SetMarkerColor(kBlue);
      grUtof->SetMarkerStyle(20);
      grScint->SetMarkerColor(kBlue);
      grScintAvg->SetMarkerColor(kBlue);
      grScintAvg->SetMarkerStyle(3);
      grbzh01Run->SetMarkerColor(kBlue);
      grbzh01Run->SetMarkerStyle(4);
      grbzh03Run->SetMarkerColor(kBlue);
      grbzh03Run->SetMarkerStyle(5);
      grbzh01RunAvg->SetMarkerColor(kBlue);
      grbzh01RunAvg->SetMarkerStyle(4);
      grbzh03RunAvg->SetMarkerColor(kBlue);
      grbzh03RunAvg->SetMarkerStyle(5);
    }
    else if (b==17) {
      grUtof->SetMarkerColor(kCyan+1);
      grUtof->SetMarkerStyle(22);
      grScint->SetMarkerColor(kCyan+1);
      grScintAvg->SetMarkerColor(kCyan+1);
      grScintAvg->SetMarkerStyle(3);
      grbzh01Run->SetMarkerColor(kCyan+1);
      grbzh01Run->SetMarkerStyle(4);
      grbzh03Run->SetMarkerColor(kCyan+1);
      grbzh03Run->SetMarkerStyle(5);
      grbzh01RunAvg->SetMarkerColor(kCyan+1);
      grbzh01RunAvg->SetMarkerStyle(4);
      grbzh03RunAvg->SetMarkerColor(kCyan+1);
      grbzh03RunAvg->SetMarkerStyle(5);
    }
    */
    else {
      grUtof->SetMarkerColor(kMagenta+1);
      grUtof->SetMarkerStyle(39);
      grScint->SetMarkerColor(kMagenta+1);
      grScintAvg->SetMarkerColor(kMagenta+1);
      grScintAvg->SetMarkerStyle(3);
      grbzh01Run->SetMarkerColor(kMagenta+1);
      grbzh01Run->SetMarkerStyle(4);
      grbzh03Run->SetMarkerColor(kMagenta+1);
      grbzh03Run->SetMarkerStyle(5);
      grbzh01RunAvg->SetMarkerColor(kMagenta+1);
      grbzh01RunAvg->SetMarkerStyle(4);
      grbzh03RunAvg->SetMarkerColor(kMagenta+1);
      grbzh03RunAvg->SetMarkerStyle(5);
    }
    
    grUtof->SetMarkerSize(2);
    mg->Add(grUtof);
    futof->Close();
    delete futof;
    fout->cd();
    grUtof->Write(Form("grUtof_%s", blockVec[b]));

    int nScint = 1;
    int nCount = 1;
    int countTotal = 0;
    int spillTotal = 0;
    vector<double> countVecT;
    vector<double> scintVecT;
    int nCountT = 0;
    double lastScintT = 300.;
    // Find appropriate scintillator values and put into graph
    // Also want to make plot of number of scintillator hits vs number of 
    // S1S2 coincidences over some timescale (say 5 minutes)
    for (int s=0; s<stampVec.size(); ++s) {
      if (stampVec[s] > utofEndStamp) break;
      if (stampVec[s] < utofStartStamp) continue;

      if (stampVec[s] > utofStartStamp && stampVec[s] < utofEndStamp) {
	grScint->SetPoint(grScint->GetN(), stampVec[s], countVec[s]);
	countTotal += countVec[s];
	spillTotal++;
	nScint++;
	nCount+=countVec[s];
	nCountT += countVec[s];
	if (stampVec[s].GetSec() > utofStartStamp.GetSec() + lastScintT) {
	  scintVecT.push_back(utofStartStamp + lastScintT);
	  lastScintT += 300.;
	  countVecT.push_back(nCountT);
	  nCountT = 0;
	}

	if (nScint % 20 == 0) {
	  grScintAvg->SetPoint(grScintAvg->GetN(), stampVec[s], (double)nCount/20.); 
	  nCount = 0;
	} // if (nScint % 20 == 0)
      } // if (stampVec[s] > utofStartStamp && stampVec[s] < utofEndStamp)
    } // for (int s=0; s<stampVec.size(); ++s)
    mg->Add(grScint);
    mgScint->Add(grScint);
    mgScintAvg->Add(grScintAvg);
    fout->cd();
    grScint->Write(Form("grScint_%s", blockVec[b]));
    grScintAvg->Write(Form("grScintAvg_%s", blockVec[b]));
    // Print out average number of scintillator counts per spill
    cout<<blockVec[b]<<" "<<(double)countTotal/(double)spillTotal<<endl;
    // Same thing for BZH01
    int bzh01Count = 1;
    double bzh01Sum = 0;
    for (int s=0; s<bzh01Vec.size(); ++s) {
      if (bzh01StampVec[s] > utofEndStamp) break;
      if (bzh01StampVec[s] < utofStartStamp) continue;

      if (bzh01StampVec[s] > utofStartStamp && bzh01StampVec[s] < utofEndStamp) {
	grbzh01Run->SetPoint(grbzh01Run->GetN(), bzh01StampVec[s], bzh01Vec[s]);
	bzh01Count++;
	bzh01Sum += bzh01Vec[s];
	if (bzh01Count % 20 == 0 && bzh01StampVec[s].GetSec() > 1500e6) {
	  grbzh01RunAvg->SetPoint(grbzh01Run->GetN(), bzh01StampVec[s], bzh01Sum/20.);
	  bzh01Sum = 0.;
	} // if (bzh01Count % 20)
      } // if (stampVec[s] > utofStartStamp && stampVec[s] < utofEndStamp)
    } // for (int s=0; s<stampVec.size(); ++s)
    mgMagnet->Add(grbzh01Run);
    mgMagnetAvg->Add(grbzh01RunAvg);
    mgbzh01Avg->Add(grbzh01RunAvg);
    mgbzh01->Add(grbzh01Run);
    fout->cd();
    grbzh01Run->Write(Form("grbzh01Run_%s", blockVec[b]));
    grbzh01RunAvg->Write(Form("grbzh01RunAvg_%s", blockVec[b]));

    // Same thing for BZH03
    int bzh03Count = 1;
    double bzh03Sum = 0;
    for (int s=0; s<bzh03Vec.size(); ++s) {
      if (bzh03StampVec[s] > utofEndStamp) break;
      if (bzh03StampVec[s] < utofStartStamp) continue;

      if (bzh03StampVec[s] > utofStartStamp && bzh03StampVec[s] < utofEndStamp) {
	grbzh03Run->SetPoint(grbzh03Run->GetN(), bzh03StampVec[s], bzh03Vec[s]);
	bzh03Count++;
	bzh03Sum += bzh03Vec[s];
	if (bzh03Count % 20 == 0  && bzh03StampVec[s].GetSec() > 1500e6) {
	  grbzh03RunAvg->SetPoint(grbzh03Run->GetN(), bzh03StampVec[s], bzh03Sum/20.);
	  bzh03Sum = 0.;
	} // if (bzh03Count % 20)
      } // if (bzh03StampVec[s] > utofStartStamp && bzh03StampVec[s] < utofEndStamp)
    } // for (int s=0; s<bzh03StampVec.size(); ++s)
    mgMagnet->Add(grbzh03Run);
    mgMagnetAvg->Add(grbzh03RunAvg);
    mgbzh03Avg->Add(grbzh03RunAvg);
    mgbzh03->Add(grbzh03Run);
    fout->cd();
    grbzh03Run->Write(Form("grbzh03Run_%s", blockVec[b]));
    grbzh03RunAvg->Write(Form("grbzh03RunAvg_%s", blockVec[b]));

    // Find dtof runs
    int runMin = -1;
    int runMax = -1;
    int startTime = startTimeUtof;
    int endTime   = endTimeUtof;
    for (int irun=890; irun<1400; irun++) {
      TFile *fin = new TFile(Form("%srun%d/DsTOFcoincidenceRun%d_tdc1.root", dstofDir, irun, irun), "read");
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
    cout << "Min and max dtof runs are " << runMin << " " << runMax << endl;
    int nDtof = 1;

    lastScintT = 300.;
    vector<double> s1s2VecT;
    int nS1S2T = 0;
    for (int irun = runMin; irun < runMax+1; irun++) {
      TFile *dbFile = new TFile(Form("%s/spillDB_run%d_run%d.root", spillDir, irun, irun), "read");
      TTree *spillTree = (TTree*)dbFile->Get("spillTree");
      double globalSpillTime;
      double ustofSpillTime;
      spillTree->SetBranchAddress("globalSpillTime", &globalSpillTime);
      spillTree->SetBranchAddress("ustofSpillTime", &ustofSpillTime);
      TFile *fdtof = new TFile(Form("%s/run%d/DsTOFtreeRun%d_tdc1.root", dstofDir, irun, irun), "read");
      RawDsTofHeader *tof = NULL;
      TTree *tofTree = (TTree*)fdtof->Get("tofTree");
      tofTree->SetBranchAddress("tof", &tof);
      int lasts = 0;
      tofTree->GetEntry(0);
      int firstTime = tof->unixTime;
      double lastS1S2Dtof = 0.;
      for (int t = 0; t < spillTree->GetEntries(); t++) {
	spillTree->GetEntry(t);
	int hits = 0;
	if (globalSpillTime >= startTime && globalSpillTime <= endTime) {

	  for (int s = lasts; s<tofTree->GetEntries(); s++) {
	    tofTree->GetEntry(s);
	    if ((tof->fakeTimeNs/1e9)+firstTime >= globalSpillTime+1.) break;
	    if ((tof->fakeTimeNs/1e9)+firstTime <= globalSpillTime) continue;

	    if ((tof->fakeTimeNs/1e9)+firstTime >= globalSpillTime &&
		(tof->fakeTimeNs/1e9)+firstTime <= globalSpillTime + 1. && 
		tof->channel == 13 && (tof->fakeTimeNs - lastS1S2Dtof) > 500.) {
	      lastS1S2Dtof = tof->fakeTimeNs;
	      hits++;
	      lasts = s;
	      nS1S2T++;
	      if ((tof->fakeTimeNs/1e9)+firstTime > startTime + lastScintT) {
		lastScintT += 300.;
		s1s2VecT.push_back(nS1S2T);
		nS1S2T = 0;
	      }
	    }
	  } // for (int s = lasts; s<tofTree->GetEntries(); s++)
	} // if (globalSpillTime >= startTime && globalSpillTime <= endTime)
	grDtof->SetPoint(grDtof->GetN(), globalSpillTime, hits);
	nDtof++;

      } // for (int t = 0; t < spillTree->GetEntries(); t++)
      fdtof->Close();
      delete fdtof;
      dbFile->Close();
      delete dbFile;
    } // for (int irun = runMin; irun < runMax+1; irun++) 
    TGraph *grScintVsS1S2 = new TGraph(s1s2VecT.size(), &s1s2VecT[0], &countVecT[0]);
    // In this case check the separate populations
    if (b == 0) {
      // Go through vector and output into the correct graph
      std::vector<double> s1s2VecT_1;
      std::vector<double> s1s2VecT_2;
      std::vector<double> countVecT_1;
      std::vector<double> countVecT_2;
      for (int s1 = 0; s1 < s1s2VecT.size(); s1++) {
	if (stampVec.at(s) < cutTime) {
	  s1s2VecT_1.push_back(s1s2VecT.at(s));
	  countVecT_1.push_back(countVecT.at(s));
	}
	else {
	  s1s2VecT_2.push_back(s1s2VecT.at(s));
	  countVecT_2.push_back(countVecT.at(s));
	}
      } // for (int s1 = 0; s1 < s1s2VecT.size(); s1++)
      TGraph *grScintVsS1S2_1 = new TGraph(s1s2VecT_1.size(), &s1s2VecT_1[0], countVecT_1[0]);
      TGraph *grScintVsS1S2_2 = new TGraph(s1s2VecT_2.size(), &s1s2VecT_2[0], countVecT_2[0]);
      fout->cd();
      grScintVsS1S2_1->Write(Form("grScintVsS1S2_%s_1", blockVec[b]));
    }
    fout->cd();
    grScintVsS1S2->Write(Form("grScintVsS1S2_%s", blockVec[b]));
  } // for (int b=0; b<blockVec.size(); b++)
  mg->Add(grDtof);
  mg->SetTitle("S1 #cap S2 dtof hits for selected 4 block utof runs; Time / s; Hits / spill");
  TCanvas *c1 = new TCanvas("c1");
  mg->Draw("AP");

  fout->cd();
  grScintAll->Write("grScintAll");
  mg->Write("mg");
  mgScint->Write("mgScint");
  mgScintAvg->Write("mgScintAvg");
  mgMagnet->Write("mgMagnet");
  mgMagnetAvg->Write("mgMagnetAvg");
  mgbzh01->Write("mgbzh01");
  mgbzh03->Write("mgbzh03");
  mgbzh01Avg->Write("mgbzh01Avg");
  mgbzh03Avg->Write("mgbzh03Avg");

  fout->Close();
  delete fout;

} // block4Times
