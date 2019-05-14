// getScintillatorValues.C

// Time format below
// "2018-08-15T00:00:00";
// "2018-09-18T19:00:00";
TTimeStamp convert_to_timestamp(TString time_string){

  TString time_to_convert = time_string;
  char* time =  (char*)time_to_convert.Data();

  //Get various time parameters from the char that the query spits out  
  UInt_t year = atoi(strtok(time, "-"));
  UInt_t month = atoi(strtok(NULL, "-"));
  UInt_t day = atoi(strtok(NULL, " "));
  UInt_t hour = atoi(strtok(NULL, ":"));
  UInt_t min = atoi(strtok(NULL, ":"));
  UInt_t sec = atoi(strtok(NULL, "")); 

  TTimeStamp time_stamp = TTimeStamp(year, month, day, hour, min, sec, 0, "kFALSE", 0);

  return time_stamp;
}

void getScintillatorValues(const char* txtFile="/scratch2/sjones/Scintillator_info.txt",
			   const char* utofDir="/zfs_home/sjones/mylinktoutof/") 
{
  const char* str4Block1 = "Data_2018_9_1_b8_800MeV_4block_bend4cm.root";
  const char* str4Block2 = "Data_2018_9_3_b2_800MeV_4block_bend4cm.root";
  const char* str4Block3 = "Data_2018_9_1_b5_800MeV_4block.root";
  const char* str4Block4 = "Data_2018_9_4_b1_800MeV_4block.root";
  const char* str4Block5 = "Data_2018_9_4_b2_800MeV_4block.root";
  const char* str4Block6 = "Data_2018_9_3_b4_800MeV_4block_1mUSB.root";
  const char* str0Block = "Data_2018_8_31_b2_800MeV_0block.root";
  const char* str1Block = "Data_2018_9_1_b4_800MeV_1block_bend4cm.root";
  const char* str2Block = "Data_2018_9_1_b2_800MeV_2block_bend4cm.root";
  const char* str3Block = "Data_2018_9_1_b3_800MeV_3block_bend4cm.root";
  std::vector<const char*> = {str4Block1, str4Block2, str4Block3, str4Block4,
			      str4Block5, str4Block6, str0Block, str1Block, 
			      str2Block, str3Block};
  std::vector<int> countVec;
  std::vector<TTimeStamp> stampVec;
  TGraph *grValues = new TGraph();
  string line;
  ifstream inFile;
  inFile.open(txtFile);
  while(!inFile.eof()) {
    getline(inFile, line);
    if (line.length()>20) {
      TString inLine(line);
      char* charLine = (char*)inLine.Data();
      cout<<line<<endl;
      // cout<<line.substr(0,4)<<" "<<line.substr(5,2)<<endl;
      stringstream ssyear(line.substr(0,4));
      stringstream ssmonth(line.substr(5,2));
      stringstream ssday(line.substr(8,2));
      stringstream sshour(line.substr(11,2));
      stringstream ssmin(line.substr(14,2));
      stringstream sssec(line.substr(17,2));
      stringstream sscount(line.substr(24));
      int year = 0;//stoi(line.substr(0, 4), &sz);
      int month = 0;//stoi(line.substr(5,2), &sz);
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
      sscount >> count;
      cout<<year<<" "<<month<<" "<<day<<" "<<hour<<" "<<min<<" "<<sec<<endl;
      TTimeStamp stamp(year, month, day, hour, min, sec, 0, "kFALSE", 7200);
      countVec.push_back(count);
      stampVec.push_back(stamp);

      if (count <20000) {
	grValues->SetPoint(grValues->GetN(), stamp, count);
      }
    }

  }
  inFile.close();

  grValues->Draw("AP*");
} // getScintillatorValues
