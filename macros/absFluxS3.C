// absFluxS3.C
// Similar to absFluxS4.C but for S3
// Look at difference between only using S1 and using S1+S2

void absFluxS3(const char* saveDir,
	       const char* ustofDir = "/zfs_home/sjones/mylinktoutof/")
{
  gROOT->SetBatch(kTRUE);

  // Define the runs to be used for varying number of blocks
  const char* str0Block = "Data_2018_8_31_b2_800MeV_0block.root";
  const char* str1Block = "Data_2018_9_1_b4_800MeV_1block_bend4cm.root";
  const char* str2Block = "Data_2018_9_1_b2_800MeV_2block_bend4cm.root";
  const char* str3Block = "Data_2018_9_1_b3_800MeV_3block_bend4cm.root";
  const char* str4Block = "Data_2018_9_1_b8_800MeV_4block_bend4cm.root";

  THStack *hsS1 = new THStack("hsS1", "Angular distribution of S3 hits with S1 trigger only; #theta / degrees; Events / spill");
  THStack *hsS1S2 = new THStack("hsS1S2", "Angular distribution of S3 hits with S1 & S2; #theta / degrees; Events / spill");

  for (int nBlocks = 0; nBlocks <=4; nBlocks++) {
    char* nustof;
    if (nBlocks == 0) nustof = Form("%s/%s", ustofDir, str0Block);
    else if (nBlocks == 1) nustof = Form("%s/%s", ustofDir, str1Block);
    else if (nBlocks == 2) nustof = Form("%s/%s", ustofDir, str2Block);
    else if (nBlocks == 3) nustof = Form("%s/%s", ustofDir, str3Block);
    else if (nBlocks == 4) nustof = Form("%s/%s", ustofDir, str4Block);
      
  } // for (int nBlocks = 0; nBlocks <=4; nBlocks++)

}
