// deadtime.C

void deadtime(const char* saveDir,
	      const char* dstofDir="/scratch0/dbrailsf/temp/mylinktodtof",
	      const char* ustofDir="/zfs_home/sjones/mylinktoutof") 
{
  gROOT->SetBatch(kTRUE);
  // Define the runs to be used for varying number of blocks
  const char* str0Block = "Data_2018_8_31_b2_800MeV_0block.root";
  const char* str1Block = "Data_2018_9_1_b4_800MeV_1block_bend4cm.root";
  const char* str2Block = "Data_2018_9_1_b2_800MeV_2block_bend4cm.root";
  const char* str3Block = "Data_2018_9_1_b3_800MeV_3block_bend4cm.root";
  const char* str4Block = "Data_2018_9_1_b8_800MeV_4block_bend4cm.root";
  // Unix timestamps for variable block moves
  // 0.8GeV/c, 0 blocks
  // 31/08/2018
  const double start0Block = 1535713289;
  const double end0Block   = 1535716132;
  // 0.8GeV/c, 1 block
  // 01/09/2018
  const double start1Block = 1535796057;
  const double end1Block   = 1535799112;
  // 0.8GeV/c, 2 blocks
  // 01/09/2018
  const double start2Block = 1535789157;
  const double end2Block   = 1535792026;
  // 0.8GeV/c, 3 block
  // 01/09/2018
  const double start3Block = 1535792404;
  const double end3Block   = 1535795300;
  // const double end3Block   = 1535798437;
  // 0.8GeV/c, 4 block
  // 4 moderator blocks with -4cm bend
  const double start4Block = 1535836129;
  const double end4Block   = 1535879634;

  
  for (int nBlocks = 0; nBlocks <=4; nBlocks++) {
    // Open files by hand since we're doing this locally
    if (nBlocks == 0) {
      
    }
    else if (nBlocks == 1) {

    }
    else if (nBlocks == 2) {

    }
    else if (nBlocks == 3) {

    }
    else if (nBlocks == 4) {

    }
  } // for (int nBlocks = 0; nBlocks <=4; nBlocks++)
  
} // deadtime
