// mcToData.C
// Take the Monte Carlo and imply the inverse of the corrections + see what we get out
// Do everything in global coordinates. It's just easier that way
#include "UsefulFunctions.C"

const char* mcstem  = "/scratch0/tnonnenm/FixedMC/Prototype-HPTPC-MC/wdir/moderator_data";
const char* wgtstem = "/scratch0/sjones/plots/angularDistS3_newSample/s3ProtonWeightTrees.root";

const vector<double> dataS3 = {1983., 1656., 1325., 899., 136.3};

void mcToData(const char* outfile)
{
  gROOT->SetBatch(kTRUE);
  TFile *fout = new TFile(outfile, "recreate");

  for (int b=0; b<5; b++) {
    // Do data stuff first
    vector<double> startTimes;
    vector<double> endTimes;
    startTimes.clear();
    endTimes.clear();
    // For spill counting normalisation
    int nSpills = 0;
    int nSpillsTrue = 0;
    double lastSpill = 0.;

    // Find the correct dstof files
    Int_t runMin=-1;
    Int_t runMax=-1;
    double startTime = 0;
    double endTime   = 0;

    if (b == 0) {
      startTimes.push_back(start0Block);
      endTimes.push_back(end0Block);
    }
    else if (b == 1) {
      startTimes.push_back(start1Block);
      endTimes.push_back(end1Block);
    }
    else if (b == 2) {
      startTimes.push_back(start2Block);
      endTimes.push_back(end2Block);
    }
    else if (b == 3) {
      startTimes.push_back(start3Block);
      endTimes.push_back(end3Block);
    }
    else if (b == 4) {
      for (int b4=0; b4<str4BlockVec.size(); b4++) {
	TFile *futofTmp = new TFile(Form("%s/%s",ustofDir, str4BlockVec.at(b4)), "read");
	TTree *treeTmp = (TTree*)futofTmp->Get("tree");
	double tS1Tmp;
	treeTmp->SetBranchAddress("tS1", &tS1Tmp);
	TNamed *start = 0;
	TNamed *end   = 0;                                                                               
	futofTmp->GetObject("start_of_run", start);
	const char* startchar = start->GetTitle();
	std::string startstr(startchar);
	std::string unixstart = startstr.substr(25,10);
	int startTime = stoi(unixstart);
	treeTmp->GetEntry(treeTmp->GetEntries() - 1);
	int endTime = startTime + (tS1Tmp/1e9);
	futofTmp->Close();
	delete futofTmp;
	startTimes.push_back(startTime);
	endTimes.push_back(endTime);
      }
    }

    TH2D *h2CosmicsEff = new TH2D(Form("h2CosmicsEff%d", b), Form("Relative efficiency, %d blocks; x / cm; Bar; Eff", b), binnumCosmics, binsCosmics, 10, 0.5, 10.5);
    setHistAttr(h2CosmicsEff);
    // Loop over subsamples
    for (int sub=0; sub<startTimes.size(); sub++) {
      startTime = startTimes.at(sub);
      endTime   = endTimes.at(sub);

      for (int irun=950; irun<1400; irun++) {
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
      cout << "Min and max runs are " << runMin << " " << runMax << endl;

      fout->cd();
      h2CosmicsEff->Write();

      // Build cosmic maps first
      vector<double> barTimeVec;
      barTimeVec.resize(10, 0.);
      for (int itdc=0; itdc<2; itdc++) {
	unsigned int lastRun = 0;
	for (int irun=runMin; irun<runMax+1; irun++) {
	  // Load input files
	  TFile *tofCoinFile = new TFile(Form("%srun%d/DsTOFcoincidenceRun%d_tdc%d.root", dstofDir, irun, irun, itdc+1));
	  TTree *tofCoinTree = (TTree*)tofCoinFile->Get("tofCoinTree");
	  RawDsTofCoincidence *tofCoin = NULL;
	  tofCoinTree->SetBranchAddress("tofCoin", &tofCoin);

	  for (int t=0; t<tofCoinTree->GetEntries(); t++) {
	    tofCoinTree->GetEntry(t);
	    if (tofCoin->unixTime[0] < startTime) continue;
	    if (tofCoin->unixTime[0] > endTime) break;
	    if (abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1]) > s4BarTime) continue;
	    double deltat = TMath::Abs(tofCoin->fakeTimeNs[0]-tofCoin->fakeTimeNs[1]);
	    double dstofHitT = min(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]) - (s4BarTime/2. - TMath::Abs(deltat) / 2. );
	    if (dstofHitT - barTimeVec.at(tofCoin->bar-1) < s4DeadtimeCut && 
		tofCoin->run==lastRun) {
	      barTimeVec.at(tofCoin->bar-1) = dstofHitT;
	      continue;
	    }
	    else if (tofCoin->run != lastRun) {
	      barTimeVec.clear();
	      barTimeVec.resize(10, 0.);
	      lastRun = tofCoin->run;
	    }
	    barTimeVec.at(tofCoin->bar-1) = dstofHitT;
	    // If the hits are not in a spill then consider them cosmics 
	    // and use them for angular efficiency
	    if (!tofCoin->inSpill && 
		abs(tofCoin->fakeTimeNs[1] - tofCoin->fakeTimeNs[0]) < s4BarTime) {	      
	      double positionX = localDtofPosition(tofCoin->fakeTimeNs[0], tofCoin->fakeTimeNs[1]);
	      h2CosmicsEff->Fill(positionX, tofCoin->bar);
	    } 
	  } // Loop over coincidences
	} // Loop over runs
      } // Loop over TDCs
      h2CosmicsEff->Scale(1., "width");
      h2CosmicsEff->Scale(1./h2CosmicsEff->GetBinContent(h2CosmicsEff->GetMaximumBin()));
    } // Loop over data subsamples

    // Now do the MC stuff
    cout<<"======================="<<endl;
    cout<<b<<" blocks"<<endl;
    cout<<"======================="<<endl;
    TFile *fmc = new TFile(Form("%s/NewGeom%dBlocks.root", mcstem, b), "read");
    TTree *treemc = (TTree*)fmc->Get("stepTree");
    UInt_t NSteps;
    float Mom[60000];
    Int_t PreVolID[60000];
    Int_t PostVolID[60000];
    Int_t PDG[60000];
    Float_t PostX[60000], PostY[60000], PostZ[60000];
    treemc->SetBranchAddress("NSteps", &NSteps);
    treemc->SetBranchAddress("PostStepMom", Mom);
    treemc->SetBranchAddress("PDG", PDG);
    treemc->SetBranchAddress("PreStepVolumeID", PreVolID);
    treemc->SetBranchAddress("PostStepVolumeID", PostVolID);
    treemc->SetBranchAddress("PostStepPosX", PostX);
    treemc->SetBranchAddress("PostStepPosY", PostY);
    treemc->SetBranchAddress("PostStepPosZ", PostZ);
    treemc->AddFriend(Form("protonTree%dBlocks", b), wgtstem);
    int spill;
    double weight;
    treemc->SetBranchAddress("spill", &spill);
    treemc->SetBranchAddress("weight", &weight);

    double nS4MC = 0;
    double nS3MC = treemc->GetEntries();

    TH1D *hMcThetaUnwgt = new TH1D(Form("hMcThetaUnwgt%d", b), Form("%d blocks", b), nBinsS4Horz, binsS4HorzLow, binsS4HorzHigh);
    setHistAttr(hMcThetaUnwgt);
    TH1D *hMcTheta = new TH1D(Form("hMcTheta%d", b), Form("%d blocks", b), nBinsS4Horz, binsS4HorzLow, binsS4HorzHigh);
    setHistAttr(hMcTheta);
    TH1D *hMcPhiUnwgt = new TH1D(Form("hMcPhiUnwgt%d", b), Form("%d blocks", b), nBinsS4Vert, binsS4VertLow, binsS4VertHigh);
    setHistAttr(hMcPhiUnwgt);
    TH1D *hMcPhi = new TH1D(Form("hMcPhi%d", b), Form("%d blocks", b), nBinsS4Vert, binsS4VertLow, binsS4VertHigh);
    setHistAttr(hMcPhi);
    TH2D *hMcLocalUnwgt = new TH2D(Form("h2McLocalUnwgt%d", b), Form("%d blocks; Bar position / cm; Bar; Events", b), 20, 0., 140., 10, 0.5, 10.5);

    cout<<"Looping over MC files"<<endl;
    for  (int e=0; e<treemc->GetEntries(); e++) {
      treemc->GetEntry(e);
      for (int n=0; n<NSteps; n++) {
	if (PreVolID[n]!=7 && PostVolID[n]==7 && Mom[n]>30. && PDG[n]==2212) {
	  nS4MC++;
	  TVector3 mcCoords(PostX[n]/1000., PostY[n]/1000., PostZ[n]/1000.);
	  TVector3 globalCoords = MCToGlobalCoords(mcCoords);
	  double localpos = localDtofFromMC(mcCoords);
	  int bar = dtofBarFromMC(mcCoords);
	  double theta = getThetaFromGlobal(globalCoords);
	  double phi   = getPhiFromGlobal(globalCoords);
	  hMcLocalUnwgt->Fill(localpos, bar);
	  hMcThetaUnwgt->Fill(theta);
	  hMcPhiUnwgt->Fill(phi);
	  double w = h2CosmicsEff->GetBinContent(localpos, bar);
	  hMcTheta->Fill(theta, w);
	  hMcPhi->Fill(phi, w);
	  break;
	} // Selection conditions for S4 protons
      } // Loop over steps
    } // Loop over MC entries

    cout<<"Unweighted MC S4/S3 = "<<nS4MC/nS3MC<<endl;

    fout->cd();
    hMcThetaUnwgt->Write();
    hMcPhiUnwgt->Write();
    hMcTheta->Write();
    hMcPhi->Write();
    hMcLocalUnwgt->Write();

    fmc->Close();
    delete fmc;
  } // Loop over blocks

  fout->Close();
  delete fout;
} // mcToData

