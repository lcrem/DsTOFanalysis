// 3BlockMcVesselTest.C
// Vary the vessel thickness slightly for the 3 block case to see what the effect is on the number
// of protons
#include "UsefulFunctions.C"

const char* stem = "/scratch0/tnonnenm/FixedMC/Prototype-HPTPC-MC/wdir/moderator_data";

vector<const char*> fileNames = {"3Blockvesselminus1mm", "3Blockvesselminus05mm", 
				 "3Blockvesselminus02mm", "NewGeom3Blocks",
				 "3Blockvesselplus02mm", "3Blockvesselplus05mm", 
				 "3Blockvesselplus1mm"};
				 
vector<const char*> descriptions = {"-1 mm", "-0.5 mm", "-0.2 mm",  
				    "Nominal", 
				    "+0.2 mm", "+0.5 mm", "+1 mm"};

vector<double> values = {-1., -0.5, -0.2, 0., 0.2, 0.5, 1.};

void ThreeBlockMcVesselTest(const char* outfile)
{
  gROOT->SetBatch(kTRUE);

  TFile *fout = new TFile(outfile, "recreate");
  fout->cd();

  TGraphErrors *gS4     = new TGraphErrors();
  TGraphErrors *gChange = new TGraphErrors();
  setHistAttr(gS4, kRed);
  setHistAttr(gChange, kRed);
  gS4->SetTitle("Protons reaching S4, 3 blocks; Change in vessel thickness / mm; # protons");
  gChange->SetTitle("Fractional change in protons reaching S4, 3 blocks; Change in vessel thickness / mm; Fractional change");

  for (int f=0; f<fileNames.size(); f++) {
    cout<<"=================="<<endl;
    cout<<f<<" sample"<<endl;
    cout<<"=================="<<endl;

    TFile *fin = new TFile(Form("%s/%s.root", stem, fileNames[f]), "read");
    TTree *t = (TTree*)fin->Get("stepTree");
    UInt_t NSteps;
    float Mom[50000];
    int PreVolID[50000];
    int PostVolID[50000];
    int PDG[50000];
    t->SetBranchAddress("NSteps", &NSteps);
    t->SetBranchAddress("PostStepMom", Mom);
    t->SetBranchAddress("PreStepVolumeID", PreVolID);
    t->SetBranchAddress("PostStepVolumeID", PostVolID);
    t->SetBranchAddress("PDG", PDG);

    double nS4 = 0;
    for(int e=0; e<t->GetEntries(); e++) { 
      t->GetEntry(e);
      for (UInt_t n=0; n<NSteps; n++) {
	if (PreVolID[n]!=7 && PostVolID[n]==7 && Mom[n] > 140. && PDG[n]==2212) {
	  nS4++;
	  break;
	} // Is an S4 proton
      } // Loop through steps
    }
    gS4->SetPoint(f, values[f], nS4);
    gS4->SetPointError(f, 0., sqrt(nS4));

    fin->Close();
    delete fin;
  } // Loop over files

  for (int p=0; p<fileNames.size(); p++) {
    double x, y;
    double xOrig, yOrig;
    gS4->GetPoint(p, x, y);
    gS4->GetPoint(3, xOrig, yOrig);
    double errSum = sqrt(yOrig + y);
    double err = abs((yOrig - y)/yOrig) * sqrt( pow(errSum/(yOrig - y), 2) + pow(sqrt(yOrig)/yOrig, 2) );
    gChange->SetPoint(p, values[p], (yOrig - y)/yOrig);
    gChange->SetPointError(p, 0., err);
  }

  fout->cd();
  gS4->Write("gS4");
  gChange->Write("gChange");

  fout->Close();
  delete fout;
}

