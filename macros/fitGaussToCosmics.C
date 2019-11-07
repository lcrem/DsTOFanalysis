// fitGaussToCosmics.C
// Attempt to fit gaussians to 1D projections of the bar by bar cosmics

void fitGaussToCosmics(const char* outFile, const char* s4Plots="/scratch0/sjones/plots/angularDistS4_newSample/withTree/includesUnweightedHists.root") 
{
  gROOt->SetBatch(kTRUE);
  TFile *fout = new TFile(outFile, "recreate");

  TFile *fin = new TFile(

  fout->Close();
  delete fout;
} // fitGaussToCosmics
