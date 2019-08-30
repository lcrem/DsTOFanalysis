const double piMass   = 0.1396;
const double proMass  = 0.9383;
const double muMass   = 0.1057;
const double kMass    = 0.4937;
const double deutMass = 1.8756;

Double_t timeCalc(const Double_t mass, const Double_t dist, const Double_t mom)
{
  double t = dist * TMath::Sqrt( pow(mass,2)/pow(mom*3e8,2) + 1/pow(3e8,2) ) * 1e9;
  return t;
}

Double_t tPiS3(Double_t *P, Double_t *par)
{
  Float_t Pp = P[0];
  Double_t t = timeCalc(piMass, 10.9, Pp);
  return t;
}

Double_t tProS3(Double_t *P, Double_t *par)
{
  Float_t Pp = P[0];
  Double_t t = timeCalc(proMass, 10.9, Pp);
  return t;
}

Double_t tMuS3(Double_t *P, Double_t *par)
{
  Float_t Pp = P[0];
  Double_t t = timeCalc(muMass, 10.9, Pp);
  return t;
}

Double_t tKS3(Double_t *P, Double_t *par)
{
  Float_t Pp = P[0];
  Double_t t = timeCalc(kMass, 10.9, Pp);
  return t;
}

Double_t tDeutS3(Double_t *P, Double_t *par)
{
  Float_t Pp = P[0];
  Double_t t = timeCalc(deutMass, 10.9, Pp);
  return t;
}

Double_t tPiS4(Double_t *P, Double_t *par)
{
  Float_t Pp = P[0];
  Double_t t = timeCalc(piMass, 12.6, Pp);
  return t;
}

Double_t tProS4(Double_t *P, Double_t *par)
{
  Float_t Pp = P[0];
  Double_t t = timeCalc(proMass, 12.6, Pp);
  return t;
}

Double_t tMuS4(Double_t *P, Double_t *par)
{
  Float_t Pp = P[0];
  Double_t t = timeCalc(muMass, 12.6, Pp);
  return t;
}

Double_t tKS4(Double_t *P, Double_t *par)
{
  Float_t Pp = P[0];
  Double_t t = timeCalc(kMass, 12.6, Pp);
  return t;
}

Double_t tDeutS4(Double_t *P, Double_t *par)
{
  Float_t Pp = P[0];
  Double_t t = timeCalc(deutMass, 12.6, Pp);
  return t;
}

void tofVsMom(const char* outDir) {
  TFile *fout = new TFile(Form("%s/tofVsMom.root", outDir), "recreate");

  TF1 *fPiS3 = new TF1("fPiS3", tPiS3, 0.1, 1.0, 0);
  fPiS3->SetLineColor(kCyan+1);
  TF1 *fProS3 = new TF1("fProS3", tProS3, 0.1, 1.0, 0);
  fProS3->SetLineColor(kRed);
  TF1 *fMuS3 = new TF1("fProS3", tMuS3, 0.1, 1.0, 0);
  fMuS3->SetLineColor(kBlue);
  TF1 *fKS3 = new TF1("fKS3", tKS3, 0.1, 1.0, 0);
  fKS3->SetLineColor(kOrange+1);
  TF1 *fDeutS3 = new TF1("fDeutS3", tDeutS3, 0.1, 1.0, 0);
  fDeutS3->SetLineColor(kBlack);
  
  TF1 *fPiS4 = new TF1("fPiS4", tPiS4, 0.1, 1.0, 0);
  fPiS4->SetLineColor(kCyan+1);
  TF1 *fProS4 = new TF1("fProS4", tProS4, 0.1, 1.0, 0);
  fProS4->SetLineColor(kRed);
  TF1 *fMuS4 = new TF1("fProS4", tMuS4, 0.1, 1.0, 0);
  fMuS4->SetLineColor(kBlue);
  TF1 *fKS4 = new TF1("fKS4", tKS4, 0.1, 1.0, 0);
  fKS4->SetLineColor(kOrange+1);
  TF1 *fDeutS4 = new TF1("fDeutS4", tDeutS4, 0.1, 1.0, 0);
  fDeutS4->SetLineColor(kBlack);

  std::vector<TF1*> s3Vec = {fPiS3, fProS3, fMuS3, fKS3, fDeutS3};
  std::vector<TF1*> s4Vec = {fPiS4, fProS4, fMuS4, fKS4, fDeutS4};

  TCanvas *c1 = new TCanvas("c1", "c1");
  TCanvas *c2 = new TCanvas("c2", "c2");
  TLegend *leg = new TLegend(0.67, 0.55, 0.9, 0.9);
  for (int i=0; i<s3Vec.size(); i++) {
    s3Vec.at(i)->SetLineWidth(2);
    s4Vec.at(i)->SetLineWidth(2);
    s3Vec.at(i)->GetHistogram()->GetYaxis()->SetRangeUser(35., 250.);
    s4Vec.at(i)->GetHistogram()->GetYaxis()->SetRangeUser(35., 250.);
    s3Vec.at(i)->GetHistogram()->GetXaxis()->SetTitleSize(0.05);
    s3Vec.at(i)->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
    s3Vec.at(i)->GetHistogram()->GetYaxis()->SetTitleSize(0.05);
    s3Vec.at(i)->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
    s4Vec.at(i)->GetHistogram()->GetXaxis()->SetTitleSize(0.05);
    s4Vec.at(i)->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
    s4Vec.at(i)->GetHistogram()->GetYaxis()->SetTitleSize(0.05);
    s4Vec.at(i)->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
    c1->cd();
    if (i==0) {
      s3Vec.at(i)->GetHistogram()->SetTitle("Predicted time of flight from S1 to S3; Particle momentum [GeV/c]; Time of flight / ns");
      s3Vec.at(i)->Draw();
    }
    else s3Vec.at(i)->Draw("same");
    c2->cd();
    if (i==0) {
      s4Vec.at(i)->GetHistogram()->SetTitle("Predicted time of flight from S2 to S4; Particle momentum [GeV/c]; Time of flight / ns");
      s4Vec.at(i)->Draw();
    }
    else s4Vec.at(i)->Draw("same");
    fout->cd();
    s3Vec.at(i)->Write();
    s4Vec.at(i)->Write();
  }

  leg->AddEntry(fPiS3, "Pion", "l");
  leg->AddEntry(fMuS3, "Muon", "l");
  leg->AddEntry(fKS3, "Kaon", "l");
  leg->AddEntry(fProS3, "Proton", "l");
  leg->AddEntry(fDeutS3, "Deuteron", "l");

  c1->cd();
  leg->Draw();
  c2->cd();
  leg->Draw();
  c1->SetLeftMargin(0.13);
  c1->SetBottomMargin(0.13);
  c2->SetLeftMargin(0.13);
  c2->SetBottomMargin(0.13);
  c1->Print(Form("%s/s1s3Times.pdf", outDir));
  c1->Print(Form("%s/s1s3Times.png", outDir));
  c1->Print(Form("%s/s1s3Times.tex", outDir));
  c2->Print(Form("%s/s2s4Times.pdf", outDir));
  c2->Print(Form("%s/s2s4Times.png", outDir));
  c2->Print(Form("%s/s2s4Times.tex", outDir));
  
  fout->Close();
  delete fout;
}
