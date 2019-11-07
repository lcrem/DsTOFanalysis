// threshRes.C
// Check how bar resolution and time difference changes with a changing threshold
// Here we are using as the trigger, the first time a signal goes over a certain point

std::vector<const char*> days = {/*"05", "06", */"09", "10", "11", "12", 
				 "13", "16",/* "17",*/ "19", "20", "23"};
std::vector<double> thresholds = {5, 10, 15, 20};
std::vector<const char*> trigCh     = {"3", "4"};
std::vector<const char*> readCh     = {"3", "4"};
std::vector<int> sourceDist   = {16, 32, 48, 64, 80, 96, 112};

void threshRes(const char *outFile, const char *fileDir="/unix/dune/tof/")
{
  TFile *fout = new TFile(outFile, "recreate");

  for (int i = 0; i<days.size(); i++) {
    void* dirp = gSystem->OpenDirectory(Form("%s/2018Jul%s", fileDir, days.at(i)));
    const char* entry;
    TString barFile;
    while ((entry = (char*)gSystem->GetDirEntry(dirp))) {
      barFile = entry; // The bar numbers
      // cout<<barFile<<endl;
      // Append bar name, source position, threshold, readout channel and trigger channel
      if (barFile.Contains("BarD")) {
	for (int t=0; t<thresholds.size(); t++) {
	  for (int tch=0; tch<trigCh.size(); tch++) { 
	    TGraph *grResDist = new TGraph();
	    TGraph *grDiffDist = new TGraph();
	    grResDist->SetName(Form("grResDist_%s_thresh%dmV_trig%s", barFile.Data(), (int)(thresholds.at(t)), trigCh.at(tch)));
	    grResDist->SetTitle(Form("%s: resolution as function of distance; Position / cm; Resolution / s", barFile.Data()));
	    grDiffDist->SetName(Form("grResDist_%s_thresh%dmV_trig%s", barFile.Data(), (int)(thresholds.at(t)), trigCh.at(tch)));
	    grDiffDist->SetTitle(Form("%s: #Deltat as function of distance; Position / cm; #Deltat / s", barFile.Data()));
	    for (int src=0; src<sourceDist.size(); src++) {
	      TFile *finch3 = new TFile(Form("%s/2018Jul%s/%s/SourceAt%dcm/TrigCh%s_thres%dmV_20mVdiv.ch3.traces.root", fileDir, days.at(i), barFile.Data(), sourceDist.at(src), trigCh.at(tch), (int)(thresholds.at(t))), "read");
	      TFile *finch4 = new TFile(Form("%s/2018Jul%s/%s/SourceAt%dcm/TrigCh%s_thres%dmV_20mVdiv.ch4.traces.root", fileDir, days.at(i), barFile.Data(), sourceDist.at(src), trigCh.at(tch), (int)(thresholds.at(t))), "read");
	      // Files are good
	      if (!finch3->IsZombie() && !finch4->IsZombie()) {
		TH1D *hDiff = new TH1D(Form("hDiff_%s_thresh%dmV_trig%s_src%d", barFile.Data(), (int)(thresholds.at(t)), trigCh.at(tch), sourceDist.at(src)), Form("#Delta t, %s, threshold %dmV, distance %dcm; Ch3 - Ch4 / s", barFile.Data(), (int)(thresholds.at(t)), sourceDist.at(src)), 100, -20e-9, 20e-9);
		// For each of the files, find the first point when the graph goes under the threshold
		int nEntries = 0;
		TKey *key;
		TIter next(finch3->GetListOfKeys());
		std::cout<<Form("%s/2018Jul%s/%s/SourceAt%dcm/TrigCh%s_thres%dmV_20mVdiv.ch4.traces.root", fileDir, days.at(i), barFile.Data(), sourceDist.at(src), trigCh.at(tch), (int)(thresholds.at(t)))<<" is good file"<<std::endl;
		while ((key = (TKey*) next())) {
		  nEntries++;

		  TGraph *grch3 = (TGraph*)finch3->Get(Form("graph%d", nEntries));
		  TGraph *grch4 = (TGraph*)finch4->Get(Form("graph%d", nEntries));
		  double t3, v3, t4, v4 = 0;
		  double trigTime3, trigTime4 = 0;
		  for (int p3=0; p3<grch3->GetN(); p3++) {
		    grch3->GetPoint(p3, t3, v3);
		    if (v3 < -0.001 * thresholds.at(t)) {
		      trigTime3 = t3;
		      break;
		    } // Is below trigger threshold
		  } // Loop over ch3 points
		  for (int p4=0; p4<grch4->GetN(); p4++) {
		    grch4->GetPoint(p4, t4, v4);
		    if (v4 < -0.001 * thresholds.at(t)) {
		      trigTime4 = t4;
		      break;
		    } // Is below trigger threshold
		  } // Loop over ch4 points
		  // Time difference 
		  if (trigTime3 != 0 && trigTime4 != 0) hDiff->Fill(trigTime3 - trigTime4);
		  delete grch3;
		  delete grch4;
		} // Loop through keys in file
		fout->cd();
		TF1 *f = new TF1("f", "gaus");
		hDiff->Fit(f, "Q");
		grResDist->SetPoint(grResDist->GetN(), sourceDist.at(src), f->GetParameter(2));
		grDiffDist->SetPoint(grDiffDist->GetN(), sourceDist.at(src), f->GetParameter(1));
		hDiff->Write();
		delete f;
		finch3->Close();
		finch4->Close();
	      } // Files are good
	      else {
		std::cout<<"One or both of "<<Form("%s/2018Jul%s/%s/SourceAt%dcm/TrigCh%s_thres%dmV_20mVdiv.ch3.traces.root", fileDir, days.at(i), barFile.Data(), sourceDist.at(src), trigCh.at(tch), (int)(thresholds.at(t)))<<" or "<<Form("%s/2018Jul%s/%s/SourceAt%dcm/TrigCh%s_thres%dmV_20mVdiv.ch3.traces.root", fileDir, days.at(i), barFile.Data(), sourceDist.at(src), trigCh.at(tch), (int)(thresholds.at(t)))<<" does not exist"<<std::endl;
	      }

	      delete finch3;
	      delete finch4;
	    } // Loop over distances
	    fout->cd();
	    grResDist->Write();
	    grDiffDist->Write();
	  } // Loop over trigger channels
	} // Loop over thresholds
      } // Correct bar format
    } // Loop through directories
  } // Loop over days
  fout->Close();
  delete fout;
} // threshRes
