// calculateDataMCTotals.C
// Calculate Data and MC ratios with all appropriate systematics

// Numbers
// Data numbers are per spill, 
const std::vector<double> nMCS4   = {1364, 4113, 4678, 973, 1076}; 
const std::vector<double> nDataS4 = {119, 181, 157.9, 38.9, 1.32};
const std::vector<double> nMCS3   = {49690, 61676, 55944, 17587, 93930};
const std::vector<double> nDataS3 = {2441, 2058, 1608, 1069, 166.5}; 
// Stats errors
const std::vector<double> nDataErrS4 = {2.1, 2.4, 1.9, 1.4, 0.15};
const std::vector<double> nDataErrS3 = {10, 8, 7, 8, 0.6};
// Fractional systematics errors
const std::vector<double> nMCS4Syst      = {0.095, 0.08, 0.085, 0.171, 0.08};
const std::vector<double> nDataS4AngSyst = {0.029, 0.015, 0.067, 0.082, 0.041};
const std::vector<double> nDataS4EffSyst = {0.10964, 0.10964, 0.10964, 0.10964, 0.10964};
const std::vector<double> nDataS4BkgSyst = {0.0018, 0.0016, 0.0114, 0.0141, 0.0811};
const std::vector<double> nDataS3EffSyst = {0.0113, 0.1143, 0.0698, 0.1143, 0.0492};

void calculateDataMCTotals()
{
  for (int b=0; b<5; b++) {
    std::cout<<"========================================="<<std::endl;
    std::cout<<b<<" blocks"<<std::endl;

    double dataRat = nDataS4[b]/nDataS3[b];
    double MCRat   = nMCS4[b]/nMCS3[b];
    double dataRatFracErr = sqrt( pow(nDataErrS4[b]/nDataS4[b],2) + pow(nDataErrS3[b]/nDataS3[b],2) + nDataS4AngSyst[b]*nDataS4AngSyst[b] + nDataS4EffSyst[b]*nDataS4EffSyst[b] + nDataS4BkgSyst[b]*nDataS4BkgSyst[b] + nDataS3EffSyst[b]*nDataS3EffSyst[b] );
    double MCRatFracErr   = sqrt( 1./nMCS4[b] + 1./nMCS3[b] + nMCS4Syst[b]*nMCS4Syst[b] );

    double dataMC        = dataRat/MCRat;
    double dataMCFracErr = sqrt(dataRatFracErr*dataRatFracErr + MCRatFracErr*MCRatFracErr);
    
    std::cout<<"S4/S3 MC = "<<MCRat<<" +- "<<MCRat*MCRatFracErr<<", Data = "<<dataRat<<" +- "<<dataRat*dataRatFracErr<<std::endl;
    std::cout<<"Data/MC = "<<dataMC<<" +- "<<dataMC*dataMCFracErr<<std::endl;

    std::cout<<"Data % error = "<<dataRatFracErr*100.<<std::endl;
  }
  
} // calculateDataMCTotals()
