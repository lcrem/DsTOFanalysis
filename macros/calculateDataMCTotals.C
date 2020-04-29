// calculateDataMCTotals.C
// Calculate Data and MC ratios with all appropriate systematics

// Numbers
// Data numbers are per spill, 
const std::vector<double> nMCS4   = {1364, 4113, 4678, 973, 1076}; 
const std::vector<double> nDataS4 = {122, 185, 167, 49, 5.8};
const std::vector<double> nMCS3   = {49690, 61676, 55944, 17587, 93930};
const std::vector<double> nDataS3 = {2227, 1881, 1473, 982, 148}; 
// Stats errors
const std::vector<double> nDataErrS4 = {2, 2, 2, 2, 0.2};
const std::vector<double> nDataErrS3 = {10, 7, 6, 7, 0.5};
// Fractional systematics errors
const std::vector<double> nMCS4Syst   = {0.095, 0.08, 0.085, 0.171, 0.08};
const std::vector<double> nDataS4Syst = {0.029, 0.015, 0.067, 0.082, 0.041};

void calculateDataMCTotals()
{
  for (int b=0; b<5; b++) {
    std::cout<<"========================================="<<std::endl;
    std::cout<<b<<" blocks"<<std::endl;

    double dataRat = nDataS4[b]/nDataS3[b];
    double MCRat   = nMCS4[b]/nMCS3[b];
    double dataRatFracErr = sqrt( pow(nDataErrS4[b]/nDataS4[b],2) + pow(nDataErrS3[b]/nDataS3[b],2) + nDataS4Syst[b]*nDataS4Syst[b] );
    double MCRatFracErr   = sqrt( 1./nMCS4[b] + 1./nMCS3[b] + nMCS4Syst[b]*nMCS4Syst[b] );

    double dataMC        = dataRat/MCRat;
    double dataMCFracErr = sqrt(dataRatFracErr*dataRatFracErr + MCRatFracErr*MCRatFracErr);
    
    std::cout<<"S4/S3: Data = "<<dataRat<<" +- "<<dataRat*dataRatFracErr<<",  MC = "<<MCRat<<" +- "<<MCRat*MCRatFracErr<<std::endl;
    std::cout<<"Data/MC = "<<dataMC<<" +- "<<dataMC*dataMCFracErr<<std::endl;
  }
  
} // calculateDataMCTotals()
