/* 
 * RawDsTOFCoincidence.h     HP TPC raw Ds TOF coincidence
 * 
 *   Stores infromation about Ds TOF hit PMTs
 */

#ifndef HPTPC_RAW_DSTOF_COINCIDENCE_H
#define HPTPC_RAW_DSTOF_COINCIDENCE_H

#include "TObject.h" 


class RawDsTofCoincidence : public TObject
{
 public:

  RawDsTofCoincidence();
  ~RawDsTofCoincidence(); 
  
  Short_t run;
  Short_t tdc; 
  Short_t bar;
  Bool_t inSpill;
  Double_t fakeTimeNs[2];
  UInt_t unixTime[2]; 
  Double_t lastRawBeamSignal;
  Double_t lastDelayedBeamSignal;
  Double_t usTofSignal;

  ClassDef(RawDsTofCoincidence, 4); 

}; 

      

#endif
