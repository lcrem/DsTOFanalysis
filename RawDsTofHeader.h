/* 
 * RawDsTOFHeader.h     HP TPC raw Ds TOF header
 * 
 *   Stores infromation about Ds TOF hit PMTs
 */

#ifndef HPTPC_RAW_DSTOF_HEADER_H
#define HPTPC_RAW_DSTOF_HEADER_H

#include "TObject.h" 


class RawDsTofHeader : public TObject
{
 public:

  RawDsTofHeader();
  ~RawDsTofHeader(); 
  
  Int_t run;
  Int_t tdc; 
  Int_t channel; 
  Int_t ticks; 
  Int_t clockCounter; 
  Int_t beamSpill; 
  UInt_t unixTime; 
  Double_t triggerTimeNs; 
  Double_t fakeTimeNs; 

  ClassDef(RawDsTofHeader, 1); 

}; 

      

#endif
