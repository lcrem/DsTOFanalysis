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
  
  Short_t run;
  Short_t tdc; 
  Short_t channel; 
  Int_t ticks; 
  Int_t count25k; 
  UInt_t unixTime; 
  Double_t fakeTimeNs; 

  ClassDef(RawDsTofHeader, 5); 

}; 

      

#endif
