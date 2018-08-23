/*
 * RawTofMini.h
 * HPTPC mini tree class with just ustof and delayed beam signal info
 *
 */

#ifndef HPTPC_RAW_TOF_MINI
#define HPTPC_RAW_TOF_MINI

#include "TObject.h"

class RawTofMini: public TObject
{
 public:
  RawTofMini();
  ~RawTofMini();

  Int_t run;
  Bool_t inSpill;
  Double_t fakeTimeNs;
  Double_t lastDelayedBeamSpillNs;

  ClassDef(RawTofMini, 3);
  
};

#endif
