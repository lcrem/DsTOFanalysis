//////////////////////////////////////////////////////////////////////////////
/////  DsTOFConventions.h        DSTOF Convetnions                       /////
/////                                                                    /////
/////  Author: Linda Cremonesi (l.cremonesi@ucl.ac.uk)                   /////
//////////////////////////////////////////////////////////////////////////////


#ifndef DSTOFCONVENTIONS_H
#define DSTOFCONVENTIONS_H
#include <string>

const std::string tdc1[10] = {"10B", "10A", "09B", "09A", "08B", "08A", "07B", "07A", "06B", "06A" }; // channels from 1 to 10
const std::string tdc2[10] = {"05B", "05A", "04B", "04A", "03B", "03A", "02B", "02A", "01B", "01A" };

const int tdc1bar[15] = {10, 10, 9, 9, 8, 8, 7, 7, 6, 6, -1, -1, -1, -1, -1 }; // channels from 1 to 1
const int tdc2bar[15] = { 5,  5, 4, 4, 3, 3, 2, 2, 1, 1, -1, -1, -1, -1, -1 };

const int tdc1pmt[15] = { 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, -1, -1, -1, -1, -1}; // 0 is PMT A and 1 is PMT B
const int tdc2pmt[15] = { 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, -1, -1, -1, -1, -1}; 

const double clockTicksNs=0.025/1.024;

const double dstofDelayAllBars[10] = { 70.2, 61.6, 61.6, 61.6, 61.6, 61.6, 61.6, 61.6, 61.6, 61.6};

const std::string whereIsMyTOFdata="/unix/dune/hptpctof/";

#endif //DSTOFCONVENTIONS_H
