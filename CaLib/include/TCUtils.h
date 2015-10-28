// SVN Info: $Id$

/*************************************************************************
 * Author: Dominik Werthmueller
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCUtils                                                              //
//                                                                      //
// CaLib utility methods namespace                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef TCUTILS_H
#define TCUTILS_H

#include "TH2.h"
#include "TMath.h"

#include "TCReadConfig.h"

// define A2 detectors
enum EA2Detector {
    kNoDetector,            // no detector
    kCBDetector,            // Crystal Ball
    kTAPSDetector           // TAPS
};
typedef EA2Detector A2Detector_t;

// define TAPS types
enum EA2TAPSType {
    kTAPS_2007,             // TAPS 2007 with 384 BaF2
    kTAPS_2008,             // TAPS 2008 with 378 BaF2 + 24 PbWO4
    kTAPS_2009              // TAPS 2009 with 366 BaF2 + 72 PbWO4
};
typedef EA2TAPSType A2TAPSType_t;


namespace TCUtils
{
    void FindBackground(TH1* h, Double_t peak, Double_t low, Double_t high,
                        Double_t* outPar0, Double_t* outPar1);
    TH1* DeriveHistogram(TH1* inH);
    void ZeroBins(TH1* inH, Double_t th = 0);
    Double_t Pi0Func(Double_t* x, Double_t* par);
    Double_t GetHistogramMinimum(TH1* h);
    Double_t GetHistogramMinimumPosition(TH1* h);
    void FormatHistogram(TH1* h, const Char_t* ident);
    Bool_t IsCBHole(Int_t elem);
    Int_t GetTAPSRing(Int_t id, A2TAPSType_t type);
    Int_t GetVetoInFrontOfElement(Int_t id, Int_t maxTAPS);
    Double_t GetDiffPercent(Double_t oldValue, Double_t newValue);
}

#endif

