#ifndef _MetCorrector_H
#define _MetCorrector_H

#include <map>
#include <string>
#include <iostream>
#include <cmath>
#include "TObject.h"
using namespace std;

class MetCorrector : public TObject {

  public:

    //constructor
    MetCorrector();

    //destructor
    ~MetCorrector();

    //default correction function for pfMET
    float metCorrect(float pfMet, int nVtx, string metType, bool isRealData );

    //correction function with user defined parameters for pfMET
    float metCorrect(float pfMet, int nVtx, float a0, float aPU, float sig0, float sigPU);

    ClassDef(MetCorrector,1);

  private:

    // default correction parameters for Qt = 25;
    static const float corStdMET_MC[4];
    static const float corZ_MC[4];
    static const float corLep_MC[4];
    static const float corFastMET_MC[4];
    static const float corRedMET1_MC[4];
    static const float corRedMET2_MC[4];
    static const float corTrack0_MC[4];
    static const float corTrack1_MC[4];
    static const float corTrack2_MC[4];
    static const float corTrackAll_MC[4];

    static const float corStdMET_Data[4];
    static const float corZ_Data[4];
    static const float corLep_Data[4];
    static const float corFastMET_Data[4];
    static const float corRedMET1_Data[4];
    static const float corRedMET2_Data[4];
    static const float corTrack0_Data[4];
    static const float corTrack1_Data[4];
    static const float corTrack2_Data[4];
    static const float corTrackAll_Data[4];

    static const float wrong[4];
};

#endif

#if !defined(__CINT__)
ClassImp(MetCorrector);
#endif
