#ifndef ZGANGLES_HH
#define ZGANGLES_HH
#define STANDALONE

#include "TLorentzVector.h"
#ifndef STANDALONE
#include "EventData.h"
#include "GenInfoStruct.h"
#endif

#define VARLIST_ZGANGLES "costheta_lm/F:costheta_lp/F:phi/F:cosTheta/F:cosThetaG/F:ptg/F:ptz/F:ptl1/F:ptl2/F:etag/F:etaz/F:etal1/F:etal2/F:mzg/F:mz/F"

class ZGAngles { 
 public:
  float costheta_lm, costheta_lp, phi, cosTheta, cosThetaG;
  float ptg, ptz, ptl1, ptl2, etag, etaz, etal1, etal2, mzg, mz;
};



class ZGLabVectors { 
 public:
  TLorentzVector veckq;
  TLorentzVector veckqbar;

  TLorentzVector veczg;
  TLorentzVector vecz;
  TLorentzVector vecg;
  TLorentzVector veclp;
  TLorentzVector veclm;
};


void getZGAngles( ZGLabVectors &l, ZGAngles &a, bool=false );
#ifndef STANDALONE
void fillZGAngles( EventData &, ZGAngles &a );
void fillZGAngles( GenInfoStruct &geninfo, ZGAngles &a, bool=false );
#endif
#endif
