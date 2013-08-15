#ifndef EvtProb_VAR
#define EvtProb_VAR

#include <TLorentzVector.h>
#include "TH2F.h"
#include "TH1F.h"

#define EBEAM 4000.00
#define fbGeV2 0.389379E12
#define smallnumber 1e-15
#define sixteen_2Pi_to_8 3.88650230418250561e+07
#define   eight_2Pi_to_5 7.83410393050320417e+04
#define    four_2Pi_to_2 39.478417604357432
class TVar{
public:
  enum VerbosityLevel {
    ERROR = 0,
    INFO = 1,
    DEBUG = 2
  };
  enum MatrixElement{
    MCFM = 0,
  };
  enum Production{
    GG = 0,
    QQB = 1,
    INDEPENDENT=2
  };

  enum Process{
    //ZZ_2e2m  =0, // eemm
    //HZZ_4l   =1, // 0+
    qqb_zgam   =0, // background
    gg_hzgam   =1, // Higgs to Z+Gamma
    Null
  };
  enum LeptonFlavor{
    Flavor_Dummy = 0, // legacy code runs on 1/2/3
    //Flavor_4e    = 1,
    //Flavor_4mu   = 2,
    //Flavor_2e2mu = 3    
    Flavor_2e    = 1,
    Flavor_2mu   = 2
  };

  //---------------------------------
  // Function
  //---------------------------------
  static TString ProcessName(int temp){ 
    //if(temp==TVar::ZZ_2e2m   ) 
    //  return TString("ZZ_2e2m");
    //else if(temp==TVar::HZZ_4l   ) 
    //  return TString("HZZ_4l");
    if(temp==TVar::qqb_zgam   )
      return TString("qqb_zgam");
    else if(temp==TVar::gg_hzgam   )
      return TString("gg_hzgam");
    else 
      return TString("UnKnown");
  };
  ClassDef(TVar,0)
};

struct branch_particle {
  int   PdgCode   ;
  int   Charge    ;
  double Px       ;
  double Py       ;
  double Pz       ;
  double E        ;
  double Eta      ;
  double Phi      ;

};
static const TString branch_format_particle =
 "PdgCode/I:"
 "Charge/I:"
 "Px/D:"
 "Py/D:"
 "Pz/D:"
 "E/D:"
 "Eta/D:"
 "Phi/D";

// in development
struct hzgamma_event_type{
  int PdgCode[3];
  TLorentzVector p[3];
  double Xsec   [10];
  double XsecErr[10];  
};

struct mcfm_event_type{
  int PdgCode[5];
  TLorentzVector p[5];
  double pswt;
};
struct event_type{
  TLorentzVector p1,p2,ep,em,nu,nb;
  double PSWeight;
};
struct anomcoup{
	   double delg1_z, delg1_g, lambda_g, lambda_z, delk_g, delk_z_,tevscale;
};

#endif
