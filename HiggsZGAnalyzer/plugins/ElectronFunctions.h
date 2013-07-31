#ifndef ELECTRON_FUNCTIONS
#define ELECTRON_FUNCTIONS

#include <utility> // to pick up std::pair

#include "TLorentzVector.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TFile.h"


// definitions
namespace CorrectionType { enum {LOWER, CENTRAL, UPPER}; } 
//-
#define SCALE_2011_EB             0.996
#define SCALE_UNCERT_2011_EB      0.000 
#define SCALE_2011_EE             0.997
#define SCALE_UNCERT_2011_EE      0.000
#define SIGMA_2011_EB             0.011  
#define SIGMA_UNCERT_2011_EB      0.011
#define SIGMA_2011_EE             0.023
#define SIGMA_UNCERT_2011_EE      0.023
//-
#define SCALE_2012_EB             1.000
#define SCALE_UNCERT_2012_EB      0.000
#define SCALE_2012_EE             1.000   
#define SCALE_UNCERT_2012_EE      0.000
#define SIGMA_2012_EB             0.016
#define SIGMA_UNCERT_2012_EB      0.016
#define SIGMA_2012_EE             0.028  
#define SIGMA_UNCERT_2012_EE      0.028 


// prototypes
void initializeEfficiencyWeights(const char * fname2011, 
				 const char * fname2012);
float getEfficiencyWeight( TLorentzVector * v, int type, int era );
std::pair<float,float> getScaleSigmaPair(int type, int era, int isEE);
float getCorrectedPt( TLorentzVector * v, int type, int era , TRandom3 * r); 

#endif
