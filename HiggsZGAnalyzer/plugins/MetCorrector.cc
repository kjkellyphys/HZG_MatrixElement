#include "MetCorrector.h" 

MetCorrector::MetCorrector(){

}


MetCorrector::~MetCorrector(){

}
const float MetCorrector::corStdMET_MC[4]   = {4.6957, 0.9424, 3.5843, 2.4206}; // MC
const float MetCorrector::corZ_MC[4]        = {4.871, 0.5698, 4.879, 3.051}; // MC
const float MetCorrector::corLep_MC[4]        = {1.3125, -0.0624, 2.9915, 3.2872}; // MC
const float MetCorrector::corFastMET_MC[4]  = {4.7373, 1.0427, 3.6755, 2.6469}; // MC
const float MetCorrector::corRedMET1_MC[4]  = {3.2368, 0.1747, 2.6703, 1.2289};
const float MetCorrector::corRedMET2_MC[4]  = {2.3425, 0.1479, 2.8386, 1.3340};
const float MetCorrector::corTrack0_MC[4]        = {4.8488, -0.0056, 4.5470, 0.000}; // MC
const float MetCorrector::corTrack1_MC[4]        = {7.5894, 0.2655, 6.7454, 1.6909}; // MC
const float MetCorrector::corTrack2_MC[4]        = {9.8626, 0.2683, 8.3771, 1.8784}; // MC
const float MetCorrector::corTrackAll_MC[4]        = {5.2967, 0.0182, 5.1775, 0.7371}; // MC

const float MetCorrector::corStdMET_Data[4]   = {6.4929, 0.7015, 4.6971, 2.2131}; //  DATA  
const float MetCorrector::corZ_Data[4]        = {6.585, 0.4831, 6.523, 2.74};  //  DATA  
const float MetCorrector::corLep_Data[4]        = {1.2000, -0.0417, 4.2684, 3.2716};  //  DATA  
const float MetCorrector::corFastMET_Data[4]  = {6.5594,  0.7747, 4.7820, 2.4107}; //  DATA
const float MetCorrector::corRedMET1_Data[4]  = {3.5726, 0.1124, 3.3772, 1.0158}; 
const float MetCorrector::corRedMET2_Data[4]  = {2.3737, 0.1157, 3.6279, 1.1336};
const float MetCorrector::corTrack0_Data[4]        = {4.6898, -0.0044, 4.7227, 0.000}; // Data
const float MetCorrector::corTrack1_Data[4]        = {8.3563, 0.2634, 7.2566, 1.6303}; // Data
const float MetCorrector::corTrack2_Data[4]        = {10.7611, 0.2734, 9.0082, 1.9471}; // Data
const float MetCorrector::corTrackAll_Data[4]        = {5.3020, 0.0174, 5.5736, 0.7667}; // Data

const float MetCorrector::wrong[4] = {-999,-999,-999,-999};

float MetCorrector::metCorrect(float pfMet, int nVtx, string metType, bool isRealData){

  const float* values;
  if (metType == "pf"){
    if (isRealData){ values = corStdMET_Data;}
    else { values = corStdMET_MC;}
  }else if (metType == "projZ"){
    if (isRealData){ values = corZ_Data;}
    else { values = corZ_MC;}
  }else if (metType == "projLep"){
    if (isRealData){ values = corLep_Data;}
    else { values = corLep_MC;}
  }else if (metType == "fast"){
    if (isRealData){ values = corFastMET_Data;}
    else { values = corFastMET_MC;}
  }else if (metType == "red1"){
    if (isRealData){ values = corRedMET1_Data;}
    else { values = corRedMET1_MC;}
  }else if (metType == "red2"){
    if (isRealData){ values = corRedMET2_Data;}
    else { values = corRedMET2_MC;}
  }else if (metType == "track0"){
    if (isRealData){ values = corTrack0_Data;}
    else { values = corTrack0_MC;}
  }else if (metType == "track1"){
    if (isRealData){ values = corTrack1_Data;}
    else { values = corTrack1_MC;}
  }else if (metType == "track2"){
    if (isRealData){ values = corTrack2_Data;}
    else { values = corTrack2_MC;}
  }else if (metType == "trackAll"){
    if (isRealData){ values = corTrackAll_Data;}
    else { values = corTrackAll_MC;}
  }else {
    values = wrong;
  }

  float corMET;
  float sigT;
  sigT     = sqrt(pow(values[2],2)+(nVtx-1)*pow(values[3],2));  
  corMET   = (values[2]/sigT)*(pfMet-values[0]-(nVtx-1)*values[1])+values[0];

  return corMET;
}



float MetCorrector::metCorrect(float pfMet, int nVtx, float a0, float aPU, float sig0, float sigPU){

  float corMET;
  float sigT;
  sigT     = sqrt(pow(sig0,2)+(nVtx-1)*pow(sigPU,2));  
  corMET   = (sig0/sigT)*(pfMet-a0-(nVtx-1)*aPU)+a0;

  return corMET;
}


