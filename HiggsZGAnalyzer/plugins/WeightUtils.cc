#include "WeightUtils.h"

WeightUtils::WeightUtils(string sampleName, string dataPeriod, string abcd, string selection, bool isRealData)
{
  _sampleName          = sampleName;
  _dataPeriod          = dataPeriod;
  _selection           = selection;
  _isRealData          = isRealData;
  _abcd                = abcd;

  Initialize();

  _inFileS6to2011      = new TFile("otherHistos/S6to2011PU.root", "OPEN");
  _inFileS6to2012      = new TFile("otherHistos/S6to2012PU.root", "OPEN");
  _inFileS7to2011      = new TFile("otherHistos/S7to2011PU.root", "OPEN");
  _inFileS7to2012      = new TFile("otherHistos/S7to2012PU.root", "OPEN");
  _inFileS6to2011obs   = new TFile("otherHistos/S6to2011PUobsPOTER.root", "OPEN");
  _inFileS7to2012obs   = new TFile("otherHistos/S7to2012PUobsPOTER.root", "OPEN");
  _inFileS7to2012obsV2 = new TFile("otherHistos/S7to2012PUobsPOTERv2.root", "OPEN");
  _inFileS10to2012ABTrue = new TFile("otherHistos/S10to2012ABTrue.root", "OPEN");
  _inFileS10to2012CDTrue = new TFile("otherHistos/S10to2012CDTrue.root", "OPEN");
  _inFileS10to2012ABCDTrue = new TFile("otherHistos/S10to2012ABCDTrue.root", "OPEN");
  _inFileS10to2012ABCDTrue_73500 = new TFile("otherHistos/S10to2012ABCDTrue_73500.root", "OPEN");
  _inFileSomeTuneto2012ABCDTrue = new TFile("otherHistos/SomeTuneTo2012ABCD.root", "OPEN");

  _MuTrig2011          = new TFile("otherHistos/Eff_HLT_Mu17_Mu8_2011_TPfit_v2.root", "OPEN");
  _kFactors            = new TFile("otherHistos/KFactors_AllScales.root", "OPEN");
  _PhoMedWP2012        = new TFile("otherHistos/PhotonSF_MediumWP.root", "OPEN");
  _PhoMedWPveto2012    = new TFile("otherHistos/PhotonSF_MediumWP_includingVetoCut_Full2012.root", "OPEN");
  _EleMoriondWP2012    = new TFile("otherHistos/efficiency_results_EleHZGammaMoriond2013WPMixed_Moriond2013.root", "OPEN");


  // PU weights
  h1_S6to2011          = (TH1F*)_inFileS6to2011->Get("pileupWeights");
  h1_S6to2012          = (TH1F*)_inFileS6to2012->Get("pileupWeights");
  h1_S7to2011          = (TH1F*)_inFileS7to2011->Get("pileupWeights");
  h1_S7to2012          = (TH1F*)_inFileS7to2012->Get("pileupWeights");
  h1_S7to2012obs       = (TH1F*)_inFileS7to2012obs->Get("pileupWeights");
  h1_S7to2012obsV2     = (TH1F*)_inFileS7to2012obsV2->Get("pileupWeights");
  h1_S6to2011obs       = (TH1F*)_inFileS6to2011obs->Get("pileupWeights");
  h1_S10to2012ABTrue     = (TH1F*)_inFileS10to2012ABTrue->Get("pileupWeights");
  h1_S10to2012CDTrue     = (TH1F*)_inFileS10to2012CDTrue->Get("pileupWeights");
  h1_S10to2012ABCDTrue     = (TH1F*)_inFileS10to2012ABCDTrue->Get("pileupWeights");
  h1_S10to2012ABCDTrue_73500     = (TH1F*)_inFileS10to2012ABCDTrue_73500->Get("pileupWeights");
  h1_SomeTuneto2012ABCDTrue     = (TH1F*)_inFileSomeTuneto2012ABCDTrue->Get("pileupWeights");

  // higgs pt weights
  kfact120_0           = (TH1D*)_kFactors->Get("kfact120_0");
  kfact125_0           = (TH1D*)_kFactors->Get("kfact125_0");
  kfact130_0           = (TH1D*)_kFactors->Get("kfact130_0");
  kfact135_0           = (TH1D*)_kFactors->Get("kfact135_0");
  kfact140_0           = (TH1D*)_kFactors->Get("kfact140_0");
  kfact145_0           = (TH1D*)_kFactors->Get("kfact145_0");
  kfact150_0           = (TH1D*)_kFactors->Get("kfact150_0");

  //trigger histo
  _HLTMu17Mu8_2011     = (TH2D*)_MuTrig2011->Get("RunAB_EtaEta_DATAoverMC");

  //photon histo
  h1_PhoMedWP2012 = (TH2F*)_PhoMedWP2012->Get("photonSF_MediumWP");
  h1_PhoMedWPveto2012 = (TH2F*)_PhoMedWPveto2012->Get("photonSF_MediumWP");

  //electron histo
  h1_EleMoriondWP2012 = (TH2F*)_EleMoriondWP2012->Get("heff_electron_selection");
  
}

  void WeightUtils::Initialize()
  {
  _puWeight          = 1.;
  _zzWeight          = 1.;
  _glugluWeight      = 1.;
  _vbfWeight         = 1.;
  _recoWeight        = 1.;
  _triggerWeight     = 1.;
  _gammaPtPuWeight   = 1.;
  _recoilLongWeight  = 1.;
  _recoilTransWeight = 1.;



}

void WeightUtils::SetDataBit(bool isRealData)
{
  _isRealData = isRealData;
}

void WeightUtils::SetDataPeriod(string dataPeriod)
{
  _dataPeriod = dataPeriod;
}

void WeightUtils::SetSampleName(string sampleName)
{
  _sampleName = sampleName;
}

void WeightUtils::SetSelection(string selection)
{
  _selection = selection;
}


float WeightUtils::PUWeight(float nPU)
{
  //cout<<nPU<<endl;
  if (!_isRealData){
    if (nPU < 100 && _dataPeriod == "2011" ){
      _puWeight = h1_S6to2011obs->GetBinContent(h1_S6to2011obs->FindBin(nPU)); 
    } else if (nPU < 100 && nPU > 27 && _dataPeriod == "2012" && _sampleName.find("pythia8")!= string::npos){
      _puWeight = h1_SomeTuneto2012ABCDTrue->GetBinContent(h1_SomeTuneto2012ABCDTrue->FindBin(nPU)); 
    } else if (nPU < 100 && _dataPeriod == "2012" && _sampleName.find("pythia8") == string::npos){
      if (_abcd == "AB") _puWeight = h1_S10to2012ABTrue->GetBinContent(h1_S10to2012ABTrue->FindBin(nPU)); 
      else if (_abcd == "CD") _puWeight = h1_S10to2012CDTrue->GetBinContent(h1_S10to2012CDTrue->FindBin(nPU)); 
      else if (_abcd == "ABCD") _puWeight = h1_S10to2012ABCDTrue->GetBinContent(h1_S10to2012ABCDTrue->FindBin(nPU)); 
      //_puWeight = h1_S10to2012ABCDTrue_73500->GetBinContent(h1_S10to2012ABCDTrue_73500->FindBin(nPU));
    } else{
      _puWeight = 1; 
    }
  } else {
    _puWeight = 1;
  }
  return _puWeight;
}


float WeightUtils::HqtWeight(TLorentzVector l1)
{
  // for now, only weight gg higgs for 7TeV

  float hqtWeight = 1.0;
  if (_isRealData) return hqtWeight;
  if (_sampleName.find("vbf") != string::npos) return hqtWeight;
  if (_dataPeriod.find("2012") != string::npos){
    return hqtWeight;
  }else{

    float hqt = l1.Pt();
    if (hqt>250) hqt=250.0;
    if (_sampleName.find("120") != string::npos){
      hqtWeight = kfact120_0->GetBinContent(kfact120_0->FindBin(hqt));
    }else if (_sampleName.find("125") != string::npos){
      hqtWeight = kfact125_0->GetBinContent(kfact125_0->FindBin(hqt));
    }else if (_sampleName.find("130") != string::npos){
      hqtWeight = kfact130_0->GetBinContent(kfact130_0->FindBin(hqt));
    }else if (_sampleName.find("135") != string::npos){
      hqtWeight = kfact135_0->GetBinContent(kfact135_0->FindBin(hqt));
    }else if (_sampleName.find("140") != string::npos){
      hqtWeight = kfact140_0->GetBinContent(kfact140_0->FindBin(hqt));
    }else if (_sampleName.find("145") != string::npos){
      hqtWeight = kfact145_0->GetBinContent(kfact145_0->FindBin(hqt));
    }else if (_sampleName.find("150") != string::npos){
      hqtWeight = kfact150_0->GetBinContent(kfact150_0->FindBin(hqt));
    }
    
  }
  return hqtWeight;



}


float WeightUtils::MuonSelectionWeight(TLorentzVector l1) 
{

  float _muonID2012[15][16] = {
    { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, },
    { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, },
    { 0.8738,0.9118,0.9528,0.9697,0.9646,0.9798,0.9886,0.9285,0.9858,0.9877,0.955,0.9946,0.9875,0.9762,0.8779,0.8578 }, 
    { 0.9034,0.9673,0.9661,0.9883,0.9541,1.004,0.9778,0.9566,0.974,0.9839,0.9999,0.9774,0.9865,0.949,0.968,0.9126 }, 
    { 0.9548,0.9622,0.9481,0.9975,0.9809,0.9961,0.9908,0.9857,0.9773,0.9913,0.9893,0.9773,0.9882,0.9727,0.9843,0.954 }, 
    { 0.9443,0.9782,0.9831,0.9894,0.9717,0.9931,0.9874,0.9677,0.9806,0.9914,0.9818,0.9817,0.9908,0.9895,0.9927,0.9769 }, 
    { 0.9662,0.9728,0.9841,0.9826,0.9807,0.9952,0.9951,0.9784,0.9841,0.9879,0.9829,0.9827,0.9891,0.997,0.9932,0.9845 }, 
    { 0.9648,0.9831,0.9858,0.9916,0.9844,0.991,0.994,0.9872,0.9841,0.9887,0.9872,0.9841,0.9972,0.9945,0.9967,0.9835 }, 
    { 0.9735,0.9832,0.9758,0.9884,0.9805,0.9917,0.9952,0.9839,0.9832,0.9951,0.9869,0.9818,0.993,0.9969,0.9912,0.9926 }, 
    { 0.9782,0.9896,0.9904,0.9952,0.9802,0.9929,0.9935,0.9821,0.9801,0.9932,0.9905,0.9862,0.9941,1,0.9962,0.9967 }, 
    { 0.9814,0.9859,0.9797,0.994,0.9866,0.9913,0.9945,0.9857,0.9851,0.9931,0.9942,0.9887,0.997,0.9962,0.9978,1.002 }, 
    { 0.9776,0.9911,0.9844,0.9954,0.9844,0.9939,0.9969,0.9816,0.9899,0.9931,0.9917,0.9892,0.9962,0.999,0.9952,1.002 }, 
    { 0.981,0.9871,0.9799,0.9943,0.9845,0.9897,0.9917,0.9799,0.9869,0.9929,0.9879,0.9837,0.9949,1.008,0.9919,1.002 }, 
    { 0.9789,0.9813,0.9828,0.993,0.9805,0.9896,0.9932,0.9854,0.9846,0.9873,0.9865,0.9802,0.9906,0.9979,0.9951,0.993 }, 
    { 0.9635,0.9833,0.9733,0.9796,0.9771,0.988,0.9895,0.9832,0.9788,0.9961,0.9839,0.9775,0.9869,1.007,0.9825,0.979 }
  };


  float _muonISO2012[15][16] = {
    { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, },
    { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, },
    { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, },
    { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, },
    { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, },
    { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, },
    { 1.027,1.009,0.9978,0.9845,0.9925,0.988,0.9939,0.9776,0.9899,0.9862,0.9877,0.9882,0.9895,0.9888,0.995,1.03 }, 
    { 1.019,1.001,1.002,0.997,1.003,0.9895,0.9915,0.9941,0.9929,0.9902,0.9898,0.9928,0.999,1.001,1.002,1.026 }, 
    { 1.015,0.999,1.003,0.9956,0.9965,0.9978,0.9995,0.9952,0.9965,0.9962,0.9993,0.9965,0.9957,0.9963,1.002,1.014 }, 
    { 1.007,0.9997,0.9996,0.9994,0.9974,0.9959,0.9954,0.9938,0.995,0.9958,1,0.9958,0.9968,0.9989,1.003,1.01 }, 
    { 1.007,1.001,0.9975,0.999,0.9978,0.9984,0.9959,0.9965,0.995,0.9971,0.9969,0.9986,0.9992,0.9984,1,1.007 }, 
    { 1.002,0.9991,0.9996,0.9986,0.9974,0.9982,0.9974,0.9963,0.9957,0.9971,0.9985,0.9992,0.9979,0.9988,1.001,1.002 }, 
    { 1.001,0.9979,0.9973,0.9983,0.9986,0.9977,0.9969,0.9947,0.9971,0.9985,0.9981,0.997,0.9986,0.9986,1,1.003 }, 
    { 1.006,1.001,0.9983,0.9988,0.9978,0.9979,0.9975,0.9944,0.9958,1,0.9973,0.9983,0.9996,0.9969,0.9999,1.003 }, 
    { 0.9967,0.9967,0.9992,0.9993,0.9988,0.9923,0.9953,0.9951,0.9988,0.9975,0.9979,0.9977,0.9985,0.9926,0.9991,0.9936 }
  };

  float _muonID2011[15][16] = {
    { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, }, 
    { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, }, 
    { 0.983,0.9641,0.9367,1.002,1,0.9625,0.988,1.003,1.003,0.9613,0.986,1.004,0.9877,0.9299,0.9619,0.9726, }, 
    { 0.9735,0.9763,0.9687,0.959,0.9951,1.009,0.9922,0.9866,0.9912,1.002,0.9982,1.004,0.993,1.012,0.9762,0.9931, }, 
    { 0.9847,0.9841,0.9668,0.9789,0.9928,0.9972,0.9961,0.9951,1.002,0.9959,0.9937,0.9897,1,0.9907,0.9767,1.001, }, 
    { 0.9929,0.9886,0.9536,0.9812,0.9911,0.9962,0.9997,0.992,0.9949,0.9982,0.9947,1.002,0.99,0.9804,0.985,0.9837, }, 
    { 0.9802,0.9811,0.9569,0.9873,0.9938,0.9996,0.9976,0.9867,0.9973,0.9988,0.9949,0.9949,0.9913,0.9899,0.9904,0.988, }, 
    { 0.9798,0.9878,0.9612,0.9831,0.9879,1,0.9987,0.9932,0.993,0.9966,0.9962,0.9951,0.9931,0.9845,0.9895,0.993, }, 
    { 0.9743,0.9807,0.9563,0.984,0.9919,0.9965,0.9978,0.9935,0.9964,0.9955,0.9975,0.9971,0.9899,0.9831,0.9934,0.991, }, 
    { 0.9664,0.985,0.9554,0.9862,0.9896,0.9987,0.9974,0.9895,0.9898,0.9975,0.9946,0.994,0.9905,0.9869,0.9923,0.9857, }, 
    { 0.9838,0.9729,0.958,0.9855,0.9892,0.9939,0.9968,0.9959,0.9884,0.9959,0.9906,0.9941,0.9912,0.9803,0.9871,0.9815, }, 
    { 0.9715,0.9802,0.9541,0.9786,0.9854,0.9943,0.9983,0.9778,0.9948,0.9974,0.9945,0.9946,0.9819,0.9722,0.9944,0.9822, }, 
    { 0.9275,0.9749,0.9554,0.9851,0.9892,0.9948,0.9926,0.9895,0.9814,1.002,0.9934,0.9961,0.9921,0.9732,0.989,0.9948, }, 
    { 1.038,0.9818,0.9379,0.9845,0.9874,1.002,0.9999,0.9944,1.007,1,0.9936,0.9901,0.9888,0.9768,0.9943,0.9304, }, 
    { 0.9111,1.024,1,0.9845,0.9993,1.004,1.005,0.9751,0.9803,0.9967,1.001,0.9895,0.9834,0.9632,0.9917,1.008, }, 
  };

  float _muonISO2011[15][16] = {
    { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, }, 
    { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, }, 
    { 0.8949,0.9304,0.76,0.8702,0.9242,0.8778,0.9139,0.8324,0.7679,0.8131,0.8504,0.9509,0.5491,0.8272,0.8651,0.946, }, 
    { 0.9837,0.9668,0.9725,0.9073,0.9462,0.9548,0.9559,0.9417,0.9068,0.9367,0.9449,0.9355,0.9436,0.9035,1.008,0.9863, }, 
    { 0.9861,0.9964,0.9772,0.9931,0.9707,0.9604,0.976,0.9711,0.9798,0.9896,0.9711,0.9915,0.9812,0.9759,0.988,0.9955, }, 
    { 0.9929,1,1.004,0.989,0.9855,0.981,0.9897,0.9875,0.9924,0.9936,1.003,0.9975,0.9927,0.9935,1.005,1.004, }, 
    { 0.9977,1.001,0.9996,1.001,1.009,0.998,0.9928,0.9939,0.9938,0.9975,0.9966,1.003,0.9973,0.9962,0.998,1.004, }, 
    { 1,0.9969,1.003,0.9984,1.002,0.9992,0.9961,0.9912,0.9949,0.9972,0.9976,0.9957,0.9987,0.9966,1,1.002, }, 
    { 1.002,1.002,0.9977,0.999,1.002,1,0.9973,0.9962,0.9975,0.9957,0.9983,0.9976,1,1.001,1.001,1.002, }, 
    { 1,0.9995,1,0.9999,0.9993,1,0.9992,0.9961,0.9975,0.9977,0.9982,0.9991,0.9991,0.9989,0.9998,0.9975, }, 
    { 0.992,1.001,0.9994,1,0.999,0.9983,0.9978,0.9966,0.9987,0.9986,0.9991,0.9994,1,0.9983,0.9996,1.006, }, 
    { 1.001,0.9945,1,0.9998,0.998,0.9995,0.9987,0.9977,0.9993,0.9987,0.9979,0.9992,0.9996,0.9995,0.9995,0.9958, }, 
    { 1.008,0.9992,1.001,0.998,0.9981,0.9973,0.9999,0.9979,1.001,0.9974,0.9979,0.9981,1,1.002,1.003,0.998, }, 
    { 1,1.003,0.9954,0.9975,0.9977,0.9999,0.9994,0.9946,1.002,0.9969,0.9993,0.9996,1.002,1.001,0.9924,0.9825, }, 
    { 1,1,1,0.9974,0.9991,1.002,0.9962,0.9981,1,0.9984,1.002,1.001,0.9997,1,0.9973,1.011, }, 
  };
  float muSF = 1.0;
  float binningPt[]  = {0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,9999999};
  float binningEta[] = {-2.4,-2.1,-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,0,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4};
  int ptBin = 0;
  int etaBin = 0;

  if (_isRealData) return muSF;

  for (int i = 0; i<16; i++){
    if (l1.Eta() > binningEta[i] && l1.Eta() < binningEta[i+1]){
      etaBin = i;
      break;
    }
  }
  for (int i = 0; i<16; i++){
    if (l1.Pt() > binningPt[i] && l1.Pt() < binningPt[i+1]){
      ptBin = i;
      break;
    }
  }

  if (_dataPeriod.find("2011") != string::npos){
    muSF = _muonID2011[ptBin][etaBin]*_muonISO2011[ptBin][etaBin];
  }else{
    muSF = _muonID2012[ptBin][etaBin]*_muonISO2012[ptBin][etaBin];
  }
  return muSF;

}

float WeightUtils::GammaSelectionWeight(TLorentzVector l1, float SCEta) 
{
  float _gammaMedID2011[4][5] = {
    //15<pt<20   20<pt<30    30<pt<40   40<pt<50   50<pt
    {1.0022,    0.9854,     0.9834,    0.9853,    0.9847},                        // 0.0<eta<0.8
    {1.0022,    0.9876,     0.9932,    0.9914,    0.9941},                        // 0.8<eta<1.4442
    {1.0189,    1.0155,     1.0048,    1.0013,    1.0033},                        // 1.566<eta<2.0
    {1.0189,    1.0056,     1.0108,    1.0057,    1.0106}                         // 2.0<eta<2.5
  };

  float _gammaMedVeto2011Barrel = 0.9971;
  float _gammaMedVeto2011Endcap = 0.9983;

  float phoSF = 1.0;
  if (_isRealData) return phoSF;

  int ptBin = 0;
  int etaBin = 0;
  float binningPt[]  = {15., 20., 30., 40., 50., 99999.};
  float binningEta[] = {0., 0.8, 1.4442, 2.0, 2.5};
  float weight = 1.;


  for (int i = 0; i < 4; ++i) {
    //if (fabs(l1.Eta()) > binningEta[i] && fabs(l1.Eta()) <= binningEta[i+1]) {
    if (fabs(SCEta) > binningEta[i] && fabs(SCEta) <= binningEta[i+1]) {
      etaBin = i;
      break;
    }
  }
  for (int i = 0; i < 5; ++i) {
    if (l1.Pt() > binningPt[i] && l1.Pt() <= binningPt[i+1]) {
      ptBin = i;
      break;
    }
  }

  if (_dataPeriod.find("2011") != string::npos){
    if (etaBin <= 1){
      phoSF = _gammaMedID2011[etaBin][ptBin]*_gammaMedVeto2011Barrel;
    }else{
      phoSF = _gammaMedID2011[etaBin][ptBin]*_gammaMedVeto2011Endcap;
    }
  }else{
    phoSF = h1_PhoMedWPveto2012->GetBinContent(h1_PhoMedWPveto2012->GetXaxis()->FindFixBin(l1.Pt()), h1_PhoMedWPveto2012->GetYaxis()->FindFixBin(fabs(SCEta)));
  }
  return phoSF;

}

float WeightUtils::ElectronSelectionWeight(TLorentzVector l1){
  float eleSF;
  int xBin;
  int yBin;

  xBin = h1_EleMoriondWP2012->GetXaxis()->FindFixBin(l1.Pt());
  if (xBin>7) xBin=7;
  yBin = h1_EleMoriondWP2012->GetYaxis()->FindFixBin(fabs(l1.Eta()));
  if (yBin>5) yBin=5;

  if (_dataPeriod.find("2012") != string::npos){
    eleSF = h1_EleMoriondWP2012->GetBinContent(xBin,yBin);
  }else{
    eleSF = 1;
  }
  return eleSF;
}


float WeightUtils::MuonTriggerWeight(TLorentzVector l1, TLorentzVector l2) 
{

  float _DataEff_HLTMu17Mu8_8Leg2012[9][3] = {  
 //|eta|<0.9 , 0.9<|eta|<1.2 , 1.2<|eta|<2.4
 {0.967172   , 0.930476      , 0.916004}, // 10<pt<20
 {0.967066   , 0.93713       , 0.920873}, // 20<pt<30
 {0.965728   , 0.935549      , 0.916849}, // 30<pt<35
 {0.965991   , 0.932407      , 0.915271}, // 35<pt<40
 {0.965568   , 0.935851      , 0.918777}, // 40<pt<50
 {0.964775   , 0.937591      , 0.917744}, // 50<pt<60
 {0.96494    , 0.933094      , 0.918446}, // 60<pt<90
 {0.960397   , 0.939106      , 0.909972}, // 90<pt<140
 {0.961868   , 0.92921       , 0.937057}  // 140<pt<500
  };

  float _MCEff_HLTMu17Mu8_8Leg2012[9][3] = {  
 //|eta|<0.9 , 0.9<|eta|<1.2 , 1.2<|eta|<2.4
    {0.969444 , 0.925885 , 0.921075} , // 10<pt<20
    {0.976136 , 0.945697 , 0.927715} , // 20<pt<30
    {0.976453 , 0.948453 , 0.926418} , // 30<pt<35
    {0.975895 , 0.944976 , 0.925758} , // 35<pt<40
    {0.976162 , 0.946528 , 0.928904} , // 40<pt<50
    {0.975507 , 0.950925 , 0.931956} , // 50<pt<60
    {0.976274 , 0.951396 , 0.926831} , // 60<pt<90
    {0.972801 , 0.951352 , 0.932385} , // 90<pt<140
    {0.971679 , 0.973072 , 0.939368}   // 140<pt<500
  };

  float _DataEff_HLTMu17Mu8_17Leg2012[9][3] = {  
 //|eta|<0.9 , 0.9<|eta|<1.2 , 1.2<|eta|<2.4
    {0.609746 , 0.496018 , 0.428991} , // 10<pt<20
    {0.964723 , 0.923791 , 0.892096} , // 20<pt<30
    {0.964065 , 0.924091 , 0.896823} , // 30<pt<35
    {0.964584 , 0.923641 , 0.898797} , // 35<pt<40
    {0.964363 , 0.928434 , 0.90573}  , // 40<pt<50
    {0.963617 , 0.930997 , 0.907169} , // 50<pt<60
    {0.963878 , 0.925745 , 0.908756} , // 60<pt<90
    {0.960051 , 0.935225 , 0.901006} , // 90<pt<140
    {0.959271 , 0.92921  , 0.937057}   // 140<pt<500
  };

  float _MCEff_HLTMu17Mu8_17Leg2012[9][3] = {  
 //|eta|<0.9 , 0.9<|eta|<1.2 , 1.2<|eta|<2.4
    {0.617508 , 0.488784 , 0.428354} , // 10<pt<20
    {0.97418  , 0.935211 , 0.893312} , // 20<pt<30
    {0.975246 , 0.93891  , 0.903676} , // 30<pt<35
    {0.974711 , 0.93787  , 0.907107} , // 35<pt<40
    {0.975291 , 0.939777 , 0.915754} , // 40<pt<50
    {0.974371 , 0.94515  , 0.920956} , // 50<pt<60
    {0.975252 , 0.946933 , 0.917094} , // 60<pt<90
    {0.972801 , 0.945771 , 0.92517}  , // 90<pt<140
    {0.971679 , 0.973072 , 0.931013}   // 140<pt<500
  };


  float _HLTMu17Mu8_2012[3][10] = {  
    //10<pt<20 , 20<pt<25 , 25<pt<30 , 30<pt<35 , 35<pt<40 , 40<pt<50 , 50<pt<60 , 60<pt<90 , 90<pt<140 , 140<pt<500
    {0.991     , 0.989    , 0.989    , 0.988    , 0.989    , 0.988    , 0.987    , 0.989    , 0.990     , 0.982}, // |eta| < 0.9
    {1.00      , 0.985    , 0.986    , 0.983    , 0.984    , 0.984    , 0.986    , 0.983    , 0.982     , 0.964}, // 0.9 < |eta| < 1.2
    {1.01      , 1.00     , 0.994    , 0.991    , 0.987    , 0.986    , 0.985    , 0.984    , 0.972     , 1.01}   // |eta| > 1.2
  };


  float muTrigSF1 = 1.0;
  float muTrigSF2 = 1.0;

  float muTrigDataA8 = 1.0;
  float muTrigDataA17 = 1.0;
  float muTrigDataB8 = 1.0;
  float muTrigDataB17 = 1.0;

  float muTrigMCA8 = 1.0;
  float muTrigMCA17 = 1.0;
  float muTrigMCB8 = 1.0;
  float muTrigMCB17 = 1.0;

  if (_isRealData) return 1.0;

  // 2011 use the 2D histogram
  if (_dataPeriod.find("2011") != string::npos){
    muTrigSF1 = _HLTMu17Mu8_2011->GetBinContent(_HLTMu17Mu8_2011->FindBin(l1.Eta(),l2.Eta()));
    if (muTrigSF1 == 0.0) muTrigSF1 = 1.0;
    return muTrigSF1;
  }else{
  // 2012 use the 2D arrays

    int ptBinV1 = 0;
    int ptBinV2 = 0;
    int etaBin = 0;
    float binningPtV1[]  = {10., 20., 25., 30., 35., 40., 50., 60., 90., 140., 500.};
    float binningPtV2[]  = {10., 20., 30., 35., 40., 50., 60., 90., 140., 500.};

    if (fabs(l1.Eta()) < 0.9) {
      etaBin = 0;
    }else if (fabs(l1.Eta()) < 1.2){
      etaBin = 1;
    }else{
      etaBin = 2;
    }
    for (int i = 0; i < 10; ++i) {
      if (l1.Pt() >= binningPtV1[i] && l1.Pt() < binningPtV1[i+1]) {
        ptBinV1 = i;
        break;
      }
    }
    for (int i = 0; i < 9; ++i) {
      if (l1.Pt() >= binningPtV2[i] && l1.Pt() < binningPtV2[i+1]) {
        ptBinV2 = i;
        break;
      }
    }
    muTrigSF1   = _HLTMu17Mu8_2012[etaBin][ptBinV1];
    muTrigDataA17 = _DataEff_HLTMu17Mu8_17Leg2012[etaBin][ptBinV2];
    muTrigDataA8  = _DataEff_HLTMu17Mu8_8Leg2012[etaBin][ptBinV2];
    muTrigMCA17 = _MCEff_HLTMu17Mu8_17Leg2012[etaBin][ptBinV2];
    muTrigMCA8  = _MCEff_HLTMu17Mu8_8Leg2012[etaBin][ptBinV2];


    if (fabs(l2.Eta()) < 0.9) {
      etaBin = 0;
    }else if (fabs(l2.Eta()) < 1.2){
      etaBin = 1;
    }else{
      etaBin = 2;
    }
    for (int i = 0; i < 10; ++i) {
      if (l2.Pt() >= binningPtV1[i] && l2.Pt() < binningPtV1[i+1]) {
        ptBinV1 = i;
        break;
      }
    }
    for (int i = 0; i < 9; ++i) {
      if (l2.Pt() >= binningPtV2[i] && l2.Pt() < binningPtV2[i+1]) {
        ptBinV2 = i;
        break;
      }
    }
    muTrigSF2   = _HLTMu17Mu8_2012[etaBin][ptBinV1];
    muTrigDataB17 = _DataEff_HLTMu17Mu8_17Leg2012[etaBin][ptBinV2];
    muTrigDataB8  = _DataEff_HLTMu17Mu8_8Leg2012[etaBin][ptBinV2];
    muTrigMCB17 = _MCEff_HLTMu17Mu8_17Leg2012[etaBin][ptBinV2];
    muTrigMCB8  = _MCEff_HLTMu17Mu8_8Leg2012[etaBin][ptBinV2];

    //return muTrigSF1*muTrigSF2;
    return (muTrigDataA8*muTrigDataB17 + muTrigDataA17*muTrigDataB8 - muTrigDataA17*muTrigDataB17)/
           (muTrigMCA8*muTrigMCB17 + muTrigMCA17*muTrigMCB8 - muTrigMCA17*muTrigMCB17);
  }


}

float WeightUtils::ElectronTriggerWeight(TLorentzVector l1, TLorentzVector l2) 
{

 ///////////////////
 // Scale Factors //
 ///////////////////
 
 float  _HLTEl17El8_8Leg2012[3][6] = {  
    //10<pt<15  15<pt<20  20<pt<30  30<pt<40  40<pt<50  50<pt
    {  0.9085 , 0.9750 , 0.9869 , 0.9908 , 0.9914 , 0.9929  },   // |eta| < 1.4442, trailing  
    {  0.9493 , 1.0459 , 1.0033 , 0.9929 , 0.9940 , 0.9944  },   // 1.4442 > |eta| > 1.566 
    {  0.9010 , 0.9715 , 0.9966 , 0.9954 , 0.9977 , 0.9979  }    // |eta| > 1.566, trailing  
  };

  float _HLTEl17El8_17Leg2012[3][4] = { 
    //20<pt<30  30<pt<40  40<pt<60  50<pt
    {  0.9863 , 0.9910 , 0.9920 , 0.9933  }, // |eta| < 1.4442, leading 
    {  0.9664 , 0.9645 , 0.9752 , 0.9868  }, // 1.4442 > |eta| > 1.566
    {  0.9892 , 0.9965 , 0.9991 , 0.9998  }  // |eta| > 1.566, leading 


  };

  float _HLTEl17El8_8Leg2011[2][6] = {  
    //10<pt<15  15<pt<20  20<pt<30  30<pt<40  40<pt<50  50<pt
    {1.0200,    1.0115,   1.0016,   0.9993,   0.9983,   0.9990}, // |eta| < 1.4442, trailing 
    {0.9905,    1.0062,   1.0039,   1.0038,   1.0020,   1.0032}  // |eta| > 1.566, trailing 
  };

  float _HLTEl17El8_17Leg2011[2][4] = { 
    //20<pt<30  30<pt<40  40<pt<60  50<pt
    {0.9953,    0.9973,   0.9971,   0.9979}, // |eta| < 1.4442, leading 
    {0.9860,    1.0014,   0.9999,   1.0017}  // |eta| > 1.566, leading 
  };

 //////////////////
 // Efficiencies //
 //////////////////
 
 float  _DataEff_HLTEl17El8_8Leg2011[2][6] = {  
    //10<pt<15  15<pt<20  20<pt<30  30<pt<40  40<pt<50  50<pt
    {0.937113,    0.967319,   0.982833,   0.990721,   0.994619,   0.995135}, // |eta| < 1.4442, trailing 
    {0.950915,    0.981737,   0.988265,   0.993117,   0.995014,   0.9956}  // |eta| > 1.566, trailing 
  };

 float  _MCEff_HLTEl17El8_8Leg2011[2][6] = {  
    //10<pt<15  15<pt<20  20<pt<30  30<pt<40  40<pt<50  50<pt
    {0.918734,    0.95629,   0.981257,   0.99138,   0.996299,   0.996177}, // |eta| < 1.4442, trailing 
    {0.96,    0.97568,   0.984404,   0.989377,   0.99303,   0.992393}  // |eta| > 1.566, trailing 
  };

  float _DataEff_HLTEl17El8_17Leg2011[2][4] = { 
    //20<pt<30  30<pt<40  40<pt<60  50<pt
    {0.975521,    0.988414,   0.99346,   0.994072}, // |eta| < 1.4442, leading 
    {0.970854,    0.991086,   0.993635,   0.994258}  // |eta| > 1.566, leading 
  };

  float _MCEff_HLTEl17El8_17Leg2011[2][4] = { 
    //20<pt<30  30<pt<40  40<pt<60  50<pt
    {0.980084,    0.991064,   0.996344,   0.996177}, // |eta| < 1.4442, leading 
    {0.984609,    0.989658,   0.993723,   0.992564}  // |eta| > 1.566, leading 
  };

  
 float  _DataEff_HLTEl17El8_8Leg2012[2][6] = {  
    //10<pt<15  15<pt<20  20<pt<30  30<pt<40  40<pt<50  50<pt
    {0.81565,    0.912922,   0.953769,   0.973848,   0.982891,   0.985804}, // |eta| < 1.4442, trailing 
    {0.791932,    0.928296,   0.967244,   0.974564,   0.979996,   0.982088}  // |eta| > 1.566, trailing 
  };

 float  _MCEff_HLTEl17El8_8Leg2012[2][6] = {  
    //10<pt<15  15<pt<20  20<pt<30  30<pt<40  40<pt<50  50<pt
    {0.768699,    0.880306,   0.949366,   0.977081,   0.989122,   0.990735}, // |eta| < 1.4442, trailing 
    {0.86336,    0.941864,   0.96837,   0.979045,   0.984052,   0.985658}  // |eta| > 1.566, trailing 
  };

  float _DataEff_HLTEl17El8_17Leg2012[2][4] = { 
    //20<pt<30  30<pt<40  40<pt<60  50<pt
    {0.954609,    0.978026,   0.987355,   0.989184}, // |eta| < 1.4442, leading 
    {0.97516,    0.984818,   0.990335,   0.990963}  // |eta| > 1.566, leading 
  };

  float _MCEff_HLTEl17El8_17Leg2012[2][4] = { 
    //20<pt<30  30<pt<40  40<pt<60  50<pt
    {0.951421,    0.979938,   0.99201,   0.99347}, // |eta| < 1.4442, leading 
    {0.977215,    0.986849,   0.992508,   0.993102}  // |eta| > 1.566, leading 
  };

 
  float elTrigSF1 = 1.0;
  float elTrigSF2 = 1.0;

  float elTrigDataA8 = 1.0;
  float elTrigDataA17 = 1.0;
  float elTrigDataB8 = 1.0;
  float elTrigDataB17 = 1.0;

  float elTrigMCA8 = 1.0;
  float elTrigMCA17 = 1.0;
  float elTrigMCB8 = 1.0;
  float elTrigMCB17 = 1.0;

  if (_isRealData) return 1.0;

  int ptBin1 = 0;
  int etaBin1 = 0;
  int ptBin2 = 0;
  int etaBin2 = 0;
  float binningPt1[] = {20., 30., 40., 50., 9999.};
  float binningPt2[] = {10., 15., 20., 30., 40., 50., 9999.};
  if (fabs(l1.Eta()) < 1.4442) {
    etaBin1 = 0;
  }else if (fabs(l1.Eta()) < 1.566){
    etaBin1 = 1;
  }else{
    etaBin1 = 2;
  }

  for (int i = 0; i < 4; ++i) {
    if (l1.Pt() > binningPt1[i] && l1.Pt() <= binningPt1[i+1]) {
      ptBin1 = i;
      break;
    }
  }

  if (fabs(l2.Eta()) < 1.4442) {
    etaBin2 = 0;
  }else if (fabs(l2.Eta()) < 1.566){
    etaBin2 = 1;
  }else{
    etaBin2 = 2;
  }

  for (int i = 0; i < 6; ++i) {
    if (l2.Pt() > binningPt2[i] && l2.Pt() <= binningPt2[i+1]) {
      ptBin2 = i;
      break;
    }
  }

  if (_dataPeriod.find("2011") != string::npos){
    elTrigSF1 = _HLTEl17El8_17Leg2011[etaBin1][ptBin1];
    elTrigSF2 = _HLTEl17El8_8Leg2011[etaBin2][ptBin2];
    //elTrigDataA17 = _DataEff_HLTEl17El8_17Leg2011[etaBin1][ptBin1];
    //elTrigDataB17 = _DataEff_HLTEl17El8_17Leg2011[etaBin2][ptBin2-2];
    //elTrigDataA8 = _DataEff_HLTEl17El8_8Leg2011[etaBin1][ptBin1+2];
    //elTrigDataB8 = _DataEff_HLTEl17El8_8Leg2011[etaBin2][ptBin2];
    //elTrigMCA17 = _MCEff_HLTEl17El8_17Leg2011[etaBin1][ptBin1];
    //elTrigMCB17 = _MCEff_HLTEl17El8_17Leg2011[etaBin2][ptBin2-2];
    //elTrigMCA8 = _MCEff_HLTEl17El8_8Leg2011[etaBin1][ptBin1+2];
    //elTrigMCB8 = _MCEff_HLTEl17El8_8Leg2011[etaBin2][ptBin2];
  }else{
    elTrigSF1 = _HLTEl17El8_17Leg2012[etaBin1][ptBin1];
    elTrigSF2 = _HLTEl17El8_8Leg2012[etaBin2][ptBin2];
    //elTrigDataA17 = _DataEff_HLTEl17El8_17Leg2012[etaBin1][ptBin1];
    //elTrigDataB17 = _DataEff_HLTEl17El8_17Leg2012[etaBin2][ptBin2-2];
    //elTrigDataA8 = _DataEff_HLTEl17El8_8Leg2012[etaBin1][ptBin1+2];
    //elTrigDataB8 = _DataEff_HLTEl17El8_8Leg2012[etaBin2][ptBin2];
    //elTrigMCA17 = _MCEff_HLTEl17El8_17Leg2012[etaBin1][ptBin1];
    //elTrigMCB17 = _MCEff_HLTEl17El8_17Leg2012[etaBin2][ptBin2-2];
    //elTrigMCA8 = _MCEff_HLTEl17El8_8Leg2012[etaBin1][ptBin1+2];
    //elTrigMCB8 = _MCEff_HLTEl17El8_8Leg2012[etaBin2][ptBin2];
  }

  //return (elTrigDataA8*elTrigDataB17 + elTrigDataA17*elTrigDataB8 - elTrigDataA17*elTrigDataB17)/
  //       (elTrigMCA8*elTrigMCB17 + elTrigMCA17*elTrigMCB8 - elTrigMCA17*elTrigMCB17);
  return elTrigSF1*elTrigSF2;
}


