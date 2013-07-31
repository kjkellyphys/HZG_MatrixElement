#include "WeightUtils.h"

WeightUtils::WeightUtils(string sampleName, string dataPeriod, string selection, bool isRealData)
{
  _sampleName = sampleName;
  _dataPeriod = dataPeriod;
  _selection  = selection;
  _isRealData = isRealData;

  Initialize();

  _inFile         = new TFile("otherHistos/gammaCombined.root", "OPEN");

  _inFileS6to2011 = new TFile("otherHistos/S6to2011PU.root", "OPEN");
  _inFileS6to2012 = new TFile("otherHistos/S6to2012PU.root", "OPEN");
  _inFileS7to2011 = new TFile("otherHistos/S7to2011PU.root", "OPEN");
  _inFileS7to2012 = new TFile("otherHistos/S7to2012PU.root", "OPEN");

  _inFilePhoA_EE = new TFile("otherHistos/runA_photonId_EEEffsSF.root", "OPEN");
  _inFilePhoA_EB = new TFile("otherHistos/runA_photonId_EBEffsSF.root", "OPEN");
  _inFilePhoB_EE = new TFile("otherHistos/runB_photonId_EEEffsSF.root", "OPEN");
  _inFilePhoB_EB = new TFile("otherHistos/runB_photonId_EBEffsSF.root", "OPEN");

  // Photon weights
  //h1_eGammaPt     = (TH1D*)_inFile->Get(("h1_eGammaPtWeight"+_dataPeriod).c_str());
  //h1_muGammaPt    = (TH1D*)_inFile->Get(("h1_muGammaPtWeight"+_dataPeriod).c_str());
  //h1_eGammaPV     = (TH1D*)_inFile->Get(("h1_eGammaPVWeight"+_dataPeriod).c_str());
  //h1_muGammaPV    = (TH1D*)_inFile->Get(("h1_muGammaPVWeight"+_dataPeriod).c_str());
  //h1_eGammaMass   = (TH1D*)_inFile->Get(("h1_diElectronMass"+_dataPeriod).c_str());
  //h1_muGammaMass  = (TH1D*)_inFile->Get(("h1_diMuonMass"+_dataPeriod).c_str());
  h1_eGammaPt     = (TH1D*)_inFile->Get("h1_eGammaPtWeightCombined");
  h1_muGammaPt    = (TH1D*)_inFile->Get("h1_muGammaPtWeightCombined");
  h1_eGammaPV     = (TH1D*)_inFile->Get("h1_eGammaPVWeightCombined");
  h1_muGammaPV    = (TH1D*)_inFile->Get("h1_muGammaPVWeightCombined");
  h1_eGammaMass   = (TH1D*)_inFile->Get("h1_diElectronMassCombined");
  h1_muGammaMass  = (TH1D*)_inFile->Get("h1_diMuonMassCombined");

  phoEBSF_A = (TH2F*)_inFilePhoA_EB->Get("BarrelSF");
  phoEBSF_B = (TH2F*)_inFilePhoB_EB->Get("BarrelSF");
  phoEESF_A = (TH2F*)_inFilePhoA_EE->Get("EndcapSF");
  phoEESF_B = (TH2F*)_inFilePhoB_EE->Get("EndcapSF");

  // PU weights
  h1_S6to2011  = (TH1F*)_inFileS6to2011->Get("pileupWeights");
  h1_S6to2012  = (TH1F*)_inFileS6to2012->Get("pileupWeights");
  h1_S7to2011  = (TH1F*)_inFileS7to2011->Get("pileupWeights");
  h1_S7to2012  = (TH1F*)_inFileS7to2012->Get("pileupWeights");

  // Recoil weights
  //h1_recoilLongMuon       = (TH1D*)_inFile->Get("h1_muonRecoilTransWeight");
  //h1_recoilTransMuon      = (TH1D*)_inFile->Get("h1_muonRecoilTransWeight");
  //h1_recoilLongElectron   = (TH1D*)_inFile->Get("h1_electronRecoilTransWeight");
  //h1_recoilTransElectron  = (TH1D*)_inFile->Get("h1_electronRecoilTransWeight");

}

void WeightUtils::Initialize()
{
  _puWeight = 1.;
  _zzWeight = 1.;
  _glugluWeight = 1.;
  _vbfWeight = 1.;
  _recoWeight = 1.;
  _triggerWeight = 1.;
  _gammaPtPuWeight = 1.;
  _recoilLongWeight = 1.;
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

float WeightUtils::GetTotalWeight(int nPV, TLorentzVector l1, TLorentzVector l2)
{
  Initialize();
  float weight = 1.;

  if (!_isRealData) {
    weight *= PUWeight(nPV);
    weight *= RecoWeight(l1, l2);
  }
  return weight;
}

float WeightUtils::PUWeight(float nPU)
{
  //cout<<nPU<<endl;
  if (!_isRealData){
    if (nPU < 60 && _dataPeriod == "2011" ){
      _puWeight = h1_S6to2011->GetBinContent(h1_S6to2011->FindBin(nPU)); 
    } else if (nPU < 60 && _dataPeriod == "2012"){
      _puWeight = h1_S7to2012->GetBinContent(h1_S7to2012->FindBin(nPU)); 
    } else{
      _puWeight = 1; 
    }
  } else {
    _puWeight = 1;
  }
  return _puWeight;
}

float WeightUtils::RecoWeight(TLorentzVector l1, TLorentzVector l2)
{
  if (_selection == "muon" || _selection == "mumuGamma") {
    _triggerWeight = GetMuTriggerEff(l1)*GetMuTriggerEff(l2);
    _recoWeight    = 1.;
  }
  if (_selection == "electron" || _selection == "eeGamma" ) {
    _triggerWeight = 0.995;
    _recoWeight    = GetElectronEff(l1)*GetElectronEff(l2);
  }
  return _triggerWeight*_recoWeight;
}

float WeightUtils::GluGluHiggsWeight(float higgsPt, int higgsMass) 
{
  int iMass = (higgsMass - 250)/50;
  _glugluWeight = h1_Higgs[iMass]->GetBinContent(h1_Higgs[iMass]->FindBin(higgsPt));
  return _glugluWeight;
}

float WeightUtils::VBFHiggsWeight(float genMass, int higgsMass)
{
  _vbfWeight = float(higgsMass)/genMass;
  return _vbfWeight;
}


float WeightUtils::GammaWeightNew(TLorentzVector p1){
  float _gammaPtEtaWeight = 1;

  float gammaPtEta2012[4][5] = {
   //15<pt<20   20<pt<30    30<pt<40   40<pt<50   50<pt
    {0.9966,    0.9960,     0.9916,    0.9921,    0.9917},                        // 0.0<eta<0.8
    {0.9966,    1.0047,     0.9985,    0.9959,    0.9982},                        // 0.8<eta<1.4442
    {1.0071,    0.9991,     1.0039,    1.0021,    1.0026},                        // 1.566<eta<2.0
    {1.0071,    1.0077,     1.0075,    1.0048,    0.9989}                         // 2.0<eta<2.5
  };


  float gammaPtEta2011[4][5] = {
   //15<pt<20   20<pt<30    30<pt<40   40<pt<50   50<pt
    {1.0064,    0.9879,     0.9913,    0.9945,    0.9944},                        // 0.0<eta<0.8
    {1.0064,    0.9958,     0.9979,    0.9982,    0.9989},                        // 0.8<eta<1.4442
    {1.0188,    1.0113,     1.0016,    1.0000,    1.0041},                        // 1.566<eta<2.0
    {1.0188,    1.0081,     1.0037,    1.0025,    1.0006}                         // 2.0<eta<2.5
  };

  int ptBin = 0;
  int etaBin = 0;
  float binningPt[]  = {15., 20., 30., 40., 50., 99999.};
  float binningEta[] = {0., 0.8, 1.4442, 2.0, 2.5};
  float weight = 1.;

  for (int i = 0; i < 4; ++i) {
    if (fabs(p1.Eta()) > binningEta[i] && fabs(p1.Eta()) <= binningEta[i+1]) {
      etaBin = i;
      break;
    }
  }
  for (int i = 0; i < 5; ++i) {
    if (p1.Pt() > binningPt[i] && p1.Pt() <= binningPt[i+1]) {
      ptBin = i;
      break;
    }
  }

  if(!_isRealData){
    if(_dataPeriod == "2011"){
      _gammaPtEtaWeight = gammaPtEta2011[etaBin][ptBin];
    }else if (_dataPeriod == "2012"){
      _gammaPtEtaWeight = gammaPtEta2012[etaBin][ptBin];
    }
  }
  return _gammaPtEtaWeight;




}


float WeightUtils::GammaWeight(int nPV, TLorentzVector p1)
{
  _gammaPtPuWeight = 1;
  if (!_isRealData && (_selection == "eeGamma" || _selection == "mumuGamma")) {

    if (_dataPeriod == "2011A"){
      phoEBSF = phoEBSF_A;
      phoEESF = phoEESF_A;
    }else if (_dataPeriod == "2011B"){
      phoEBSF = phoEBSF_B;
      phoEESF = phoEESF_B;
    }

    if (fabs(p1.Eta()) < 1.45){
      int puMaxBin = phoEBSF->GetXaxis()->GetNbins();
      int puBin = phoEBSF->GetXaxis()->FindBin(nPV);
      if (puBin > puMaxBin) puBin = puMaxBin;
      if (puBin<1) puBin=1;
      int ptMaxBin = phoEBSF->GetYaxis()->GetNbins();
      int ptBin = phoEBSF->GetYaxis()->FindBin(p1.Pt());
      if (ptBin > ptMaxBin) ptBin = ptMaxBin;
      if (ptBin<1) ptBin=1;
      //cout << "max pu bin \t" << puMaxBin << "\t nVtx \t" << nGoodVtx      << "\t picked pu bin \t" << puBin << endl;
      //cout << "max pt bin \t" << ptMaxBin << "\t pt   \t" << phoEt[phoind] << "\t picked pt bin \t" << ptBin << endl;
      _gammaPtPuWeight = phoEBSF->GetBinContent(puBin,ptBin);
      //cout << "SF applied based on picked bin \t" << phoEBSF->GetBinContent(puBin,ptBin) << endl;

      //pixel seed veto section, barrel:
      if (_dataPeriod == "2011A") _gammaPtPuWeight *= .994;
      else                        _gammaPtPuWeight *= .990;
    }
    else{
      int puMaxBin = phoEESF->GetXaxis()->GetNbins();
      int puBin = phoEESF->GetXaxis()->FindBin(nPV);
      if (puBin > puMaxBin) puBin = puMaxBin;
      if (puBin<1) puBin=1;
      int ptMaxBin = phoEESF->GetYaxis()->GetNbins();
      int ptBin = phoEESF->GetYaxis()->FindBin(p1.Pt());
      if (ptBin > ptMaxBin) ptBin = ptMaxBin;
      if (ptBin<1) ptBin=1;
      //cout << puMaxBin << "\t" << nGoodVtx << "\t" << puBin << endl;
      //cout << ptMaxBin << "\t" << phoEt[phoind] << "\t" << ptBin << endl;
      _gammaPtPuWeight = phoEESF->GetBinContent(puBin,ptBin);
      //cout << phoEESF->GetBinContent(puBin,ptBin) << endl;

      //pixel seed veto section, endcap:
      if (_dataPeriod == "2011A") _gammaPtPuWeight *= .989;
      else                        _gammaPtPuWeight *= .978;
    }
  }
  return _gammaPtPuWeight;
}

float WeightUtils::GetMuTriggerEff(TLorentzVector l1) const
{
  //float doubleMu7Scale[4][4] = {
  //    {0.98, 0.974, 0.978, 0.978},
  //    {0.947, 0.947, 0.948, 0.947},
  //    {0.959, 0.958, 0.958, 0.956},
  //    {0.926, 0.925, 0.913, 0.921}
  //};

  float mu13_mu8Scale[4][4] = {
    {0.9753, 0.9751, 0.9746, 0.9745},
    {0.9627, 0.9568, 0.9563, 0.9575},
    {0.9490, 0.9482, 0.9503, 0.9462},
    {0.8594, 0.8614, 0.867, 0.8529}
  };

  int ptBin = 0;
  int etaBin = 0;
  float binningPt[]  = {20., 30., 40., 50., 99999.};
  float binningEta[] = {0., 0.8, 1.2, 2.1, 2.4};
  float weight = 1.;

  for (int i = 0; i < 4; ++i) {
    if (fabs(l1.Eta()) > binningEta[i] && fabs(l1.Eta()) <= binningEta[i+1]) {
      etaBin = i;
      break;
    }
  }
  for (int i = 0; i < 4; ++i) {
    if (l1.Pt() > binningPt[i] && l1.Pt() <= binningPt[i+1]) {
      ptBin = i;
      break;
    }
  }
  weight = mu13_mu8Scale[etaBin][ptBin];
  return weight;
}

float WeightUtils::GetElectronEff(TLorentzVector l1) const
{
  /*
  ////old shit
  float eleScale[2][5] = {
  {1.0, 0.938, 0.971, 0.980, 0.987},
  {1.0, 0.967, 0.97, 0.989, 0.989}
  };
  */

  //WP80 Data/MC scalefactors
  float eleScaleNew[5][5] = {
    {1.05, 1.02, 1.00, 0.98, 0.98},
    {1.19, 1.04, 0.97, 0.98, 0.98},
    {1.28, 1.11, 1.07, 1.05, 1.00},
    {1.30, 1.18, 1.07, 1.00, 0.98},
    {1.09, 1.09, 1.04, 1.00, 1.00}
  };
  int ptBin = 0;
  int etaBin = 0;
  float binningPt[] = {10., 15., 20., 25., 50., 9999.};
  float binningEta[] = {0., 1., 1.4442, 1.566, 2., 2.5};
  float weight = 1.;
  /*
     if (fabs(l1.Eta()) < 1.479) {
     etaBin = 0;
     } else {
     etaBin = 1;
     }
     */
  for (int i = 0; i < 5; ++i) {
    if (l1.Pt() > binningPt[i] && l1.Pt() <= binningPt[i+1]) {
      ptBin = i;
      break;
    }
  }
  for (int i = 0; i < 5; ++i) {
    if (fabs(l1.Eta()) >= binningEta[i] && fabs(l1.Eta()) <= binningEta[i+1]) {
      etaBin = i;
      break;
    }
  }

  weight = eleScaleNew[etaBin][ptBin];
  //weight = 1.;
  return weight;
}

float WeightUtils::GetPhotonMass() const
{
  float photonMass = 91.2;
  if (_selection == "eGamma") {
    photonMass = h1_eGammaMass->GetRandom();
  }
  if (_selection == "muGamma") {
    photonMass = h1_muGammaMass->GetRandom();
  }
  return photonMass;
}

