/*
   Utilities for retrieving weights for PU,etc.
 */

#ifndef _WeightUtils_H
#define _WeightUtils_H

#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "TROOT.h"
#include "TObject.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLorentzVector.h"
//#include "../../src/TCJet.h"

using namespace std;

class WeightUtils: public TObject {
    public:
        WeightUtils() {};
        virtual ~WeightUtils() {};
        WeightUtils(string sampleName, string dataPeriod, string abcd, string selection, bool isRealData);
        void  Initialize();
        void  SetDataBit(bool isRealData);
        void  SetDataPeriod(string dataPeriod);
        void  SetSampleName(string sampleName);
        void  SetSelection(string selection);

        float PUWeight(float nPU);
        float MuonSelectionWeight(TLorentzVector l1);
        float ElectronSelectionWeight(TLorentzVector l1);
        float GammaSelectionWeight(TLorentzVector l1, float SCEta);
        float HqtWeight(TLorentzVector l1);

        float ElectronTriggerWeight(TLorentzVector l1, TLorentzVector l2);
        float MuonTriggerWeight(TLorentzVector l1, TLorentzVector l2);

        ClassDef(WeightUtils, 0);

    private:
        //input parameters
        string _dataPeriod;
        string _abcd;
        string _sampleName;
        string _selection;
        bool   _isRealData;

        //sources
        TFile * _inFileS6to2011;
        TFile * _inFileS6to2012;       
        TFile * _inFileS7to2011;
        TFile * _inFileS7to2012;
        TFile * _inFileS6to2011obs;
        TFile * _inFileS7to2012obs;
        TFile * _inFileS7to2012obsV2;
        TFile * _inFileS10to2012ABTrue;
        TFile * _inFileS10to2012CDTrue;
        TFile * _inFileS10to2012ABCDTrue;
        TFile * _inFileS10to2012ABCDTrue_73500;
        TFile * _inFileSomeTuneto2012ABCDTrue;

        TFile *_MuTrig2011;

        TFile *_kFactors;

        TFile *_PhoMedWP2012;
        TFile *_PhoMedWPveto2012;

        TFile *_EleMoriondWP2012;



        TH2D * _HLTMu17Mu8_2011;

        TH1F * h1_S6to2011;
        TH1F * h1_S6to2012; 
        TH1F * h1_S7to2011; 
        TH1F * h1_S7to2012; 
        TH1F * h1_S7to2012obs; 
        TH1F * h1_S7to2012obsV2; 
        TH1F * h1_S6to2011obs;
        TH1F * h1_S10to2012ABTrue;
        TH1F * h1_S10to2012CDTrue;
        TH1F * h1_S10to2012ABCDTrue;
        TH1F * h1_S10to2012ABCDTrue_73500;
        TH1F * h1_SomeTuneto2012ABCDTrue;

        TH1D * kfact120_0;
        TH1D * kfact125_0;
        TH1D * kfact130_0;
        TH1D * kfact135_0;
        TH1D * kfact140_0;
        TH1D * kfact145_0;
        TH1D * kfact150_0;

        TH2F * h1_PhoMedWP2012;
        TH2F * h1_PhoMedWPveto2012;

        TH2F * h1_EleMoriondWP2012;


        //weights
        float _puWeight;
        float _zzWeight;
        float _glugluWeight;
        float _vbfWeight;
        float _recoWeight;
        float _triggerWeight;
        float _recoilLongWeight;
        float _recoilTransWeight;
        float _gammaPtPuWeight;

};

#endif

#if !defined(__CINT__)
ClassImp(WeightUtils);
#endif
