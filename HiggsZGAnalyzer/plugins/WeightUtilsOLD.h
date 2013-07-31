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
#include "../../src/TCJet.h"

using namespace std;

class WeightUtils: public TObject {
    public:
        WeightUtils() {};
        virtual ~WeightUtils() {};
        WeightUtils(string sampleName, string dataPeriod, string selection, bool isRealData);
        void  Initialize();
        void  SetDataBit(bool isRealData);
        void  SetDataPeriod(string dataPeriod);
        void  SetSampleName(string sampleName);
        void  SetSelection(string selection);

        float PUWeight(float nPU);
        float RecoWeight(TLorentzVector l1, TLorentzVector l2);
        float GammaWeight(int nPV, TLorentzVector p1);
        float GammaWeightNew(TLorentzVector p1);
        float GluGluHiggsWeight(float higgsPt, int higgsMass);
        float VBFHiggsWeight(float genMass, int higgsMass);
        float GetTotalWeight(int nPV, TLorentzVector l1, TLorentzVector l2);

        float GetElectronEff(TLorentzVector l1) const;
        float GetMuTriggerEff(TLorentzVector l1) const;
        float GetPhotonMass() const;

        float GetGammaPtPuWeight() {return _gammaPtPuWeight;}

        ClassDef(WeightUtils, 0);

    private:
        //input parameters
        string _dataPeriod;
        string _sampleName;
        string _selection;
        bool   _isRealData;

        //sources
        TFile * _inFile;

        TFile * _inFileS6to2011;
        TFile * _inFileS6to2012;       
        TFile * _inFileS7to2011;
        TFile * _inFileS7to2012;

        TFile *_inFilePhoA_EE;
        TFile *_inFilePhoA_EB;
        TFile *_inFilePhoB_EE;
        TFile *_inFilePhoB_EB;

        TH1F * h1_S6to2011;
        TH1F * h1_S6to2012; 
        TH1F * h1_S7to2011; 
        TH1F * h1_S7to2012; 


        TH1D  *h1_eGammaPt;   
        TH1D  *h1_muGammaPt;  
        TH1D  *h1_eGammaPV;   
        TH1D  *h1_muGammaPV;  
        TH1D  *h1_eGammaMass; 
        TH1D  *h1_muGammaMass;
        TH1F  *fixed;

        TH1D  *h1_Higgs[8];

        TH2F  *phoEBSF_A;
        TH2F  *phoEBSF_B;
        TH2F  *phoEESF_A;
        TH2F  *phoEESF_B;

        TH2F  *phoEESF;
        TH2F  *phoEBSF;

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
