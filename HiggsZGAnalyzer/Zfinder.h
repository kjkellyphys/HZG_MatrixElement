//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 17 13:59:30 2012 by ROOT version 5.32/00
// from TTree eventTree/eventTree
// found on file: /eos/uscms/store/user/devildog/V01_10_8TeV/DYJetsToLLV2/nuTuple_100_1_4DF.root
//////////////////////////////////////////////////////////

#ifndef Zfinder_h
#define Zfinder_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.
#include <TClonesArray.h>
#include <TVector3.h>
#include <vector>
#include "plugins/HistManager.h"

#include "../src/TCJet.h"
#include "../src/TCMET.h"
#include "../src/TCElectron.h"
#include "../src/TCMuon.h"
#include "../src/TCTau.h"
#include "../src/TCPhoton.h"
#include "../src/TCGenJet.h"
#include "../src/TCPrimaryVtx.h"
#include "../src/TCTriggerObject.h"
#include "../src/TCGenParticle.h"
// Fixed size dimensions of array or collections stored in the TTree if any.

class Zfinder : public TSelector {
  public :
    TFile* histoFile;
    HistManager* hm;

    float         unskimmedEvents;
    float unskimmedEventsTotal;
    int fileCount;
    TTree *thisTree;
    TFile *file0;
    TH1F * h1_numOfEvents;
    TVector3 *pvPosition;

    TTree          *fChain;   //!pointer to the analyzed TTree or TChain

    // Declaration of leaf types
    TClonesArray    *recoJets;
    TClonesArray    *recoJPT;
    TClonesArray    *recoElectrons;
    TClonesArray    *recoMuons;
    TClonesArray    *pfMuons;
    TClonesArray    *recoTaus;
    TClonesArray    *recoPhotons;
    TClonesArray    *pfPhotons;
    TCMET           *recoMET;
    TCMET           *recoMETNoPU;
    TClonesArray    *triggerObjects;
    TClonesArray    *genJets;
    TClonesArray    *genParticles;
    TClonesArray    *primaryVtx;
    TVector3        *beamSpot;
    Int_t           nPUVertices;
    Float_t         nPUVerticesTrue;
    Bool_t          isRealData;
    ULong64_t       runNumber;
    ULong64_t       eventNumber;
    UInt_t          lumiSection;
    UInt_t          bunchCross;
    Bool_t          isScraping;
    Bool_t          isNoiseHcal;
    Bool_t          isDeadEcalCluster;
    Bool_t          isCSCTightHalo;
    Bool_t          isCSCLooseHalo;
    Float_t         ptHat;
    Float_t         qScale;
    Float_t         evtWeight;
    Float_t         rhoFactor;
    Float_t         rhoFactorMu;
    ULong64_t       triggerStatus;
    UInt_t          hltPrescale[64];
    vector<string>  *mapCuts;

    // List of branches
    TBranch        *b_recoJets;   //!
    TBranch        *b_recoJPT;   //!
    TBranch        *b_recoElectrons;   //!
    TBranch        *b_recoMuons;   //!
    TBranch        *b_pfMuons;   //!
    TBranch        *b_recoTaus;   //!
    TBranch        *b_recoPhotons;   //!
    TBranch        *b_pfPhotons;   //!
    TBranch        *b_recoMET;   //!
    TBranch        *b_recoMETNoPU;   //!
    TBranch        *b_triggerObjects;   //!
    TBranch        *b_genJets;   //!
    TBranch        *b_genParticles;   //!
    TBranch        *b_primaryVtx;   //!
    TBranch        *b_beamSpot;   //!
    TBranch        *b_nPUVertices;   //!
    TBranch        *b_nPUVerticesTrue;   //!
    TBranch        *b_isRealData;   //!
    TBranch        *b_runNumber;   //!
    TBranch        *b_eventNumber;   //!
    TBranch        *b_lumiSection;   //!
    TBranch        *b_bunchCross;   //!
    TBranch        *b_isScraping;   //!
    TBranch        *b_isNoiseHcal;   //!
    TBranch        *b_isDeadEcalCluster;   //!
    TBranch        *b_isCSCTightHalo;   //!
    TBranch        *b_isCSCLooseHalo;   //!
    TBranch        *b_ptHat;   //!
    TBranch        *b_qScale;   //!
    TBranch        *b_evtWeight;   //!
    TBranch        *b_rhoFactor;   //!
    TBranch        *b_rhoFactorMu;   //!
    TBranch        *b_triggerStatus;   //!
    TBranch        *b_hltPrescale;   //!
    TBranch        *b_mapCuts;   //!

    Zfinder(TTree * /*tree*/ =0) : fChain(0) { }
    virtual ~Zfinder() { }
    virtual Int_t   Version() const { return 2; }
    virtual void    Begin(TTree *tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Notify();
    virtual Bool_t  Process(Long64_t entry);
    virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
    virtual void    SetOption(const char *option) { fOption = option; }
    virtual void    SetObject(TObject *obj) { fObject = obj; }
    virtual void    SetInputList(TList *input) { fInput = input; }
    virtual TList  *GetOutputList() const { return fOutput; }
    virtual void    SlaveTerminate();
    virtual void    Terminate();

    float EAMu[6]; 

    struct muIDCuts{
      float IsPF;
      float IsGLB;
      float NormalizedChi2;
      float NumberOfValidMuonHits;
      float NumberOfMatchedStations;
      float NumberOfValidPixelHits;
      float TrackLayersWithMeasurement;
      float dxy;
      float dz;
      string cutName;
    } tightMuID;

    struct muIsoCuts{
      float chIso04;
      float nhIso04;
      float phIso04;
      float relCombIso04;
      string cutName;
    } tightMuIso;

    virtual bool      PassMuonID(TCMuon *mu, muIDCuts cutLevel);
    virtual bool      PassMuonIso(TCMuon *mu, muIsoCuts cutLevel);

    ClassDef(Zfinder,0);
};

#endif

#ifdef Zfinder_cxx
void Zfinder::Init(TTree *tree)
{
  fileCount = 0;
  unskimmedEvents = 0;
  unskimmedEventsTotal = 0;

  // Set object pointer
  recoJets = 0;
  recoJPT = 0;
  recoElectrons = 0;
  recoMuons = 0;
  pfMuons = 0;
  recoTaus = 0;
  recoPhotons = 0;
  pfPhotons = 0;
  recoMET = 0;
  recoMETNoPU = 0;
  triggerObjects = 0;
  genJets = 0;
  genParticles = 0;
  primaryVtx = 0;
  beamSpot = 0;
  mapCuts = 0;
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  thisTree = tree;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("recoJets", &recoJets, &b_recoJets);
  fChain->SetBranchAddress("recoJPT", &recoJPT, &b_recoJPT);
  fChain->SetBranchAddress("recoElectrons", &recoElectrons, &b_recoElectrons);
  fChain->SetBranchAddress("recoMuons", &recoMuons, &b_recoMuons);
  fChain->SetBranchAddress("pfMuons", &pfMuons, &b_pfMuons);
  fChain->SetBranchAddress("recoTaus", &recoTaus, &b_recoTaus);
  fChain->SetBranchAddress("recoPhotons", &recoPhotons, &b_recoPhotons);
  fChain->SetBranchAddress("pfPhotons", &pfPhotons, &b_pfPhotons);
  fChain->SetBranchAddress("recoMET", &recoMET, &b_recoMET);
  fChain->SetBranchAddress("recoMETNoPU", &recoMETNoPU, &b_recoMETNoPU);
  fChain->SetBranchAddress("triggerObjects", &triggerObjects, &b_triggerObjects);
  fChain->SetBranchAddress("genJets", &genJets, &b_genJets);
  fChain->SetBranchAddress("genParticles", &genParticles, &b_genParticles);
  fChain->SetBranchAddress("primaryVtx", &primaryVtx, &b_primaryVtx);
  fChain->SetBranchAddress("beamSpot", &beamSpot, &b_beamSpot);
  fChain->SetBranchAddress("nPUVertices", &nPUVertices, &b_nPUVertices);
  fChain->SetBranchAddress("nPUVerticesTrue", &nPUVerticesTrue, &b_nPUVerticesTrue);
  fChain->SetBranchAddress("isRealData", &isRealData, &b_isRealData);
  fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
  fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
  fChain->SetBranchAddress("lumiSection", &lumiSection, &b_lumiSection);
  fChain->SetBranchAddress("bunchCross", &bunchCross, &b_bunchCross);
  fChain->SetBranchAddress("isScraping", &isScraping, &b_isScraping);
  fChain->SetBranchAddress("isNoiseHcal", &isNoiseHcal, &b_isNoiseHcal);
  fChain->SetBranchAddress("isDeadEcalCluster", &isDeadEcalCluster, &b_isDeadEcalCluster);
  fChain->SetBranchAddress("isCSCTightHalo", &isCSCTightHalo, &b_isCSCTightHalo);
  fChain->SetBranchAddress("isCSCLooseHalo", &isCSCLooseHalo, &b_isCSCLooseHalo);
  fChain->SetBranchAddress("ptHat", &ptHat, &b_ptHat);
  fChain->SetBranchAddress("qScale", &qScale, &b_qScale);
  fChain->SetBranchAddress("evtWeight", &evtWeight, &b_evtWeight);
  fChain->SetBranchAddress("rhoFactor", &rhoFactor, &b_rhoFactor);
  fChain->SetBranchAddress("rhoFactorMu", &rhoFactorMu, &b_rhoFactorMu);
  fChain->SetBranchAddress("triggerStatus", &triggerStatus, &b_triggerStatus);
  fChain->SetBranchAddress("hltPrescale", hltPrescale, &b_hltPrescale);
  fChain->SetBranchAddress("mapCuts", &mapCuts, &b_mapCuts);

  ///////////////
  // Cuts Init //
  ///////////////


  tightMuID.cutName =                     "tightMuID";
  tightMuID.IsPF =                        1;
  tightMuID.IsGLB =                       1;
  tightMuID.NormalizedChi2 =              10;
  tightMuID.NumberOfValidMuonHits =       0;
  tightMuID.NumberOfMatchedStations =     1;
  tightMuID.NumberOfValidPixelHits =      0;
  tightMuID.TrackLayersWithMeasurement =  5;
  tightMuID.dxy =                         0.2;
  tightMuID.dz =                          0.5;

  tightMuIso.cutName =                    "tightMuIso";
  tightMuIso.chIso04 =                    99999;
  tightMuIso.nhIso04 =                    99999;
  tightMuIso.phIso04 =                    99999;
  tightMuIso.relCombIso04 =               0.12;

}

Bool_t Zfinder::Notify()
{
  fileCount+= 1;
  file0 = thisTree->GetCurrentFile();
  h1_numOfEvents = (TH1F*) file0->Get("ntupleProducer/numOfEvents");
  unskimmedEvents = h1_numOfEvents->GetBinContent(1);
  cout<<"THIS IS FILE NUMBER: "<<fileCount<<" and it has this many events: "<<unskimmedEvents<<endl;
  unskimmedEventsTotal += unskimmedEvents;
  return kTRUE;
}

#endif // #ifdef Zfinder_cxx
