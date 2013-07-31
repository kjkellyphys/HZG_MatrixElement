//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May 30 06:57:44 2013 by ROOT version 5.32/00
// from TTree eventTree/eventTree
// found on file: /eos/uscms/store/user/bpollack/V02_01_8TeV/ggHZG_M125/nuTuple_10_2_xWc.root
//////////////////////////////////////////////////////////

#ifndef simple_v2_h
#define simple_v2_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.
#include <TClonesArray.h>
#include <TVector3.h>
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
#include "plugins/HistManager.h"
#include "plugins/ZGAngles.h"

// Fixed size dimensions of array or collections stored in the TTree if any.

class simple_v2 : public TSelector {
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
    float weight;


    TTree          *fChain;   //!pointer to the analyzed TTree or TChain

    // Declaration of leaf types
    TClonesArray    *recoJets;
    TClonesArray    *recoJPT;
    TClonesArray    *recoElectrons;
    TClonesArray    *recoMuons;
    TClonesArray    *recoTaus;
    TClonesArray    *recoPhotons;
    TCMET           *recoMET;
    TClonesArray    *triggerObjects;
    TClonesArray    *genJets;
    TClonesArray    *genParticles;
    TClonesArray    *primaryVtx;
    TVector3        *beamSpot;
    Int_t           nPUVertices;
    Float_t         nPUVerticesTrue;
    Bool_t          isRealData;
    UInt_t          runNumber;
    ULong64_t       eventNumber;
    UInt_t          lumiSection;
    UInt_t          bunchCross;
    Float_t         ptHat;
    Float_t         qScale;
    Float_t         evtWeight;
    Float_t         rhoFactor;
    Float_t         rho25Factor;
    Float_t         rhoMuFactor;
    ULong64_t       triggerStatus;
    UInt_t          hltPrescale[64];
    Bool_t          NoiseFilters_isScraping;
    Bool_t          NoiseFilters_isNoiseHcalHBHE;
    Bool_t          NoiseFilters_isNoiseHcalLaser;
    Bool_t          NoiseFilters_isNoiseEcalTP;
    Bool_t          NoiseFilters_isNoiseEcalBE;
    Bool_t          NoiseFilters_isCSCTightHalo;
    Bool_t          NoiseFilters_isCSCLooseHalo;

    // List of branches
    TBranch        *b_recoJets;   //!
    TBranch        *b_recoJPT;   //!
    TBranch        *b_recoElectrons;   //!
    TBranch        *b_recoMuons;   //!
    TBranch        *b_recoTaus;   //!
    TBranch        *b_recoPhotons;   //!
    TBranch        *b_recoMET;   //!
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
    TBranch        *b_ptHat;   //!
    TBranch        *b_qScale;   //!
    TBranch        *b_evtWeight;   //!
    TBranch        *b_rhoFactor;   //!
    TBranch        *b_rho25Factor;   //!
    TBranch        *b_rhoMuFactor;   //!
    TBranch        *b_triggerStatus;   //!
    TBranch        *b_hltPrescale;   //!
    TBranch        *b_NoiseFilters;   //!

    simple_v2(TTree * /*tree*/ =0) : fChain(0) { }
    virtual ~simple_v2() { }
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

    static bool P4SortCondition(const TLorentzVector& p1, const TLorentzVector& p2) {return (p1.Pt() > p2.Pt());} 

    //////////////////////
    // Cut and Iso Defs //
    //////////////////////

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
    } tightMuIso, looseMuIso;

    virtual bool      PassMuonID(TCMuon *mu, muIDCuts cutLevel);
    virtual bool      PassMuonIso(TCMuon *mu, muIsoCuts cutLevel);

    struct genHZGParticles{
      TCGenParticle* h;
      TCGenParticle* z;
      TCGenParticle* w;
      TCGenParticle* g;
      TCGenParticle* lp;
      TCGenParticle* lm;
    }genHZG;

    void FindGenParticles(TClonesArray* genParticles, string selection, vector<TCGenParticle>& vetoPhotons, genHZGParticles& _genHZG);
    void CleanUpGen(genHZGParticles& _genHZG);
    virtual void    AnglePlots(ZGAngles &zga, float eventWeight);

    ClassDef(simple_v2,0);
};

#endif

#ifdef simple_v2_cxx
void simple_v2::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer
  recoJets = 0;
  recoJPT = 0;
  recoElectrons = 0;
  recoMuons = 0;
  recoTaus = 0;
  recoPhotons = 0;
  recoMET = 0;
  triggerObjects = 0;
  genJets = 0;
  genParticles = 0;
  primaryVtx = 0;
  beamSpot = 0;
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  thisTree = tree;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("recoJets", &recoJets, &b_recoJets);
  fChain->SetBranchAddress("recoJPT", &recoJPT, &b_recoJPT);
  fChain->SetBranchAddress("recoElectrons", &recoElectrons, &b_recoElectrons);
  fChain->SetBranchAddress("recoMuons", &recoMuons, &b_recoMuons);
  fChain->SetBranchAddress("recoTaus", &recoTaus, &b_recoTaus);
  fChain->SetBranchAddress("recoPhotons", &recoPhotons, &b_recoPhotons);
  fChain->SetBranchAddress("recoMET", &recoMET, &b_recoMET);
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
  fChain->SetBranchAddress("ptHat", &ptHat, &b_ptHat);
  fChain->SetBranchAddress("qScale", &qScale, &b_qScale);
  fChain->SetBranchAddress("evtWeight", &evtWeight, &b_evtWeight);
  fChain->SetBranchAddress("rhoFactor", &rhoFactor, &b_rhoFactor);
  fChain->SetBranchAddress("rho25Factor", &rho25Factor, &b_rho25Factor);
  fChain->SetBranchAddress("rhoMuFactor", &rhoMuFactor, &b_rhoMuFactor);
  fChain->SetBranchAddress("triggerStatus", &triggerStatus, &b_triggerStatus);
  fChain->SetBranchAddress("hltPrescale", hltPrescale, &b_hltPrescale);
  fChain->SetBranchAddress("NoiseFilters", &NoiseFilters_isScraping, &b_NoiseFilters);


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

  looseMuIso.cutName =                    "looseMuIso";
  looseMuIso.chIso04 =                    99999;
  looseMuIso.nhIso04 =                    99999;
  looseMuIso.phIso04 =                    99999;
  looseMuIso.relCombIso04 =               0.4;

}

Bool_t simple_v2::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  fileCount+= 1;
  file0 = thisTree->GetCurrentFile();
  h1_numOfEvents = (TH1F*) file0->Get("ntupleProducer/numOfEvents");
  unskimmedEvents = h1_numOfEvents->GetBinContent(1);
  cout<<"THIS IS FILE NUMBER: "<<fileCount<<" and it has this many events: "<<unskimmedEvents<<endl;
  unskimmedEventsTotal += unskimmedEvents;
  return kTRUE;
}



#endif // #ifdef simple_v2_cxx
