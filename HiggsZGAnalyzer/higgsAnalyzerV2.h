#ifndef higgsAnalyzerV2_h
#define higgsAnalyzerV2_h

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <set>
#include <map>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH2.h>
#include <TStyle.h>

#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

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
#include "plugins/MetCorrector.h"
#include "plugins/TriggerSelector.h"
#include "plugins/WeightUtils.h"
#include "plugins/ElectronFunctions.h"
#include "plugins/rochcor_2011.h"
#include "plugins/rochcor2012v2.h"
#include "plugins/PhosphorCorrectorFunctor.hh"
#include "plugins/LeptonScaleCorrections.h"
#include "plugins/EGammaMvaEleEstimator.cc"
#include "plugins/ZGAngles.h"

#ifdef __MAKECINT__
#pragma link C++ class vector<string>+;
#endif

class higgsAnalyzerV2 : public TSelector {

	private:

		TFile* histoFile;
		TFile* trainingFile;
		TFile* sampleFile;
		TFile* higgsFile;
    TFile* eleSmearFile;
    TFile* eleIDISOFile;
    TFile* m_llgFile;
    HistManager* hm;
    HistManager* hmHiggs;
    HistManager* hmEleIDISO;
    HistManager* hmEleSmear;

    ofstream muDump1;
    ofstream elDump2;
    ofstream elDumpMVA;
    ofstream phDump1;
    ofstream phDump2;

    ofstream muDumpFinal;
    ofstream elDumpFinal;

    ofstream dataDump;

    ofstream finalDump;

    bool electronDump;
    bool muonDump;
    
    float unskimmedEventsTotal;
    int fileCount;

    TVector3 *pvPosition;

    // Random number generator
    TRandom3* rnGenerator;
    TRandom3* rEl;
    TRandom3* rMuRun;

    // Selectors
    //GenParticleSelector genParticleSelector;
    TriggerSelector *triggerSelector;

    // Utilities
    WeightUtils *weighter;
    rochcor_2011 * rmcor2011;
    rochcor2012 * rmcor2012;
    zgamma::PhosphorCorrectionFunctor* phoCorrector;

    //Unskimmed events
    TTree *thisTree;
    TFile *file0;
    TH1F * h1_numOfEvents;

    //R9
    float R9;
    float R9Cor;

    //ElectronMVA selection
    EGammaMvaEleEstimator* myMVATrig;

    TMVA::Reader             *myTMVAReader;

    //ZGAngles MVA shit



  public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    TTree          *sampleChain; 
    TTree          *trainingChain; 
    TTree          *m_llgChain;

    //MVA Branches
    float          diffZGscalar;
    float          diffZGvector;
    float          threeBodyMass;
    float          threeBodyPt;
    float          divPt;
    float          cosZ;
    float          cosG;
    float          METdivQt;
    float          GPt;
    float          ZPt;
    float          DPhi;
    float          diffPlaneMVA;
    float          vtxVariable;
    float          dr1;
    float          dr2;
    float          M12;
    float          scaleFactor;

    //M_llg branches
    float          m_llg;
    float          m_llgCAT1;
    float          m_llgCAT2;
    float          m_llgCAT3;
    float          m_llgCAT4;
    float          unBinnedWeight;
    float          unBinnedLumiXS;

		// Declaration of leaf types
		TClonesArray  *recoJets;
		TClonesArray  *recoMET;
		TClonesArray  *genJets;
		TClonesArray  *genParticles;
		TClonesArray  *recoMuons;
		TClonesArray  *recoElectrons;
		TClonesArray  *recoPhotons;
		TClonesArray  *primaryVtx;
		TClonesArray  *triggerObjects;
    float         unskimmedEvents;
		int           runNumber;
		int           eventNumber;
		int           lumiSection;
		int           nPUVertices;
		float         nPUVerticesTrue;
		double        ptHat;
		float         rhoFactor;
		ULong64_t     triggerStatus;
		unsigned int           hltPrescale[64];
		bool          isRealData;

		// List of branches
		TBranch        *b_recoJets;   //!
		TBranch        *b_recoMET;   //!
		TBranch        *b_genJets;   //!
		TBranch        *b_genParticles;   //!
		TBranch        *b_recoMuons;     //!
		TBranch        *b_recoElectrons;  //!
		TBranch        *b_recoPhotons;    //!
		TBranch        *b_primaryVtx;   //!
		TBranch        *b_triggerObjects;   //!
		TBranch        *b_runNumber;   //!
		TBranch        *b_eventNumber;   //!
		TBranch        *b_lumiSection;   //!
		TBranch        *b_rhoFactor;   //!
		TBranch        *b_ptHat;   //!
		TBranch        *b_triggerStatus;   //!
		TBranch        *b_hltPrescale;   //!
		TBranch        *b_nPUVertices;   //!
		TBranch        *b_nPUVerticesTrue;   //!
		TBranch        *b_isRealData;   //!


		//For counting events
    int          nEvents[100];
    int          genAccept[2];

		higgsAnalyzerV2(TTree * /*tree*/ =0): fChain(0){ }
		virtual ~higgsAnalyzerV2() { }
		virtual int     Version() const { return 2; }
		virtual void    Begin(TTree *tree);
		//virtual void    SlaveBegin(TTree *tree) { TString option = GetOption();};
		virtual void    Init(TTree *tree);
		virtual bool    Notify();
		virtual bool    Process(Long64_t entry);
		virtual int     GetEntry(Long64_t entry, int getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
		virtual void    SetOption(const char *option) { fOption = option; }
		virtual void    SetObject(TObject *obj) { fObject = obj; }
		virtual void    SetInputList(TList *input) { fInput = input; }
		virtual TList  *GetOutputList() const { return fOutput; }
		virtual void    SlaveTerminate() {};
		virtual void    Terminate();

    static bool P4SortCondition(const TLorentzVector& p1, const TLorentzVector& p2) {return (p1.Pt() > p2.Pt());} 
    static bool VertexSortCondition(const TCPrimaryVtx& pv1, const TCPrimaryVtx& pv2) {return (pv1.SumPt2Trks() > pv2.SumPt2Trks());}

    virtual float   Dz(TVector3 objVtx, TLorentzVector objP4, TVector3 vtx);
    virtual float   Dxy(TVector3 objVtx, TLorentzVector objP4, TVector3 vtx);
		virtual void    MetPlusZPlots(TLorentzVector metP4, TLorentzVector ZP4, float evtWeight);
		virtual void    MetPlusLeptonPlots(TLorentzVector metP4, TLorentzVector p1, TLorentzVector p2, float evtWeight);
		virtual void    LeptonBasicPlots(TLorentzVector p1, TLorentzVector p2, float evtWeight);
		virtual void    GenPlots(vector<TCGenParticle> Zs, vector<TCGenParticle> leps, vector<TCGenParticle> phots, vector<TCGenParticle> Hs, TLorentzVector ZP4,TLorentzVector GP4, float evtWeight); 
    virtual void    StandardPlots(TLorentzVector p1, TLorentzVector p2, TLorentzVector gamma, float evtWeight,string tag, string folder);
    virtual void    AnglePlots(ZGAngles &zga, float eventWeight);
		virtual void    DileptonBasicPlots(TLorentzVector ZP4, float evtWeight);
		virtual bool    CosmicMuonFilter(TCMuon muon1, TCMuon muon2);
		virtual float   CalculateTransMass(TLorentzVector p1, TLorentzVector p2);
		virtual float   CalculateTransMassAlt(TLorentzVector p1, TLorentzVector p2);
    virtual float   GetPhotonMass();

    virtual float          CalculateM12sqrd(TLorentzVector p1, TLorentzVector p2);
    virtual float          CalculateX1(TLorentzVector p1,TLorentzVector p2);
    virtual float          CalculateX2(TLorentzVector p1, TLorentzVector p2);
    virtual float          Zeppenfeld(TLorentzVector p, TLorentzVector pj1, TLorentzVector pj2);

    virtual bool           FindGoodZElectron(vector<TCElectron> electronList, TLorentzVector* lepton1, TLorentzVector* lepton2, TLorentzVector* ZP4, float* eta1, float* eta2, int* int1, int* int2); 
    virtual bool           FindGoodZMuon(vector<TCMuon> muonList, TLorentzVector* lepton1, TLorentzVector* lepton2, TLorentzVector* ZP4, int* int1, int* int2 ); 
    virtual bool           FindGoodZElectron(vector<TCElectron> electronList, vector<TCElectron> uncorElectronList, TLorentzVector* lepton1, TLorentzVector* lepton2, TLorentzVector* uncorLepton1, TLorentzVector* uncorLepton2, TLorentzVector* ZP4, float* eta1, float* eta2, int* int1, int* int2); 
    virtual bool           FindGoodZMuon(vector<TCMuon> muonList, vector<TCMuon> uncorMuonList, TLorentzVector* lepton1, TLorentzVector* lepton2, TLorentzVector* uncorLepton1, TLorentzVector* uncorLepton2, TLorentzVector* ZP4, int* int1, int* int2); 

    /////////////////////
    // Gen Level Stuff //
    /////////////////////

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

    //////////////////////
    // Cut and Iso Defs //
    //////////////////////

    float EAMu[6]; 
    float EAEle[7]; 
    float EAPho[7][3];

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

    struct elIDCuts{
      //broken into [0] barrel and [1] endcap
      float dEtaIn[2];
      float dPhiIn[2];
      float sigmaIetaIeta[2];
      float HadOverEm[2];
      float dxy[2];
      float dz[2];
      float fabsEPDiff[2];
      float ConversionMissHits[2];
      float PassedConversionProb[2];
      string cutName;
    } vetoElID,looseElID;

    struct elIsoCuts{
      float chIso04;
      float nhIso04;
      float phIso04;
      float relCombIso04;
      string cutName;
    } looseElIso;

    struct phIDCuts{
      //broken into [0] barrel and [1] endcap
      float PassedEleSafeVeto[2];
      float HadOverEm[2];
      float sigmaIetaIeta[2];
      string cutName;
    } loosePhID, mediumPhID;

    struct phIsoCuts{
      float chIso03[2];
      float nhIso03[2];
      float phIso03[2];
      float relCombIso03[2];
      string cutName;
    } loosePhIso, mediumPhIso;

    virtual bool      PassMuonID(TCMuon *mu, muIDCuts cutLevel);
    virtual bool      PassMuonIso(TCMuon *mu, muIsoCuts cutLevel);
    virtual bool      PassElectronID(TCElectron *el, elIDCuts cutLevel);
    virtual bool      PassElectronIso(TCElectron *el, elIsoCuts cutLevel);
    virtual bool      PassPhotonID(TCPhoton *ph, phIDCuts cutLevel);
    virtual bool      PassPhotonIso(TCPhoton *ph, phIsoCuts cutLevel);

    virtual void      LumiXSWeight(float *lumiXS);


    /////////////////////
    // Debugging Dumps //
    /////////////////////

    virtual void  ElectronDump(TCElectron *el, elIDCuts cutLevelID, elIsoCuts cutLevelIso, ofstream & dump);
    virtual void  MVADumper(TCElectron *ele, EGammaMvaEleEstimator* mvaMaker, double rhoFactor, ofstream & dump);
    virtual void  MuonDump(TCMuon *mu, muIDCuts cutLevelID, muIsoCuts cutLevelIso, ofstream & dump);
    virtual void  PhotonDump(TCPhoton *ph, phIDCuts cutLevelID, phIsoCuts cutLevelIso, ofstream & dump);
    virtual void  PhotonDump2(TCPhoton *ph, phIDCuts cutLevelID, phIsoCuts cutLevelIso, TLorentzVector lepton1, TLorentzVector lepton2, ofstream & dump);
    virtual void  DataDumper(TLorentzVector* lepton1, TLorentzVector* lepton2, TLorentzVector* gamma, float R9, float SCEta, ofstream & dump, float eta1, float eta2);
    virtual void  DataDumper(TLorentzVector* lepton1, TLorentzVector* lepton2, TLorentzVector* gamma,
        TLorentzVector* uncorLepton1, TLorentzVector* uncorLepton2, TLorentzVector* uncorGamma, float R9, float SCEta, ofstream & dump, float eta1, float eta2);
    virtual void  FinalDumper(TLorentzVector* lepton1, TLorentzVector* lepton2, TLorentzVector* gamma, int catNum, ofstream & dump);


    ///////////////////////
    // ZGAngles MVA defs //
    ///////////////////////

    // TMVA weights directory

#define N_DISCR_METHODS 3
#define N_HIGGS_MASSES 1

    // here we will use only BDTG... but will keep the structure 
    enum DISCR_METHOD {
      MLPBNN, BDTG, BDT
    };

    enum DISCR_TYPE {
      D_ZJets
    };


    TMVA::Reader* tmvaReader;

    struct mvaInitStruct{
      TString weightsDir; 
      TString discrMethodName[3];
      TString discrSampleName; 
      Int_t mvaHiggsMassPoint[N_HIGGS_MASSES];
      Float_t bdtCut[N_HIGGS_MASSES];
    } mvaInits;

    struct mvaVarStruct{
      Float_t _diffZGvector;
      Float_t _threeBodyPt;
      Float_t _GPt;
      Float_t _cosZ;
      Float_t _diffPlaneMVA;
      Float_t _vtxVariable;
      Float_t _threeBodyMass;
      Float_t _dr1;
      Float_t _dr2;
      Float_t _M12;
    } mvaVars; 

    virtual TMVA::Reader*   MVAInitializer(mvaVarStruct vars, mvaInitStruct inits);
    virtual void MVACalulator (mvaVarStruct vars, mvaInitStruct inits, TMVA::Reader* tmvaReader);




    ClassDef(higgsAnalyzerV2,0);
};

#endif

#ifdef higgsAnalyzerV2_cxx
void higgsAnalyzerV2::Init(TTree *tree)
{
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the branch addresses and branch
	// pointers of the tree will be set.
	// It is normally not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).

	// Set object pointer
  fileCount = 0;
  unskimmedEvents = 0;
  unskimmedEventsTotal = 0;
	recoJets = 0;
	recoMET = 0;
	genJets = 0;
	genParticles = 0;
	primaryVtx = 0;
	recoMuons = 0;
	recoElectrons = 0;
	recoPhotons = 0;
  triggerObjects = 0;
	// Set branch addresses and branch pointers
	if (!tree) return;
  thisTree = tree;
  //TTree *cloneTree = new TTree(&tree->CloneTree());
	//fChain = tree->CloneTree();
	fChain = tree;
	fChain->SetMakeClass(1);
  cout<<"init"<<endl;

	fChain->SetBranchAddress("recoJets", &recoJets, &b_recoJets);
	fChain->SetBranchAddress("recoElectrons", &recoElectrons, &b_recoElectrons);
	fChain->SetBranchAddress("recoMuons", &recoMuons, &b_recoMuons);
	fChain->SetBranchAddress("recoPhotons", &recoPhotons, &b_recoPhotons);
	fChain->SetBranchAddress("recoMET", &recoMET, &b_recoMET);
	fChain->SetBranchAddress("genJets", &genJets, &b_genJets);
	fChain->SetBranchAddress("genParticles", &genParticles, &b_genParticles);
	fChain->SetBranchAddress("primaryVtx", &primaryVtx, &b_primaryVtx);
	fChain->SetBranchAddress("triggerObjects", &triggerObjects, &b_triggerObjects);

	fChain->SetBranchAddress("isRealData", &isRealData, &b_isRealData);
  fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
	fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
	fChain->SetBranchAddress("lumiSection", &lumiSection, &b_lumiSection);
	fChain->SetBranchAddress("rhoFactor", &rhoFactor, &b_rhoFactor);
	fChain->SetBranchAddress("ptHat", &ptHat, &b_ptHat);
	fChain->SetBranchAddress("nPUVertices", &nPUVertices, &b_nPUVertices);
	fChain->SetBranchAddress("nPUVerticesTrue", &nPUVerticesTrue, &b_nPUVerticesTrue);
	fChain->SetBranchAddress("triggerStatus", &triggerStatus, &b_triggerStatus);
	fChain->SetBranchAddress("hltPrescale",&hltPrescale, &b_hltPrescale);

  //jobTree->SetBranchAddress("nEvents", &unskimmedEvents, &b_unskimmedEvents);

	// ADD THESE!
	//eventTree->Branch("recoTaus",&recoTaus, 6400, 0);
	//eventTree->Branch("beamSpot", &beamSpot, 6400, 0);
	//eventTree->Branch("bunchCross",&bunchCross, "bunchCross/i");
	//eventTree->Branch("qScale", &qScale, "qScale/F");
	//eventTree->Branch("evtWeight", &evtWeight, "evtWeight/F");

  ///////////////
  // Cuts Init //
  ///////////////


  vetoElID.cutName =                     "vetoElID";
  vetoElID.dEtaIn[0] =                   0.007;
  vetoElID.dPhiIn[0] =                   0.8;
  vetoElID.sigmaIetaIeta[0] =            0.01;
  vetoElID.HadOverEm[0] =                0.15;
  vetoElID.dxy[0] =                      0.04;
  vetoElID.dz[0] =                       0.2;
  vetoElID.fabsEPDiff[0] =               99999;
  vetoElID.ConversionMissHits[0] =       99999;
  vetoElID.PassedConversionProb[0] =     1;

  vetoElID.dEtaIn[1] =                   0.01;
  vetoElID.dPhiIn[1] =                   0.7;
  vetoElID.sigmaIetaIeta[1] =            0.03;
  vetoElID.HadOverEm[1] =                999999;
  vetoElID.dxy[1] =                      0.04;
  vetoElID.dz[1] =                       0.2;
  vetoElID.fabsEPDiff[1] =               999999;
  vetoElID.ConversionMissHits[1] =       999999;
  vetoElID.PassedConversionProb[1] =     1;

  looseElID.cutName =                     "looseElID";
  looseElID.dEtaIn[0] =                   0.007;
  looseElID.dPhiIn[0] =                   0.15;
  looseElID.sigmaIetaIeta[0] =            0.01;
  looseElID.HadOverEm[0] =                0.12;
  looseElID.dxy[0] =                      0.02;
  looseElID.dz[0] =                       0.2;
  looseElID.fabsEPDiff[0] =               0.05;
  looseElID.ConversionMissHits[0] =       1;
  looseElID.PassedConversionProb[0] =     1;

  looseElID.dEtaIn[1] =                   0.009;
  looseElID.dPhiIn[1] =                   0.10;
  looseElID.sigmaIetaIeta[1] =            0.03;
  looseElID.HadOverEm[1] =                0.10;
  looseElID.dxy[1] =                      0.02;
  looseElID.dz[1] =                       0.2;
  looseElID.fabsEPDiff[1] =               0.05;
  looseElID.ConversionMissHits[1] =       1;
  looseElID.PassedConversionProb[1] =     1;
  
  looseElIso.cutName =                    "looseElIso";
  looseElIso.chIso04 =                    99999;
  looseElIso.nhIso04 =                    99999;
  looseElIso.phIso04 =                    99999;
  looseElIso.relCombIso04 =               0.4;

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

  loosePhID.cutName =                     "loosePhID";
  loosePhID.PassedEleSafeVeto[0] =        1;
  loosePhID.HadOverEm[0] =                0.05;
  loosePhID.sigmaIetaIeta[0] =            0.012;

  loosePhID.PassedEleSafeVeto[1] =        1;
  loosePhID.HadOverEm[1] =                0.05;
  loosePhID.sigmaIetaIeta[1] =            0.034;

  loosePhIso.cutName =                    "loosePhIso";
  loosePhIso.chIso03[0] =                 2.6;
  loosePhIso.nhIso03[0] =                 3.5;
  loosePhIso.phIso03[0] =                 1.3;

  loosePhIso.chIso03[1] =                 2.3;
  loosePhIso.nhIso03[1] =                 2.9;
  loosePhIso.phIso03[1] =                 99999;

  mediumPhID.cutName =                     "mediumPhID";
  mediumPhID.PassedEleSafeVeto[0] =        1;
  mediumPhID.HadOverEm[0] =                0.05;
  mediumPhID.sigmaIetaIeta[0] =            0.011;

  mediumPhID.PassedEleSafeVeto[1] =        1;
  mediumPhID.HadOverEm[1] =                0.05;
  mediumPhID.sigmaIetaIeta[1] =            0.033;

  mediumPhIso.cutName =                    "mediumPhIso";
  mediumPhIso.chIso03[0] =                 1.5;
  mediumPhIso.nhIso03[0] =                 1.0;
  mediumPhIso.phIso03[0] =                 0.7;

  mediumPhIso.chIso03[1] =                 1.2;
  mediumPhIso.nhIso03[1] =                 1.5;
  mediumPhIso.phIso03[1] =                 1.0;


  ///////////////////////
  // ZGAngles MVA init //
  ///////////////////////

  // TMVA weights directory
  mvaInits.weightsDir = "/uscms_data/d2/bpollack/CMSSW_4_4_2/src/MVACodes/mvaExamples_brian/weights/";

  // here we will use only BDTG... but will keep the structure 
  mvaInits.discrMethodName[0] = "MLPBNN";
  mvaInits.discrMethodName[1] = "BDTG";
  mvaInits.discrMethodName[2] = "BDT";

  mvaInits.discrSampleName = "MVA_ZJets";


  mvaInits.mvaHiggsMassPoint[0] = 125;

  mvaInits.bdtCut[0] = 0.75;

}

bool higgsAnalyzerV2::Notify()
{
  fileCount+= 1;
  file0 = thisTree->GetCurrentFile();
  h1_numOfEvents = (TH1F*) file0->Get("ntupleProducer/numOfEvents");
  unskimmedEvents = h1_numOfEvents->GetBinContent(1);
  cout<<"THIS IS FILE NUMBER: "<<fileCount<<" and it has this many events: "<<unskimmedEvents<<endl;
  unskimmedEventsTotal += unskimmedEvents;
	return kTRUE;
}

#endif 

