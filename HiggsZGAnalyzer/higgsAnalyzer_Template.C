#define higgsAnalyzer_cxx

#include "higgsAnalyzer.h"

using namespace std;

/////////////////////////////
//Specify parameters here. //
/////////////////////////////

string  selection      = "SELECTION";
string  period         = "PERIOD";
int     JC_LVL         = 0;
string  abcd           = "ABCD";
string  suffix         = "SUFFIX";


/////////////////
//Analysis cuts//
/////////////////

//nominal selection
//float jetPtCut[] = {30,30}, muPtCut[] = {20,10}, elePtCut[] = {20,10}, gammaPtCut[] = {15./110.,15.}, zMassCut[] = {50,99999.}, metCut[] = {-9999.,9999.},
//  metByQtCut[] = {-999.,999.}, bJetVeto = 9999, qtCut = 0., nJets[] = {0,999}, MassZG[] = {100.,180.}, SumEt = 99999., AngCut = 999,
//  M12Cut = 999999999999999999999, drCut = 0.4, MZpMZGcut = 185, VBFdeltaEta = 3.5, VBFdeltaPhi = 2.4, diJetMassCut = 500, Zeppenfeld = 2.5, R9Cut = 0.94; 

float jetPtCut[] = {30,30}, muPtCut[] = {20,10}, elePtCut[] = {20,10}, gammaPtCut[] = {15./110.,15.}, zMassCut[] = {50,99999.}, metCut[] = {-9999.,9999.},
  metByQtCut[] = {-999.,999.}, bJetVeto = 9999, qtCut = 0., nJets[] = {0,999}, MassZG[] = {100.,190.}, SumEt = 99999., AngCut = 999,
  M12Cut = 999999999999999999999, drCut = 0.4, MZpMZGcut = 185, VBFdeltaEta = 3.5, VBFdeltaPhi = 2.4, diJetMassCut = 500, Zeppenfeld = 2.5, R9Cut = 0.94; 


bool    VBFcuts         = false;

bool    DYGammaVeto    = true;
bool    customPhotoID  = false;
bool    spikeVeto      = true;

bool    R9switch       = false;

bool    doEleMVA       = true;

bool    doLooseMuIso   = true;
//#define USE_MVA


///// debugging dumps /////
bool dumps = false;
bool dataDumps = false;
int EVENTNUMBER = -999;

//// energy corrections ////
bool engCor = true;
bool doR9Cor = true;
bool doEleReg  = true;

//// Scale Factors ////
bool doScaleFactors = true;
bool doLumiXS= false;

///////////////////////////
//Resources for weighting//
///////////////////////////

// Do something about these: should just have one sort condition function
bool P4SortCondition(const TLorentzVector& p1, const TLorentzVector& p2) {return (p1.Pt() > p2.Pt());} 
bool MuonSortCondition(const TCMuon& m1, const TCMuon& m2) {return (m1.Pt() > m2.Pt());}
bool ElectronSortCondition(const TCElectron& e1, const TCElectron& e2) {return (e1.Pt() > e2.Pt());}
bool PhotonSortCondition(const TCPhoton& g1, const TCPhoton& g2) {return (g1.Pt() > g2.Pt());}
bool JetSortCondition(const TCJet& j1, const TCJet& j2) {return (j1.Pt() > j2.Pt());}
bool GenSortCondition(const TCGenParticle& ge1, const TCGenParticle& ge2) {return (ge1.Pt() > ge2.Pt());}
bool VertexSortCondition(const TCPrimaryVtx& pv1, const TCPrimaryVtx& pv2) {return (pv1.SumPt2Trks() > pv2.SumPt2Trks());}


// stuff for BDT selection ----------------------------------------------------------------------------
// inclusions made with the 'define'
// to isolate it from the rest of the code


#ifdef USE_MVA

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

// TMVA weights directory
TString weightsDir = "/uscms_data/d2/bpollack/CMSSW_4_4_2/src/MVACodes/mvaExamples_brian/weights/";

#define N_DISCR_METHODS 3
// here we will use only BDTG... but will keep the structure 
enum DISCR_METHOD {
  MLPBNN, BDTG, BDT
};
TString discrMethodName[3] = {
  "MLPBNN", "BDTG", "BDT"
};

enum DISCR_TYPE {
  D_ZJets
};
TString discrSampleName = "MVA_ZJets";


#define N_HIGGS_MASSES 1
Int_t mvaHiggsMassPoint[N_HIGGS_MASSES] = {125};

//      [mass point]
Float_t bdtCut[N_HIGGS_MASSES] = {0.75};
//Float_t bdtCut[N_HIGGS_MASSES] = {0.1};
// Float_t bdtCut[N_HIGGS_MASSES] = {0.5}; //NEURAL NET

//      [sample discr]
TMVA::Reader* tmvaReader;

// variables for the mva tree 

Float_t _diffZGvector, _threeBodyPt, _GPt, _cosZ, _diffPlaneMVA, _vtxVariable, _threeBodyMass, _dr1, _dr2, _M12;



#endif
/// end of MVA related stuff


void higgsAnalyzer::Begin(TTree * tree)
{

  for (int i =0; i<100; i++){
    nEvents[i] = 0;
  }
  for (int i =0; i<2; i++){
    genAccept[i] = 0;
  }

  // Get trigger names from jobTree
  vector<string>* triggerNames = 0;
  TFile   *inFile         = tree->GetCurrentFile();
  TTree   *jobTree        = (TTree*)inFile->Get("ntupleProducer/jobTree");

  jobTree->SetBranchAddress("triggerNames", &triggerNames);
  jobTree->GetEntry();

  // Initialize utilities and selectors here //
  weighter        = new WeightUtils(suffix, period, abcd, selection, isRealData);
  triggerSelector = new TriggerSelector(selection, period, *triggerNames);
  rmcor2011       = new rochcor_2011(229);
  rmcor2012       = new rochcor2012(229);
  rEl             = new TRandom3(1234);
  rMuRun          = new TRandom3(187);
  phoCorrector    = new zgamma::PhosphorCorrectionFunctor("plugins/PHOSPHOR_NUMBERS_EXPFIT_ERRORS.txt", true);

  // Random numbers! //
  rnGenerator = new TRandom3();

  TString option = GetOption();
  TH1::SetDefaultSumw2(kTRUE);
  TH2::SetDefaultSumw2(kTRUE);

  histoFile = new TFile("higgsHistograms_DATANAME_COUNT.root", "RECREATE");
  trainingFile = new TFile("higgsTraining_DATANAME_COUNT.root", "RECREATE");
  sampleFile = new TFile("higgsSample_DATANAME_COUNT.root", "RECREATE");
  higgsFile = new TFile("higgsFile_DATANAME_COUNT.root", "RECREATE");
  eleSmearFile = new TFile("eleSmearFile_DATANAME_COUNT.root", "RECREATE");
  eleIDISOFile = new TFile("eleIDISOFile_DATANAME_COUNT.root", "RECREATE");
  m_llgFile = new TFile("m_llgFile_DATANAME_COUNT.root","RECREATE");

  trainingFile->cd();
  trainingChain = new TTree("varMVA","hey everyone it's the training tree");
  sampleFile->cd();
  sampleChain   = new TTree("varMVA","ZOMG the sample tree!");

  m_llgFile->cd();
  m_llgChain = new TTree("m_llg_SUFFIX","three body mass values");

  initializeEfficiencyWeights( "otherHistos/elsf2011.root", "otherHistos/elsf2012.root");
  eleSmearFile->cd();
  hmEleSmear = new HistManager(eleSmearFile);
  eleIDISOFile->cd();
  hmEleIDISO = new HistManager(eleIDISOFile);

  higgsFile->cd();
  hmHiggs   = new HistManager(higgsFile);
  histoFile->cd();
  hm        = new HistManager(histoFile);
  histoFile->mkdir("Misc", "Misc");
  histoFile->mkdir("Lepton", "Lepton");
  histoFile->mkdir("Jet", "Jet");
  histoFile->mkdir("TESTS", "TESTS");
  histoFile->mkdir("ZGamma", "ZGamma");
  histoFile->mkdir("pT-Eta-Phi", "pT-Eta-Phi");
  histoFile->mkdir("Vtx", "Vtx");
  histoFile->mkdir("PreSel", "PreSel");
  histoFile->mkdir("PhotonIsoDR04", "PhotonIsoDR04");
  histoFile->mkdir("MuonIso", "MuonIso");
  histoFile->mkdir("GenLvl", "GenLvl");
  histoFile->mkdir("GenPlots", "GenPlots");
  histoFile->mkdir("Jet-Photons", "Jet-Photons");
  histoFile->mkdir("CAT1", "CAT1");
  histoFile->mkdir("CAT2", "CAT2");
  histoFile->mkdir("CAT3", "CAT3");
  histoFile->mkdir("CAT4", "CAT4");
  histoFile->mkdir("PreGen", "PreGen");
  histoFile->mkdir("PostGen", "PostGen");
  histoFile->mkdir("ZGAngles", "ZGAngles");

  if (period.find("2012") != string::npos){
    //combined rel ISO, 2012 Data, 0.5 GeV
    EAMu[0] =   0.674; //         eta < 1.0
    EAMu[1] =   0.565; // 1.0   < eta < 1.5
    EAMu[2] =   0.442; // 1.5   < eta < 2.0
    EAMu[3] =   0.515; // 2.0   < eta < 2.2
    EAMu[4] =   0.821; // 2.2   < eta < 2.3
    EAMu[5] =   0.660; // 2.3   < eta < 2.4

    EAEle[0] =   0.208; //         eta < 1.0
    EAEle[1] =   0.209; // 1.0   < eta < 1.5
    EAEle[2] =   0.115; // 1.5   < eta < 2.0
    EAEle[3] =   0.143; // 2.0   < eta < 2.2
    EAEle[4] =   0.183; // 2.2   < eta < 2.3
    EAEle[5] =   0.194; // 2.3   < eta < 2.4
    EAEle[6] =   0.261; // 2.4   < eta


  }else if (period.find("2011") != string::npos){
    //combined rel ISO, 2011 Data, 0.5 GeV
    EAMu[0] =   0.132; //         eta < 1.0
    EAMu[1] =   0.120; // 1.0   < eta < 1.5
    EAMu[2] =   0.114; // 1.5   < eta < 2.0
    EAMu[3] =   0.139; // 2.0   < eta < 2.2
    EAMu[4] =   0.168; // 2.2   < eta < 2.3
    EAMu[5] =   0.189; // 2.3   < eta < 2.4


  }else{
    Abort("period != 2011 OR 2012, figure your shit out");
  }

  // ch      nh       ph
  float EAPhoTemp[7][3] = {
    {0.012,  0.030,   0.148}, //         eta < 1.0
    {0.010,  0.057,   0.130}, // 1.0   < eta < 1.479
    {0.014,  0.039,   0.112}, // 1.479 < eta < 2.0
    {0.012,  0.015,   0.216}, // 2.0   < eta < 2.2
    {0.016,  0.024,   0.262}, // 2.2   < eta < 2.3
    {0.020,  0.039,   0.260}, // 2.3   < eta < 2.4
    {0.012,  0.072,   0.266}  // 2.4   < eta 
  };

  for (unsigned int i =0; i<7; i++){
    for (unsigned int j =0; j<3; j++){
      EAPho[i][j] = EAPhoTemp[i][j];
    }
  }

  muonDump = false;
  electronDump = false;
  if (dumps){
    if (selection == "mumuGamma") muonDump = true;
    if (selection == "eeGamma") electronDump= true;
  }

  if (electronDump){
    elDump2.open("dumps/electronDump2_SUFFIX.txt");
    if (!elDump2.good()) cout << "ERROR: can't open file for writing." << endl;
    elDump2<<"run"             <<" "<<"event"       <<" "<<"pt"
      <<" "<<"eta"             <<" "<<"dEtaIn"      <<" "<<"dPhiIn"
      <<" "<<"sigmaIetaIeta"   <<" "<<"HoE"         <<" "<<"ooemoop"
      <<" "<<"conPass"         <<" "<<"mHits"       <<" "<<"d0"
      <<" "<<"dz"              <<" "<<"chIso"       <<" "<<"nhIso"
      <<" "<<"phIso"           <<" "<<"combIso"     <<" "<<"rho"
      <<" "<<"EA"              <<" "<<"passID"      <<" "<<"passIso"
      <<endl<<endl;
    phDump1.open("dumps/photonDump1El_SUFFIX.txt");
    if (!phDump1.good()) cout << "ERROR: can't open file for writing." << endl;
    phDump2.open("dumps/photonDump2El_SUFFIX.txt");
    if (!phDump2.good()) cout << "ERROR: can't open file for writing." << endl;

    elDumpFinal.open("dumps/electronDumpFinal_SUFFIX.txt");
    if (!elDumpFinal.good()) cout << "ERROR: can't open file for writing." << endl;
    elDumpFinal<<"run"             <<" "<<"event"       <<" "<<"pt"
      <<" "<<"eta"             <<" "<<"dEtaIn"      <<" "<<"dPhiIn"
      <<" "<<"sigmaIetaIeta"   <<" "<<"HoE"         <<" "<<"ooemoop"
      <<" "<<"conPass"         <<" "<<"mHits"       <<" "<<"d0"
      <<" "<<"dz"              <<" "<<"chIso"       <<" "<<"nhIso"
      <<" "<<"phIso"           <<" "<<"combIso"     <<" "<<"rho"
      <<" "<<"EA"              <<" "<<"passID"      <<" "<<"passIso"
      <<endl<<endl;

    if(doEleMVA){
      elDumpMVA.open("dumps/electronDumpMVA_SUFFIX.txt");
      if (!elDumpMVA.good()) cout << "ERROR: can't open file for writing." << endl;
    }
  }
  if (muonDump){
    muDump1.open("dumps/muonDump1_SUFFIX.txt");
    if (!muDump1.good()) cout << "ERROR: can't open file for writing." << endl;
    muDump1<<"run"       <<" "<<"event"       <<" "<<"pt"
      <<" "<<"eta"       <<" "<<"isGlobal"    <<" "<<"isPf"
      <<" "<<"rChi2"     <<" "<<"nHits"       <<" "<<"nMatch"
      <<" "<<"d0"        <<" "<<"dz"          <<" "<<"nPixel"
      <<" "<<"nLayer"    <<" "<<"chIso"       <<" "<<"nhIso"
      <<" "<<"phIso"     <<" "<<"combIso"     <<" "<<"dB"
      <<" "<<"passID"    <<" "<<"passIso"
      <<endl<<endl;
    phDump1.open("dumps/photonDump1Mu_SUFFIX.txt");
    if (!phDump1.good()) cout << "ERROR: can't open file for writing." << endl;
    phDump2.open("dumps/photonDump2Mu_SUFFIX.txt");
    if (!phDump2.good()) cout << "ERROR: can't open file for writing." << endl;

    muDumpFinal.open("dumps/muonDumpFinal_SUFFIX.txt");
    if (!muDumpFinal.good()) cout << "ERROR: can't open file for writing." << endl;
    muDumpFinal<<"run"       <<" "<<"event"       <<" "<<"pt"
      <<" "<<"eta"       <<" "<<"isGlobal"    <<" "<<"isPf"
      <<" "<<"rChi2"     <<" "<<"nHits"       <<" "<<"nMatch"
      <<" "<<"d0"        <<" "<<"dz"          <<" "<<"nPixel"
      <<" "<<"nLayer"    <<" "<<"chIso"       <<" "<<"nhIso"
      <<" "<<"phIso"     <<" "<<"combIso"     <<" "<<"dB"
      <<" "<<"passID"    <<" "<<"passIso"
      <<endl<<endl;
  }

  if (dumps){
    phDump1<<"run"             <<" "<<"event"       <<" "<<"pt"
      <<" "<<"eta"              <<" "<<"CSElVeto"    <<" "<<"HoE"
      <<" "<<"sigmaIetaIeta"    <<" "<<"chIso"       <<" "<<"nhIso"
      <<" "<<"phIso"            <<" "<<"rho"         <<" "<<"chEA"
      <<" "<<"nhEA"             <<" "<<"phEA"        <<" "<<"passID"
      <<" "<<"passIso"
      <<endl<<endl;

    phDump2<<"run"             <<" "<<"event"       <<" "<<"pt"
      <<" "<<"eta"              <<" "<<"CSElVeto"    <<" "<<"HoE"
      <<" "<<"sigmaIetaIeta"    <<" "<<"chIso"       <<" "<<"nhIso"
      <<" "<<"phIso"            <<" "<<"rho"         <<" "<<"chEA"
      <<" "<<"nhEA"             <<" "<<"phEA"        <<" "<<"passID"
      <<" "<<"passIso"          <<" "<<"dr1"         <<" "<<"dr2"
      <<" "<<"pt/Mllg"
      <<endl<<endl;

    finalDump.open("dumps/finalDump_SELECTION_PERIOD_SUFFIX.txt");
  }

  if (dataDumps && suffix == "DATA"){
    dataDump.open("dataDump_DATANAME_COUNT.txt");
    if (!dataDump.good()) cout << "ERROR: can't open file for writing." << endl;
  }


  sampleChain->Branch("diffZGscalar",&diffZGscalar,"diffZGscalar/F");
  trainingChain->Branch("diffZGscalar",&diffZGscalar,"diffZGscalar/F");
  sampleChain->Branch("diffZGvector",&diffZGvector,"diffZGvector/F");
  trainingChain->Branch("diffZGvector",&diffZGvector,"diffZGvector/F");
  sampleChain->Branch("threeBodyMass",&threeBodyMass,"threeBodyMass/F");
  trainingChain->Branch("threeBodyMass",&threeBodyMass,"threeBodyMass/F");
  sampleChain->Branch("threeBodyPt",&threeBodyPt,"threeBodyPt/F");
  trainingChain->Branch("threeBodyPt",&threeBodyPt,"threeBodyPt/F");
  sampleChain->Branch("divPt",&divPt,"divPt/F");
  trainingChain->Branch("divPt",&divPt,"divPt/F");
  sampleChain->Branch("cosZ",&cosZ,"cosZ/F");
  trainingChain->Branch("cosZ",&cosZ,"cosZ/F");
  sampleChain->Branch("cosG",&cosG,"cosG/F");
  trainingChain->Branch("cosG",&cosG,"cosG/F");
  sampleChain->Branch("METdivQt",&METdivQt,"METdivQt/F");
  trainingChain->Branch("METdivQt",&METdivQt,"METdivQt/F");
  sampleChain->Branch("GPt",&GPt,"GPt/F");
  trainingChain->Branch("GPt",&GPt,"GPt/F");
  sampleChain->Branch("ZPt",&ZPt,"ZPt/F");
  trainingChain->Branch("ZPt",&ZPt,"ZPt/F");
  sampleChain->Branch("DPhi",&DPhi,"DPhi/F");
  trainingChain->Branch("DPhi",&DPhi,"DPhi/F");
  sampleChain->Branch("diffPlaneMVA",&diffPlaneMVA,"diffPlaneMVA/F");
  trainingChain->Branch("diffPlaneMVA",&diffPlaneMVA,"diffPlaneMVA/F");
  sampleChain->Branch("vtxVariable",&vtxVariable,"vtxVariable/F");
  trainingChain->Branch("vtxVariable",&vtxVariable,"vtxVariable/F");
  sampleChain->Branch("dr1",&dr1,"dr1/F");
  trainingChain->Branch("dr1",&dr1,"dr1/F");
  sampleChain->Branch("dr2",&dr2,"dr2/F");
  trainingChain->Branch("dr2",&dr2,"dr2/F");
  sampleChain->Branch("M12",&M12,"M12/F");
  trainingChain->Branch("M12",&M12,"M12/F");
  sampleChain->Branch("scaleFactor",&scaleFactor,"scaleFactor/F");
  trainingChain->Branch("scaleFactor",&scaleFactor,"scaleFactor/F");

  m_llgChain->Branch("m_llg_SUFFIX", &m_llg, "m_llg/F");
  m_llgChain->Branch("m_llgCAT1_SUFFIX", &m_llgCAT1, "m_llgCAT1/F");
  m_llgChain->Branch("m_llgCAT2_SUFFIX", &m_llgCAT2, "m_llgCAT2/F");
  m_llgChain->Branch("m_llgCAT3_SUFFIX", &m_llgCAT3, "m_llgCAT3/F");
  m_llgChain->Branch("m_llgCAT4_SUFFIX", &m_llgCAT4, "m_llgCAT4/F");
  m_llgChain->Branch("unBinnedWeight_SUFFIX", &unBinnedWeight, "unBinnedWeight/F");
  m_llgChain->Branch("unBinnedLumiXS_SUFFIX", &unBinnedLumiXS, "unBinnedLumiXS/F");

  // -----------------  MVA stuff -----------------------------------------------------------
#ifdef USE_MVA

  // set up the TMVA readers (the input variables have been set already as global)

  tmvaReader = new TMVA::Reader("!Color:!Silent");

  //add  variables... some are exclusive to particular sample/jet multi discriminators
  tmvaReader->AddVariable("diffZGvector", &_diffZGvector);
  //tmvaReader->AddVariable("threeBodyMass", &_threeBodyMass);
  tmvaReader->AddVariable("threeBodyPt", &_threeBodyPt);
  tmvaReader->AddVariable("GPt", &_GPt);
  tmvaReader->AddVariable("cosZ", &_cosZ);
  tmvaReader->AddVariable("diffPlaneMVA", &_diffPlaneMVA);
  tmvaReader->AddVariable("vtxVariable", &_vtxVariable);
  tmvaReader->AddVariable("dr1", &_dr1);
  tmvaReader->AddVariable("dr2", &_dr2);
  tmvaReader->AddVariable("M12", &_M12);

  // Book the methods
  // for testing we will set only the BDT and hotwire this loop

  int discr = BDTG;
  //int discr = MLPBNN;

  for (int mh = 0; mh < N_HIGGS_MASSES; ++mh) {

    TString label = TString::Format("%s_%s_MVA_HZG%i", discrMethodName[discr].Data(), discrSampleName.Data(),
        mvaHiggsMassPoint[mh]);

    //discr_MVA_ZJets_MVA_HZG125___BDTG.weights.xml
    TString weightFile = TString::Format("%sdiscr_%s_MVA_HZG%i___%s.weights.xml",
        weightsDir.Data(), discrSampleName.Data(), mvaHiggsMassPoint[mh], discrMethodName[discr].Data());

    tmvaReader->BookMVA(label.Data(), weightFile.Data());

  }


#endif
  // ------------------ End of MVA stuff --------------------------------------------------------------------------


  cout<<"load and initialize MVA"<<endl;
  //ElectronMVA Selection
  std::vector<std::string> myManualCatWeigthsTrig;

  myManualCatWeigthsTrig.push_back("plugins/MVAWeights/Electrons_BDTG_TrigV0_Cat1.weights.xml");
  myManualCatWeigthsTrig.push_back("plugins/MVAWeights/Electrons_BDTG_TrigV0_Cat2.weights.xml");
  myManualCatWeigthsTrig.push_back("plugins/MVAWeights/Electrons_BDTG_TrigV0_Cat3.weights.xml");
  myManualCatWeigthsTrig.push_back("plugins/MVAWeights/Electrons_BDTG_TrigV0_Cat4.weights.xml");
  myManualCatWeigthsTrig.push_back("plugins/MVAWeights/Electrons_BDTG_TrigV0_Cat5.weights.xml");
  myManualCatWeigthsTrig.push_back("plugins/MVAWeights/Electrons_BDTG_TrigV0_Cat6.weights.xml");

  myMVATrig = new EGammaMvaEleEstimator();
  myMVATrig->initialize("BDT",
      EGammaMvaEleEstimator::kTrig,
      true,
      myManualCatWeigthsTrig);

  cout<<"MVA Loaded"<<endl;
}


Bool_t higgsAnalyzer::Process(Long64_t entry)
{

  GetEntry(entry,1);
  if(eventNumber == EVENTNUMBER) cout<<EVENTNUMBER<<endl;

  // 2011 bad electron run veto
  if(selection == "eeGamma" && isRealData && period.find("2011") != string::npos){
    if( (runNumber == 171050 && (lumiSection >= 47 && lumiSection <= 92 )) ||
        (runNumber == 171156 && (lumiSection >= 42 && lumiSection <= 211)) ||
        (runNumber == 171219 && (lumiSection >= 48 && lumiSection <= 163)) ||
        (runNumber == 171274 && (lumiSection >= 88 && lumiSection <= 148)) ||
        (runNumber == 171282 && (lumiSection >= 1  && lumiSection <= 172)) ||
        (runNumber == 171315 && (lumiSection >= 53 && lumiSection <= 226)) ||
        (runNumber == 171369 && (lumiSection >= 42 && lumiSection <= 162)) ||
        (runNumber == 171446 && (lumiSection >= 58 && lumiSection <= 374)) ||
        (runNumber == 171484 && (lumiSection >= 61 && lumiSection <= 358)) ||
        (runNumber == 171578 && (lumiSection >= 47 && lumiSection <= 347))
      ) return kTRUE;
  }
    

  hm->fill1DHist(1,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  hm->fill1DHist(1,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  ++nEvents[0];

  if (nEvents[0] == 1) weighter->SetDataBit(isRealData);

  if (nEvents[0] % (int)1e5 == 0) cout<<nEvents[17]<<" events passed of "<<nEvents[0]<<" checked!"<<endl;

  m_llg = m_llgCAT1 = m_llgCAT2 = m_llgCAT3 = m_llgCAT4 = -1;
  unBinnedWeight = unBinnedLumiXS = 1;

  ///////////////////////////
  // Electron Preselection //
  ///////////////////////////

  /*
  if (selection == "eeGamma"){
    if (recoElectrons->GetSize() < 2) return kTRUE;  

    bool leadingEl = false;
    int leadingPos = -1;
    for (int i = 0; i <  recoElectrons->GetSize(); ++i) {
      TCElectron* thisElec = (TCElectron*) recoElectrons->At(i);    
      if (thisElec->Pt() > 20){
        leadingEl = true;
        leadingPos = i;
        break;
      }
    }
    if (!leadingEl) return kTRUE;

    bool trailingEl = false;
    for (int i = 0; i <  recoElectrons->GetSize(); ++i) {
      TCElectron* thisElec = (TCElectron*) recoElectrons->At(i);    
      if ((i != leadingPos) && (thisElec->Pt() > 10)){
        trailingEl = true;
        break;
      }
    }
    if (!trailingEl) return kTRUE;
  }
*/

  hm->fill1DHist(2,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  hm->fill1DHist(2,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  if (!isRealData) hm->fill1DHist(nPUVerticesTrue,"h1_simVertexMultTrueUltraFine_SUFFIX", "Multiplicity of simulated vertices true", 500, 0, 100,1,"Misc");
  ++nEvents[1];



  //////////////////
  //Trigger status//
  //////////////////

  for(int i = 0; i < 64; ++i) {
    unsigned long iHLT = 0x0; 
    iHLT = 0x01 << i;  
    if ((triggerStatus & iHLT) == iHLT) hm->fill1DHist(i+1,"h1_triggerStatus_SUFFIX", "Triggers", 64, 0.5, 64.5,1,"Misc");  
  } 

  bool triggerPass   = triggerSelector->SelectTrigger(triggerStatus, hltPrescale);
  //cout<<"triggerStatus: "<<triggerStatus<<" hltPrescale: "<<hltPrescale<<" triggerPass: "<<triggerPass<<endl;
  //int  eventPrescale = triggerSelector->GetEventPrescale();
  //cout<<eventPrescale<<endl;

  /*
  cout<<"new trig event"<<endl;
  for (int i = 0; i <  triggerObjects->GetSize(); ++i) {
    TCTriggerObject* thisTrigObj = (TCTriggerObject*) triggerObjects->At(i);    
    cout<<" HLT: "<<thisTrigObj->GetHLTName()<<endl;
    cout<<" module: "<<thisTrigObj->GetModuleName()<<endl;
    cout<<" id: "<<thisTrigObj->GetId()<<endl;
    thisTrigObj->Print();
    cout<<endl;
  }
  */
  if (!triggerPass) return kTRUE;
  hm->fill1DHist(3,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  hm->fill1DHist(3,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  ++nEvents[2];




  ////////////////////////////
  //Check the event vertices//
  ////////////////////////////


  if (!isRealData) hm->fill1DHist(nPUVertices,"h1_simVertexMultStoyan_SUFFIX", "Multiplicity of simulated vertices", 60, 0, 60, 1, "Misc");
  if (!isRealData) hm->fill1DHist(nPUVertices,"h1_simVertexMultPoter_SUFFIX", "Multiplicity of simulated vertices", 100, 0, 100, 1, "Misc");
  if (!isRealData) hm->fill1DHist(nPUVertices,"h1_simVertexMult_SUFFIX", "Multiplicity of simulated vertices", 50, -0.5, 49.5,1,"Misc");
  if (!isRealData) hm->fill1DHist(nPUVerticesTrue,"h1_simVertexMultTrue_SUFFIX", "Multiplicity of simulated vertices true", 50, -0.5, 49.5,1,"Misc");
  if (!isRealData) hm->fill1DHist(nPUVerticesTrue,"h1_simVertexMultTrueFine_SUFFIX", "Multiplicity of simulated vertices true", 200, 0, 100,1,"Misc");


  vector<TVector3> goodVertices;
  vector<TCPrimaryVtx> pvVect;
  for (int i = 0; i < primaryVtx->GetSize(); ++i) {
    TCPrimaryVtx* pVtx = (TCPrimaryVtx*) primaryVtx->At(i);
    if (
        !pVtx->IsFake() 
        && pVtx->NDof() > 4.
        && fabs(pVtx->z()) <= 24.
        && fabs(pVtx->Perp()) <= 2.
       )

      goodVertices.push_back(*pVtx);

  }
  if (goodVertices.size() < 1) return kTRUE;
  hm->fill1DHist(4,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  hm->fill1DHist(4,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  ++nEvents[3];

  pvPosition = new TVector3();
  *pvPosition = goodVertices[0];



  //////////////////
  // Data quality //
  //////////////////

  //if (isRealData && (isDeadEcalCluster || isScraping || isCSCTightHalo)) return kTRUE;
  //if (isRealData && (isNoiseHcal || isDeadEcalCluster || isScraping || isCSCTightHalo)) return kTRUE;
  hm->fill1DHist(5,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  hm->fill1DHist(5,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  ++nEvents[4];

  ///////////////////
  // gen Particles //
  ///////////////////

  vector<TCGenParticle> genElectrons;
  vector<TCGenParticle> genMuons;
  vector<TCGenParticle> genZs;
  vector<TCGenParticle> genWs;
  vector<TCGenParticle> genHs;
  vector<TCGenParticle> genPhotons;
  bool isMuMuGamma = false;
  bool isEEGamma = false;

  for (int i = 0; i < genParticles->GetSize(); ++i) {
    TCGenParticle* thisGen = (TCGenParticle*) genParticles->At(i);    
    //   cout<<thisGen->GetPDGId()<<endl;
    if (abs(thisGen->GetPDGId()) == 11){
      genElectrons.push_back(*thisGen);
      if (abs(thisGen->Mother())==23) isEEGamma = true;
    }else if (abs(thisGen->GetPDGId()) == 13){
      genMuons.push_back(*thisGen);
      if (abs(thisGen->Mother())==23) isMuMuGamma = true;
    }else if (abs(thisGen->GetPDGId()) == 23) genZs.push_back(*thisGen);
    else if (abs(thisGen->GetPDGId()) == 24) genWs.push_back(*thisGen);
    else if (abs(thisGen->GetPDGId()) == 22) genPhotons.push_back(*thisGen);
    else if (abs(thisGen->GetPDGId()) == 25) genHs.push_back(*thisGen);
  }
  ///////// whzh decomposition /////////////////

  if (suffix.find("zh") != string::npos && genWs.size() >  0) return kTRUE;
  if (suffix.find("wh") != string::npos && genWs.size() == 0) return kTRUE;


  ///////// sort gen particles by PT ///////////

  sort(genElectrons.begin(), genElectrons.end(), GenSortCondition);
  sort(genMuons.begin(), genMuons.end(), GenSortCondition);
  sort(genZs.begin(), genZs.end(), GenSortCondition);
  sort(genPhotons.begin(), genPhotons.end(), GenSortCondition);
  sort(genHs.begin(), genHs.end(), GenSortCondition);

  ZGLabVectors genLevelInputs;
  ZGAngles     genLevelOutputs;
  TCGenParticle gPhot;
  TCGenParticle gMuon1;
  TCGenParticle gMuon2;
  TCGenParticle gElectron1;
  TCGenParticle gElectron2;
  TLorentzVector genAccL1;
  TLorentzVector genAccL2;
  TLorentzVector genAccG;
  TLorentzVector genZ;
  TLorentzVector genZG;
  bool gFirstMu = false;
  bool gBothMu = false;
  bool gFirstEl = false;
  bool gBothEl = false;
  vector<TCGenParticle>::iterator testIt;

  if (genPhotons.size() > 0){
      for (testIt=genPhotons.begin(); testIt<genPhotons.end(); testIt++){
        //cout<<"mother: "<<testIt->Mother()<<"\tstatus: "<<testIt->GetStatus()<<endl;
        if (testIt->Mother() == 25) genAccG = *testIt; break;
      }
  }

  int gCharge = 0;
  if (selection == "mumuGamma" && isMuMuGamma && genMuons.size() > 1 && genPhotons.size() > 0){
    if (genMuons.size() > 1){
      for (testIt=genMuons.begin(); testIt<genMuons.end(); testIt++){
        if(testIt->Mother()== 23 && !gFirstMu){
          gMuon1 = *testIt;
          gFirstMu = true;

          gCharge = testIt->Charge();
        }else if(testIt->Mother()== 23 &&  gFirstMu && testIt->Charge()!=gCharge){
          gMuon2 = *testIt;
          gBothMu = true;
          break;
        }
      }
    }
    if (gMuon1.Charge() == 1){
      genAccL1 = gMuon1;
      genAccL2 = gMuon2;
    }else{
      genAccL1 = gMuon2;
      genAccL2 = gMuon1;
    }
    //genZ = genAccL1+genAccL2;
    //genZG = genAccL1+genAccL2+genAccG;

    if(genAccL1.Pt() > 0.1 && genAccL2.Pt() > 0.1 && genAccG.Pt() > 0.1){
      StandardPlots(genAccL1,genAccL2,genAccG,1,"PreGen", "PreGen");
      genAccG.Print();
      genAccL1.Print();
      genAccL2.Print();
      genLevelInputs.veczg = genAccL1+genAccL2+genAccG; 
      cout<<"higgsMass: "<<genLevelInputs.veczg.M()<<endl;
      genLevelInputs.vecz = genAccL1+genAccL2; 
      genLevelInputs.vecg = genAccG;

      genLevelInputs.veclp = genAccL1;
      genLevelInputs.veclm = genAccL2;

      getZGAngles(genLevelInputs,genLevelOutputs);
      AnglePlots(genLevelOutputs,1);
      //cout<<"costheta_lm: "<<genLevelOutputs.costheta_lm<<"\tcostheta_lp: "<<genLevelOutputs.costheta_lp<<"\tphi: "<<genLevelOutputs.phi<<"\tcosTheta: "<<genLevelOutputs.cosTheta<<"\tcosThetaG: "<<genLevelOutputs.cosThetaG<<endl;
    }
  }else if (selection == "eeGamma" && isEEGamma && genElectrons.size() > 1 && genPhotons.size() > 0){
    if (genElectrons.size() > 1){
      for (testIt=genElectrons.begin(); testIt<genElectrons.end(); testIt++){
        if(testIt->Mother()== 23 && !gFirstEl){
          gElectron1 = *testIt;
          gFirstEl = true;
          gCharge = testIt->Charge();
        }else if(testIt->Mother()== 23 &&  gFirstEl && testIt->Charge()!=gCharge){
          gElectron2 = *testIt;
          gBothEl = true;
          break;
        }
      }
    }
    if (gElectron1.Charge() == 1){
      genAccL1 = gElectron1;
      genAccL2 = gElectron2;
    }else{
      genAccL1 = gElectron1;
      genAccL2 = gElectron2;
    }
    if(genAccL2.Pt() > 1 && genAccL1.Pt() > 1 && genAccG.Pt() > 1){
      StandardPlots(genAccL1,genAccL2,genAccG,1,"PreGen", "PreGen");
      genLevelInputs.veczg = (genAccL1+genAccL2+genAccG);
      genLevelInputs.vecz = (genAccL1+genAccL2);
      genLevelInputs.vecg = genAccG;

      genLevelInputs.veclp = genAccL1;
      genLevelInputs.veclm = genAccL2;

      getZGAngles(genLevelInputs,genLevelOutputs);
      AnglePlots(genLevelOutputs,1);
    }
  }



  ////////// clean the gen photons (fiducial) /////////
  vector<TCGenParticle> genPhotonsEtaClean;

  if (genHs.size()>0)      hm->fill1DHist(genHs[0].M(),"h1_genHM_SUFFIX", "Gen: Higgs Mass; M (GeV); N_{evts}",40,124,126,1,"GenPlots");


  bool isRealPho = false;
  TCGenParticle realPho;
  if (genHs.size()>0 && genPhotons.size()>0){
    for (unsigned int i = 0; i<genPhotons.size(); i++){
      if (genPhotons[i].GetStatus() == 1){
        TLorentzVector GGenBoost(genPhotons[i]);
        hm->fill1DHist(genPhotons[i].Pt(),"h1_realGenPhotonPt_SUFFIX", "Gen: real Photon pT; pT; N_{evts}", 150,0,150,1,"GenPlots");
        hm->fill1DHist(genPhotons[i].E(),"h1_realGenPhotonE_SUFFIX", "Gen: real Photon E; E; N_{evts}", 300,0,300,1,"GenPlots");
        GGenBoost.Boost(-genHs[0].BoostVector());
        hm->fill1DHist(GGenBoost.E(),"h1_realGenPhotonCoME_SUFFIX", "Gen: real Photon (higgs CoM) E; E; N_{evts}", 40,50,70,1,"GenPlots");
        realPho = genPhotons[i];
        isRealPho = true;
        break;
      }

  }
  if (selection == "eeGamma"){
    if (genElectrons.size() > 1 && genHs.size()>0){
      int elcharge = genElectrons[0].Charge();
      for (unsigned int i =1; i<genElectrons.size(); i++){
        if (genElectrons[i].Charge() != elcharge){
          TLorentzVector SGenBoost((genElectrons[0]+genElectrons[i]));
          SGenBoost.Boost(-genHs[0].BoostVector());
          hm->fill1DHist((genElectrons[0]+genElectrons[i]).Pt(),"h1_starGenPhotonPt_SUFFIX", "Gen: star Photon pT; pT; N_{evts}", 150,0,150,1,"GenPlots");
          hm->fill1DHist((genElectrons[0]+genElectrons[i]).E(),"h1_starGenPhotonE_SUFFIX", "Gen: star Photon E; E; N_{evts}", 300,0,300,1,"GenPlots");
          hm->fill1DHist(SGenBoost.E(),"h1_starGenPhotonCoME_SUFFIX", "Gen: star Photon (higgs CoM) E; E; N_{evts}", 40,50,70,1,"GenPlots");
          hm->fill1DHist((genElectrons[0]+genElectrons[i]).M(),"h1_genDiLepMass_SUFFIX", "Gen: diLepton Mass; M (GeV); N_{evts}", 50,-5,5,1,"GenPlots");
          hm->fill1DHist(genElectrons[0].Pt(), "h1_genLeadLepPt_SUFFIX","Gen: leading lepton pt; pT (GeV); N_{evts}", 100,0,100,1,"GenPlots");
          hm->fill1DHist(genElectrons[i].Pt(), "h1_genTrailLepPt_SUFFIX","Gen: trailing lepton pt; pT (GeV); N_{evts}", 100,0,100,1,"GenPlots");
          if (isRealPho){
            hm->fill1DHist((genElectrons[0]+genElectrons[i]+realPho).M(),"h1_gen3BodyMass_SUFFIX", "Gen: 3body Mass; M (GeV); N_{evts}", 40,124,126,1,"GenPlots");
          }
        }
      }
      //cout<<"electron mother: "<<genElectrons[0].Mother()<<" "<<genElectrons[1].Mother()<<" "<<genElectrons[2].Mother()<<" "<<genElectrons[3].Mother()<<endl;
    }
  } else if (selection == "mumuGamma"){
    for (unsigned int i = 0; i<genMuons.size(); i++){
      //if (genMuons[i].Mother() == 22) cout<<"Muon mother is Z "<<endl;
    }
    if (genMuons.size() > 1 && genHs.size() > 0){
      int mucharge = genMuons[0].Charge();
      for (unsigned int i =1; i<genMuons.size(); i++){
        if (genMuons[i].Charge() != mucharge){
          TLorentzVector SGenBoost((genMuons[0]+genMuons[i]));
          SGenBoost.Boost(-genHs[0].BoostVector());
          hm->fill1DHist((genMuons[0]+genMuons[i]).Pt(),"h1_starGenPhotonPt_SUFFIX", "Gen: star Photon pT; pT; N_{evts}", 150,0,150,1,"GenPlots");
          hm->fill1DHist((genMuons[0]+genMuons[i]).E(),"h1_starGenPhotonE_SUFFIX", "Gen: star Photon E; E; N_{evts}", 300,0,300,1,"GenPlots");
          hm->fill1DHist(SGenBoost.E(),"h1_starGenPhotonCoME_SUFFIX", "Gen: star Photon (higgs CoM) E; E; N_{evts}", 40,50,70,1,"GenPlots");
          hm->fill1DHist((genMuons[0]+genMuons[i]).M(),"h1_genDiLepMass_SUFFIX", "Gen: diLepton Mass; M (GeV); N_{evts}", 50,-5,5,1,"GenPlots");
          hm->fill1DHist(genMuons[0].Pt(), "h1_genLeadLepPt_SUFFIX","Gen: leading lepton pt; pT (GeV); N_{evts}", 100,0,100,1,"GenPlots");
          hm->fill1DHist(genMuons[i].Pt(), "h1_genTrailLepPt_SUFFIX","Gen: trailing lepton pt; pT (GeV); N_{evts}", 100,0,100,1,"GenPlots");
          if (isRealPho){
            hm->fill1DHist((genMuons[0]+genMuons[i]+realPho).M(),"h1_gen3BodyMass_SUFFIX", "Gen: 3body Mass; M (GeV); N_{evts}", 40,124,126,1,"GenPlots");
          }
        }
      }
    }
  }
}

  if (genPhotons.size() > 0){
    //cout<<"photon mother ID:"<<endl;
    for (testIt=genPhotons.begin(); testIt<genPhotons.end(); testIt++){
      //cout<<" "<<testIt->Mother()<<" "<<testIt->Grandmother()<<endl;
      if(testIt->Mother()== 22 || testIt->Grandmother() == 25 || testIt->Mother() == 23) gPhot = *testIt;

      ///////////////////////
      // DYJets Gamma Veto //
      ///////////////////////

      if (DYGammaVeto && (suffix.find("ZJets") != string::npos)){
        // if the photon's mother and grandmother is a lepton, kill it
        if ( (abs(testIt->Mother()) == 11 || abs(testIt->Mother()) == 13 || abs(testIt->Mother()) == 15)
           && (abs(testIt->Grandmother()) == 11 || abs(testIt->Grandmother()) == 13 || abs(testIt->Grandmother()) == 15) ) return kTRUE; 
        // if the photon's mother is a lepton, and the grandmother is a Z or W, destroy it
        if ( (abs(testIt->Mother()) == 11 || abs(testIt->Mother()) == 13 || abs(testIt->Mother()) == 15)
           && (abs(testIt->Grandmother()) == 23 || abs(testIt->Grandmother()) == 24) ) return kTRUE;
        // if the photon's mother is a photon, and the grandmother is an electron, raze it to the ground
        if ( abs(testIt->Mother()) == 22 && abs(testIt->Grandmother()) == 11) return kTRUE;
        // if the photon's mother is a gluon (?!) or a quark, eliminate it with prejudice
        if ( abs(testIt->Mother()) == 21 || abs(testIt->Mother()) < 7) return kTRUE;
      }


    }
  }

  ////////////////////
  // gen acceptance //
  ////////////////////

if(!isRealData){
  if(
      (selection == "mumuGamma" && gMuon2.Pt() > 0 && gMuon1.Pt() > 0)
      && gBothMu

      && gMuon1.Pt()               > muPtCut[0]
      && gMuon2.Pt()               > muPtCut[1]
      && fabs(gMuon1.Eta())   < 2.4
      && fabs(gMuon2.Eta())   < 2.4

    ){
    ++genAccept[0]; //cout<<"event passes fuck yeah!"<<endl;
  }
  if(
      (selection == "eeGamma" && gElectron2.Pt() > 0 && gElectron1.Pt() > 0)
      && gBothEl

      && gElectron1.Pt()               > elePtCut[0]
      && gElectron2.Pt()               > elePtCut[1]
      && fabs(gElectron1.Eta())   < 2.4
      && fabs(gElectron2.Eta())   < 2.4

    ){
    ++genAccept[0]; //cout<<"event passes fuck yeah!"<<endl;
  }
  if(
      (selection == "mumuGamma" && gMuon2.Pt() > 0 && gMuon1.Pt() > 0 && gPhot.Pt() > 0)
      && gBothMu

      && gMuon1.Pt()                            > muPtCut[0]
      && gMuon2.Pt()                            > muPtCut[1]
      && fabs(gMuon1.Eta())                     < 2.4
      && fabs(gMuon2.Eta())                     < 2.4

      && (gMuon1+gMuon2).M()                    > zMassCut[0]
      && (gMuon1+gMuon2).M()                    < zMassCut[1]
      && gMuon1.DeltaR(gPhot)                   > 0.4
      && gMuon2.DeltaR(gPhot)                   > 0.4
      && (gMuon1+gMuon2+gPhot).M()              > MassZG[0]
      && (gMuon1+gMuon2+gPhot).M()              < MassZG[1]
      && gPhot.Pt()/((gMuon1+gMuon2+gPhot).M()) > gammaPtCut[0]
      && fabs(gPhot.Eta())                      < 2.5
      && (fabs(gPhot.Eta())                     < 1.4442 || fabs(gPhot.Eta()) > 1.566)

    ){
    ++genAccept[1]; //cout<<"event passes fuck yeah!"<<endl;
    StandardPlots(genAccL1,genAccL2,genAccG,1,"PostGen", "PostGen");
    
  }
  if(
      (selection == "eeGamma" && gElectron2.Pt() > 0 && gElectron1.Pt() > 0 && gPhot.Pt() > 0)
      && gBothEl

      && gElectron1.Pt()                                > elePtCut[0]
      && gElectron2.Pt()                                > elePtCut[1]
      && fabs(gElectron1.Eta())                         < 2.4
      && fabs(gElectron2.Eta())                         < 2.4

      && (gElectron1+gElectron2).M()                    > zMassCut[0]
      && (gElectron1+gElectron2).M()                    < zMassCut[1]
      && gElectron1.DeltaR(gPhot)                       > 0.4
      && gElectron2.DeltaR(gPhot)                       > 0.4
      && (gElectron1+gElectron2+gPhot).M()              > MassZG[0]
      && (gElectron1+gElectron2+gPhot).M()              < MassZG[1]
      && gPhot.Pt()/((gElectron1+gElectron2+gPhot).M()) > gammaPtCut[0]
      && fabs(gPhot.Eta())                              < 2.5
      && (fabs(gPhot.Eta())                             < 1.4442 || fabs(gPhot.Eta()) > 1.566)

    ){
    ++genAccept[1]; //cout<<"event passes fuck yeah!"<<endl;
    StandardPlots(genAccL1,genAccL2,genAccG,1,"PostGen", "PostGen");
  }
}

  //////////////////////
  // Gen Level Angles //
  //////////////////////





  

  /*
     cout<<endl;
     genHs[0].P4().Print();
     (genPhotons[0].P4()+genZs[0].P4()).Print();
     cout<<genMuons.size()<<endl;
     if(genMuons.size()>0) (genPhotons[0].P4()+genMuons[0].P4()+genMuons[1].P4()).Print();
     */
  ///////////////
  // electrons //
  ///////////////

  vector<TLorentzVector> extraLeptons;
  vector<TCElectron> electronsID;
  vector<TCElectron> electronsIDIso;
  vector<TCElectron> electronsIDIsoUnCor;
  TCElectron * cloneElectron;
  //int eleCount = 0;


  for (int i = 0; i <  recoElectrons->GetSize(); ++i) {
    TCElectron* thisElec = (TCElectron*) recoElectrons->At(i);    



    // Section for electron energy/momentum corrections.  NOTE: this will change the pt and thus ID/ISO of electron

    thisElec->SetPtEtaPhiM(thisElec->Pt(),thisElec->Eta(),thisElec->Phi(),0.000511);
    

    // loose cuts on electron for diLepton requirements

    

    if(doEleMVA){
      bool passAll = false;

      if(electronDump) MVADumper(thisElec, myMVATrig, rhoFactor, elDumpMVA);

      if (thisElec->IdMap("preSelPassV1")) electronsID.push_back(*thisElec);			
      double tmpMVAValue = myMVATrig->mvaValue(
          thisElec->IdMap("fbrem"),    
          thisElec->IdMap("kfChi2"),                            
          thisElec->IdMap("kfNLayers"),                            
          thisElec->IdMap("gsfChi2"),                           
          thisElec->IdMap("dEta"),
          thisElec->IdMap("dPhi"),            
          thisElec->IdMap("dEtaAtCalo"),
          thisElec->SigmaIEtaIEta(), 
          thisElec->IdMap("SigmaIPhiIPhi"),
          thisElec->IdMap("SCEtaWidth"),
          thisElec->IdMap("SCPhiWidth"),
          thisElec->IdMap("ome1x5oe5x5"),
          thisElec->IdMap("R9"),
          thisElec->HadOverEm(),
          thisElec->IdMap("EoP"),
          thisElec->IdMap("ooemoopV1"), 
          thisElec->IdMap("eopOut"), 
          thisElec->IdMap("preShowerORaw"), 
          thisElec->IdMap("d0"),
          thisElec->IdMap("ip3d"), 
          thisElec->SCEta(),
          thisElec->Pt(),
          false);                
                                                                  
      /// inner barrel
      if (thisElec->IdMap("preSelPassV1") && thisElec->Pt() > 20 && tmpMVAValue > -0.5 && PassElectronIso(thisElec, looseElIso)){
        passAll = true;
      /// outer barrel
      }else if (thisElec->IdMap("preSelPassV1") && thisElec->Pt() < 20 && tmpMVAValue > -0.90 && PassElectronIso(thisElec, looseElIso)){
        passAll = true;
      }

      if (passAll){
        cloneElectron = thisElec;
        electronsIDIso.push_back(*thisElec);			
        if (EVENTNUMBER == eventNumber) cout<<"regE: "<<thisElec->IdMap("EnergyRegression")<<" regE-p: "<<thisElec->RegressionMomCombP4().E()<<endl;
        if (engCor || doEleReg) electronsIDIsoUnCor.push_back(*cloneElectron);
      }
                                                                  

    }else{
      if(engCor){
        float energyElCor;
        if ( period.find("2011") != string::npos ){
          energyElCor = correctedElectronEnergy( thisElec->E(), thisElec->SCEta(), thisElec->R9(), runNumber, 0, "2011", !isRealData, rEl );
        }else{
          energyElCor = correctedElectronEnergy( thisElec->E(), thisElec->SCEta(), thisElec->R9(), runNumber, 0, "HCP2012", !isRealData, rEl );
        }
        float newPt = sqrt((pow(energyElCor,2)-pow(0.000511,2))/pow(cosh(thisElec->Eta()),2));
        cloneElectron = thisElec;
        thisElec->SetPtEtaPhiM(newPt,thisElec->Eta(),thisElec->Phi(),0.000511);
      }
      if (electronDump) ElectronDump(thisElec, looseElID, looseElIso, elDump2);
      if (PassElectronID(thisElec, looseElID)) electronsID.push_back(*thisElec);			
      if (PassElectronID(thisElec, looseElID) && PassElectronIso(thisElec, looseElIso)){
        electronsIDIso.push_back(*thisElec);			
        if (engCor) electronsIDIsoUnCor.push_back(*cloneElectron);
      }
    } 

    // 3rd lepton veto

  }

  sort(electronsID.begin(), electronsID.end(), ElectronSortCondition);
  sort(electronsIDIso.begin(), electronsIDIso.end(), ElectronSortCondition);


  ///////////
  // muons //
  ///////////

  vector<TCMuon> muonsID;
  vector<TCMuon> muonsIDIso;
  vector<TCMuon> muonsIDIsoUnCor;
  TLorentzVector tmpMuCor;
  TCMuon* cloneMuon;

  /*


     float EAMu[6] = {     //combined rel ISO, 2012 Data, 1.0 GeV
     0.340, //         eta < 1.0
     0.310, // 1.0   < eta < 1.5
     0.315, // 1.5   < eta < 2.0
     0.415, // 2.0   < eta < 2.2
     0.658, // 2.2   < eta < 2.3
     0.405, // 2.3   < eta < 2.4
     };
     */



  for (int i = 0; i < recoMuons->GetSize(); ++ i)
    //for (int i = 0; i < pfMuons->GetSize(); ++ i)
  {
    TCMuon* thisMuon = (TCMuon*) recoMuons->At(i);    
    /*
    if(eventNumber==EVENTNUMBER){
      float ptErrMuTmp = 1.0;
      thisMuon->Print(); cout<<" charge: "<<thisMuon->Charge()<<endl;
      tmpMuCor.SetXYZT(11.398083,-24.497328,38.084204,46.695322);
      rmcor2012->momcor_mc(tmpMuCor,thisMuon->Charge(),0,0,ptErrMuTmp);
      cout<<" cor val: "<<endl;
      tmpMuCor.Print();

    }
    */


    /*
       thisMuon->LoadMap(*mapCuts);
       map<string, float> muCut = thisMuon->MuonMap();
       cout<<"new event muons"<<endl;
       for( map<string,float>::iterator ii=muCut.begin(); ii!=muCut.end(); ++ii)
       {                                                                                                        
       cout << (*ii).first << ": " << (*ii).second << endl;                                                   
       } 
       */

    // Section for muon energy/momentum corrections.  NOTE: this will change the pt and thus ID/ISO of muon

    if(engCor){
      tmpMuCor = *thisMuon;
      if ( period.find("2011") != string::npos ){
        if (isRealData){
          if ( period.find("A") != string::npos ){
            rmcor2011->momcor_data(tmpMuCor,thisMuon->Charge(),0,0);
          }else{
            rmcor2011->momcor_data(tmpMuCor,thisMuon->Charge(),0,1);
          }
        }
        if (!isRealData){
          if (rMuRun->Rndm() < 0.458){
            rmcor2011->momcor_mc(tmpMuCor,thisMuon->Charge(),0,0);
          }else{
            rmcor2011->momcor_mc(tmpMuCor,thisMuon->Charge(),0,1);
          }
        }
      }else{
        float ptErrMu = 1.0;
        if (isRealData){
          if (abcd.find("D") != string::npos ){
            rmcor2012->momcor_data(tmpMuCor,thisMuon->Charge(),1,ptErrMu);
          }else{
            rmcor2012->momcor_data(tmpMuCor,thisMuon->Charge(),0,ptErrMu);
          }
        }
         
        if (!isRealData) rmcor2012->momcor_mc(tmpMuCor,thisMuon->Charge(),0,ptErrMu);
      }
      cloneMuon = thisMuon;
      thisMuon->SetPtEtaPhiM(tmpMuCor.Pt(), tmpMuCor.Eta(), tmpMuCor.Phi(), tmpMuCor.M());
    }


    // tight muon id

    if (doLooseMuIso){
      if(muonDump) MuonDump(thisMuon, tightMuID, looseMuIso, muDump1);
    }else{
      if(muonDump) MuonDump(thisMuon, tightMuID, tightMuIso, muDump1);
    }

    if (PassMuonID(thisMuon, tightMuID)) muonsID.push_back(*thisMuon);

    //tight ID and Iso

    if (doLooseMuIso){
      if (PassMuonID(thisMuon, tightMuID) && PassMuonIso(thisMuon, looseMuIso)){
        muonsIDIso.push_back(*thisMuon);
        if (engCor) muonsIDIsoUnCor.push_back(*cloneMuon);
      }
    }else{
      if (PassMuonID(thisMuon, tightMuID) && PassMuonIso(thisMuon, tightMuIso)){
        muonsIDIso.push_back(*thisMuon);
        if (engCor) muonsIDIsoUnCor.push_back(*cloneMuon);
      }
    }

  }

  sort(muonsID.begin(), muonsID.end(), MuonSortCondition);
  sort(muonsIDIso.begin(), muonsIDIso.end(), MuonSortCondition);
  sort(extraLeptons.begin(), extraLeptons.begin(), P4SortCondition);


  /////////////
  // photons //
  /////////////

  vector<TCPhoton> photonsID; 
  vector<TCPhoton> photonsIDIso; 
  vector<TCPhoton> photonsIDIsoUnCor; 
  TCPhoton* clonePhoton;


  if (selection == "mumuGamma" || selection == "eeGamma") {
    for (Int_t i = 0; i < recoPhotons->GetSize(); ++i) {
      //cout<<endl;
      //cout<<"new photon!!!!!!!"<<endl;
      vector<float> TrkIsoVec;
      TCPhoton* thisPhoton = (TCPhoton*) recoPhotons->At(i);

      if (spikeVeto && (period == "2012A_Jul13" || period == "2012A_Aug06rec" || period == "2012B_Jul13")){
      //  cout<<"SCeta: "<<thisPhoton->SCEta()<<" SCphi: "<<thisPhoton->SCPhi()<<endl;
        float dEta = thisPhoton->SCEta()-(-1.76);
        float dPhi = TVector2::Phi_mpi_pi(thisPhoton->SCPhi()-1.37);
        float spikeDR = sqrt(dEta*dEta+dPhi*dPhi);
        //cout<<"dEta: "<<dEta<<" dPhi: "<<dPhi<<" spikeDr: "<<spikeDR<<endl;
        if (spikeDR < 0.05) continue;
        dEta = thisPhoton->SCEta()-(2.37);
        dPhi = TVector2::Phi_mpi_pi(thisPhoton->SCPhi()-2.69);
        spikeDR = sqrt((dEta*dEta)+(dPhi*dPhi));
        if (spikeDR < 0.05) continue;
      }

      // Section for photon energy/momentum corrections.  NOTE: this will change the pt and thus ID/ISO of photon
      
      if(engCor){
        R9Cor = thisPhoton->R9();
        if (doR9Cor){
          if (suffix != "DATA"){
            if (period == "2011"){
              if (fabs(thisPhoton->SCEta()) < 1.479){
                R9Cor = thisPhoton->R9()*1.0048;
              }else{
                R9Cor = thisPhoton->R9()*1.00492;
              }
            } else if (period == "2012"){
              if (fabs(thisPhoton->SCEta()) < 1.479){
                R9Cor = thisPhoton->R9()*1.0045 + 0.0010;
              }else{
                R9Cor = thisPhoton->R9()*1.0086 - 0.0007;
              }
            }
          }
        }
        clonePhoton = thisPhoton;
        float corrPhoPt = -1;
        int periodNum   = -1;
        if (period.find("2011") != string::npos) periodNum = 2011;
        if (period.find("2012") != string::npos) periodNum = 2012;
        if(!isRealData && thisPhoton->Pt()>10.){
          TCGenParticle goodGenPhoton;
          float testDr = 9999;
          for (UInt_t j = 0; j<genPhotons.size(); j++){
            if(thisPhoton->DeltaR(genPhotons[j]) < 0.2 && genPhotons[j].GetStatus()==1 && fabs(genPhotons[j].Mother()) == 22){
              if(thisPhoton->DeltaR(genPhotons[j])<testDr){
                goodGenPhoton = genPhotons[j];
                testDr = thisPhoton->DeltaR(genPhotons[j]);
              }
            }
          }
          if (testDr < 0.2){
            corrPhoPt = phoCorrector->GetCorrEtMC(R9Cor, periodNum, thisPhoton->Pt(), thisPhoton->Eta(), goodGenPhoton.E());
            //cout<<"uncor pt: "<<thisPhoton->Pt()<<" cor pt: "<<corrPhoPt<<endl;
            thisPhoton->SetPtEtaPhiM(corrPhoPt,thisPhoton->Eta(),thisPhoton->Phi(),0.0);
          }
          //else cout<<" no match, pt: "<<thisPhoton->Pt()<<endl;
        }else if (isRealData && thisPhoton->Pt()>10.){
          corrPhoPt = phoCorrector->GetCorrEtData(thisPhoton->R9(), periodNum, thisPhoton->Pt(), thisPhoton->Eta());
          thisPhoton->SetPtEtaPhiM(corrPhoPt,thisPhoton->Eta(),thisPhoton->Phi(),0.0);
        }
       
      }


      ////// Currently Using Cut-Based Photon ID, 2012

      if (dumps) PhotonDump(thisPhoton, mediumPhID, mediumPhIso, phDump1);
      if (PassPhotonID(thisPhoton, mediumPhID)) photonsID.push_back(*thisPhoton);
      if (PassPhotonID(thisPhoton, mediumPhID) && PassPhotonIso(thisPhoton, mediumPhIso)){
        photonsIDIso.push_back(*thisPhoton);
        if (engCor) photonsIDIsoUnCor.push_back(*clonePhoton);
      }




    }
    sort(photonsID.begin(), photonsID.end(), PhotonSortCondition);
    sort(photonsIDIso.begin(), photonsIDIso.end(), PhotonSortCondition);
  }

  ////////////////////////
  // Analysis selection //
  ////////////////////////

  TLorentzVector ZP4;
  TLorentzVector ZP4Gen;
  TLorentzVector GP4;
  TLorentzVector uncorGP4;
  TLorentzVector GP4Gen;
  TLorentzVector JetOne, JetTwo;
  float GP4scEta = -9999;
  float eventWeight = 1;
  float eventWeightPU = 1;
  float eventWeightPho = 1;
  float eventWeightLep = 1;
  float eventWeightTrig = 1;

  //float deltaPhiJetMET = 2*TMath::Pi();
  //if (jetP4.size() > 0) deltaPhiJetMET= DeltaPhiJetMET(metP4, jetP4, eventWeight);

  TLorentzVector lepton1;
  int            lepton1int =-1;
  TLorentzVector lepton2;
  int            lepton2int =-1;
  TLorentzVector uncorLepton1;
  TLorentzVector uncorLepton2;

  if (selection == "electron" || selection == "eeGamma") {

    //////////////////////
    // 2 good electrons //
    //////////////////////
    if (electronsID.size() < 2) return kTRUE;

    bool firstEl = false;
    bool bothEls = false;
    for (UInt_t i = 0; i<electronsID.size(); i++){
      if (!firstEl && (electronsID[i].Pt() >elePtCut[1])){
        firstEl = true;
      }else if (firstEl && (electronsID[i].Pt() >elePtCut[1])){
        bothEls = true;
        break;
      }
    }
    if (!bothEls) return kTRUE;

    hm->fill1DHist(6,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5, eventWeight,"Misc");
    hm->fill1DHist(6,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
    ++nEvents[5];

    if (electronsIDIso.size() < 2) return kTRUE;
    firstEl = false;
    bothEls = false;
    for (UInt_t i = 0; i<electronsIDIso.size(); i++){
      if (!firstEl && (electronsIDIso[i].Pt() >elePtCut[1])){
        firstEl = true;
      }else if (firstEl && (electronsIDIso[i].Pt() >elePtCut[1])){
        bothEls = true;
        break;
      }
    }
    if (!bothEls) return kTRUE;






    //////////////////
    // 2 good muons //
    //////////////////

  } else if (selection == "muon" || selection =="mumuGamma") {
    if (eventNumber == EVENTNUMBER) cout<<"two ID muons?: "<<muonsID.size()<<endl;
    if (muonsID.size() < 2) return kTRUE;

    bool firstMu = false;
    bool bothMus = false;
    for (UInt_t i = 0; i<muonsID.size(); i++){
      if (!firstMu && (muonsID[i].Pt() >muPtCut[1])){
        firstMu = true;
      }else if (firstMu && (muonsID[i].Pt() >muPtCut[1])){
        bothMus = true;
        break;
      }
    }
    if (!bothMus) return kTRUE;

    hm->fill1DHist(6,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5, eventWeight,"Misc");
    hm->fill1DHist(6,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
    ++nEvents[5];

    if (eventNumber == EVENTNUMBER) cout<<"two ISO muons?: "<<muonsIDIso.size()<<endl;
    if (muonsIDIso.size() < 2) return kTRUE;
    firstMu = false;
    bothMus = false;
    for (UInt_t i = 0; i<muonsIDIso.size(); i++){
      if (!firstMu && (muonsIDIso[i].Pt() >muPtCut[1])){
        firstMu = true;
      }else if (firstMu && (muonsIDIso[i].Pt() >muPtCut[1])){
        bothMus = true;
        break;
      }
    }
    if (!bothMus) return kTRUE;

  }

  if (period.find("2011") != string::npos){
    if (doScaleFactors) eventWeight   *= weighter->PUWeight(nPUVertices);
    eventWeightPU   *= weighter->PUWeight(nPUVertices);
  }else if (period.find("2012") != string::npos){
    if (doScaleFactors) eventWeight   *= weighter->PUWeight(nPUVerticesTrue);
    eventWeightPU   *= weighter->PUWeight(nPUVerticesTrue);
  }


  ///////////////////////
  // Cosmics rejection //
  ///////////////////////

  //if (selection == "mumuGamma"){
  //if (muons.size() > 1 && CosmicMuonFilter(muons[0], muons[1])) return kTRUE;
  //}
  hm->fill1DHist(7,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5, eventWeight,"Misc");
  hm->fill1DHist(7,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  ++nEvents[6];


  //////////////////////
  // Good Z selection //
  //////////////////////

  bool goodZ = false;
  float SCetaEl1 = -99999;
  float SCetaEl2 = -99999;
  if(selection == "eeGamma"){
    if (eventNumber == EVENTNUMBER) cout<<goodZ<<endl;
    if (engCor || doEleReg) goodZ = FindGoodZElectron(electronsIDIso,electronsIDIsoUnCor,&lepton1,&lepton2,&uncorLepton1,&uncorLepton2,&ZP4,&SCetaEl1,&SCetaEl2,&lepton1int,&lepton2int);
    else goodZ = FindGoodZElectron(electronsIDIso,&lepton1,&lepton2,&ZP4,&SCetaEl1,&SCetaEl2,&lepton1int,&lepton2int);
    if (eventNumber == EVENTNUMBER) cout<<"el num: "<<electronsIDIso.size()<<" goodZ: "<<goodZ<<endl;
    if (!goodZ) return kTRUE;
    if (!isRealData){ 
      if (doScaleFactors){
        eventWeight   *= weighter->ElectronTriggerWeight(lepton1, lepton2);
        if (period.find("2011") != string::npos){
          eventWeight   *= getEfficiencyWeight( &lepton1, CorrectionType::CENTRAL, atoi(period.c_str()));
          eventWeight   *= getEfficiencyWeight( &lepton2, CorrectionType::CENTRAL, atoi(period.c_str()));
        }else if (period.find("2012") != string::npos){
          eventWeight   *= weighter->ElectronSelectionWeight(lepton1);
          eventWeight   *= weighter->ElectronSelectionWeight(lepton2);
        }
      }
      eventWeightTrig   *= weighter->ElectronTriggerWeight(lepton1, lepton2);
      if (period.find("2011") != string::npos){
        eventWeightLep   *= getEfficiencyWeight( &lepton1, CorrectionType::CENTRAL, atoi(period.c_str()));
        eventWeightLep   *= getEfficiencyWeight( &lepton2, CorrectionType::CENTRAL, atoi(period.c_str()));
      }else if (period.find("2012") != string::npos){
        eventWeightLep   *= weighter->ElectronSelectionWeight(lepton1);
        eventWeightLep   *= weighter->ElectronSelectionWeight(lepton2);
      }

    }
    //cout<<eventWeight<<endl;
  }else{
    if (engCor) goodZ = FindGoodZMuon(muonsIDIso,muonsIDIsoUnCor,&lepton1,&lepton2,&uncorLepton1,&uncorLepton2,&ZP4,&lepton1int,&lepton2int);
    else goodZ = FindGoodZMuon(muonsIDIso,&lepton1,&lepton2,&ZP4,&lepton1int,&lepton2int);
    if (eventNumber == EVENTNUMBER) cout<<"goodZ?: "<<goodZ<<endl;
    if (!goodZ) return kTRUE;
    if (eventNumber == EVENTNUMBER) cout<<"goodZ?: "<<goodZ<<endl;
    if (!isRealData){ 
      if (doScaleFactors){
        eventWeight   *= weighter->MuonTriggerWeight(lepton1, lepton2);
        eventWeight   *= weighter->MuonSelectionWeight(lepton1);
        eventWeight   *= weighter->MuonSelectionWeight(lepton2);
      }
      eventWeightTrig   *= weighter->MuonTriggerWeight(lepton1, lepton2);
      eventWeightLep   *= weighter->MuonSelectionWeight(lepton1);
      eventWeightLep   *= weighter->MuonSelectionWeight(lepton2);
    }
    //cout<<eventWeight<<endl;
  }

  hm->fill1DHist(8,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5, eventWeight,"Misc");
  hm->fill1DHist(8,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  ++nEvents[7];

  //////////////////////////////////////////////////////////////
  // lepton, gamma, Z plots directly after diLepton selection //
  //////////////////////////////////////////////////////////////



  if (photonsIDIso.size()>0 && ZP4.M() > 50 ){
    StandardPlots(lepton1,lepton2,photonsIDIso[0],eventWeight,"PreSel", "PreSel");
  }else if (ZP4.M() > zMassCut[0]){
    TLorentzVector noPhoton;
    noPhoton.SetPtEtaPhiE(0.00001,0,0,0.00001);
    StandardPlots(lepton1,lepton2,noPhoton,eventWeight,"PreSel", "PreSel");
  }

  if (ZP4.M() > zMassCut[0]) hm->fill1DHist(primaryVtx->GetSize(),"h1_pvMultPreSel_SUFFIX", "Multiplicity of PVs;N_{PV};N_{evts}", 25, 0.5, 25.5, eventWeight,"PreSel");


  ////////////
  // Z mass //
  ////////////

  if (eventNumber == EVENTNUMBER) cout<<ZP4.M()<<endl;
  if ((ZP4.M() < zMassCut[0] || ZP4.M() > zMassCut[1])) return kTRUE;  
  hm->fill1DHist(9,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5, eventWeight,"Misc");
  hm->fill1DHist(9,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  ++nEvents[8];
  hm->fill1DHist(ZP4.M(),"h1_diLeptonMassPreSel_SUFFIX", "M_{ll}; M_{ll};N_{evts}", 40, 70., 110.,eventWeight,"PreSel");

  if (dumps){
    for (Int_t i = 0; i < photonsIDIso.size(); ++i) {
      PhotonDump2(&photonsIDIso[i], mediumPhID, mediumPhIso, lepton1, lepton2, phDump2);
    }
  }





  ////////////////////////////
  //**ZGamma** Gamma Energy //
  ////////////////////////////

  if(selection == "mumuGamma" || selection == "eeGamma")
  {

    if (photonsIDIso.size() < 1) return kTRUE;
    //if (photonsIDIso.size() > 1) cout<<"more than one good photon?! no fuckin way man"<<endl;
    //if (R9switch && (R9Cor < R9Cut)) return kTRUE;

    bool goodPhoton = false;
    for (UInt_t i = 0; i<photonsIDIso.size(); i++){

      GP4 = photonsIDIso[i];    // define GP4
      GP4scEta = photonsIDIso[i].SCEta();
      if (engCor) uncorGP4 = photonsIDIsoUnCor[i]; 
      if (eventNumber == EVENTNUMBER) cout<<"dr1: "<<GP4.DeltaR(lepton1)<<" dr2: "<<GP4.DeltaR(lepton2)<<endl<<" pt/M: "<<GP4.Pt()/(GP4+ZP4).M()<<endl;
      if ((GP4.DeltaR(lepton1)<drCut || GP4.DeltaR(lepton2)<drCut)) continue;
      if (GP4.Pt()/(GP4+ZP4).M() > gammaPtCut[0] && GP4.Pt() > gammaPtCut[1]) goodPhoton = true;
      if(goodPhoton){
        R9 = photonsIDIso[i].R9();
        R9Cor = R9;
        if (doR9Cor){
          if (suffix != "DATA"){
            if (period == "2011"){
              if (fabs(photonsIDIso[i].SCEta()) < 1.479){
                R9Cor = R9*1.0048;
              }else{
                R9Cor = R9*1.00492;
              }
            } else if (period == "2012"){
              if (fabs(photonsIDIso[i].SCEta()) < 1.479){
                R9Cor = R9*1.0045 + 0.0010;
              }else{
                R9Cor = R9*1.0086 - 0.0007;
              }
            }
          }
          break;
        }
      }
    }
    if(!goodPhoton) return kTRUE;
    //cout<<GP4.Pt()<<endl;

    if (doScaleFactors) eventWeight   *= weighter->GammaSelectionWeight(GP4, GP4scEta);
    eventWeightPho   *= weighter->GammaSelectionWeight(GP4, GP4scEta);
  }

  hm->fill1DHist(10,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5, eventWeight,"Misc");
  hm->fill1DHist(10,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  ++nEvents[9];

  ///////////////////////////////
  //**ZGamma** Photon DR Stuff //
  ///////////////////////////////

  if(selection == "mumuGamma" || selection == "eeGamma")
  {
    // Make an FSR plot for 1337 normalization of the ZGToMuMuG set

    if (ZP4.M() > zMassCut[0] && (GP4.DeltaR(lepton1)<drCut || GP4.DeltaR(lepton2)<drCut)){
      hm->fill1DHist(ZP4.M(),"h1_diLeptonMassPreSelFSR_SUFFIX", "M_{ll}; M_{ll};N_{evts}", 30, 50., 80.,eventWeight,"PreSel");
    }


  //  if ((GP4.DeltaR(lepton1)<drCut || GP4.DeltaR(lepton2)<drCut)) return kTRUE;
    //if (GP4.DeltaR(muonsIDIso[0].P4())>0.6 && GP4.DeltaR(muonsIDIso[1].P4())>0.6) return kTRUE;

  }
  hm->fill1DHist(11,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5, eventWeight,"Misc");
  hm->fill1DHist(11,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  ++nEvents[10];

  /////////////////////
  // 3rd lepton veto //
  /////////////////////

  //if (extraLeptons.size() > 0) return kTRUE;
  //if ((electronsIDIso.size() > 0 || muonsIDIso.size() > 2) && (selection == "muon" || selection == "mumuGamma")) return kTRUE;
  //if ((muonsIDIso.size() > 0 || electronsIDIso.size() > 2) && (selection == "electron" || selection == "eeGamma")) return kTRUE;
  //if ((muonsIDIso.size() > 0 || electronsIDIso.size() > 0) && (selection == "eGamma" || selection == "muGamma" || selection == "gamma")) return kTRUE;
  hm->fill1DHist(12,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5, eventWeight,"Misc");
  hm->fill1DHist(12,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  ++nEvents[11];

  ////////////////
  // b-jet veto //
  ////////////////

  hm->fill1DHist(13,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5, eventWeight,"Misc");
  hm->fill1DHist(13,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  ++nEvents[12];

  /////////
  // MET //
  /////////

  hm->fill1DHist(14,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5, eventWeight,"Misc");
  hm->fill1DHist(14,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  ++nEvents[13];

  ////////////
  // MET/QT //
  ////////////

  hm->fill1DHist(15,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5, eventWeight,"Misc");
  hm->fill1DHist(15,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  ++nEvents[14];

  //////////////////////
  // Jet Multiplicity //
  //////////////////////

  hm->fill1DHist(16,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5, eventWeight,"Misc");
  hm->fill1DHist(16,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  ++nEvents[15];

  //////////////////////////
  //**ZGamma**    MZG Low //
  //////////////////////////

  if (selection == "mumuGamma"|| selection == "eeGamma" || selection == "gammaJets")  if ((ZP4+GP4).M() < MassZG[0]) return kTRUE;
  if (doScaleFactors) eventWeight *= weighter->HqtWeight((ZP4+GP4));
  hm->fill1DHist(17,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5, eventWeight,"Misc");
  hm->fill1DHist(17,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  ++nEvents[16];

  ///////////////////////////
  //**ZGamma**    MZG High //
  ///////////////////////////

  if (selection == "mumuGamma"|| selection == "eeGamma" || selection == "gammaJets")  if ((ZP4+GP4).M() > MassZG[1]  ) return kTRUE;
  hm->fill1DHist(18,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5, eventWeight,"Misc");
  hm->fill1DHist(18,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  ++nEvents[17];

  //////////////////////////
  //**ZGamma** MZ+MZG Cut //
  //////////////////////////

  if(selection == "mumuGamma" || selection == "eeGamma") if(ZP4.M()+(ZP4+GP4).M() < MZpMZGcut) return kTRUE;
  hm->fill1DHist(19,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5, eventWeight,"Misc");
  hm->fill1DHist(19,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  ++nEvents[18];



  ////////////////////////
  //Angular Correlations//
  ////////////////////////

  //Define required angles
  double polarZBoost= -999.;
  double polarGBoost= -999.;
  double cospolarZBoost= -999.;
  double cospolarGBoost= -999.;
  double angleGLepNeg = -999;
  double angleGLepPos = -999;
  double azimuthLepPlane = -999;
  double diffPlane = -999;

  if(selection == "mumuGamma" || selection == "eeGamma")
  {

    TLorentzVector HP4 = GP4 + ZP4; //Define 4-vector for the Higgs
    TLorentzVector ZBoost(ZP4);
    TLorentzVector GBoost(GP4);
    TLorentzVector Lep1Boost;
    TLorentzVector Lep2Boost;
    if (selection == "mumuGamma"){
      Lep1Boost = muonsIDIso[0];
      Lep2Boost = muonsIDIso[1];
    }else{
      Lep1Boost = electronsIDIso[0];
      Lep2Boost = electronsIDIso[1];
    }

    TVector3 v_z(0.,0.,1.);
    TVector3 v_y(0.,1.,0.);
    TVector3 v_x(1.,0.,0.);
    TVector3 v_1;
    TVector3 v_2;

    //// Boost everything that's important into the Higgs CoM
    ZBoost.Boost(-HP4.BoostVector());
    GBoost.Boost(-HP4.BoostVector());
    Lep1Boost.Boost(-HP4.BoostVector());
    Lep2Boost.Boost(-HP4.BoostVector());

    //Rotate gamma and Z into the Z-X plane
    Lep1Boost.RotateZ(-ZBoost.Phi());
    Lep2Boost.RotateZ(-ZBoost.Phi());
    GBoost.RotateZ(-ZBoost.Phi());
    ZBoost.RotateZ(-ZBoost.Phi());


    //Sanity Checks
    hm->fill1DHist(ZBoost.P()-GBoost.P(),"h1_HiggsCoMPSanity_SUFFIX","Higgs mom in CoM (sanity)",100,-50,50,eventWeight,"ZGamma");
    hm->fill1DHist(ZBoost.Phi(),"h1_PhiOfZSanity_SUFFIX","Z Phi in CoM (sanity)",100,-2*TMath::Pi(),2*TMath::Pi(),eventWeight,"ZGamma");
    hm->fill1DHist(GBoost.Phi(),"h1_PhiOfGSanity_SUFFIX","gamma Phi in CoM (sanity)",100,-2*TMath::Pi(),2*TMath::Pi(),eventWeight,"ZGamma");

    //Define angles using above boosted vectors
    polarZBoost= ZBoost.Theta();
    polarGBoost= GBoost.Theta();
    cospolarZBoost = ZBoost.CosTheta();
    cospolarGBoost = GBoost.CosTheta();


    //Define angles between gamma and leptons
    if ((selection == "mumuGamma" && muonsIDIso[0].Charge() < 0) || (selection == "eeGamma" && electronsIDIso[0].Charge() <0 )){
      angleGLepNeg = GBoost.Angle(Lep1Boost.Vect());
      angleGLepPos = GBoost.Angle(Lep2Boost.Vect());
    } else {
      angleGLepNeg = GBoost.Angle(Lep2Boost.Vect());
      angleGLepPos = GBoost.Angle(Lep1Boost.Vect());
    }

    //Correct for negative or >Pi angles
    if (fabs(angleGLepNeg)>TMath::Pi()) {
      angleGLepNeg =2*TMath::Pi()-fabs(angleGLepNeg);
    } else {
      angleGLepNeg = fabs(angleGLepNeg);
    }
    if (fabs(angleGLepPos)>TMath::Pi()) {
      angleGLepPos =2*TMath::Pi()-fabs(angleGLepPos);
    } else {
      angleGLepPos = fabs(angleGLepPos);
    }

    //Define ZGamma to be in x-z plane
    if (fabs(ZBoost.Angle(v_z)) > 0 && fabs(ZBoost.Angle(v_z)) < TMath::Pi()) {
      v_1 = ZBoost.Vect().Cross(v_z);
    } else {
      v_1 = ZBoost.Vect().Cross(v_x);
    }

    //Define diLepton plane
    v_2 = Lep1Boost.Vect().Cross(Lep2Boost.Vect());
    azimuthLepPlane = v_2.Phi();

    //Define angle between diLepton and ZGamma planes
    diffPlane = v_1.Angle(v_2);
    if (fabs(diffPlane)>TMath::Pi()) {
      diffPlane =2*TMath::Pi()-fabs(diffPlane);
    } else {
      diffPlane = fabs(diffPlane);
    }

  }


  //////////////////
  // Angular Cuts //
  //////////////////

  if (fabs(cospolarZBoost) > AngCut) return kTRUE;
  hm->fill1DHist(20,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5, eventWeight,"Misc");
  hm->fill1DHist(20,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  ++nEvents[19];

  //////////////
  // M12 Cuts //
  //////////////

  if (CalculateM12sqrd(ZP4,GP4) > M12Cut) return kTRUE;
  hm->fill1DHist(21,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5, eventWeight,"Misc");
  hm->fill1DHist(21,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  ++nEvents[20];


  // -------------------------- MVA stuff -------------------------------------------
#ifdef USE_MVA
  // the sequence of cuts is a bit different for the pre-selection
  // thi si the latest where we can plug in the BDT for now

  // *additional* cuts for the BDT pre-selection on top of the ones already applied
  // assign values to the variables used in the BDT


  TLorentzVector HP4tmp = ZP4+GP4;
  TLorentzVector totalTrackP4tmp = sumJetP4; 
  if (selection == "mumuGamma"){
    for (int j = 0; j < (int)muonsIDIso.size(); ++j) totalTrackP4tmp += muonsIDIso[j];
    for (int j = 0; j < (int)extraLeptons.size(); ++j) totalTrackP4tmp += extraLeptons[j];
  }else if (selection == "eeGamma"){
    for (int j = 0; j < (int)electronsIDIso.size(); ++j) totalTrackP4tmp += electronsIDIso[j];
    for (int j = 0; j < (int)extraLeptons.size(); ++j) totalTrackP4tmp += extraLeptons[j];
  }
  float vtxPtRat2tmp = (totalTrackP4tmp.Pt()-HP4tmp.Pt())/(totalTrackP4tmp.Pt() + HP4tmp.Pt());

  _diffZGvector       = (ZP4-GP4).Pt();
  _threeBodyPt        = (ZP4+GP4).Pt();
  _GPt                = GP4.Pt();
  _cosZ               = cospolarZBoost;
  _diffPlaneMVA       = diffPlane;
  _vtxVariable        = vtxPtRat2tmp;
  _dr1                = lepton1.DeltaR(GP4);
  _dr2                = lepton2.DeltaR(GP4);
  _M12                = CalculateM12sqrd(ZP4,GP4);
  _threeBodyMass      = (ZP4+GP4).M();

  // get the MVA discriminators for the considered methods

  //             [higgs mass points]
  Bool_t passBdtCut[N_HIGGS_MASSES] = {kFALSE};
  Bool_t passAllBdtCuts[N_HIGGS_MASSES] = {kTRUE};


  //                    [mva method][higgs mass point]
  Float_t tmvaValue[N_DISCR_METHODS][N_HIGGS_MASSES] = {{0.0}};

  int discr = BDTG; // use only this one for now
  //int discr = MLPBNN; // use only this one for now

  for (int mh = 0; mh<N_HIGGS_MASSES; ++mh) {

    TString label = TString::Format("%s_%s_MVA_HZG%i", discrMethodName[discr].Data(), discrSampleName.Data(),
        mvaHiggsMassPoint[mh]);
    tmvaValue[discr][mh] = tmvaReader->EvaluateMVA(label.Data());

    if (tmvaValue[discr][mh] > bdtCut[mh]) passBdtCut[mh] = kTRUE;
    passAllBdtCuts[mh] = (passAllBdtCuts[mh] && passBdtCut[mh]);

  }



  // here we can count events, fill histograms etc

  for (int mh = 0; mh<N_HIGGS_MASSES; ++mh) { // loop over mass points -> at the moment I put only one!

    if (!passAllBdtCuts[mh]) return kTRUE;
    if (nEvents[0]%2 == 0 && suffix != "DATA") return kTRUE;
    if (suffix != "DATA") eventWeight = eventWeight*2;
    hm->fill1DHist(tmvaValue[discr][mh],"h1_mvaValue_SUFFIX", "MVA value;;N_{evts}", 20, -1, 1, eventWeight,"ZGamma");
  }

  // we can also produce histograms and count events based on passing a subset of th discriminator cuts


#endif
  // ------------------------------------- End of MVA stuff -------------------------------------

  ///////////////////////
  // MVA Event Counter //
  ///////////////////////

  hm->fill1DHist(22,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5, eventWeight,"Misc");
  hm->fill1DHist(22,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  ++nEvents[21];

  //////////////
  // VBF Cuts //
  //////////////

  if (VBFcuts){

  }
  hm->fill1DHist(23,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5, eventWeight,"Misc");
  hm->fill1DHist(23,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  ++nEvents[22];

  //////////////
  // Non-Cuts //
  //////////////

  // Cross Section Weighting //
  LumiXSWeight(&unBinnedLumiXS);
  if(doLumiXS) eventWeight *= unBinnedLumiXS;

  if ( (GP4+ZP4).M() > 100){
    hm->fill1DHist(50,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5,eventWeight,"Misc");
    hm->fill1DHist(50,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
    hm->fill1DHist(51,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5,eventWeightPU,"Misc");
    hm->fill1DHist(51,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
    hm->fill1DHist(52,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5,eventWeightPho,"Misc");
    hm->fill1DHist(52,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
    hm->fill1DHist(53,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5,eventWeightTrig,"Misc");
    hm->fill1DHist(53,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
    hm->fill1DHist(54,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5,eventWeightLep,"Misc");
    hm->fill1DHist(54,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
    ++nEvents[51];
  }

  if ( (GP4+ZP4).M() > 120 && (GP4+ZP4).M() < 150){
    hm->fill1DHist(70,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5,eventWeight,"Misc");
    hm->fill1DHist(70,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
  }

  ////////////////////////////
  // Fill lepton histograms //
  ////////////////////////////

  if (selection == "muon" || selection == "mumuGamma") {

    LeptonBasicPlots(muonsIDIso[0], muonsIDIso[1], eventWeight);
    hm->fill1DHist(muonsIDIso[0].Charge(),"h1_leadLeptonCharge_SUFFIX", "Charge leading lepton;Charge;N_{evts}", 5, -2.5, 2.5, eventWeight,"Lepton");     
    hm->fill1DHist(muonsIDIso[1].Charge(),"h1_trailingLeptonCharge_SUFFIX", "Charge trailing lepton;Charge;N_{evts}", 5, -2.5, 2.5, eventWeight,"Lepton");     

  } else if (selection == "electron" || selection == "eeGamma") {

    LeptonBasicPlots(electronsIDIso[0], electronsIDIso[1], eventWeight);
    hm->fill1DHist(electronsIDIso[0].Charge(),"h1_leadLeptonCharge_SUFFIX", "Charge leading lepton;Charge;N_{evts}", 5, -2.5, 2.5, eventWeight,"Lepton");     
    hm->fill1DHist(electronsIDIso[1].Charge(),"h1_trailingLeptonCharge_SUFFIX", "Charge trailing lepton;Charge;N_{evts}", 5, -2.5, 2.5, eventWeight,"Lepton");     
  } 

  DileptonBasicPlots(ZP4, eventWeight);

  //////////////////////////////
  // Scale Factor Corrections //
  //////////////////////////////

  // Up to this point, eventWeight should only include pileup reweighting.  Now we will add
  // scalefactors on every selection, while still saving the uncorrected versions.  We will
  // also save the total up and down errors on the scalefactors to be used for systematic
  // shape analysis.

  //float phoEffWeight; 
  //float phoEffWeightUp; 
  //float phoEffWeightDown; 

  if (selection == "eeGamma"){
    //cout<<"ele thing firing off"<<endl;
    TRandom3 *rl1 = new TRandom3(0);
    TRandom3 *rl2 = new TRandom3(0);

    float wgt1        = getEfficiencyWeight( &lepton1, CorrectionType::CENTRAL, atoi(period.c_str()));
    float wgt_hi1     = getEfficiencyWeight( &lepton1, CorrectionType::UPPER, atoi(period.c_str()));
    float wgt_lo1     = getEfficiencyWeight( &lepton1, CorrectionType::LOWER, atoi(period.c_str()));
    float ptcorr1     = getCorrectedPt( &lepton1, CorrectionType::CENTRAL, atoi(period.c_str()),rl1);
    float ptcorr_hi1  = getCorrectedPt( &lepton1, CorrectionType::UPPER, atoi(period.c_str()),rl1);
    float ptcorr_lo1  = getCorrectedPt( &lepton1, CorrectionType::LOWER, atoi(period.c_str()),rl1);
    float wgt2        = getEfficiencyWeight( &lepton2, CorrectionType::CENTRAL, atoi(period.c_str()));
    float wgt_hi2     = getEfficiencyWeight( &lepton2, CorrectionType::UPPER, atoi(period.c_str()));
    float wgt_lo2     = getEfficiencyWeight( &lepton2, CorrectionType::LOWER, atoi(period.c_str()));
    float ptcorr2     = getCorrectedPt( &lepton2, CorrectionType::CENTRAL, atoi(period.c_str()),rl2);
    float ptcorr_hi2  = getCorrectedPt( &lepton2, CorrectionType::UPPER, atoi(period.c_str()),rl2);
    float ptcorr_lo2  = getCorrectedPt( &lepton2, CorrectionType::LOWER, atoi(period.c_str()),rl2);

    //cout<<"mid: "<<ptcorr1<<" plus: "<<ptcorr_hi1-ptcorr1<<" minus: "<<ptcorr1-ptcorr_lo1<<endl;

    TLorentzVector el1High(lepton1);
    TLorentzVector el1Mid(lepton1);
    TLorentzVector el1Low(lepton1);
    TLorentzVector el2High(lepton2);
    TLorentzVector el2Mid(lepton2);
    TLorentzVector el2Low(lepton2);

    el1High.SetPtEtaPhiM(ptcorr_hi1,lepton1.Eta(),lepton1.Phi(),0.000511);
    el1Mid.SetPtEtaPhiM(ptcorr1,lepton1.Eta(),lepton1.Phi(),0.000511);
    el1Low.SetPtEtaPhiM(ptcorr_lo1,lepton1.Eta(),lepton1.Phi(),0.000511);
    el2High.SetPtEtaPhiM(ptcorr_hi2,lepton2.Eta(),lepton2.Phi(),0.000511);
    el2Mid.SetPtEtaPhiM(ptcorr2,lepton2.Eta(),lepton2.Phi(),0.000511);
    el2Low.SetPtEtaPhiM(ptcorr_lo2,lepton2.Eta(),lepton2.Phi(),0.000511);

    bool goodEle = false;
    goodEle = (fabs(lepton1.Eta()) < 1.4442 && fabs(lepton2.Eta()) < 1.4442);
    string catName = "";
    if (goodEle && (fabs(GP4.Eta()) < 1.4442) && R9Cor >= R9Cut ){       catName = "CAT1";
    } else if (goodEle && (fabs(GP4.Eta()) < 1.4442) && R9Cor < R9Cut ){ catName = "CAT4";
    } else if ( fabs(GP4.Eta()) < 1.4442){                               catName = "CAT2";
    } else {                                                             catName = "CAT3";
    }
    

    hmEleSmear->fill1DHist((el1Mid+el2Mid+GP4).M(),"h1_middleElSmear_"+catName+"_SUFFIX","middle Smear;M (GeV);Entries",65,115,190,eventWeight);
    hmEleSmear->fill1DHist((el1High+el2High+GP4).M(),"h1_highElSmear_"+catName+"_SUFFIX","high Smear;M (GeV);Entries",65,115,190,eventWeight);
    hmEleSmear->fill1DHist((el1Low+el2Low+GP4).M(),"h1_lowElSmear_"+catName+"_SUFFIX","low Smear;M (GeV);Entries",65,115,190,eventWeight);
    hmEleSmear->fill1DHist((ZP4+GP4).M(),"h1_noCorSmear_"+catName+"_SUFFIX","no cor;M (GeV);Entries",65,115,190,eventWeight);
    hmEleIDISO->fill1DHist((lepton1+lepton2+GP4).M(),"h1_noCorIDISO_"+catName+"_SUFFIX","no weight cor;M (GeV);Entries",65,115,190,eventWeight);

    // make a bunch of toy histos 
    
    float elEffMean1 = wgt1;
    float elEffSig1  = max(fabs(wgt1-wgt_hi1),fabs(wgt1-wgt_lo1));
    float elEffMean2 = wgt2;
    float elEffSig2  = max(fabs(wgt2-wgt_hi2),fabs(wgt2-wgt_lo2));
    TRandom3 *g1 = new TRandom3(0);
    TRandom3 *g2 = new TRandom3(0);
    float tw1,tw2 = 0;
    stringstream itoa;
    for (int i =0; i<100; i++){
      tw1 = g1->Gaus(elEffMean1,elEffSig1);
      tw2 = g2->Gaus(elEffMean2,elEffSig2);

      itoa << i;

      hmEleIDISO->fill1DHist((lepton1+lepton2+GP4).M(),"h1_toy"+itoa.str()+"_ElIDISO_"+catName+"_SUFFIX","low IDISO;M (GeV);Entries",65,115,190,eventWeight*tw2*tw1);
      itoa.str("");
    }


  }

  ////////////////////////////
  // Fill ZGamma histograms //
  ////////////////////////////



  if(selection == "mumuGamma" || selection == "eeGamma")
  {

    if (MassZG[0] < 90){
      hm->fill1DHist((ZP4+GP4).M(),"h1_InvariantMass_SUFFIX","Invariant Mass (Z#gamma);Mass (GeV);Entries",100,0,100,eventWeight,"ZGamma");
    }else{
      hm->fill1DHist((ZP4+GP4).M(),"h1_InvariantMass_SUFFIX","Invariant Mass (Z#gamma);Mass (GeV);Entries",100,100,200,eventWeight,"ZGamma");
    }
    hm->fill1DHist((ZP4+GP4).M(),"h1_InvariantMassCut_SUFFIX","Invariant Mass (H->Z#gamma);Mass (GeV);Entries",65,115,190,eventWeight,"ZGamma");


    ///// OK Lets make some stoyan plots /////

    hmHiggs->fill1DHist((ZP4+GP4).M(),"h1_InvariantMassRecoStoyan1GevFULLRANGE_SUFFIX","Invariant Mass (H->Z#gamma);Mass (GeV);Entries",90,90,190,eventWeight);
    hmHiggs->fill1DHist((ZP4+GP4).M(),"h1_InvariantMassRecoStoyan1Gev_SUFFIX","Invariant Mass (H->Z#gamma);Mass (GeV);Entries",65,115,190,eventWeight);
    hmHiggs->fill1DHist(R9, "h1_R9FullStoyan_SUFFIX","R9;R9;Entries",100,0,1,eventWeight);
    hmHiggs->fill1DHist(R9Cor, "h1_R9CorFullStoyan_SUFFIX","R9 Corrected;R9;Entries",100,0,1,eventWeight);
    if (genHs.size() > 0) hmHiggs->fill1DHist(genHs[0].M()-(ZP4+GP4).M(),"h1_genHiggsMassResFull_SUFFIX", "Gen-Reco M_{ll#gamma} Full; #Delta M_{ll};N_{evts}", 40, -20., 20., eventWeight);
    



    bool goodLep = false;
    if (selection == "mumuGamma"){
      goodLep = (fabs(lepton1.Eta()) < 0.9 || fabs(lepton2.Eta()) < 0.9) && (fabs(lepton1.Eta()) < 2.1 && fabs(lepton2.Eta()) < 2.1);
    }else if (selection == "eeGamma"){
      goodLep = (fabs(SCetaEl1) < 1.4442 && fabs(SCetaEl2) < 1.4442);
    }

    /*
    float PhoRegWeight;
    if (selection == "mumuGamma"){
      if (suffix.find("Signal") != string::npos){
        PhoRegWeight = 0.219;
      }else PhoRegWeight = 0.587;
    }else{
      if (suffix.find("Signal") != string::npos){
        PhoRegWeight = 0.219;
      }else PhoRegWeight = 0.515;
    }
*/

    m_llg = (GP4+ZP4).M();

    if (dumps){
      if (selection == "mumuGamma"){
        if (doLooseMuIso){
          MuonDump(&muonsIDIso[lepton1int], tightMuID, looseMuIso, muDumpFinal);
          MuonDump(&muonsIDIso[lepton2int], tightMuID, looseMuIso, muDumpFinal);
        }else{
          MuonDump(&muonsIDIso[lepton1int], tightMuID, tightMuIso, muDumpFinal);
          MuonDump(&muonsIDIso[lepton2int], tightMuID, tightMuIso, muDumpFinal);
        }

      } else if (selection == "eeGamma"){
        ElectronDump(&electronsIDIso[lepton1int], looseElID, looseElIso, elDumpFinal);
        ElectronDump(&electronsIDIso[lepton2int], looseElID, looseElIso, elDumpFinal);
      }
    }

    int catNum = -1;
    // If lep1 or lep2 is in 0.9 and both are in 2.1 //
    if (goodLep && (fabs(GP4scEta) < 1.4442) ){
      if (R9Cor >= R9Cut){
        //CAT 1
        catNum = 1;
        hmHiggs->fill1DHist((ZP4+GP4).M(),"h1_InvariantMassReco1GevCAT1FULLRANGE_SUFFIX","Invariant Mass (H->Z#gamma);Mass (GeV);Entries",90,90,190,eventWeight);
        if (genHs.size() > 0) hmHiggs->fill1DHist(genHs[0].M()-(ZP4+GP4).M(),"h1_genHiggsMassResCAT1_SUFFIX", "Gen-Reco M_{ll#gamma} CAT1; #Delta M_{ll};N_{evts}", 40, -20., 20., eventWeight);
        hm->fill1DHist((ZP4+GP4).M(),"h1_InvariantMassReco1GevCAT1FULLRANGE_SUFFIX","Invariant Mass (H->Z#gamma);Mass (GeV);Entries",90,90,190,eventWeight);
        StandardPlots(lepton1,lepton2,GP4,eventWeight,"CAT1", "CAT1");
        hmHiggs->fill1DHist(R9, "h1_R9CAT1_SUFFIX","R9;R9;Entries",100,0,1,eventWeight);
        hm->fill1DHist(R9, "h1_R9CAT1_SUFFIX","R9;R9;Entries",100,0,1,eventWeight);
        hmHiggs->fill1DHist(R9Cor, "h1_R9CorCAT1_SUFFIX","R9Cor;R9Cor;Entries",100,0,1,eventWeight);
        hm->fill1DHist(R9Cor, "h1_R9CorCAT1_SUFFIX","R9Cor;R9Cor;Entries",100,0,1,eventWeight);
        hm->fill2DHist(lepton1.Eta(),lepton2.Eta(),"h2_dilepEtaCAT1_SUFFIX","Dilepton Eta CAT1; Eta (leading); Eta (trailing)", 50,-2.5,2.5,50,-2.5,2.5,eventWeight,"CAT1");
        m_llgCAT1 = (GP4+ZP4).M();
        hm->fill1DHist(60,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5,eventWeight,"Misc");
        hm->fill1DHist(60,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
        ++nEvents[61];
      }else{
        //CAT 4
        catNum = 4;
        hmHiggs->fill1DHist((ZP4+GP4).M(),"h1_InvariantMassReco1GevCAT4FULLRANGE_SUFFIX","Invariant Mass (H->Z#gamma);Mass (GeV);Entries",90,90,190,eventWeight);
        if (genHs.size() > 0) hmHiggs->fill1DHist(genHs[0].M()-(ZP4+GP4).M(),"h1_genHiggsMassResCAT4_SUFFIX", "Gen-Reco M_{ll#gamma} CAT4; #Delta M_{ll};N_{evts}", 40, -20., 20., eventWeight);
        hm->fill1DHist((ZP4+GP4).M(),"h1_InvariantMassReco1GevCAT4FULLRANGE_SUFFIX","Invariant Mass (H->Z#gamma);Mass (GeV);Entries",90,90,190,eventWeight);
        StandardPlots(lepton1,lepton2,GP4,eventWeight,"CAT4", "CAT4");
        hmHiggs->fill1DHist(R9, "h1_R9CAT4_SUFFIX","R9;R9;Entries",100,0,1,eventWeight);
        hm->fill1DHist(R9, "h1_R9CAT4_SUFFIX","R9;R9;Entries",100,0,1,eventWeight);
        hmHiggs->fill1DHist(R9Cor, "h1_R9CorCAT4_SUFFIX","R9Cor;R9Cor;Entries",100,0,1,eventWeight);
        hm->fill1DHist(R9Cor, "h1_R9CorCAT4_SUFFIX","R9Cor;R9Cor;Entries",100,0,1,eventWeight);
        hm->fill2DHist(lepton1.Eta(),lepton2.Eta(),"h2_dilepEtaCAT4_SUFFIX","Dilepton Eta CAT4; Eta (leading); Eta (trailing)", 50,-2.5,2.5,50,-2.5,2.5,eventWeight,"CAT4");
        m_llgCAT4 = (GP4+ZP4).M();
        hm->fill1DHist(63,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5,eventWeight,"Misc");
        hm->fill1DHist(63,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
        ++nEvents[64];
      }
    } else if ( fabs(GP4scEta) < 1.4442){
      // CAT 2
      catNum = 2;
      hmHiggs->fill1DHist((ZP4+GP4).M(),"h1_InvariantMassReco1GevCAT2FULLRANGE_SUFFIX","Invariant Mass (H->Z#gamma);Mass (GeV);Entries",90,90,190,eventWeight);
      if (genHs.size() > 0) hmHiggs->fill1DHist(genHs[0].M()-(ZP4+GP4).M(),"h1_genHiggsMassResCAT2_SUFFIX", "Gen-Reco M_{ll#gamma} CAT2; #Delta M_{ll};N_{evts}", 40, -20., 20., eventWeight);
      hm->fill1DHist((ZP4+GP4).M(),"h1_InvariantMassReco1GevCAT2FULLRANGE_SUFFIX","Invariant Mass (H->Z#gamma);Mass (GeV);Entries",90,90,190,eventWeight);
      StandardPlots(lepton1,lepton2,GP4,eventWeight,"CAT2", "CAT2");
      hmHiggs->fill1DHist(R9, "h1_R9CAT2_SUFFIX","R9;R9;Entries",100,0,1,eventWeight);
      hm->fill1DHist(R9, "h1_R9CAT2_SUFFIX","R9;R9;Entries",100,0,1,eventWeight);
      hmHiggs->fill1DHist(R9Cor, "h1_R9CorCAT2_SUFFIX","R9Cor;R9Cor;Entries",100,0,1,eventWeight);
      hm->fill1DHist(R9Cor, "h1_R9CorCAT2_SUFFIX","R9Cor;R9Cor;Entries",100,0,1,eventWeight);
      hm->fill2DHist(lepton1.Eta(),lepton2.Eta(),"h2_dilepEtaCAT2_SUFFIX","Dilepton Eta CAT2; Eta (leading); Eta (trailing)", 50,-2.5,2.5,50,-2.5,2.5,eventWeight,"CAT2");
      m_llgCAT2 = (GP4+ZP4).M();
      hm->fill1DHist(61,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5,eventWeight,"Misc");
      hm->fill1DHist(61,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
      ++nEvents[62];
    } else {
      // CAT 3
      catNum = 3;
      hm->fill1DHist((ZP4+GP4).M(),"h1_InvariantMassReco1GevCAT3FULLRANGE_SUFFIX","Invariant Mass (H->Z#gamma);Mass (GeV);Entries",90,90,190,eventWeight);
      if (genHs.size() > 0) hmHiggs->fill1DHist(genHs[0].M()-(ZP4+GP4).M(),"h1_genHiggsMassResCAT3_SUFFIX", "Gen-Reco M_{ll#gamma} CAT3; #Delta M_{ll};N_{evts}", 40, -20., 20., eventWeight);
      hmHiggs->fill1DHist((ZP4+GP4).M(),"h1_InvariantMassReco1GevCAT3FULLRANGE_SUFFIX","Invariant Mass (H->Z#gamma);Mass (GeV);Entries",90,90,190,eventWeight);
      hmHiggs->fill1DHist(R9, "h1_R9CAT3_SUFFIX","R9;R9;Entries",100,0,1,eventWeight);
      hm->fill1DHist(R9, "h1_R9CAT3_SUFFIX","R9;R9;Entries",100,0,1,eventWeight);
      hmHiggs->fill1DHist(R9Cor, "h1_R9CorCAT3_SUFFIX","R9Cor;R9Cor;Entries",100,0,1,eventWeight);
      hm->fill1DHist(R9Cor, "h1_R9CorCAT3_SUFFIX","R9Cor;R9Cor;Entries",100,0,1,eventWeight);
      StandardPlots(lepton1,lepton2,GP4,eventWeight,"CAT3", "CAT3");
      hm->fill2DHist(lepton1.Eta(),lepton2.Eta(),"h2_dilepEtaCAT3_SUFFIX","Dilepton Eta CAT3; Eta (leading); Eta (trailing)", 50,-2.5,2.5,50,-2.5,2.5,eventWeight,"CAT3");
      hm->fill2DHist(GP4.Eta(),GP4.Phi(),"h2_EtaVsPhiCAT3_SUFFIX","Photon Eta vs Phi CAT3; Eta Photon; Phi Photon", 50,-2.5,2.5,50,-TMath::Pi(),TMath::Pi(),eventWeight,"CAT3");
      if (GP4.Eta() > -2.0 && GP4.Eta() < -1.5 && GP4.Phi() > 1.0 && GP4.Phi() < 2){
        hm->fill2DHist(GP4.Eta(),GP4.Phi(),"h2_EtaVsPhiCAT3Zoom_SUFFIX","Photon Eta vs Phi CAT3; Eta Photon; Phi Photon", 50,-2.0,-1.5,100,1.0,2.0,eventWeight,"CAT3");
      }
      if (GP4.Eta() > -1.80 && GP4.Eta() < -1.75 && GP4.Phi() > 1.35 && GP4.Phi() < 1.40){
        hm->fill1DHist(runNumber,"h1_spikeVsRunCAT3_SUFFIX","runNumber for eta/phi spike CAT3; runNumber; entries", 6075,190456,196531,eventWeight,"CAT3");
      }
      m_llgCAT3 = (GP4+ZP4).M();
      hm->fill1DHist(62,"h1_acceptanceByCut_SUFFIX", "Weighted number of events passing cuts by cut; cut; N_{evts}", 100, 0.5, 100.5,eventWeight,"Misc");
      hm->fill1DHist(62,"h1_acceptanceByCutRaw_SUFFIX", "Raw number of events passing cuts; cut; N_{evts}", 100, 0.5, 100.5,1,"Misc");
      ++nEvents[63];

      
    }

    //High R9
    if (R9Cor>R9Cut){
      hmHiggs->fill1DHist((ZP4+GP4).M(),"h1_InvariantMassReco1GevHighR9FULLRANGE_SUFFIX","Invariant Mass (H->Z#gamma);Mass (GeV);Entries",90,90,190,eventWeight);
      if (genHs.size() > 0) hmHiggs->fill1DHist(genHs[0].M()-(ZP4+GP4).M(),"h1_genHiggsMassResHighR9_SUFFIX", "Gen-Reco M_{ll#gamma} HighR9; #Delta M_{ll};N_{evts}", 40, -20., 20., eventWeight);
      hmHiggs->fill1DHist(R9, "h1_R9HighR9_SUFFIX","R9;R9;Entries",100,0,1,eventWeight);
      hmHiggs->fill1DHist(R9Cor, "h1_R9CorHighR9Cor_SUFFIX","R9Cor;R9Cor;Entries",100,0,1,eventWeight);
    } else {
      hmHiggs->fill1DHist((ZP4+GP4).M(),"h1_InvariantMassReco1GevLowR9FULLRANGE_SUFFIX","Invariant Mass (H->Z#gamma);Mass (GeV);Entries",90,90,190,eventWeight);
      if (genHs.size() > 0) hmHiggs->fill1DHist(genHs[0].M()-(ZP4+GP4).M(),"h1_genHiggsMassResLowR9_SUFFIX", "Gen-Reco M_{ll#gamma} LowR9; #Delta M_{ll};N_{evts}", 40, -20., 20., eventWeight);
      hmHiggs->fill1DHist(R9, "h1_R9LowR9_SUFFIX","R9;R9;Entries",100,0,1,eventWeight);
      hmHiggs->fill1DHist(R9Cor, "h1_R9CorLowR9Cor_SUFFIX","R9Cor;R9Cor;Entries",100,0,1,eventWeight);
    }

    if (dataDumps && isRealData){
      //if(engCor) DataDumper(&lepton1, &lepton2, &GP4, &uncorLepton1, &uncorLepton2, &uncorGP4, R9Cor, GP4scEta, dataDump,SCetaEl1,SCetaEl2);
      DataDumper(&lepton1, &lepton2, &GP4, R9Cor, GP4scEta, dataDump,SCetaEl1,SCetaEl2);
    }

    if (dumps){
      FinalDumper(&lepton1, &lepton2, &GP4, catNum, finalDump);
    }

    unBinnedWeight = eventWeight;



    hm->fill1DHist(R9, "h1_R9Full_SUFFIX","R9;R9;Entries",100,0,1,eventWeight);
    hm->fill1DHist(R9Cor, "h1_R9CorFull_SUFFIX","R9Cor;R9Cor;Entries",100,0,1,eventWeight);
    hm->fill1DHist(lepton1.DeltaR(GP4),"h1_DeltaRLeading_SUFFIX","DeltaR leading lepton vs photon;#Delta R;Entries",16,0,4,eventWeight,"ZGamma");
    hm->fill1DHist(lepton2.DeltaR(GP4),"h1_DeltaRTrailing_SUFFIX","DeltaR trailing lepton vs photon;#Delta R;Entries",16,0,4,eventWeight,"ZGamma");
    hm->fill1DHist(CalculateM12sqrd(ZP4,GP4),"h1_M12sqrd_SUFFIX","M12^{2} CORRECTED (H->Z#gamma);Mass^{2} (GeV^{2});Entries",25,3600,34400,eventWeight,"ZGamma");
    hm->fill1DHist(CalculateX1(ZP4,GP4),"h1_X1_SUFFIX","x1;x1;Entries",50,0,0.4,eventWeight,"ZGamma");
    hm->fill1DHist(CalculateX2(ZP4,GP4),"h1_X2_SUFFIX","x2;x2;Entries",50,0,0.4,eventWeight,"ZGamma");
    hm->fill1DHist(ZP4.Pt()-GP4.Pt(),"h1_PtDiffZG_SUFFIX","Zp_{T}-Gp_{T} (scalar);#Delta p_{T} (GeV);Entries",50,-100,100,eventWeight,"ZGamma");
    hm->fill1DHist((ZP4-GP4).Pt(),"h1_PtVecDiffZG_SUFFIX","Zp_{T}-Gp_{T} (vector);#Delta p_{T} (GeV);Entries",50,0,200,eventWeight,"ZGamma");
    hm->fill1DHist(ZP4.Pt()/GP4.Pt(),"h1_PtRatZG_SUFFIX","Zp_{T}/Gp_{T};#Ratio p_{T} (GeV);Entries",20,0,10,eventWeight,"ZGamma");
    hm->fill1DHist((ZP4+GP4).Pt(),"h1_PtSumZG_SUFFIX","3-Body p_{T};p_{T};Entries",25,0,70,eventWeight,"ZGamma");
    hm->fill1DHist(GP4.E(),"h1_GammaEnergy_SUFFIX","Gamma Energy;E (GeV);Entries",30,0,150,eventWeight,"ZGamma");
    hm->fill1DHist(GP4.Pt(),"h1_GammaPt_SUFFIX","Gamma p_{T};p_{T} (GeV);Entries",12,0,60,eventWeight,"ZGamma");
    hm->fill2DHist(CalculateX1(ZP4,GP4),CalculateX2(ZP4,GP4),"h2_X1X2_SUFFIX","x1 vs x2;x1;x2",50,0,0.4,50,0,0.4,eventWeight,"ZGamma");
    hm->fill2DHist(ZP4.M(),(ZP4+GP4).M(),"h2_InvariantMasses_SUFFIX","2 Body vs 3 Body Invariant Mass; Z (GeV); Z#gamma (GeV)",60,60,120,100,60,190,eventWeight,"ZGamma");
    hm->fill2DHist((ZP4+GP4).M(),GP4.Pt(),"h2_3BodyVsGPt_SUFFIX","3-Body Mass vs #gamma p_{T};Z#gamma;#gamma p_{T}",20,100,190,20,0,60,eventWeight,"ZGamma");
    hm->fill2DHist((ZP4+GP4).M(),(ZP4.Pt()-GP4.Pt()),"h2_3BodyVsDiffPt_SUFFIX","3-Body Mass vs diff p_{T};Z#gamma M;Z p_{T} - #gamma p_{T}",50,100,190,50,-100,100,eventWeight,"ZGamma");
    hm->fill2DHist((ZP4+GP4).M(),(ZP4.Pt()/GP4.Pt()),"h2_3BodyVsRatPt_SUFFIX","3-Body Mass vs ratio p_{T};Z#gamma M;Ratio p_{T}",20,100,190,20,0,10,eventWeight,"ZGamma");
    hm->fill2DHist(ZP4.M(),GP4.Pt(),"h2_2BodyVsGPt_SUFFIX","Dilepton Mass vs #gamma p_{T};Z#gamma;#gamma p_{T}",20,80,105,20,0,60,eventWeight,"ZGamma");
    hm->fill2DHist((ZP4+GP4).M(),ZP4.Pt(),"h2_3BodyVsZPt_SUFFIX","3-Body Mass vs #Z p_{T};Z#gamma;Z p_{T}",20,100,190,20,0,60,eventWeight,"ZGamma");
    hm->fillProfile(CalculateX1(ZP4,GP4),CalculateX2(ZP4,GP4),"profile_X1X2Pro_SUFFIX","x1 vs x2;x1;x2",50,0,0.4,0,0.4,eventWeight,"ZGamma");
    hm->fill1DHist(lepton1.Pt()/lepton2.Pt(),"h1_MuonRatPt_SUFFIX","leading/trailing p_{T};p_{T} Ratio;Entries",15,0,5,eventWeight,"ZGamma");
    hm->fill1DHist(photonsIDIso.size(),"h1_photonMult_SUFFIX","Photon Multiplicity;nPhoton;Entries",10,0.5,10.5,eventWeight,"ZGamma");
    hm->fill1DHist(fabs(ZP4.DeltaPhi(GP4)),"h1_DeltaPhiZG_SUFFIX","#Delta#phi (Z,#gamma);#Delta#phi (rad);Entries",20,0,TMath::Pi()+0.5,eventWeight,"ZGamma");
    if (suffix.find("GammaStar") !=string::npos){
      hm->fill1DHist(ZP4.M(),"h1_gammaStarMll1_SUFFIX","gammaStar M_ll full range; M_ll (GeV);Entries",120,-5,115,eventWeight,"ZGamma");
      hm->fill1DHist(ZP4.M(),"h1_gammaStarMll2_SUFFIX","gammaStar M_ll low range; M_ll (GeV);Entries",50,-5,45,eventWeight,"ZGamma");
    }

    diffZGscalar    = ZP4.Pt()-GP4.Pt();
    diffZGvector    = (ZP4-GP4).Pt();
    threeBodyMass   = (ZP4+GP4).M();
    threeBodyPt     = (ZP4+GP4).Pt();
    divPt           = ZP4.Pt()/GP4.Pt();
    cosZ            = cospolarZBoost;
    cosG            = cospolarGBoost;
    GPt             = GP4.Pt();
    ZPt             = ZP4.Pt();
    DPhi            = fabs(ZP4.DeltaPhi(GP4));
    diffPlaneMVA    = diffPlane;
    dr1             = lepton1.DeltaR(GP4);
    dr2             = lepton2.DeltaR(GP4);
    M12             = CalculateM12sqrd(ZP4,GP4);
    scaleFactor     = eventWeight;
    if (suffix == "HZG125Signal") scaleFactor *= 4.98/124532.73;
    if (suffix == "ZJets") scaleFactor *= (4.98*2.5)/11.900;
    if (suffix == "ZGMuMu") scaleFactor *= 4.98/7.11;
    if (suffix == "ZGEE") scaleFactor *= 4.98/7.11;
    if (suffix == "ZZJets") scaleFactor *= 4.98/6175;
    if (suffix == "WZJets") scaleFactor *= 4.98/1426.5;
    if (suffix == "WWJets") scaleFactor *= 4.98/250.4;

    //// angular correlation parts /////
    hm->fill1DHist(polarZBoost,"h1_polarZBoost_SUFFIX","#theta of Z in Higgs CoM;#theta (rad);Entries",15,0,TMath::Pi(),eventWeight,"ZGamma");
    hm->fill1DHist(polarGBoost,"h1_polarGBoost_SUFFIX","#theta of #gamma in Higgs CoM;#theta (rad);Entries",15,0,TMath::Pi(),eventWeight,"ZGamma");
    hm->fill1DHist(cospolarZBoost,"h1_cospolarZBoost_SUFFIX","cos(#theta) of Z in Higgs CoM;cos(#theta);Entries",10,-1,1,eventWeight,"ZGamma");
    hm->fill1DHist(cospolarGBoost,"h1_cospolarGBoost_SUFFIX","cos(#theta) of #gamma in Higgs CoM;cos(#theta);Entries",10,-1,1,eventWeight,"ZGamma");
    hm->fill1DHist(angleGLepNeg,"h1_angleGLepNeg_SUFFIX","angle between #gamma and neg lepton, Higgs CoM;#Delta angle (rad);Entries",15,0,TMath::Pi(),eventWeight,"ZGamma");
    hm->fill1DHist(angleGLepPos,"h1_angleGLepPos_SUFFIX","angle between #gamma and pos lepton, Higgs CoM;#Delta angle (rad);Entries",15,0,TMath::Pi(),eventWeight,"ZGamma");
    hm->fill1DHist(azimuthLepPlane,"h1_azimuthLepPlane_SUFFIX","#phi of diLepton plane in Higgs CoM;#phi (rad);Entries",20,-TMath::Pi(),TMath::Pi(),eventWeight,"ZGamma");
    hm->fill1DHist(diffPlane,"h1_diffPlane_SUFFIX","angle between Z#gamma plane and diLepton plane, Higgs CoM;#Delta angle (rad);Entries",15,0,TMath::Pi(),eventWeight,"ZGamma");


  }


  /////////////////////////
  // Fill std histograms //
  /////////////////////////

  StandardPlots(lepton1,lepton2,GP4,eventWeight,"", "pT-Eta-Phi");

  //////////////////////////////
  // Fill Vtx variable histos //
  //////////////////////////////


  hm->fill1DHist(primaryVtx->GetSize(),"h1_pvMultCopy_SUFFIX", "Multiplicity of PVs", 50, 0.5, 50.5, eventWeight,"Vtx");
  if (isRealData) {
    hm->fillProfile(runNumber,primaryVtx->GetSize(),"p1_nVtcsCopy", "Average number of vertices per run; Run Number; nVertices", 8700.0, 135000.0, 144200.0, 0.0, 6.0, 1,"Vtx");
  }

  //////////////////////////////
  // Fill Gen Particle Histos //
  //////////////////////////////


  GenPlots(genZs,genMuons,genPhotons,genHs,ZP4,GP4,eventWeight);

  //////////
  // misc //
  //////////

  hm->fill1DHist(eventWeight,"h1_eventWeight_SUFFIX", "event weight", 100, 0., 2.,1,"Misc");
  hm->fill1DHist(primaryVtx->GetSize(),"h1_pvMult_SUFFIX", "Multiplicity of PVs", 25, 0.5, 25.5, eventWeight,"Misc");
  if (!isRealData) hm->fill1DHist(ptHat,"h1_ptHat_SUFFIX","ptHat",37, 15.0, 200.0, eventWeight,"Misc");

  if (isRealData) {
    hm->fill1DHist(runNumber,"h1_goodRuns_SUFFIX", "Number of events passing selection cuts by run;Run number; nEvents", 9200, 160000., 175000.,1,"Misc");
    hm->fillProfile(runNumber,primaryVtx->GetSize(),"p1_nVtcs", "Average number of vertices per run; Run Number; nVertices", 8700.0, 135000.0, 144200.0, 0.0, 6.0, 1,"Misc");
  }
  pvPosition->Delete();

  if (nEvents[0]%2 == 0){
    trainingChain->Fill();
  }else{
    sampleChain->Fill();
  }

  m_llgChain->Fill();

   return kTRUE;
}

void higgsAnalyzer::Terminate()
{

  TH1F* eventHisto = (TH1F*)histoFile->GetDirectory("Misc")->Get("h1_acceptanceByCut_SUFFIX");
  

  cout<<"\nRunning over "<<suffix<<" dataset with "<<selection<<" selection."<<"\n"<<endl;
  cout << "| CUT DESCRIPTION                    |\t" << "\t|"                         << endl;
  cout << "| Unskimmed events:                  |\t" << unskimmedEventsTotal          << "\t|" << endl;
  cout << "| Initial number of events:          |\t" << nEvents[0]                    << "\t|" << endl;
  //cout << "| Electrion Preselection:            |\t" << nEvents[1]                    << "\t|" << endl;
  cout << "| Pass HLT selection:                |\t" << nEvents[2]                    << "\t|" << endl;
  cout << "| Good PV:                           |\t" << nEvents[3]                    << "\t|" << endl;
  cout << "| Data quality bits:                 |\t" << nEvents[4]                    << "\t|" << endl;
  cout << "| Two opp charged leptons (ID):      |\t" << nEvents[5]                    << "\t|" << endl;
  cout << "| Two opp charged leptons (ID/Iso):  |\t" << nEvents[6]                    << "\t|" << endl;
  cout << "| Good Z:                            |\t" << nEvents[7]                    << "\t|" << endl;
  cout << "| Z mass window:                     |\t" << nEvents[8]                    << "\t|" << endl;
  cout << "| Gamma Acc and Energy:              |\t" << nEvents[9]                    << "\t|" << endl;
  cout << "| Photon DR Cuts:                    |\t" << nEvents[10]                   << "\t|" << endl;
  //cout << "| Third lepton veto:                 |\t" << nEvents[11]                   << "\t|" << endl;
  //cout << "| b-jet veto:                        |\t" << nEvents[12]                   << "\t|" << endl;
  //cout << "| MET:                               |\t" << nEvents[13]                   << "\t|" << endl;
  //cout << "| MET/QT:                            |\t" << nEvents[14]                   << "\t|" << endl;
  //cout << "| Jet Multiplicity:                  |\t" << nEvents[15]                   << "\t|" << endl;
  cout << "| Higgs Mass Low:                    |\t" << nEvents[16]                   << "\t|" << endl;
  cout << "| Higgs Mass High:                   |\t" << nEvents[17]                   << "\t|" << endl;
  cout << "| MZ + MZG Cut:                      |\t" << nEvents[18]                   << "\t|" << endl;
  //cout << "| Angular Cut:                       |\t" << nEvents[19]                   << "\t|" << endl;
  //cout << "| M12 Cut:                           |\t" << nEvents[20]                   << "\t|" << endl;
  //cout << "| MVA Cuts:                          |\t" << nEvents[21]                   << "\t|" << endl;
  cout << "| VBF Cuts:                          |\t" << nEvents[22]                   << "\t|" << endl;
  cout << "| 100-190 range:                     |\t" << nEvents[51]                   << "\t|" << endl;
  cout << "| 100-190 weighted (all):            |\t" << eventHisto->Integral(50,50) << "\t|" << endl;
  cout << "| 100-190 weighted (PU):             |\t" << eventHisto->Integral(51,51) << "\t|" << endl;
  cout << "| 100-190 weighted (pho):            |\t" << eventHisto->Integral(52,52) << "\t|" << endl;
  cout << "| 100-190 weighted (trig):           |\t" << eventHisto->Integral(53,53) << "\t|" << endl;
  cout << "| 100-190 weighted (lep):            |\t" << eventHisto->Integral(54,54) << "\t|" << endl;
  cout << "| CAT 1:                             |\t" << nEvents[61]                   << "\t|" << endl;
  cout << "| CAT 2:                             |\t" << nEvents[62]                   << "\t|" << endl;
  cout << "| CAT 3:                             |\t" << nEvents[63]                   << "\t|" << endl;
  cout << "| CAT 4:                             |\t" << nEvents[64]                   << "\t|" << endl;
  cout << "| GEN ACCEPTANCE Leptons:            |\t" << genAccept[0]                  << "\t|" << endl;
  cout << "| GEN ACCEPTANCE Total:              |\t" << genAccept[1]                  << "\t|" << endl;

  //hm->writeHists();
  trainingFile->Write();
  sampleFile->Write();
  histoFile->Write();
  higgsFile->Write();
  eleSmearFile->Write();
  eleIDISOFile->Write();
  m_llgFile->Write();
  trainingFile->Close();
  sampleFile->Close();
  histoFile->Close();  
  higgsFile->Close();
  eleSmearFile->Close();
  eleIDISOFile->Close();
  m_llgFile->Close();

  if(electronDump){
    elDump2.close();
    elDumpFinal.close();
  }
  if(muonDump){
    muDump1.close();
    muDumpFinal.close();
  }
  if(dumps){
    phDump1.close();
    phDump2.close();
    finalDump.close();
  }

  if(dataDumps && isRealData){
    dataDump.close();
  }
}

////////////////////////////
// User Defined Functions //
////////////////////////////

void higgsAnalyzer::MetPlusZPlots(TLorentzVector metP4, TLorentzVector ZP4, float eventWeight)
{

  hm->fill1DHist(metP4.Pt()/ZP4.Pt(),"h1_MetOverQt_SUFFIX", "MET/Q_{T,ll};MET/q_{T};N_{evts}", 45, 0., 9., eventWeight,"MET+Lepton");
  hm->fill1DHist((metP4 + ZP4).Pt(),"h1_MetPlusQtMagnitude_SUFFIX", "|MET + q_{T}|;|MET + Q_{T,ll}|;N_{evts}", 60, 0., 300., eventWeight,"MET+Lepton");
  hm->fill1DHist(fabs(metP4.DeltaPhi(ZP4)),"h1_MetQtDeltaPhi_SUFFIX", "#Delta#phi(q_{T}, MET);#Delta#phi;N_{evts}", 36, 0., TMath::Pi(), eventWeight,"MET+Lepton");
  hm->fill2DHist(metP4.Pt(), metP4.Pt()/ZP4.Pt(),"h2_MetByMetOverQt_SUFFIX", "MET/q_{T} vs. MET; MET; MET/q_{T}", 40, 0., 200., 40, 0., 4., eventWeight,"2D");
  hm->fill2DHist(metP4.Pt(), (metP4 + ZP4).Pt(),"h2_MetByMetPlusQtMag_SUFFIX", "|MET + q_{T}| vs. MET; MET; MET/q_{T}", 40, 0., 200., 60, 0., 300., eventWeight,"2D");

  //MT
  hm->fill1DHist(CalculateTransMass(metP4, ZP4),"h1_MetDileptonMT_SUFFIX", "M_{T};M_{T};N_{evts}", 90, 50., 500., eventWeight,"MET+Lepton");

  //ProjectedMET
  hm->fill1DHist(metP4.Pt()*cos(metP4.DeltaPhi(ZP4)),"h1_ProjMetByQt_SUFFIX", "Longitudinal MET by q_{T};MET_{longitudinal};N_{evts}", 52, -190., 80., eventWeight,"MET+Lepton");
  hm->fill1DHist(metP4.Pt()*sin(metP4.DeltaPhi(ZP4)),"h1_OrthoMetByQt_SUFFIX", "Transverse MET by q_{T};MET_{Transverse};N_{evts}", 64, -100., 80., eventWeight,"MET+Lepton");
  if (fabs(metP4.DeltaPhi(ZP4)) < TMath::Pi()/2) {
    hm->fill1DHist(metP4.Pt() * sin(fabs(metP4.DeltaPhi(ZP4))),"h1_ProjectedMet_SUFFIX", "ProjMET; ProjMET; N_{evts}", 50, 0., 250., eventWeight,"MET");
  } else {
    hm->fill1DHist(metP4.Pt(),"h1_ProjectedMet_SUFFIX", "ProjMET; ProjMET; N_{evts}", 50, 0., 250., eventWeight,"MET");
  }
}

void higgsAnalyzer::MetPlusLeptonPlots(TLorentzVector metP4, TLorentzVector p1, TLorentzVector p2, float eventWeight)
{
  hm->fill1DHist(CalculateTransMassAlt(metP4, p1),"h1_MetLeadLeptonMT_SUFFIX", "M_{T,l1};M_{T,l1};N_{evts}", 50, 0., 250., eventWeight,"MET+Lepton");
  hm->fill1DHist(CalculateTransMassAlt(metP4, p1),"h1_MetTrailingLeptonMT_SUFFIX", "M_{T,l2};M_{T,l2};N_{evts}", 50, 0., 250., eventWeight,"MET+Lepton");
  hm->fill2DHist(CalculateTransMassAlt(metP4, p1), CalculateTransMassAlt(metP4, p2),"h2_MetLeptonMT_SUFFIX", "MT(MET,lepton);MT(MET, l1);MT(MET, l2)", 80, 0., 400., 80, 0., 400., eventWeight,"2D");
  hm->fill1DHist(fabs(metP4.DeltaPhi(p1)),"h1_MetLeadLeptonDeltaPhi_SUFFIX", "#Delta#phi(l1, MET);#Delta#phi;N_{evts}", 36, 0., TMath::Pi(), eventWeight,"MET+Lepton");
  hm->fill1DHist(fabs(metP4.DeltaPhi(p2)),"h1_MetTrailingLeptonDeltaPhi_SUFFIX", "#Delta#phi(l2, MET);#Delta#phi;N_{evts}", 36, 0., TMath::Pi(), eventWeight,"MET+Lepton");
}

void higgsAnalyzer::LeptonBasicPlots(TLorentzVector p1, TLorentzVector p2, float eventWeight)
{
  hm->fill1DHist(p1.Pt(),"h1_leadLeptonPt_SUFFIX", "p_{T} leading lepton;p_{T};N_{evts}", 48, 10., 250., eventWeight,"Lepton");     
  hm->fill1DHist(p1.Eta(),"h1_leadLeptonEta_SUFFIX", "#eta leading lepton;#eta;N_{evts}", 25, -2.5, 2.5, eventWeight,"Lepton");    
  hm->fill1DHist(p1.Phi(),"h1_leadLeptonPhi_SUFFIX", "#phi leading lepton;#phi;N_{evts}", 36, -TMath::Pi(), TMath::Pi(), eventWeight,"Lepton");    

  hm->fill1DHist(p2.Pt(),"h1_trailingLeptonPt_SUFFIX", "p_{T} trailing lepton;p_{T};N_{evts}", 38, 10., 200., eventWeight,"Lepton"); 
  hm->fill1DHist(p2.Eta(),"h1_trailingLeptonEta_SUFFIX", "#eta trailing lepton;#eta;N_{evts}", 25, -2.5, 2.5, eventWeight,"Lepton");
  hm->fill1DHist(p2.Phi(),"h1_trailingLeptonPhi_SUFFIX", "#phi trailing lepton;#phi;N_{evts}", 36, -TMath::Pi(), TMath::Pi(), eventWeight,"Lepton");

  hm->fill2DHist(p2.Pt(),fabs(p2.Eta()),"h2_trailingLeptonPtEta_SUFFIX","Trailing Lepton p_{T} vs. #||{#eta}",5,0,100,5,0,2.5,eventWeight,"Lepton");

  hm->fill1DHist(p2.Pt()/p1.Pt(),"h1_diLeptonPtRatio_SUFFIX", "dilepton p_{T} ratio;p_{T,2}/p_{T,1};N_{evts}", 25, 0., 1., eventWeight,"Lepton"); 
  hm->fill1DHist(fabs(p2.Eta() - p1.Eta()),"h1_diLeptonDeltaEta_SUFFIX", "dilepton #Delta#eta;#Delta#eta;N_{evts}", 40, 0., 8., eventWeight,"Lepton");
  hm->fill1DHist(fabs(p2.DeltaPhi(p1)),"h1_diLeptonDeltaPhi_SUFFIX", "dilepton #Delta#phi;#Delta#phi;N_{evts}", 18, 0., TMath::Pi(), eventWeight,"Lepton");
  hm->fill1DHist(p2.DeltaR(p1),"h1_diLeptonDeltaR_SUFFIX", "dilepton #Delta R;#Delta R;N_{evts}", 18, 0., 4.5, eventWeight,"Lepton");  
}

void higgsAnalyzer::StandardPlots(TLorentzVector p1, TLorentzVector p2, TLorentzVector gamma, float eventWeight,string tag, string folder)
{
  TLorentzVector diLep = p1+p2;
  TLorentzVector threeBody = p1+p2+gamma;

  hm->fill1DHist(p1.Pt(),"h1_leadLeptonStdPt"+tag+"_SUFFIX", "p_{T} leading lepton;p_{T};N_{evts}", 28, 0., 140., eventWeight,folder);     
  hm->fill1DHist(p1.Eta(),"h1_leadLeptonStdEta"+tag+"_SUFFIX", "#eta leading lepton;#eta;N_{evts}", 20, -2.5, 2.5, eventWeight,folder);    
  hm->fill1DHist(p1.Phi(),"h1_leadLeptonStdPhi"+tag+"_SUFFIX", "#phi leading lepton;#phi;N_{evts}", 20, -TMath::Pi(), TMath::Pi(), eventWeight,folder);    

  hm->fill1DHist(p2.Pt(),"h1_trailingLeptonStdPt"+tag+"_SUFFIX", "p_{T} trailing lepton;p_{T};N_{evts}", 20, 0., 100., eventWeight,folder); 
  hm->fill1DHist(p2.Eta(),"h1_trailingLeptonStdEta"+tag+"_SUFFIX", "#eta trailing lepton;#eta;N_{evts}", 20, -2.5, 2.5, eventWeight,folder);
  hm->fill1DHist(p2.Phi(),"h1_trailingLeptonStdPhi"+tag+"_SUFFIX", "#phi trailing lepton;#phi;N_{evts}", 20, -TMath::Pi(), TMath::Pi(), eventWeight,folder);

  hm->fill1DHist(gamma.Pt(),"h1_photonPt"+tag+"_SUFFIX", "p_{T} gamma;p_{T};N_{evts}", 16, 0., 80., eventWeight,folder); 
  hm->fill1DHist(gamma.Eta(),"h1_photonEta"+tag+"_SUFFIX", "#eta gamma;#eta;N_{evts}", 20, -2.5, 2.5, eventWeight,folder);//20
  hm->fill1DHist(gamma.Phi(),"h1_photonPhi"+tag+"_SUFFIX", "#phi gamma;#phi;N_{evts}", 20, -TMath::Pi(), TMath::Pi(), eventWeight,folder);//18

  hm->fill1DHist(diLep.Pt(),"h1_diLepPt"+tag+"_SUFFIX", "p_{T} Z;p_{T};N_{evts}", 20, 0., 100., eventWeight,folder);     
  hm->fill1DHist(diLep.Eta(),"h1_diLepEta"+tag+"_SUFFIX", "#eta Z;#eta;N_{evts}", 18, -2.5, 2.5, eventWeight,folder);    
  hm->fill1DHist(diLep.Phi(),"h1_diLepPhi"+tag+"_SUFFIX", "#phi Z;#phi;N_{evts}", 18, -TMath::Pi(), TMath::Pi(), eventWeight,folder);    
  hm->fill1DHist(diLep.M(),"h1_diLepMass"+tag+"_SUFFIX", "M_{Z};M (GeV);N_{evts}", 25, 66, 116, eventWeight,folder);    

  hm->fill1DHist(threeBody.Pt(),"h1_threeBodyPt"+tag+"_SUFFIX", "p_{T} 3body;p_{T};N_{evts}", 20, 0., 100., eventWeight,folder);     
  hm->fill1DHist(threeBody.Eta(),"h1_threeBodyEta"+tag+"_SUFFIX", "#eta 3body;#eta;N_{evts}", 18, -2.5, 2.5, eventWeight,folder);    
  hm->fill1DHist(threeBody.Phi(),"h1_threeBodyPhi"+tag+"_SUFFIX", "#phi 3body;#phi;N_{evts}", 18, -TMath::Pi(), TMath::Pi(), eventWeight,folder);    
  hm->fill1DHist(threeBody.M(),"h1_threeBodyMass"+tag+"_SUFFIX", "M_{3body};M (GeV);N_{evts}", 75, 50, 200, eventWeight,folder);    


  hm->fill1DHist(primaryVtx->GetSize(),"h1_pvMult"+tag+"_SUFFIX", "Multiplicity of PVs;N_{PV};N_{evts}", 25, 0.5, 25.5, eventWeight, folder);
  for(int i = 0; i < 64; ++i) {
    unsigned long iHLT = 0x0; 
    iHLT = 0x01 << i;  
    if ((triggerStatus & iHLT) == iHLT) hm->fill1DHist(i+1,"h1_triggerStatus"+tag+"_SUFFIX", "Triggers (no weight)", 64, 0.5, 64.5,1,"Misc");  
  } 
}

void higgsAnalyzer::AnglePlots(ZGAngles &zga,float eventWeight)
{
  hm->fill1DHist(zga.costheta_lp,"h1_costhetaLP_SUFFIX", "Cos(#theta) positive lepton;cos(#theta);N_{evts}", 20, -1., 1., eventWeight,"ZGAngles");     
  hm->fill1DHist(zga.costheta_lm,"h1_costhetaLM_SUFFIX", "Cos(#theta) negative lepton;cos(#theta);N_{evts}", 20, -1., 1., eventWeight,"ZGAngles");     
  hm->fill1DHist(zga.phi,"h1_phi_SUFFIX", "#phi positive lepton;#phi;N_{evts}", 20, -TMath::Pi(), TMath::Pi(), eventWeight,"ZGAngles");     
  hm->fill1DHist(zga.cosTheta,"h1_costhetaZG_SUFFIX", "Cos(#Theta) ZG system;cos(#Theta);N_{evts}", 20, -1., 1., eventWeight,"ZGAngles");     
}

void higgsAnalyzer::DileptonBasicPlots(TLorentzVector ZP4, float eventWeight)
{
  hm->fill1DHist(ZP4.Mt(),"h1_diLeptonTransMass_SUFFIX", "M_{T,ll};M_{T,ll};N_{evts}", 100, 55., 255., eventWeight,"Lepton");
  hm->fill1DHist(ZP4.M(),"h1_diLeptonMass_SUFFIX", "M_{ll}; M_{ll};N_{evts}", 40, 70., 110., eventWeight,"Lepton");     
  hm->fill1DHist(ZP4.Pt(),"h1_diLeptonQt_SUFFIX", "q_{T};Q_{T};N_{evts}", 50, 0., 500., eventWeight,"Lepton");       
}

void higgsAnalyzer::GenPlots(vector<TCGenParticle> Zs, vector<TCGenParticle> leps, vector<TCGenParticle> phots, vector<TCGenParticle> Hs, TLorentzVector ZP4, TLorentzVector GP4, float eventWeight)
{
  //cout<<"gen function has activated"<<endl;
  if (Zs.size()>0){
    //cout<<"yarp there's some Zs"<<endl;
    hm->fill1DHist(Zs[0].M(),"h1_genZMass_SUFFIX", "GEN M_{Z}; M_{Z};N_{evts}", 40, 70., 110., eventWeight,"GenLvl");     
    hm->fill1DHist(Zs[0].Pt(),"h1_genZQt_SUFFIX", "q_{T};Q_{T};N_{evts}", 50, 0., 500., eventWeight,"GenLvl");       
  }
  if (Hs.size()>0){
    //cout<<"yarp there's some Zs"<<endl;
    hmHiggs->fill1DHist(Hs[0].M(),"h1_genHMass_SUFFIX", "GEN M_{H}; M_{H};N_{evts}", 90, 90., 190., eventWeight);     
    hm->fill1DHist(Hs[0].M(),"h1_genHMass_SUFFIX", "GEN M_{H}; M_{H};N_{evts}", 90, 90., 190., eventWeight,"GenLvl");     
    hm->fill1DHist(Hs[0].Pt(),"h1_genHQt_SUFFIX", "H p_{T};Q_{T};N_{evts}", 50, 0., 500., eventWeight,"GenLvl");       
  }

  vector<TCGenParticle>::iterator lepIt;
  TCGenParticle zMuon1;
  TCGenParticle zMuon2;
  bool oneSet = false;

  if (leps.size() > 1){
    for (lepIt=leps.begin(); lepIt<leps.end(); lepIt++){
      if (oneSet && abs((*lepIt).Mother()) == 23){
        zMuon2 = *lepIt;
        break;
      }
      if (!oneSet && abs((*lepIt).Mother()) == 23){
        zMuon1 = *lepIt;
        oneSet = true;
      }
    }
    //cout<<"yarp there's some MUOOOOONs"<<endl;
    hm->fill1DHist((zMuon1+zMuon2).M(),"h1_genDiLeptonMass_SUFFIX", "M_{ll} (pre FSR); M_{ll};N_{evts}", 40, 70., 110., eventWeight,"GenLvl");
    hm->fill1DHist((zMuon1+zMuon2).Pt(),"h1_genDiLeptonPt_SUFFIX", "pT_{ll} (pre FSR); M_{ll};N_{evts}", 50, 0., 500., eventWeight,"GenLvl");
    hm->fill1DHist((zMuon1+zMuon2).M()-ZP4.M(),"h1_genDiLeptonMassRes_SUFFIX", "Gen-Reco No Gamma M_{ll}; #Delta M_{ll};N_{evts}", 40, -10., 10., eventWeight,"GenLvl");

    /// combinatorix time

    if(suffix == "HZG125Signal"){
      //cout<<"SIGNAL GEN SHIT"<<endl;
      TCGenParticle hPhot;
      vector<TCGenParticle>::iterator photIt;
      for (photIt=leps.begin(); photIt<leps.end(); photIt++){
        //cout<<(*photIt).Mother()<<endl;
        if ((*photIt).Mother()==23){
          hPhot = *photIt;
          break;
        }
      }
      //cout<<"End PhotoLoop"<<endl;
      //cout<<"Gen: "<<(Hs[0]).M()<<endl;
      //cout<<"Reco: "<<(ZP4+GP4).M()<<endl;
      hm->fill1DHist((Hs[0]).M()-(ZP4+GP4).M(),"h1_genGamDiLeptonMassRes_SUFFIX", "Gen-Reco M_{ll#gamma} Total Acc; #Delta M_{ll};N_{evts}", 40, -20., 20., eventWeight,"GenLvl");
      // gamma barrel, Z barrel
      if (fabs(GP4.Eta())<1.4442 && fabs(ZP4.Eta())<1.4442){
        hm->fill1DHist((Hs[0]).M()-(ZP4+GP4).M(),"h1_genGamDiLeptonMassResBB_SUFFIX", "Gen-Reco M_{ll#gamma} Barrel Only; #Delta M_{ll};N_{evts}", 40, -20., 20., eventWeight,"GenLvl");
      }
      // gamma barrel, Z endcap
      if (fabs(GP4.Eta())<1.4442 && fabs(ZP4.Eta())>1.566){
        hm->fill1DHist((Hs[0]).M()-(ZP4+GP4).M(),"h1_genGamDiLeptonMassResBE_SUFFIX", "Gen-Reco M_{ll#gamma} #gamma Barrel, diLep Endcap; #Delta M_{ll};N_{evts}", 40, -20., 20., eventWeight,"GenLvl");
      }
      // gamma endcap, Z barrel
      if (fabs(GP4.Eta())>1.566 && fabs(ZP4.Eta())<1.4442){
        hm->fill1DHist((Hs[0]).M()-(ZP4+GP4).M(),"h1_genGamDiLeptonMassResEB_SUFFIX", "Gen-Reco M_{ll#gamma} #gamma Endcap, diLep Barrel; #Delta M_{ll};N_{evts}", 40, -20., 20., eventWeight,"GenLvl");
      }
      // gamma endcap, Z endcap
      if (fabs(GP4.Eta())>1.566 && fabs(ZP4.Eta())>1.566){
        hm->fill1DHist((Hs[0]).M()-(ZP4+GP4).M(),"h1_genGamDiLeptonMassResEE_SUFFIX", "Gen-Reco M_{ll#gamma} Endcap Only; #Delta M_{ll};N_{evts}", 40, -20., 20., eventWeight,"GenLvl");
      }
      // gamma endcap, Z both
      if (fabs(GP4.Eta())>1.566){
        hm->fill1DHist((Hs[0]).M()-(ZP4+GP4).M(),"h1_genGamDiLeptonMassResE_SUFFIX", "Gen-Reco M_{ll#gamma} #gamma Endcap, diLep Both; #Delta M_{ll};N_{evts}", 40, -20., 20., eventWeight,"GenLvl");
      }
      // gamma barrel, Z both
      if (fabs(GP4.Eta())<1.4442){
        hm->fill1DHist((Hs[0]).M()-(ZP4+GP4).M(),"h1_genGamDiLeptonMassResB_SUFFIX", "Gen-Reco M_{ll#gamma} #gamma Barrel, diLep Both; #Delta M_{ll};N_{evts}", 40, -20., 20., eventWeight,"GenLvl");
      }
    }



    else{
      hm->fill1DHist((zMuon1+zMuon2).M()-(ZP4+GP4).M(),"h1_genGamDiLeptonMassRes_SUFFIX", "Gen-Reco M_{ll#gamma} Total Acc; #Delta M_{ll};N_{evts}", 40, -10., 10., eventWeight,"GenLvl");
      // gamma barrel, Z barrel
      if (fabs(GP4.Eta())<1.4442 && fabs(ZP4.Eta())<1.4442){
        hm->fill1DHist((zMuon1+zMuon2).M()-(ZP4+GP4).M(),"h1_genGamDiLeptonMassResBB_SUFFIX", "Gen-Reco M_{ll#gamma} Barrel Only; #Delta M_{ll};N_{evts}", 40, -10., 10., eventWeight,"GenLvl");
      }
      // gamma barrel, Z endcap
      if (fabs(GP4.Eta())<1.4442 && fabs(ZP4.Eta())>1.566){
        hm->fill1DHist((zMuon1+zMuon2).M()-(ZP4+GP4).M(),"h1_genGamDiLeptonMassResBE_SUFFIX", "Gen-Reco M_{ll#gamma} #gamma Barrel, diLep Endcap; #Delta M_{ll};N_{evts}", 40, -10., 10., eventWeight,"GenLvl");
      }
      // gamma endcap, Z barrel
      if (fabs(GP4.Eta())>1.566 && fabs(ZP4.Eta())<1.4442){
        hm->fill1DHist((zMuon1+zMuon2).M()-(ZP4+GP4).M(),"h1_genGamDiLeptonMassResEB_SUFFIX", "Gen-Reco M_{ll#gamma} #gamma Endcap, diLep Barrel; #Delta M_{ll};N_{evts}", 40, -10., 10., eventWeight,"GenLvl");
      }
      // gamma endcap, Z endcap
      if (fabs(GP4.Eta())>1.566 && fabs(ZP4.Eta())>1.566){
        hm->fill1DHist((zMuon1+zMuon2).M()-(ZP4+GP4).M(),"h1_genGamDiLeptonMassResEE_SUFFIX", "Gen-Reco M_{ll#gamma} Endcap Only; #Delta M_{ll};N_{evts}", 40, -10., 10., eventWeight,"GenLvl");
      }
      // gamma endcap, Z both
      if (fabs(GP4.Eta())>1.566){
        hm->fill1DHist((zMuon1+zMuon2).M()-(ZP4+GP4).M(),"h1_genGamDiLeptonMassResE_SUFFIX", "Gen-Reco M_{ll#gamma} #gamma Endcap, diLep Both; #Delta M_{ll};N_{evts}", 40, -10., 10., eventWeight,"GenLvl");
      }
      // gamma barrel, Z both
      if (fabs(GP4.Eta())<1.4442){
        hm->fill1DHist((zMuon1+zMuon2).M()-(ZP4+GP4).M(),"h1_genGamDiLeptonMassResB_SUFFIX", "Gen-Reco M_{ll#gamma} #gamma Barrel, diLep Both; #Delta M_{ll};N_{evts}", 40, -10., 10., eventWeight,"GenLvl");
      }
    }

  }

  if (phots.size() > 1){
    hm->fill1DHist(phots[0].Pt(),"h1_genLeadPhotonPt_SUFFIX", "p_{T} gamma;p_{T};N_{evts}", 16, 0., 80., eventWeight,"GenLvl");
  }

}


float higgsAnalyzer::CalculateTransMass(TLorentzVector p1, TLorentzVector p2)
{
  //float transE    = sqrt(p1.Pt()*p1.Pt() + pow(91.2,2)) + sqrt(p2.Pt()*p2.Pt() + pow(91.2,2)); //p2.M()*p2.M());
  float transE    = sqrt(p1.Pt()*p1.Pt() + p2.M()*p2.M()) + sqrt(p2.Pt()*p2.Pt() +p2.M()*p2.M());
  float transPt   = (p1 + p2).Pt();
  float transMass = sqrt(transE*transE - transPt*transPt);

  return transMass;
}

float higgsAnalyzer::CalculateTransMassAlt(TLorentzVector p1, TLorentzVector p2)
{
  float transMass = sqrt(2*p2.Pt()*p1.Pt() * (1 - cos(fabs(p2.DeltaPhi(p1)))));
  return transMass;
}

bool higgsAnalyzer::CosmicMuonFilter(TCMuon muon1, TCMuon muon2)
{
  float dimuonAngle = muon1.Angle(muon2.Vect());
  if (fabs(TMath::Pi() - dimuonAngle)  < 0.05)
    return true;
  else
    return false;
}


float higgsAnalyzer::Dz(TVector3 objVtx, TLorentzVector objP4, TVector3 vtx)
{
  float vx  = objVtx.x(), vy = objVtx.y(), vz = objVtx.z();
  float px  = objP4.Px(), py = objP4.Py();
  float pz  = objP4.Pz(), pt = objP4.Pt();
  float pvx = vtx.x(), pvy = vtx.y(), pvz = vtx.z();

  float dZ =  (vz-pvz)-((vx-pvx)*px +(vy-pvy)*py)/pt*(pz/pt);
  return dZ;
}

float higgsAnalyzer::Dxy(TVector3 objVtx, TLorentzVector objP4, TVector3 vtx)
{
  float vx  = objVtx.x(), vy = objVtx.y();
  float px  = objP4.Px(), py = objP4.Py(), pt = objP4.Pt();
  float pvx = vtx.x(), pvy = vtx.y();

  double dXY =  (-(vx-pvx)*py + (vy-pvy)*px)/pt;
  return dXY;
}

float higgsAnalyzer::GetPhotonMass()
{
  float photonMass = 91.2;
  if (selection == "eGamma") {
    //        TH1D  *h_mass = (TH1D*)reweightFile.Get("h1_diElectronMass");
    //        photonMass = h_mass->GetRandom();
  }
  if (selection == "muGamma") {
    //        TH1D  *h_mass = (TH1D*)reweightFile.Get("h1_diMuonMass");
    //        photonMass = h_mass->GetRandom();
  }
  return photonMass;
}

//////////////////////////////
// ZGamma Functions !!!!!   //
//////////////////////////////


float higgsAnalyzer::CalculateX2(TLorentzVector p1, TLorentzVector p2)
{
  float Et = sqrt(p1.Pt()*p1.Pt() + p1.M()*p1.M()) + sqrt(p2.Pt()*p2.Pt() + p2.M()*p2.M()); //same points as above.
  float X2 = (Et/7000.)*(exp(-p1.Eta())+exp(-p2.Eta()));
  return X2;
}

float higgsAnalyzer::CalculateX1(TLorentzVector p1,TLorentzVector p2)
{
  float Et = sqrt(p1.Pt()*p1.Pt() + p1.M()*p1.M() + sqrt(p2.Pt()*p2.Pt() + p2.M()*p2.M())); //uses p.M() in both cases - should i change this for the photon, also uses theoretical z mass.
  float X1 = (Et/7000.)*(exp(p1.Eta())+exp(p2.Eta()));
  return X1;
}


float higgsAnalyzer::CalculateM12sqrd(TLorentzVector p1, TLorentzVector p2)
{
  float Et1Et2 =p1.Pt()*p2.Pt();
  float M12sqrd = 2*Et1Et2*(cosh(p1.Eta()-p2.Eta())-cos(p1.DeltaPhi(p2)));
  return M12sqrd;
}

float higgsAnalyzer::Zeppenfeld(TLorentzVector p, TLorentzVector pj1, TLorentzVector pj2)
{
  float zep = p.Eta()-(pj1.Eta()+pj2.Eta())/2.;
  return zep;
}


//////////////////////////////////
// ID Cut and Iso Cut Functions //
//////////////////////////////////


bool higgsAnalyzer::PassMuonID(TCMuon *mu, muIDCuts cutLevel){

  bool muPass = false;

  if (suffix.find("2011") != string::npos){
    if (
        fabs(mu->Eta()) < 2.4
        && mu->IsGLB()                         == cutLevel.IsGLB
        && mu->NormalizedChi2()                < cutLevel.NormalizedChi2
        && mu->NumberOfValidMuonHits()         > cutLevel.NumberOfValidMuonHits
        && mu->NumberOfMatchedStations()       > cutLevel.NumberOfMatchedStations
        && mu->NumberOfValidPixelHits()        > cutLevel.NumberOfValidPixelHits
        && mu->TrackLayersWithMeasurement()    > cutLevel.TrackLayersWithMeasurement
        && fabs(mu->Dxy(pvPosition))           < cutLevel.dxy
        && fabs(mu->Dz(pvPosition))            < cutLevel.dz
       ) muPass = true;
  }else{
    if (
        fabs(mu->Eta()) < 2.4
        && mu->IsPF()                          == cutLevel.IsPF
        && mu->IsGLB()                         == cutLevel.IsGLB
        && mu->NormalizedChi2()                < cutLevel.NormalizedChi2
        && mu->NumberOfValidMuonHits()         > cutLevel.NumberOfValidMuonHits
        && mu->NumberOfMatchedStations()       > cutLevel.NumberOfMatchedStations
        && mu->NumberOfValidPixelHits()        > cutLevel.NumberOfValidPixelHits
        && mu->TrackLayersWithMeasurement()    > cutLevel.TrackLayersWithMeasurement
        && fabs(mu->Dxy(pvPosition))           < cutLevel.dxy
        && fabs(mu->Dz(pvPosition))            < cutLevel.dz
       ) muPass = true;
  }
  return muPass;
}

bool higgsAnalyzer::PassMuonIso(TCMuon *mu, muIsoCuts cutLevel){

  float combIso;

  combIso = (mu->IsoMap("pfChargedHadronPt_R04")
    + max(0.,(double)mu->IsoMap("pfNeutralHadronEt_R04") + mu->IsoMap("pfPhotonEt_R04") - 0.5*mu->IsoMap("pfPUPt_R04")));

  bool isoPass = false;
  if (combIso/mu->Pt() < cutLevel.relCombIso04) isoPass = true;
  return isoPass;
}


bool higgsAnalyzer::PassElectronID(TCElectron *el, elIDCuts cutLevel)
{
  bool elPass = false;
  if (fabs(el->SCEta()) > 2.5) return elPass;
  if (fabs(el->SCEta()) > 1.4442 && fabs(el->SCEta()) < 1.566) return elPass;
  if (
    (fabs(el->Eta()) < 1.566
      && fabs(el->DetaSuperCluster())    < cutLevel.dEtaIn[0]
      && fabs(el->DphiSuperCluster())    < cutLevel.dPhiIn[0]
      && el->SigmaIEtaIEta()             < cutLevel.sigmaIetaIeta[0]
      && el->HadOverEm()                 < cutLevel.HadOverEm[0]
      && fabs(el->Dxy(pvPosition))       < cutLevel.dxy[0]
      && fabs(el->Dz(pvPosition))        < cutLevel.dz[0]
      && el->IdMap("fabsEPDiff")         < cutLevel.fabsEPDiff[0]
      && el->ConversionMissHits()        <= cutLevel.ConversionMissHits[0]
      && el->ConversionVeto()            == cutLevel.PassedConversionProb[0]
      ) ||
    (fabs(el->Eta()) > 1.566  
      && fabs(el->DetaSuperCluster())    < cutLevel.dEtaIn[1]
      && fabs(el->DphiSuperCluster())    < cutLevel.dPhiIn[1]
      && el->SigmaIEtaIEta()             < cutLevel.sigmaIetaIeta[1]
      && el->HadOverEm()                 < cutLevel.HadOverEm[1]
      && fabs(el->Dxy(pvPosition))       < cutLevel.dxy[1]
      && fabs(el->Dz(pvPosition))        < cutLevel.dz[1]
      && el->IdMap("fabsEPDiff")         < cutLevel.fabsEPDiff[1]
      && el->ConversionMissHits()        <= cutLevel.ConversionMissHits[1]
      && el->ConversionVeto()            == cutLevel.PassedConversionProb[1]
      )
      ) elPass = true; 
    //   cout<<"evt: "<<eventNumber<<" muon num: "<<recoMuons->GetSize()<<endl;
       for (int j = 0; j < recoMuons->GetSize(); ++ j)
       {
         TCMuon* thisMuon = (TCMuon*) recoMuons->At(j);    
    //     if (eventNumber==11944 || eventNumber==1780) cout<<thisMuon->DeltaR(*el)<<endl;
         if (thisMuon->DeltaR(*el) < 0.05){
           //cout<<"event: "<<eventNumber<<endl;
           //cout<<"unclean"<<endl;
           //el->Print();
           elPass = false;
           break;
         }
       }
       return elPass;
}

bool higgsAnalyzer::PassElectronIso(TCElectron *el, elIsoCuts cutLevel){
  float thisEA = 0;
  if (fabs(el->Eta())     <  1.0) thisEA = EAEle[0];
  else if (fabs(el->Eta())     <  1.5) thisEA = EAEle[1];
  else if (fabs(el->Eta())     <  2.0) thisEA = EAEle[2];
  else if (fabs(el->Eta())     <  2.2) thisEA = EAEle[3];
  else if (fabs(el->Eta())     <  2.3) thisEA = EAEle[4];
  else if (fabs(el->Eta())     <  2.4) thisEA = EAEle[5];
  else if (fabs(el->Eta())     >  2.4) thisEA = EAEle[6];

  float combIso = (el->IsoMap("pfChIso_R04")
    + max(0.,(double)el->IsoMap("pfNeuIso_R04") + el->IsoMap("pfPhoIso_R04") - rhoFactor*thisEA));
  bool isoPass = false;
  if (combIso/el->Pt() < cutLevel.relCombIso04) isoPass = true;
  return isoPass;
}

bool higgsAnalyzer::PassPhotonID(TCPhoton *ph, phIDCuts cutLevel){
  float tmpEta;
  bool phoPass = false;
  tmpEta = ph->SCEta();
  if (fabs(tmpEta) > 2.5) return phoPass;
  if (fabs(tmpEta) > 1.4442 && fabs(tmpEta) < 1.566) return phoPass;
  if(
      (fabs(tmpEta)  < 1.4442
       && ph->ConversionVeto()       == cutLevel.PassedEleSafeVeto[0]
       && ph->HadOverEm()               < cutLevel.HadOverEm[0]
       && ph->SigmaIEtaIEta()           < cutLevel.sigmaIetaIeta[0]
      ) ||
      (fabs(tmpEta)  > 1.566
       && ph->ConversionVeto()       == cutLevel.PassedEleSafeVeto[1]
       && ph->HadOverEm()               < cutLevel.HadOverEm[1]
       && ph->SigmaIEtaIEta()           < cutLevel.sigmaIetaIeta[1]
      )
    ) phoPass = true;
  return phoPass;
}

bool higgsAnalyzer::PassPhotonIso(TCPhoton *ph, phIsoCuts cutLevel){
  float chEA,nhEA,phEA,chIsoCor,nhIsoCor,phIsoCor,tmpEta;
  bool isoPass = false;
  tmpEta = ph->SCEta();

  if(fabs(tmpEta) > 2.5) return isoPass;

  if (fabs(tmpEta) < 1.0){
    chEA = EAPho[0][0];
    nhEA = EAPho[0][1];
    phEA = EAPho[0][2];
  }else if (fabs(tmpEta) < 1.479){
    chEA = EAPho[1][0];
    nhEA = EAPho[1][1];
    phEA = EAPho[1][2];
  }else if (fabs(tmpEta) < 2.0){
    chEA = EAPho[2][0];
    nhEA = EAPho[2][1];
    phEA = EAPho[2][2];
  }else if (fabs(tmpEta) < 2.2){
    chEA = EAPho[3][0];
    nhEA = EAPho[3][1];
    phEA = EAPho[3][2];
  }else if (fabs(tmpEta) < 2.3){ 
    chEA = EAPho[4][0];
    nhEA = EAPho[4][1];
    phEA = EAPho[4][2];
  }else if (fabs(tmpEta) < 2.4){
    chEA = EAPho[5][0];
    nhEA = EAPho[5][1];
    phEA = EAPho[5][2];
  }else{                                  
    chEA = EAPho[6][0];
    nhEA = EAPho[6][1];
    phEA = EAPho[6][2];
  }

  chIsoCor = ph->IsoMap("chIso03")-rhoFactor*chEA;
  nhIsoCor = ph->IsoMap("nhIso03")-rhoFactor*nhEA;
  phIsoCor = ph->IsoMap("phIso03")-rhoFactor*phEA;

  if (cutLevel.cutName == "loosePhIso"){
    if (
        (fabs(tmpEta) < 1.4442
        //(fabs(ph->Eta())  < 1.566
         && max((double)chIsoCor,0.)          < cutLevel.chIso03[0]
         && max((double)nhIsoCor,0.)          < cutLevel.nhIso03[0] + 0.04*ph->Pt()
         && max((double)phIsoCor,0.)          < cutLevel.phIso03[0] + 0.005*ph->Pt() 
        ) ||
        (fabs(tmpEta) > 1.566
        //(fabs(ph->Eta())  > 1.566
         && max((double)chIsoCor,0.)          < cutLevel.chIso03[1]
         && max((double)nhIsoCor,0.)          < cutLevel.nhIso03[1] + 0.04*ph->Pt() 
         //&& phoCut["phIso03"]/ph->Pt() < nuthin
        )
       ) isoPass = true;
  } else {
    if (
        //(fabs(ph->Eta())  < 1.566
        (fabs(tmpEta) < 1.4442
         && max((double)chIsoCor,0.)          < cutLevel.chIso03[0]
         && max((double)nhIsoCor,0.)          < cutLevel.nhIso03[0] + 0.04*ph->Pt()
         && max((double)phIsoCor,0.)          < cutLevel.phIso03[0] + 0.005*ph->Pt() 
        ) ||
        //(fabs(ph->Eta())  > 1.566
        (fabs(tmpEta) > 1.566
         && max((double)chIsoCor,0.)          < cutLevel.chIso03[1]
         && max((double)nhIsoCor,0.)          < cutLevel.nhIso03[1] + 0.04*ph->Pt() 
         && max((double)phIsoCor,0.)          < cutLevel.phIso03[1] + 0.005*ph->Pt() 
        )
       ) isoPass = true;
  }
  return isoPass;
}
//////////////////////////////
// Debugging Dump Functions //
//////////////////////////////


void higgsAnalyzer::ElectronDump(TCElectron *el, elIDCuts cutLevelID, elIsoCuts cutLevelIso, ofstream & dump)
{
  float thisEA = 0;
  if (fabs(el->Eta())     <  1.0) thisEA = EAEle[0];
  else if (fabs(el->Eta())     <  1.5) thisEA = EAEle[1];
  else if (fabs(el->Eta())     <  2.0) thisEA = EAEle[2];
  else if (fabs(el->Eta())     <  2.2) thisEA = EAEle[3];
  else if (fabs(el->Eta())     <  2.3) thisEA = EAEle[4];
  else if (fabs(el->Eta())     <  2.4) thisEA = EAEle[5];
  else if (fabs(el->Eta())     >  2.4) thisEA = EAEle[6];

  bool idPass = PassElectronID(el, cutLevelID);
  bool isoPass = PassElectronIso(el, cutLevelIso);

  float combIso = (el->IsoMap("pfChIso_R04")
    + max(0.,(double)el->IsoMap("pfNeuIso_R04") + el->IsoMap("pfPhoIso_R04") - rhoFactor*el->IsoMap("EffArea_R04")));

  dump << runNumber << " "                          << eventNumber << " "                                 << el->Pt()
       << " "       << el->Eta()                    << " "         << el->DetaSuperCluster()              << " "      << el->DphiSuperCluster()
       << " "       << el->SigmaIEtaIEta()          << " "         << el->HadOverEm()                     << " "      << el->IdMap("fabsEPDiff")
       << " "       << el->ConversionVeto()         << " "         << el->ConversionMissHits()            << " "      << el->Dxy(pvPosition)
       << " "       << el->Dz(pvPosition)           << " "         << el->IsoMap("pfChIso_R04")           << " "      << el->IsoMap("pfNeuIso_R04")
       << " "       << el->IsoMap("pfPhoIso_R04")   << " "         << combIso                             << " "      << rhoFactor
       << " "       << thisEA                       << " "         << idPass                              << " "      << isoPass
       << endl;
}

void higgsAnalyzer::MuonDump(TCMuon *mu, muIDCuts cutLevelID, muIsoCuts cutLevelIso, ofstream & dump)
{
  float thisEA = 0;
  if (fabs(mu->Eta())     <  1.0) thisEA = EAMu[0];
  else if (fabs(mu->Eta())     <  1.5) thisEA = EAMu[1];
  else if (fabs(mu->Eta())     <  2.0) thisEA = EAMu[2];
  else if (fabs(mu->Eta())     <  2.2) thisEA = EAMu[3];
  else if (fabs(mu->Eta())     <  2.3) thisEA = EAMu[4];
  else if (fabs(mu->Eta())     <  2.4) thisEA = EAMu[5];

  bool idPass = PassMuonID(mu, cutLevelID);
  bool isoPass = PassMuonIso(mu, cutLevelIso);

  float combIso; 
  combIso = (mu->IsoMap("pfChargedHadronPt_R04")
    + max(0.,(double)mu->IsoMap("pfNeutralHadronEt_R04") + mu->IsoMap("pfPhotonEt_R04") - 0.5*mu->IsoMap("pfPUPt_R04")));

  dump << runNumber << " "                              << eventNumber << " "                                 << mu->Pt()
       << " "       << mu->Eta()                        << " "         << mu->IsGLB()                         << " "      << mu->IsPF()
       << " "       << mu->NormalizedChi2()             << " "         << mu->NumberOfValidMuonHits()         << " "      << mu->NumberOfMatchedStations()
       << " "       << mu->Dxy(pvPosition)              << " "         << mu->Dz(pvPosition)                  << " "      << mu->NumberOfValidPixelHits()
       << " "       << mu->TrackLayersWithMeasurement() << " "         << mu->IsoMap("pfChargedHadronPt_R04") << " "      << mu->IsoMap("pfNeutralHadronEt_R04")
       << " "       << mu->IsoMap("pfPhotonEt_R04")     << " "         << combIso                             << " "      << mu->IsoMap("pfPUPt_R04")
       << " "         << idPass                              << " "      << isoPass
       << endl;
}

void higgsAnalyzer::PhotonDump(TCPhoton *ph, phIDCuts cutLevelID, phIsoCuts cutLevelIso, ofstream & dump)
{
  float chEA,nhEA,phEA,tmpEta;
  tmpEta = ph->SCEta();

  if (fabs(tmpEta) < 1.0){
    chEA = EAPho[0][0];
    nhEA = EAPho[0][1];
    phEA = EAPho[0][2];
  }else if (fabs(tmpEta) < 1.479){
    chEA = EAPho[1][0];
    nhEA = EAPho[1][1];
    phEA = EAPho[1][2];
  }else if (fabs(tmpEta) < 2.0){
    chEA = EAPho[2][0];
    nhEA = EAPho[2][1];
    phEA = EAPho[2][2];
  }else if (fabs(tmpEta) < 2.2){
    chEA = EAPho[3][0];
    nhEA = EAPho[3][1];
    phEA = EAPho[3][2];
  }else if (fabs(tmpEta) < 2.3){ 
    chEA = EAPho[4][0];
    nhEA = EAPho[4][1];
    phEA = EAPho[4][2];
  }else if (fabs(tmpEta) < 2.4){
    chEA = EAPho[5][0];
    nhEA = EAPho[5][1];
    phEA = EAPho[5][2];
  }else{                                  
    chEA = EAPho[6][0];
    nhEA = EAPho[6][1];
    phEA = EAPho[6][2];
  }

  bool idPass = PassPhotonID(ph, cutLevelID);
  bool isoPass = PassPhotonIso(ph, cutLevelIso);

  dump << runNumber << " "                   << eventNumber << " "                   << ph->Pt()
       << " "       << tmpEta                << " "         << ph->ConversionVeto()  << " "      << ph->HadOverEm()
       << " "       << ph->SigmaIEtaIEta()   << " "         << ph->IsoMap("chIso03") << " "      << ph->IsoMap("nhIso03")
       << " "       << ph->IsoMap("phIso03") << " "         << rhoFactor             << " "      << chEA
       << " "       << nhEA                  << " "         << phEA                  << " "      << idPass
       << " "       << isoPass
       << endl;
}

void higgsAnalyzer::PhotonDump2(TCPhoton *ph, phIDCuts cutLevelID, phIsoCuts cutLevelIso, TLorentzVector lepton1, TLorentzVector lepton2, ofstream & dump)
{
  float chEA,nhEA,phEA,tmpEta;
  tmpEta = ph->SCEta();

  if (fabs(tmpEta) < 1.0){
    chEA = EAPho[0][0];
    nhEA = EAPho[0][1];
    phEA = EAPho[0][2];
  }else if (fabs(tmpEta) < 1.479){
    chEA = EAPho[1][0];
    nhEA = EAPho[1][1];
    phEA = EAPho[1][2];
  }else if (fabs(tmpEta) < 2.0){
    chEA = EAPho[2][0];
    nhEA = EAPho[2][1];
    phEA = EAPho[2][2];
  }else if (fabs(tmpEta) < 2.2){
    chEA = EAPho[3][0];
    nhEA = EAPho[3][1];
    phEA = EAPho[3][2];
  }else if (fabs(tmpEta) < 2.3){ 
    chEA = EAPho[4][0];
    nhEA = EAPho[4][1];
    phEA = EAPho[4][2];
  }else if (fabs(tmpEta) < 2.4){
    chEA = EAPho[5][0];
    nhEA = EAPho[5][1];
    phEA = EAPho[5][2];
  }else{                                  
    chEA = EAPho[6][0];
    nhEA = EAPho[6][1];
    phEA = EAPho[6][2];
  }

  bool idPass = PassPhotonID(ph, cutLevelID);
  bool isoPass = PassPhotonIso(ph, cutLevelIso);

  float dr1 = ph->DeltaR(lepton1);
  float dr2 = ph->DeltaR(lepton2);

  float scaledPt = ph->Pt()/(*ph+lepton1+lepton2).M();

  dump << runNumber << " "                   << eventNumber << " "                   << ph->Pt()
       << " "       << tmpEta                << " "         << ph->ConversionVeto()  << " "      << ph->HadOverEm()
       << " "       << ph->SigmaIEtaIEta()   << " "         << ph->IsoMap("chIso03") << " "      << ph->IsoMap("nhIso03")
       << " "       << ph->IsoMap("phIso03") << " "         << rhoFactor             << " "      << chEA
       << " "       << nhEA                  << " "         << phEA                  << " "      << idPass
       << " "       << isoPass               << " "         << dr1                   << " "      << dr2
       << " "       << scaledPt              << endl;
}

void higgsAnalyzer::DataDumper(TLorentzVector* lepton1, TLorentzVector* lepton2, TLorentzVector* gamma, float R9input, float SCEta, ofstream & dump, float eta1, float eta2){
  int cat = -1;
  bool goodLep = false;
  if (selection == "mumuGamma"){
    goodLep = (fabs(lepton1->Eta()) < 0.9 || fabs(lepton2->Eta()) < 0.9) && (fabs(lepton1->Eta()) < 2.1 && fabs(lepton2->Eta()) < 2.1);
  }else if (selection == "eeGamma"){
    goodLep = (fabs(eta1) < 1.4442 && fabs(eta2) < 1.4442);
  }

  if (goodLep && (fabs(SCEta) < 1.4442) ){
    if (R9input>=R9Cut) cat = 0;
    else cat = 3;
  }else if (!goodLep && (fabs(SCEta) < 1.4442) ){
    cat = 1;
  } else {
    cat = 2;
  }
  float mll = (*lepton1 + *lepton2).M();
  float mllg = (*lepton1 + *lepton2 + *gamma).M();

  dump << "cat: "     << setw(2)  << cat         << " run: "   << setw(7)  << runNumber      << " evt: "   << setw(10) << eventNumber
       << " lumi: "   << setw(5)  << lumiSection << " mll: "   << setw(10) << mll            << " phopt: " << setw(10) << gamma->Pt()
       << " mllg: "   << setw(10) << mllg        << " l1eta: " << setw(10) << lepton1->Eta() << " l2eta: " << setw(10) << lepton2->Eta()
       << " pSCeta: " << setw(10) << SCEta       << " r9: "    << setw(10) << R9input        << endl;

}

void higgsAnalyzer::DataDumper(TLorentzVector* lepton1, TLorentzVector* lepton2, TLorentzVector* gamma, TLorentzVector* uncorLepton1, TLorentzVector* uncorLepton2, TLorentzVector* uncorGamma, float R9input, float SCEta, ofstream & dump, float eta1, float eta2){
  int cat = -1;
  bool goodLep = false;
  if (selection == "mumuGamma"){
    goodLep = (fabs(lepton1->Eta()) < 0.9 || fabs(lepton2->Eta()) < 0.9) && (fabs(lepton1->Eta()) < 2.1 && fabs(lepton2->Eta()) < 2.1);
  }else if (selection == "eeGamma"){
    goodLep = (fabs(eta1) < 1.4442 && fabs(eta2) < 1.4442);
  }

  if (goodLep && (fabs(SCEta) < 1.4442) ){
    if (R9input>=R9Cut) cat = 0;
    else cat = 3;
  }else if (!goodLep && (fabs(SCEta) < 1.4442) ){
    cat = 1;
  } else {
    cat = 2;
  }
  float mll = (*lepton1 + *lepton2).M();
  float mllg = (*lepton1 + *lepton2 + *gamma).M();
  float mllUncor = (*uncorLepton1 + *uncorLepton2).M();
  float mllgUncor = (*uncorLepton1 + *uncorLepton2 + *uncorGamma).M();

  dump << "cat: "      << setw(2)  << cat           << " run: "      << setw(7)  << runNumber           << " evt: "     << setw(10) << eventNumber
       << " lumi: "    << setw(5)  << lumiSection   << " mll: "      << setw(10) << mllUncor            << " phopt: "   << setw(10) << uncorGamma->Pt()
       << " mllg: "    << setw(10) << mllgUncor     << " l1eta: "    << setw(10) << uncorLepton1->Eta() << " l2eta: "   << setw(10) << uncorLepton2->Eta()
       << " pSCeta: "  << setw(10) << SCEta

       << " mllCor: "  << setw(10) << mll           << " phoptCor: " << setw(10) << gamma->Pt()         << " l1PtCor: " << setw(10) << lepton1->Pt()
       << " l2PtCor: " << setw(10) << lepton2->Pt() << " mllgCor: "  << setw(10) << mllg
       << " r9: "      << setw(10) << R9input       << endl;
}

void higgsAnalyzer::FinalDumper(TLorentzVector* lepton1, TLorentzVector* lepton2, TLorentzVector* gamma, int catNum, ofstream & dump){
  float mll = (*lepton1 + *lepton2).M();
  float mllg = (*lepton1 + *lepton2 + *gamma).M();

  if (doEleReg){
    dump << "cat: "   << setw(2)  << catNum        << " evt: "   << setw(10) << eventNumber << " l1RegE-P: " << setw(10) << lepton1->E() << " l1pt: " << setw(10) << lepton1->Pt()
         << " l2RegE-P: " << setw(10) << lepton2->E() << " l2pt: " << setw(10) << lepton2->Pt() << " l1eta: " << setw(10) << lepton1->Eta() << " l2eta: " << setw(10) << lepton2->Eta() << " phopt: " << setw(10) << gamma->Pt() << " mll: "  << setw(10) << mll
         << " mllg: " << setw(10) << mllg          << endl;
  }else{
    dump << "cat: "   << setw(2)  << catNum        << " evt: "   << setw(10) << eventNumber << " l1pt: " << setw(10) << lepton1->Pt()
         << " l2pt: " << setw(10) << lepton2->Pt() << " l1eta: " << setw(10) << lepton1->Eta() << " l2eta: " << setw(10) << lepton2->Eta() << " phopt: " << setw(10) << gamma->Pt() << " mll: "  << setw(10) << mll
         << " mllg: " << setw(10) << mllg          << endl;
  }
}

void higgsAnalyzer::MVADumper(TCElectron* ele, EGammaMvaEleEstimator* mvaMaker, double rhoFactor, ofstream & dump){
  bool passPreSel = false;
  bool passMVA = false;
  if (ele->IdMap("preSelPassV1")) passPreSel = true;			
  double tmpMVAValue = myMVATrig->mvaValue(
      ele->IdMap("fbrem"),    
      ele->IdMap("kfChi2"),                            
      ele->IdMap("kfNLayers"),                            
      ele->IdMap("gsfChi2"),                           
      ele->IdMap("dEta"),
      ele->IdMap("dPhi"),            
      ele->IdMap("dEtaAtCalo"),
      ele->SigmaIEtaIEta(), 
      ele->IdMap("SigmaIPhiIPhi"),
      ele->IdMap("SCEtaWidth"),
      ele->IdMap("SCPhiWidth"),
      ele->IdMap("ome1x5oe5x5"),
      ele->IdMap("R9"),
      ele->HadOverEm(),
      ele->IdMap("EoP"),
      ele->IdMap("ooemoopV1"), 
      ele->IdMap("eopOut"), 
      ele->IdMap("preShowerORaw"), 
      ele->IdMap("d0"),
      ele->IdMap("ip3d"), 
      ele->SCEta(),
      ele->Pt(),
      false);                

  if (fabs(ele->SCEta()) < 0.8 && tmpMVAValue > 0.5){
    passMVA = true;
  }else if (fabs(ele->SCEta()) > 0.8 && fabs(ele->SCEta()) < 1.479 && tmpMVAValue > 0.12){
    passMVA = true;
  }else if (fabs(ele->SCEta()) > 1.479 && tmpMVAValue > 0.6){
    passMVA = true;
  }

  bool passIso = false;
  passIso = PassElectronIso(ele, looseElIso);

dump << " run: "                   << setw(7)  << runNumber                            << " evt: "                   << setw(10) << eventNumber                          << " lumi: "                  << setw(5)  << lumiSection
     << " pt: "                    << setw(10) << ele->Pt()                            << " SCeta: "                 << setw(10) << ele->SCEta()                         << " fbrem: "                 << setw(10) << ele->IdMap("fbrem")
     << " kfChi2: "                << setw(10) << ele->IdMap("kfChi2")                 << " kfHits: "                << setw(10) << ele->IdMap("kfNLayers")              << " gsfChi2: "               << setw(10) << ele->IdMap("gsfChi2")
     << " dEta: "                  << setw(10) << ele->IdMap("dEta")                   << " dPhi: "                  << setw(10) << ele->IdMap("dPhi")                   << " dEtaCalo: "              << setw(10) << ele->IdMap("dEtaAtCalo")
     << " sieie: "                 << setw(10) << ele->SigmaIEtaIEta()                 << " sipip: "                 << setw(10) << ele->IdMap("SigmaIPhiIPhi")          << " etaWidth: "              << setw(10) << ele->IdMap("SCEtaWidth")
     << " phiWidth: "              << setw(10) << ele->IdMap("SCPhiWidth")             << " ome1x5oe5x5: "           << setw(10) << ele->IdMap("ome1x5oe5x5")            << " R9: "                    << setw(10) << ele->IdMap("R9")
     << " HoE: "                   << setw(10) << ele->HadOverEm()                     << " EoP: "                   << setw(10) << ele->IdMap("EoP")                    << " ooemoop: "               << setw(10) << ele->IdMap("ooemoopV1")
     << " EoPout: "                << setw(10) << ele->IdMap("eopOut")                 << " preShowerOverRaw: "      << setw(10) << ele->IdMap("preShowerORaw")          << " d0: "                    << setw(10) << ele->IdMap("d0")
     << " ip3d: "                  << setw(10) << ele->IdMap("ip3d")                   << " ChargedIso_DR0p0To0p1: " << setw(10) << ele->IsoMap("ChargedIso_DR0p0To0p1") << " ChargedIso_DR0p1To0p2: " << setw(10) << ele->IsoMap("ChargedIso_DR0p1To0p2")
     << " ChargedIso_DR0p2To0p3: " << setw(10) << ele->IsoMap("ChargedIso_DR0p2To0p3") << " ChargedIso_DR0p3To0p4: " << setw(10) << ele->IsoMap("ChargedIso_DR0p3To0p4") << " ChargedIso_DR0p4To0p5: " << setw(10) << ele->IsoMap("ChargedIso_DR0p4To0p5")
     << " GammaIso_DR0p0To0p1: "   << setw(10) << ele->IsoMap("GammaIso_DR0p0To0p1")   << " GammaIso_DR0p1To0p2: "   << setw(10) << ele->IsoMap("GammaIso_DR0p1To0p2")   << " GammaIso_DR0p2To0p3: "   << setw(10) << ele->IsoMap("GammaIso_DR0p2To0p3")
     << " GammaIso_DR0p3To0p4: "   << setw(10) << ele->IsoMap("GammaIso_DR0p3To0p4")   << " GammaIso_DR0p4To0p5: "   << setw(10) << ele->IsoMap("GammaIso_DR0p4To0p5")
     << " NeutralHadronIso_DR0p0To0p1: " << setw(10) << ele->IsoMap("NeutralHadronIso_DR0p0To0p1") << " NeutralHadronIso_DR0p1To0p2: " << setw(10) << ele->IsoMap("NeutralHadronIso_DR0p1To0p2")
     << " NeutralHadronIso_DR0p2To0p3: " << setw(10) << ele->IsoMap("NeutralHadronIso_DR0p2To0p3") << " NeutralHadronIso_DR0p3To0p4: " << setw(10) << ele->IsoMap("NeutralHadronIso_DR0p3To0p4")
     << " NeutralHadronIso_DR0p4To0p5: " << setw(10) << ele->IsoMap("NeutralHadronIso_DR0p4To0p5") << " rhoFactor: "                   << setw(10) << rhoFactor
     << " mva Value: " << setw(10) << tmpMVAValue << " pass PreSel: " << setw(5) << passPreSel << " pass MVA: " << setw(5) << passMVA << " pass ISO: " << setw(5) << passIso
     << endl;
}

bool higgsAnalyzer::FindGoodZElectron(vector<TCElectron> electronList, TLorentzVector* lepton1, TLorentzVector* lepton2, TLorentzVector* ZP4,float* eta1, float* eta2, int* int1, int* int2){
  TLorentzVector tmpZ;
  bool goodZ = false;
  float ZmassDiff=99999;
  for(unsigned int i =0; i<electronList.size(); i++){
    if (electronList[i].Pt() > elePtCut[0]){
      for(unsigned int j =1; j<electronList.size(); j++){
        if (electronList[j].Pt() > elePtCut[1] && electronList[j].Charge() != electronList[i].Charge()){
          goodZ = true;
          tmpZ = (electronList[i]+electronList[j]);
          if(fabs(91.1876-tmpZ.M()) < ZmassDiff){
            *ZP4 = (electronList[i]+electronList[j]);
            *lepton1 = electronList[i];
            *int1 = i;
            *lepton2 = electronList[j];
            *int2 = j;
            *eta1 = electronList[i].SCEta();
            *eta2 = electronList[j].SCEta();
            ZmassDiff = fabs(91.1876-tmpZ.M());
          }
        }
      }
    }
  }
  return goodZ;
}
          
bool higgsAnalyzer::FindGoodZElectron(vector<TCElectron> electronList, vector<TCElectron> uncorElectronList, TLorentzVector* lepton1, TLorentzVector* lepton2, TLorentzVector* uncorLepton1, TLorentzVector* uncorLepton2, TLorentzVector* ZP4, float* eta1, float* eta2, int* int1, int* int2){
  TLorentzVector tmpZ;
  bool goodZ = false;
  float ZmassDiff=99999;
  for(unsigned int i =0; i<electronList.size(); i++){
    if (electronList[i].Pt() > elePtCut[0]){
      for(unsigned int j =1; j<electronList.size(); j++){
        if (electronList[j].Pt() > elePtCut[1] && electronList[j].Charge() != electronList[i].Charge()){
          goodZ = true;
          /// do regression
          if (doEleReg && period.find("2012") != string::npos){
            if (electronList[i].RegressionMomCombP4().E() != 0){
              electronList[i].SetPtEtaPhiE(electronList[i].RegressionMomCombP4().E()/cosh(electronList[i].RegressionMomCombP4().Eta()),electronList[i].RegressionMomCombP4().Eta(),electronList[i].RegressionMomCombP4().Phi(),electronList[i].RegressionMomCombP4().E());
              //electronList[i].SetPtEtaPhiE(electronList[i].RegressionMomCombP4().Pt(),electronList[i].RegressionMomCombP4().Eta(),electronList[i].RegressionMomCombP4().Phi(),electronList[i].RegressionMomCombP4().E());
            }
            if (electronList[j].RegressionMomCombP4().E() != 0){
              electronList[j].SetPtEtaPhiE(electronList[j].RegressionMomCombP4().E()/cosh(electronList[j].RegressionMomCombP4().Eta()),electronList[j].RegressionMomCombP4().Eta(),electronList[j].RegressionMomCombP4().Phi(),electronList[j].RegressionMomCombP4().E());
              //electronList[j].SetPtEtaPhiE(electronList[j].RegressionMomCombP4().Pt(),electronList[j].RegressionMomCombP4().Eta(),electronList[j].RegressionMomCombP4().Phi(),electronList[j].RegressionMomCombP4().E());
            }
          }
          /// do eng cor
          if(engCor){
            float energyElCorI;
            float energyElCorJ;
            if ( period.find("2011") != string::npos ){
              energyElCorI = correctedElectronEnergy( electronList[i].E(), electronList[i].SCEta(), electronList[i].R9(), runNumber, 0, "2011", !isRealData, rEl );
              energyElCorJ = correctedElectronEnergy( electronList[j].E(), electronList[j].SCEta(), electronList[j].R9(), runNumber, 0, "2011", !isRealData, rEl );
            }else{
              if(doEleReg && electronList[i].RegressionMomCombP4().E() != 0){
                energyElCorI = correctedElectronEnergy( electronList[i].E(), electronList[i].SCEta(), electronList[i].R9(), runNumber, 1, "Moriond2013", !isRealData, rEl );
              }else{
                energyElCorI = correctedElectronEnergy( electronList[i].E(), electronList[i].SCEta(), electronList[i].R9(), runNumber, 0, "HCP2012", !isRealData, rEl );
              }
              if(doEleReg && electronList[j].RegressionMomCombP4().E() != 0){
                energyElCorJ = correctedElectronEnergy( electronList[j].E(), electronList[j].SCEta(), electronList[j].R9(), runNumber, 1, "Moriond2013", !isRealData, rEl );
              }else{
                energyElCorJ = correctedElectronEnergy( electronList[j].E(), electronList[j].SCEta(), electronList[j].R9(), runNumber, 0, "HCP2012", !isRealData, rEl );
              }
            }
            //float newPtI = sqrt((pow(energyElCorI,2)-pow(0.000511,2))/pow(cosh(electronList[i].Eta()),2));
            //float newPtJ = sqrt((pow(energyElCorJ,2)-pow(0.000511,2))/pow(cosh(electronList[j].Eta()),2));
            float newPtI = energyElCorI/cosh(electronList[i].Eta());
            float newPtJ = energyElCorJ/cosh(electronList[j].Eta());
            electronList[i].SetPtEtaPhiM(newPtI,electronList[i].Eta(),electronList[i].Phi(),0.000511);
            electronList[j].SetPtEtaPhiM(newPtJ,electronList[j].Eta(),electronList[j].Phi(),0.000511);
          }
          tmpZ = (electronList[i]+electronList[j]);
          if(fabs(91.1876-tmpZ.M()) < ZmassDiff){
            *ZP4 = (electronList[i]+electronList[j]);
            *lepton1 = electronList[i];
            *lepton2 = electronList[j];
            *int1 = i;
            *int2 = j;
            *eta1 = electronList[i].SCEta();
            *eta2 = electronList[j].SCEta();
            *uncorLepton1 = uncorElectronList[i];
            *uncorLepton2 = uncorElectronList[j];
            ZmassDiff = fabs(91.1876-tmpZ.M());
          }
        }
      }
    }
  }
  return goodZ;
}
          
bool higgsAnalyzer::FindGoodZMuon(vector<TCMuon> muonList, TLorentzVector* lepton1, TLorentzVector* lepton2, TLorentzVector* ZP4, int* int1, int* int2){
  TLorentzVector tmpZ;
  bool goodZ = false;
  float ZmassDiff=99999;
  for(unsigned int i =0; i<muonList.size(); i++){
    if (muonList[i].Pt() > muPtCut[0]){
      for(unsigned int j =1; j<muonList.size(); j++){
        if (muonList[j].Pt() > muPtCut[1] && muonList[j].Charge() != muonList[i].Charge()){
          goodZ = true;
          tmpZ = (muonList[i]+muonList[j]);
          if(fabs(91.1876-tmpZ.M()) < ZmassDiff){
            *ZP4 = (muonList[i]+muonList[j]);
            *lepton1 = muonList[i];
            *lepton2 = muonList[j];
            *int1 = i;
            *int2 = j;
            ZmassDiff = fabs(91.1876-tmpZ.M());
          }
        }
      }
    }
  }
  return goodZ;
}

bool higgsAnalyzer::FindGoodZMuon(vector<TCMuon> muonList, vector<TCMuon> uncorMuonList, TLorentzVector* lepton1, TLorentzVector* lepton2, TLorentzVector* uncorLepton1, TLorentzVector* uncorLepton2, TLorentzVector* ZP4, int* int1, int* int2){
  TLorentzVector tmpZ;
  bool goodZ = false;
  float ZmassDiff=99999;
  for(unsigned int i =0; i<muonList.size(); i++){
    if (muonList[i].Pt() > muPtCut[0]){
      for(unsigned int j =1; j<muonList.size(); j++){
        if (muonList[j].Pt() > muPtCut[1] && muonList[j].Charge() != muonList[i].Charge()){
          goodZ = true;
          tmpZ = (muonList[i]+muonList[j]);
          if(fabs(91.1876-tmpZ.M()) < ZmassDiff){
            *ZP4 = (muonList[i]+muonList[j]);
            *lepton1 = muonList[i];
            *lepton2 = muonList[j];
            *uncorLepton1 = uncorMuonList[i];
            *uncorLepton2 = uncorMuonList[j];
            *int1 = i;
            *int2 = j;
            ZmassDiff = fabs(91.1876-tmpZ.M());
          }
        }
      }
    }
  }
  return goodZ;
}


void higgsAnalyzer::LumiXSWeight(float * LumiXSWeight){
  if (!isRealData){
    if(period.find("2011") != string::npos){

      if(suffix.find("gg") != string::npos){
        if(suffix.find("120")!=string::npos) *LumiXSWeight = 99990/(16.65*0.00111*0.10098*1000);
        if(suffix.find("125")!=string::npos) *LumiXSWeight = 95193/(15.32*0.00154*0.10098*1000);
        if(suffix.find("130")!=string::npos) *LumiXSWeight = 99992/(14.16*0.00195*0.10098*1000);
        if(suffix.find("135")!=string::npos) *LumiXSWeight = 99995/(13.11*0.00227*0.10098*1000);
        if(suffix.find("140")!=string::npos) *LumiXSWeight = 99988/(12.18*0.00246*0.10098*1000);
        if(suffix.find("145")!=string::npos) *LumiXSWeight = 99997/(11.33*0.00248*0.10098*1000);
        if(suffix.find("150")!=string::npos) *LumiXSWeight = 99998/(10.58*0.00231*0.10098*1000);

      }else if(suffix.find("vbf") != string::npos){
        if(suffix.find("120")!=string::npos) *LumiXSWeight = 99885/(1.269*0.00111*0.10098*1000);
        if(suffix.find("125")!=string::npos) *LumiXSWeight = 99890/(1.211*0.00154*0.10098*1000);
        if(suffix.find("130")!=string::npos) *LumiXSWeight = 99899/(1.154*0.00195*0.10098*1000);
        if(suffix.find("135")!=string::npos) *LumiXSWeight = 99913/(1.100*0.00227*0.10098*1000);
        if(suffix.find("140")!=string::npos) *LumiXSWeight = 99882/(1.052*0.00246*0.10098*1000);
        if(suffix.find("145")!=string::npos) *LumiXSWeight = 99893/(1.023*0.00248*0.10098*1000);
        if(suffix.find("150")!=string::npos) *LumiXSWeight = 99892/(0.9800*0.00231*0.10098*1000);

      }else if(suffix.find("tth") != string::npos){
        if(suffix.find("120")!=string::npos) *LumiXSWeight = 99999/(0.0976*0.00111*0.10098*1000);
        if(suffix.find("125")!=string::npos) *LumiXSWeight = 99999/(0.0863*0.00154*0.10098*1000);
        if(suffix.find("130")!=string::npos) *LumiXSWeight = 100076/(0.0766*0.00195*0.10098*1000);
        if(suffix.find("135")!=string::npos) *LumiXSWeight = 99999/(0.0681*0.00227*0.10098*1000);
        if(suffix.find("140")!=string::npos) *LumiXSWeight = 99999/(0.0607*0.00246*0.10098*1000);
        if(suffix.find("145")!=string::npos) *LumiXSWeight = 100284/(0.0544*0.00248*0.10098*1000);
        if(suffix.find("150")!=string::npos) *LumiXSWeight = 100000/(0.0487*0.00231*0.10098*1000);

      }else if(suffix.find("wh") != string::npos){
        if(suffix.find("120")!=string::npos) *LumiXSWeight = 65152/(0.6561*0.00111*0.10098*1000);
        if(suffix.find("125")!=string::npos) *LumiXSWeight = 65507/(0.5729*0.00154*0.10098*1000);
        if(suffix.find("130")!=string::npos) *LumiXSWeight = 65659/(0.5008*0.00195*0.10098*1000);
        if(suffix.find("135")!=string::npos) *LumiXSWeight = 65461/(0.4390*0.00227*0.10098*1000);
        if(suffix.find("140")!=string::npos) *LumiXSWeight = 65491/(0.3857*0.00246*0.10098*1000);
        if(suffix.find("145")!=string::npos) *LumiXSWeight = 65491/(0.3406*0.00248*0.10098*1000);
        if(suffix.find("150")!=string::npos) *LumiXSWeight = 65205/(0.3001*0.00231*0.10098*1000);

      }else if(suffix.find("zh") != string::npos){
        if(suffix.find("120")!=string::npos) *LumiXSWeight = 34294/(0.3598*0.00111*0.10098*1000);
        if(suffix.find("125")!=string::npos) *LumiXSWeight = 34492/(0.3158*0.00154*0.10098*1000);
        if(suffix.find("130")!=string::npos) *LumiXSWeight = 34340/(0.2778*0.00195*0.10098*1000);
        if(suffix.find("135")!=string::npos) *LumiXSWeight = 34538/(0.2453*0.00227*0.10098*1000);
        if(suffix.find("140")!=string::npos) *LumiXSWeight = 34049/(0.2172*0.00246*0.10098*1000);
        if(suffix.find("145")!=string::npos) *LumiXSWeight = 34653/(0.1930*0.00248*0.10098*1000);
        if(suffix.find("150")!=string::npos) *LumiXSWeight = 34794/(0.1713*0.00231*0.10098*1000);
      }
      if (selection=="mumuGamma") *LumiXSWeight = 5.05/(*LumiXSWeight);
      if (selection=="eeGamma") *LumiXSWeight = 4.98/(*LumiXSWeight);

    }else if(period.find("2012") != string::npos){

      if(suffix.find("gg") != string::npos){
        if(suffix.find("120")!=string::npos) *LumiXSWeight = 99992/(21.13*0.00111*0.10098*1000);
        if(suffix.find("125")!=string::npos) *LumiXSWeight = 99991/(19.52*0.00154*0.100974*1000);
        if(suffix.find("130")!=string::npos) *LumiXSWeight = 99991/(18.07*0.00195*0.10098*1000);
        if(suffix.find("135")!=string::npos) *LumiXSWeight = 99996/(16.79*0.00227*0.10098*1000);
        if(suffix.find("140")!=string::npos) *LumiXSWeight = 96994/(15.63*0.00246*0.10098*1000);
        if(suffix.find("145")!=string::npos) *LumiXSWeight = 99990/(14.59*0.00248*0.10098*1000);
        if(suffix.find("150")!=string::npos) *LumiXSWeight = 99991/(13.65*0.00231*0.10098*1000);

      }else if(suffix.find("vbf") != string::npos){
        if(suffix.find("120")!=string::npos) *LumiXSWeight = 99886/(1.649*0.00111*0.10098*1000);
        if(suffix.find("125")!=string::npos) *LumiXSWeight = 99885/(1.578*0.00154*0.10098*1000);
        if(suffix.find("130")!=string::npos) *LumiXSWeight = 99899/(1.511*0.00195*0.10098*1000);
        if(suffix.find("135")!=string::npos) *LumiXSWeight = 99004/(1.448*0.00227*0.10098*1000);
        if(suffix.find("140")!=string::npos) *LumiXSWeight = 99890/(1.389*0.00246*0.10098*1000);
        if(suffix.find("145")!=string::npos) *LumiXSWeight = 99888/(1.333*0.00248*0.10098*1000);
        if(suffix.find("150")!=string::npos) *LumiXSWeight = 99893/(1.280*0.00231*0.10098*1000);

      }else if(suffix.find("tth") != string::npos){
        if(suffix.find("120")!=string::npos) *LumiXSWeight = 100080/(0.1470*0.00111*0.10098*1000);
        if(suffix.find("125")!=string::npos) *LumiXSWeight = 100048/(0.1302*0.00154*0.10098*1000);
        if(suffix.find("130")!=string::npos) *LumiXSWeight = 100368/(0.1157*0.00195*0.10098*1000);
        if(suffix.find("135")!=string::npos) *LumiXSWeight = 100130/(0.1031*0.00227*0.10098*1000);
        if(suffix.find("140")!=string::npos) *LumiXSWeight = 100340/(0.09207*0.00246*0.10098*1000);
        if(suffix.find("145")!=string::npos) *LumiXSWeight = 97869/(0.08246*0.00248*0.10098*1000);
        if(suffix.find("150")!=string::npos) *LumiXSWeight = 100048/(0.07403*0.00231*0.10098*1000);

      }else if(suffix.find("wh") != string::npos){
        if(suffix.find("120")!=string::npos) *LumiXSWeight = 656153/(0.7966*0.00111*0.10098*1000);
        if(suffix.find("125")!=string::npos) *LumiXSWeight = 656101/(0.6966*0.00154*0.10098*1000);
        if(suffix.find("130")!=string::npos) *LumiXSWeight = 614920/(0.6095*0.00195*0.10098*1000);
        if(suffix.find("135")!=string::npos) *LumiXSWeight = 328690/(0.5351*0.00227*0.10098*1000);
        if(suffix.find("140")!=string::npos) *LumiXSWeight = 653556/(0.4713*0.00246*0.10098*1000);
        if(suffix.find("145")!=string::npos) *LumiXSWeight = 652944/(0.4164*0.00248*0.10098*1000);
        if(suffix.find("150")!=string::npos) *LumiXSWeight = 652973/(0.3681*0.00231*0.10098*1000);

      }else if(suffix.find("zh") != string::npos){
        if(suffix.find("120")!=string::npos) *LumiXSWeight = 344220/(0.4483*0.00111*0.10098*1000);
        if(suffix.find("125")!=string::npos) *LumiXSWeight = 344143/(0.3943*0.00154*0.10098*1000);
        if(suffix.find("130")!=string::npos) *LumiXSWeight = 323352/(0.3473*0.00195*0.10098*1000);
        if(suffix.find("135")!=string::npos) *LumiXSWeight = 173822/(0.3074*0.00227*0.10098*1000);
        if(suffix.find("140")!=string::npos) *LumiXSWeight = 346650/(0.2728*0.00246*0.10098*1000);
        if(suffix.find("145")!=string::npos) *LumiXSWeight = 347056/(0.2424*0.00248*0.10098*1000);
        if(suffix.find("150")!=string::npos) *LumiXSWeight = 347459/(0.2159*0.00231*0.10098*1000);
      }
      if (selection=="mumuGamma"){
        if (abcd == "AB") *LumiXSWeight = 5.3185/(*LumiXSWeight);
        else if (abcd == "CD") *LumiXSWeight = 14.299/(*LumiXSWeight);
        else if (abcd == "ABCD") *LumiXSWeight = 19.6175/(*LumiXSWeight);
        else cout<<"abcd not specified"<<endl;
      }else if (selection=="eeGamma"){
        if (abcd == "AB") *LumiXSWeight = 5.3185/(*LumiXSWeight);
        else if (abcd == "CD") *LumiXSWeight = 14.301/(*LumiXSWeight);
        else if (abcd == "ABCD") *LumiXSWeight = 19.6195/(*LumiXSWeight);
        else cout<<"abcd not specified"<<endl;
      }
    }
  }
}


