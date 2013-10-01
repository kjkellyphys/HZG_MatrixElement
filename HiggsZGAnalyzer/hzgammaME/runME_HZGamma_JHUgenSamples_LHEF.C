// 
// This code calculates the ME
// 

#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TObject.h"
#include "TProfile.h"
#include <iostream>
#include "Math/LorentzVector.h"
#include "TLorentzRotation.h"
#include "Math/VectorUtil.h"
#include "TParticle.h"
#include "ExRootAnalysis/src/ExRootClasses.cc"
#include "TClonesArray.h"

// ME related
#include "TVar.hh"
#include "TEvtProb.hh"
#include "math.h"

// ZG Angles from Brian/Kristian
#include "../plugins/ZGAngles.cc"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector; 

float ERRORthreshold=1.0;
using namespace std;

void xseccalc(TString inputDir, TString fileName, TString outputDir, int maxevt, TVar::VerbosityLevel verbosity);


//###################
//# main function
//###################
void runME_hzgamma_JHUgenSamples(TString inputDir, TString fileName, TString outputDir, int maxevt, TVar::VerbosityLevel verbosity=TVar::INFO){

  if (verbosity >= TVar::INFO) cout <<"=== Calculating differential cross-section ==========" <<endl;  
  xseccalc(inputDir, fileName, outputDir, maxevt, verbosity); 
 
}

void xseccalc(TString inputDir, TString fileName, TString outputDir, int maxevt, TVar::VerbosityLevel verbosity){

  if (verbosity >= TVar::INFO) cout << "Input File: " << fileName << " \n";
  
  bool verbose;
  if (verbosity > TVar::INFO)
    verbose = true;
  else
    verbose = false;

  TFile* fin = new TFile(inputDir+fileName);

  ofstream KKTxt;
  bool Signal = true;
  KKTxt.open ("Disc_125.txt");

  ofstream CutTxt;
  CutTxt.open("CutFile.txt");

  //Output file 1: Matrix Element Calculations (Stored in Branches)
  TString outFileName = outputDir+fileName;
  outFileName.ReplaceAll(".root","_ME.root");
  cout << outFileName <<endl;
  TFile *newfile = new TFile(outFileName, "recreate");
  
  //Output file 2: Histograms (weighted from MCFM)
  TString histFileName = outputDir+fileName;
  histFileName.ReplaceAll(".root","_AnglePlots.root");
  cout << histFileName <<endl;
  TFile *histfile = new TFile(histFileName, "RECREATE");

  newfile->cd();
  TTree* ch=(TTree*)fin->Get("LHEF"); 
  if (ch==0x0) return; 
  TTree* evt_tree=(TTree*) ch->CloneTree(0, "fast");
  evt_tree->SetName("newTree");

  //------------------------------
  //       Input Variables
  //------------------------------

  //Weight from MCFM-Produced Files
  float weight = 0.;
  weight = 1.0; //For uniform weight.
  
  TClonesArray *Event = new TClonesArray("TRootLHEFEvent");
  TClonesArray *Particle = new TClonesArray("TRootLHEFParticle"); 
  std::vector<TRootLHEFParticle*> PartCollec(10);

  ch->SetBranchAddress("Event", &Event);
  ch->SetBranchAddress("Particle", &Particle);

  //------------------------
  //      Output Files
  //------------------------

  // output variables (in _ME.root file)
  float dXsec_ZGam_MCFM = 0.;
  float dXsec_HZGam_MCFM = 0.;
  float Discriminant = 0.;
  float logBkg(0.), logSig(0.);
  
  evt_tree->Branch("dxSec_ZGam_MCFM"   , &dXsec_ZGam_MCFM ,"dXsec_ZGam_MCFM/F");
  evt_tree->Branch("dxSec_HZGam_MCFM"   , &dXsec_HZGam_MCFM ,"dXsec_HZGam_MCFM/F");
  evt_tree->Branch("logBkg", &logBkg, "logBkg/F");
  evt_tree->Branch("logSig", &logSig, "logSig/F");
  evt_tree->Branch("Discriminant", &Discriminant, "Discriminant/F");
  
  //User-defined Histograms (_AnglePlots.root file)
  histfile->cd();

  TH2F* DvMlla = new TH2F("DvMlla", "Discriminant v. 3-Body Mass;D;M_{ll#gamma}", 100, 0, 0.5, 100, 90, 190);
  TH2F* PsvMlla = new TH2F("PsvMlla", "Discriminant v. 3-Body Mass;Ps;M_{ll#gamma}", 100, 0, 0.0001, 100, 90, 190);
  TH1F* Mlla_PostDCut = new TH1F("Mlla_PostDCut", "Mass of 3-Body", 100, 90, 190);

  TH1F* me_pg = new TH1F("me_pg", "Mass of Positron/Photon", 100, 0, 150);
  TH1F* me_ng = new TH1F("me_ng", "Mass of Electron/Photon", 100, 0, 150);
  TH1F* mm_pg = new TH1F("mm_pg", "Mass of AntiMuon/Photon", 100, 0, 150);
  TH1F* mm_ng = new TH1F("mm_ng", "Mass of Muon/Photon", 100, 0, 150);
  
  TH2F* megVmee = new TH2F("megVmee", "Mass of Electron/Photon pair vs. Mass of Di-Electron Pair", 100, 0, 150, 100, 0, 150);
  TH2F* mmgVmmm = new TH2F("mmgVmmm", "Mass of Muon/Photon pair vs. Mass of Di-Muon Pair", 100, 0, 150, 100, 0, 150);

  TH2F* meegVmeg = new TH2F("meegVmeg", "Mass of 3-Body vs. Mass of Electron/Photon;M_{ee#gamma};M_{e#gamma}", 50, 90, 190, 100, 0, 150);
  TH2F* mmmgVmmg = new TH2F("mmmgVmmg", "Mass of 3-Body vs. Mass of Muon/Photon;M_{#mu#mu#gamma};M_{#mu#gamma}", 50, 90, 190, 100, 0, 150);

  TH1F* CosT_lp = new TH1F("CosT_lp", "Cos(#theta) positive lepton;cos(#theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* CosT_lm = new TH1F("CosT_lm", "Cos(#theta) negative lepton;cos(#theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* Phi_lp = new TH1F("Phi_lp", "#phi positive lepton;#phi;N_{evts}", 50, -3.2, 3.2);
  TH1F* CosTZG = new TH1F("CosTZG", "Cos(#Theta) ZG system;cos(#Theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* Pt_g = new TH1F("Pt_g", "Pt of photon in ZG System;Pt;N_{evts}", 50, 0., 200.);
  TH1F* Pt_lp = new TH1F("Pt_lp", "Pt of positive lepton;Pt;N_{evts}", 50, 0., 200.);
  TH1F* Pt_lm = new TH1F("Pt_lm", "Pt of negative lepton;Pt;N_{evts}", 50, 0., 200.);
  TH1F* eta_g = new TH1F("eta_g", "Eta of Photon", 50, -5.0, 5.0);
  TH1F* eta_lp = new TH1F("eta_lp", "Eta of positive lepton;#eta;N_{evts}", 50, -5.0, 5.0);
  TH1F* eta_lm = new TH1F("eta_lm", "Eta of negative lepton;#eta;N_{evts}", 50, -5.0, 5.0);
  TH1F* mzg = new TH1F("mzg", "Mass of Z-Gamma System;Mass;N_{evts}", 100, 90, 190);
  TH1F* mz = new TH1F("mz", "Mass of Z Boson;Mass;N_{evts}", 100, 0, 150);

  TH1F* CosT_lpA = new TH1F("CosT_lpA", "Cos(#theta) positive lepton;cos(#theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* CosT_lmA = new TH1F("CosT_lmA", "Cos(#theta) negative lepton;cos(#theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* Phi_lpA = new TH1F("Phi_lpA", "#phi positive lepton;#phi;N_{evts}", 50, -3.2, 3.2);
  TH1F* CosTZGA = new TH1F("CosTZGA", "Cos(#Theta) ZG system;cos(#Theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* Pt_gA = new TH1F("Pt_gA", "Pt of photon in ZG System;Pt;N_{evts}", 50, 0., 200.);
  TH1F* Pt_lpA = new TH1F("Pt_lpA", "Pt of positive lepton;Pt;N_{evts}", 50, 0., 200.);
  TH1F* Pt_lmA = new TH1F("Pt_lmA", "Pt of negative lepton;Pt;N_{evts}", 50, 0., 200.);
  TH1F* eta_gA = new TH1F("eta_gA", "Eta of Photon", 50, -5.0, 5.0);
  TH1F* eta_lpA = new TH1F("eta_lpA", "Eta of positive lepton;#eta;N_{evts}", 50, -5.0, 5.0);
  TH1F* eta_lmA = new TH1F("eta_lmA", "Eta of negative lepton;#eta;N_{evts}", 50, -5.0, 5.0);
  TH1F* mllgA = new TH1F("mllgA", "Mass of ll-Gamma System;Mass;N_{evts}", 100, 90, 150);
  TH1F* mllA = new TH1F("mllA", "Mass of Lepton Pair;Mass;N_{evts}", 100, 0, 150);

  TH1F* CosT_lpB = new TH1F("CosT_lpB", "Cos(#theta) positive lepton;cos(#theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* CosT_lmB = new TH1F("CosT_lmB", "Cos(#theta) negative lepton;cos(#theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* Phi_lpB = new TH1F("Phi_lpB", "#phi positive lepton;#phi;N_{evts}", 50, -3.2, 3.2);
  TH1F* CosTZGB = new TH1F("CosTZGB", "Cos(#Theta) ZG system;cos(#Theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* Pt_gB = new TH1F("Pt_gB", "Pt of photon in ZG System;Pt;N_{evts}", 50, 0., 200.);
  TH1F* Pt_lpB = new TH1F("Pt_lpB", "Pt of positive lepton;Pt;N_{evts}", 50, 0., 200.);
  TH1F* Pt_lmB = new TH1F("Pt_lmB", "Pt of negative lepton;Pt;N_{evts}", 50, 0., 200.);
  TH1F* eta_gB = new TH1F("eta_gB", "Eta of Photon", 50, -5.0, 5.0);
  TH1F* eta_lpB = new TH1F("eta_lpB", "Eta of positive lepton;#eta;N_{evts}", 50, -5.0, 5.0);
  TH1F* eta_lmB = new TH1F("eta_lmB", "Eta of negative lepton;#eta;N_{evts}", 50, -5.0, 5.0);
  TH1F* mllgB = new TH1F("mllgB", "Mass of ll-Gamma System;Mass;N_{evts}", 100, 90, 150);
  TH1F* mllB = new TH1F("mllB", "Mass of Lepton Pair;Mass;N_{evts}", 100, 0, 150);

  TH1F* CosT_lpC = new TH1F("CosT_lpC", "Cos(#theta) positive lepton;cos(#theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* CosT_lmC = new TH1F("CosT_lmC", "Cos(#theta) negative lepton;cos(#theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* Phi_lpC = new TH1F("Phi_lpC", "#phi positive lepton;#phi;N_{evts}", 50, -3.2, 3.2);
  TH1F* CosTZGC = new TH1F("CosTZGC", "Cos(#Theta) ZG system;cos(#Theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* Pt_gC = new TH1F("Pt_gC", "Pt of photon in ZG System;Pt;N_{evts}", 50, 0., 200.);
  TH1F* Pt_lpC = new TH1F("Pt_lpC", "Pt of positive lepton;Pt;N_{evts}", 50, 0., 200.);
  TH1F* Pt_lmC = new TH1F("Pt_lmC", "Pt of negative lepton;Pt;N_{evts}", 50, 0., 200.);
  TH1F* eta_gC = new TH1F("eta_gC", "Eta of Photon", 50, -5.0, 5.0);
  TH1F* eta_lpC = new TH1F("eta_lpC", "Eta of positive lepton;#eta;N_{evts}", 50, -5.0, 5.0);
  TH1F* eta_lmC = new TH1F("eta_lmC", "Eta of negative lepton;#eta;N_{evts}", 50, -5.0, 5.0);
  TH1F* mllgC = new TH1F("mllgC", "Mass of ll-Gamma System;Mass;N_{evts}", 100, 90, 150);
  TH1F* mllC = new TH1F("mllC", "Mass of Lepton Pair;Mass;N_{evts}", 100, 0, 150);

  TH1F* CosT_lpD = new TH1F("CosT_lpD", "Cos(#theta) positive lepton;cos(#theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* CosT_lmD = new TH1F("CosT_lmD", "Cos(#theta) negative lepton;cos(#theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* Phi_lpD = new TH1F("Phi_lpD", "#phi positive lepton;#phi;N_{evts}", 50, -3.2, 3.2);
  TH1F* CosTZGD = new TH1F("CosTZGD", "Cos(#Theta) ZG system;cos(#Theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* Pt_gD = new TH1F("Pt_gD", "Pt of photon in ZG System;Pt;N_{evts}", 50, 0., 200.);
  TH1F* Pt_lpD = new TH1F("Pt_lpD", "Pt of positive lepton;Pt;N_{evts}", 50, 0., 200.);
  TH1F* Pt_lmD = new TH1F("Pt_lmD", "Pt of negative lepton;Pt;N_{evts}", 50, 0., 200.);
  TH1F* eta_gD = new TH1F("eta_gD", "Eta of Photon", 50, -5.0, 5.0);
  TH1F* eta_lpD = new TH1F("eta_lpD", "Eta of positive lepton;#eta;N_{evts}", 50, -5.0, 5.0);
  TH1F* eta_lmD = new TH1F("eta_lmD", "Eta of negative lepton;#eta;N_{evts}", 50, -5.0, 5.0);
  TH1F* mllgD = new TH1F("mllgD", "Mass of ll-Gamma System;Mass;N_{evts}", 100, 90, 150);
  TH1F* mllD = new TH1F("mllD", "Mass of Lepton Pair;Mass;N_{evts}", 100, 0, 150);

  TH1F* CosT_lpE = new TH1F("CosT_lpE", "Cos(#theta) positive lepton;cos(#theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* CosT_lmE = new TH1F("CosT_lmE", "Cos(#theta) negative lepton;cos(#theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* Phi_lpE = new TH1F("Phi_lpE", "#phi positive lepton;#phi;N_{evts}", 50, -3.2, 3.2);
  TH1F* CosTZGE = new TH1F("CosTZGE", "Cos(#Theta) ZG system;cos(#Theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* Pt_gE = new TH1F("Pt_gE", "Pt of photon in ZG System;Pt;N_{evts}", 50, 0., 200.);
  TH1F* Pt_lpE = new TH1F("Pt_lpE", "Pt of positive lepton;Pt;N_{evts}", 50, 0., 200.);
  TH1F* Pt_lmE = new TH1F("Pt_lmE", "Pt of negative lepton;Pt;N_{evts}", 50, 0., 200.);
  TH1F* eta_gE = new TH1F("eta_gE", "Eta of Photon", 50, -5.0, 5.0);
  TH1F* eta_lpE = new TH1F("eta_lpE", "Eta of positive lepton;#eta;N_{evts}", 50, -5.0, 5.0);
  TH1F* eta_lmE = new TH1F("eta_lmE", "Eta of negative lepton;#eta;N_{evts}", 50, -5.0, 5.0);
  TH1F* mllgE = new TH1F("mllgE", "Mass of ll-Gamma System;Mass;N_{evts}", 100, 90, 150);
  TH1F* mllE = new TH1F("mllE", "Mass of Lepton Pair;Mass;N_{evts}", 100, 0, 150);

  TH1F* CosT_lpF = new TH1F("CosT_lpF", "Cos(#theta) positive lepton;cos(#theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* CosT_lmF = new TH1F("CosT_lmF", "Cos(#theta) negative lepton;cos(#theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* Phi_lpF = new TH1F("Phi_lpF", "#phi positive lepton;#phi;N_{evts}", 50, -3.2, 3.2);
  TH1F* CosTZGF = new TH1F("CosTZGF", "Cos(#Theta) ZG system;cos(#Theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* Pt_gF = new TH1F("Pt_gF", "Pt of photon in ZG System;Pt;N_{evts}", 50, 0., 200.);
  TH1F* Pt_lpF = new TH1F("Pt_lpF", "Pt of positive lepton;Pt;N_{evts}", 50, 0., 200.);
  TH1F* Pt_lmF = new TH1F("Pt_lmF", "Pt of negative lepton;Pt;N_{evts}", 50, 0., 200.);
  TH1F* eta_gF = new TH1F("eta_gF", "Eta of Photon", 50, -5.0, 5.0);
  TH1F* eta_lpF = new TH1F("eta_lpF", "Eta of positive lepton;#eta;N_{evts}", 50, -5.0, 5.0);
  TH1F* eta_lmF = new TH1F("eta_lmF", "Eta of negative lepton;#eta;N_{evts}", 50, -5.0, 5.0);
  TH1F* mllgF = new TH1F("mllgF", "Mass of ll-Gamma System;Mass;N_{evts}", 100, 90, 150);
  TH1F* mllF = new TH1F("mllF", "Mass of Lepton Pair;Mass;N_{evts}", 100, 0, 150);

  TH1F* CosT_lpG = new TH1F("CosT_lpG", "Cos(#theta) positive lepton;cos(#theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* CosT_lmG = new TH1F("CosT_lmG", "Cos(#theta) negative lepton;cos(#theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* Phi_lpG = new TH1F("Phi_lpG", "#phi positive lepton;#phi;N_{evts}", 50, -3.2, 3.2);
  TH1F* CosTZGG = new TH1F("CosTZGG", "Cos(#Theta) ZG system;cos(#Theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* Pt_gG = new TH1F("Pt_gG", "Pt of photon in ZG System;Pt;N_{evts}", 50, 0., 200.);
  TH1F* Pt_lpG = new TH1F("Pt_lpG", "Pt of positive lepton;Pt;N_{evts}", 50, 0., 200.);
  TH1F* Pt_lmG = new TH1F("Pt_lmG", "Pt of negative lepton;Pt;N_{evts}", 50, 0., 200.);
  TH1F* eta_gG = new TH1F("eta_gG", "Eta of Photon", 50, -5.0, 5.0);
  TH1F* eta_lpG = new TH1F("eta_lpG", "Eta of positive lepton;#eta;N_{evts}", 50, -5.0, 5.0);
  TH1F* eta_lmG = new TH1F("eta_lmG", "Eta of negative lepton;#eta;N_{evts}", 50, -5.0, 5.0);
  TH1F* mllgG = new TH1F("mllgG", "Mass of ll-Gamma System;Mass;N_{evts}", 100, 90, 150);
  TH1F* mllG = new TH1F("mllG", "Mass of Lepton Pair;Mass;N_{evts}", 100, 0, 150);

  TH1F* CosT_lpH = new TH1F("CosT_lpH", "Cos(#theta) positive lepton;cos(#theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* CosT_lmH = new TH1F("CosT_lmH", "Cos(#theta) negative lepton;cos(#theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* Phi_lpH = new TH1F("Phi_lpH", "#phi positive lepton;#phi;N_{evts}", 50, -3.2, 3.2);
  TH1F* CosTZGH = new TH1F("CosTZGH", "Cos(#Theta) ZG system;cos(#Theta);N_{evts}", 50, -1.1, 1.1);
  TH1F* Pt_gH = new TH1F("Pt_gH", "Pt of photon in ZG System;Pt;N_{evts}", 50, 0., 200.);
  TH1F* Pt_lpH = new TH1F("Pt_lpH", "Pt of positive lepton;Pt;N_{evts}", 50, 0., 200.);
  TH1F* Pt_lmH = new TH1F("Pt_lmH", "Pt of negative lepton;Pt;N_{evts}", 50, 0., 200.);
  TH1F* eta_gH = new TH1F("eta_gH", "Eta of Photon", 50, -5.0, 5.0);
  TH1F* eta_lpH = new TH1F("eta_lpH", "Eta of positive lepton;#eta;N_{evts}", 50, -5.0, 5.0);
  TH1F* eta_lmH = new TH1F("eta_lmH", "Eta of negative lepton;#eta;N_{evts}", 50, -5.0, 5.0);
  TH1F* mllgH = new TH1F("mllgH", "Mass of ll-Gamma System;Mass;N_{evts}", 100, 90, 150);
  TH1F* mllH = new TH1F("mllH", "Mass of Lepton Pair;Mass;N_{evts}", 100, 0, 150);

  TH1F* mz_precut = new TH1F("mz_precut", "Mass of Z Boson (before cut);N_{evts}", 100, 0, 150);
  TH1F* Pb = new TH1F("Pb", "Background Probability (weighted);dXsec_ZGam;N_{evts}", 100, 0, .0001);
  TH1F* logPb = new TH1F("logPb", "Log of Background Probability (weighted);-logPb;N_{evts}", 110, -1.0, 10.0); 
  TH1F* Ps = new TH1F("Ps", "Signal Probability (weighted);dXsec_HZGam;N_{evts}", 100, 0, .0001);
  TH1F* logPs = new TH1F("logPs", "Log of Signal Probability (weighted);-logPs;N_{evts}", 110, -1.0, 10.0);
  TH1F* WD = new TH1F("WD", "Weighted Discriminant;D;N_{evts}", 1000, 0, 1.5);
  TH1F* WD_A = new TH1F("WD_A", "Weighted Discriminant;D;N_{evts}", 1000, 0, 1.5);
  TH1F* WD_B = new TH1F("WD_B", "Weighted Discriminant;D;N_{evts}", 1000, 0, 1.5);
  TH1F* WD_C = new TH1F("WD_C", "Weighted Discriminant;D;N_{evts}", 1000, 0, 1.5);
  TH1F* WD_D = new TH1F("WD_D", "Weighted Discriminant;D;N_{evts}", 1000, 0, 1.5);
  TH1F* mzg_cut1 = new TH1F("mzg_cut1", "Mass of Z-Gamma System;Mass;N_{evts}", 100, 90, 190);
  TH1F* mzg_cut2 = new TH1F("mzg_cut2", "Mass of Z-Gamma System;Mass;N_{evts}", 100, 90, 190);

  TH1F* D_EtaLow = new TH1F("D_EtaLow", "Discriminant;D;N_{evts}", 1000, 0, 1.5);
  TH1F* D_EtaHigh = new TH1F("D_EtaHigh", "Discriminant;D;N_{evts}", 1000, 0, 1.5);

  // Create the instance of TEvtProb to calculate the differential cross-section
  TEvtProb Xcal2;  
  hzgamma_event_type hzgamma_event;

  //==========================================
  // Loop All Events
  //==========================================

  int Ntot = maxevt > ch->GetEntries() ? ch->GetEntries() : maxevt; 
  if ( maxevt < 0. ) Ntot =  ch->GetEntries();
  
  pair<float,float> prob;
  
  if (verbosity >= TVar::INFO) printf("Total number of events = %d\n", Ntot);

  for(int ievt = 0; ievt < Ntot; ievt++){
    if (verbosity >= TVar::INFO && (ievt % 1000 == 0)) 
      std::cout << "Doing Event: " << ievt << std::endl;
    
    //-------------------------------------------
    // initialize the differential cross-sections
    // and other output variables
    //-------------------------------------------
    dXsec_ZGam_MCFM = 0.;
    dXsec_HZGam_MCFM = 0.;
    Discriminant = 0;
    PartCollec.clear();
    TRootLHEFEvent *EvtAccess;

    int N_Neg(0), N_Pos(0), N_Gam(0);

    int PDGID_0(0), PDGID_1(0);

    TLorentzVector p0, p1, p2, psum;
    TVector3 bv;
 
    ch->GetEvent(ievt);

    int Num = Event->GetEntriesFast();
    for (int i = 0; i < Num; i++){
      EvtAccess = (TRootLHEFEvent*)Event->At(i);
    }

    weight = EvtAccess->Weight;

    //Read in the number of particles per Event:
    for (int i = 0; i< EvtAccess->Nparticles; i++)
      {
	PartCollec[i] = (TRootLHEFParticle*)Particle->At(i);
	
	if (PartCollec[i]->PID == 11 || PartCollec[i]->PID == 13){
	  p0.SetPxPyPzE(PartCollec[i]->Px, PartCollec[i]->Py, PartCollec[i]->Pz, PartCollec[i]->E);
	  N_Neg++;
	  PDGID_0 = PartCollec[i]->PID;}
	if (PartCollec[i]->PID == -11 || PartCollec[i]->PID == -13){
	  p1.SetPxPyPzE(PartCollec[i]->Px, PartCollec[i]->Py, PartCollec[i]->Pz, PartCollec[i]->E);
	  N_Pos++;
	  PDGID_1 = PartCollec[i]->PID;}
	if (PartCollec[i]->PID == 22){
	  p2.SetPxPyPzE(PartCollec[i]->Px, PartCollec[i]->Py, PartCollec[i]->Pz, PartCollec[i]->E);
	  N_Gam++;}
      }

    if (N_Neg != 1 || N_Pos != 1 || N_Gam != 1)
      continue;
    if (PDGID_0 != -PDGID_1)
      continue;
    
    //--------------------
    //Calculating ZGAngles
    //--------------------
    ZGLabVectors genLevelInputs;
    ZGAngles genLevelOutputs;

    genLevelInputs.veczg = p0 + p1 + p2;
    genLevelInputs.vecz = p0 + p1;
    genLevelInputs.vecg = p2;
    genLevelInputs.veclm = p0;
    genLevelInputs.veclp = p1;
    
    getZGAngles(genLevelInputs,genLevelOutputs, false);

    mz_precut->Fill(genLevelOutputs.mz, weight);
    
    //----------------------------------------------------
    //Making cuts from H->ZG Analysis Group
    //----------------------------------------------------

    CutTxt << genLevelOutputs.costheta_lp << " " << genLevelOutputs.phi << " " << genLevelOutputs.cosTheta << endl;
    
    if ((genLevelOutputs.ptg < 15.0) ||
        (std::abs(genLevelOutputs.etal1) > 2.5) ||
        (std::abs(genLevelOutputs.etal2) > 2.5) ||
        (std::abs(genLevelOutputs.etag) > 2.5) ||
        (genLevelOutputs.mz < 50.0) ||
        (genLevelOutputs.mzg < 100.0) ||
        (genLevelOutputs.mzg + genLevelOutputs.mz < 185.0) ||
        (std::abs(p0.DeltaR(p2)) < 0.4) ||
        (std::abs(p1.DeltaR(p2)) < 0.4) ||
        ((genLevelOutputs.ptl1 > genLevelOutputs.ptl2) && (genLevelOutputs.ptl1 < 20.0 || genLevelOutputs.ptl2 < 10.0)) ||
        ((genLevelOutputs.ptl2 > genLevelOutputs.ptl1) && (genLevelOutputs.ptl2 < 20.0 || genLevelOutputs.ptl1 < 10.0)))
      {
	if (verbosity >= TVar::DEBUG){
	  std::cout << "Event does not pass cuts:" << endl;

	  CutTxt << "Event #" << ievt << " fails: " << endl;

	  if (genLevelOutputs.ptg < 15.0)
	    CutTxt << "Event fails photon-Pt Cut: Photon Pt = " << genLevelOutputs.ptg << endl;
	  if (std::abs(genLevelOutputs.etal1) > 2.5)
	    CutTxt << "Event fails lepton-1 Eta Cut: Eta = " << genLevelOutputs.etal1 << endl;
	  if (std::abs(genLevelOutputs.etal2) > 2.5)
	    CutTxt << "Event fails lepton-2 Eta Cut: Eta = " << genLevelOutputs.etal2 << endl;
	  if (std::abs(genLevelOutputs.etag) > 2.5)
	    CutTxt << "Event fails Photon Eta Cut: Eta = " << genLevelOutputs.etag << endl;
	  if (genLevelOutputs.mz < 50.0)
	    CutTxt << "Event fails z-mass cut: mz = " << genLevelOutputs.mz << endl;
	  if (genLevelOutputs.mzg < 100.0)
	    CutTxt << "Event fails H-mass cut: mzg = " << genLevelOutputs.mzg << endl;
	  if (genLevelOutputs.mzg + genLevelOutputs.mz < 185.0)
	    CutTxt << "Event fails H+z-mass cut: mzg + mz = " << genLevelOutputs.mzg + genLevelOutputs.mz << endl;
	  if (std::abs(p0.DeltaR(p2)) < 0.4)
	    CutTxt << "Event fails lepton-1/gamma DR Cut: DR = " << p0.DeltaR(p2) << endl;
	  if (std::abs(p1.DeltaR(p2)) < 0.4)
	    CutTxt << "Event fails lepton-2/gamma DR Cut: DR = " << p1.DeltaR(p2) << endl;
	  if ((genLevelOutputs.ptl1 > genLevelOutputs.ptl2) && (genLevelOutputs.ptl1 < 20.0 || genLevelOutputs.ptl2 < 10.0))
	    CutTxt << "Event fails Lepton-Pt Cut: Pt1 = " << genLevelOutputs.ptl1 << " Pt2 = " << genLevelOutputs.ptl2 << endl;
	  if ((genLevelOutputs.ptl2 > genLevelOutputs.ptl1) && (genLevelOutputs.ptl2 < 20.0 || genLevelOutputs.ptl1 < 10.0))
	    CutTxt << "Event fails Lepton-Pt Cut: Pt1 = " << genLevelOutputs.ptl1 << " Pt2 = " << genLevelOutputs.ptl2 << endl;
	}
	continue;
      }

    //-------------------------------------------------
    //        Filling User-Defined Histograms
    //-------------------------------------------------
    CosT_lp->Fill(genLevelOutputs.costheta_lp, weight);
    CosT_lm->Fill(genLevelOutputs.costheta_lm, weight);
    Phi_lp->Fill(genLevelOutputs.phi, weight);
    CosTZG->Fill(genLevelOutputs.cosTheta, weight);
    Pt_g->Fill(genLevelOutputs.ptg, weight);
    Pt_lp->Fill(genLevelOutputs.ptl1, weight);
    Pt_lm->Fill(genLevelOutputs.ptl2, weight);
    eta_g->Fill(genLevelOutputs.etag, weight);
    eta_lp->Fill(genLevelOutputs.etal1, weight);
    eta_lm->Fill(genLevelOutputs.etal2, weight);
    mzg->Fill(genLevelOutputs.mzg, weight);
    mz->Fill(genLevelOutputs.mz, weight);

    if (genLevelOutputs.mzg > 120.0 && genLevelOutputs.mzg < 125.0){
      CosT_lpA->Fill(genLevelOutputs.costheta_lp, weight);
      CosT_lmA->Fill(genLevelOutputs.costheta_lm, weight);
      Phi_lpA->Fill(genLevelOutputs.phi, weight);
      CosTZGA->Fill(genLevelOutputs.cosTheta, weight);
      Pt_gA->Fill(genLevelOutputs.ptg, weight);
      Pt_lpA->Fill(genLevelOutputs.ptl1, weight);
      Pt_lmA->Fill(genLevelOutputs.ptl2, weight);
      eta_gA->Fill(genLevelOutputs.etag, weight);
      eta_lpA->Fill(genLevelOutputs.etal1, weight);
      eta_lmA->Fill(genLevelOutputs.etal2, weight);
      mllgA->Fill(genLevelOutputs.mzg, weight);
      mllA->Fill(genLevelOutputs.mz, weight);
    }

    if (genLevelOutputs.mzg > 125.0 && genLevelOutputs.mzg < 130.0){
      CosT_lpB->Fill(genLevelOutputs.costheta_lp, weight);
      CosT_lmB->Fill(genLevelOutputs.costheta_lm, weight);
      Phi_lpB->Fill(genLevelOutputs.phi, weight);
      CosTZGB->Fill(genLevelOutputs.cosTheta, weight);
      Pt_gB->Fill(genLevelOutputs.ptg, weight);
      Pt_lpB->Fill(genLevelOutputs.ptl1, weight);
      Pt_lmB->Fill(genLevelOutputs.ptl2, weight);
      eta_gB->Fill(genLevelOutputs.etag, weight);
      eta_lpB->Fill(genLevelOutputs.etal1, weight);
      eta_lmB->Fill(genLevelOutputs.etal2, weight);
      mllgB->Fill(genLevelOutputs.mzg, weight);
      mllB->Fill(genLevelOutputs.mz, weight);
    }

    if (genLevelOutputs.mzg > 130.0 && genLevelOutputs.mzg < 135.0){
      CosT_lpC->Fill(genLevelOutputs.costheta_lp, weight);
      CosT_lmC->Fill(genLevelOutputs.costheta_lm, weight);
      Phi_lpC->Fill(genLevelOutputs.phi, weight);
      CosTZGC->Fill(genLevelOutputs.cosTheta, weight);
      Pt_gC->Fill(genLevelOutputs.ptg, weight);
      Pt_lpC->Fill(genLevelOutputs.ptl1, weight);
      Pt_lmC->Fill(genLevelOutputs.ptl2, weight);
      eta_gC->Fill(genLevelOutputs.etag, weight);
      eta_lpC->Fill(genLevelOutputs.etal1, weight);
      eta_lmC->Fill(genLevelOutputs.etal2, weight);
      mllgC->Fill(genLevelOutputs.mzg, weight);
      mllC->Fill(genLevelOutputs.mz, weight);
    }

    if (genLevelOutputs.mzg > 135.0 && genLevelOutputs.mzg < 140.0){
      CosT_lpD->Fill(genLevelOutputs.costheta_lp, weight);
      CosT_lmD->Fill(genLevelOutputs.costheta_lm, weight);
      Phi_lpD->Fill(genLevelOutputs.phi, weight);
      CosTZGD->Fill(genLevelOutputs.cosTheta, weight);
      Pt_gD->Fill(genLevelOutputs.ptg, weight);
      Pt_lpD->Fill(genLevelOutputs.ptl1, weight);
      Pt_lmD->Fill(genLevelOutputs.ptl2, weight);
      eta_gD->Fill(genLevelOutputs.etag, weight);
      eta_lpD->Fill(genLevelOutputs.etal1, weight);
      eta_lmD->Fill(genLevelOutputs.etal2, weight);
      mllgD->Fill(genLevelOutputs.mzg, weight);
      mllD->Fill(genLevelOutputs.mz, weight);
    }

    if (genLevelOutputs.mzg > 140.0 && genLevelOutputs.mzg < 145.0){
      CosT_lpE->Fill(genLevelOutputs.costheta_lp, weight);
      CosT_lmE->Fill(genLevelOutputs.costheta_lm, weight);
      Phi_lpE->Fill(genLevelOutputs.phi, weight);
      CosTZGE->Fill(genLevelOutputs.cosTheta, weight);
      Pt_gE->Fill(genLevelOutputs.ptg, weight);
      Pt_lpE->Fill(genLevelOutputs.ptl1, weight);
      Pt_lmE->Fill(genLevelOutputs.ptl2, weight);
      eta_gE->Fill(genLevelOutputs.etag, weight);
      eta_lpE->Fill(genLevelOutputs.etal1, weight);
      eta_lmE->Fill(genLevelOutputs.etal2, weight);
      mllgE->Fill(genLevelOutputs.mzg, weight);
      mllE->Fill(genLevelOutputs.mz, weight);
    }
    if (genLevelOutputs.mzg > 145.0 && genLevelOutputs.mzg < 150){
      CosT_lpF->Fill(genLevelOutputs.costheta_lp, weight);
      CosT_lmF->Fill(genLevelOutputs.costheta_lm, weight);
      Phi_lpF->Fill(genLevelOutputs.phi, weight);
      CosTZGF->Fill(genLevelOutputs.cosTheta, weight);
      Pt_gF->Fill(genLevelOutputs.ptg, weight);
      Pt_lpF->Fill(genLevelOutputs.ptl1, weight);
      Pt_lmF->Fill(genLevelOutputs.ptl2, weight);
      eta_gF->Fill(genLevelOutputs.etag, weight);
      eta_lpF->Fill(genLevelOutputs.etal1, weight);
      eta_lmF->Fill(genLevelOutputs.etal2, weight);
      mllgF->Fill(genLevelOutputs.mzg, weight);
      mllF->Fill(genLevelOutputs.mz, weight);
    }
    if (genLevelOutputs.mzg > 150.0 && genLevelOutputs.mzg < 155){
      CosT_lpG->Fill(genLevelOutputs.costheta_lp, weight);
      CosT_lmG->Fill(genLevelOutputs.costheta_lm, weight);
      Phi_lpG->Fill(genLevelOutputs.phi, weight);
      CosTZGG->Fill(genLevelOutputs.cosTheta, weight);
      Pt_gG->Fill(genLevelOutputs.ptg, weight);
      Pt_lpG->Fill(genLevelOutputs.ptl1, weight);
      Pt_lmG->Fill(genLevelOutputs.ptl2, weight);
      eta_gG->Fill(genLevelOutputs.etag, weight);
      eta_lpG->Fill(genLevelOutputs.etal1, weight);
      eta_lmG->Fill(genLevelOutputs.etal2, weight);
      mllgG->Fill(genLevelOutputs.mzg, weight);
      mllG->Fill(genLevelOutputs.mz, weight);
    }
    if (genLevelOutputs.mzg > 155.0 && genLevelOutputs.mzg < 160){
      CosT_lpH->Fill(genLevelOutputs.costheta_lp, weight);
      CosT_lmH->Fill(genLevelOutputs.costheta_lm, weight);
      Phi_lpH->Fill(genLevelOutputs.phi, weight);
      CosTZGH->Fill(genLevelOutputs.cosTheta, weight);
      Pt_gH->Fill(genLevelOutputs.ptg, weight);
      Pt_lpH->Fill(genLevelOutputs.ptl1, weight);
      Pt_lmH->Fill(genLevelOutputs.ptl2, weight);
      eta_gH->Fill(genLevelOutputs.etag, weight);
      eta_lpH->Fill(genLevelOutputs.etal1, weight);
      eta_lmH->Fill(genLevelOutputs.etal2, weight);
      mllgH->Fill(genLevelOutputs.mzg, weight);
      mllH->Fill(genLevelOutputs.mz, weight);
    }


    //----------------------------------------------------------------------
    //     Setting Up hzgamma_event for MCFM Matrix Element Calculation
    //----------------------------------------------------------------------

    //Boost TLorentzVectors to ZG-Restframe for MCFM
    psum = p0 + p1 + p2;
    bv = -psum.BoostVector();
    p0.Boost(bv);
    p1.Boost(bv);
    p2.Boost(bv);

    TLorentzVector Z_minus = p0;
    TLorentzVector Z_plus  = p1;
    TLorentzVector Gamma = p2;

    hzgamma_event.p[0].SetPxPyPzE(Z_minus.Px(), Z_minus.Py(), Z_minus.Pz(), Z_minus.Energy());
    hzgamma_event.p[1].SetPxPyPzE(Z_plus.Px(), Z_plus.Py(), Z_plus.Pz(), Z_plus.Energy());
    hzgamma_event.p[2].SetPxPyPzE(Gamma.Px(), Gamma.Py(), Gamma.Pz(), Gamma.Energy());

    hzgamma_event.PdgCode[0] = PDGID_0;
    hzgamma_event.PdgCode[1] = PDGID_1;
    hzgamma_event.PdgCode[2] = 22;

    if (PDGID_0 == 11 || PDGID_0 == -11){
      if (PDGID_0 == -11){
	me_pg->Fill((hzgamma_event.p[0] + hzgamma_event.p[2]).M());
	me_ng->Fill((hzgamma_event.p[1] + hzgamma_event.p[2]).M());
      }
      else{
	me_pg->Fill((hzgamma_event.p[1] + hzgamma_event.p[2]).M());
	me_ng->Fill((hzgamma_event.p[0] + hzgamma_event.p[2]).M());
      }
      megVmee->Fill((hzgamma_event.p[0]+hzgamma_event.p[2]).M(),(hzgamma_event.p[0]+hzgamma_event.p[1]).M());
      megVmee->Fill((hzgamma_event.p[1]+hzgamma_event.p[2]).M(),(hzgamma_event.p[0]+hzgamma_event.p[1]).M());
      meegVmeg->Fill((hzgamma_event.p[0]+hzgamma_event.p[1]+hzgamma_event.p[2]).M(), (hzgamma_event.p[0]+hzgamma_event.p[2]).M());
      meegVmeg->Fill((hzgamma_event.p[0]+hzgamma_event.p[1]+hzgamma_event.p[2]).M(), (hzgamma_event.p[1]+hzgamma_event.p[2]).M());
    }
    if (PDGID_0 == 13 || PDGID_0 == -13){
      if (PDGID_0 == -13){
	mm_pg->Fill((hzgamma_event.p[0] + hzgamma_event.p[2]).M());
	mm_ng->Fill((hzgamma_event.p[1] + hzgamma_event.p[2]).M());
      }
      else{
	mm_pg->Fill((hzgamma_event.p[1] + hzgamma_event.p[2]).M());
	mm_ng->Fill((hzgamma_event.p[0] + hzgamma_event.p[2]).M());
      }
      mmgVmmm->Fill((hzgamma_event.p[0]+hzgamma_event.p[2]).M(),(hzgamma_event.p[0]+hzgamma_event.p[1]).M());
      mmgVmmm->Fill((hzgamma_event.p[1]+hzgamma_event.p[2]).M(),(hzgamma_event.p[0]+hzgamma_event.p[1]).M());
      mmmgVmmg->Fill((hzgamma_event.p[0]+hzgamma_event.p[1]+hzgamma_event.p[2]).M(), (hzgamma_event.p[0]+hzgamma_event.p[2]).M());
      mmmgVmmg->Fill((hzgamma_event.p[0]+hzgamma_event.p[1]+hzgamma_event.p[2]).M(), (hzgamma_event.p[1]+hzgamma_event.p[2]).M());
    }

    
    float zmass = (hzgamma_event.p[0]+hzgamma_event.p[1]).M();
    float gammass = (hzgamma_event.p[2]).M();
    float zgammass = (hzgamma_event.p[0]+hzgamma_event.p[1]+hzgamma_event.p[2]).M();    

    if (verbosity >= TVar::DEBUG) {
      cout << "\n=========================================================\n";
      cout << "Entry: " << ievt << "\n";
      cout << "Input: ==================================================" <<endl;
      printf("lep1 (Px, Py, Pz, E) = (%4.4f, %4.4f, %4.4f, %4.4f)\n",  p0.Px(), p0.Py(), p0.Pz(), p0.E());
      printf("lep2 (Px, Py, Pz, E) = (%4.4f, %4.4f, %4.4f, %4.4f)\n",  p1.Px(), p1.Py(), p1.Pz(), p1.E()); 
      printf("gam (Px, Py, Pz, E) = (%4.4f, %4.4f, %4.4f, %4.4f)\n",  p2.Px(), p2.Py(), p2.Pz(), p2.E());
      std::cout << "ZGam system (pX, pY, pZ, E, mass) = ( " 
		<< (hzgamma_event.p[0]+hzgamma_event.p[1]+hzgamma_event.p[2]).Px() << ", "
		<< (hzgamma_event.p[0]+hzgamma_event.p[1]+hzgamma_event.p[2]).Py() << ", "
		<< (hzgamma_event.p[0]+hzgamma_event.p[1]+hzgamma_event.p[2]).Pz() << ", "
		<< (hzgamma_event.p[0]+hzgamma_event.p[1]+hzgamma_event.p[2]).Energy()  << ", "
		<< zgammass << ")\n";
      std::cout << "Z mass = " << zmass << "\tgammass = " << gammass << "\n";
      cout << "=========================================================\n";
    } 
    // finish loading event information
    
    // ==== Begin the differential cross-section calculation ====
    Xcal2.SetHiggsMass(zgammass);
    Xcal2.SetMatrixElement(TVar::MCFM);
    dXsec_ZGam_MCFM = Xcal2.XsecCalc(TVar::qqb_zgam, TVar::QQB, hzgamma_event,verbose);
    dXsec_HZGam_MCFM = Xcal2.XsecCalc(TVar::gg_hzgam, TVar::GG, hzgamma_event,verbose);
    Pb->Fill(dXsec_ZGam_MCFM, weight);
    logPb->Fill(-log10(dXsec_ZGam_MCFM),weight);
    Ps->Fill(dXsec_HZGam_MCFM, weight);
    logPs->Fill(-log10(dXsec_HZGam_MCFM),weight);

    logBkg = -log10(dXsec_ZGam_MCFM);
    logSig = -log10(dXsec_HZGam_MCFM);
    Discriminant = -log(dXsec_ZGam_MCFM/(dXsec_ZGam_MCFM+dXsec_HZGam_MCFM));
    WD->Fill(Discriminant, weight);

    if (std::abs(genLevelOutputs.etag) < 1.44 && std::abs(genLevelOutputs.etag) != 0)
      D_EtaLow->Fill(Discriminant, weight);
    else if (std::abs(genLevelOutputs.etag) > 1.57 && std::abs(genLevelOutputs.etag) < 2.5)
      D_EtaHigh->Fill(Discriminant, weight);

    DvMlla->Fill(Discriminant, genLevelOutputs.mzg);
    PsvMlla->Fill(dXsec_HZGam_MCFM, genLevelOutputs.mzg);
    if (Discriminant > 0.02)
      Mlla_PostDCut->Fill(genLevelOutputs.mzg, weight);

    if (genLevelOutputs.mzg >= 120.0 && genLevelOutputs.mzg < 125.0)
      WD_A->Fill(Discriminant, weight);
    if (genLevelOutputs.mzg >= 125.0 && genLevelOutputs.mzg < 130.0)
      WD_B->Fill(Discriminant, weight);
    if (genLevelOutputs.mzg >= 130.0 && genLevelOutputs.mzg < 135.0)
      WD_C->Fill(Discriminant, weight);
    if (genLevelOutputs.mzg >= 135.0 && genLevelOutputs.mzg < 140.0)
      WD_D->Fill(Discriminant, weight);

    if (Discriminant > 0.02)
      mzg_cut1->Fill(genLevelOutputs.mzg, weight);
    if (Discriminant > 0.063)
      mzg_cut2->Fill(genLevelOutputs.mzg, weight);

    KKTxt << Discriminant;
    if (Signal)
      KKTxt << " 3" << endl;
    else
      KKTxt << " 4" << endl;

    evt_tree->Fill();
    
  }//nevent
  
  cout << "TotCalculatedEvts: " << evt_tree->GetEntries() <<endl 
       << "Fraction Cut On: " << (double)evt_tree->GetEntries()/Ntot;
  newfile->cd();
  evt_tree->Write(); 
  newfile->Close();
  histfile->cd();
  histfile->Write();
  histfile->Close();
  KKTxt.close();
  CutTxt.close();
}  
