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
  KKTxt.open ("Disc_Sig.txt");

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
  float Discriminant2 = 0.;
  float PreBoostMass = 0.;
  float PostBoostMass = 0.;
  float logBkg(0.), logSig(0.);
  
  evt_tree->Branch("dxSec_ZGam_MCFM"   , &dXsec_ZGam_MCFM ,"dXsec_ZGam_MCFM/F");
  evt_tree->Branch("dxSec_HZGam_MCFM"   , &dXsec_HZGam_MCFM ,"dXsec_HZGam_MCFM/F");
  evt_tree->Branch("logBkg", &logBkg, "logBkg/F");
  evt_tree->Branch("logSig", &logSig, "logSig/F");
  evt_tree->Branch("Discriminant", &Discriminant, "Discriminant/F");
  evt_tree->Branch("Discriminant2", &Discriminant2, "Discriminant2/F");  
  evt_tree->Branch("PreBoostMass", &PreBoostMass, "PreBoostMass/F");
  evt_tree->Branch("PostBoostMass", &PostBoostMass, "PostBoostMass/F");
  
  //User-defined Histograms (_AnglePlots.root file)
  histfile->cd();
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
  TH1F* mllg = new TH1F("mllg", "Mass of ll-Gamma System;Mass;N_{evts}", 100, 90, 150);
  TH1F* mll = new TH1F("mll", "Mass of Lepton Pair;Mass;N_{evts}", 100, 0, 150);
  TH1F* mll_precut = new TH1F("mll_precut", "Mass of Lepton Pair (before cut);N_{evts}", 100, 0, 150);
  TH1F* Pb = new TH1F("Pb", "Background Probability (weighted);dXsec_ZGam;N_{evts}", 100, 0, .0001);
  TH1F* logPb = new TH1F("logPb", "Log of Background Probability (weighted);-logPb;N_{evts}", 110, -1.0, 10.0); 
  TH1F* Ps = new TH1F("Ps", "Signal Probability (weighted);dXsec_HZGam;N_{evts}", 100, 0, .0001);
  TH1F* logPs = new TH1F("logPs", "Log of Signal Probability (weighted);-logPs;N_{evts}", 110, -1.0, 10.0);
  TH1F* WD = new TH1F("WD", "Weighted Discriminant;D;N_{evts}", 1000, 0, 5);

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
    Discriminant2 = 0;
    PreBoostMass = 0;
    PostBoostMass = 0;
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

    mll_precut->Fill(genLevelOutputs.mz, weight);

    if (genLevelOutputs.mz > 50.0)
      continue;
    
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
    mllg->Fill(genLevelOutputs.mzg, weight);
    mll->Fill(genLevelOutputs.mz, weight);

    //----------------------------------------------------------------------
    //     Setting Up hzgamma_event for MCFM Matrix Element Calculation
    //----------------------------------------------------------------------

    //Boost TLorentzVectors to ZG-Restframe for MCFM
    psum = p0 + p1 + p2;
    PreBoostMass = psum.M();
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
    
    float llmass = (hzgamma_event.p[0]+hzgamma_event.p[1]).M();
    float gammass = (hzgamma_event.p[2]).M();
    float llgammass = (hzgamma_event.p[0]+hzgamma_event.p[1]+hzgamma_event.p[2]).M();    
    PostBoostMass = llgammass;

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
		<< llgammass << ")\n";
      std::cout << "Z mass = " << llmass << "\tgammass = " << gammass << "\n";
      cout << "=========================================================\n";
    } 
    // finish loading event information
    
    // ==== Begin the differential cross-section calculation ====
    Xcal2.SetHiggsMass(llgammass);
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
    Discriminant2 = -log(dXsec_HZGam_MCFM/(dXsec_ZGam_MCFM+dXsec_HZGam_MCFM));
    WD->Fill(Discriminant, weight);

    KKTxt << Discriminant;
    if (Signal)
      KKTxt << " 0" << endl;
    else
      KKTxt << " 1" << endl;

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
}  
