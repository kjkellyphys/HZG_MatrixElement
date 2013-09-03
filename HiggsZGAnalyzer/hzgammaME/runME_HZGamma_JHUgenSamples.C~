/// 
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
  bool Signal = false;
  KKTxt.open ("./TxtOutput/Discriminant.txt");

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
  TTree* ch=(TTree*)fin->Get("h300"); 
  if (ch==0x0) return; 
  TTree* evt_tree=(TTree*) ch->CloneTree(0, "fast");
  evt_tree->SetName("newTree");

  //------------------------------
  //       Input Variables
  //------------------------------

  //Weight from MCFM-Produced Files
  float weight = 0.;
  weight = 1.0; //For uniform weight.
  ch->SetBranchAddress("wt", &weight);
  
  //For Reading in components of TLorentzVectors individually.
  /*
  float l_minus_px(0.), l_minus_py(0.), l_minus_pz(0.), l_minus_e(0.);
  float l_plus_px(0.), l_plus_py(0.), l_plus_pz(0.), l_plus_e(0.);
  float gamma_px(0.), gamma_py(0.), gamma_pz(0.), gamma_e(0.);

  ch->SetBranchAddress("px3", &l_minus_px);
  ch->SetBranchAddress("py3", &l_minus_py);
  ch->SetBranchAddress("pz3", &l_minus_pz);
  ch->SetBranchAddress("E3", &l_minus_e);
  ch->SetBranchAddress("px4", &l_plus_px);
  ch->SetBranchAddress("py4", &l_plus_py);
  ch->SetBranchAddress("pz4", &l_plus_pz);
  ch->SetBranchAddress("E4", &l_plus_e);
  ch->SetBranchAddress("px5", &gamma_px);
  ch->SetBranchAddress("py5", &gamma_py);
  ch->SetBranchAddress("pz5", &gamma_pz);
  ch->SetBranchAddress("E5", &gamma_e);
  */
  //For Reading in TLorentzVectors as a whole
  
  TBranch *b_l_minus; TBranch *b_l_plus; TBranch *b_gamma;
  TLorentzVector *l_minus_event, *l_plus_event, *gamma_event;
  l_minus_event = 0; l_plus_event = 0; gamma_event = 0;

  ch->SetBranchAddress("l_minus", &l_minus_event, &b_l_minus);
  ch->SetBranchAddress("l_plus", &l_plus_event, &b_l_plus);
  ch->SetBranchAddress("gamma", &gamma_event, &b_gamma);
  

  //------------------------
  //      Output Files
  //------------------------

  // output variables (in _ME.root file)
  float dXsec_ZGam_MCFM = 0.;
  float dXsec_HZGam_MCFM = 0.;
  float Discriminant = 0.;
  float PreBoostMass = 0.;
  float PostBoostMass = 0.;
  float logBkg(0.), logSig(0.);
  
  evt_tree->Branch("dxSec_ZGam_MCFM"   , &dXsec_ZGam_MCFM ,"dXsec_ZGam_MCFM/F");
  evt_tree->Branch("dxSec_HZGam_MCFM"   , &dXsec_HZGam_MCFM ,"dXsec_HZGam_MCFM/F");
  evt_tree->Branch("logBkg", &logBkg, "logBkg/F");
  evt_tree->Branch("logSig", &logSig, "logSig/F");
  evt_tree->Branch("Discriminant", &Discriminant, "Discriminant/F");
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
  TH1F* mzg = new TH1F("mzg", "Mass of Z-Gamma System;Mass;N_{evts}", 100, 90, 150);
  TH1F* mz = new TH1F("mz", "Mass of Z Boson;Mass;N_{evts}", 100, 0, 150);
  TH1F* mz_precut = new TH1F("mz_precut", "Mass of Z Boson (before cut);N_{evts}", 100, 0, 150);
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
    PreBoostMass = 0;
    PostBoostMass = 0;

    ch->GetEntry(ievt);

    TLorentzVector p0, p1, p2, psum;
    TVector3 bv;

    //When Pulled from Components
    //p0.SetPxPyPzE(l_minus_px, l_minus_py, l_minus_pz, l_minus_e);
    //p1.SetPxPyPzE(l_plus_px, l_plus_py, l_plus_pz, l_plus_e);
    //p2.SetPxPyPzE(gamma_px, gamma_py, gamma_pz, gamma_e);

    //When Pulled from TLorentzVector
    p0.SetPxPyPzE(l_minus_event->Px(), l_minus_event->Py(), l_minus_event->Pz(), l_minus_event->Energy());
    p1.SetPxPyPzE(l_plus_event->Px(), l_plus_event->Py(), l_plus_event->Pz(), l_plus_event->Energy());
    p2.SetPxPyPzE(gamma_event->Px(), gamma_event->Py(), gamma_event->Pz(), gamma_event->Energy());
    
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

    hzgamma_event.PdgCode[0] = 11;
    hzgamma_event.PdgCode[1] = -11;
    hzgamma_event.PdgCode[2] = 22;
    
    float zmass = (hzgamma_event.p[0]+hzgamma_event.p[1]).M();
    float gammass = (hzgamma_event.p[2]).M();
    float zgammass = (hzgamma_event.p[0]+hzgamma_event.p[1]+hzgamma_event.p[2]).M();    
    PostBoostMass = zgammass;

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

    KKTxt << Discriminant;
    if (Signal)
      KKTxt << " 0" << endl;
    else
      KKTxt << " 1" << endl;

    evt_tree->Fill();
    
  }//nevent
  
  cout << "TotCalculatedEvts: " << evt_tree->GetEntries() <<endl; 
  cout << "Fraction Cut On  : " << (double)evt_tree->GetEntries()/Ntot;
  newfile->cd();
  evt_tree->Write(); 
  newfile->Close();
  histfile->cd();
  histfile->Write();
  histfile->Close();
  KKTxt.close();
  CutTxt.close();
}  
