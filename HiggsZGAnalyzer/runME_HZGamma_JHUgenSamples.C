// 
// This code calculates the ME
// 

#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
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
#include "plugins/ZGAngles.cc"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector; 

float ERRORthreshold=1.0;
using namespace std;

void xseccalc(TString inputDir, TString fileName, TString outputDir, int maxevt, TVar::VerbosityLevel verbosity);

vector<TLorentzVector> Calculate4Momentum(float Mx,float M1,float M2,float theta,float theta1,float theta2,float Phi1,float Phi)
{
    float phi1,phi2;
    phi1=TMath::Pi()-Phi1;
    phi2=Phi1+Phi;
    
    
    float gamma1,gamma2,beta1,beta2;
    
    gamma1=(Mx*Mx+M1*M1-M2*M2)/(2*Mx*M1);
    gamma2=(Mx*Mx-M1*M1+M2*M2)/(2*Mx*M2);
    beta1=sqrt(1-1/(gamma1*gamma1));
    beta2=sqrt(1-1/(gamma2*gamma2));
    
    
    //incoming particle 4 vectors
    TLorentzVector p1CM(0,0,Mx/2,Mx/2);
    TLorentzVector p2CM(0,0,-Mx/2,Mx/2);
    
    //vector boson 4 vectors
    TLorentzVector kZ(gamma1*M1*sin(theta)*beta1,0, gamma1*M1*cos(theta)*beta1,gamma1*M1*1);   
    TLorentzVector kG(-gamma2*M2*sin(theta)*beta2,0, -gamma2*M2*cos(theta)*beta2,gamma2*M2*1);
    
    //Rotation and Boost matrices. Note gamma1*beta1*M1=gamma2*beta2*M2.
    
    TLorentzRotation Z1ToZ,Z2ToZ;
    
    Z1ToZ.Boost(0,0,beta1);
    Z2ToZ.Boost(0,0,beta2);
    Z1ToZ.RotateY(theta);
    Z2ToZ.RotateY(TMath::Pi()+theta);
    
    
    //fermons 4 vectors in vector boson rest frame
    
    TLorentzVector p3Z((M1/2)*sin(theta1)*cos(phi1),(M1/2)*sin(theta1)*sin(phi1),(M1/2)*cos(theta1),(M1/2)*1);
       
    TLorentzVector p4Z(-(M1/2)*sin(theta1)*cos(phi1),-(M1/2)*sin(theta1)*sin(phi1),-(M1/2)*cos(theta1),(M1/2)*1);
      
    TLorentzVector p5G((M2)*sin(theta2)*cos(phi2),(M2)*sin(theta2)*sin(phi2),(M2)*cos(theta2),(M2)*1);
    
    //TLorentzVector p6Z2(-(M2/2)*sin(theta2)*cos(phi2),-(M2/2)*sin(theta2)*sin(phi2),-(M2/2)*cos(theta2),(M2/2)*1);
      
    
    // fermions 4 vectors in CM frame
    
    TLorentzVector p3CM,p4CM,p5CM,p6CM;
    
    p3CM=Z1ToZ*p3Z;
    p4CM=Z1ToZ*p4Z;
    p5CM=Z2ToZ*p5G;
    //p6CM=Z2ToZ*p6Z2;
    
    vector<TLorentzVector> p;
    
    p.push_back(p3CM);
    p.push_back(p4CM);
    p.push_back(p5CM);
    //p.push_back(p6CM);http://www.chiralcomp.com/support/mixing_f77_c_cpp/

    return p;
}



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

  TString outFileName = outputDir+fileName;
  outFileName.ReplaceAll(".root","_ME.root");
  cout << outFileName <<endl;

  TFile *newfile = new TFile(outFileName,"recreate");

  TTree* ch=(TTree*)fin->Get("h300"); 
  if (ch==0x0) return; 
  TTree* evt_tree=(TTree*) ch->CloneTree(0, "fast");
  evt_tree->SetName("newTree");
  //For Reading in components of TLorentzVectors individually.

  float l_minus_px(0.), l_minus_py(0.), l_minus_pz(0.), l_minus_e(0.);
  float l_plus_px(0.), l_plus_py(0.), l_plus_pz(0.), l_plus_e(0.);
  float gamma_px(0.), gamma_py(0.), gamma_pz(0.), gamma_e(0.);

  ch->SetBranchAddress("l_minus_px", &l_minus_px);
  ch->SetBranchAddress("l_minus_py", &l_minus_py);
  ch->SetBranchAddress("l_minus_pz", &l_minus_pz);
  ch->SetBranchAddress("l_minus_e", &l_minus_e);
  ch->SetBranchAddress("l_plus_px", &l_plus_px);
  ch->SetBranchAddress("l_plus_py", &l_plus_py);
  ch->SetBranchAddress("l_plus_pz", &l_plus_pz);
  ch->SetBranchAddress("l_plus_e", &l_plus_e);
  ch->SetBranchAddress("gamma_px", &gamma_px);
  ch->SetBranchAddress("gamma_py", &gamma_py);
  ch->SetBranchAddress("gamma_pz", &gamma_pz);
  ch->SetBranchAddress("gamma_e", &gamma_e);

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

  

  //For Reading in TLorentzVectors as a whole
  /*
  TBranch *b_l_minus; TBranch *b_l_plus; TBranch *b_gamma;
  TLorentzVector *l_minus_event, *l_plus_event, *gamma_event;
  l_minus_event = 0; l_plus_event = 0; gamma_event = 0;

  ch->SetBranchAddress("l_minus", &l_minus_event, &b_l_minus);
  ch->SetBranchAddress("l_plus", &l_plus_event, &b_l_plus);
  ch->SetBranchAddress("gamma", &gamma_event, &b_gamma);
  */
  // output variables
  //float dXsec_ZZ_MCFM = 0.;
  //float dXsec_HZZ_MCFM = 0.;
  float dXsec_ZGam_MCFM = 0.;
  float dXsec_HZGam_MCFM = 0.;
  float Discriminant = 0.;
  float PreBoostMass = 0.;
  float PostBoostMass = 0.;
  float logBkg(0.), logSig(0.);
  float CosT_lp = 0.;
  float CosT_lm = 0.;
  float Phi_lp = 0.;
  float CosTZG = 0.;
  float PtG = 0.;
  float Eta_lp = 0.;
  float Eta_lm = 0.;
  float Eta_g = 0.;

  //evt_tree->Branch("dXsec_ZZ_MCFM"   , &dXsec_ZZ_MCFM   ,"dXsec_ZZ_MCFM/F");
  //evt_tree->Branch("dXsec_HZZ_MCFM"  , &dXsec_HZZ_MCFM   ,"dXsec_HZZ_MCFM/F");
  evt_tree->Branch("dxSec_ZGam_MCFM"   , &dXsec_ZGam_MCFM ,"dXsec_ZGam_MCFM/F");
  evt_tree->Branch("dxSec_HZGam_MCFM"   , &dXsec_HZGam_MCFM ,"dXsec_HZGam_MCFM/F");
  evt_tree->Branch("logBkg", &logBkg, "logBkg/F");
  evt_tree->Branch("logSig", &logSig, "logSig/F");
  evt_tree->Branch("Discriminant", &Discriminant, "Discriminant/F");
  evt_tree->Branch("PreBoostMass", &PreBoostMass, "PreBoostMass/F");
  evt_tree->Branch("PostBoostMass", &PostBoostMass, "PostBoostMass/F");
  evt_tree->Branch("CosT_lp", &CosT_lp, "Cos(#theta) Positive Lepton/F");
  evt_tree->Branch("CosT_lm", &CosT_lm, "Cos(#theta) Negative Lepton/F");
  evt_tree->Branch("Phi_lp", &Phi_lp, "#phi of Positive Lepton/F");
  evt_tree->Branch("CosTZG", &CosTZG, "Cos(#Theta) of ZG System/F");
  evt_tree->Branch("PtG", &PtG, "Pt of Photon/F");
  evt_tree->Branch("Eta_lp", &Eta_lp, "Eta of Positive Lepton/F");
  evt_tree->Branch("Eta_lm", &Eta_lm, "Eta of Negative Lepton/F");
  evt_tree->Branch("Eta_g", &Eta_g, "Eta of Photon/F");
  

// input variables 
  /*
  float m1,m2,h1,h2,hs,phi,phi1,mzg;  
  ch->SetBranchAddress( "ZMass"        , &m1      );   
  ch->SetBranchAddress( "GamMass"        , &m2      );   
  ch->SetBranchAddress( "helcosthetaZ" , &h1      );   
  ch->SetBranchAddress( "helcosthetaGam" , &h2      );   
  ch->SetBranchAddress( "costhetastar"  , &hs      );   
  ch->SetBranchAddress( "helphi"        , &phi     );   
  ch->SetBranchAddress( "phistarZ"     , &phi1    );   
  ch->SetBranchAddress( "ZGamMass"        , &mzg     );   
  */
  /*
  ch->SetBranchAddress( "Z1Mass"        , &m1      );   
  ch->SetBranchAddress( "Z2Mass"        , &m2      );   
  ch->SetBranchAddress( "helcosthetaZ1" , &h1      );   
  ch->SetBranchAddress( "helcosthetaZ2" , &h2      );   
  ch->SetBranchAddress( "costhetastar"  , &hs      );   
  ch->SetBranchAddress( "helphi"        , &phi     );   
  ch->SetBranchAddress( "phistarZ1"     , &phi1    );   
  ch->SetBranchAddress( "ZZMass"        , &mzg     );   
  */

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
    
    // 
    // initialise the differential cross-sections
    // 
    dXsec_ZGam_MCFM = 0.;
    dXsec_HZGam_MCFM = 0.;
    Discriminant = 0;
    PreBoostMass = 0;
    PostBoostMass = 0;
    CosT_lp = 0; CosT_lm = 0; Phi_lp = 0; CosTZG = 0; PtG = 0; Eta_lp = 0; Eta_lm = 0; Eta_g = 0;

    ch->GetEntry(ievt);

    // set four momenta
    //vector<TLorentzVector> p;
    //p=Calculate4Momentum(mzg,m1,m2,acos(hs),acos(h1),acos(h2),phi1,phi);

    TLorentzVector p0, p1, p2, psum;
    TVector3 bv;
    //When Pulled from Components
    p0.SetPxPyPzE(l_minus_px, l_minus_py, l_minus_pz, l_minus_e);
    p1.SetPxPyPzE(l_plus_px, l_plus_py, l_plus_pz, l_plus_e);
    p2.SetPxPyPzE(gamma_px, gamma_py, gamma_pz, gamma_e);

    //When Pulled from TLorentzVector
    //p0.SetPxPyPzE(l_minus_event->Px(), l_minus_event->Py(), l_minus_event->Pz(), l_minus_event->Energy());
    //p1.SetPxPyPzE(l_plus_event->Px(), l_plus_event->Py(), l_plus_event->Pz(), l_plus_event->Energy());
    //p2.SetPxPyPzE(gamma_event->Px(), gamma_event->Py(), gamma_event->Pz(), gamma_event->Energy());
    
    psum = p0 + p1 + p2;
    PreBoostMass = psum.M();
    bv = -psum.BoostVector();
    p0.Boost(bv);
    p1.Boost(bv);
    p2.Boost(bv);

    //Calculating ZGAngles
    ZGLabVectors genLevelInputs;
    ZGAngles genLevelOutputs;

    genLevelInputs.veczg = p0 + p1 + p2;
    genLevelInputs.vecz = p0 + p1;
    genLevelInputs.vecg = p2;
    genLevelInputs.veclm = p0;
    genLevelInputs.veclp = p1;
    
    getZGAngles(genLevelInputs,genLevelOutputs, false);

    CosT_lp = genLevelOutputs.costheta_lp;
    CosT_lm = genLevelOutputs.costheta_lm;
    Phi_lp = genLevelOutputs.phi;
    CosTZG = genLevelOutputs.cosTheta;
    PtG = genLevelOutputs.ptg;
    Eta_lp = genLevelOutputs.etal1;
    Eta_lm = genLevelOutputs.etal2;
    Eta_g = genLevelOutputs.etag;

    //Making cuts from Campbell/Ellis/Giele/Williams Paper
    if ((genLevelOutputs.ptg < 15.0) ||
	(std::abs(genLevelOutputs.etal1) > 2.0) ||
	(std::abs(genLevelOutputs.etal2) > 2.0) ||
	(std::abs(genLevelOutputs.etag) > 2.0) ||
	(genLevelOutputs.ptl1 < 20.0) ||
	(genLevelOutputs.ptl2 < 20.0) ||
	(genLevelOutputs.mz < 60.0) ||
	(genLevelOutputs.mz > 120.0) ||
	(genLevelOutputs.mzg < 115.0) ||
	(genLevelOutputs.mzg > 135.0) ||
	(std::abs(p0.DeltaR(p2)) < 0.7) ||
	(std::abs(p1.DeltaR(p2)) < 0.7)){
      if (verbosity >= TVar::DEBUG)
	std::cout << "Event does not pass cuts" << endl;
      continue;
    }

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
    
    // ==== Begin the differential cross-section calculation
    Xcal2.SetHiggsMass(zgammass);
    Xcal2.SetMatrixElement(TVar::MCFM);
    //dXsec_ZZ_MCFM = Xcal2.XsecCalc(TVar::ZZ_2e2m, TVar::GG, hzz4l_event,verbosity);
    //dXsec_HZZ_MCFM = Xcal2.XsecCalc(TVar::HZZ_4l, TVar::GG, hzz4l_event,verbosity);
    dXsec_ZGam_MCFM = Xcal2.XsecCalc(TVar::qqb_zgam, TVar::QQB, hzgamma_event,verbose);
    dXsec_HZGam_MCFM = Xcal2.XsecCalc(TVar::gg_hzgam, TVar::GG, hzgamma_event,verbose);
    logBkg = -log10(dXsec_ZGam_MCFM);
    logSig = -log10(dXsec_HZGam_MCFM);
    Discriminant = -log(dXsec_ZGam_MCFM/(dXsec_ZGam_MCFM+dXsec_HZGam_MCFM));
    if (dXsec_HZGam_MCFM > 10000000000)
      {
	std::cout << "HUGE XSEC" << endl;
	continue;
      }

    evt_tree->Fill();
    
  }//nevent
  
  if (verbosity >= TVar::INFO) cout << "TotCalculatedEvts: " << evt_tree->GetEntries() <<endl; 
  
  newfile->cd(); 
  evt_tree->Write(); 
  newfile->Close();
  
}  
