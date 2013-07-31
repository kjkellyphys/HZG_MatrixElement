#include "ElectronFunctions.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>

#define VARLIST_KINEMATICS "l1type/I:l2type/I:l3type/I:l4type/I:l1pt/F:l2pt/F:l3pt/F:l4pt/F:l1eta/F:l2eta/F:l3eta/F:l4eta/F:l1phi/F:l2phi/F:l3phi/F:l4phi/F:Z1pt/F:Z2pt/F:ZZpt/F:ZZdotZ1/F:ZZdotZ2/F:ZZptCosDphiZ1pt/F:ZZptCosDphiZ2pt/F:Z1eta/F:Z2eta/F:ZZy/F:mZ1/F:mZ2/F:m4l/F:channel/I:mZ1err/F:mZ2err/F:m4lerr/F"

typedef struct { 
  int l1type,l2type,l3type,l4type;
  float l1pt,l2pt,l3pt,l4pt;
  float l1eta,l2eta,l3eta,l4eta;
  float l1phi,l2phi,l3phi,l4phi;
  float Z1pt, Z2pt, ZZpt;
  float ZZdotZ1, ZZdotZ2, ZZptZ1ptCosDphi, ZZptZ2ptCosDphi; 
  float Z1eta, Z2eta, ZZy;
  float mZ1, mZ2, m4l;
  int channel;
  float mZ1err, mZ2err, m4lerr;
} KinematicsStruct;


int main() { 

  TFile * f = new TFile("f11.angles_b.root");
  TTree * t = (TTree*)(f->Get("zznt"));

  initializeEfficiencyWeights( "sf2011.root", "sf2012.root");

  KinematicsStruct k;
  t->SetBranchAddress("kinematics", &k);
  TLorentzVector * v = new TLorentzVector();

  for( int i=0; i<t->GetEntries()&&i<10; i++ ) {
    t->GetEntry(i);
    v->SetPtEtaPhiM(k.l1pt,k.l1eta,k.l1phi,0.);

    float wgt        = getEfficiencyWeight( v, CorrectionType::CENTRAL, 2011);
    float wgt_hi     = getEfficiencyWeight( v, CorrectionType::UPPER, 2011);
    float wgt_lo     = getEfficiencyWeight( v, CorrectionType::LOWER, 2011);
    float ptcorr     = getCorrectedPt( v, CorrectionType::CENTRAL, 2011);
    float ptcorr_hi  = getCorrectedPt( v, CorrectionType::UPPER, 2011);
    float ptcorr_lo  = getCorrectedPt( v, CorrectionType::LOWER, 2011);

    std::cout << "pt: " << v->Pt() << "\t"
	      << "eta: " << v->Eta() << "\t"
	      << "wgt: " << wgt << " (" << wgt_lo << " - " << wgt_hi << ")\t"
	      << "ptcorr: " << ptcorr << " (" << ptcorr_lo << " - " << ptcorr_hi << ")\t"
	      << std::endl;
  } 

}
