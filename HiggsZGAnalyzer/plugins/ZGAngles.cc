// ======================================================================================================
//
// angles following conventions in Gainer, et. al.  
// arXiv:1112.1405v2 [hep-ph], 6 May 2012
// see also : Gainer, et. al. 
// arVix:1108.2274v1 [hep-ph], 10 Aug 2011
//
// ======================================================================================================

#include "ZGAngles.h"
#include <iostream>
using namespace std;

#ifndef STANDALONE 
#include "CommonDefs.h"

// --------------------------------------------------------------------------------------------------
// 
// fill angle struct from ZZ-type GenInfo objects
//
void 
fillZGAngles( GenInfoStruct &geninfo, ZGAngles &a, bool debug )  
// --------------------------------------------------------------------------------------------------
{
  ZGLabVectors l;

  l.veckq = geninfo.x_1*TLorentzVector( 0,0,3500, 3500);
  l.veckqbar = geninfo.x_2*TLorentzVector( 0,0,-3500, 3500);

  if( geninfo.id_1_a < 0 ) {
    l.veclp.SetPtEtaPhiM( geninfo.pt_1_a, geninfo.eta_1_a, geninfo.phi_1_a, 
        (abs(geninfo.id_1_a) == 11 ? ELECTRON_MASS : MUON_MASS)  );
    l.veclm.SetPtEtaPhiM( geninfo.pt_2_a, geninfo.eta_2_a, geninfo.phi_2_a, 
        (abs(geninfo.id_2_a) == 11 ? ELECTRON_MASS : MUON_MASS)  );
  } else { 
    l.veclp.SetPtEtaPhiM( geninfo.pt_2_a, geninfo.eta_2_a, geninfo.phi_2_a, 
        (abs(geninfo.id_2_a) == 11 ? ELECTRON_MASS : MUON_MASS)  );
    l.veclm.SetPtEtaPhiM( geninfo.pt_1_a, geninfo.eta_1_a, geninfo.phi_1_a, 
        (abs(geninfo.id_1_a) == 11 ? ELECTRON_MASS : MUON_MASS)  );
 }

  
  //  l.vecz = (l.veclp+l.veclm);
  l.vecg.SetPtEtaPhiM( geninfo.vpt_b, geninfo.veta_b, geninfo.vphi_b, geninfo.vmass_b);
  l.vecz.SetPtEtaPhiM( geninfo.vpt_a, geninfo.veta_a, geninfo.vphi_a, geninfo.vmass_a);
  l.veczg = (l.vecz+l.vecg);
  
  getZGAngles( l, a, debug );  
}

// --------------------------------------------------------------------------------------------------
//
// fill angle struct from ZZ-type EventData
//
void 
fillZGAngles( EventData &data, ZGAngles &a ) 
// --------------------------------------------------------------------------------------------------
{ 
  ZGLabVectors l;
  
  if( data.Z1leptons[0].charge > 0 ) {
    l.veclp = data.Z1leptons[0].vec;  
    l.veclm = data.Z1leptons[1].vec;  
    l.vecg  = data.Z2leptons[0].vec;  
  } else { 
    l.veclp = data.Z1leptons[1].vec;  
    l.veclm = data.Z1leptons[0].vec;  
    l.vecg = data.Z2leptons[0].vec;  
 }

  l.vecz = (l.veclp+l.veclm);
  l.veczg = (l.vecz+l.vecg);
  
  getZGAngles( l, a );
}

#endif



// --------------------------------------------------------------------------------------------------
//
// set angles from a struct of lab-frame 4vecs 
//
void 
getZGAngles( ZGLabVectors &l, ZGAngles &a, bool debug ) 
// --------------------------------------------------------------------------------------------------
{
  if( debug ) { 
    cout << endl << endl;
    cout << "veczg(x:y:z) :: " << l.veczg.X() <<":"<< l.veczg.Y() <<":"<< l.veczg.Z() << endl;
    cout << "\tpt: " << l.veczg.Pt() << "\teta: " << l.veczg.Eta() << "\tphi:" << l.veczg.Phi() << endl;
    cout << "vecz(x:y:z) :: "  << l.vecz.X() <<":"<< l.vecz.Y() <<":"<< l.vecz.Z() << endl;
    cout << "\tpt: " << l.vecz.Pt() << "\teta: " << l.vecz.Eta() << "\tphi:" << l.vecz.Phi() << endl;
    cout << "veclp(x:y:z) :: "  << l.veclp.X() <<":"<< l.veclp.Y() <<":"<< l.veclp.Z() << endl;
    cout << "\tpt: " << l.veclp.Pt() << "\teta: " << l.veclp.Eta() << "\tphi:" << l.veclp.Phi() << endl;
    cout << "veclm(x:y:z) :: "  << l.veclm.X() <<":"<< l.veclm.Y() <<":"<< l.veclm.Z() << endl;
    cout << "\tpt: " << l.veclm.Pt() << "\teta: " << l.veclm.Eta() << "\tphi:" << l.veclm.Phi() << endl;
    cout << "vecg(x:y:z) :: "  << l.vecg.X() <<":"<< l.vecg.Y() <<":"<< l.vecg.Z() << endl;
    cout << "\tpt: " << l.vecg.Pt() << "\teta: " << l.vecg.Eta() << "\tphi:" << l.vecg.Phi() << endl;
  }



  //
  // define the frames & coordinate systems
  // --------------------------------------------------------------------------------------------------

  TVector3 Xframe  = l.veczg.BoostVector();
  TVector3 Z1frame = l.vecz.BoostVector();

  if( debug ) { 
    cout << "Xboost(x:y:z) :: " << Xframe.X() <<":"<< Xframe.Y() <<":"<< Xframe.Z() << endl;
    cout << "Zboost(x:y:z) :: " << Z1frame.X() <<":"<< Z1frame.Y() <<":"<< Z1frame.Z() << endl;
  }

  // "partons" 
  TLorentzVector  kq( 0, 0, (l.veczg.E()+l.veczg.Pz())/2, (l.veczg.E()+l.veczg.Pz())/2 );
  TLorentzVector  kqbar( 0, 0, (l.veczg.Pz()-l.veczg.E())/2, (l.veczg.E()-l.veczg.Pz())/2 );
  TLorentzVector veckq_in_Xframe(kq);
  TLorentzVector veckqbar_in_Xframe(kqbar);
  veckq_in_Xframe.Boost(-1*Xframe);
  veckqbar_in_Xframe.Boost(-1*Xframe);
  if( debug ) { 
    cout << "kq_X(x:y:z) :: " << veckq_in_Xframe.X() <<":"<< veckq_in_Xframe.Y() <<":"<< veckq_in_Xframe.Z() << endl;
    cout << "kqb_X(x:y:z) :: " << veckqbar_in_Xframe.X() <<":"<< veckqbar_in_Xframe.Y() <<":"<< veckqbar_in_Xframe.Z() << endl;
  }

  // Z vectors
  if( debug ) { cout << "boosting Z,g ..." << endl;}
  TLorentzVector vecz_in_Xframe (l.vecz); 
  TLorentzVector vecg_in_Xframe (l.vecg); 
  TLorentzVector vecz_in_Z1frame (l.vecz); 
  vecz_in_Xframe.Boost(-1*Xframe);
  vecg_in_Xframe.Boost(-1*Xframe);
  vecz_in_Z1frame.Boost(-1*Z1frame);

  if( debug ) { cout << "rotating ..." << endl;}
  // coord system in the CM frame
  TVector3 uz_in_Xframe = vecz_in_Xframe.Vect().Unit();
  TVector3 uy_in_Xframe = (veckq_in_Xframe.Vect().Unit().Cross(uz_in_Xframe.Unit())).Unit();
  TVector3 ux_in_Xframe = (uy_in_Xframe.Unit().Cross(uz_in_Xframe.Unit())).Unit();
  TRotation rotation;
  rotation = rotation.RotateAxes( ux_in_Xframe,uy_in_Xframe,uz_in_Xframe ).Inverse();

  //
  // for going to the Z frames from the CM frame,
  // boost after transform
  //
  if( debug ) { cout << "xform z ..." << endl; }
  TLorentzVector vecz_in_Xframe_newcoords(vecz_in_Xframe);
  vecz_in_Xframe_newcoords.Transform(rotation);
  TVector3 Z1frame_from_Xframe_newcoords = vecz_in_Xframe_newcoords.BoostVector();
  if( debug ) { 
    cout << "vecz_in_Xframe_newcoords (x:y:z) :: " 
   << vecz_in_Xframe_newcoords.X() <<":"
   << vecz_in_Xframe_newcoords.Y() <<":"
   << vecz_in_Xframe_newcoords.Z() << endl;
    
    cout << "Z1frame_from_Xframe_newcoords (x:y:z) :: " 
   << Z1frame_from_Xframe_newcoords.X() <<":"
   << Z1frame_from_Xframe_newcoords.Y() <<":"
   << Z1frame_from_Xframe_newcoords.Z() << endl;
  }


  //  
  // theta(lm), phi(lm) in Z1 frame
  // --------------------------------------------------------------------------------------------------
  TLorentzVector veclm_in_Z1frame(l.veclm);
  TLorentzVector veclp_in_Z1frame(l.veclp);

  // first boost to CM, then redefine coords  
  if( debug ) { cout << "xform lm..." << endl; }
  veclm_in_Z1frame.Boost(-1*Xframe);
  veclm_in_Z1frame.Transform(rotation);
  veclp_in_Z1frame.Boost(-1*Xframe);
  veclp_in_Z1frame.Transform(rotation);
  
  // then boost to Z1
  if ( debug ) { cout << "boost lm ..." << endl;}
  veclm_in_Z1frame.Boost(-1*Z1frame_from_Xframe_newcoords);
  veclp_in_Z1frame.Boost(-1*Z1frame_from_Xframe_newcoords);

  // now get angles
  a.phi      = veclm_in_Z1frame.Phi();
  a.costheta_lm = veclm_in_Z1frame.CosTheta();
  a.costheta_lp = veclp_in_Z1frame.CosTheta();
  
  //
  // or, just go directly to Z frame, skip CM frame
  //
  /*
  // direct
  TLorentzVector veclm_in_Z1frame_direct(l.veclm);
  veclm_in_Z1frame_direct.Boost(-1*Z1frame);
  TLorentzVector veclp_in_Z1frame_direct(l.veclp);
  veclp_in_Z1frame_direct.Boost(-1*Z1frame);
  a.costheta_lm = veclm_in_Z1frame_direct.Vect().Unit().Dot((-1*l.vecz.Vect()).Unit());
  a.costheta_lp = veclp_in_Z1frame_direct.Vect().Unit().Dot((-1*l.vecz.Vect()).Unit());
  */
    
  if( debug ) { 
    TVector3 unitz(0,0,1);
    cout << "compare :: ct: " << a.costheta_lm
   << "\tv.ct: " << veclm_in_Z1frame.Vect().CosTheta()
   << "\tc(t): " << TMath::Cos(veclm_in_Z1frame.Theta())
   << "\tc(v.t): " << TMath::Cos(veclm_in_Z1frame.Theta())
   << "\tdot: " << veclm_in_Z1frame.Vect().Unit().Dot(unitz)
   << endl;
    cout << "costheta lp vs lm : " << veclp_in_Z1frame.CosTheta() << " vs " << veclm_in_Z1frame.CosTheta() << endl;
  }


  //  
  // Theta(Z1,Zg-lab) in X frame 
  // --------------------------------------------------------------------------------------------------
  TLorentzVector veczg_in_Xframe(l.veczg);
  //  veczg_in_Xframe.Boost(-1*Xframe);
  veczg_in_Xframe.Transform(rotation);

  if( debug ) { 
    cout << "veczg_X (x:y:z) ::  " << veczg_in_Xframe.X() 
   <<":"<< veczg_in_Xframe.Y()  
   <<":"<< veczg_in_Xframe.Z()  
   << endl;
    cout << "vecz_X (x:y:z) ::  " << vecz_in_Xframe.X() 
   <<":"<< vecz_in_Xframe.Y()  
   <<":"<< vecz_in_Xframe.Z()  
   << endl;
    cout << "vecg_X (x:y:z) ::  " << vecg_in_Xframe.X() 
   <<":"<< vecg_in_Xframe.Y()  
   <<":"<< vecg_in_Xframe.Z()  
   << endl;

    // check
    TLorentzVector ztmp(l.vecz);
    ztmp.Boost(-1*Xframe);
    ztmp.Transform(rotation);
    
    TLorentzVector gtmp(l.vecg);
    gtmp.Boost(-1*Xframe);
    gtmp.Transform(rotation);
    
    cout << "ztmp.CosTheta(): " << ztmp.CosTheta() << endl;
    cout << "gtmp.CosTheta(): " << gtmp.CosTheta() << endl;
  }
  

  TLorentzVector veczg_in_Xframe_newcoords(l.veczg);
  //  veczg_in_Xframe_newcoords.Boost(-1*Xframe);
  veczg_in_Xframe_newcoords.Transform(rotation);
  a.cosTheta = (-1*veczg_in_Xframe_newcoords.Vect()).CosTheta();

  if( debug ) { cout << "phi: " << a.phi << "\tct: " << a.costheta_lm << "\tcT:" << a.cosTheta << endl;}

  a.ptg   = l.vecg.Pt();
  a.etag  = l.vecg.Eta();
  a.ptl1  = l.veclp.Pt();
  a.etal1 = l.veclp.Eta();
  a.ptl2  = l.veclm.Pt();
  a.etal2 = l.veclm.Eta();
  a.ptz   = (l.veclp+l.veclm).Pt();
  a.etaz  = (l.veclp+l.veclm).Eta();
  a.mzg   = l.veczg.M();
  a.mz    = l.vecz.M();
};
