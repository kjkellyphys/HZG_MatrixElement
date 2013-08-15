//-----------------------------------------------------------------------------
//
// Class EventProb Module
//
//   EventProb Module
//
// March 21 2011
// S. Jindariani (sergo@fnal.gov)
// Y. Gao (ygao@fnal.gov)
// K. Burkett (burkett@fnal.gov)
//-----------------------------------------------------------------------------

#include "TEvtProb.hh"
#include "TVar.hh"


ClassImp(TEvtProb)

    using namespace std;

    //-----------------------------------------------------------------------------
    // Constructors and Destructor
    //-----------------------------------------------------------------------------
TEvtProb::TEvtProb() {
  mcfm_init_("input.DAT", "./");
  SetEwkCoupligParameters();
  coupling_();
  myCSW_ = new HiggsCSandWidth("../txtFiles");
}


TEvtProb::~TEvtProb() {
    delete myCSW_;
}

//
// Directly calculate the ZZ->4l differential cross-section 
// WARNING: in development
// 
double TEvtProb::XsecCalc(TVar::Process proc, TVar::Production production, const hzgamma_event_type &hzgamma_event,
			  bool verbose){

    //Initialize Process

    SetProduction(production);

    if ( _matrixElement == TVar::MCFM) 
      My_choose(proc);
    
    //constants
    double sqrts = 2.*EBEAM;
    double W=sqrts*sqrts;
    
    //Weight calculation
    double flux=1.;
    double dXsec=0.;
    
    mcfm_event_type mcfm_event; 
    // assign the right initial momentum
    // assumes the events are boosted to have 0 transverse momenta
    double sysPz= ( hzgamma_event.p[0] + hzgamma_event.p[1] + hzgamma_event.p[2]).Pz();
    double sysE = ( hzgamma_event.p[0] + hzgamma_event.p[1] + hzgamma_event.p[2]).Energy();
    double pz0 = (sysE+sysPz)/2.; 
    double pz1 = -(sysE-sysPz)/2.;
    mcfm_event.p[0].SetPxPyPzE   (0., 0., pz0, TMath::Abs(pz0));
    mcfm_event.p[1].SetPxPyPzE   (0., 0., pz1, TMath::Abs(pz1));
    mcfm_event.p[2].SetPxPyPzE   (hzgamma_event.p[0].Px(), hzgamma_event.p[0].Py(), hzgamma_event.p[0].Pz(), hzgamma_event.p[0].Energy());
    mcfm_event.p[3].SetPxPyPzE   (hzgamma_event.p[1].Px(), hzgamma_event.p[1].Py(), hzgamma_event.p[1].Pz(), hzgamma_event.p[1].Energy());
    mcfm_event.p[4].SetPxPyPzE   (hzgamma_event.p[2].Px(), hzgamma_event.p[2].Py(), hzgamma_event.p[2].Pz(), hzgamma_event.p[2].Energy());
    //mcfm_event.p[5].SetPxPyPzE   (hzgamma_event.p[3].Px(), hzgamma_event.p[3].Py(), hzgamma_event.p[3].Pz(), hzgamma_event.p[3].Energy());

    if (production == TVar::QQB){
      mcfm_event.PdgCode[0] = 2;
      mcfm_event.PdgCode[1] = -2;
    }
    else if (production == TVar::GG){
      mcfm_event.PdgCode[0] = 21;
      mcfm_event.PdgCode[1] = 21;
    }
    else{
      mcfm_event.PdgCode[0] = 21;
      mcfm_event.PdgCode[1] = 21;
    }
    mcfm_event.PdgCode[2] = hzgamma_event.PdgCode[0];
    mcfm_event.PdgCode[3] = hzgamma_event.PdgCode[1];
    mcfm_event.PdgCode[4] = hzgamma_event.PdgCode[2];
    //mcfm_event.PdgCode[5] = hzgamma_event.PdgCode[3];

    /*
    //for ( int i = 0; i < 6; i++ ) {
    for (int i = 0; i < 5; i++) {  
    std::cout << "Particle " << i << " (Px, Py, Pz, E): " <<  mcfm_event.p[i].Px() << ", " << mcfm_event.p[i].Py() 
		<< ", " << mcfm_event.p[i].Pz() << ", " << mcfm_event.p[i].Energy() <<  "\n";
    }
    */
    //Matrix Element evaluation in qX=qY=0 frame
    //Evaluate f(x1)f(x2)|M(q)|/x1/x2 
    // 
    double qX=mcfm_event.p[0].Px()+mcfm_event.p[1].Px();
    double qY=mcfm_event.p[0].Py()+mcfm_event.p[1].Py();
    
    if((qX*qX+qY*qY)>0){
      double qE = mcfm_event.p[0].Energy()+mcfm_event.p[1].Energy();
      TVector3 boostV(qX/qE,qY/qE,0);
      for(int ipt=0;ipt<6;ipt++) mcfm_event.p[ipt].Boost(-boostV);
    }
    //event selections in Lab Frame
    double flavor_msq[nmsq][nmsq];
    double msqjk=0; 
    if ( _matrixElement == TVar::MCFM ) 
      msqjk = SumMatrixElementPDF(proc, &mcfm_event, flavor_msq, &flux);
    
    if(msqjk<=0){ mcfm_event.pswt=0; }
    
    flux=fbGeV2/(mcfm_event.p[0].Energy()*mcfm_event.p[1].Energy())	/(4*W);
    //    dXsec=msqjk*flux;
    dXsec=msqjk;
    
    if (verbose)
      {
	cout <<"Process " << TVar::ProcessName(proc) << 
	  " TEvtProb::XsecCalc(): dXsec=" << dXsec
	     <<" Msq="<<msqjk 
	     <<" flux="<<flux 
	     <<endl;
      }

    return dXsec;

}

// this appears to be some kind of 
// way of setting MCFM parameters through
// an interface defined in TMCFM.hh
void TEvtProb::SetHiggsMass(double mass){
    masses_mcfm_.hmass=mass;
    masses_mcfm_.hwidth=myCSW_->HiggsWidth(0, min(mass,1000.) );
    _hmass = mass;
    _hwidth = myCSW_->HiggsWidth(0, min(mass,1000.) );
    /*
    //
    // get higgs width for 125 and 250 GeV
    // 
    std::cout << "H125 width " << myCSW_->HiggsWidth(0, 125);
    std::cout << "H250 width " << myCSW_->HiggsWidth(0, 250);
    */
}
