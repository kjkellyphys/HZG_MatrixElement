#define simple_v2_cxx
// The class definition in simple_v2.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("simple_v2.C")
// Root > T->Process("simple_v2.C","some options")
// Root > T->Process("simple_v2.C+")
//

#include "simple_v2.h"
#include <TH2.h>
#include <TStyle.h>

#include "TVar.hh"
#include "TMCFM.hh"
#include "TUtil.cc"
#include "TEvtProb.cc"

string  selection      = "eeGamma";
string  period         = "2012";
string  abcd           = "ABCD";
string  suffix         = "SUFFIX";
bool verbose           = false;
int quickCount         = 0;

bool MuonSortCondition(const TCMuon& m1, const TCMuon& m2) {return (m1.Pt() > m2.Pt());}

void simple_v2::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  histoFile = new TFile("simpleHistograms_SUFFIX.root", "RECREATE");
  histoFile->cd();
  hm        = new HistManager(histoFile);
  TString option = GetOption();

  genHZG = {};

  //KK
  /*
  newfile = new TFile("simpleHistograms_KK.root", "RECREATE");
  newfile->cd();
  ch = new TTree("K", "4Vectors of Events");

  p_lminus = new TLorentzVector();
  p_lplus = new TLorentzVector();
  p_gamma = new TLorentzVector();

  float dXsec_ZGam_MCFM = 0.;
  float dXsec_HZGam_MCFM = 0.;
  float Discriminant = 0.;
  float PreBoostMass = 0.;
  float PostBoostMass = 0.;
  float logBkg(0.), logSig(0.);

  ch->Branch("dxSec_ZGam_MCFM"   , &dXsec_ZGam_MCFM ,"dXsec_ZGam_MCFM/F");
  ch->Branch("dxSec_HZGam_MCFM"   , &dXsec_HZGam_MCFM ,"dXsec_HZGam_MCFM/F");
  ch->Branch("logBkg", &logBkg, "logBkg/F");
  ch->Branch("logSig", &logSig, "logSig/F");
  ch->Branch("Discriminant", &Discriminant, "Discriminant/F");
  ch->Branch("PreBoostMass", &PreBoostMass, "PreBoostMass/F");
  ch->Branch("PostBoostMass", &PostBoostMass, "PostBoostMass/F");
  */
}

void simple_v2::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

}

Bool_t simple_v2::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either simple_v2::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  GetEntry(entry);

  ////////////////////
  // Find Gen Stuff //
  ////////////////////

  //KK
  /*
  float dXsec_ZGam_MCFM = 0.;
  float dXsec_HZGam_MCFM = 0.;
  float Discriminant = 0;
  float PreBoostMass = 0;
  float PostBoostMass = 0;
  float logBkg = 0;
  float logSig = 0;

  TEvtProb Xcal2;
  hzgamma_event_type hzgamma_event;
  */
  vector<TCGenParticle> vetoPhotons;
  CleanUpGen(genHZG);
  genHZG = {};
  if(!isRealData){
    ///////// load all the relevent particles into a struct /////////
    FindGenParticles(genParticles, selection, vetoPhotons, genHZG);

    ///////// gen angles, plots before any kinematic/fiducial cleaning //////////////

    ZGLabVectors genLevelInputs;
    ZGAngles     genLevelOutputs;


    if(genHZG.lp && genHZG.lm && genHZG.g){
      //genHZG.lp->Print();
      //genHZG.lm->Print();
      //genHZG.g->Print();
      genLevelInputs.veczg = *genHZG.lp+*genHZG.lm+*genHZG.g;
      genLevelInputs.vecz = *genHZG.lp+*genHZG.lm; 
      genLevelInputs.vecg = *genHZG.g;

      genLevelInputs.veclp = *genHZG.lp;
      genLevelInputs.veclm = *genHZG.lm;

      if (verbose){
        cout<<"event: "<<eventNumber<<endl;
        cout<<"lp:"<<endl;
        genHZG.lp->Print();
        cout<<"lm:"<<endl;
        genHZG.lm->Print();
        cout<<"g:"<<endl;
        genHZG.g->Print();
      }


      if (fabs(genLevelInputs.veczg.M()-125.0)>2.0){
        cout<<"odd mass"<<endl;
        cout<<"higgsMass: "<<genLevelInputs.veczg.M()<<endl;
        cout<<"zMass (ll): "<<genLevelInputs.vecz.M()<<endl;
        cout<<"zMass (ll): "<<genLevelInputs.vecz.M()<<endl;
        cout<<endl;
        return kTRUE;
      }
      getZGAngles(genLevelInputs,genLevelOutputs, false);
      AnglePlots(genLevelOutputs,1);
      quickCount += 1;
      //cout<<"costheta_lm: "<<genLevelOutputs.costheta_lm<<"\tcostheta_lp: "<<genLevelOutputs.costheta_lp<<"\tphi: "<<genLevelOutputs.phi<<"\tcosTheta: "<<genLevelOutputs.cosTheta<<"\tcosThetaG: "<<genLevelOutputs.cosThetaG<<endl;

      //KK
      /*
      p_lminus->SetPxPyPzE(genHZG.lm->Px(), genHZG.lm->Py(), genHZG.lm->Pz(), genHZG.lm->Energy());
      p_lplus->SetPxPyPzE(genHZG.lp->Px(), genHZG.lp->Py(), genHZG.lp->Pz(), genHZG.lp->Energy());
      p_gamma->SetPxPyPzE(genHZG.g->Px(), genHZG.g->Py(), genHZG.g->Pz(), genHZG.g->Energy());

      TLorentzVector p0, p1, p2, psum;
      TVector3 bv;
      p0.SetPxPyPzE(genHZG.lm->Px(), genHZG.lm->Py(), genHZG.lm->Pz(), genHZG.lm->Energy());
      p1.SetPxPyPzE(genHZG.lp->Px(), genHZG.lp->Py(), genHZG.lp->Pz(), genHZG.lp->Energy());
      p2.SetPxPyPzE(genHZG.g->Px(), genHZG.g->Py(), genHZG.g->Pz(), genHZG.g->Energy());
      psum = p0 + p1 + p2;
      PreBoostMass = psum.M();
      bv = -psum.BoostVector();
      p0.Boost(bv);
      p1.Boost(bv);
      p2.Boost(bv);

      hzgamma_event.p[0].SetPxPyPzE(p0.Px(), p0.Py(), p0.Pz(), p0.Energy());
      hzgamma_event.p[1].SetPxPyPzE(p1.Px(), p1.Py(), p1.Pz(), p1.Energy());
      hzgamma_event.p[2].SetPxPyPzE(p2.Px(), p2.Py(), p2.Pz(), p2.Energy());

      hzgamma_event.PdgCode[0] = 11;
      hzgamma_event.PdgCode[1] = -11;
      hzgamma_event.PdgCode[2] = 22;

      float zmass = (hzgamma_event.p[0]+hzgamma_event.p[1]).M();
      float gammass = (hzgamma_event.p[2]).M();
      float zgammass = (hzgamma_event.p[0]+hzgamma_event.p[1]+hzgamma_event.p[2]).M();
      PostBoostMass = zgammass;

      if (verbose) {
	cout << "\n=========================================================\n";
	cout << "Entry: " << entry << "\n";
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
      Xcal2.SetHiggsMass(zgammass);
      Xcal2.SetMatrixElement(TVar::MCFM);
      dXsec_ZGam_MCFM = Xcal2.XsecCalc(TVar::qqb_zgam, TVar::QQB, hzgamma_event,verbose);
      dXsec_HZGam_MCFM = Xcal2.XsecCalc(TVar::gg_hzgam, TVar::GG, hzgamma_event,verbose);
      logBkg = -log10(dXsec_ZGam_MCFM);
      logSig = -log10(dXsec_HZGam_MCFM);
      Discriminant = -log(dXsec_ZGam_MCFM/(dXsec_ZGam_MCFM+dXsec_HZGam_MCFM));


      ch->SetBranchAddress("dxSec_ZGam_MCFM", &dXsec_ZGam_MCFM);
      ch->SetBranchAddress("dxSec_HZGam_MCFM", &dXsec_HZGam_MCFM);
      ch->SetBranchAddress("logBkg", &logBkg);
      ch->SetBranchAddress("logSig", &logSig);
      ch->SetBranchAddress("Discriminant", &Discriminant);
      ch->SetBranchAddress("PreBoostMass", &PreBoostMass);
      ch->SetBranchAddress("PostBoostMass", &PostBoostMass);

      ch->Fill();
      */
    }
  }
  return kTRUE;
}

void simple_v2::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

}

void simple_v2::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  cout<<"unskimmed total: "<<unskimmedEventsTotal<<endl;
  cout<<"quickCount: "<<quickCount<<endl;
  histoFile->Write();
  histoFile->Close();  
  //KK
  //newfile->Write();
  //newfile->Close();
}

bool simple_v2::PassMuonID(TCMuon *mu, muIDCuts cutLevel){

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

bool simple_v2::PassMuonIso(TCMuon *mu, muIsoCuts cutLevel){

  float combIso;

  combIso = (mu->IsoMap("pfChargedHadronPt_R04")
    + max(0.,(double)mu->IsoMap("pfNeutralHadronEt_R04") + mu->IsoMap("pfPhotonEt_R04") - 0.5*mu->IsoMap("pfPUPt_R04")));

  bool isoPass = false;
  if (combIso/mu->Pt() < cutLevel.relCombIso04) isoPass = true;
  return isoPass;
}

void simple_v2::AnglePlots(ZGAngles &zga,float eventWeight)
{
  hm->fill1DHist(zga.costheta_lp,"h1_costhetaLP_SUFFIX", "Cos(#theta) positive lepton;cos(#theta);N_{evts}", 50, -1.1, 1.1, eventWeight);     
  hm->fill1DHist(zga.costheta_lm,"h1_costhetaLM_SUFFIX", "Cos(#theta) negative lepton;cos(#theta);N_{evts}", 50, -1.1, 1.1, eventWeight);     
  hm->fill1DHist(zga.phi,"h1_phi_SUFFIX", "#phi positive lepton;#phi;N_{evts}", 50, -3.2, 3.2, eventWeight);     
  hm->fill1DHist(zga.cosTheta,"h1_costhetaZG_SUFFIX", "Cos(#Theta) ZG system;cos(#Theta);N_{evts}", 50, -1.1, 1.1, eventWeight);     
  hm->fill1DHist(zga.ptg,"h1_ptgamma_SUFFIX", "Pt of photon in ZG System;Pt;N_{evts}",50,0,200,eventWeight);
  hm->fill1DHist(zga.etal1,"h1_etaLP_SUFFIX", "Eta of positive lepton;#eta;N_{evts}",50,-5.0,5.0,eventWeight);
  hm->fill1DHist(zga.etal2,"h1_etaLM_SUFFIX", "Eta of negative lepton;#eta;N_{evts}",50,-5.0,5.0,eventWeight);
  hm->fill1DHist(zga.etag,"h1_etagamma_SUFFIX", "Eta of Gamma;#eta;N_{evts}",50,-5.0,5.0,eventWeight);
}

void  simple_v2::FindGenParticles(TClonesArray *genParticles, string selection, vector<TCGenParticle>& vetoPhotons, genHZGParticles& _genHZG){
  vector<TCGenParticle> genElectrons;
  vector<TCGenParticle> genMuons;
  vector<TCGenParticle> genZs;
  vector<TCGenParticle> genWs;
  vector<TCGenParticle> genHs;
  vector<TCGenParticle> genPhotons;
  vector<TCGenParticle> genLeptons;
  bool isMuMuGamma = false;
  bool isEEGamma = false;
  bool goodPhot = false;

  for (int i = 0; i < genParticles->GetSize(); ++i) {
    TCGenParticle* thisGen = (TCGenParticle*) genParticles->At(i);    
  //  cout<<thisGen->GetPDGId()<<endl;
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
  ///////// sort gen particles by PT ///////////

  sort(genElectrons.begin(), genElectrons.end(), P4SortCondition);
  sort(genMuons.begin(), genMuons.end(), P4SortCondition);
  sort(genZs.begin(), genZs.end(), P4SortCondition);
  sort(genWs.begin(), genWs.end(), P4SortCondition);
  sort(genPhotons.begin(), genPhotons.end(), P4SortCondition);
  sort(genHs.begin(), genHs.end(), P4SortCondition);

  vetoPhotons = genPhotons;

  if (isMuMuGamma && (selection == "mumuGamma")) genLeptons = genMuons;
  else if (isEEGamma && (selection == "eeGamma")) genLeptons = genElectrons;
  
  if (_genHZG.lp){
    cout<<"well here's your fucking problem"<<endl;
    _genHZG.lp->Print();
    cout<<endl;
  }
  
  bool posLep = false;
  bool negLep = false;
  vector<TCGenParticle>::iterator testIt;

  if (genLeptons.size() > 1){
    for (testIt=genLeptons.begin(); testIt<genLeptons.end(); testIt++){
      if(testIt->Mother() == 23 && testIt->Charge() == 1 ){
        _genHZG.lp = new TCGenParticle(*testIt);
        posLep = true;
      }else if(testIt->Mother()== 23 && testIt->Charge() == -1){
        _genHZG.lm = new TCGenParticle((*testIt));
        negLep = true;
      }
      if (posLep && negLep) break;
    }
  }else { return;}

  if (genPhotons.size() > 0 && posLep && negLep){
      for (testIt=genPhotons.begin(); testIt<genPhotons.end(); testIt++){
        //cout<<"mother: "<<testIt->Mother()<<"\tstatus: "<<testIt->GetStatus()<<endl;
        if (testIt->Mother() == 25 && fabs((*testIt+*_genHZG.lm+*_genHZG.lp).M()-125.0) < 0.1) _genHZG.g = new TCGenParticle(*testIt); goodPhot = true; break;
      }
      if (!goodPhot) return;
    //_genHZG.g = new TCGenParticle(genPhotons.front());
  }else{ return;}


  if (genZs.size() > 0) _genHZG.z = new TCGenParticle(genZs.front());
  if (genWs.size() > 0) _genHZG.w = new TCGenParticle(genWs.front());
  if (genHs.size() > 0) _genHZG.h = new TCGenParticle(genHs.front());

  return;
}

void simple_v2::CleanUpGen(genHZGParticles& _genHZG){
  if (_genHZG.lp) delete _genHZG.lp;
  if (_genHZG.lm) delete _genHZG.lm;
  if (_genHZG.g) delete _genHZG.g;
  if (_genHZG.w) delete _genHZG.w;
  if (_genHZG.z) delete _genHZG.z;
  if (_genHZG.h) delete _genHZG.h;
}
