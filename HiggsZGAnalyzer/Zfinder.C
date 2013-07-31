#define Zfinder_cxx

#include "Zfinder.h"
#include <TH2.h>
#include <TStyle.h>
//
//PARAMETERS

string  selection      = "mumuGamma";
string  period         = "2012";
int     JC_LVL         = 0;
int     trigger[]      = {3,4};
string  suffix         = "SUFFIX";

using namespace std;


void Zfinder::Begin(TTree * /*tree*/)
{
  histoFile = new TFile("higgsHistograms.root", "RECREATE");
  histoFile->cd();
  hm        = new HistManager(histoFile);
  TString option = GetOption();

}

void Zfinder::SlaveBegin(TTree * /*tree*/)
{
  TString option = GetOption();

}

bool MuonSortCondition(const TCMuon& m1, const TCMuon& m2) {return (m1.Pt() > m2.Pt());}

Bool_t Zfinder::Process(Long64_t entry)
{

  GetEntry(entry);

  //////////////////
  // PV selection //
  //////////////////

  vector<TVector3> goodVertices;
  vector<TCPrimaryVtx> pvVect;

  if (primaryVtx->GetSize() < 1) return kTRUE;
  for (int i = 0; i < primaryVtx->GetSize(); ++i) {
    TCPrimaryVtx* pVtx = (TCPrimaryVtx*) primaryVtx->At(i);
    if (
        !pVtx->IsFake() 
        && pVtx->NDof() > 4.
        && fabs(pVtx->Position().z()) <= 24.
        && fabs(pVtx->Position().Perp()) <= 2.
       )

      goodVertices.push_back(pVtx->Position());
  }

  if (goodVertices.size() < 1) return kTRUE;
  pvPosition = new TVector3();
  *pvPosition = goodVertices[0];

  ////////////////////
  // Muon Selection //
  ////////////////////

  vector<TCMuon> muonsIDIso;
  for (int i = 0; i < recoMuons->GetSize(); ++ i)
  {
    TCMuon* thisMuon = (TCMuon*) recoMuons->At(i);    
    if (PassMuonID(thisMuon, tightMuID) && PassMuonIso(thisMuon, tightMuIso)) muonsIDIso.push_back(*thisMuon);

  }
  sort(muonsIDIso.begin(), muonsIDIso.end(), MuonSortCondition);
  /*
     for (int i = 0; i < genParticles->GetSize(); ++i) {
     TCGenParticle* thisGen = (TCGenParticle*) genParticles->At(i);    
     if (abs(thisGen->GetPDGId()) == 13){
     cout<<"hey stoyan i got some gen muon pt for you WOOOOOOOO: "<<thisGen->Pt()<<endl;

     hm->fill1DHist(thisGen->Pt(),"muonPt", "#mu_{pT};pT;N_{evts}", 100, 0., 100., weight);
     }

     }
     */

  if (selection =="mumuGamma") {
    if (muonsIDIso.size() < 2) return kTRUE;
    if (muonsIDIso[0].Charge() == muonsIDIso[1].Charge()) return kTRUE;
    if ((muonsIDIso[0].Pt() < 18) || (muonsIDIso[1].Pt() < 9)) return kTRUE;

    TLorentzVector ZP4          = muonsIDIso[0] + muonsIDIso[1];

    if (ZP4.M() > 50) hm->fill1DHist(ZP4.M(),"h1_diLeptonMassPreSel_SUFFIX", "M_{ll}; M_{ll};N_{evts}", 40, 70., 110.,1); 

  }

  return kTRUE;
}

void Zfinder::SlaveTerminate()
{

}

void Zfinder::Terminate()
{
  cout<<unskimmedEventsTotal<<endl;
  histoFile->Write();
  histoFile->Close();  
}

bool Zfinder::PassMuonID(TCMuon *mu, muIDCuts cutLevel){

  bool muPass = false;

  mu->LoadMap(*mapCuts);
  map<string, float> muCut = mu->MuonMap();

  if (suffix.find("2011") != string::npos){
    if (
        fabs(mu->Eta()) < 2.4
        && muCut["IsGLB"]                         == cutLevel.IsGLB
        && muCut["NormalizedChi2"]                < cutLevel.NormalizedChi2
        && muCut["NumberOfValidMuonHits"]         > cutLevel.NumberOfValidMuonHits
        && muCut["NumberOfMatchedStations"]       > cutLevel.NumberOfMatchedStations
        && muCut["NumberOfValidPixelHits"]        > cutLevel.NumberOfValidPixelHits
        && muCut["TrackLayersWithMeasurement"]    > cutLevel.TrackLayersWithMeasurement
        && fabs(mu->Dxy(pvPosition))              < cutLevel.dxy
        && fabs(mu->Dz(pvPosition))               < cutLevel.dz
       ) muPass = true;
  }else{
    if (
        fabs(mu->Eta()) < 2.4
        && muCut["IsPF"]                          == cutLevel.IsPF
        && muCut["IsGLB"]                         == cutLevel.IsGLB
        && muCut["NormalizedChi2"]                < cutLevel.NormalizedChi2
        && muCut["NumberOfValidMuonHits"]         > cutLevel.NumberOfValidMuonHits
        && muCut["NumberOfMatchedStations"]       > cutLevel.NumberOfMatchedStations
        && muCut["NumberOfValidPixelHits"]        > cutLevel.NumberOfValidPixelHits
        && muCut["TrackLayersWithMeasurement"]    > cutLevel.TrackLayersWithMeasurement
        && fabs(mu->Dxy(pvPosition))              < cutLevel.dxy
        && fabs(mu->Dz(pvPosition))               < cutLevel.dz
       ) muPass = true;
  }
  return muPass;
}

bool Zfinder::PassMuonIso(TCMuon *mu, muIsoCuts cutLevel){
  float thisEA = 0;
  if (fabs(mu->Eta())     <  1.0) thisEA = EAMu[0];
  else if (fabs(mu->Eta())     <  1.5) thisEA = EAMu[1];
  else if (fabs(mu->Eta())     <  2.0) thisEA = EAMu[2];
  else if (fabs(mu->Eta())     <  2.2) thisEA = EAMu[3];
  else if (fabs(mu->Eta())     <  2.3) thisEA = EAMu[4];
  else if (fabs(mu->Eta())     <  2.4) thisEA = EAMu[5];

  mu->LoadMap(*mapCuts);
  map<string, float> muCut = mu->MuonMap();
  float combIso;
  if ( (suffix.find("Signal") != string::npos) || (period.find("2011") != string::npos) ){
    combIso = (muCut["chIso04"] + max(0.,(double)muCut["nhIso04"] + muCut["phIso04"] - rhoFactor*thisEA));
  }else{
    combIso = (muCut["chIso04"] + max(0.,(double)muCut["nhIso04"] + muCut["phIso04"] - rhoFactorMu*thisEA));
  }
  bool isoPass = false;
  if (combIso/mu->Pt() < cutLevel.relCombIso04) isoPass = true;
  return isoPass;
}
