#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>


using namespace std;

void runLocal(int nEvents = 999999999) 
{

  gROOT->LoadMacro("../src/TCJet.cc+");
  gROOT->LoadMacro("../src/TCMET.cc+");
  gROOT->LoadMacro("../src/TCElectron.cc+");
  gROOT->LoadMacro("../src/TCMuon.cc+");
  gROOT->LoadMacro("../src/TCTau.cc+");
  gROOT->LoadMacro("../src/TCPhoton.cc+");
  gROOT->LoadMacro("../src/TCGenJet.cc+");
  gROOT->LoadMacro("../src/TCPrimaryVtx.cc+");
  gROOT->LoadMacro("../src/TCTriggerObject.cc+");

  TChain* fChain = new TChain("ntupleProducer/eventTree");

  // add files to be analyzed
  fChain->Add("/uscms/home/bpollack/nobackup/Data/DoubleMu_Run2011A_May10_v1/nuTuple_100_1_TaM.root");      

  TStopwatch timer;
  timer.Start();

  fChain->Process("higgsAnalyzer.C+", "", nEvents, 0);

  cout << "\n\nDone!" << endl;
  cout << "CPU Time : " << timer.CpuTime() << endl;
  cout << "RealTime : " << timer.RealTime() << endl;
  cout << "\n";
}
