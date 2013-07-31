#!/bin/csh

cp simple_Template.C simple_v2.C

sed -i "s/SUFFIX/$1/g" simple_v2.C

cat > run.C << +EOF
    
	#include <iostream>
	#include <fstream>
	#include <string>
	#include <vector>
	#include <cstdlib>

	using namespace std;

	void run() {

     gROOT->LoadMacro("../src/TCPhysObject.cc+");
     gROOT->LoadMacro("../src/TCJet.cc+");
     gROOT->LoadMacro("../src/TCMET.cc+");
     gROOT->LoadMacro("../src/TCElectron.cc+");
     gROOT->LoadMacro("../src/TCMuon.cc+");
     gROOT->LoadMacro("../src/TCTau.cc+");
     gROOT->LoadMacro("../src/TCPhoton.cc+");
     gROOT->LoadMacro("../src/TCGenJet.cc+");
     gROOT->LoadMacro("../src/TCGenParticle.cc+");
     gROOT->LoadMacro("../src/TCPrimaryVtx.cc+");
     gROOT->LoadMacro("../src/TCTriggerObject.cc+");

     gROOT->LoadMacro("plugins/HistManager.cc+");
     gROOT->LoadMacro("plugins/ZGAngles.cc+");

	  TChain* fChain = new TChain("ntupleProducer/eventTree");

	  ifstream sourceFiles("sourceFiles/$1.txt");
	  char line[128];
     int  count = 0;
	  cout<<"Adding files from $1 to chain..."<<endl;

	  while (sourceFiles >> line) {
	    fChain->Add(line);      
       ++count;
	  }
	  cout<<count<<" files added!"<<endl;
	  sourceFiles.close();

	  TStopwatch timer;
	  timer.Start();

	  fChain->Process("simple_v2.C+");

	  cout << "\n\nDone!" << endl;
	  cout << "CPU Time : " << timer.CpuTime() << endl;
	  cout << "RealTime : " << timer.RealTime() << endl;
     cout << "\n";
	}
                                          
+EOF

root -l -b -q run.C

rm run.C
rm simple_v2.C
