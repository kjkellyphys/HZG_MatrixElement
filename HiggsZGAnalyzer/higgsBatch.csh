#!/bin/csh
source /uscmst1/prod/sw/cms/cshrc prod
scram pro CMSSW CMSSW_4_2_7
cd CMSSW_4_2_7/src
cmsenv 

set dir = /uscms/home/bpollack/work/CMSSW_4_2_8/src/Higgs/UserCode/NatesAnalyzer/HiggsZGAnalyzer

cp $dir/higgsAnalyzer_Template.C ./higgsAnalyzer.C
cp $dir/higgsAnalyzer.h .
cp $dir/../src/TC*.cc .
cp $dir/../src/TC*.h .
cp $dir/sourceFilesOct09/$3.txt .
mkdir data
cp -r $dir/data/*root data/.
mkdir plugins
cp -r $dir/plugins/ReducedMETComputer.* plugins/.
cp -r $dir/plugins/MetCorrector.* plugins/.
cp -r $dir/plugins/HistManager.* plugins/.

sed -i "s/SUFFIX/$1/g" higgsAnalyzer.C
sed -i "s/TRIGGER/$2/g" higgsAnalyzer.C
sed -i "s/SELECTION/$4/g" higgsAnalyzer.C

cat > run.C << +EOF
    
	#include <iostream>
	#include <fstream>
	#include <string>
	#include <vector>
	#include <cstdlib>

	using namespace std;

	void run() {

	  gROOT->LoadMacro("TCJet.cc+");
	  gROOT->LoadMacro("TCMET.cc+");
	  gROOT->LoadMacro("TCGenJet.cc+");
	  gROOT->LoadMacro("TCElectron.cc+");
	  gROOT->LoadMacro("TCMuon.cc+");
      gROOT->LoadMacro("TCTau.cc+");
      gROOT->LoadMacro("TCPhoton.cc+");
	  gROOT->LoadMacro("TCPrimaryVtx.cc+");
	  gROOT->LoadMacro("TCTriggerObject.cc+");
	  gROOT->LoadMacro("plugins/MetCorrector.cc+");
	  gROOT->LoadMacro("plugins/HistManager.cc+");

	  TChain* fChain = new TChain("ntupleProducer/eventTree");

	  ifstream sourceFilesOct09("$3.txt");
	  char line[100];
	  while (sourceFilesOct09 >> line) {
	    fChain->Add(line);      
	  }
	  sourceFilesOct09.close();

	  TStopwatch timer;
	  timer.Start();

	  fChain->Process("higgsAnalyzer.C+");
	}
                                          
+EOF

root -l -b -q run.C
cp higgsHistograms.root $dir/histos/batchHists/$4/higgsHistograms_$3.root
rm -r ./*
