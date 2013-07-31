#!/bin/csh

set dir=`echo $1 | cut -d _ -f 1 `

cp $6 higgsAnalyzerV2.C

sed -i "s/SUFFIX/$1/g" higgsAnalyzerV2.C
sed -i "s/ABCD/$2/g" higgsAnalyzerV2.C
sed -i "s/SELECTION/$4/g" higgsAnalyzerV2.C
sed -i "s/PERIOD/$5/g" higgsAnalyzerV2.C
sed -i "s/DATANAME/$3/g" higgsAnalyzerV2.C
sed -i "s/COUNT/local/g" higgsAnalyzerV2.C

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
     gROOT->LoadMacro("plugins/MetCorrector.cc+");
     gROOT->LoadMacro("plugins/HistManager.cc+");
     gROOT->LoadMacro("plugins/WeightUtils.cc+");
     gROOT->LoadMacro("plugins/TriggerSelector.cc+");
     gROOT->LoadMacro("plugins/ElectronFunctions.cc+");
     gROOT->LoadMacro("plugins/rochcor_2011.cc+");
     gROOT->LoadMacro("plugins/rochcor2012v2.C+");
     gROOT->LoadMacro("plugins/PhosphorCorrectorFunctor.cc+");
     gROOT->LoadMacro("plugins/LeptonScaleCorrections.h+");
     gROOT->LoadMacro("plugins/EGammaMvaEleEstimator.cc+");
     gROOT->LoadMacro("plugins/ZGAngles.cc+");

	  TChain* fChain = new TChain("ntupleProducer/eventTree");

	  ifstream sourceFiles("sourceFiles/$3.txt");
	  string line;
     int  count = 0;
	  cout<<"Adding files from $3 to chain..."<<endl;
    
    while (sourceFiles >> line) {
      if (count == 0 && line.find("dcache")==string::npos){
        float rhoFactor;
		    TBranch        *b_rhoFactor;   //!
        TFile fixFile(line.c_str(),"open");
        TTree *fixTree = (TTree*)fixFile.Get("ntupleProducer/eventTree");
        fixTree->SetBranchAddress("rhoFactor",&rhoFactor,&b_rhoFactor);
        for(int i =0; i<fixTree->GetEntries();i++){
          fixTree->GetEntry(i);
        }
        delete fixTree;

      }
	    fChain->Add(line.c_str());      
       ++count;
	  }
	  cout<<count<<" files added!"<<endl;
	  sourceFiles.close();

	  TStopwatch timer;
	  timer.Start();

	  fChain->Process("higgsAnalyzerV2.C+");

	  cout << "\n\nDone!" << endl;
	  cout << "CPU Time : " << timer.CpuTime() << endl;
	  cout << "RealTime : " << timer.RealTime() << endl;
     cout << "\n";
	}
                                          
+EOF

root -l -b -q run.C

rm run.C
mv *local.root localHistos/.
