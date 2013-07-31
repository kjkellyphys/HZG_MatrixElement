#!/bin/csh
source /uscmst1/prod/sw/cms/cshrc prod
scram pro CMSSW CMSSW_5_3_6
cd CMSSW_5_3_6/src
cmsenv 
cd ${_CONDOR_SCRATCH_DIR}
#### Leave this blank #######

#############################
set count     = $1
set dataName  = $2

set suffix    = $3
set abcd      = $4
set selection = $5
set period    = $6

#mkdir otherHistos
#mkdir plugins
#mkdir printouts


tar -zxf stageball.tar.gz

cp higgsAnalyzer_Template.C higgsAnalyzer.C
#cp higgsAnalyzer.h .
#cp TC*.cc .
#cp TC*.h .
#cp -r otherHistos/*root otherHistos/.
#cp -r plugins/* plugins/.

#cp $srcDir/higgsAnalyzer_Template.C ./higgsAnalyzer.C
#cp $srcDir/higgsAnalyzer.h .
#cp $srcDir/TC*.cc .
#cp $srcDir/TC*.h .
#cp -r $srcDir/otherHistos/*root otherHistos/.
#cp -r $srcDir/plugins/* plugins/.

sed -i "s/SUFFIX/$suffix/g" higgsAnalyzer.C
sed -i "s/ABCD/$abcd/g" higgsAnalyzer.C
sed -i "s/SELECTION/$selection/g" higgsAnalyzer.C
sed -i "s/PERIOD/$period/g" higgsAnalyzer.C
sed -i "s/COUNT/$count/g" higgsAnalyzer.C
sed -i "s/DATANAME/$dataName/g" higgsAnalyzer.C
sed -i "s/\.\.\/src/src/g" higgsAnalyzer.h

ls
cat > run.C << +EOF
    
  #include <iostream>
  #include <fstream>
  #include <string>
  #include <vector>
  #include <cstdlib>

  using namespace std;

  void run() {

    gROOT->LoadMacro("src/TCPhysObject.cc+");
    gROOT->LoadMacro("src/TCJet.cc+");
    gROOT->LoadMacro("src/TCMET.cc+");
    gROOT->LoadMacro("src/TCElectron.cc+");
    gROOT->LoadMacro("src/TCMuon.cc+");
    gROOT->LoadMacro("src/TCTau.cc+");
    gROOT->LoadMacro("src/TCPhoton.cc+");
    gROOT->LoadMacro("src/TCGenJet.cc+");
    gROOT->LoadMacro("src/TCGenParticle.cc+");
    gROOT->LoadMacro("src/TCPrimaryVtx.cc+");
    gROOT->LoadMacro("src/TCTriggerObject.cc+");
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

    TChain* fChain = new TChain("ntupleProducer/eventTree");

    ifstream sourceFiles("input.txt");
    string line;
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
    }
    sourceFiles.close();

    TStopwatch timer;
    timer.Start();

    fChain->Process("higgsAnalyzer.C+");
  }
                                          
+EOF

root -l -b -q run.C


rm higgsAnalyzer*
rm input.txt 
rm run.C
