#!/bin/csh
source /uscmst1/prod/sw/cms/cshrc prod
scram pro CMSSW CMSSW_4_2_8
cd CMSSW_4_2_8/src
cmsenv 
#### Leave this blank #######

#############################
set srcDir    = $1
set outDir    = $2
set count     = $3
set dataName  = $4

set suffix    = $5
set trigger   = $6
set selection = $7
set period    = $8

mkdir data
mkdir plugins
mkdir printouts

cp $srcDir/zLumi_Template.C ./higgsAnalyzer.C
cp $srcDir/higgsAnalyzer.h .
cp $srcDir/../src/TC*.cc .
cp $srcDir/../src/TC*.h .
cp -r $srcDir/data/*root data/.
cp -r $srcDir/plugins/*.* plugins/.
#cp -r $srcDir/mvaWeights .

sed -i "s/SUFFIX/$suffix/g" higgsAnalyzer.C
sed -i "s/TRIGGER/$trigger/g" higgsAnalyzer.C
sed -i "s/SELECTION/$selection/g" higgsAnalyzer.C
sed -i "s/PERIOD/$period/g" higgsAnalyzer.C

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
    gROOT->LoadMacro("TCGenParticle.cc+");
    gROOT->LoadMacro("TCElectron.cc+");
    gROOT->LoadMacro("TCMuon.cc+");
    gROOT->LoadMacro("TCTau.cc+");
    gROOT->LoadMacro("TCPhoton.cc+");
    gROOT->LoadMacro("TCPrimaryVtx.cc+");
    gROOT->LoadMacro("TCTriggerObject.cc+");
	  gROOT->LoadMacro("plugins/MetCorrector.cc+");
	  gROOT->LoadMacro("plugins/HistManager.cc+");
    gROOT->LoadMacro("plugins/WeightUtils.cc+");
    gROOT->LoadMacro("plugins/TriggerSelector.cc+");

    TChain* fChain = new TChain("ntupleProducer/eventTree");

    ifstream sourceFiles("input.txt");
    char line[1000];
    while (sourceFiles >> line) {
        //cout<<line<<endl;
      fChain->Add(line);      
    }
    sourceFiles.close();

    TStopwatch timer;
    timer.Start();

    fChain->Process("higgsAnalyzer.C+");
  }
                                          
+EOF

root -l -b -q run.C

cp higgsHistograms.root $outDir/higgsHistograms_${dataName}_$count.root
rm -r ./*


