
#include <iostream>
#include <map>

#include "TROOT.h"
#include "TApplication.h"

#include "TFile.h"
#include "TChain.h"
#include "TString.h"

#include "TH2.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLorentzVector.h"

#include "LHEF.h"

#include "ExRootAnalysis/ExRootClasses.h"

#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"

#include "ExRootAnalysis/ExRootUtilities.h"
#include "ExRootAnalysis/ExRootProgressBar.h"

using namespace std;

//---------------------------------------------------------------------------

void AnalyseEvent(LHEF::Reader *reader, ExRootTreeBranch *branch, Long64_t eventNumber)
{
  const LHEF::HEPEUP &hepeup = reader->hepeup;

  TRootLHEFEvent *element;

  element = (TRootLHEFEvent*) branch->NewEntry();


  element->Number = eventNumber;
  element->Nparticles = hepeup.NUP;
  element->ProcessID = hepeup.IDPRUP;
  element->Weight = hepeup.XWGTUP;
  element->ScalePDF = hepeup.SCALUP;
  element->CouplingQED = hepeup.AQEDUP;
  element->CouplingQCD = hepeup.AQCDUP;
}

//---------------------------------------------------------------------------

void AnalyseParticles(LHEF::Reader *reader, ExRootTreeBranch *branch)
{
  const LHEF::HEPEUP &hepeup = reader->hepeup;

  Int_t particle;
  Double_t signPz, cosTheta;

  TLorentzVector momentum;

  TRootLHEFParticle *element;

  for(particle = 0; particle < hepeup.NUP; ++particle)
  {
    element = (TRootLHEFParticle*) branch->NewEntry();

    element->PID = hepeup.IDUP[particle];
    element->Status = hepeup.ISTUP[particle];
    element->Mother1 = hepeup.MOTHUP[particle].first;
    element->Mother2 = hepeup.MOTHUP[particle].second;
    element->ColorLine1 = hepeup.ICOLUP[particle].first;
    element->ColorLine2 = hepeup.ICOLUP[particle].second;
    element->Px = hepeup.PUP[particle][0];
    element->Py = hepeup.PUP[particle][1];
    element->Pz = hepeup.PUP[particle][2];
    element->E = hepeup.PUP[particle][3];
    element->M = hepeup.PUP[particle][4];

    momentum.SetPxPyPzE(element->Px, element->Py, element->Pz, element->E);

    cosTheta = TMath::Abs(momentum.CosTheta());
    signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;

    element->PT = momentum.Perp();
    element->Eta = (cosTheta == 1.0 ? signPz*999.9 : momentum.Eta());
    element->Phi = (cosTheta == 1.0 ? signPz*999.9 : momentum.Rapidity());
    element->Rapidity = (cosTheta == 1.0 ? signPz*999.9 : momentum.Rapidity());

    element->LifeTime = hepeup.VTIMUP[particle];
    element->Spin = hepeup.SPINUP[particle];
  }
}

//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char appName[] = "ExRootLHEFConverter";

  if(argc != 3)
  {
    cout << " Usage: " << appName << " input_file" << " output_file" << endl;
    cout << " input_file - input file in LHEF format," << endl;
    cout << " output_file - output file in ROOT format." << endl;
    return 1;
  }

  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);

  // Open a stream connected to an event file:
  ifstream inputFileStream(argv[1]);

  // Create the Reader object:
  LHEF::Reader *inputReader = new LHEF::Reader(inputFileStream);

  TFile *outputFile = TFile::Open(argv[2], "RECREATE");
  ExRootTreeWriter *treeWriter = new ExRootTreeWriter(outputFile, "LHEF");

  // generated event from LHEF
  ExRootTreeBranch *branchEvent = treeWriter->NewBranch("Event", TRootLHEFEvent::Class());

  // generated partons from LHEF
  ExRootTreeBranch *branchParticle = treeWriter->NewBranch("Particle", TRootLHEFParticle::Class());

  cout << "** Calculating number of events to process. Please wait..." << endl;
  Long64_t allEntries = inputReader->getNumberOfEvents();
  cout << "** Input file contains " << allEntries << " events" << endl;

  if(allEntries > 0)
  {
    ExRootProgressBar progressBar(allEntries);
    
    // Loop over all events
    Long64_t entry = 0;
    while(inputReader->readEvent())
    {
      treeWriter->Clear();

      AnalyseEvent(inputReader, branchEvent, entry + 1);
      AnalyseParticles(inputReader, branchParticle);

      treeWriter->Fill();

      progressBar.Update(entry);

      ++entry;
    }

    progressBar.Finish();
  }

  treeWriter->Write();

  cout << "** Exiting..." << endl;

  delete treeWriter;
  delete outputFile;
  delete inputReader;
}



