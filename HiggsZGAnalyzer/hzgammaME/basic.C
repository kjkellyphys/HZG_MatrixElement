void basic(){
TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
dir.ReplaceAll("basic2.C","");
dir.ReplaceAll("/./","/");
   
TFile *f = new TFile("basic2.root","RECREATE");
TH1F *h1 = new TH1F("h1","x distribution",100,-4,4);
TTree *T = new TTree("ntuple","data from ascii file");
Long64_t nlines = T->ReadFile(Form("DiscCuts.txt",dir.Data()),"Cut:N_S:N_B:SB:Ov:Sq");
printf(" found %lld points\n",nlines);
T->Draw("Sq:Cut");
T->Write();
}
