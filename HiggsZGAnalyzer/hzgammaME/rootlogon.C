{
  gSystem->SetIncludePath("-I$ROOTSYS/include -I../../txtFiles/Higgs_CS_and_Width/include");
  gSystem->Load("libgfortran.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libEG.so");
  gSystem->Load("./libmcfm_6p6.so");
  gSystem->Load("./libME.so");
  gSystem->AddIncludePath("-I$ROOFITSYS/include/");
}
