void LoadLib(){
 // load dielectron lib
 gSystem->Load("libPWGDQdielectron.so");
 gROOT->ProcessLine(".include $ALICE_PHYSICS/include/");
}


