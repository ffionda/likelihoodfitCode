#include "FitCDFLikelihoodPbPb.C+"

void doLikelihoodFitMVA(TString candType,Double_t ptmin, Double_t ptmax, Double_t invMmin,Double_t invMmax, Bool_t isQuadratic, Double_t centmin, Double_t centmax, Bool_t fixFsig=kFALSE, Int_t mBkgOpt=0 /*Bool_t useMixedEvent=kFALSE*/){
 // load lib and execute likelihood fit using MVA approach
 //gSystem->Load("libPWGDQdielectron.so");
 gROOT->ProcessLine(".include $ALICE_PHYSICS/include/");
 FitCDFLikelihoodMVA(candType,ptmin,ptmax,invMmin,invMmax,isQuadratic,centmin,centmax,fixFsig,mBkgOpt);
 return;
}
