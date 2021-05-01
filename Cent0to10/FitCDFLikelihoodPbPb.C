#include "fstream"
#include "TLatex.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TSystem.h"
#include "TMath.h"
#include "TFile.h"
#include "AliDielectronBtoJPSItoEle.h"
#include "AliDielectronBtoJPSItoEleCDFfitFCN.h"
#include "AliDielectronBtoJPSItoEleCDFfitHandler.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TFitter.h"
#include "TMinuit.h"
#include "TLegend.h"

// NEW PID SETTINGS
// Input files containing TNtuple for both data and MC. 
// The Ntuple can be produced using the method 
// MakeNtupleFromDstFilteredTree() and running on the 
// standard output of the AliReducedAnalysisFilterTrees
//

//Use the tigther PID cuts
TString inputDistr = "./InputRootFiles/NtuplePP5TeV_DATA_OS_tightPID.root"; //
TString inputDistrMC = "./InputRootFiles/NtuplePP5TeV_MC_OS_tightPID.root"; // only prompt J/psi

enum fitParameters {
// x background
kWResolution=0, kWPosL, kWNegL, kWSymL, kLPos, kLNeg, kLSym, 
// fb ,fsig
kFb, kFsig,
// signal invariant mass 
kMean, kNexp, kSigma, kAlpha, kNormSig, 
// background invariant mass
kNormBkg, kMeanBkg, kSlopeBkg, kConstBkg, 
// resolution FF
kWG1FF, kWG2FF, kMeanG1FF, kResG1FF, kMeanG2FF, kResG2FF, kAlphaResFF, kLambdaResFF, kNormPowerLawFF, 
// resolution FS
kWG1FS, kWG2FS, kMeanG1FS, kResG1FS, kMeanG2FS, kResG2FS, kAlphaResFS, kLambdaResFS, kNormPowerLawFS, 
// resolution SS
kWG1SS, kWG2SS, kMeanG1SS, kResG1SS, kMeanG2SS, kResG2SS, kAlphaResSS, kLambdaResSS, kNormPowerLawSS, 
// additional parameters background 
kLSym1, kWSym1L, 
// additional parameters mass 
kx4PolMass, kx5PolMass,
// total number of parameters
kNumPar
};

AliDielectronBtoJPSItoEleCDFfitHandler* fithandler=0x0;
const int numParam = kNumPar;

Double_t paramInputValues[numParam];

Bool_t useExpoBkgMass = kFALSE; // use exponential function for background (if false polynomial function is used)
Int_t kPolynOrd = 3; // polyn. order for invariant mass bkg (3, 4, or 5)

Bool_t saveBkgParametersForPtIntegratedCase = kTRUE;

TString fitParameterNames[numParam]={
// x background
"WeightResolution", "WPosLm", "WNegLm", "WSymLm", "LPos", "LNeg", "LSym",
// fb ,fsig
"Fb", "Fsig",
// signal invariant mass 
"Mean", "Nexp", "Sigma", "Alpha", "NormSig",
// background invariant mass
useExpoBkgMass ? "NormBkg" : "xPolynOrd0",useExpoBkgMass ? "MeanBkg" : "xPolynOrd1", useExpoBkgMass ? "SlopeBkg" : "xPolynOrd2", useExpoBkgMass ? "ConstBkg" : "xPolynOrd3",
// resolution FF
"WG1FF", "WG2FF", "MeanG1FF", "ResG1FF", "MeanG2FF", "ResG2FF", "AlphaResFF", "LambdaResFF", "NormPowerLawFF",
// resolution FS
"WG1FS", "WG2FS", "MeanG1FS", "ResG1FS", "MeanG2FS", "ResG2FS", "AlphaResFS", "LambdaResFS", "NormPowerLawFS",
// resolution SS
"WG1SS", "WG2SS", "MeanG1SS", "ResG1SS", "MeanG2SS", "ResG2SS", "AlphaResSS", "LambdaResSS", "NormPowerLawSS",
// additional parameters background (extra symmetric part, if needed)
"LSym1", "WSym1L",
// additional parameters mass (polynomial bkg) 
"xPolynOrd4", "xPolynOrd5",
};


Int_t invMassSignalParam[]={kMean, kNexp, kSigma, kAlpha, kNormSig}; // inv mass signal parameters
Int_t invMassBkgExpoParam[]={kNormBkg, kMeanBkg, kSlopeBkg, kConstBkg, kx4PolMass, kx5PolMass}; // inv mass bkg parameters

//Int_t psProperBkgParam[]={kWPosL, kWNegL, kWSymL, kLPos, kLNeg, kLSym}; // x-background parameters to be fitted
Int_t psProperBkgParam[]={kWPosL, kWNegL, kWSymL, kWSym1L, kLPos, kLNeg, kLSym, kLSym1}; // x-background parameters to be fitted [+additional symmetric part]

Int_t resolutionParamFF[]={kWG1FF, kMeanG1FF, kResG1FF, kWG2FF, kMeanG2FF, kResG2FF, kAlphaResFF, kLambdaResFF, kNormPowerLawFF}; // resol parameters (FF)
Int_t resolutionParamFS[]={kWG1FS, kMeanG1FS, kResG1FS, kWG2FS, kMeanG2FS, kResG2FS, kAlphaResFS, kLambdaResFS, kNormPowerLawFS}; // resol parameters (FS)
Int_t resolutionParamSS[]={kWG1SS, kMeanG1SS, kResG1SS, kWG2SS, kMeanG2SS, kResG2SS, kAlphaResSS, kLambdaResSS, kNormPowerLawSS}; // resol parameters (SS)

// all parameters
Int_t allParam[]={kWResolution, kWPosL, kWNegL, kWSymL, kLPos, kLNeg, kLSym, 
kFb, kFsig, kMean, kNexp, kSigma, kAlpha, kNormSig, kNormBkg, kMeanBkg, kSlopeBkg, kConstBkg,
kWG1FF, kWG2FF, kMeanG1FF, kResG1FF, kMeanG2FF, kResG2FF, kAlphaResFF, kLambdaResFF, kNormPowerLawFF,
kWG1FS, kWG2FS, kMeanG1FS, kResG1FS, kMeanG2FS, kResG2FS, kAlphaResFS, kLambdaResFS, kNormPowerLawFS,
kWG1SS, kWG2SS, kMeanG1SS, kResG1SS, kMeanG2SS, kResG2SS, kAlphaResSS, kLambdaResSS, kNormPowerLawSS,
kWSym1L, kLSym1,kx4PolMass, kx5PolMass};

// all parameters (x-backgroud)
Int_t psProperBkgParamAll[]={kWPosL, kWNegL, kWSymL, kWSym1L, kLPos, kLNeg, kLSym, kLSym1}; // x-background parameters (all)

//
enum fitMode {kExtrFb=0, kFitChi2XResolFF, kFitChi2XResolFS, kFitChi2XResolSS, kFitXBkg, kFitChi2SigMass, kFitChi2BkgMass, kFitBkgMass, kFitXBkgMVA, kFitXBkgMVAChi2};
// fit mode option used in the method FitCDFLikelihoodPbPb:
// kExtrFb -> likelihood fit to extract Fsig / Fb
// kFitChi2XResolFF (FS, SS) -> chi2 fit of x-resolution functions (MC) for different candidate's types
// kFitXBkg -> likelihood fit for pseudo-proper decay length background determination
// kFitXBkgMVA (use unbinned likelihood) -> likelihood fit for pseudo-proper decay length background determination separately for candidate's type - pt - mass range (for MVA)
// kFitXBkgMVAChi2 (use chi2) -> likelihood fit for pseudo-proper decay length background determination separately for candidate's type - pt - mass range (for MVA)
// kFitChi2SigMass -> chi2 fit to fix invariant mass signal shape (MC)
// kFitChi2BkgMass -> chi2 fit to fix invariant mass background shape (signal fixed from MC)
// kFitBkgMass -> likelihood fit to fix invariant mass background shape (signal fixed from MC)
//
// all fitted parameters are saved in a file which can be re-used as input to set starting parameters.
//

Double_t Fb =0.0;
Double_t Fsig =0.2;

Double_t weightType[] = {0.,0.,0.};
Double_t weightTypeSignal[] = {0.,0.,0.};

void AddNtupla(TNtuple *ntToFill, TNtuple *ntNew, Double_t ptmin=0., Double_t pmax=200., Double_t mMin = 2., Double_t mMax = 6., Double_t centmin=-1, Double_t centmax=-1);
void SaveParameters(TString filename, Int_t param[], Int_t size);
void SaveParameters(TString filename, Int_t param[], Int_t size, Double_t paramValues[]);
void SetStartingParameters(TString filename);
void SetReleasedParameters(Int_t param[], Int_t npar);

Double_t CDFResolutionFunction(Double_t *x, Double_t *par);
TF1 *SetupResolutionFunction(TString name, Int_t resolutionParam[], Int_t type);
Double_t InvariantMassSignalFunction(Double_t *x, Double_t *par);
TF1* SetupInvariantMassSignalFunction(TString name, Double_t lowEd, Double_t upEd);
Double_t InvariantMassSignalPlusBkgExpoFunction(Double_t *x, Double_t *par);
TF1 *SetupInvariantMassSignalPlusBkgFunction(TString name, Double_t lowEd, Double_t upEd, Bool_t isExpo);
TF1 *SetupBackgroundFunction(TString name, Double_t bkgParam[],TString type);
void SaveFunctions(TString resType, Double_t ptMin, Double_t ptMax, Double_t bandLow, Double_t bandUp);
Double_t CDFxBackgroundFunction(Double_t *x, Double_t *par);
Double_t GetAverageFunctions(Double_t *x, Double_t *par);
Double_t GetAverageFunctionsSignal(Double_t *x, Double_t *par);

void GetParameters(TString filename, Float_t parameters[]);
Double_t GetMeanMass(TNtuple *ntCp, Double_t massMin, Double_t massMax, Double_t ptMin, Double_t ptMax,Int_t type);
Double_t*** ComputeWeightsForInterpolation(TNtuple *ntOrig, Bool_t kQuadratic, Bool_t excludeExt);
Double_t ***LoadResolutionParameters(Double_t ptLimits[]);
void LoadBackgroundFunctionsAndWeights(AliDielectronBtoJPSItoEleCDFfitFCN *likeObj,TString filename, TNtuple *ntCand, Bool_t Quadratic, Bool_t excludeExt);

void FitCDFLikelihoodMVA(TString resType, Double_t ptMin, Double_t ptMax, Double_t bandLow, Double_t bandUp, Bool_t kQuadratic, Double_t centmin=-1, Double_t centmax=-1, Bool_t fixFsig=kFALSE, Int_t mBkgOpt = 0 );
void FitCDFLikelihoodPbPb(Int_t fitmode, Double_t ptMin, Double_t ptMax, Double_t bandLow=2.2, Double_t bandUp=4.0, TString resType="FF;FS",Double_t  centmin=-1, Double_t centmax = -1, Int_t mBkgOpt = 0);
// global var
AliDielectronBtoJPSItoEleCDFfitFCN *likely_obj = 0x0;
AliDielectronBtoJPSItoEleCDFfitFCN *likely_obj_proj = 0x0;
Bool_t computeIntMass = kTRUE; 
const Int_t kNbinsMax = 10;

Double_t ptEdges[]={1.5,3.0};  // pt range(s) //CHANGED, uncomment if pt-independent integrated
//Double_t ptEdges[]={1.5,3.0,5.,10.};  // pt range(s)
Double_t massBins[]={2.6,2.7,2.8, 3.2,3.4,3.6}; // low mass (0) - interp. region (1) - high mass (2)
//Double_t massBins[]={2.70, 2.75, 2.82, 3.2, 3.26, 3.4}; // low mass (0) - interp. region (1) - high mass (2)
Int_t extrRegion = 2;  // interp. region
Int_t kPtBins = sizeof(ptEdges)/sizeof(Double_t)-1;
Int_t kMassBins = sizeof(massBins)/sizeof(Double_t)-1;
Int_t kTypes = 3;
Float_t nCandPtMassWnd[kNbinsMax][kNbinsMax][kNbinsMax];
Float_t nCandPtMassWndSignal[kNbinsMax][kNbinsMax];

void  FitCDFLikelihoodPbPb(Int_t fitmode, Double_t ptMin, Double_t ptMax, Double_t bandLow, Double_t bandUp, TString resType, Double_t  centmin, Double_t centmax, Int_t mBkgOpt){

  ///////////////////////////////////////////////////////////////////
  //
  // Example macro to read local N-Tuples of JPSI 
  // and bkg candidates and perform log-likelihood 
  // minimization using the minimization handler class
  // origin: fiorella.fionda@cern.ch
  //
  ///////////////////////////////////////////////////////////////////

  // Set inv mass / pt limits
  //inv mass functions normalized to  1 between bandLow - bandUp 
  // 
  // projection of pseudoproper decay length functions drawn in the signal region 
  Double_t bandLowSignal=2.92;
  Double_t bandUpSignal=3.16;
  if(kPtBins == 1) {ptEdges[0] = ptMin; ptEdges[1] = ptMax;}
  ifstream filenameEdges(Form("inputFiles_%1.1f_%1.1f/massEdges.txt",ptMin,ptMax));
  if(filenameEdges) { for(int im=0; im<kMassBins+1; im++) { filenameEdges >> massBins[im]; /*printf("mass %f \n",massBins[im]);*/} }
  //
 TNtuple *nt=new TNtuple("ntuplaSigna_new","NtuplaSignal","Xdecaytime:Mass:Type:Pt");
 TString inputFileName = inputDistr; 
  //if(fitmode == kFitChi2XResolFF || fitmode == kFitChi2XResolFS || fitmode==kFitChi2XResolSS || fitmode == kFitChi2SigMass) { inputFileName = inputDistrMC; centmin = -1; centmax = -1; }
  if(fitmode == kFitChi2XResolFF || fitmode == kFitChi2XResolFS || fitmode==kFitChi2XResolSS || fitmode == kFitChi2SigMass) {
      inputFileName = inputDistrMC;
      if(fitmode == kFitChi2SigMass) {centmin = -1; centmax = -1;} 
      
}
  TFile f(inputFileName);
  TFile ftemplate(Form("inputFiles_%1.1f_%1.1f/XtemplateNonPromptJpsi.root",kPtBins > 1 ? ptEdges[0] : ptMin, kPtBins > 1 ? ptEdges[kPtBins] : ptMax)); //template non-prompt J/psi
  TH1F *hCsiMCEvtGen = (TH1F*)ftemplate.Get("psTemplate");

  Double_t integral = 0;
  for(int i=1;i<hCsiMCEvtGen->GetNbinsX()+1; i++) integral += (hCsiMCEvtGen->GetBinContent(i)*hCsiMCEvtGen->GetBinWidth(i));
  hCsiMCEvtGen->Scale(1./integral);
  
  Double_t* x=0x0; Double_t* m=0x0; Double_t *pt =0x0; Int_t* type=0; Int_t n=0; 
  AliDielectronBtoJPSItoEle* aBtoJPSItoEle =new AliDielectronBtoJPSItoEle();
  // set all starting parameters 
  SetStartingParameters(Form("inputFiles_%1.1f_%1.1f/startingParameters.root",kPtBins > 1 ? ptEdges[0] : ptMin, kPtBins > 1 ? ptEdges[kPtBins] : ptMax));
  // set inv mass signal parameters 
  if(fitmode != kFitChi2SigMass) SetStartingParameters(Form("inputFiles_%1.1f_%1.1f/InvariantMassSignalMC.root",kPtBins > 1 ? ptEdges[0]: ptMin, kPtBins > 1 ? ptEdges[kPtBins] : ptMax));
  // set inv mass signal parameters
  if(fitmode < kFitChi2SigMass || fitmode > kFitBkgMass) { 
      TString resTypeSt = resType; 
      resTypeSt.ReplaceAll(";","_");
      if(mBkgOpt == 1 ) SetStartingParameters(Form("inputFiles_%1.1f_%1.1f/InvariantMassBkgME%s_.root",kPtBins > 1 ? ptEdges[0]: ptMin, kPtBins > 1 ? ptEdges[kPtBins] : ptMax,resTypeSt.Data()));
      else if(mBkgOpt == 2 ) SetStartingParameters(Form("inputFiles_%1.1f_%1.1f/InvariantMassBkgLS%s_.root",kPtBins > 1 ? ptEdges[0]: ptMin, kPtBins > 1 ? ptEdges[kPtBins] : ptMax,resTypeSt.Data()));
      else SetStartingParameters(Form("inputFiles_%1.1f_%1.1f/InvariantMassBkgExpoChi2Fit.root",kPtBins > 1 ? ptEdges[0]: ptMin, kPtBins > 1 ? ptEdges[kPtBins] : ptMax));
  } 
  
  //if(fitmode < kFitChi2SigMass) SetStartingParameters(Form("inputFiles_%1.1f_%1.1f/InvariantMassBkgExpoChi2Fit.root",kPtBins > 1 ? ptEdges[0]: ptMin, kPtBins > 1 ? ptEdges[kPtBins] : ptMax));
  //if(fitmode < kFitChi2SigMass) SetStartingParameters(Form("inputFiles_%1.1f_%1.1f/InvariantMassBkgME.root",kPtBins > 1 ? ptEdges[0]: ptMin, kPtBins > 1 ? ptEdges[kPtBins] : ptMax));


  // set resol parameters -> overwrite only resolution parameters
  // 
  if(fitmode == kFitXBkg || fitmode == kFitXBkgMVA || fitmode == kExtrFb || fitmode == kFitXBkgMVAChi2){
  SetStartingParameters(Form("inputFiles_%1.1f_%1.1f/XResolParametersFF_%1.1f_%1.1f.root",kPtBins > 1 ? ptEdges[0]: ptMin, kPtBins > 1 ? ptEdges[kPtBins] : ptMax, ptMin, ptMax));
  SetStartingParameters(Form("inputFiles_%1.1f_%1.1f/XResolParametersFS_%1.1f_%1.1f.root",kPtBins > 1 ? ptEdges[0]: ptMin, kPtBins > 1 ? ptEdges[kPtBins] : ptMax, ptMin, ptMax));
  SetStartingParameters(Form("inputFiles_%1.1f_%1.1f/XResolParametersSS_%1.1f_%1.1f.root",kPtBins > 1 ? ptEdges[0]: ptMin, kPtBins > 1 ? ptEdges[kPtBins] : ptMax, ptMin, ptMax)); 
  //
  if( !(ptMin==1.5 && ptMax==10.0) ){ 
      if(!resType.Contains(";")) SetStartingParameters(Form("inputFiles_1.5_10.0/bkgPtIntegrated/XBkgParameters_mass%1.1f_%1.1f_resType%s_pt1.5_10.0.root",bandLow,bandUp,resType.Data()));
      else SetStartingParameters(Form("inputFiles_1.5_10.0/bkgPtIntegrated/XBkgParameters_mass%1.1f_%1.1f_resTypeFF_pt1.5_10.0.root",bandLow,bandUp));
      // if(resType.Contains("FS")) SetStartingParameters(Form("inputFiles_3.0_5.0/XBkgParameters_mass%1.1f_%1.1f_resTypeFF_pt3.0_5.0.root",bandLow,bandUp));
            }else{
        // starting parameters for centrality (10-30%)
        if(centmin==10 || centmin==0){
        TString massR;
        if(bandLow < 3.0) { massR = "L";
        if(bandLow == 2.72) massR.Append("1"); else massR.Append("2");  }
        if(bandLow > 3.0) { massR = "H";
        if(bandUp == 3.40) massR.Append("2"); else massR.Append("1");  }
        SetStartingParameters(Form("startingParametersBkg/XBkgParameters_mass%s_%s_pt1.5_10.0.root",massR.Data(),resType.Data()));
        }
      
}
  }

  TNtuple *ntOrig=0x0;
  ntOrig=(TNtuple*)f.Get("fNtupleJPSI");
  if(fitmode==kFitXBkg) {AddNtupla(nt,ntOrig,ptMin,ptMax,2.2,2.6,centmin,centmax);  AddNtupla(nt,ntOrig,ptMin,ptMax,3.2,4.0,centmin,centmax); } // add candidates within low+high mass sidebands (pp)
  else AddNtupla(nt,ntOrig,ptMin,ptMax,bandLow,bandUp,centmin,centmax);  // add candidates in [bandLow, bandUp]

  // if resolution function is fitted the corresponding candidate's type is considered
  if(fitmode == kFitChi2XResolFF) resType="FF"; 
  else if(fitmode == kFitChi2XResolFS) resType="FS"; 
  else if(fitmode == kFitChi2XResolSS) resType="SS"; 
  aBtoJPSItoEle->SetResTypeAnalysis(resType);
   
  aBtoJPSItoEle->ReadCandidates(nt,x,m,pt,type,n,bandLow,bandUp); // read N-Tuples
  printf("+++\n+++ Number of total candidates (prim J/psi + secondary J/psi + bkg) ---> %d candidates \n+++\n",n);

  aBtoJPSItoEle->SetFitHandler(x,m,pt,type,n); // Set the fit handler with given values of x, m, # of candidates 

  aBtoJPSItoEle->CloneMCtemplate(hCsiMCEvtGen);    // clone MC template and copy internally
                                                   // in this way any model can be setted from outside

  aBtoJPSItoEle->SetCsiMC(); // Pass the MC template to the CDF fit function


  //AliDielectronBtoJPSItoEleCDFfitHandler* 
  fithandler = aBtoJPSItoEle->GetCDFFitHandler(); // Get the fit handler
  //
  // Set some fit options through the handler class
  //

  if(fitmode==kExtrFb || fitmode==kFitXBkg || fitmode == kFitXBkgMVA || fitmode==kFitBkgMass || fitmode == kFitXBkgMVAChi2) fithandler->SetPrintStatus(kTRUE);
  
  fithandler->SetCrystalBallFunction(kTRUE);
  fithandler->SetExponentialFunction(useExpoBkgMass); // kFALSE - polynomial function used
  
  Double_t massLow = (TDatabasePDG::Instance()->GetParticle(443)->Mass()) - bandLow;
  Double_t massHigh = bandUp - (TDatabasePDG::Instance()->GetParticle(443)->Mass());
  fithandler->SetMassWndLow(massLow);
  fithandler->SetMassWndHigh(massHigh);

   switch(fitmode)
   {
   case kExtrFb:
   // 2D likelihood fit for fB extraction
   paramInputValues[kFb] = 0.10; paramInputValues[kFsig] = 0.010;
   paramInputValues[kWSym1L]=0.; paramInputValues[kLSym1]=0.; //bkg Only
   fithandler->SetParamStartValues(paramInputValues);
   for(int ipar=0; ipar<kNumPar; ipar++) { if(ipar != kFb && ipar !=kFsig) fithandler->FixParam(ipar,kTRUE);  }
   break;
   case kFitXBkg:
   case kFitXBkgMVA:
   case kFitXBkgMVAChi2:
   // likelihood fit for fixing x-bkg parameters (common to all candidate's types)
        paramInputValues[kFsig]=0.; paramInputValues[kFb]=0.; //bkg Only
        // if(ptMin > 3.) { paramInputValues[kWSym1L]=0.; paramInputValues[kLSym1]=0.; } //bkg Only 
        fithandler->SetParamStartValues(paramInputValues);
        SetReleasedParameters(psProperBkgParam,sizeof(psProperBkgParam)/sizeof(Int_t)); 
	if(ptMin > 3.){
	fithandler->FixParam(kLPos,kTRUE);
	fithandler->FixParam(kLNeg,kTRUE);
	//fithandler->FixParam(kLSym,kTRUE);
	//fithandler->FixParam(kLSym1,kTRUE);
	//fithandler->FixParam(kWSym1L,kTRUE); 
	}
   break;
   case kFitBkgMass:
   // likelihood fit to fix inv. mass background shape 
        paramInputValues[kFsig]=0.2; paramInputValues[kFb]=0.2; //sig+bkg
        fithandler->SetParamStartValues(paramInputValues);
        Int_t size = sizeof(invMassBkgExpoParam)/sizeof(Int_t);
        invMassBkgExpoParam[size] = kFsig;
        SetReleasedParameters(invMassBkgExpoParam,size+1);
	if(kPolynOrd == 3){
	/// fix higher polyn. x^4 and x^5 
        fithandler->FixParam(kx4PolMass,kTRUE);
        fithandler->FixParam(kx5PolMass,kTRUE);
	printf("3th polyn. for bkg \n");
	}else if(kPolynOrd == 4){
        fithandler->FixParam(kx5PolMass,kTRUE); printf("4th polyn. for bkg \n");
	}else {printf("5th polyn. for bkg \n");}
   break;
   }

  //fill histos from Ntupla
  Int_t nbinsMass = (Int_t)((bandUp - bandLow)/0.04);  
  TH1F *histMass = new TH1F("histMass","Invariant Mass; InvMass[GeV]; Entries/40MeV",nbinsMass,bandLow,bandUp);
  histMass->SetLineColor(1);
  histMass->SetMarkerColor(1);
  histMass->SetMarkerStyle(20);
  histMass->SetMarkerSize(0.7);

  Double_t bwidth = 40.; Int_t nbinsX = 300;
  if(fitmode == kFitChi2XResolFF || fitmode == kFitChi2XResolFS || fitmode == kFitChi2XResolSS) { bwidth = 10.; nbinsX = 1200; }
  TH1F *histpsproper = new TH1F("psproper_decay_length",Form("psproper_decay_length_distrib(%1.1f < M < %1.1f GeV/c^{2});pseudoproper decay length [#mum];Entries/%1.0f#mum",bandLow,bandUp,bwidth),nbinsX,-6000.,6000.);
  histpsproper->SetLineColor(1);
  histpsproper->SetMarkerColor(1);
  histpsproper->SetMarkerStyle(20);
  histpsproper->SetMarkerSize(0.7);

   Double_t maxWd = 2*6000.; Double_t nBins = 2*300.;

   TH1F *histpsproperFF = new TH1F("psproper_decay_length_FF",Form("psproper_decay_length_distrib_FF(%1.1f < M < %1.1f GeV/c^{2});pseudoproper decay length for FF [#mum];Entries/%1.0f#mum",bandLow,bandUp,2.*maxWd/nBins),(Int_t)nBins,-1.*maxWd,maxWd);
  histpsproperFF->SetLineColor(1);
  histpsproperFF->SetMarkerColor(1);
  histpsproperFF->SetMarkerStyle(20);
  histpsproperFF->SetMarkerSize(0.7);


    TH1F *histpsproperFS = new TH1F("psproper_decay_length_FS",Form("psproper_decay_length_distrib_FS(%1.1f < M < %1.1f GeV/c^{2});pseudoproper decay length for FS [#mum];Entries/%1.0f#mum",bandLow,bandUp,2.*maxWd/nBins),(Int_t)nBins,-1.*maxWd,maxWd);
  histpsproperFS->SetLineColor(1);
  histpsproperFS->SetMarkerColor(1);
  histpsproperFS->SetMarkerStyle(20);
  histpsproperFS->SetMarkerSize(0.7);

    TH1F *histpsproperSS = new TH1F("psproper_decay_length_SS",Form("psproper_decay_length_distrib_SS(%1.1f < M < %1.1f GeV/c^{2});pseudoproper decay length for SS [#mum];Entries/%1.0f#mum",bandLow,bandUp,2.*maxWd/nBins),(Int_t)nBins,-1.*maxWd,maxWd);
  histpsproperSS->SetLineColor(1);
  histpsproperSS->SetMarkerColor(1);
  histpsproperSS->SetMarkerStyle(20);
  histpsproperSS->SetMarkerSize(0.7);

  TH1F *histpsproperSignal = new TH1F("psproper_decay_length_signal",Form("psproper_decay_length_distrib(%1.2f < M < %1.2f GeV/c^{2});pseudoproper decay length [#mum];Entries/%1.0f#mum",bandLowSignal,bandUpSignal,bwidth),nbinsX,-6000.,6000.);
  histpsproperSignal->SetLineColor(1);
  histpsproperSignal->SetMarkerColor(1);
  histpsproperSignal->SetMarkerStyle(20);
  histpsproperSignal->SetMarkerSize(0.7);

  Float_t mass =0.; Float_t psproper = 0.; Float_t typeCand = 0.; Int_t nb = 0;
  Float_t ptCand = 0.;
  TString arrType[]={"SS","FS","FF"};
  nt->SetBranchAddress("Xdecaytime",&psproper);
  nt->SetBranchAddress("Mass",&mass);
  nt->SetBranchAddress("Type",&typeCand);
  nt->SetBranchAddress("Pt",&ptCand);
  Int_t fNcurrent=0; Double_t nCandSel = 0 ; Double_t nCandSelSignal=0.;
  nb = (Int_t)nt->GetEvent(fNcurrent);
  cout << nb << endl;

  for(Int_t iev=0; iev<(nt->GetEntries()); iev++){
   if(resType.Contains(arrType[(Int_t)typeCand]) && mass > bandLow && mass < bandUp){
   nCandSel += 1;
   weightType[(Int_t)typeCand] += 1.;
   histMass->Fill(mass);
   histpsproper->Fill(psproper);
   if(mass > bandLowSignal && mass < bandUpSignal) {
   histpsproperSignal->Fill(psproper);
   weightTypeSignal[(Int_t)typeCand] += 1.;
   nCandSelSignal += 1;
   }
   }
   fNcurrent++;
   nb = (Int_t) nt->GetEvent(fNcurrent);
   }

  likely_obj = fithandler->LikelihoodPointer();
  likely_obj->SetAllParameters(paramInputValues);

  likely_obj->ComputeMassIntegral();
  if(fitmode==kExtrFb || fitmode==kFitXBkg || fitmode == kFitXBkgMVA || fitmode==kFitBkgMass || fitmode == kFitXBkgMVAChi2) likely_obj->PrintStatus(); 

  if(fitmode==kFitChi2XResolFF || fitmode==kFitChi2XResolFS || fitmode==kFitChi2XResolSS){
   // fit x-resolution using chi2. Parameters are saved in a file that can be 
   // used for 2D likelihood fit. 
   TF1 *resolutionFuncFit = 0x0; 
   if(fitmode==kFitChi2XResolFF) resolutionFuncFit = SetupResolutionFunction("resolutionFuncFF",resolutionParamFF,2);
   else if(fitmode==kFitChi2XResolFS) resolutionFuncFit = SetupResolutionFunction("resolutionFuncFS",resolutionParamFF,1);
   else if(fitmode==kFitChi2XResolSS) resolutionFuncFit = SetupResolutionFunction("resolutionFuncSS",resolutionParamFS,0);
   ///
   resolutionFuncFit->SetParameter(10,histpsproper->GetEntries()*histpsproper->GetBinWidth(1));
   // resolutionFuncFit->Print();
   TCanvas *cResolFitChi2 = new TCanvas("XresolutionFit","XresolutionFit");
   cResolFitChi2->SetLogy();
   Double_t rangeResFit = 6000.; 
   if(fitmode==kFitChi2XResolSS){
   resolutionFuncFit->FixParameter(1,resolutionFuncFit->GetParameter(1));
   //resolutionFuncFit->FixParameter(2,resolutionFuncFit->GetParameter(2));
   //resolutionFuncFit->FixParameter(6,resolutionFuncFit->GetParameter(6));
   //resolutionFuncFit->FixParameter(7,resolutionFuncFit->GetParameter(7));
   }
   if(ptMin > 3.){
   resolutionFuncFit->FixParameter(1,resolutionFuncFit->GetParameter(1));
   resolutionFuncFit->FixParameter(4,resolutionFuncFit->GetParameter(4));
   resolutionFuncFit->FixParameter(6,resolutionFuncFit->GetParameter(6));
   //resolutionFuncFit->FixParameter(7,resolutionFuncFit->GetParameter(7));
	  }
   TFitResultPtr rFitRes = histpsproper->Fit(resolutionFuncFit->GetName(),"S0L","L",-1.*rangeResFit,rangeResFit);
   histpsproper->GetYaxis()->SetRangeUser(0.01,histpsproper->GetMaximum()*1.3);
   histpsproper->DrawCopy("E");   
   resolutionFuncFit->Draw("same"); 
   TLatex *latex = new TLatex();
   latex->SetNDC();
   latex->SetTextFont(42);
   latex->SetTextSize(0.035);
   latex->SetLineWidth(2);
   latex->DrawLatex(0.57, 0.85, Form("#chi^{2}/dof = %4.3f ",(rFitRes->Chi2()/(Double_t)rFitRes->Ndf())));
   latex->DrawLatex(0.57, 0.75, Form("%1.1f - %1.1f GeV/c",ptMin,ptMax));
   if(fitmode==kFitChi2XResolFF) latex->DrawLatex(0.60, 0.80, "FF");
   else if(fitmode==kFitChi2XResolFS) latex->DrawLatex(0.60, 0.80, "FS");
   else if(fitmode==kFitChi2XResolSS) latex->DrawLatex(0.60, 0.80, "SS");
   cResolFitChi2->SaveAs(Form("inputFiles_%1.1f_%1.1f/resolutionFit%s_%1.1f_%1.1f.pdf",kPtBins > 1 ? ptEdges[0]: ptMin, kPtBins > 1 ? ptEdges[kPtBins] : ptMax, resType.Data(),ptMin,ptMax)); 
 
   // save parameters
   if(fitmode==kFitChi2XResolFF) SaveParameters(Form("inputFiles_%1.1f_%1.1f/XResolParametersFF_%1.1f_%1.1f",kPtBins > 1 ? ptEdges[0]: ptMin, kPtBins > 1 ? ptEdges[kPtBins] : ptMax, ptMin, ptMax),resolutionParamFF,sizeof(resolutionParamFF)/sizeof(Int_t),resolutionFuncFit->GetParameters());
   if(fitmode==kFitChi2XResolFS) SaveParameters(Form("inputFiles_%1.1f_%1.1f/XResolParametersFS_%1.1f_%1.1f",kPtBins > 1 ? ptEdges[0]: ptMin, kPtBins > 1 ? ptEdges[kPtBins] : ptMax, ptMin, ptMax),resolutionParamFS,sizeof(resolutionParamFS)/sizeof(Int_t),resolutionFuncFit->GetParameters());
   if(fitmode==kFitChi2XResolSS) SaveParameters(Form("inputFiles_%1.1f_%1.1f/XResolParametersSS_%1.1f_%1.1f",kPtBins > 1 ? ptEdges[0]: ptMin, kPtBins > 1 ? ptEdges[kPtBins] : ptMax, ptMin, ptMax),resolutionParamSS,sizeof(resolutionParamSS)/sizeof(Int_t),resolutionFuncFit->GetParameters());
   return;
  }

  if(fitmode==kFitChi2SigMass){
   TF1 *invMassSignFuncFit = SetupInvariantMassSignalFunction("invMassSignalFunc", bandLow, bandUp);
   TCanvas *cMassSigFitChi2 = new TCanvas("invMassSigFit","invMassSigFit");
   invMassSignFuncFit->SetParameter(4, histMass->GetEntries());
   histMass->Fit(invMassSignFuncFit->GetName(),"0","",bandLow,bandUp);
   histMass->DrawCopy("E");
   invMassSignFuncFit->Draw("same");   
   cMassSigFitChi2->SaveAs(Form("inputFiles_%1.1f_%1.1f/invariantMassSignalMC.pdf",kPtBins > 1 ? ptEdges[0]: ptMin, kPtBins > 1 ? ptEdges[kPtBins] : ptMax));
   SaveParameters(Form("inputFiles_%1.1f_%1.1f/InvariantMassSignalMC",kPtBins > 1 ? ptEdges[0]: ptMin, kPtBins > 1 ? ptEdges[kPtBins] : ptMax),invMassSignalParam,sizeof(invMassSignalParam)/sizeof(Int_t),invMassSignFuncFit->GetParameters());
   return;
  }

  if(fitmode==kFitChi2BkgMass){
  Int_t sizeSig = sizeof(invMassSignalParam)/sizeof(Int_t);
  Int_t sizeBkg = sizeof(invMassBkgExpoParam)/sizeof(Int_t); 
  // inv mass sig+bkg fit (signal fixed from MC)
  TF1 *invMassSigPlusBkg = SetupInvariantMassSignalPlusBkgFunction("SignalPlusBkg", bandLow, bandUp, useExpoBkgMass);
  invMassSigPlusBkg->SetParameter(sizeSig+sizeBkg+1, histMass->GetEntries()*histMass->GetBinWidth(1));
  if(kPolynOrd == 3){
  invMassSigPlusBkg->FixParameter(sizeSig+sizeBkg, 0.); // fix x^5 term 
  invMassSigPlusBkg->FixParameter(sizeSig+sizeBkg-1, 0.); // fix x^4 term 
  printf("3th polyn. for bkg \n");
  }else if(kPolynOrd == 4){
  invMassSigPlusBkg->FixParameter(sizeSig+sizeBkg, 0.); // fix x^5 term 
  printf("4th polyn. for bkg \n");
  } else{
  printf("5th polyn. for bkg \n");	 
  }
  TCanvas *cMassSigPlusBkgFitChi2 = new TCanvas("invMassSigPlusBkgFit","invMassSigPlusBkgFit");
  TFitResultPtr rPsproperBackMass = histMass->Fit(invMassSigPlusBkg->GetName(),"S0","",bandLow,bandUp);
  histMass->DrawCopy("E");
  computeIntMass = kFALSE; // not needed for drawing 
  invMassSigPlusBkg->SetLineColor(1);
  // draw signal part
  TF1 *invMassSigOnly = SetupInvariantMassSignalPlusBkgFunction("Signal", bandLow, bandUp, useExpoBkgMass);
  for(int ipar=0; ipar<sizeSig+sizeBkg+2;ipar++) { 
     if( (ipar < sizeSig+1) || (ipar == sizeSig+sizeBkg+1) ) invMassSigOnly->SetParameter(ipar, invMassSigPlusBkg->GetParameter(ipar));
     else invMassSigOnly->SetParameter(ipar, 0.);
  }
  invMassSigOnly->SetLineColor(2); 
  invMassSigOnly->Draw("same");
  // draw bkg part
  TF1 *invMassBkgOnly = SetupInvariantMassSignalPlusBkgFunction("Background", bandLow, bandUp, useExpoBkgMass);
  for(int ipar=0; ipar<sizeSig+sizeBkg+2;ipar++) {
     if(ipar < sizeSig) invMassBkgOnly->SetParameter(ipar, 0.);
     else  invMassBkgOnly->SetParameter(ipar, invMassSigPlusBkg->GetParameter(ipar)); 
  }
  invMassBkgOnly->SetLineColor(8);
  invMassBkgOnly->Draw("same");
  invMassSigPlusBkg->Draw("same");
  TLatex *tex2 = new TLatex();
  tex2->SetNDC();
  tex2->SetTextFont(42);
  tex2->SetTextSize(0.035);
  tex2->SetLineWidth(2);
  tex2->DrawLatex(0.14, 0.72, Form("#chi^{2}/dof = %4.3f ",(rPsproperBackMass->Chi2()/(Double_t)rPsproperBackMass->Ndf())));
  tex2->DrawLatex(0.14, 0.77, Form("Fsig = %4.3f +- %4.3f ",invMassSigPlusBkg->GetParameter(sizeSig), invMassSigPlusBkg->GetParError(sizeSig)));
  tex2->DrawLatex(0.14, 0.82, Form("S[%1.2f-%1.2f GeV/c^{2}] = %4.1f +- %4.1f ",bandLowSignal,bandUpSignal,invMassSigOnly->Integral(bandLowSignal,bandUpSignal)/histMass->GetBinWidth(1), TMath::Sqrt(invMassSigOnly->Integral(bandLowSignal,bandUpSignal)/histMass->GetBinWidth(1))));
  cMassSigPlusBkgFitChi2->SaveAs(Form("inputFiles_%1.1f_%1.1f/invariantMassOS_data.pdf",kPtBins > 1 ? ptEdges[0]: ptMin, kPtBins > 1 ? ptEdges[kPtBins] : ptMax));
  // save bkg parameters 
  Double_t bkgParamInvMassFromFit[sizeBkg]; 
  for(int ipar=0; ipar<sizeBkg; ipar++) bkgParamInvMassFromFit[ipar] = invMassSigPlusBkg->GetParameter(fitParameterNames[invMassBkgExpoParam[ipar]]); 
  SaveParameters(Form("inputFiles_%1.1f_%1.1f/InvariantMassBkgExpoChi2Fit",kPtBins > 1 ? ptEdges[0]: ptMin, kPtBins > 1 ? ptEdges[kPtBins] : ptMax),invMassBkgExpoParam,sizeBkg,bkgParamInvMassFromFit); 

  std::ofstream ofs;
  TString stringmass=Form("%1.2f_%1.2f",bandLow,bandUp);
  TString resTypeS = resType; resTypeS.ReplaceAll(";","_");
  stringmass.Append(Form("_type%s",resTypeS.Data()));
  //ofs.open (Form("fSIGPt%1.2f_%1.2f_%s.txt",ptMin,ptMax,stringmass.Data()), std::ofstream::out | std::ofstream::app);
  ofs.open (Form("fSIGPt%1.2f_%1.2f_%s.txt",ptMin,ptMax,stringmass.Data()), std::ofstream::out);
  ofs << Form("%f %f \n",invMassSigPlusBkg->GetParameter(sizeSig), invMassSigPlusBkg->GetParError(sizeSig)); 
  ofs.close();

  std::ofstream ofsSig2;
  stringmass=Form("%1.2f_%1.2f",2.92,3.16);
  stringmass.Append(Form("_type%s",resTypeS.Data()));
  ofsSig2.open (Form("inputFiles_%1.1f_%1.1f/fSIGPt_%s.txt",ptMin,ptMax,stringmass.Data()), std::ofstream::out);
  ofsSig2 << Form("%f \n",invMassSigOnly->Integral(2.92,3.16)/(invMassSigOnly->Integral(2.92,3.16)+invMassBkgOnly->Integral(2.92,3.16)));
  ofsSig2.close();

  std::ofstream ofsSig3;
  stringmass=Form("%1.2f_%1.2f",bandLow,bandUp);
  stringmass.Append(Form("_type%s",resTypeS.Data()));
  ofsSig3.open (Form("inputFiles_%1.1f_%1.1f/fSIGPt_%s.txt",ptMin,ptMax,stringmass.Data()), std::ofstream::out);
  ofsSig3 << Form("%f \n",invMassSigOnly->Integral(bandLow,bandUp)/(invMassSigOnly->Integral(bandLow,bandUp)+invMassBkgOnly->Integral(bandLow,bandUp)));
  ofsSig3.close();


  return;
  }


  if(fitmode == kFitXBkgMVAChi2){
   Int_t sizeBkg = sizeof(psProperBkgParamAll)/sizeof(Int_t);
   Double_t bkgStartingParameters[sizeBkg];
   for(int iparam=0; iparam<sizeBkg; iparam++) bkgStartingParameters[iparam+1] = paramInputValues[psProperBkgParamAll[iparam]];
   TF1 *psProperBkgFit = SetupBackgroundFunction("XBackground",bkgStartingParameters, resType);
   psProperBkgFit->FixParameter(0,10.);
   
      if(!(ptMin == 1.5 && ptMax == 10.0)){
     if(ptMin > 1.0){
     // 
     psProperBkgFit->FixParameter(5,psProperBkgFit->GetParameter(5));
     psProperBkgFit->FixParameter(6,psProperBkgFit->GetParameter(6));
     psProperBkgFit->FixParameter(8,psProperBkgFit->GetParameter(8));
     if(resType.Contains("FF") && ptMin == 1.5 && bandLow == 3.20) psProperBkgFit->FixParameter(7,psProperBkgFit->GetParameter(7));
    }
     
     }else{
    
	psProperBkgFit->SetParLimits(5, 0.0008,0.005);
        psProperBkgFit->SetParLimits(6, 0.0008,0.005);
        psProperBkgFit->SetParLimits(8, 0.00005,0.001);
        if(bandLow == 3.3 && resType.Contains("FF")){ 
		psProperBkgFit->SetParLimits(7, 0.001,0.1);
        	psProperBkgFit->SetParLimits(8, 0.00009,0.0005);
	}
        if( (bandLow == 2.80 || bandLow == 3.3) && resType.Contains("FS") )
        {
                psProperBkgFit->SetParLimits(7, 0.001,0.05);
                psProperBkgFit->SetParLimits(5, 0.0001,0.004);
                psProperBkgFit->SetParLimits(6, 0.0001,0.004);
        }
        if( bandLow == 2.72 && resType.Contains("FS")){  
	psProperBkgFit->SetParLimits(7, 0.001,0.1);
        psProperBkgFit->SetParLimits(5, 0.0001,0.008);
        psProperBkgFit->SetParLimits(6, 0.0001,0.008);
        }
    }
   //
   if(resType.Contains("SS")) {
     psProperBkgFit->FixParameter(5,psProperBkgFit->GetParameter(5));
     psProperBkgFit->FixParameter(6,psProperBkgFit->GetParameter(6));
   }
   psProperBkgFit->SetParameter(10,histpsproper->GetEntries()/histpsproper->GetBinWidth(1)); 
   psProperBkgFit->Print();
   TCanvas *psTotCan = new TCanvas("pseudoProperDecayLengthBkgFitChi2","pseudoProperDecayLengthBkgFitChi2");
   psTotCan->cd()->SetLogy();
   // histpsproper->Fit("XBackground","","L0");
   histpsproper->Fit("XBackground","L0");
   for(int iparam=0; iparam<sizeBkg; iparam++) { paramInputValues[psProperBkgParamAll[iparam]] = psProperBkgFit->GetParameter(iparam+1); printf("param bkg %d -> %f \n",psProperBkgParamAll[iparam],paramInputValues[psProperBkgParamAll[iparam]]); }
   fithandler->SetParamStartValues(paramInputValues);
   likely_obj->SetAllParameters(paramInputValues);
   for(int ipar=0; ipar<kNumPar; ipar++) fithandler->FixParam(ipar,kTRUE); 
 }

  // likelihood fit minimization 
  aBtoJPSItoEle->DoMinimization();

  Double_t FsigFromFit = fithandler->GetParameter(kFsig);
  Double_t FbFromFit = fithandler->GetParameter(kFb);
  Double_t FsigErr = fithandler->GetParameterError(kFsig);
  Double_t FbErr = fithandler->GetParameterError(kFb);  
 
  printf("fb = %f #pm %f \n",FbFromFit,FbErr);
  printf("fsig = %f #pm %f \n",FsigFromFit,FsigErr);
 
  likely_obj->SetWeightType(weightType[2]/nCandSel,weightType[1]/nCandSel,weightType[0]/nCandSel);
  printf("FF %f - FS %f - SS %f ncand %f \n",weightType[2], weightType[1], weightType[0], nCandSel);
 
  // draw psproper total
  TLegend *leg=new TLegend(0.17,0.72,0.42,0.88);
  if(fitmode != kFitChi2SigMass && fitmode != kFitBkgMass){
  TCanvas *psTotCan = new TCanvas("pseudoProperDecayLength","pseudoProperDecayLength");
  psTotCan->SetLogy();
  Double_t maximum = histpsproper->GetMaximum();
  histpsproper->GetYaxis()->SetRangeUser(0.5*maximum/histpsproper->GetEntries(), maximum*1.50);
  histpsproper->DrawCopy("E");
  TLatex *tex = 0x0;
  tex = new TLatex(0.54,0.876,Form("%2.1f < M(e^{+}e^{-}) < %2.1f GeV/c^{2}",bandLow,bandUp));
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.035);
  tex->SetLineWidth(2);
  tex->Draw();
  Double_t normTot = ((Double_t)histpsproper->GetEntries())*histpsproper->GetBinWidth(1);
  TF1 *psproperTot = likely_obj->GetEvaluateCDFDecayTimeTotalDistrAllTypes(-1.e+04, 1.e+04,normTot);
  psproperTot->SetLineColor(kBlack);
  TFitResultPtr rPsproper = histpsproper->Fit(psproperTot->GetName(),"S0Q");
   
  //legend
  leg->SetBorderSize(0); leg->SetFillColor(0); leg->SetTextFont(42);
  leg->SetFillStyle(0); leg->SetMargin(0.25); 
  leg->SetEntrySeparation(0.15);
 
  if(fitmode == kExtrFb){
  TLatex *lat=new TLatex;
  lat->SetNDC(kTRUE);
  lat->SetTextColor(1);lat->SetTextFont(42);lat->SetTextSize(.035);
  lat->DrawLatex(0.53, 0.82, Form("#chi^{2}/dof = %4.3f ",(rPsproper->Chi2()/(Double_t)rPsproper->Ndf())));
  lat->DrawLatex(0.53, 0.72, Form("F_{Sig}[%1.2f - %1.2f] = %4.3f #pm %4.3f", bandLow, bandUp, FsigFromFit, FsigErr));
  lat->DrawLatex(0.53, 0.62, Form("F_{B}^{uncorrected} = %4.3f #pm %4.3f", FbFromFit, FbErr)); 
  leg->AddEntry(psproperTot, "fit, all","l");
  
  //prompt jpsi
  Double_t normPrompt = (rPsproper->Parameter(0))*FsigFromFit*(1 - FbFromFit);
  TF1 *prompt = likely_obj->GetResolutionFuncAllTypes(-1.e+04,1.e+04,normPrompt);
  prompt->SetLineColor(2);
  prompt->Draw("same");
  leg->AddEntry(prompt, "fit, prompt J/#psi","l");
  
  Double_t normSec = (rPsproper->Parameter(0))*FsigFromFit*FbFromFit;
  TF1 *templateMC = likely_obj->GetFunBAllTypes(-1.e+04,1.e+04,normSec);
  templateMC->SetLineColor(6);
  templateMC->SetFillColor(6);
  templateMC->SetFillStyle(3005);
  templateMC->Draw("same");
  leg->AddEntry(templateMC, "fit, secondary J/#psi","l");
  
  Double_t normBkg =  (rPsproper->Parameter(0))*(1 - FsigFromFit);
  TF1 *psProperBack = (TF1*)(likely_obj->GetEvaluateCDFDecayTimeBkgDistrAllTypes(-1.e+04,1.e+04,normBkg));
  psProperBack->SetName(Form("backFuncTotal_pt%1.1f_%1.1f_mass%1.1f_%1.1f",ptMin,ptMax,bandLow,bandUp));
  psProperBack->SetLineColor(3);
  psProperBack->SetLineWidth(2);
  psProperBack->Draw("same");
  leg->AddEntry(psProperBack, "fit, bkg","l");
  }
 
  if(fitmode == kFitXBkg || fitmode == kFitXBkgMVA || fitmode == kFitXBkgMVAChi2){
  //bkg
  // Double_t normBkg = (rPsproper->Parameter(0))*(1 - FsigFromFit);
  Double_t normBkg =  histpsproper->GetEntries()*histpsproper->GetBinWidth(1);
  likely_obj->SetFMinus(fithandler->GetParameter(kWNegL));
  likely_obj->SetResWeight(fithandler->GetParameter(kWResolution));
  likely_obj->SetFPlus(fithandler->GetParameter(kWPosL));
  likely_obj->SetFSym(fithandler->GetParameter(kWSymL));
  likely_obj->SetFSym1(fithandler->GetParameter(kWSym1L));
  likely_obj->SetLamPlus(fithandler->GetParameter(kLPos));
  likely_obj->SetLamMinus(fithandler->GetParameter(kLNeg));
  likely_obj->SetLamSym(fithandler->GetParameter(kLSym));
  likely_obj->SetLamSym1(fithandler->GetParameter(kLSym1));
  TF1 *psProperBack = likely_obj->GetEvaluateCDFDecayTimeBkgDistrAllTypes(-1.e+04,1.e+04,normBkg);
  TFitResultPtr rPsproperBack = histpsproper->Fit(psProperBack->GetName(),"S0Q");
  tex->DrawLatex(0.53, 0.82, Form("#chi^{2}/dof = %4.3f ",(rPsproperBack->Chi2()/(Double_t)rPsproperBack->Ndf())));  

  std::ofstream ofsbkg;
  TString stringmass=Form("_mass%1.2f_%1.2f",bandLow,bandUp);
  //for(int im=0; im<kMassBins+1; im++) stringmass.Append(Form("_m%d%1.2f",im+1,massBins[im]));
  TString resTypeS = resType; resTypeS.ReplaceAll(";","_");
  stringmass.Append(Form("_type%s",resTypeS.Data()));
  ofsbkg.open (Form("xbkgChi2_Pt%1.2f_%1.2f%s.txt",kPtBins > 1 ? ptEdges[0]: ptMin, kPtBins > 1 ? ptEdges[kPtBins] : ptMax,stringmass.Data()), std::ofstream::out);
  ofsbkg << Form("%f \n",(rPsproperBack->Chi2()/(Double_t)rPsproperBack->Ndf()));
  ofsbkg.close();

  Double_t normBkgPos = fithandler->GetParameter(kWPosL)/(fithandler->GetParameter(kWResolution)+fithandler->GetParameter(kWPosL)+fithandler->GetParameter(kWNegL)+fithandler->GetParameter(kWSymL)+fithandler->GetParameter(kWSym1L))*normBkg;
  likely_obj->SetResWeight(0.);  
  likely_obj->SetFMinus(0.);  
  likely_obj->SetFSym(0.);  
  likely_obj->SetFSym1(0.);  
  TF1 *psProperBackPos = likely_obj->GetEvaluateCDFDecayTimeBkgDistrAllTypes(-1.e+04,1.e+04,normBkgPos);
  psProperBackPos->SetName(Form("backFuncPos_pt%1.1f_%1.1f_mass%1.1f_%1.1f",ptMin,ptMax,bandLow,bandUp));
  psProperBackPos->SetLineColor(4);
  psProperBackPos->SetLineWidth(2);
  psProperBackPos->Draw("same");

  Double_t normBkgNeg = fithandler->GetParameter(kWNegL)/(fithandler->GetParameter(kWResolution)+fithandler->GetParameter(kWPosL)+fithandler->GetParameter(kWNegL)+fithandler->GetParameter(kWSymL)+fithandler->GetParameter(kWSym1L))*normBkg;
  likely_obj->SetFMinus(fithandler->GetParameter(kWNegL));  
  likely_obj->SetResWeight(0.);
  likely_obj->SetFPlus(0.);
  likely_obj->SetFSym(0.);
  likely_obj->SetFSym1(0.);  
  TF1 *psProperBackNeg = likely_obj->GetEvaluateCDFDecayTimeBkgDistrAllTypes(-1.e+04,1.e+04,normBkgNeg);
  psProperBackNeg->SetName(Form("backFuncNeg_pt%1.1f_%1.1f_mass%1.1f_%1.1f",ptMin,ptMax,bandLow,bandUp));
  psProperBackNeg->SetLineColor(kMagenta);
  psProperBackNeg->SetLineWidth(2);
  psProperBackNeg->Draw("same");

  Double_t normBkgSym = fithandler->GetParameter(kWSymL)/(fithandler->GetParameter(kWResolution)+fithandler->GetParameter(kWPosL)+fithandler->GetParameter(kWNegL)+fithandler->GetParameter(kWSymL)+fithandler->GetParameter(kWSym1L))*normBkg;
  likely_obj->SetFSym(fithandler->GetParameter(kWSymL));  
  likely_obj->SetResWeight(0.);
  likely_obj->SetFPlus(0.);
  likely_obj->SetFMinus(0.);
  likely_obj->SetFSym1(0.);  
  TF1 *psProperBackSym = likely_obj->GetEvaluateCDFDecayTimeBkgDistrAllTypes(-1.e+04,1.e+04,normBkgSym);
  psProperBackSym->SetName(Form("backFuncSym_pt%1.1f_%1.1f_mass%1.1f_%1.1f",ptMin,ptMax,bandLow,bandUp));
  psProperBackSym->SetLineColor(kOrange);
  psProperBackSym->SetLineWidth(2);
  psProperBackSym->Draw("same");

  Double_t normBkgSym1 = fithandler->GetParameter(kWSym1L)/(fithandler->GetParameter(kWResolution)+fithandler->GetParameter(kWPosL)+fithandler->GetParameter(kWNegL)+fithandler->GetParameter(kWSymL)+fithandler->GetParameter(kWSym1L))*normBkg;
  likely_obj->SetFSym1(fithandler->GetParameter(kWSym1L));
  likely_obj->SetResWeight(0.);
  likely_obj->SetFPlus(0.);
  likely_obj->SetFMinus(0.);
  likely_obj->SetFSym(0.);
  TF1 *psProperBackSym1 = likely_obj->GetEvaluateCDFDecayTimeBkgDistrAllTypes(-1.e+04,1.e+04,normBkgSym1);
  psProperBackSym1->SetName(Form("backFuncSym_pt%1.1f_%1.1f_mass%1.1f_%1.1f",ptMin,ptMax,bandLow,bandUp));
  psProperBackSym1->SetLineColor(kGray);
  psProperBackSym1->SetLineWidth(2);
  psProperBackSym1->Draw("same");

  Double_t normBkgResolution = fithandler->GetParameter(kWResolution)/(fithandler->GetParameter(kWResolution)+fithandler->GetParameter(kWPosL)+fithandler->GetParameter(kWNegL)+fithandler->GetParameter(kWSymL)+fithandler->GetParameter(kWSym1L))*normBkg;
  likely_obj->SetResWeight(fithandler->GetParameter(kWResolution));
  likely_obj->SetFPlus(0.);
  likely_obj->SetFMinus(0.);
  likely_obj->SetFSym(0.);
  likely_obj->SetFSym1(0.);  
  TF1 *psProperBackResol = likely_obj->GetEvaluateCDFDecayTimeBkgDistrAllTypes(-1.e+04,1.e+04,normBkgResolution);
  psProperBackResol->SetName(Form("backFuncResolution_pt%1.1f_%1.1f_mass%1.1f_%1.1f",ptMin,ptMax,bandLow,bandUp));
  psProperBackResol->SetLineColor(kRed);
  psProperBackResol->SetLineWidth(2);
  psProperBackResol->Draw("same");
  
  // }
  psProperBack->SetName(Form("backFuncTotal_pt%1.1f_%1.1f_mass%1.1f_%1.1f",ptMin,ptMax,bandLow,bandUp));
  psProperBack->SetLineColor(3);
  psProperBack->SetLineWidth(2);
  psProperBack->Draw("same");
  leg->AddEntry(psProperBack, "fit, bkg","l");
        
  	// TString resTypeS = resType; resTypeS.ReplaceAll(";","_"); 
	if(fitmode == kFitXBkgMVA || fitmode == kFitXBkgMVAChi2) SaveParameters(Form("inputFiles_%1.1f_%1.1f/XBkgParameters_mass%1.1f_%1.1f_resType%s_pt%1.1f_%1.1f",ptEdges[0],ptEdges[kPtBins],bandLow,bandUp,resTypeS.Data(),ptMin,ptMax),psProperBkgParam,sizeof(psProperBkgParam)/sizeof(Int_t));
  	else SaveParameters(Form("inputFiles_%1.1f_%1.1f/XBkgParameters",ptMin,ptMax),psProperBkgParam,sizeof(psProperBkgParam)/sizeof(Int_t));
  	psTotCan->SaveAs(Form("inputFiles_%1.1f_%1.1f/XBkgProjection_mass%1.1f_%1.1f_types%s_pt%1.1f_%1.1f.pdf",kPtBins > 1 ? ptEdges[0]: ptMin, kPtBins > 1 ? ptEdges[kPtBins] : ptMax,bandLow,bandUp,resTypeS.Data(),ptMin,ptMax));
  	psTotCan->SaveAs(Form("inputFiles_%1.1f_%1.1f/XBkgProjection_mass%1.1f_%1.1f_types%s_pt%1.1f_%1.1f.root",kPtBins > 1 ? ptEdges[0]: ptMin, kPtBins > 1 ? ptEdges[kPtBins] : ptMax,bandLow,bandUp,resTypeS.Data(),ptMin,ptMax));
  	// save bkg function 
        if(fitmode == kFitXBkgMVA || fitmode == kFitXBkgMVAChi2) SaveFunctions(resType, ptMin, ptMax, bandLow, bandUp);  // save functions for likelihood fit
     }
   
  leg->Draw("same");
  if(fitmode == kExtrFb) { psproperTot->Draw("same");  psTotCan->SaveAs(Form("inputFiles_%1.1f_%1.1f/LikelihoodFitProjWholeMass.root",ptMin,ptMax)); }
  }

  Double_t FsigNew=0.; Double_t FsigNewErr=0.; 
  if(fitmode == kFitBkgMass || fitmode == kExtrFb){
  // draw invariant mass 
  Double_t normMass =((Double_t)histMass->GetEntries())*histMass->GetBinWidth(1); 
  TF1 *invMassFunc = likely_obj->GetEvaluateCDFInvMassTotalDistr(bandLow,bandUp,normMass);
  invMassFunc->SetLineColor(kBlack);
  TCanvas *invMassCan = new TCanvas("invMassCanvas","invMassCanvas");
  invMassCan->cd();
  TFitResultPtr rMass = histMass->Fit(invMassFunc->GetName(),"S0Q");
  histMass->DrawCopy("E"); 
  invMassFunc->SetLineColor(1);
  Double_t intTot = invMassFunc->Integral(bandLow,bandUp); 
  TLatex *lat=new TLatex;
  lat->SetNDC(kTRUE);
  lat->SetTextColor(1);lat->SetTextFont(42);lat->SetTextSize(.035);
  if(fitmode == kFitBkgMass || fitmode == kExtrFb) lat->DrawLatex(0.53, 0.82, Form("#chi^{2}/dof = %4.2f ",(rMass->Chi2()/(Double_t)rMass->Ndf())));
  invMassFunc->SetChisquare(rMass->Ndf());
  invMassFunc->SetNDF(rMass->Ndf());
  if(fitmode == kFitBkgMass || fitmode == kExtrFb) lat->DrawLatex(0.53, 0.72, Form("F_{Sig}[%1.2f - %1.2f] = %4.3f #pm %4.3f", bandLow, bandUp, FsigFromFit, FsigErr));
  TLegend *legMass=new TLegend(0.17,0.72,0.42,0.88);
  legMass->SetBorderSize(0); legMass->SetFillColor(0); legMass->SetTextFont(42);
  legMass->SetFillStyle(0); legMass->SetMargin(0.25); 
  legMass->SetEntrySeparation(0.15);
  TF1 *invMassSig = likely_obj->GetEvaluateCDFInvMassSigDistr(bandLow,bandUp,intTot*FsigFromFit);
  TF1 *invMassBkg = likely_obj->GetEvaluateCDFInvMassBkgDistr(bandLow,bandUp,intTot*(1.-FsigFromFit));
  invMassSig->SetLineColor(4);
  invMassBkg->SetLineColor(3);
  invMassSig->Draw("same");
  legMass->AddEntry(invMassSig, "fit, signal","l");
  if(fitmode == kFitBkgMass || fitmode == kExtrFb ){
  invMassBkg->Draw("same");
  invMassFunc->Draw("same");
  legMass->AddEntry(invMassBkg, "fit, background","l");
  legMass->AddEntry(invMassFunc, "fit, all","l");
  }
  legMass->Draw("same");
  Double_t integSig = invMassSig->Integral(bandLowSignal,bandUpSignal);
  Double_t integBkg = invMassBkg->Integral(bandLowSignal,bandUpSignal);
  FsigNew = integSig/(integSig+integBkg);
  FsigNewErr = (FsigNew*FsigErr)/FsigFromFit;
  if(fitmode == kFitBkgMass) SaveParameters(Form("inputFiles_%1.1f_%1.1f/InvMassBkgExpoLikelihoodFit",ptMin,ptMax),invMassBkgExpoParam,sizeof(invMassBkgExpoParam)/sizeof(Int_t));
  if(fitmode == kExtrFb) invMassCan->SaveAs(Form("inputFiles_%1.1f_%1.1f/LikelihoodFitInvariantMass.root",ptMin,ptMax));
  }
   // save all parameters -> can be re-used as input 
   if(fitmode == kExtrFb) SaveParameters(Form("inputFiles_%1.1f_%1.1f/startingParametersLast",ptMin,ptMax),allParam,sizeof(allParam)/sizeof(Int_t)); 

   if(fitmode == kExtrFb){
   // draw pseudo-proper decay length projections under the signal region
   AliDielectronBtoJPSItoEleCDFfitFCN *likely_obj_proj = 0x0;
   AliDielectronBtoJPSItoEle* aBtoJPSItoEle_proj =new AliDielectronBtoJPSItoEle();
   aBtoJPSItoEle_proj->SetResTypeAnalysis(resType);
   for(int j =0; j < kNumPar; j++) {paramInputValues[j] = fithandler->GetParameter(j);}
   paramInputValues[kFsig] = FsigNew;
   Double_t* xx=0x0; Double_t* mm=0x0; Double_t *ppt=0; Int_t*tt=0; Int_t nn=0;
   aBtoJPSItoEle_proj->ReadCandidates(nt,xx,mm,ppt,tt,nn,bandLowSignal,bandUpSignal,ptMin,ptMax);
   aBtoJPSItoEle_proj->SetFitHandler(xx,mm,ppt,tt,nn);
   aBtoJPSItoEle_proj->CloneMCtemplate(hCsiMCEvtGen);  
   aBtoJPSItoEle_proj->SetCsiMC(); 

   AliDielectronBtoJPSItoEleCDFfitHandler* fithandler_proj = aBtoJPSItoEle_proj->GetCDFFitHandler();
   fithandler_proj->SetPrintStatus(kTRUE);
   fithandler_proj->SetParamStartValues(paramInputValues);
   fithandler_proj->SetCrystalBallFunction(kTRUE);
   fithandler_proj->SetMassWndLow((TDatabasePDG::Instance()->GetParticle(443)->Mass()) - bandLowSignal);
   fithandler_proj->SetMassWndHigh(bandUpSignal - (TDatabasePDG::Instance()->GetParticle(443)->Mass()));
   likely_obj_proj = fithandler_proj->LikelihoodPointer();
   likely_obj_proj->SetAllParameters(paramInputValues);
   likely_obj_proj->SetWeightType(weightTypeSignal[2]/nCandSelSignal,weightTypeSignal[1]/nCandSelSignal,weightTypeSignal[0]/nCandSelSignal);
   likely_obj_proj->ComputeMassIntegral();
   // likely_obj_proj->PrintStatus();
   // 
   TCanvas *psSignalCan = new TCanvas("pseudoProperDecayLengthSignal","pseudoProperDecayLengthSignal");
   psSignalCan->SetLogy();
   Double_t maximum = histpsproperSignal->GetMaximum();
   histpsproperSignal->GetYaxis()->SetRangeUser(0.5*maximum/histpsproperSignal->GetEntries(), maximum*1.50);
   histpsproperSignal->GetXaxis()->SetRangeUser(-3000.,3000.);
   histpsproperSignal->DrawCopy("E");
   TLatex *tex = 0x0;
   tex = new TLatex(0.54,0.876,Form("%2.1f < M(e^{+}e^{-}) < %2.1f GeV/c^{2}",bandLowSignal,bandUpSignal));
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();
 
   Double_t normPsSignal = ((Double_t)histpsproperSignal->GetEntries())*histpsproperSignal->GetBinWidth(1);
   TF1 *psproperTot_sig = likely_obj_proj->GetEvaluateCDFDecayTimeTotalDistrAllTypes(-1.e+04, 1.e+04, normPsSignal);
   psproperTot_sig->SetLineColor(kBlack);
   psproperTot_sig->SetName("psProperTotal_sig");
   TFitResultPtr rPsproper_sig = histpsproperSignal->Fit(psproperTot_sig->GetName(),"S0Q");
   tex->SetNDC(kTRUE);
   tex->SetTextColor(1);tex->SetTextFont(42);tex->SetTextSize(.035);
   tex->DrawLatex(0.53, 0.82, Form("#chi^{2}/dof = %4.3f ",(rPsproper_sig->Chi2()/(Double_t)rPsproper_sig->Ndf())));
   tex->DrawLatex(0.53, 0.72, Form("F_{Sig}(scaled)[%1.2f-%1.2f] = %4.3f #pm %4.3f", bandLowSignal, bandUpSignal, FsigNew,FsigNewErr));
   tex->DrawLatex(0.53, 0.62, Form("F_{B}^{uncorrected} = %4.3f #pm %4.3f", FbFromFit, FbErr));
   ///
   Double_t normPromptSig = (normPsSignal)*FsigNew*(1 - FbFromFit);
   TF1 *prompt_sig = likely_obj_proj->GetResolutionFuncAllTypes(-1.e+04,1.e+04,normPromptSig);
   prompt_sig->SetLineColor(2);
   prompt_sig->Draw("same");
   
   Double_t normSecSig = (normPsSignal)*FsigNew*FbFromFit;
   TF1 *templateMC_sig = likely_obj_proj->GetFunBAllTypes(-1.e+04,1.e+04,normSecSig);
   templateMC_sig->SetLineColor(6);
   templateMC_sig->SetFillColor(6);
   templateMC_sig->SetFillStyle(3005);
   templateMC_sig->Draw("same");

   Double_t normBkgSig = (normPsSignal)*(1 - FsigNew);
   TF1 *psProperBack_sig = likely_obj_proj->GetEvaluateCDFDecayTimeBkgDistrAllTypes(-1.e+04,1.e+04,normBkgSig);
   psProperBack_sig->SetLineColor(3);
   psProperBack_sig->Draw("same");
   psproperTot_sig->Draw("same");
   leg->Draw("same");
   psSignalCan->SaveAs(Form("inputFiles_%1.1f_%1.1f/LikelihoodFitProjSignalRegionMass.root",ptMin,ptMax)); 
   }

return;  
}


void FitCDFLikelihoodMVA(TString resType, Double_t ptMin, Double_t ptMax, Double_t bandLow, Double_t bandUp, Bool_t kQuadratic,  Double_t centmin, Double_t centmax, Bool_t fixFsig, Int_t mBkgOpt){
///////////////////////////////////////////////////////////////////
// Set inv mass / pt limits
  //inv mass functions normalized to  1 between bandLow - bandUp 
  //
  // projection of pseudoproper decay length functions drawn in the signal region 
  Double_t bandLowSignal=2.92;
  Double_t bandUpSignal=3.16; 
  Double_t fsigFrom1D = -1.;
  TString resTypeS = resType; resTypeS.ReplaceAll(";","_");
  if(kPtBins == 1) {ptEdges[0] = ptMin; ptEdges[1] = ptMax;}
  ifstream filenameEdges(Form("inputFiles_%1.1f_%1.1f/massEdges.txt",ptMin,ptMax));
  if(filenameEdges) { for(int im=0; im<kMassBins+1; im++) { filenameEdges >> massBins[im]; /*printf("mass %f \n",massBins[im]);*/} }
  if(fixFsig){
      ifstream filenameFsig;
      if(mBkgOpt == 1 ) filenameFsig.open(Form("inputFiles_%1.1f_%1.1f/FractionOfSignalMassME%s_.txt",ptMin,ptMax,resTypeS.Data()));    
      else if(mBkgOpt == 2 ) filenameFsig.open(Form("inputFiles_%1.1f_%1.1f/FractionOfSignalMassLS%s_.txt",ptMin,ptMax,resTypeS.Data()));    
      else filenameFsig.open(Form("inputFiles_%1.1f_%1.1f/fSIGPt_%1.2f_%1.2f_type%s.txt",ptMin,ptMax,bandLow,bandUp,resTypeS.Data()));  

  if(filenameFsig) {filenameFsig >> fsigFrom1D; printf("Fsig from 1D in %1.2f-%1.2f GeV/c^{2} -> %1.2f",bandLow,bandUp,fsigFrom1D); }     
}
  // 
  TFile f(inputDistr); // input: data
  TFile ftemplate(Form("inputFiles_%1.1f_%1.1f/XtemplateNonPromptJpsi.root",ptMin,ptMax));
  TH1F *hCsiMCEvtGen = (TH1F*)ftemplate.Get("psTemplate");
  Double_t integral = 0;
  for(int i=1;i<hCsiMCEvtGen->GetNbinsX()+1; i++) integral += (hCsiMCEvtGen->GetBinContent(i)*hCsiMCEvtGen->GetBinWidth(i));
  hCsiMCEvtGen->Scale(1./integral);
  
  Double_t* x=0x0; Double_t* m=0x0; Double_t *pt =0x0; Int_t* type=0; Int_t n=0;
  AliDielectronBtoJPSItoEle* aBtoJPSItoEle =new AliDielectronBtoJPSItoEle();

  TNtuple *ntOrig=0x0;
  ntOrig=(TNtuple*)f.Get("fNtupleJPSI");
  TNtuple *nt=new TNtuple("ntuplaSigna_new","NtuplaSignal","Xdecaytime:Mass:Type:Pt",10000000);
  AddNtupla(nt,ntOrig,ptMin,ptMax,bandLow,bandUp,centmin,centmax);  // add candidates in [bandLow, bandUp]

  // set all starting parameters 
  SetStartingParameters(Form("inputFiles_%1.1f_%1.1f/startingParameters.root",ptMin,ptMax));
  // set inv mass signal parameters 
  SetStartingParameters(Form("inputFiles_%1.1f_%1.1f/InvariantMassSignalMC.root",ptMin,ptMax));
  // set inv mass signal parameters
  //SetStartingParameters(Form("inputFiles_%1.1f_%1.1f/InvariantMassBkgExpoChi2Fit.root",ptMin,ptMax));
  if(mBkgOpt == 1 ) SetStartingParameters(Form("inputFiles_%1.1f_%1.1f/InvariantMassBkgME%s_.root",ptMin,ptMax,resTypeS.Data()));
  else if(mBkgOpt == 2 ) SetStartingParameters(Form("inputFiles_%1.1f_%1.1f/InvariantMassBkgLS%s_.root",ptMin,ptMax,resTypeS.Data()));
   else SetStartingParameters(Form("inputFiles_%1.1f_%1.1f/InvariantMassBkgExpoChi2Fit.root",ptMin,ptMax));  
  
  
  aBtoJPSItoEle->SetResTypeAnalysis(resType);
  aBtoJPSItoEle->ReadCandidates(nt,x,m,pt,type,n,bandLow,bandUp); // read N-Tuples
  printf("+++\n+++ Number of total candidates (prim J/psi + secondary J/psi + bkg) ---> %d candidates \n+++\n",n);

  aBtoJPSItoEle->SetFitHandler(x,m,pt,type,n); // Set the fit handler with given values of x, m, # of candidates 

  aBtoJPSItoEle->CloneMCtemplate(hCsiMCEvtGen);    // clone MC template and copy internally
                                                   // in this way any model can be setted from outside
  aBtoJPSItoEle->SetCsiMC(); // Pass the MC template to the CDF fit function

  fithandler = aBtoJPSItoEle->GetCDFFitHandler(); // Get the fit handler
  //
  // Set some fit options through the handler class
  //
  fithandler->SetPrintStatus(kTRUE);
  fithandler->SetCrystalBallFunction(kTRUE);
  fithandler->SetExponentialFunction(useExpoBkgMass); // kFALSE - polynomial function used

  Double_t massLow = (TDatabasePDG::Instance()->GetParticle(443)->Mass()) - bandLow;
  Double_t massHigh = bandUp - (TDatabasePDG::Instance()->GetParticle(443)->Mass());
  fithandler->SetMassWndLow(massLow);
  fithandler->SetMassWndHigh(massHigh);
  // 2D likelihood fit for fB extraction
  paramInputValues[kFb] = 0.10; paramInputValues[kFsig] = 0.010;
  if(fixFsig) paramInputValues[kFsig] = fsigFrom1D;
  fithandler->SetParamStartValues(paramInputValues);
  for(int ipar=0; ipar<kNumPar; ipar++) { if(ipar != kFb && ipar !=kFsig) fithandler->FixParam(ipar,kTRUE);  }
  if(fixFsig) fithandler->FixParam(kFsig,kTRUE);
  //for(int ipar=0; ipar<kNumPar; ipar++) { if(ipar != kFb) fithandler->FixParam(ipar,kTRUE);  }
  
  likely_obj = fithandler->LikelihoodPointer();
  likely_obj->SetAllParameters(paramInputValues);

  likely_obj->ComputeMassIntegral();
  likely_obj->PrintStatus();

   // MVA settings 
  likely_obj->SetMultivariateFit(kTRUE);
  likely_obj->InitializeFunctions(kPtBins,kMassBins);

  // masswindows
  TArrayD *massWind =new TArrayD(kMassBins+1); // mass windows to define signal and adiacent regions
  for(int im=0; im<kMassBins+1; im++)  massWind->AddAt(massBins[im],im);
  likely_obj->SetMassWindows(massWind);

  TArrayD *ptWind = new TArrayD(kPtBins+1); // pt windows
  for(int ipt=0; ipt<kPtBins+1; ipt++) ptWind->AddAt(ptEdges[ipt],ipt);
  likely_obj->SetPtWindows(ptWind);

  // set interpolation region
  likely_obj->SetExtrapolationRegion(extrRegion);

  //
  //fill histos from Ntupla
  //Int_t nbinsMass = (Int_t)((bandUp - bandLow)/0.04);  
  Int_t nbinsMass = (Int_t)((bandUp - bandLow)/0.02);  
  TH1F *histMass = new TH1F("histMass","Invariant Mass; InvMass[GeV]; Entries/40MeV",nbinsMass,bandLow,bandUp);
  histMass->SetLineColor(1);
  histMass->SetMarkerColor(1);
  histMass->SetMarkerStyle(20);
  histMass->SetMarkerSize(0.7);

  TH1F *histpsproper = new TH1F("psproper_decay_length",Form("psproper_decay_length_distrib(%1.1f < M < %1.1f GeV/c^{2});pseudoproper decay length [#mum];Entries/40#mum",bandLow,bandUp),300,-6000.,6000.);
  histpsproper->SetLineColor(1);
  histpsproper->SetMarkerColor(1);
  histpsproper->SetMarkerStyle(20);
  histpsproper->SetMarkerSize(0.7);

  TH1F *histpsproperSignal = new TH1F("psproper_decay_length_signal",Form("psproper_decay_length_distrib(%1.2f < M < %1.2f GeV/c^{2});pseudoproper decay length [#mum];Entries/40#mum",bandLowSignal,bandUpSignal),300,-6000.,6000.);
  histpsproperSignal->SetLineColor(1);
  histpsproperSignal->SetMarkerColor(1);
  histpsproperSignal->SetMarkerStyle(20);
  histpsproperSignal->SetMarkerSize(0.7);

  Float_t mass =0.; Float_t psproper = 0.; Float_t typeCand = 0.; Int_t nb = 0;
  Float_t ptCand = 0.;
  TString arrType[]={"SS","FS","FF"};
  nt->SetBranchAddress("Xdecaytime",&psproper);
  nt->SetBranchAddress("Mass",&mass);
  nt->SetBranchAddress("Type",&typeCand);
  nt->SetBranchAddress("Pt",&ptCand);
  Int_t fNcurrent=0; Double_t nCandSel = 0 ; Double_t nCandSelSignal=0.;
  nb = (Int_t)nt->GetEvent(fNcurrent);
  cout << nb << endl;

  for(Int_t iev=0; iev<(nt->GetEntries()); iev++){
   if(resType.Contains(arrType[(Int_t)typeCand]) && mass > bandLow && mass < bandUp){
   	nCandSel += 1;
   	weightType[(Int_t)typeCand] += 1.;
   	histMass->Fill(mass);
   	histpsproper->Fill(psproper);
   	if(mass > bandLowSignal && mass < bandUpSignal) {
   	histpsproperSignal->Fill(psproper);
   	weightTypeSignal[(Int_t)typeCand] += 1.;
   	nCandSelSignal += 1;
   	}
      }
   /////
   //  
   if(resType.Contains(arrType[(Int_t)typeCand])){ 
   for(int ipt=0; ipt<kPtBins; ipt++){
        for(int im=0; im<kMassBins; im++){
         if((ptWind->At(ipt) < ptCand) && (ptCand < ptWind->At(ipt+1)) && ((massWind->At(im) < mass) && (mass < massWind->At(im+1)) ))
          nCandPtMassWnd[im][ipt][(Int_t)typeCand]++;
          }
          if((ptWind->At(ipt) < ptCand) && (ptCand < ptWind->At(ipt+1)) && (mass>bandLowSignal && mass<bandUpSignal)) nCandPtMassWndSignal[ipt][(Int_t)typeCand]++; 
       } 
    }
   fNcurrent++;
   nb = (Int_t) nt->GetEvent(fNcurrent);
   }
 
  //// load pt-dep. resolution parameters
  Double_t ***resolParamAll = LoadResolutionParameters(ptEdges);
  likely_obj->SetResParams(resolParamAll);
 
  // load background functions
  LoadBackgroundFunctionsAndWeights(likely_obj,Form("inputFiles_%1.1f_%1.1f/XBkgFunctions.root",ptMin,ptMax),ntOrig,kQuadratic, kMassBins > 3 ? kTRUE : kFALSE);
  //LoadBackgroundFunctionsAndWeights(likely_obj,Form("inputFiles_%1.1f_%1.1f/XBkgFunctions.root",ptMin,ptMax),ntOrig,kQuadratic, kFALSE);

  // do minimization 
  aBtoJPSItoEle->DoMinimization();
 
  Double_t FsigFromFit = fithandler->GetParameter(kFsig);
  Double_t FbFromFit = fithandler->GetParameter(kFb);
  Double_t FsigErr = fithandler->GetParameterError(kFsig);
  Double_t FbErr = fithandler->GetParameterError(kFb);

  printf("fb = %f #pm %f \n",FbFromFit,FbErr);
  printf("fsig = %f #pm %f \n",FsigFromFit,FsigErr);

  likely_obj->SetWeightType(weightType[2]/nCandSel,weightType[1]/nCandSel,weightType[0]/nCandSel);
  printf("FF %f - FS %f - SS %f ncand %f \n",weightType[2], weightType[1], weightType[0], nCandSel);
/* 
  std::ofstream ofs;
  TString stringmass="";
  for(int im=0; im<kMassBins+1; im++) stringmass.Append(Form("_m%d%1.2f",im+1,massBins[im]));  
  TString resTypeS = resType; resTypeS.ReplaceAll(";","_");
  stringmass.Append(Form("_type%s",resTypeS.Data()));
  if(kQuadratic) stringmass.Append("_quadratic");
  if(fixFsig) stringmass.Append("_FsigFixed");
  ofs.open (Form("fbPt%1.2f_%1.2f%s.txt",ptMin,ptMax,stringmass.Data()), std::ofstream::out | std::ofstream::app);
  ofs << Form("%f %f %f %f \n",FbFromFit,FbErr,FsigFromFit,FsigErr);
  */
  ///
  // draw psproper total
  /*TCanvas *dPsProp = new TCanvas("pseudoproperDL","pseudoproperDL");
  dPsProp->SetLogy();
  Double_t norm1 = ((Double_t)histpsproper->GetEntries())*histpsproper->GetBinWidth(1);

  TF1 *psproperTot = new TF1("averageTotal",GetAverageFunctions,-8000.,8000.,3);
  psproperTot->SetParameter(0,1.);
  psproperTot->FixParameter(1,FsigFromFit);
  psproperTot->FixParameter(2,FbFromFit);
  psproperTot->SetNpx(2500);
  psproperTot->SetLineColor(kBlack);
  
  TFitResultPtr rPsproper = histpsproper->Fit(psproperTot->GetName(),"S0");
  TLatex *lat=new TLatex;
  lat->SetNDC(kTRUE);
  lat->SetTextColor(1);lat->SetTextFont(42);lat->SetTextSize(.035);

  histpsproper->SetLineColor(1);
  histpsproper->SetMarkerColor(1);
  histpsproper->SetMarkerStyle(20);
  histpsproper->SetMarkerSize(0.7);
  Double_t maximum = histpsproper->GetMaximum();
  histpsproper->GetYaxis()->SetRangeUser(0.5*maximum/histpsproperSignal->GetEntries(), maximum*1.50);
  histpsproper->DrawCopy("E");
  lat->DrawLatex(0.53, 0.82, Form("#chi^{2}/dof = %4.3f ",(rPsproper->Chi2()/(Double_t)rPsproper->Ndf())));
  lat->DrawLatex(0.53, 0.72, Form("F_{Sig}[%1.2f - %1.2f] = %4.3f #pm %4.3f", bandLow, bandUp, FsigFromFit, FsigErr));
  lat->DrawLatex(0.53, 0.62, Form("F_{B}^{uncorrected} = %4.3f #pm %4.3f", FbFromFit, FbErr));

   TLatex *tex = 0x0;
   tex = new TLatex(0.54,0.876,Form("%2.1f < M(e^{+}e^{-}) < %2.1f GeV/c^{2}",bandLow,bandUp));
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();
  
   //prompt jpsi
  Double_t normPrompt = (rPsproper->Parameter(0))*FsigFromFit*(1 - FbFromFit);
  TF1 *prompt = new TF1("averageResolution",GetAverageFunctions,-3000.,3000.,3);
  prompt->SetParameter(0,normPrompt);
  prompt->SetParameter(1,1.);
  prompt->SetParameter(2,0.);
  prompt->SetNpx(500);

  prompt->SetLineColor(2);
  prompt->Draw("same");

  Double_t normSec = (rPsproper->Parameter(0))*FsigFromFit*FbFromFit;
  TF1 *templateMC = new TF1("averageNonPrompt",GetAverageFunctions,-5000.,5000.,3);
  templateMC->SetParameter(0,normSec);
  templateMC->SetParameter(1,1.);
  templateMC->SetParameter(2,1.);
  templateMC->SetNpx(500);

  templateMC->SetLineColor(6);
  templateMC->SetFillColor(6);
  templateMC->SetFillStyle(3005);
  templateMC->Draw("same");

  //bkg
  Double_t normBkg = (rPsproper->Parameter(0))*(1 - FsigFromFit);
  TF1 *psProperBack= new TF1("averageBkg",GetAverageFunctions,-8000.,8000.,3);
  psProperBack->FixParameter(0,normBkg);
  psProperBack->FixParameter(1,0.);
  psProperBack->FixParameter(2,0.);
  psProperBack->SetNpx(500);
  psProperBack->SetLineColor(3);
  psProperBack->SetLineWidth(2);
  psProperBack->Draw("same");
  psproperTot->SetLineColor(kBlack);
  psproperTot->Draw("same");
  
  TLegend *leg=new TLegend(0.17,0.72,0.42,0.88);
  leg->SetBorderSize(0); leg->SetFillColor(0); leg->SetTextFont(42);
  leg->SetFillStyle(0); leg->SetMargin(0.25); 
  leg->SetEntrySeparation(0.15);
  leg->AddEntry(psproperTot, "fit, all","l");
  leg->AddEntry(prompt, "fit, prompt J/#psi","l");
  leg->AddEntry(templateMC, "fit, secondary J/#psi","l");
  leg->AddEntry(psProperBack, "fit, bkg","l");
  leg->Draw("same");

  psproperTot->Draw("same");
  dPsProp->SaveAs(Form("inputFiles_%1.1f_%1.1f/LikelihoodFitProjWholeMass.root",ptMin,ptMax));
  */
  // draw invariant mass 
  Double_t normMass =((Double_t)histMass->GetEntries())*histMass->GetBinWidth(1);
  TF1 *invMassFunc = likely_obj->GetEvaluateCDFInvMassTotalDistr(bandLow,bandUp,normMass);
  invMassFunc->SetLineColor(kBlack);
  TCanvas *invMassCan = new TCanvas("invMassCanvas","invMassCanvas");
  invMassCan->cd();
  TFitResultPtr rMass = histMass->Fit(invMassFunc->GetName(),"S0Q");
  histMass->SetTitle(" ");
  histMass->GetYaxis()->SetRangeUser(0.,histMass->GetMaximum()*1.4);
  histMass->DrawCopy("E");
  invMassFunc->SetLineColor(1);
  Double_t intTot = invMassFunc->Integral(bandLow,bandUp,1.e-06);
  TLatex *lat=new TLatex;
  lat->SetNDC(kTRUE);
  lat->SetTextColor(1);lat->SetTextFont(42);lat->SetTextSize(.035);
  lat->DrawLatex(0.53, 0.82, Form("#chi^{2}/dof = %4.2f ",(rMass->Chi2()/(Double_t)rMass->Ndf())));
  invMassFunc->SetChisquare(rMass->Ndf());
  invMassFunc->SetNDF(rMass->Ndf());
  lat->DrawLatex(0.53, 0.72, Form("F_{Sig}[%1.2f - %1.2f] = %4.3f #pm %4.3f", bandLow, bandUp, FsigFromFit, FsigErr));
  TLegend *legMass=new TLegend(0.17,0.72,0.42,0.88);
  legMass->SetBorderSize(0); legMass->SetFillColor(0); legMass->SetTextFont(42);
  legMass->SetFillStyle(0); legMass->SetMargin(0.25); 
  legMass->SetEntrySeparation(0.15);
  TF1 *invMassSig = likely_obj->GetEvaluateCDFInvMassSigDistr(bandLow,bandUp,intTot*FsigFromFit);
  TF1 *invMassBkg = likely_obj->GetEvaluateCDFInvMassBkgDistr(bandLow,bandUp,intTot*(1.-FsigFromFit));
  invMassSig->SetLineColor(4);
  invMassBkg->SetLineColor(3);
  invMassSig->Draw("same");
  legMass->AddEntry(invMassSig, "fit, signal","l");
  invMassBkg->Draw("same");
  invMassFunc->Draw("same");
  legMass->AddEntry(invMassBkg, "fit, background","l");
  legMass->AddEntry(invMassFunc, "fit, all","l");
  legMass->Draw("same"); 
  TString resTypeSlikeFit = resType; resTypeSlikeFit.ReplaceAll(";","_");
  if(kQuadratic) resTypeSlikeFit.Append("_quadraticWeights");
  if(fixFsig) resTypeSlikeFit.Append("_fSigFixed");
  if(mBkgOpt == 1 ) resTypeSlikeFit.Append("_ME");
  else if(mBkgOpt == 2 ) resTypeSlikeFit.Append("_LS");
  invMassCan->SaveAs(Form("inputFiles_%1.1f_%1.1f/LikelihoodFitInvariantMass%s.root",ptMin,ptMax,resTypeSlikeFit.Data()));
  invMassCan->SaveAs(Form("inputFiles_%1.1f_%1.1f/LikelihoodFitInvariantMass%s.pdf",ptMin,ptMax,resTypeSlikeFit.Data()));
  Double_t integSig = invMassSig->Integral(bandLowSignal,bandUpSignal,1.e-06);
  Double_t integBkg = invMassBkg->Integral(bandLowSignal,bandUpSignal,1.e-06);
  Double_t FsigNew = integSig/(integSig+integBkg);
  Double_t FsigNewErr = (FsigNew*FsigErr)/FsigFromFit;
  //// 
 
  // draw pseudo-proper decay length projections under the signal region
   TCanvas *psSignalCan = new TCanvas("pseudoProperDecayLengthSignal","pseudoProperDecayLengthSignal");
   psSignalCan->SetLogy();
   histpsproperSignal->SetTitle(" ");
   Double_t maximum = histpsproperSignal->GetMaximum();
   histpsproperSignal->GetYaxis()->SetRangeUser(0.5*maximum/histpsproperSignal->GetEntries(), maximum*1.50);
   histpsproperSignal->GetXaxis()->SetRangeUser(-3000.,3000.);
   histpsproperSignal->DrawCopy("E");
   TLatex *tex = new TLatex(0.54,0.876,Form("%2.1f < M(e^{+}e^{-}) < %2.1f GeV/c^{2}",bandLowSignal,bandUpSignal));
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();

   Double_t normPsSignal = ((Double_t)histpsproperSignal->GetEntries())*histpsproperSignal->GetBinWidth(1);
   TF1 *psproperTot_sig = new TF1("averageTotalSig",GetAverageFunctionsSignal,-8000.,8000.,3);
   psproperTot_sig->SetParameter(0,normPsSignal);
   psproperTot_sig->FixParameter(1,FsigNew);
   psproperTot_sig->FixParameter(2,FbFromFit);
   //psproperTot_sig->SetNpx(2500);
   psproperTot_sig->SetNpx(1000);
   psproperTot_sig->SetLineColor(kBlack);
   psproperTot_sig->SetName("psProperTotal_sig");
   TFitResultPtr rPsproper_sig = histpsproperSignal->Fit(psproperTot_sig->GetName(),"S0Q");
   //TFitResultPtr rPsproper_sig = histpsproperSignal->Fit(psproperTot_sig->GetName(),"S");
   tex->SetNDC(kTRUE);
   tex->SetTextColor(1);tex->SetTextFont(42);tex->SetTextSize(.035);
   tex->DrawLatex(0.53, 0.82, Form("#chi^{2}/dof = %4.3f ",(rPsproper_sig->Chi2()/(Double_t)rPsproper_sig->Ndf())));
   tex->DrawLatex(0.53, 0.72, Form("F_{Sig}(scaled)[%1.2f-%1.2f] = %4.3f #pm %4.3f", bandLowSignal, bandUpSignal, FsigNew,FsigNewErr));
   tex->DrawLatex(0.53, 0.62, Form("F_{B}^{uncorrected} = %4.3f #pm %4.3f", FbFromFit, FbErr));
   //
   std::ofstream ofs;
   TString stringmass="";
   for(int im=0; im<kMassBins+1; im++) stringmass.Append(Form("_m%d%1.2f",im+1,massBins[im]));
   //TString resTypeS = resType; resTypeS.ReplaceAll(";","_");
   stringmass.Append(Form("_type%s",resTypeS.Data()));
   if(kQuadratic) stringmass.Append("_quadratic");
   if(fixFsig) stringmass.Append("_FsigFixed");
   if(mBkgOpt == 1 ) stringmass.Append("_ME");
   else if(mBkgOpt == 2 ) stringmass.Append("_LS");
   ofs.open (Form("fbPt%1.2f_%1.2f%s.txt",ptMin,ptMax,stringmass.Data()), std::ofstream::out | std::ofstream::app);
   ofs << Form("%f %f %f %f %f %f \n",FbFromFit,FbErr,FsigFromFit,FsigErr,(rPsproper_sig->Chi2()/(Double_t)rPsproper_sig->Ndf()),(rMass->Chi2()/(Double_t)rMass->Ndf()));
   ofs.close();

   ///
   //prompt jpsi
   Double_t normPromptSig = (rPsproper_sig->Parameter(0))*FsigNew*(1 - FbFromFit);
   TF1 *promptSig = new TF1("averageResolutionSignal",GetAverageFunctionsSignal,-3000.,3000.,3);
   promptSig->SetParameter(0,normPromptSig);
   promptSig->SetParameter(1,1.);
   promptSig->SetParameter(2,0.);
   promptSig->SetNpx(500);
   promptSig->SetLineColor(2);
   promptSig->Draw("same");
 
   //non prompt J/psi
   Double_t normSecSig = (rPsproper_sig->Parameter(0))*FsigNew*FbFromFit;
   TF1 *templateMCSig = new TF1("averageNonPromptSignal",GetAverageFunctionsSignal,-5000.,5000.,3);
   templateMCSig->SetParameter(0,normSecSig);
   templateMCSig->SetParameter(1,1.);
   templateMCSig->SetParameter(2,1.);
   templateMCSig->SetNpx(500);
   templateMCSig->SetLineColor(6);
   templateMCSig->SetFillColor(6);
   templateMCSig->SetFillStyle(3005);
   templateMCSig->Draw("same");

   //bkg
   Double_t normBkgSigregion = (rPsproper_sig->Parameter(0))*(1 - FsigNew);
   TF1 *psProperBackSigregion= new TF1("averageBkgSignalRegion",GetAverageFunctionsSignal,-8000.,8000.,3);
   psProperBackSigregion->FixParameter(0,normBkgSigregion);
   psProperBackSigregion->FixParameter(1,0.);
   psProperBackSigregion->FixParameter(2,0.);
   psProperBackSigregion->SetNpx(500);
   psProperBackSigregion->SetLineColor(3);
   psProperBackSigregion->SetLineWidth(2);
   psProperBackSigregion->Draw("same");
   psproperTot_sig->SetLineColor(kBlack);
   psproperTot_sig->Draw("same");
   //
   TLegend *leg=new TLegend(0.17,0.72,0.42,0.88);
   leg->SetBorderSize(0); leg->SetFillColor(0); leg->SetTextFont(42);
   leg->SetFillStyle(0); leg->SetMargin(0.25); 
   leg->SetEntrySeparation(0.15);
   leg->AddEntry(psproperTot_sig, "fit, all","l");
   leg->AddEntry(promptSig, "fit, prompt J/#psi","l");
   leg->AddEntry(templateMCSig, "fit, secondary J/#psi","l");
   leg->AddEntry(psProperBackSigregion, "fit, bkg","l");
   leg->Draw("same");
   psSignalCan->SaveAs(Form("inputFiles_%1.1f_%1.1f/LikelihoodFitProjSignalRegionMass%s.root",ptMin,ptMax,resTypeSlikeFit.Data())); 
   psSignalCan->SaveAs(Form("inputFiles_%1.1f_%1.1f/LikelihoodFitProjSignalRegionMass%s.pdf",ptMin,ptMax,resTypeSlikeFit.Data())); 

   // save all parameters -> can be re-used as input 
   SaveParameters(Form("inputFiles_%1.1f_%1.1f/startingParametersLast",ptMin,ptMax),allParam,sizeof(allParam)/sizeof(Int_t));
   return;
}


void AddNtupla(TNtuple *ntToFill, TNtuple *ntNew, Double_t ptmin, Double_t ptmax, Double_t mMin, Double_t mMax,Double_t centmin, Double_t centmax){
  // add mc ntupla
  ntNew->ResetBranchAddresses();
  Float_t mm , xx, tt, pt, cent;
  ntNew->SetBranchAddress("Xdecaytime",&xx);
  ntNew->SetBranchAddress("Mass",&mm);
  ntNew->SetBranchAddress("Type",&tt);
  ntNew->SetBranchAddress("Pt",&pt);
  ntNew->SetBranchAddress("Centrality",&cent);
  Int_t fNcurrent=0;
  Int_t nb = (Int_t)ntNew->GetEvent(fNcurrent);
  //printf("new ntpula orig %d \n",ntNew->GetEntries()); getchar();

  for (Int_t iev=0; iev<(ntNew->GetEntries()); iev++){
   if(pt < ptmin || pt > ptmax) {fNcurrent++; nb = (Int_t)ntNew->GetEvent(fNcurrent); continue;}
   if(mm < mMin || mm > mMax || mm == mMin || mm == mMax) {fNcurrent++; nb = (Int_t)ntNew->GetEvent(fNcurrent); continue;}
   if(centmin>0 && (cent < centmin || cent >= centmax)) {fNcurrent++; nb = (Int_t)ntNew->GetEvent(fNcurrent); continue;}
   ntToFill->Fill(xx,mm,tt,pt);
   fNcurrent++;
   nb = (Int_t)ntNew->GetEvent(fNcurrent);
  }
 cout <<  nb << endl;
 //printf("new ntpula entries %d \n",ntToFill->GetEntries()); getchar();
 return;
}


void SaveParameters(TString filename, Int_t param[], Int_t size){

   Float_t val[size];
   TFile foutFile(Form("%s.root",filename.Data()),"RECREATE");
   TString parSaved=""; for(int ipar=0; ipar<size-1; ipar++) parSaved += Form("%s:",fitParameterNames[param[ipar]].Data());
   parSaved += Form("%s",fitParameterNames[param[size-1]].Data());
   TNtuple *ntPar = new TNtuple("parameters",Form("%s parameters",filename.Data()),parSaved.Data());
   for(int ipar=0; ipar<size; ipar++) val[ipar] = fithandler->GetParameter(param[ipar]);
   ntPar->Fill(val);
   foutFile.cd();
   ntPar->Write();
   return;
}

void SaveParameters(TString filename, Int_t param[], Int_t size, Double_t paramValues[]){

   Float_t val[size];
   TFile foutFile(Form("%s.root",filename.Data()),"RECREATE");
   TString parSaved=""; for(int ipar=0; ipar<size-1; ipar++) parSaved += Form("%s:",fitParameterNames[param[ipar]].Data());
   parSaved += Form("%s",fitParameterNames[param[size-1]].Data());
   TNtuple *ntPar = new TNtuple("parameters",Form("%s parameters",filename.Data()),parSaved.Data());
   for(int ipar=0; ipar<size; ipar++) val[ipar] = paramValues[ipar];//fithandler->GetParameter(param[ipar]);
   ntPar->Fill(val);
   foutFile.cd();
   ntPar->Write();
   return;
}

void SetParameter(TString parName, Double_t parValue){
        // printf("parName 1 %s \n",parName.Data());
	for(int ipar=0; ipar<kNumPar; ipar++) {
	if(fitParameterNames[ipar] == parName) paramInputValues[ipar] = parValue;
	}
	return;
}

void SetStartingParameters(TString filename){
  
  TFile f(filename);
  TNtuple *ntParam = (TNtuple*)f.Get("parameters");
  Int_t nBranch = ntParam->GetNbranches();
  Float_t parVal[nBranch]; 
  TString parName;
  for(int ipar=0; ipar<nBranch; ipar++) {
  parName = ntParam->GetListOfBranches()->At(ipar)->GetName();
  ntParam->SetBranchAddress(parName, &parVal[ipar]); 
  }
  ntParam->GetEvent(0);
  for(int ipar=0; ipar<nBranch; ipar++) {
  parName = ntParam->GetListOfBranches()->At(ipar)->GetName();
  // printf("par name %s",parName.Data());
  SetParameter(parName, parVal[ipar]);
  }

return;
}

void GetParameters(TString filename, Float_t parameters[]){
//
  TFile f(filename);
  TNtuple *ntParam = (TNtuple*)f.Get("parameters");
  Int_t nBranch = ntParam->GetNbranches();
  Float_t parVal[nBranch];
  TString parName;
  for(int ipar=0; ipar<nBranch; ipar++) {
  parName = ntParam->GetListOfBranches()->At(ipar)->GetName();
  ntParam->SetBranchAddress(parName, &parameters[ipar]);
  }
  ntParam->GetEvent(0);
 return;  
}

void SetReleasedParameters(Int_t param[], Int_t npar){

   Bool_t toBeReleased=kFALSE;
   for(int ipar=0; ipar<kNumPar; ipar++) {
          toBeReleased=kFALSE;
          for(int im=0; im<npar; im++) { if(ipar == param[im]) toBeReleased=kTRUE; }
          if(!toBeReleased) fithandler->FixParam(ipar,kTRUE);
        }
   return;
}


TF1 *SetupResolutionFunction(TString name, Int_t resolutionParam[], Int_t type){
  // set up resolution function for chi2 fit

   Int_t sizeResol = sizeof(resolutionParamFF)/sizeof(Int_t) + 2;
   TF1 *resolutionFuncFit = new TF1(name,CDFResolutionFunction,-6000.,6000.,sizeResol);
   Double_t paramInputCDFResol[sizeResol];
   for(int ipar=0; ipar<sizeResol-2; ipar++) paramInputCDFResol[ipar] = paramInputValues[resolutionParam[ipar]];
   resolutionFuncFit->SetParameters(paramInputCDFResol);
   resolutionFuncFit->FixParameter(9,(Double_t)type); // type  
   resolutionFuncFit->SetParLimits(1, 0., 100.); // mean G1
   resolutionFuncFit->SetParLimits(4, 0., 100.); // mean G2
   for(int ipar=0; ipar<sizeResol-2; ipar++) resolutionFuncFit->SetParName(ipar, fitParameterNames[resolutionParam[ipar]]);
   resolutionFuncFit->SetParName(9, "SPDtype");
   resolutionFuncFit->SetParName(10,"Normalization");
   resolutionFuncFit->SetNpx(5000);
   return resolutionFuncFit;
}

Double_t CDFResolutionFunction(Double_t *x, Double_t *par)
{
 // resolution function
 Double_t normalization = par[10];
 Double_t type = par[9];
 likely_obj->SetResolutionConstants(par, (Int_t)type);
 return normalization*likely_obj->ResolutionFunc(x[0],-1.,(Int_t)type);
}



Double_t InvariantMassSignalFunction(Double_t *x, Double_t *par){
 // invariant mass signal function (Crystal Ball)
 likely_obj->SetCrystalBallMmean(par[0]); 
 likely_obj->SetCrystalBallNexp(par[1]);
 likely_obj->SetCrystalBallSigma(par[2]);
 likely_obj->SetCrystalBallAlpha(par[3]);
 likely_obj->SetCrystalBallNorm(par[4]);
 return likely_obj->EvaluateCDFInvMassSigDistr(x[0]);
}

TF1* SetupInvariantMassSignalFunction(TString name, Double_t lowEd, Double_t upEd){
 // set up invariant mass signal function for chi2 fit 
 Int_t sizeSig = sizeof(invMassSignalParam)/sizeof(Int_t);
 TF1 *invMassSignalFunc = new TF1(name,InvariantMassSignalFunction,lowEd,upEd,sizeSig);
 invMassSignalFunc->SetNpx(5000);
 //Double_t paramInputInvMass[sizeSig];
 for(int ipar=0; ipar<sizeSig; ipar++) 
 {
 //paramInputInvMass[ipar] = paramInputValues[invMassSignalParam[ipar]]; 
 invMassSignalFunc->SetParameter(ipar,paramInputValues[invMassSignalParam[ipar]]);
 invMassSignalFunc->SetParName(ipar, fitParameterNames[invMassSignalParam[ipar]]);
 }
 return invMassSignalFunc;
}

Double_t InvariantMassSignalPlusBkgExpoFunction(Double_t *x, Double_t *par){
 //  set up invariant mass signal+background function for chi2 fit 
 // invariant mass signal function (Crystal Ball)
 likely_obj->SetCrystalBallMmean(par[0]);
 likely_obj->SetCrystalBallNexp(par[1]);
 likely_obj->SetCrystalBallSigma(par[2]);
 likely_obj->SetCrystalBallAlpha(par[3]);
 likely_obj->SetCrystalBallNorm(par[4]);
 // Fsig
 Double_t FsigChi2 = par[5];
// invariant mass bkg function  [expo]
 likely_obj->SetBkgInvMassNorm(par[6]);
 likely_obj->SetBkgInvMassMean(par[7]);
 likely_obj->SetBkgInvMassSlope(par[8]);
 likely_obj->SetBkgInvMassConst(par[9]);
 likely_obj->SetBkgInvMassPolyn4(par[10]);
 likely_obj->SetBkgInvMassPolyn5(par[11]);

  
 Double_t kGlobNorm = par[12];
 if(computeIntMass) likely_obj->ComputeMassIntegral(); // Sig(x) and Bkg(x) should be normalized to 1 in bandLow - bandUp
 // InvMass(x) = Fsig*Sig(x) + (1-Fsig)*Bkg(x)  -> normalized to 1 in bandLow - bandUp 
 return (FsigChi2*likely_obj->EvaluateCDFInvMassSigDistr(x[0])/likely_obj->GetIntegralMassSig() + (1-FsigChi2)*likely_obj->EvaluateCDFInvMassBkgDistr(x[0])/likely_obj->GetIntegralMassBkg())*kGlobNorm;
}


TF1 *SetupInvariantMassSignalPlusBkgFunction(TString name, Double_t lowEd, Double_t upEd, Bool_t isExpo){
 //  set up invariant mass signal+background function for chi2 fit 
 // this fit is used to determine inv mass bkg parameters
 // (signal fixed from MC)
 Int_t sizeSig = sizeof(invMassSignalParam)/sizeof(Int_t);
 Int_t sizeBkg = sizeof(invMassBkgExpoParam)/sizeof(Int_t);
 TF1 *invMassSignalPlusBkgFunc = new TF1(name,InvariantMassSignalPlusBkgExpoFunction,lowEd,upEd,sizeSig+sizeBkg+2);
 invMassSignalPlusBkgFunc->SetNpx(5000);
 for(int ipar=0; ipar<sizeSig; ipar++){ 
 invMassSignalPlusBkgFunc->FixParameter(ipar,paramInputValues[invMassSignalParam[ipar]]);
 invMassSignalPlusBkgFunc->SetParName(ipar, fitParameterNames[invMassSignalParam[ipar]]);
 }
 invMassSignalPlusBkgFunc->SetParameter(sizeSig, 0.2); // Set Fsig starting value
 invMassSignalPlusBkgFunc->SetParLimits(sizeSig, 0.,1.); // Set Fsig starting value
 invMassSignalPlusBkgFunc->SetParName(sizeSig,"Fsig");
 for(int ipar=0; ipar<sizeBkg; ipar++){
 invMassSignalPlusBkgFunc->SetParameter(sizeSig+1+ipar,paramInputValues[invMassBkgExpoParam[ipar]]);
 invMassSignalPlusBkgFunc->SetParName(sizeSig+1+ipar, fitParameterNames[invMassBkgExpoParam[ipar]]);
 if(!isExpo) invMassSignalPlusBkgFunc->SetParName(sizeSig+1+ipar, Form("xPolynOrd%d",ipar));
 }
 invMassSignalPlusBkgFunc->SetParName(sizeSig+sizeBkg+1,"Normalization");
 
 return invMassSignalPlusBkgFunc;

}

Double_t GetMeanMass(TNtuple *ntCp, Double_t massMin, Double_t massMax, Double_t ptMin, Double_t ptMax,Int_t type){
  // return average of invariant mass for candidates in the specified mass / pt / type interval
  ntCp->ResetBranchAddresses();
  Float_t mm , xx, tt, pt;
  ntCp->SetBranchAddress("Xdecaytime",&xx);
  ntCp->SetBranchAddress("Mass",&mm);
  ntCp->SetBranchAddress("Type",&tt);
  ntCp->SetBranchAddress("Pt",&pt);
  Int_t fNcurrent=0;
  Int_t nb = (Int_t)ntCp->GetEvent(fNcurrent);
  Double_t meanMass = 0.;
  Double_t candSelected = 0.;


  for (Int_t iev=0; iev<(ntCp->GetEntries()); iev++){

  if((mm>massMin && mm<massMax) && (pt>ptMin && pt<ptMax) && ((Int_t)tt == type) ){
   meanMass += mm;
   candSelected += 1.;
   }

   fNcurrent++;
   nb = (Int_t)ntCp->GetEvent(fNcurrent);
  }

 if(candSelected == 0) {printf("mMin = %f - mMax = %f ------ ptMin = %f - ptMax = %f \n",massMin,massMax,ptMin,ptMax); return 0.;}
 meanMass = meanMass/candSelected;
 //cout <<  nb << endl;
 return meanMass;
}

Double_t*** ComputeWeightsForInterpolation(TNtuple *ntOrig, Bool_t kQuadratic, Bool_t excludeExt){

  // compute weights to be used in the interpolation procedure
  Double_t ***weights = new Double_t**[kMassBins-1];
   for(int im = 0; im <(kMassBins-1); im++) {
   weights[im] = new Double_t*[kPtBins];
   for(int ipt=0; ipt<(kPtBins); ipt++){
    weights[im][ipt] = new Double_t[kTypes];
   }
  }

  // mass mean 
  Double_t means[kMassBins-1][kPtBins][kTypes];
  // mass mean signal region
  Double_t meanMassSignal[kPtBins][kTypes];

  Int_t indMass = 0;
   for(int ipt = 0; ipt<kPtBins;ipt++){
   for(int im=0; im < (kMassBins-1);im++){
      for(int itype=0; itype<kTypes; itype++) {
    if(im == extrRegion) meanMassSignal[ipt][itype] = GetMeanMass(ntOrig,massBins[im],massBins[im+1],ptEdges[ipt],ptEdges[ipt+1],itype); 
    if(im >=extrRegion) indMass = im+1; else indMass = im;
       means[im][ipt][itype] = GetMeanMass(ntOrig, massBins[indMass],massBins[indMass+1],ptEdges[ipt],ptEdges[ipt+1],itype);
       }// type
      printf("\n");
     } // mass
   } // pt

  // compute weights
   Double_t sum = 0.;
   for(int ipt=0; ipt<(kPtBins); ipt++){
      for(int itype=0; itype<kTypes; itype++) {
        sum = 0.;
        for(int im=0; im<(kMassBins-1);im++)
        {
         if((TMath::Abs(means[im][ipt][itype] - meanMassSignal[ipt][itype]))>1.e-08) // generally always true...
         weights[im][ipt][itype] = 1./TMath::Abs(means[im][ipt][itype] - meanMassSignal[ipt][itype]);
         else weights[im][ipt][itype] = 0.;
         if(kQuadratic) weights[im][ipt][itype] = weights[im][ipt][itype]*weights[im][ipt][itype]; 
         if(excludeExt && (im == 0 || im == kMassBins-2)) weights[im][ipt][itype] = 0.;  
         sum += weights[im][ipt][itype];
        }
        printf("\n");

        for(int sm=0; sm<(kMassBins-1);sm++) {
                if(sum != 0) weights[sm][ipt][itype] = weights[sm][ipt][itype]/sum;
                else weights[sm][ipt][itype] = 0.;
                printf("weights[mass:%d][pt:%d][type:%d] = %f ",sm,ipt,itype,weights[sm][ipt][itype]);
                }
        }
         printf("\n");
   }
  return weights;
}

Double_t ***LoadResolutionParameters(Double_t ptLimits[]){
///
const Int_t sizeResol = sizeof(resolutionParamFF)/sizeof(Int_t);
  Double_t ***resolParamAll = new Double_t**[kPtBins]; 
  TString arrType[]={"SS","FS","FF"};
   for(int ipt=0; ipt<(kPtBins); ipt++){
      resolParamAll[ipt] = new Double_t*[kTypes];
      for(int itype=0; itype<kTypes; itype++) resolParamAll[ipt][itype] = new Double_t[sizeResol];
     }

  Float_t resParam[sizeResol];
  // load resolutions parameters
  for(int ipt=0; ipt<kPtBins; ipt++){
        for(int itype=0; itype<kTypes; itype++){
                // GetParameters(Form("inputFiles_%1.1f_%1.1f/XResolParameters%s.root",ptLimits[ipt],ptLimits[ipt+1],arrType[itype].Data()),resParam); 
                GetParameters(Form("inputFiles_%1.1f_%1.1f/XResolParameters%s_%1.1f_%1.1f.root",kPtBins > 1 ? ptLimits[0]: ptLimits[ipt], kPtBins > 1 ? ptLimits[kPtBins] : ptLimits[ipt+1], arrType[itype].Data(), ptLimits[ipt],ptLimits[ipt+1]),resParam); 
                for(int iparam=0; iparam<sizeResol; iparam++) { resolParamAll[ipt][itype][iparam] = (Double_t)resParam[iparam]; /*printf("resol type %d - pt %d - par %f \n", itype,ipt,resolParamAll[ipt][itype][iparam]); */  }
        }
  }
return resolParamAll;
}

void LoadBackgroundFunctionsAndWeights(AliDielectronBtoJPSItoEleCDFfitFCN *likeObj, TString filename, TNtuple *ntCand,Bool_t kQuadratic, Bool_t excludeExt){
   // load background functions 
   TString arrType[]={"SS","FS","FF"};
   likeObj->SetLoadFunction(kTRUE);
   TFile fFunctions(filename);
   // for(int itype = 0; itype<kTypes;itype++){
   for(int itype = 1; itype<kTypes;itype++){
         for(int ipt =0; ipt<(kPtBins);ipt++){
           for(int ims = 0; ims<kMassBins;ims++){
                if(ims != extrRegion) { likeObj->SetBkgFunction(ims,itype,ipt,(TF1*)fFunctions.Get(Form("xbkgFunc%s_pt%1.1f_%1.1f_mass%1.2f_%1.2f",arrType[itype].Data(),ptEdges[ipt],ptEdges[ipt+1],massBins[ims],massBins[ims+1])));
                   printf("%s \n",Form("xbkgFunc%s_pt%1.1f_%1.1f_mass%1.2f_%1.2f",arrType[itype].Data(),ptEdges[ipt],ptEdges[ipt+1],massBins[ims],massBins[ims+1]));
		 }   
	}
         }
   }

     //load SS (dummy), not used 
     for(int ipt =0; ipt<(kPtBins);ipt++){
     int itype=0;
           for(int ims = 0; ims<kMassBins;ims++){
                if(ims != extrRegion)  likeObj->SetBkgFunction(ims,itype,ipt,(TF1*)fFunctions.Get(Form("xbkgFuncFF_pt%1.1f_%1.1f_mass%1.2f_%1.2f",ptEdges[ipt],ptEdges[ipt+1],massBins[ims],massBins[ims+1])));
           }
         }

 
  // compute mass weights 
  Double_t ***weights = ComputeWeightsForInterpolation(ntCand, kQuadratic,excludeExt);
  likeObj->SetBkgWeights(weights);
  likeObj->SetFunctionsSaved(15000,10000,10000.,20000.,extrRegion);
  
  return;
}


Double_t GetAverageFunctions(Double_t *x, Double_t *par){
//return averaged functions (used for drawing only)
Double_t normalization = par[0];
Double_t fSignal = par[1];
Double_t fB = par[2];

Double_t val = 0.;
Double_t valRes = 0.;
Double_t valB = 0.;
Double_t valBkg = 0.;
for(int itype=0; itype<kTypes; itype++)  {
	for(int ipt=0; ipt<kPtBins; ipt++) {
        Double_t massTotalWeight = 0.;
        for(int ims=0; ims <kMassBins; ims++) {massTotalWeight += nCandPtMassWnd[ims][ipt][itype];
         }  
       if(massTotalWeight > 0.){
           if(fSignal > 0. && fB < 1.) valRes += ((TF1*)(likely_obj->GetResolutionFunc(-3000.,3000.,massTotalWeight,(ptEdges[ipt+1]+ptEdges[ipt])/2.,itype)))->Eval(x[0]);
           if(fSignal > 0. && fB > 0.) valB += massTotalWeight*((TF1*)likely_obj->GetFunBFunction(ipt, itype))->Eval(x[0]); 

   if(fSignal != 1.){
    	for(int ims=0; ims<kMassBins;ims++){
    	if(nCandPtMassWnd[ims][ipt][itype] > 0.)  valBkg += nCandPtMassWnd[ims][ipt][itype]*((TF1*)likely_obj->GetBkgFunction(ims,ipt,itype))->Eval(x[0]);
    	}
       }
     } // close massTotalWeight > 0
   } // close pt loop
  } // close type loop 

val =( 1-fSignal)*valBkg + fSignal*((1-fB)*valRes + fB*valB);

return val*normalization;
}

Double_t GetAverageFunctionsSignal(Double_t *x, Double_t *par){
//return averaged functions (used for drawing only)
Double_t normalization = par[0];
Double_t fSignal = par[1];
Double_t fB = par[2];

Double_t val = 0.;
Double_t valRes = 0.;
Double_t valB = 0.;
Double_t valBkg = 0.;
for(int itype=0; itype<kTypes; itype++)  {
for(int ipt=0; ipt<kPtBins; ipt++) {

  Double_t massTotalWeight = nCandPtMassWndSignal[ipt][itype];

if(massTotalWeight > 0.){
    if(fSignal > 0. && fB < 1.) valRes += ((TF1*)(likely_obj->GetResolutionFunc(-3000.,3000.,massTotalWeight,(ptEdges[ipt+1]+ptEdges[ipt])/2.,itype)))->Eval(x[0]);
    if(fSignal > 0. && fB > 0.) valB += massTotalWeight*((TF1*)likely_obj->GetFunBFunction(ipt, itype))->Eval(x[0]); 

   if(fSignal != 1){
     if(nCandPtMassWndSignal[ipt][itype] > 0.)  valBkg += nCandPtMassWndSignal[ipt][itype]*((TF1*)likely_obj->GetBkgFunction(extrRegion,ipt,itype))->Eval(x[0]);
    }
   }
  }
 }

val =( 1-fSignal)*valBkg + fSignal*((1-fB)*valRes + fB*valB);

return val*normalization;
}

TF1 *SetupBackgroundFunction(TString name, Double_t bkgParam[], TString type){
  // set up background function for drawing or saving x-background
   Int_t sizeBkg = sizeof(psProperBkgParamAll)/sizeof(Int_t) + 3;
   TF1 *backgrondFuncFit = new TF1(name,CDFxBackgroundFunction,-20000.,20000.,sizeBkg);
   bkgParam[0] = 10.; // resolution weight  
   if(type=="FF") bkgParam[sizeBkg-2] = 2;
   if(type=="FS") bkgParam[sizeBkg-2] = 1;
   if(type=="SS") bkgParam[sizeBkg-2] = 0;
   bkgParam[sizeBkg-1] = 1000.; // normalization
   backgrondFuncFit->SetParameters(bkgParam);
   for(int iparam=0; iparam<(sizeBkg-3); iparam++) backgrondFuncFit->SetParName(iparam+1,fitParameterNames[psProperBkgParamAll[iparam]]);
   backgrondFuncFit->SetParName(0, "fResWeight");
   backgrondFuncFit->SetParName(sizeBkg-2, "SPDType");
   backgrondFuncFit->SetParName(sizeBkg-1, "normalization");

   for(int ij=0; ij<5; ij++) backgrondFuncFit->SetParLimits(ij,0.,1.e+06);
   backgrondFuncFit->SetParLimits(5,0.,0.01);
   backgrondFuncFit->SetParLimits(6,0.,0.01);
   backgrondFuncFit->SetParLimits(7,0.,0.02); 
   backgrondFuncFit->SetParLimits(8,0.,0.02); 
   backgrondFuncFit->FixParameter(sizeBkg-2, bkgParam[sizeBkg-2]);
   backgrondFuncFit->SetNpx(10000);
   return backgrondFuncFit;
}


Double_t CDFxBackgroundFunction(Double_t *x, Double_t *par){
 //x-background function 
 likely_obj->SetResWeight(par[0]);
 likely_obj->SetFPlus(par[1]);
 likely_obj->SetFMinus(par[2]);
 likely_obj->SetFSym(par[3]);
 likely_obj->SetFSym1(par[4]);  // not used in pp / pPb 
 likely_obj->SetLamPlus(par[5]);
 likely_obj->SetLamMinus(par[6]);
 likely_obj->SetLamSym(par[7]);
 likely_obj->SetLamSym1(par[8]);  // not used in pp / pPb
 Double_t type = par[9];
 Double_t normalization = par[10];
 return normalization*likely_obj->EvaluateCDFDecayTimeBkgDistr(x[0],(Int_t)type);
}


void SaveFunctions(TString resType, Double_t ptMin, Double_t ptMax, Double_t bandLow, Double_t bandUp){
        // Save functions for x-bkg used for the likelihood fit
        TString arrType[]={"SS","FS","FF"};
        Int_t sizeBkg = sizeof(psProperBkgParamAll)/sizeof(Int_t);
        Double_t bkgParamAll[sizeBkg+3];  
        for(int iparam=0; iparam<sizeBkg; iparam++) bkgParamAll[iparam+1] = fithandler->GetParameter(psProperBkgParamAll[iparam]); 
        bkgParamAll[0] = fithandler->GetParameter(kWResolution); // resolution weight 

        AliDielectronBtoJPSItoEleCDFfitFCN* likely_obj_new = new AliDielectronBtoJPSItoEleCDFfitFCN();//fithandler->LikelihoodPointer(); 
        likely_obj_new->SetMultivariateFit(kTRUE);
        likely_obj_new->InitializeFunctions(kPtBins,kMassBins);

        // masswindows
        TArrayD *massWind =new TArrayD(kMassBins+1); // mass windows to define signal and adiacent regions
        for(int im=0; im<kMassBins+1; im++)  massWind->AddAt(massBins[im],im);
        likely_obj_new->SetMassWindows(massWind);

        TArrayD *ptWind = new TArrayD(kPtBins+1); // pt windows
        for(int ipt=0; ipt<kPtBins+1; ipt++) ptWind->AddAt(ptEdges[ipt],ipt);
        likely_obj_new->SetPtWindows(ptWind);
 
        // set interpolation region
        likely_obj_new->SetExtrapolationRegion(extrRegion);
        //// load pt-dep. resolution parameters
        Double_t ***resolParamAll = LoadResolutionParameters(ptEdges);
        likely_obj_new->SetResParams(resolParamAll);
        // set bkg parameters
        Float_t ****bkgParametersAll = 0x0; 
        bkgParametersAll = new Float_t***[kPtBins];
          for(int ii=0; ii<kPtBins; ii++) {
            bkgParametersAll[ii] = new Float_t**[kMassBins-1];
            for(int ij=0; ij<(kMassBins-1); ij++) {bkgParametersAll[ii][ij] = new Float_t*[3];
               for(int ik=0;ik<3;ik++) bkgParametersAll[ii][ij][ik] = new Float_t[9];
            }
          }
      
         for(int ipt =0;ipt<kPtBins;ipt++){ // exclude extrapolation region
           for(int km =0;km< (kMassBins-1);km++){ // exclude extrapolation region
             for(int ktype =0;ktype<3;ktype++){
               for(int ipar=0; ipar < 9; ipar++) {
                   bkgParametersAll[ipt][km][ktype][ipar] = bkgParamAll[ipar]; 
                   //printf("ipar %d %f \n",ipar,bkgParamAll[ipar]);
               }
             }
           }
         }
        likely_obj_new->SetBkgParams(bkgParametersAll);

        TFile fbkgOut(Form("inputFiles_%1.1f_%1.1f/XBkgFunctions.root",kPtBins > 1 ? ptEdges[0]: ptMin, kPtBins > 1 ? ptEdges[kPtBins] : ptMax),"UPDATE"); 
        for(int itype=0; itype<3;  itype++){
           if(resType.Contains(arrType[(Int_t)itype])) {
                likely_obj_new->SetBackgroundSpecificParameters( 0, 0, itype);  
                printf("%f %f %f %f %f %f %f %f %f \n \n",likely_obj_new->GetResWeight(), likely_obj_new->GetFPlus(), likely_obj_new->GetFMinus(), likely_obj_new->GetFSym(), likely_obj_new->GetFSym1(), likely_obj_new->GetLamPlus(), likely_obj_new->GetLamMinus(), likely_obj_new->GetLamSym(), likely_obj_new->GetLamSym1()); 
                TF1 *bkgFunc = (TF1*)likely_obj_new->GetEvaluateCDFDecayTimeBkgDistr(-40000., 40000., 1., itype, (bandLow+bandUp)/2., (ptMax+ptMin)/2.,28000);
                bkgFunc->SetName(Form("xbkgFunc%s_pt%1.1f_%1.1f_mass%1.2f_%1.2f",arrType[itype].Data(),ptMin,ptMax,bandLow,bandUp));
                fbkgOut.cd(); 
                cout << bkgFunc << endl;
                
                 if(fbkgOut.Get(Form("xbkgFunc%s_pt%1.1f_%1.1f_mass%1.2f_%1.2f;1",arrType[itype].Data(),ptMin,ptMax,bandLow,bandUp))) 
                 {
                     fbkgOut.Delete(Form("xbkgFunc%s_pt%1.1f_%1.1f_mass%1.2f_%1.2f;1",arrType[itype].Data(),ptMin,ptMax,bandLow,bandUp));
                 }
                
                bkgFunc->Write();

		if(saveBkgParametersForPtIntegratedCase){
        	TFile fbkgOutPtIntegr(Form("inputFiles_1.5_10.0/XBkgFunctions.root"),"UPDATE"); 
                bkgFunc->Write();
	        fbkgOutPtIntegr.Close();	
		
		}

                }
          }
return;
}
/*
void MakeNtupleFromDstFilteredTree(TString name = "DATA"){
 // create Ntuple used for the fit from dstTreeFileterd.root
 Double_t kConvFromCmtoMicron = 10000.;
 //
 TString filenameOutput= Form("NtuplePP13TeV_%s.root",name.Data());
 TString filenameInput = Form("dstTreeFiltered_%s.root",name.Data());
 TFile fout(filenameOutput,"RECREATE");
 TNtuple* alljpsi  = new TNtuple("fNtupleJPSI","Ntupla JPSI","Xdecaytime:Mass:Type:Pt");
 TFile finput(filenameInput);

 TTree *tree = (TTree*)finput.Get("DstTree");
 tree->Print();

AliReducedEventInfo *ev=new AliReducedEventInfo(); // to be used accordingly  
//AliReducedBaseEvent *ev=new AliReducedBaseEvent();
tree->SetBranchAddress("Event",&ev);
   Int_t fNcurrentNt=0;
   Int_t nbNt = (Int_t)tree->GetEvent(fNcurrentNt);
///
for(int iev=0; iev<(tree->GetEntries()); iev++){
///
if ( iev % ( tree->GetEntries() / 10 ) == 0 ) printf(" At Event %d out of %d \n",iev, (int)tree->GetEntries());
 TClonesArray *clArr=ev->GetPairs();
 for(int ipair=0; ipair<clArr->GetEntriesFast(); ipair++){

 AliReducedPairInfo *pair = (AliReducedPairInfo*)clArr->At(ipair);
 if(!pair) continue;
 if(pair->PairType() == 1) alljpsi->Fill(pair->PsProper()*kConvFromCmtoMicron, pair->Mass(), pair->PairTypeSPD(), pair->Pt());
 }
 fNcurrentNt++;
 nbNt = (Int_t) tree->GetEvent(fNcurrentNt);
}

fout.cd();
alljpsi->Write();

}

*/

