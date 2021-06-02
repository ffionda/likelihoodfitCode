

   AliResonanceFits gFits;

TH1* gOldPointer = 0x0;
Double_t sig;
Double_t sigErr;
Double_t sigOvB;
Double_t sigOvBErr;
Double_t significance;
Double_t chi2;
Double_t yield;

void ExtractJpsiSignal(const Char_t* filename,
                       Double_t ptMin, Double_t ptMax,
                       Double_t centMin, Double_t centMax, Int_t pairTypeMin, Int_t pairTypeMax, Bool_t useLS=kFALSE);
void ConfigureSignalExtractionHybrid(Double_t ptMin, Double_t ptMax,
                                     Double_t centMin, Double_t centMax,
                                    Double_t minFit, Double_t maxFit, Int_t pairTypeMin, Int_t pairTypeMax);
void Draw(AliResonanceFits* fits, Double_t massSignalMin, Double_t massSignalMax,
          Bool_t savePicture=kFALSE, const Char_t* outputDir=".", TString pictureName="standard", TString centralityString = "", TString ptString = "", Bool_t preliminaryPlot=kFALSE);

void ConfigureSignalExtractionFit(Double_t ptMin, Double_t ptMax,
                                  Double_t centMin, Double_t centMax,
                                  Double_t minFit, Double_t maxFit);

void ConfigureSignalExtractionME(Double_t ptMin, Double_t ptMax,
                                 Double_t centMin, Double_t centMax,
                     Double_t minFit, Double_t maxFit,
                                 Double_t minExcl, Double_t maxExcl);

void BeutifyTAxis(TAxis* ax, 
                const Char_t* title, Bool_t centerTitle, Double_t titleSize, Double_t titleOffset, Int_t titleFont,
                Double_t labelSize, Int_t labelFont, Int_t nDivisions);

void BeutifyTH1(TH1* h, const Char_t* title, Double_t lineWidth, Int_t lineColor,
                Int_t markerStyle=1, Int_t markerColor=1);

void ConfigureSignalExtractionLS(Double_t ptMin, Double_t ptMax, 
                                 Double_t centMin, Double_t centMax,
                                 Double_t minFit, Double_t maxFit, 
                                 Double_t minExcl, Double_t maxExcl,Int_t pairTypeMin, Int_t pairTypeMax);

void FillSignalHist(TH1F *ntToFill, TNtuple *ntNew, Double_t ptmin, Double_t ptmax, Int_t pairTypeMin, Int_t pairTypeMax,Double_t mMin, Double_t mMax,Double_t centmin, Double_t centmax);
//void getSignals(const Char_t* filename, Double_t ptMin, Double_t ptMax, Double_t centMin, Double_t centMax);


//I would check different PID configurations for pT > 1 GeV/c, pT > 1.5 GeV/c, pT > 2 GeV/c for centralities 0-10%, 10-30%, 30-50%, 50-90%
// and same for some pT bins, in central semicentral collisions i.e. 1-3 GeV/c, 1.5-3 GeV/c, 2-3 GeV/c, 3-5 GeV/c, 5-10 GeV/c in 0-10% and 30-50%


//void getSignalVariations(const Char_t* filename, TString outputDirectory);
/*void getAllVariations(const Char_t* filename) {
    
    const Int_t kNptBins = 9;
    //Double_t ptBinLims[kNptBins+1] = {
    //   0.0, 1.0, 2.0, 3.0, 4.0,
    //   5.0, 6.0, 8.0, 10.0, 20.0
    // };
    
    Double_t ptBinLims[kNptBins+1] = {
       0.0, 1.0, 2.0, 3.0, 4.0,
       5.0, 6.0, 8.0, 10.0, 20.0
    };
    
    Double_t centBins[10] = {
        
    };
    
    for(Int_t i=1; i<=kNptBins; ++i) {
       Double_t yield = -1.0; Double_t yieldErr = -1.0;
       getSignals(filename, ptBinLims[i-1], ptBinLims[i], 0., 10., "standardCorrected");
       jpsiYields->SetBinContent(i, yield);
       jpsiYields->SetBinError(i, yieldErr);
    }
    
}*/

/*void getSignalVariations(const Char_t* filename, TString outputDirectory) {

    Double_t ptbins[9] = {0.,1.,2.,3.,4.,5.,7.,10.,15.};
    Double_t centralities[9] = {0.0, 5.0, 10.0, 20.,30.,40.,50.,70.,90.};
       for(Int_t icent=0; icent<8; ++icent) {
          for(Int_t ipt=0; ipt<8; ++ipt) {

              ExtractJpsiSignal(filename, ptbins[ipt], ptbins[ipt+1], centralities[icent], centralities[icent+1],"electron30_prot35_pion35");
              printf("pt: %f-%f centralities: %f-%f %% \n",ptbins[ipt],ptbins[ipt+1],centralities[icent],centralities[icent+1]);
              printf("sig +- sigErr: %f +- %f\n",sig,sigErr);
              printf("Signal over background: %f+-%f\n",sigOvB,sigOvBErr);
              printf("Significance: %f\n",significance);
              printf("chi2 over ndf: %f\n",chi2);

              
              
       }  // end loop over proton rej settings
    }  // end loop over electron inclusion setting
}
*/

/*void getSignals(const Char_t* filename, Double_t ptMin = 0.0, Double_t ptMax = 10.0, Double_t centMin = 0.0, Double_t centMax = 10.0) {
    
    Double_t eleCuts[1] = {3.0};
    //Double_t protCuts[3] = {3.0, 3.5, 4.0};
    //Double_t pionCuts[3] = {3.0, 3.5, 4.0};
    Double_t protCuts[1] = {3.5};
    Double_t pionCuts[1] = {3.5};
   // TH1D *SignalExtrHist = new TH1D( "MassSignal", "PIDconf;Signal;SignalExtr",9,0,9);
    //TH1D *SOvBExtrHist = new TH1D( "MassSignal", "PIDconf;SOvB;SOvBExtr",9,0,9);
    //TH1D *SignificanceExtrHist = new TH1D( "MassSignal", "PIDconf;Significance;SignificanceExtr",9,0,9);
    Int_t icount = 0 ;
    for(Int_t iele=0; iele < 1; ++iele) {
       for(Int_t ipro=0; ipro<1; ++ipro) { //Set to 3 when over loop
          for(Int_t ipio=0; ipio<1; ++ipio) { //Set to 3 when over loop

              ExtractJpsiSignal(filename, ptMin, ptMax, centMin, centMax, Form("electron%.0f_prot%.0f_pion%.0f",
                                                               eleCuts[iele]*10, protCuts[ipro]*10, pionCuts[ipio]*10));
              printf("sig +- sigErr: %f +- %f\n",sig,sigErr);
              printf("Signal over background: %f+-%f\n",sigOvB,sigOvBErr);
              printf("Significance: %f\n",significance);
              printf("chi2 over ndf: %f\n",chi2);
              printf("%s\n",Form("electron%.0f_prot%.0f_pion%.0f",
                                 eleCuts[iele]*10, protCuts[ipro]*10, pionCuts[ipio]*10));
              printf("Yield: %f\n",yield);
              SignalExtrHist->GetXaxis()->SetBinLabel(icount+1, Form("electron%.0f_prot%.0f_pion%.0f",
                                                                     eleCuts[iele]*10, protCuts[ipro]*10, pionCuts[ipio]*10) );
              SignalExtrHist->SetBinContent(icount+1, sig);
              SignalExtrHist->SetBinError(icount+1,sigErr);
              SOvBExtrHist->GetXaxis()->SetBinLabel(icount+1, Form("electron%.0f_prot%.0f_pion%.0f",
                                                                     eleCuts[iele]*10, protCuts[ipro]*10, pionCuts[ipio]*10) );
              SOvBExtrHist->SetBinContent(icount+1, sigOvB);
              SOvBExtrHist->SetBinError(icount+1,sigOvBErr);
              
              SignificanceExtrHist->GetXaxis()->SetBinLabel(icount+1, Form("electron%.0f_prot%.0f_pion%.0f",
                                                                     eleCuts[iele]*10, protCuts[ipro]*10, pionCuts[ipio]*10) );
                                                                     
              SignificanceExtrHist->SetBinContent(icount+1,significance); 
              icount++;
              
          }  // end loop over pion rej settings
       }  // end loop over proton rej settings
    }  // end loop over electron inclusion setting */
   
   // When not loop, commented
   /*
   TCanvas *plotSigVsPID = new TCanvas("Signal", "Signal",2000,1000);
    SignalExtrHist->SetStats(kFALSE);
    SignalExtrHist->DrawCopy();
    plotSigVsPID->SaveAs("picstest3/SigVsPID.pdf");
    
    TCanvas *plotSOvBVsPID = new TCanvas("Signal over background", "Signal over background",2000,1000);
    SOvBExtrHist->SetStats(kFALSE);
    SOvBExtrHist->DrawCopy();
    plotSOvBVsPID->cd();
    plotSOvBVsPID->SaveAs("picstest3/SOvB.pdf");
    
    TCanvas *plotSignificanceVsPID = new TCanvas("Significance", "Significance",2000,1000);
    SignificanceExtrHist->SetStats(kFALSE);
    SignificanceExtrHist->DrawCopy();
    plotSignificanceVsPID->cd();
    plotSignificanceVsPID->SaveAs("picstest3/significance.pdf"); 
} */


void GetPtDependentYield(const Char_t* filename) {
   //
   //
   //
   const Int_t kNptBins = 9;
   Double_t ptBinLims[kNptBins+1] = {
      0.0, 1.0, 2.0, 3.0, 4.0, 
      5.0, 6.0, 8.0, 10.0, 20.0
   };
   TH1D* jpsiYields = new TH1D("jpsiYield", "", kNptBins, ptBinLims);
   jpsiYields->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   jpsiYields->GetYaxis()->SetTitle("Observed counts");
   
   for(Int_t i=1; i<=kNptBins; ++i) {
      Double_t yield = -1.0; Double_t yieldErr = -1.0;
      ExtractJpsiSignal(filename, ptBinLims[i-1], ptBinLims[i], 0., 10., -1.,-1./*"standardCorrected"*/);
      jpsiYields->SetBinContent(i, yield);
      jpsiYields->SetBinError(i, yieldErr);
   }
   
   TCanvas* cv=new TCanvas();
   TF1* fpl = new TF1("fpl","[0]*x/TMath::Power((1+TMath::Power(x/[1],[3])),[2])",0.,30.);
   fpl->SetParameters(1.0,3.67987,3.00625,2.0);
   fpl->FixParameter(3, 2.0);
   //jpsiYields->Fit(fpl, "QME", "", 0.0, 10.0);
   jpsiYields->Draw();
}

void ExtractJpsiSignal(const Char_t* filename, 
                       Double_t ptMin=-1.0, Double_t ptMax=-1,
                       Double_t centMin=-1.0, Double_t centMax=-1.0, Int_t pairTypeMin = -1.0, Int_t pairTypeMax = -1.0, Bool_t useLS = kFALSE) {
   
    //SET PAIRTYPE MIN AND PAIRTYPE MAX   
    //FF: 1.5 - 2.5
    //FS 0.5 - 1.5
    //SS -0.5-1.5
   // Double_t pairTypeMin = 0.51;
    //Double_t pairTypeMax = 2.49;
    if(useLS) gSystem->Exec(Form("mkdir -p LSFiles"));
    gSystem->Exec(Form("mkdir -p MixedEventFiles"));
    const Char_t* setting="standard";
    
   AliHistogramManager* histManData = new AliHistogramManager();
    std::string s = setting;
   histManData->InitFile(filename, "jpsi2eeHistos");
   //histManData->InitFile(filename, "jpsi2eeHistos_normalEvents");
   THnF* seos = (THnF*)histManData->GetHistogram(Form("PairSEPM_%s",setting),"PairInvMass");
   THnF* meos = (THnF*)histManData->GetHistogram(Form("PairMEPM_%s",setting),"PairInvMass");
   THnF* sels1 = (THnF*)histManData->GetHistogram(Form("PairSEPP_%s",setting),"PairInvMass");
   THnF* sels2 = (THnF*)histManData->GetHistogram(Form("PairSEMM_%s",setting),"PairInvMass");
   THnF* mels1 = (THnF*)histManData->GetHistogram(Form("PairMEPP_%s",setting),"PairInvMass");
   THnF* mels2 = (THnF*)histManData->GetHistogram(Form("PairMEMM_%s",setting),"PairInvMass");
   
   TH1F* events = (TH1F*)histManData->GetHistogram("Event_AfterCuts", "CentVZERO");
   cout << "Events in " << centMin << " - " << centMax << " interval: ";
   cout << events->Integral(events->GetXaxis()->FindBin(centMin+0.001), events->GetXaxis()->FindBin(centMax-0.001)) << endl;
   
   gFits.SetSEOSHistogram(seos);
   gFits.SetMEOSHistogram(meos);
   gFits.SetSELSHistograms(sels1, sels2);
   gFits.SetMELSHistograms(mels1, mels2);
   
   
   TFile* mcShapeFile=TFile::Open("LHC15o_mcShape.root");
   TH1F* mcShape = (TH1F*)mcShapeFile->Get("Mass");
   mcShape->Reset();
   TFile* mcNtupla=TFile::Open("./InputRootFiles/NtuplePP5TeV_MC_all.root"); 
   TNtuple *ntupla = (TNtuple*)mcNtupla->Get("fNtupleJPSI");
   FillSignalHist(mcShape,ntupla,ptMin,ptMax,pairTypeMin,pairTypeMax,0.,5.,centMin,centMax); 
   
   
   //TFile* mcNtupla=TFile::Open("./InputRootFiles/NtuplePP5TeV_MC_OS_tightPID.root"); 
   //TH1F* mcShape = (TH1F*)mcNtupla->Get("fInvMassAll");
   //TFile* mcShapeFile=TFile::Open("LHC15o_mcShape.root");
   //TH1F* mcShape = (TH1F*)mcShapeFile->Get("Mass");
   gFits.SetSignalMCshape(mcShape);
   
   Double_t fitPtMin = 1.8;
   Double_t fitPtMax = 4.2;
   if(ptMax<1.1 && ptMax>0.9) {
      fitPtMin = 2.1; fitPtMax = 3.9;
   }
   if(ptMax>1.1 && ptMax<2.1) {
      fitPtMin = 2.0; fitPtMax = 4.0;
   }
   //ConfigureSignalExtractionME(ptMin, ptMax, centMin, centMax, 1.5,4.5,2.5,3.2);
   if(useLS) ConfigureSignalExtractionLS(ptMin, ptMax, centMin, centMax, 1.5, 4.5, 2.5, 3.2, pairTypeMin, pairTypeMax);
   //ConfigureSignalExtractionFit(ptMin, ptMax, centMin, centMax, 1.8, 4.2);
   else ConfigureSignalExtractionHybrid(ptMin, ptMax, centMin, centMax, 2.0, 3.7, pairTypeMin, pairTypeMax);
   //ConfigureSignalExtractionHybrid(ptMin, ptMax, centMin, centMax, 2.0, 4.0);
   
   gFits.Process();
   Double_t* fitVals = gFits.ComputeOutputValues(2.92, 3.16);
   //yield = fitVals[AliResonanceFits::kChisqMCTotal];
   //err = fitVals[AliResonanceFits::kSignifErr];
   
   
   TString outputFileName = Form("MixedEventFiles/MixedEventResultsPt%.1f%.1fCent%.1f%.1fCand%.1d%.1d.root",ptMin,ptMax,centMin,centMax,pairTypeMin,pairTypeMax);
   if(useLS) outputFileName = Form("LSFiles/LSResultsPt%.1f%.1fCent%.1f%.1fCand%.1d%.1d.root",ptMin,ptMax,centMin,centMax,pairTypeMin,pairTypeMax);
   TFile foutput(outputFileName,"RECREATE"); 
   TH1 *mixEvHist = (TH1*)gFits.GetBkg();
   TH1 *hsignal = (TH1*)gFits.GetSignal();
   //TH1 *mixEvHist = (TH1*)gFits.GetBkgCombinatorial();
   TF1 *resFit = (TF1*)gFits.GetBkgFitFunction();
   foutput.cd();
   hsignal->SetName("fSig");
   mixEvHist->SetName("fBkg");
   hsignal->Write();
   mixEvHist->Write();
   if(!useLS){
   TF1 *resFit = (TF1*)gFits.GetBkgFitFunction();
   resFit->Write();
   TH1* SplusB = (TH1*)gFits.GetSplusB();
   SplusB->Write();
   }
    sig = fitVals[AliResonanceFits::kSig];
    sigErr = fitVals[AliResonanceFits::kSigErr];
    sigOvB = fitVals[AliResonanceFits::kSoverB];
    sigOvBErr = fitVals[AliResonanceFits::kSoverBerr];
    significance = fitVals[AliResonanceFits::kSignif];
    chi2 = fitVals[AliResonanceFits::kChisqMCTotal];
    yield = sig/events->Integral(events->GetXaxis()->FindBin(centMin+0.001), events->GetXaxis()->FindBin(centMax-0.001));
    TString name = Form("%s",setting);
    printf("name: %s",setting);
   //cout << "sig = " << fitVals[AliResonanceFits::kSig] << endl;
   //cout << "err = " << fitVals[AliResonanceFits::kSigErr]    << endl;
   
   TString centralityString = Form("Cent_%0.0f_%0.0f", centMin, centMax);
   TString ptString = Form("pt_%0.2f_%0.0f", ptMin, ptMax); 
    
   if(useLS) Draw(&gFits, 2.92, 3.16, kTRUE, "./LSFiles",s, ptString, centralityString);
   else Draw(&gFits, 2.92, 3.16, kTRUE, "./MixedEventFiles",s, ptString, centralityString);
}


//____________________________________________________________________________________________________
void ConfigureSignalExtractionME(Double_t ptMin = -1.0, Double_t ptMax = -1.0,
                                 Double_t centMin = -1.0, Double_t centMax = -1.0,
                     Double_t minFit = 1.5, Double_t maxFit = 4.72, 
                     Double_t minExcl = 2.5, Double_t maxExcl = 3.72) {
   // 
   //  configure the AliResonanceFits object
   //
   gFits.AddVariable(AliReducedVarManager::kCentVZERO, 2);
   //gFits.AddVariable(AliReducedVarManager::kSPDntracklets, 3);
   gFits.AddVariable(AliReducedVarManager::kPt, 1);
   gFits.AddVariable(AliReducedVarManager::kMass, 0);
   gFits.SetVarRange(AliReducedVarManager::kMass, 1.0, 5.0);
   if (ptMin >= 0.0 || ptMax >= 0.0)
      gFits.SetVarRange(AliReducedVarManager::kPt, ptMin+1.0e-6, ptMax-1.0e-6);
   if (centMin >= 0.0 || centMax >= 0.0)
      gFits.SetVarRange(AliReducedVarManager::kCentVZERO, centMin+1.0e-6, centMax-1.0e-6);

   gFits.SetBkgMethod(AliResonanceFits::kBkgMixedEvent);
   gFits.SetMEMatchingMethod(AliResonanceFits::kMatchSEOS); // kMatchSEOS, kMatchSELS
   gFits.SetLSmethod(AliResonanceFits::kLSArithmeticMean);  // kLSArithmeticMean, kLSGeometricMean
   gFits.SetUseRfactorCorrection(kTRUE);
   gFits.SetMassFitRange(minFit,maxFit);
   gFits.SetMassExclusionRange(minExcl,maxExcl);
   gFits.SetScaleSummedBkg(kTRUE);
   //gFits.SetUseSignificantZero(kFALSE);
   gFits.SetScalingOption(AliResonanceFits::kScaleEntries); // kScaleEntries, kScaleFit, kScaleWeightedAverage
   gFits.SetWeightedAveragePower(2.0);                      // 2 - is default, but it can be modified in some specific situations
   // kMinuitMethodChi2, kMinuitMethodLikelihood
   gFits.SetMinuitFitOption(AliResonanceFits::kMinuitMethodChi2);
   
}

//____________________________________________________________________________________________________
void ConfigureSignalExtractionLS(Double_t ptMin = -1.0, Double_t ptMax = -1.0, 
                                 Double_t centMin = -1.0, Double_t centMax = -1.0,
                                 Double_t minFit = 1.5, Double_t maxFit = 4.72, 
                                 Double_t minExcl = 2.5, Double_t maxExcl = 3.72,Int_t pairTypeMin = -1, Int_t pairTypeMax = -1) {
      // 
      //  configure the AliResonanceFits object
      //
      gFits.AddVariable(AliReducedVarManager::kCentVZERO, 2);
      //gFits.AddVariable(AliReducedVarManager::kSPDntracklets, 3);
      gFits.AddVariable(AliReducedVarManager::kPairTypeSPD, 3);
      gFits.SetVarRange(AliReducedVarManager::kPairTypeSPD, pairTypeMin-0.3,pairTypeMax+0.3);
      gFits.AddVariable(AliReducedVarManager::kPt, 1);
      gFits.AddVariable(AliReducedVarManager::kMass, 0);
      gFits.SetVarRange(AliReducedVarManager::kMass, 1.0, 5.0);
      if (ptMin >= 0.0 || ptMax >= 0.0)
         gFits.SetVarRange(AliReducedVarManager::kPt, ptMin+1.0e-6, ptMax-1.0e-6);
      if (centMin >= 0.0 || centMax >= 0.0)
         gFits.SetVarRange(AliReducedVarManager::kCentVZERO, centMin+1.0e-6, centMax-1.0e-6);
      
      gFits.SetBkgMethod(AliResonanceFits::kBkgLikeSign);
      gFits.SetLSmethod(AliResonanceFits::kLSArithmeticMean);  // kLSArithmeticMean, kLSGeometricMean
      gFits.SetUseRfactorCorrection(kTRUE);
      gFits.SetMassFitRange(minFit,maxFit);
      gFits.SetScaleSummedBkg(kTRUE);
      gFits.SetMassExclusionRange(minExcl,maxExcl);
      //gFits.SetUseSignificantZero(kFALSE);
      gFits.SetScalingOption(AliResonanceFits::kScaleEntries); // kScaleEntries, kScaleFit, kScaleWeightedAverage
      gFits.SetWeightedAveragePower(2.0);                      // 2 - is default, but it can be modified in some specific situations
      // kMinuitMethodChi2, kMinuitMethodLikelihood
      gFits.SetMinuitFitOption(AliResonanceFits::kMinuitMethodChi2);

}

//____________________________________________________________________________________________________
void ConfigureSignalExtractionHybrid(Double_t ptMin = -1.0, Double_t ptMax = -1.0,
                                     Double_t centMin = -1.0, Double_t centMax = -1.0,
                                    Double_t minFit = 1.5, Double_t maxFit = 4.72, Int_t pairTypeMin = -1, Int_t pairTypeMax = -1) {
      // 
      //  configure the AliResonanceFits object
      //
      printf("%f and %f\n",centMin, centMax);
      gFits.AddVariable(AliReducedVarManager::kCentVZERO, 2);
      //gFits.AddVariable(AliReducedVarManager::kVZERORP+6*1+1, 3);
      //gFits.AddVariable(AliReducedVarManager::kTPCpileupZAC, 3);
      //gFits.AddVariable(AliReducedVarManager::kTPCpileupContributorsAC, 4);
      gFits.AddVariable(AliReducedVarManager::kPairTypeSPD, 3);
      gFits.SetVarRange(AliReducedVarManager::kPairTypeSPD, pairTypeMin-0.3,pairTypeMax+0.3);
      gFits.AddVariable(AliReducedVarManager::kPt, 1);
      gFits.AddVariable(AliReducedVarManager::kMass, 0);
      gFits.SetVarRange(AliReducedVarManager::kMass, 0.0, 5.0);
      //gFits.SetVarRange(AliReducedVarManager::kTPCpileupZAC, -200.0, 0.0);
      //gFits.SetVarRange(AliReducedVarManager::kTPCpileupContributorsAC, 0.0, 1000.0);
      if (ptMin >= 0.0 && ptMax >= 0.0)
         gFits.SetVarRange(AliReducedVarManager::kPt, ptMin+1.0e-6, ptMax-1.0e-6);
      Int_t test = 0;
      if (centMin >= 0.0 || centMax >= 0.0) test = 1;
      TString testString = Form("Testnummer: %d\n\n",test);
      printf("%s",testString.Data());
      if (centMin >= 0.0 || centMax >= 0.0)
         gFits.SetVarRange(AliReducedVarManager::kCentVZERO, centMin+1.0e-6, centMax-1.0e-6);
      //gFits.SetVarRange(AliReducedVarManager::kNTPCclustersFromPileup, -5.0e+5, 1.0e+5);
      
      gFits.SetBkgMethod(AliResonanceFits::kBkgMixedEventAndResidualFit);
      gFits.SetMEMatchingMethod(AliResonanceFits::kMatchSELS); // kMatchSEOS, kMatchSELS
      gFits.SetLSmethod(AliResonanceFits::kLSArithmeticMean);  // kLSArithmeticMean, kLSGeometricMean
      gFits.SetUseRfactorCorrection(kTRUE);
      if (TMath::Abs(ptMin-0.0)<1.0e-5 && ptMax<2.0 && minFit<1.8)
         gFits.SetMassFitRange(1.8,maxFit);
      else
         gFits.SetMassFitRange(minFit,maxFit);
      gFits.SetMassExclusionRange(2.5,3.2);                 //   used for an initial guess of the bkg function parameters
      gFits.SetUseSignificantZero(kTRUE);
      gFits.SetScalingOption(AliResonanceFits::kScaleEntries); // kScaleEntries, kScaleFit, kScaleWeightedAverage
      gFits.SetWeightedAveragePower(2.0);                      // 2 - is default, but it can be modified in some specific situations
      // kMinuitMethodChi2, kMinuitMethodLikelihood
      gFits.SetMinuitFitOption(AliResonanceFits::kMinuitMethodChi2);
      //TF1* fitFunc = new TF1("fitFunc", "[0]*exp([1]*x)", 0., 5.0);
      //TF1* fitFunc = new TF1("fitFunc", "expo", 0., 5.0);
      //fitFunc->SetParameters(0.1, 0.1);
      TF1* fitFunc = (TF1*)gROOT->FindObject("fitFunc");
      if(fitFunc) delete fitFunc;
      fitFunc = 0x0;
      fitFunc = new TF1("fitFunc", "pol1", 0., 5.0);
      /*if(ptMax<3.1 || Double_t(ptMax-5.0)<0.1)
         fitFunc = new TF1("fitFunc", "pol3", 0., 5.0);
      else
         fitFunc = new TF1("fitFunc", "expo", 0., 5.0);*/
      gFits.SetBkgFitFunction(fitFunc);
}

//____________________________________________________________________________________________________
void ConfigureSignalExtractionFit(Double_t ptMin = -1.0, Double_t ptMax = -1.0, 
                                  Double_t centMin = -1.0, Double_t centMax = -1.0,
                                 Double_t minFit = 1.5, Double_t maxFit = 4.72) {
      // 
      //  configure the AliResonanceFits object
      //
      gFits.AddVariable(AliReducedVarManager::kCentVZERO, 2);
      //gFits.AddVariable(AliReducedVarManager::kSPDntracklets, 3);
      gFits.AddVariable(AliReducedVarManager::kPt, 1);
      gFits.AddVariable(AliReducedVarManager::kMass, 0);
      gFits.SetVarRange(AliReducedVarManager::kMass, 1.0, 5.0);
      if (ptMin >= 0.0 || ptMax >= 0.0)
         gFits.SetVarRange(AliReducedVarManager::kPt, ptMin+1.0e-6, ptMax-1.0e-6);
      if (centMin >= 0.0 || centMax >= 0.0)
         gFits.SetVarRange(AliReducedVarManager::kCentVZERO, centMin+1.0e-6, centMax-1.0e-6);
      
      gFits.SetBkgMethod(AliResonanceFits::kBkgFitFunction);
      gFits.SetMassFitRange(minFit,maxFit);
      gFits.SetUseSignificantZero(kTRUE);
      gFits.SetScalingOption(AliResonanceFits::kScaleEntries);
      gFits.SetMinuitFitOption(AliResonanceFits::kMinuitMethodChi2);
      TF1* fitFunc = new TF1("fitFunc", "pol2/pol3(3)", 1.5, 5.0);
      gFits.SetBkgFitFunction(fitFunc);
}

//________________________________________________________________________________________
void Draw(AliResonanceFits* fits, Double_t massSignalMin, Double_t massSignalMax, 
          Bool_t savePicture=kFALSE, const Char_t* outputDir=".", TString pictureName="standard", TString ptString="",TString centralityString="", Bool_t preliminaryPlot=kFALSE) {
   //
   // Make signal extraction plot using the AliResonanceFits output
   //
   TH1* histSplusB = fits->GetSplusB(); if(!histSplusB) return;
   TH1* histSig = fits->GetSignal();  if(!histSig) return;
   TH1* histBkg = fits->GetBkg(); // if(!histBkg) return;
   TH1* histSplusResidualBkg = fits->GetSplusResidualBkg();
   TH1* histSigMC = fits->GetSignalMC();
   TF1* bkgFit = ((fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEventAndResidualFit || 
   fits->GetBkgMethod()==AliResonanceFits::kBkgFitFunction) ? fits->GetBkgFitFunction() : 0x0);
   TF1* globalFit = ((fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEventAndResidualFit || 
   fits->GetBkgMethod()==AliResonanceFits::kBkgFitFunction) ? fits->GetGlobalFitFunction() : 0x0);
   
   Bool_t histogramsAreUnchanged = kFALSE;
   if(histSplusB==gOldPointer) histogramsAreUnchanged = kTRUE;
   else gOldPointer = histSplusB;
   
   TLatex* lat=new TLatex();
   lat->SetNDC();
   lat->SetTextSize(0.06);
   lat->SetTextColor(1);
   
   TCanvas* oldCanvas = (TCanvas*)gROOT->FindObject("AliResonanceFitsCanvas");
   if(oldCanvas) delete oldCanvas; 
   TCanvas* c1=new TCanvas("AliResonanceFitsCanvas", "AliResonanceFits canvas", 980, 1200);
   c1->SetTopMargin(0.01);
   c1->SetRightMargin(0.005);
   c1->SetBottomMargin(0.01);
   c1->SetLeftMargin(0.005);
   c1->Divide(1,2,0.0,0.0);
   
   // Draw the top pad ================================================
   TVirtualPad* pad = c1->cd(1);
   pad->SetTopMargin(0.01);
   pad->SetRightMargin(0.02);
   pad->SetBottomMargin(0.0);
   pad->SetLeftMargin(0.15);
   pad->SetTicks(1,1);
   
   histSplusB->SetStats(kFALSE);
   if(!histogramsAreUnchanged)
      histSplusB->GetYaxis()->SetRangeUser(0.01, 1.3*histSplusB->GetMaximum());
   //histSplusB->GetXaxis()->SetRangeUser(fits->GetMassFitRange()[0]-0.1, fits->GetMassFitRange()[1]+0.1);
   BeutifyTH1(histSplusB, "", 3.0, 1);
   BeutifyTAxis(histSplusB->GetXaxis(), "", kFALSE, 0.04, 1., 42, 0.04, 42, 510);
   BeutifyTAxis(histSplusB->GetYaxis(), Form("entries per %.0f MeV/#it{c}^{2}", 1000.*histSplusB->GetXaxis()->GetBinWidth(1)), kTRUE, 0.08, 0.84, 42, 0.065, 42, 507);
   histSplusB->GetXaxis()->SetRangeUser(2.0,3.7);
   histSplusB->Draw("EXY"); 
   if(histBkg) {
      histBkg->SetStats(kFALSE);
      BeutifyTH1(histBkg, "", 2.0, 2);
      //histBkg->GetXaxis()->SetRangeUser(fits->GetMassFitRange()[0]-0.1, fits->GetMassFitRange()[1]+0.1);
      histBkg->Draw("sameE");
   }
   
   if(bkgFit && globalFit && fits->GetBkgMethod()==AliResonanceFits::kBkgFitFunction) {
      bkgFit->SetLineColor(4); bkgFit->SetLineWidth(1);
      globalFit->SetLineColor(6); globalFit->SetLineWidth(1);
      globalFit->Draw("same");
      bkgFit->Draw("same");
   }
   
   /*
    *  cout << "Draw 3.2" << endl;
    *  TBox* boxL = new TBox(fits->GetMassFitRange()[0], histSplusB->GetMaximum()*0.1, fits->GetMassExclusionRange()[0], histSplusB->GetMaximum()*0.15);
    *  boxL->SetFillColorAlpha(4, 0.35);
    *  boxL->SetFillStyle(3003);
    *  boxL->Draw();
    *  
    *  cout << "Draw 3.3" << endl;
    *  TBox* boxR = new TBox(fits->GetMassExclusionRange()[1], histSplusB->GetMaximum()*0.1, fits->GetMassFitRange()[1], histSplusB->GetMaximum()*0.15);
    *  boxR->SetFillColorAlpha(4, 0.35);
    *  boxR->SetFillStyle(3003);
    *  boxR->Draw();
    *  cout << "Draw 4" << endl;
    */
   
   TLegend* legend1 = new TLegend(0.18, 0.55, 0.42, 0.73);
   legend1->SetTextSize(0.06);
   legend1->SetTextFont(42);
   legend1->SetFillColor(0);
   legend1->SetBorderSize(0);
   legend1->AddEntry(histSplusB, "Same event", "l");
   if(fits->GetBkgMethod()==AliResonanceFits::kBkgFitFunction) {
      legend1->AddEntry(globalFit, "Global fit", "l");
      legend1->AddEntry(bkgFit, "Bkg fit", "l");
   }
   else if(fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEventAndResidualFit) {
      legend1->AddEntry(histBkg, "ME like sign", "l");
   }
   else 
      legend1->AddEntry(histBkg, (fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEvent ? "Mixed event" : "Like sign"), "l");
   legend1->Draw();
   
   TLine line;
   line.SetLineColor(1);
   line.SetLineStyle(1);
   line.SetLineWidth(1.0);
   line.DrawLine(massSignalMin, 0.0, massSignalMin, histSplusB->GetMaximum()*0.75);
   line.DrawLine(massSignalMax, 0.0, massSignalMax, histSplusB->GetMaximum()*0.75);
   lat->SetTextFont(42);
   lat->SetTextColor(1);
   lat->SetTextSize(0.060);
   
   if(preliminaryPlot) {
      lat->DrawLatex(0.2, 0.90, "ALICE Preliminary");
      lat->DrawLatex(0.2, 0.82, "0#minus90% Xe#minusXe #sqrt{#it{s}_{NN}} = 5.44 TeV");
   }
   
   if(!preliminaryPlot && (fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEvent || fits->GetBkgMethod()==AliResonanceFits::kBkgLikeSign))
      lat->DrawLatex(0.56, 0.82, Form("#chi^{2}/NDF (%.1f, %.1f) = %.2f", fits->GetMassFitRange()[0], fits->GetMassFitRange()[1], fits->GetFitValues()[AliResonanceFits::kChisqSideBands]));
   //if(fEventVsCent) lat->DrawLatex((noYlabels ? 0.52 : 0.56), 0.72, Form("# events = %.2e", fEventVsCent->Integral(fEventVsCent->GetXaxis()->FindBin(fCentralityRange[0]+0.001), fEventVsCent->GetXaxis()->FindBin(fCentralityRange[1]-0.001))));
   /*lat->DrawLatex(0.63, 0.66, "Fit ranges:");
    *  lat->DrawLatex(0.63, 0.60, Form("%.2f<m_{ee}<%.2f GeV/c^{2}", fFitRange[0], fFitRange[1]));
    *  lat->DrawLatex(0.63, 0.54, Form("%.2f<p_{T}<%.2f GeV/c", fPtFitRange[0], fPtFitRange[1]));
    *  lat->DrawLatex(0.63, 0.45, "Signal:");
    *  lat->DrawLatex(0.63, 0.39, Form("%.2f<m_{ee}<%.2f GeV/c^{2}", fSignalRange[0], fSignalRange[1]));
    *  lat->DrawLatex(0.63, 0.33, Form("%.2f<p_{T}<%.2f GeV/c", fPtRange[0], fPtRange[1])); */
   
   // Draw the bottom pad ================================================
   pad = c1->cd(2);
   pad->SetTopMargin(0.0);
   pad->SetRightMargin(0.02);
   pad->SetBottomMargin(0.17);
   pad->SetLeftMargin(0.15);
   pad->SetTicks(1,1);
   
   if(fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEvent ||
      fits->GetBkgMethod()==AliResonanceFits::kBkgLikeSign ||
      fits->GetBkgMethod()==AliResonanceFits::kBkgFitFunction || 
      fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEventAndResidualFit) {
      histSig->SetStats(kFALSE);
   BeutifyTH1(histSig, "", 2.0, 2, 20, 2);
   if(!histogramsAreUnchanged)
      histSig->GetYaxis()->SetRangeUser(histSig->GetMinimum()*1.4, histSig->GetMaximum()*1.8);
   //histSig->GetXaxis()->SetRangeUser(fits->GetMassFitRange()[0]-0.1, fits->GetMassFitRange()[1]+0.1);
   BeutifyTAxis(histSig->GetYaxis(), "", kFALSE, 0.08, 1.2, 42, 0.065, 42, 508);
   BeutifyTAxis(histSig->GetXaxis(), "#it{m}_{e^{#plus}e^{#minus}} (GeV/#it{c}^{2})", kTRUE, 0.08, 0.84, 42, 0.065, 42, 508);
   histSig->GetXaxis()->SetRangeUser(2.0,3.7);
   histSig->Draw("XY");
   if(histSigMC) {
      //histSigMC->GetXaxis()->SetRangeUser(fits->GetMassFitRange()[0]-0.1, fits->GetMassFitRange()[1]+0.1);
      BeutifyTH1(histSigMC, "", 2.0, 1);
      TH1* histSigMCclone = (TH1*)histSigMC->Clone("histSigMCclone");
      if(fits->GetBkgMethod()==AliResonanceFits::kBkgFitFunction) 
         histSigMCclone->Scale(globalFit->GetParameter(0));
      histSigMCclone->Draw("sameHISTC"); //smooth MC shape drawing
   }
      }
      if(fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEventAndResidualFit) {
         histSplusResidualBkg->SetStats(kFALSE);
         BeutifyTH1(histSplusResidualBkg, "", 2.0, 2, 20, 2);
         if(!histogramsAreUnchanged)
            histSplusResidualBkg->GetYaxis()->SetRangeUser(histSplusResidualBkg->GetMinimum()*1.4, histSplusResidualBkg->GetMaximum()*1.8);
         BeutifyTAxis(histSplusResidualBkg->GetYaxis(), "", kFALSE, 0.08, 1.2, 42, 0.065, 42, 508);
         BeutifyTAxis(histSplusResidualBkg->GetXaxis(), "#it{m}_{e^{#plus}e^{#minus}} (GeV/#it{c}^{2})", kTRUE, 0.08, 0.84, 42, 0.065, 42, 508);
         histSplusResidualBkg->GetXaxis()->SetRangeUser(2.,3.7);
         histSplusResidualBkg->Draw("XY");
         bkgFit->SetLineColor(4); bkgFit->SetLineWidth(2);
         globalFit->SetLineColor(6); globalFit->SetLineWidth(2);
         
         globalFit->Draw("same");
         bkgFit->Draw("same");
      }
      line.DrawLine(fits->GetMassFitRange()[0],0.0,fits->GetMassFitRange()[1],0.0);
      
      const Double_t* fitValues = fits->GetFitValues();
      
      line.DrawLine(massSignalMin, histSig->GetMinimum()*0.9, massSignalMin, histSig->GetMaximum()*0.8);
      line.DrawLine(massSignalMax, histSig->GetMinimum()*0.9, massSignalMax, histSig->GetMaximum()*0.8);
      lat->DrawLatex(0.18, 0.92, Form("Signal:    %.0f #pm %.0f", fitValues[AliResonanceFits::kSig], fitValues[AliResonanceFits::kSigErr]));
      lat->DrawLatex(0.18, 0.84, Form("#it{S}/#it{B}:        %.3f #pm %.3f", fitValues[AliResonanceFits::kSoverB], fitValues[AliResonanceFits::kSoverBerr]));
      if(fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEvent || fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEventAndResidualFit) 
         lat->DrawLatex(0.18, 0.76, Form("#it{S}/#sqrt{#it{S}+#it{B}}: %.1f", fitValues[AliResonanceFits::kSignif]));
      if(fits->GetBkgMethod()==AliResonanceFits::kBkgLikeSign) 
         lat->DrawLatex(0.18, 0.76, Form("#it{S}/#sqrt{#it{S}+2#it{B}}: %.1f", fitValues[AliResonanceFits::kSignif]));
      if(fits->GetBkgMethod()==AliResonanceFits::kBkgFitFunction) 
         lat->DrawLatex(0.18, 0.76, Form("#it{S}/#delta_{#it{S}}: %.1f", fitValues[AliResonanceFits::kSignif]));
      if(histSigMC && !preliminaryPlot) {
         lat->DrawLatex(0.18, 0.68, Form("#chi^{2}/NDF (%.1f, %.1f) = %.2f", fits->GetMassFitRange()[0], fits->GetMassFitRange()[1], fitValues[AliResonanceFits::kChisqMCTotal]));
         lat->DrawLatex(0.18, 0.60, Form("Probability = %.2f", fitValues[AliResonanceFits::kFitProbability]));
        //lat->DrawLatex(0.18, 0.36, Form("%s", pictureName.Data()));
         lat->DrawLatex(0.18, 0.44, Form("%s GeV/c",ptString.Data()));
         lat->DrawLatex(0.18, 0.52, Form("%s %%", centralityString.Data()));
         lat->DrawLatex(0.74,0.60, Form("Yield: %f",yield));
         //lat->DrawLatex(0.56, 0.56, Form("MC yield fraction  = %.2f", fitValues[AliResonanceFits::kMCYieldFraction])); 
      }
      TLegend* legend2 = new TLegend(0.70, 0.83, 0.96, 0.97);
      legend2->SetFillColor(0);
      legend2->SetBorderSize(0);
      legend2->SetTextFont(42);
      legend2->SetTextSize(0.06);
      legend2->SetMargin(0.15);
      if(fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEventAndResidualFit)
         legend2->AddEntry(histSplusResidualBkg, "Data", "lp");
      else
         legend2->AddEntry(histSig, "Data", "lp");
      if(fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEventAndResidualFit) {
         legend2->AddEntry(globalFit, "global function", "l");
         legend2->AddEntry(bkgFit, "bkg function", "l");
      }
      else
         if(histSigMC) legend2->AddEntry(histSigMC, "Monte Carlo", "l");
         legend2->Draw();  
      
      if(savePicture) {
         TString pictureDetailedName = pictureName;
          printf("pictureName %s\n", pictureName.Data());
         if(fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEvent) pictureDetailedName += "_BkgME";
         if(fits->GetBkgMethod()==AliResonanceFits::kBkgLikeSign) pictureDetailedName += "_BkgLS";
         if(fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEventAndResidualFit) pictureDetailedName += "_BkgMEresidual";
         if(fits->GetBkgMethod()==AliResonanceFits::kBkgFitFunction) pictureDetailedName += "_BkgFit";
         if(fits->GetScalingOption() == AliResonanceFits::kScaleEntries) pictureDetailedName += "_ScaleEnt";
         if(fits->GetScalingOption() == AliResonanceFits::kScaleWeightedAverage) pictureDetailedName += "_ScaleWav";
         if(fits->GetScalingOption() == AliResonanceFits::kScaleFit) pictureDetailedName += "_ScaleFit";
         if(fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEvent) {
            if(fits->GetMEMatchingMethod()==AliResonanceFits::kMatchSEOS) pictureDetailedName += "_MatchSEOS";
            if(fits->GetMEMatchingMethod()==AliResonanceFits::kMatchSELS) pictureDetailedName += "_MatchSELS";
         }
         if(fits->GetScalingOption() == AliResonanceFits::kScaleFit) {
            if(fits->GetMinuitFitOption() == AliResonanceFits::kMinuitMethodChi2) pictureDetailedName += "_FitChi2";
            if(fits->GetMinuitFitOption() == AliResonanceFits::kMinuitMethodLikelihood) pictureDetailedName += "_FitLikelihood";
         }
         
         c1->SaveAs(Form("%s/%s%s%s.png", outputDir, pictureDetailedName.Data(),ptString.Data(),centralityString.Data()));
         //c1->SaveAs(Form("%s/%s.eps", outputDir, pictureDetailedName.Data()));
      }
}

//________________________________________________________________________________________
void BeutifyTH1(TH1* h, const Char_t* title, Double_t lineWidth, Int_t lineColor,
                Int_t markerStyle=1, Int_t markerColor=1) {
   //
   // set drawing options for a TH1
   //
   if(!h) return;
   h->SetTitle(title);
   h->SetLineWidth(lineWidth); h->SetLineColor(lineColor);
   h->SetMarkerStyle(markerStyle);
   h->SetMarkerColor(markerColor);
}
                
//________________________________________________________________________________________
void BeutifyTAxis(TAxis* ax, 
                const Char_t* title, Bool_t centerTitle, Double_t titleSize, Double_t titleOffset, Int_t titleFont,
                Double_t labelSize, Int_t labelFont, Int_t nDivisions) {
//
// set drawing options for TAxis
//
   if(!ax) return;
   ax->SetTitle(title);
   if(centerTitle) ax->CenterTitle();
   ax->SetTitleSize(titleSize);
   ax->SetTitleOffset(titleOffset);
   ax->SetTitleFont(titleFont);
   ax->SetLabelSize(labelSize);
   ax->SetLabelFont(labelFont);
   ax->SetNdivisions(nDivisions);
}


void FillSignalHist(TH1F *ntToFill, TNtuple *ntNew, Double_t ptmin, Double_t ptmax, Int_t pairTypeMin, Int_t pairTypeMax,Double_t mMin, Double_t mMax,Double_t centmin, Double_t centmax){
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

  for (Int_t iev=0; iev<(ntNew->GetEntries()); iev++){
   if(pt < ptmin || pt > ptmax) {fNcurrent++; nb = (Int_t)ntNew->GetEvent(fNcurrent); continue;}
   if(mm < mMin || mm > mMax || mm == mMin || mm == mMax) {fNcurrent++; nb = (Int_t)ntNew->GetEvent(fNcurrent); continue;}
   if(centmin>0 && (cent < centmin || cent >= centmax)) {fNcurrent++; nb = (Int_t)ntNew->GetEvent(fNcurrent); continue;}
   if(pairTypeMin>0 && (tt < pairTypeMin || tt > pairTypeMax)) {fNcurrent++; nb = (Int_t)ntNew->GetEvent(fNcurrent); continue;}
   //ntToFill->Fill(xx,mm,tt,pt);
   ntToFill->Fill(mm);
   fNcurrent++;
   nb = (Int_t)ntNew->GetEvent(fNcurrent);
  }
 cout <<  nb << endl;
 //printf("new ntpula entries %d \n",ntToFill->GetEntries()); getchar();
 TFile foutSig("InvMass_FF.root","RECREATE");
 ntToFill->Write();
 return;
}


