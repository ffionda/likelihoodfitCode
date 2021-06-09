void composeMEBackgroundforLikelihoodFit(Double_t ptmin = 1.5, Double_t ptmax = 10., Double_t mMin = 2.5, Double_t mMax = 3.6, Double_t centMin = 30., Double_t centMax = 50., Int_t pairTypeMin = 1 , Int_t pairTypeMax = 2, Int_t polynOrd = 4){

TFile *fmix = new TFile(Form("./MixedEventFiles/MixedEventResultsPt%.1f%.1fCent%.1f%.1fCand%.1d%.1d.root",ptmin, ptmax,centMin,centMax,pairTypeMin,pairTypeMax),"READ");

//TH1 *mixScaled = (TH1*)fmix->Get("fBkg_0.744305");
fmix->ls();
TH1 *mixScaled = (TH1*)fmix->Get("fBkg");
TH1 *resBkg = (TH1*)mixScaled->Clone("residualBkg");
TH1 *mixScaledOrig = (TH1*)mixScaled->Clone("mixEvOrig");
//TH1 *sig = (TH1*)fmix->Get("fSig_0.560628");
TH1 *hsignal = (TH1*)fmix->Get("fSig");
resBkg->Reset();
TF1 *func = (TF1*)fmix->Get("BkgFitFunction"); 
TCanvas *cMix = new TCanvas();
mixScaled->SetMarkerStyle(20);
mixScaled->SetMarkerColor(1);
mixScaled->SetLineColor(1);
mixScaledOrig->SetMarkerStyle(24);
mixScaledOrig->SetMarkerColor(1);
mixScaledOrig->SetLineColor(1);

// add residual bkg to mixed event
Double_t edg1, edg2, integral;
for(int ij=1; ij<mixScaled->GetNbinsX()+1; ij++){
edg1 = mixScaled->GetBinLowEdge(ij);
edg2 = mixScaled->GetBinLowEdge(ij+1);
integral = func->Integral(edg1,edg2);
//printf("integral %f \n",integral/(edg2-edg1));
for(int icount=0; icount < (integral/(edg2-edg1)); icount++) resBkg->Fill((edg2+edg1)/2.);
}
mixScaled->Add(resBkg);

//// fit total distribution
mixScaled->Fit(Form("pol%d",polynOrd),"","",mMin,mMax);
TF1 *pol3total = (TF1*)(mixScaled->GetListOfFunctions()->At(0));
mixScaled->GetXaxis()->SetRangeUser(mMin-0.1, mMax+0.1);
mixScaled->DrawCopy();
resBkg->DrawCopy("same");
mixScaledOrig->DrawCopy("same"); // original mixed event

//// save pol3 parameters -> to be used in the likelihood fit
// start from the old ntupla containig bkg parameters
TFile *finvmassparamLike = new TFile(Form("inputFiles_%1.1f_%1.1f/InvariantMassBkgExpoChi2Fit.root",ptmin,ptmax),"READ"); 
TNtuple *invMbkgParam = (TNtuple*)finvmassparamLike->Get("parameters");
TString branches = ""; 
for(int ib=0; ib<invMbkgParam->GetNbranches()-1; ib++) branches.Append(Form("%s:",invMbkgParam->GetListOfBranches()->At(ib)->GetName()));
branches.Append(Form("%s",invMbkgParam->GetListOfBranches()->At(invMbkgParam->GetNbranches()-1)->GetName()));
// printf("branches -> %s \n",branches.Data());

TString pairString = "";
TString arrTypes[]={"SS","FS","FF"}; 
for(int i=pairTypeMin; i<=pairTypeMax; i++) pairString.Prepend(Form("%s_",arrTypes[i].Data())); 
//printf("%s",pairString);

// create a new ntupla for storing ME parameters
TFile *finvmassparamLikeME = new TFile(Form("inputFiles_%1.1f_%1.1f/InvariantMassBkgME%s.root",ptmin,ptmax,pairString.Data()),"RECREATE");
TNtuple *invMbkgME=new TNtuple("parameters","parameters",branches);

// fill ntupla with parameters from ME + residual bkg fit
Float_t val[6]={0.,0.,0.,0.,0.,0.}; 
for(int ij=0; ij<polynOrd+1;ij++) {val[ij] = (Float_t)(pol3total->GetParameter(ij)/pol3total->Integral(mMin,mMax)); printf("param %f \n",val[ij]);}
invMbkgME->Fill(val);
// save
finvmassparamLikeME->cd();
invMbkgME->Write();

Int_t b1s  = hsignal->FindBin(mMin);
Int_t b2s  = hsignal->FindBin(mMax);
Int_t b1bkg  = mixScaled->FindBin(mMin);
Int_t b2bkg  = mixScaled->FindBin(mMax);
Double_t sigSum =0.;
Double_t bkgSum = 0.;
for(int i=b1s; i <=b2s; i++) 
{
   sigSum += hsignal->GetBinContent(i)*hsignal->GetBinWidth(i);
};

for(int i=b1bkg; i <=b2bkg; i++) 
{
   bkgSum += mixScaled->GetBinContent(i)*mixScaled->GetBinWidth(i);
};

Double_t sigFrac = sigSum/(sigSum+bkgSum);

std::ofstream signalFracFile;
  //for(int im=0; im<kMassBins+1; im++) stringmass.Append(Form("_m%d%1.2f",im+1,massBins[im]));
  signalFracFile.open (Form("inputFiles_%1.1f_%1.1f/FractionOfSignalMass%s.txt",ptmin,ptmax,pairString.Data()), std::ofstream::out);
  signalFracFile << Form("%f \n",sigFrac);
  signalFracFile.close();

    cMix->SaveAs(Form("MixedEventFiles/MEbackground_pt_%1.1f_%1.1f_cand_%s.png",ptmin,ptmax,pairString.Data()));
return;

}


