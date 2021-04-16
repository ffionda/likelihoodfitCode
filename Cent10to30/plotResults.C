#include "TMath.h"
Int_t nptbins = 3;
Double_t ptbins[] = {1.5,3.0,5.0,10.0};
TString types[] = {"FF","FF","FF","FF_FS"};
Int_t nmass = 6;
Int_t nmassHyp = 9;
Double_t mbins[9][6] = {
        {2.60, 2.70, 2.80, 3.20, 3.40, 3.60},
        {2.60, 2.72, 2.80, 3.20, 3.38, 3.60},
        {2.60, 2.68, 2.80, 3.20, 3.42, 3.60},
        {2.60, 2.70, 2.78, 3.22, 3.40, 3.60},
        {2.60, 2.70, 2.82, 3.18, 3.40, 3.60},
        {2.60, 2.72, 2.82, 3.22, 3.42, 3.60},
        {2.60, 2.68, 2.78, 3.18, 3.38, 3.60},
        {2.60, 2.72, 2.82, 3.18, 3.38, 3.60},
        {2.60, 2.68, 2.78, 3.22, 3.42, 3.60}
        };
Int_t mbinsSize =  nmassHyp;// sizeof(mbins)/sizeof(*mbins);
//nmassHyp = mbinsSize;        

TGraphErrors *BuildGraphFb(Double_t massBins[], TString type, Bool_t isQuadr, Bool_t fSigFixed, Bool_t mixedEvents, Int_t color, Int_t style);
TGraphErrors *BuildGraphFbCompTypes(Double_t massBins[], TString types[], Bool_t isQuadr, Bool_t fSigFixed, Int_t bkgType,TString nameGraph, Int_t color, Int_t style);
TGraphErrors *BuildGraphFSig(Double_t massBins[], TString type, Bool_t isQuadr, Bool_t fSigFixed, Int_t color, Int_t style);
TGraphErrors *BuildGraphChi2M(Double_t massBins[], TString type, Bool_t isQuadr, Bool_t fSigFixed, Int_t color, Int_t style); 
TGraphErrors *BuildGraphChi2X(Double_t massBins[], TString type, Bool_t isQuadr, Bool_t fSigFixed, Int_t color, Int_t style); 
TGraphErrors *BuildGraphFSig1D(Double_t massBins[], TString type, Int_t color, Int_t style);
TGraphErrors *BuildGraphChi2XProjCompTypes(Double_t massBins[], TString types[], Bool_t isQuadr, Bool_t fSigFixed, Bool_t mixedEvents,TString nameGraph, Int_t color, Int_t style);
TGraphErrors *BuildGraphChi2MProjCompTypes(Double_t massBins[], TString types[], Bool_t isQuadr, Bool_t fSigFixed, Bool_t mixedEvents,TString nameGraph, Int_t color, Int_t style);
void readFbResults(Double_t massBins[], Bool_t pdf);
void plotVsMassEdges(TString nameGraph="",Bool_t isQuadr=kFALSE,Double_t ymin=-1, Double_t ymax=-1);
void plotResults(Int_t num, Bool_t createPdf);
void plotVsAllMassEdges();
void plotVsMassHypothesis(Int_t num);
void plotVsAllMassHypothesis();
void plotAllResults(Bool_t createPdf);
TGraphErrors* computeRMSandAveragesForMassEdges(Int_t massHyp, Bool_t draw);
void computeRMSandAverages();
void computeSystInvMassBkg();


void computeSystInvMassBkg(){

        Double_t averageUnc[nptbins+1]; for(int ipt=0; ipt<nptbins+1; ipt++) averageUnc[ipt] = 0.;
        TGraphErrors *grSyst = new TGraphErrors(nptbins+1);
        TCanvas *avgSystMass = new TCanvas("allM","allM");
        Int_t colors[] = {1,2,4,6,7,8,9,3,15};
        Int_t styles[] = {20,24,21,32,28,31,33,35,22};
        Double_t xx, yy;
        TLegend *legend=new TLegend(0.37,0.30,0.62,0.70);
        legend->SetBorderSize(0); legend->SetFillColor(0); legend->SetTextFont(42);
        legend->SetFillStyle(0); legend->SetMargin(0.25);
        legend->SetEntrySeparation(0.15);
        legend->SetTextSize(0.045);

        for(int ij=0; ij<nmassHyp; ij++) {
        TString path = Form("mass_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f",mbins[ij][0],mbins[ij][1],mbins[ij][2],mbins[ij][3],mbins[ij][4],mbins[ij][5]);
        TString pathLegend = path;
        pathLegend.ReplaceAll("_"," ");
                TGraphErrors *UncertPerMassEdge = computeRMSandAveragesForMassEdges(ij,kFALSE);

                for(int ibin=0; ibin<nptbins+1; ibin++) {

                        UncertPerMassEdge->GetPoint(ibin,xx,yy);
                        UncertPerMassEdge->SetPoint(ibin,xx,UncertPerMassEdge->GetErrorY(ibin));
                        averageUnc[ibin] += UncertPerMassEdge->GetErrorY(ibin)/nmassHyp;
                        // printf("bin %d err %f \n",ibin, UncertPerMassEdge->GetErrorY(ibin));
                        UncertPerMassEdge->SetPointError(ibin,UncertPerMassEdge->GetErrorX(ibin),0.);
                	if(ij == nmassHyp-1) {
                         grSyst->SetPoint(ibin, xx, 1.);
                         grSyst->SetPointError(ibin, UncertPerMassEdge->GetErrorX(ibin), averageUnc[ibin]);
                        }
                }
                UncertPerMassEdge->SetLineColor(colors[ij]);
                UncertPerMassEdge->SetMarkerColor(colors[ij]);
                UncertPerMassEdge->SetMarkerStyle(styles[ij]);
                UncertPerMassEdge->GetXaxis()->SetTitle("p_{T}");
                UncertPerMassEdge->GetYaxis()->SetTitle("rel. diff. inv. mass bkg");
                UncertPerMassEdge->GetYaxis()->SetRangeUser(0.,1.);
                if(ij == 0) UncertPerMassEdge->Draw("AP");
                else UncertPerMassEdge->Draw("P same");

               legend->AddEntry(UncertPerMassEdge,Form("%s GeV/c^{2}",path.Data()),"p");
        }
 legend->Draw();
 TCanvas *cSystMass = new TCanvas();
 grSyst->SetTitle(" ");
 grSyst->GetXaxis()->SetTitle("p_{T}");
 grSyst->GetYaxis()->SetTitle("average rel. diff. inv. mass bkg");
 grSyst->SetLineColor(1);
 grSyst->SetLineWidth(2);
 grSyst->Draw("APL");

 return;
}


TGraphErrors* computeRMSandAveragesForMassEdges(Int_t massHyp, Bool_t draw){
TString nameGraphs[] = {"fbVsPtCompLS","fbVsPtCompME"};
TString nameGraphs2[] = {"fbVsPtCombinedQuadrWeightsLS","fbVsPtCombinedQuadrWeightsME"};
//TString nameGraphs[] = {"fbVsPtComp","fbVsPtCompME","fbVsPtCompFsigFixed","fbVsPtCompFsigFixedME"};
//TString nameGraphs[] = {"fbVsPtComp","fbVsPtCompME","fbVsPtCombinedQuadrWeights","fbVsPtCombinedQuadrWeightsME"};
TString titlePtString[] = {"1.5-10","1.5-3","3-5","5-10"};
TH1F *fbVal[nptbins+1];
    TGraphErrors *grAverage = 0x0;
    TGraphErrors *gRelativeUncertainties = 0x0;
        TString path = "mass_";
        for(int i=0; i<nmass-1; i++) path.Append(Form("%1.2f_",mbins[massHyp][i]));
        path.Append(Form("%1.2f",mbins[massHyp][nmass-1]));

        //TString path = Form("mass_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f",mbins[0][0],mbins[0][1],mbins[0][2],mbins[0][3],mbins[0][4],mbins[0][5]);
        TFile *fin = new TFile(Form("%s/fbResultsComb.root",path.Data()),"READ");
        TFile *finquadratic = new TFile(Form("%s/fbResultsCombQuadrWeights.root",path.Data()),"READ");
        
        Double_t xx, yy;
for(int ipt=0; ipt<nptbins+1; ipt++) fbVal[ipt] = new TH1F(Form("Pt bin:%s",titlePtString[ipt].Data()), "", 100, 0.,1.);

            Double_t minValue[4] = {100,100,100,100};
            Double_t maxValue[4] = {0.,0.,0.,0.};
         for(int imeth=0;imeth<2;imeth++){
             TGraphErrors *gr = (TGraphErrors*)fin->Get(nameGraphs[imeth]);
              TGraphErrors *grquadratic = (TGraphErrors*)finquadratic->Get(nameGraphs2[imeth]);
        if(imeth==0) grAverage = (TGraphErrors*)gr->Clone("grAverage");
         if(imeth==0) gRelativeUncertainties = (TGraphErrors*)gr->Clone("gRelativeUncertainties");
         for(int ipt=0; ipt<nptbins+1;ipt++){
            
             gr->GetPoint(ipt, xx, yy);
             if(yy < minValue[ipt]) minValue[ipt] = yy;
             if(yy > maxValue[ipt]) maxValue[ipt] = yy;
             fbVal[ipt]->Fill(yy);
             grquadratic->GetPoint(ipt, xx, yy);
             fbVal[ipt]->Fill(yy);
	    }
         }

	for(int ipt=0; ipt<nptbins+1; ipt++){
		grAverage->GetPoint(ipt,xx,yy);
		grAverage->SetPoint(ipt,xx, fbVal[ipt]->GetMean());
		grAverage->SetPointError(ipt,grAverage->GetErrorX(ipt), maxValue[ipt]-minValue[ipt]);
		gRelativeUncertainties->SetPoint(ipt,xx,1);
		gRelativeUncertainties->SetPointError(ipt,grAverage->GetErrorX(ipt),(maxValue[ipt]-minValue[ipt])/fbVal[ipt]->GetMean());
	}

if(draw){
TCanvas *cfbVal = new TCanvas();
cfbVal->Divide(2,2);
for(int ipt=0; ipt<nptbins+1; ipt++){
cfbVal->cd(ipt+1);
fbVal[ipt]->GetXaxis()->SetTitle("f_{B}");
fbVal[ipt]->GetYaxis()->SetTitle("Entries");
fbVal[ipt]->DrawCopy();
}
cfbVal->Draw();
TCanvas *cAvRMS = new TCanvas();
cAvRMS->cd();
gRelativeUncertainties->GetYaxis()->SetTitle("Rel. syst. uncertainty - inv. mass bkg (%)");
gRelativeUncertainties->Draw();
cAvRMS->Draw();

TCanvas *cfbVsPtUnc = new TCanvas();
grAverage->GetYaxis()->SetTitle("Average +- RMS");
grAverage->Draw();
cAvRMS->SaveAs("ComparisonGraphs/RelativUncertaintiesOnInvMassBkg.pdf");
}

return gRelativeUncertainties;
}

void computeRMSandAverages(){
    TString methodName="fbVsPtCompME";
    TString methodName2="fbVsPtCombinedQuadrWeightsME";
    TH1F *fbVal[nptbins+1]; 
    TString titlePtString[] = {"1.5-10","1.5-3","3-5","5-10"};
    

    
    TGraphErrors *grAverage = 0x0;
    TGraphErrors *gRelativeUncertainties = 0x0;
    TGraphErrors *gr[nmassHyp];
    TGraphErrors *grquadratic[nmassHyp];
for(int ipt=0; ipt<nptbins+1; ipt++) fbVal[ipt] = new TH1F(Form("Pt bin:%s",titlePtString[ipt].Data()), "", 100, 0.,1.);
    Double_t xx, yy;
for(int imass=0;imass<nmassHyp;imass++){
    TString path = "mass_";
    for(int i=0; i<nmass-1; i++) path.Append(Form("%1.2f_",mbins[imass][i]));
    path.Append(Form("%1.2f",mbins[imass][nmass-1]));
    //TString path = Form("mass_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f",mbins[imass][0],mbins[imass][1],mbins[imass][2],mbins[imass][3],mbins[imass][4],mbins[imass][5]);
         TFile *fin = new TFile(Form("%s/fbResultsComb.root",path.Data()),"READ");
         gr[imass] = (TGraphErrors*)fin->Get(methodName);
         if(imass==0) grAverage = (TGraphErrors*)gr[imass]->Clone("grAverage");
         if(imass==0) gRelativeUncertainties = (TGraphErrors*)gr[imass]->Clone("gRelativeUncertainties");
         TFile *finquadratic = new TFile(Form("%s/fbResultsCombQuadrWeights.root",path.Data()),"READ");
         grquadratic[imass] = (TGraphErrors*)finquadratic->Get(methodName2);
for(int ipt=0; ipt<nptbins+1;ipt++){
gr[imass]->GetPoint(ipt, xx, yy);
fbVal[ipt]->Fill(yy);
grquadratic[imass]->GetPoint(ipt, xx, yy);
fbVal[ipt]->Fill(yy);
}
}

TCanvas *cfbVal = new TCanvas();
cfbVal->Divide(2,2);
for(int ipt=0; ipt<nptbins+1; ipt++){
cfbVal->cd(ipt+1);
fbVal[ipt]->GetXaxis()->SetTitle("f_{B}");
    fbVal[ipt]->GetYaxis()->SetTitle("Entries");
fbVal[ipt]->DrawCopy();
grAverage->GetPoint(ipt,xx,yy);
grAverage->SetPoint(ipt,xx, fbVal[ipt]->GetMean());
grAverage->SetPointError(ipt,grAverage->GetErrorX(ipt), fbVal[ipt]->GetRMS());
gRelativeUncertainties->SetPoint(ipt,xx,1);
gRelativeUncertainties->SetPointError(ipt,grAverage->GetErrorX(ipt),fbVal[ipt]->GetRMS()/fbVal[ipt]->GetMean());
}
cfbVal->Draw();
TCanvas *cAvRMS = new TCanvas();
cAvRMS->cd();
grAverage->GetYaxis()->SetTitle("F_{B} (Average +- RMS)");
gRelativeUncertainties->GetYaxis()->SetTitle("Rel. syst. uncertainty on x-bkg (%)");
gRelativeUncertainties->Draw();
cAvRMS->Draw();

TCanvas *cfbVsPtUnc = new TCanvas();
grAverage->Draw();

cAvRMS->SaveAs("ComparisonGraphs/RelativUncertaintiesOnXBkg.pdf");
cfbVal->SaveAs("ComparisonGraphs/FbValuesXBkg.pdf");
cfbVsPtUnc->SaveAs("ComparisonGraphs/AverageFbWithUncertaintyXBkg.pdf");

//Compute averages and uncertainties
Double_t averageVal[nptbins+1];
Double_t uncertaintyVal[nptbins+1];
Double_t nTries = ((Double_t)nmass)*2;
printf("%f",nTries);
for(int ipt = 0; ipt<nptbins+1; ipt++){
    Double_t y1,y2;
    averageVal[ipt] = 0.;
    uncertaintyVal[ipt] = 0.;
    for(int imass = 0; imass < nmass; imass++) {
        gr[imass]->GetPoint(ipt,xx,y1);
        grquadratic[imass]->GetPoint(ipt,xx,y2);
        averageVal[ipt] += (1/nTries)*(y1+y2);
    for(int jmass = 0; jmass < nmass; jmass++){
        uncertaintyVal[ipt] +=  (1/nTries)*(1/nTries)*gr[imass]->GetErrorY(ipt)*grquadratic[jmass]->GetErrorY(ipt)+(1/nTries)*(1/nTries)*grquadratic[imass]->GetErrorY(ipt)*grquadratic[jmass]->GetErrorY(ipt)+(1/nTries)*(1/nTries)*gr[imass]->GetErrorY(ipt)*gr[jmass]->GetErrorY(ipt);
        printf("Uncertainty: %1.5f \n",uncertaintyVal[ipt]);
        
}}
    uncertaintyVal[ipt] = TMath::Sqrt(uncertaintyVal[ipt]);
    
}


for(int i = 0; i<4; i++){
 printf("Val: %1.5f and uncert: %1.5f \n\n",averageVal[i],uncertaintyVal[i]);
    }

TCanvas *c1 = new TCanvas("c1","A simple",200,10,700,500);
Double_t xVal[4] = {0.725,2.25,4.,7.5};
Double_t xValErr[4] = {0.725,0.725,1.,2.5};
c1->SetFillColor(42);
c1->SetGrid();
c1->GetFrame()->SetFillColor(21);
c1->GetFrame()->SetBorderSize(12);
TGraphErrors *grc1 = new TGraphErrors(4,xVal,averageVal,xValErr,uncertaintyVal);
grc1->SetMarkerColor(4);
grc1->SetMarkerStyle(21);
grc1->SetTitle("Averages with uncertainties");
grc1->GetXaxis()->SetTitle("p_{T}(GeV/c)");
grc1->GetYaxis()->SetTitle("f_{B}");
grc1->GetYaxis()->SetRangeUser(0.,0.5);
grc1->Draw("ALP");
TGraphErrors *grc1syst = new TGraphErrors(4,xVal,averageVal,xValErr,uncertaintyVal);
for(int ipt=0; ipt<4; ipt++)
grc1syst->SetPointError(ipt, gRelativeUncertainties->GetErrorX(ipt),gRelativeUncertainties->GetErrorY(ipt)*averageVal[ipt]);
grc1syst->SetFillStyle(0);
grc1syst->Draw(" E2 SAME");
/////
	TLegend *legend=new TLegend(0.37,0.10,0.62,0.30);
        legend->SetBorderSize(0); legend->SetFillColor(0); legend->SetTextFont(42);
        legend->SetFillStyle(0); legend->SetMargin(0.25);
        legend->SetEntrySeparation(0.15);
        legend->SetTextSize(0.045);
 	legend->AddEntry(grc1,"stat. uncertainty","pel");
 	legend->AddEntry(grc1syst,"syst. uncertainty from x-bkg","f");
	legend->Draw();
c1->SaveAs("ComparisonGraphs/AveragesAndUncertainties.pdf");
}

void plotVsAllMassHypothesis(){
 
    for(int i=0;i<nmassHyp;i++) plotVsMassHypothesis(i);
    return;
        
        
}
void plotAllResults(Bool_t createPdf=kFALSE){
 for(int i=0; i<mbinsSize; i++) plotResults(i,createPdf);
 return;
}

void plotResults(Int_t num, Bool_t createPdf=kFALSE){
    printf("%d",mbinsSize);
	readFbResults(mbins[num],createPdf);
/////
}

void plotVsAllMassEdges(){
 
    
    TString nameGraphs[] = {"fbVsPtComp","fbVsPtCompME","fbVsPtCompFsigFixed","fbVsPtCompFsigFixedME"};
    TString nameGraphs2[] = {"fbVsPtCombinedQuadrWeights","fbVsPtCombinedQuadrWeightsFsigFixed","fbVsPtCombinedQuadrWeightsFsigFixedME","fbVsPtCombinedQuadrWeightsME"};
    for(int i=0; i<4; i++){
     plotVsMassEdges(nameGraphs[i],kFALSE,0.,0.5);   
     plotVsMassEdges(nameGraphs2[i],kTRUE,0.,0.5);  
    };
}


void plotChi2andFbVsAllMassEdgesME(){


    TString nameGraphs[] = {"fbVsPtCompME","Chi2XprojVsPtCompME","Chi2MprojVsPtCompME"};
    Double_t ymin[] = {0., 0., 0.};
    Double_t ymax[] = {0.5, 5., 5.};
    for(int i=0; i<3; i++){
     plotVsMassEdges(nameGraphs[i],kFALSE,ymin[i],ymax[i]);
    };
  return;
}


void plotVsMassHypothesis(Int_t num){
        TCanvas *canv = new TCanvas();
        gSystem->Exec(Form("mkdir -p ComparisonGraphs"));
            TLegend *legend=new TLegend(0.37,0.10,0.62,0.30);
  	legend->SetBorderSize(0); legend->SetFillColor(0); legend->SetTextFont(42);
  	legend->SetFillStyle(0); legend->SetMargin(0.25);
  	legend->SetEntrySeparation(0.15);
  	legend->SetTextSize(0.045);
       //TString path = Form("mass_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f",mbins[0][0],mbins[0][1],mbins[0][2],mbins[0][3],mbins[0][4],mbins[0][5]);
        TString path = "mass_";
        for(int i=0; i<nmass-1; i++) path.Append(Form("%1.2f_",mbins[num][i]));
        path.Append(Form("%1.2f",mbins[num][nmass-1]));
	 TString pathLegend = path;
    pathLegend.ReplaceAll("_"," ");
    legend->SetHeader(Form("%s",pathLegend.Data()));
        
	Int_t colors[] = {1,2,4,6,7};
	Int_t styles[] = {20,24,21,32,28};
    TString nameGraphs[] = {"fbVsPtComp","fbVsPtCompME","fbVsPtCompFsigFixed","fbVsPtCompFsigFixedME"};
        TFile *fin = new TFile(Form("%s/fbResultsComb.root",path.Data()),"READ");
        fin->ls();
        for(int i=0; i<4; i++){
            TGraphErrors *gr = (TGraphErrors*)fin->Get(nameGraphs[i].Data());
            	gr->SetLineColor(colors[i]);
            gr->SetMarkerColor(colors[i]);
            gr->SetMarkerStyle(styles[i]);
            gr->GetYaxis()->SetRangeUser(0,0.3);
            legend->AddEntry(gr,Form("%s",nameGraphs[i].Data()),"p"); 
            if(i==0) gr->Draw("AP");
            else gr->Draw("P SAME");
            
        }
        fin->Close();
        legend->Draw();
        canv->SaveAs(Form("ComparisonGraphs/%s.pdf",path.Data()));
        return;
    
}

void plotVsMassEdges(TString nameGraph, Bool_t isQuadr, Double_t ymin, Double_t ymax){
  
        TString quadrString = "";
        if(isQuadr) quadrString.Append("QuadrWeights");
gSystem->Exec(Form("mkdir -p ComparisonGraphs"));
        TLegend *legend=new TLegend(0.37,0.05,0.62,0.45);
  	legend->SetBorderSize(0); legend->SetFillColor(0); legend->SetTextFont(42);
  	legend->SetFillStyle(0); legend->SetMargin(0.25);
  	legend->SetEntrySeparation(0.15);
  	legend->SetTextSize(0.045);
    legend->SetHeader(Form("%s",nameGraph.Data()));

	Int_t colors[] = {1,2,4,6,7,8,9,3,15};
        Int_t styles[] = {20,24,21,32,28,31,33,35,22};
        TCanvas *fbout = new TCanvas();
	for(int im=0; im<nmassHyp; im++){
        TString path = "mass_";
        for(int i=0; i<nmass-1; i++) path.Append(Form("%1.2f_",mbins[im][i]));
        path.Append(Form("%1.2f",mbins[im][nmass-1]));
        
        TFile *fin = new TFile(Form("%s/fbResultsComb%s.root",path.Data(),quadrString.Data()),"READ");
        fin->ls();
        TGraphErrors *gr = (TGraphErrors*)fin->Get(nameGraph);
	gr->SetLineColor(colors[im]);
	gr->SetMarkerColor(colors[im]);
	gr->SetMarkerStyle(styles[im]);
    	gr->GetYaxis()->SetRangeUser(ymin,ymax);
	path.ReplaceAll("_"," ");
  	legend->AddEntry(gr,Form("%s GeV/c^{2}",path.Data()),"p"); 
        if(im==0) gr->Draw("AP");
	else gr->Draw("P SAME"); 
	fin->Close();
	}
	legend->Draw();
    fbout->SaveAs(Form("ComparisonGraphs/%s%s.pdf",nameGraph.Data(),quadrString.Data()));
	return;
}

void readFbResults(Double_t massBins[], Bool_t pdf){

TString path = "mass_";
for(int i=0; i<nmass-1; i++) path.Append(Form("%1.2f_",massBins[i]));
path.Append(Form("%1.2f",massBins[nmass-1]));

TGraphErrors *fbVsPt =  BuildGraphFb(massBins,"FF_FS",kFALSE,kFALSE,kFALSE,1,20);
TGraphErrors *fbVsPtFsigFixed = BuildGraphFb(massBins,"FF_FS",kFALSE,kTRUE,kFALSE,4,20);
TGraphErrors *fbVsPtFF =  BuildGraphFb(massBins,"FF",kFALSE,kFALSE,kFALSE,2,24);
TGraphErrors *fbVsPtFFFsigFixed =  BuildGraphFb(massBins,"FF",kFALSE,kTRUE,kFALSE,6,24);
TGraphErrors *fbVsPtQuadr = BuildGraphFb(massBins,"FF_FS",kTRUE,kFALSE,kFALSE,1,21);
TGraphErrors *fbVsPtFFQuadr = BuildGraphFb(massBins,"FF",kTRUE,kFALSE,kFALSE,2,21);
TGraphErrors *fbVsPtFsigFixedQuadr = BuildGraphFb(massBins,"FF_FS",kTRUE,kTRUE,kFALSE,4,21);
TGraphErrors *fbVsPtFFFsigFixedQuadr = BuildGraphFb(massBins,"FF",kTRUE,kTRUE,kFALSE,6,21);
    
    
TGraphErrors *fbVsPtME =  BuildGraphFb(massBins,"FF_FS",kFALSE,kFALSE,kTRUE,1,20);
TGraphErrors *fbVsPtFsigFixedME = BuildGraphFb(massBins,"FF_FS",kFALSE,kTRUE,kTRUE,4,20);
TGraphErrors *fbVsPtFFME =  BuildGraphFb(massBins,"FF",kFALSE,kFALSE,kTRUE,2,24);
TGraphErrors *fbVsPtFFFsigFixedME =  BuildGraphFb(massBins,"FF",kFALSE,kTRUE,kTRUE,6,24);
TGraphErrors *fbVsPtQuadrME = BuildGraphFb(massBins,"FF_FS",kTRUE,kFALSE,kTRUE,1,21);
TGraphErrors *fbVsPtFFQuadrME = BuildGraphFb(massBins,"FF",kTRUE,kFALSE,kTRUE,2,21);
TGraphErrors *fbVsPtFsigFixedQuadrME = BuildGraphFb(massBins,"FF_FS",kTRUE,kTRUE,kTRUE,4,21);
TGraphErrors *fbVsPtFFFsigFixedQuadrME = BuildGraphFb(massBins,"FF",kTRUE,kTRUE,kTRUE,6,21);
//
TGraphErrors *fbVsPtComp = BuildGraphFbCompTypes(massBins,types,kFALSE,kFALSE,0,"fbVsPtComp",1,20);
TGraphErrors *fbVsPtCompFsigFixed = BuildGraphFbCompTypes(massBins,types,kFALSE,kTRUE,0,"fbVsPtCompFsigFixed",2,24);
TGraphErrors *fbVsPtCompME = BuildGraphFbCompTypes(massBins,types,kFALSE,kFALSE,1,"fbVsPtCompME",4,21);
TGraphErrors *fbVsPtCompFsigFixedME = BuildGraphFbCompTypes(massBins,types,kFALSE,kTRUE,1,"fbVsPtCompFsigFixedME",6,5);
TGraphErrors *fbVsPtCompLS = BuildGraphFbCompTypes(massBins,types,kFALSE,kFALSE,2,"fbVsPtCompLS",8,26);
TGraphErrors *fbVsPtCompFsigFixedLS = BuildGraphFbCompTypes(massBins,types,kFALSE,kTRUE,2,"fbVsPtCompFsigFixedLS",7,32);
TGraphErrors *fbVsPtCompQuadr = BuildGraphFbCompTypes(massBins,types,kTRUE,kFALSE,0,"fbVsPtCombinedQuadrWeights",1,20);
TGraphErrors *fbVsPtCompQuadrFsigFixed = BuildGraphFbCompTypes(massBins,types,kTRUE,kTRUE,0,"fbVsPtCombinedQuadrWeightsFsigFixed",2,24);
TGraphErrors *fbVsPtCompQuadrFsigFixedME = BuildGraphFbCompTypes(massBins,types,kTRUE,kTRUE,1,"fbVsPtCombinedQuadrWeightsFsigFixedME",4,21);
TGraphErrors *fbVsPtCompQuadrME = BuildGraphFbCompTypes(massBins,types,kTRUE,kFALSE,1,"fbVsPtCombinedQuadrWeightsME",6,5);
TGraphErrors *fbVsPtCompQuadrFsigFixedLS = BuildGraphFbCompTypes(massBins,types,kTRUE,kTRUE,1,"fbVsPtCombinedQuadrWeightsFsigFixedLS",8,26);
TGraphErrors *fbVsPtCompQuadrLS = BuildGraphFbCompTypes(massBins,types,kTRUE,kFALSE,1,"fbVsPtCombinedQuadrWeightsLS",7,32);
///
TGraphErrors *Chi2XVsPt = BuildGraphChi2X(massBins,"FF_FS",kFALSE,kFALSE,1,20); 
TGraphErrors *Chi2XVsPtFF = BuildGraphChi2X(massBins,"FF",kFALSE,kFALSE,2,20); 
TGraphErrors *Chi2XVsPtFsigFixed = BuildGraphChi2X(massBins,"FF_FS",kFALSE,kTRUE,1,20);
TGraphErrors *Chi2XVsPtFFFsigFixed = BuildGraphChi2X(massBins,"FF",kFALSE,kTRUE,2,20);
TGraphErrors *Chi2MVsPt = BuildGraphChi2M(massBins,"FF_FS",kFALSE,kFALSE,1,20); 
TGraphErrors *Chi2MVsPtFF = BuildGraphChi2M(massBins,"FF",kFALSE,kFALSE,2,20);
TGraphErrors *Chi2MVsPtFsigFixed = BuildGraphChi2M(massBins,"FF_FS",kFALSE,kTRUE,1,20);
TGraphErrors *Chi2MVsPtFFFsigFixed = BuildGraphChi2M(massBins,"FF",kFALSE,kTRUE,2,20);

TGraphErrors *Chi2XVsPtCombME = BuildGraphChi2XProjCompTypes(massBins,types,kFALSE,kFALSE,kTRUE,"Chi2XprojVsPtCompME",4,21);
TGraphErrors *Chi2XVsPtQuadrCombME = BuildGraphChi2XProjCompTypes(massBins,types,kTRUE,kFALSE,kTRUE,"Chi2XprojVsPtQuadrCompME",4,21);
TGraphErrors *Chi2MVsPtCombME = BuildGraphChi2MProjCompTypes(massBins,types,kFALSE,kFALSE,kTRUE,"Chi2MprojVsPtCompME",4,21);
TGraphErrors *Chi2MVsPtQuadrCombME = BuildGraphChi2MProjCompTypes(massBins,types,kTRUE,kFALSE,kTRUE,"Chi2MprojVsPtQuadrCompME",4,21);
///
TGraphErrors *fSigVsPt = BuildGraphFSig(massBins,"FF_FS",kFALSE,kFALSE,1,24); 
TGraphErrors *fSigVsPtFF = BuildGraphFSig(massBins,"FF",kFALSE,kFALSE,1,24); 
TGraphErrors *fSigVsPt1D = BuildGraphFSig1D(massBins,"FF_FS",2,20); 
TGraphErrors *fSigVsPt1DFF = BuildGraphFSig1D(massBins,"FF",2,20);


TString massfilename="mass_";
TString massfilename1="";
for(int ij=0;ij<nmass; ij++) { 
	if(ij != nmass-1) { massfilename.Append(Form("%1.2f_",massBins[ij])); massfilename1.Append(Form("m%d%1.2f_",ij+1,massBins[ij]));  }
	else {massfilename.Append(Form("%1.2f",massBins[ij])); massfilename1.Append(Form("m%d%1.2f",ij+1,massBins[ij]));}
}



TCanvas *resFb = new TCanvas("fbvsPt","",1200,600);
resFb->Divide(2,1);
resFb->cd(1);
fbVsPt->GetYaxis()->SetRangeUser(0.,0.5);
fbVsPt->Draw("AP");
fbVsPtFF->Draw("P SAME");
fbVsPtQuadr->Draw("P SAME");
fbVsPtFFQuadr->Draw("P SAME");
TLegend *legend=new TLegend(0.17,0.72,0.42,0.88);
  legend->SetBorderSize(0); legend->SetFillColor(0); legend->SetTextFont(42);
  legend->SetFillStyle(0); legend->SetMargin(0.25);
  legend->SetEntrySeparation(0.15);
  legend->SetTextSize(0.045);
  legend->SetHeader("Fsig free in the likelihood fit");
  legend->AddEntry(fbVsPt,"FF+FS","p");
  legend->AddEntry(fbVsPtFF,"FF","p"); 
  legend->AddEntry(fbVsPtQuadr,"FF+FS (quadratic weights)","p");
  legend->AddEntry(fbVsPtFFQuadr,"FF (quadratic weights)","p"); 
  legend->Draw();
resFb->cd(2);
fbVsPtFsigFixed->GetYaxis()->SetRangeUser(0.,0.5);
fbVsPtFsigFixed->Draw("AP");
fbVsPtFFFsigFixed->Draw("P SAME");
fbVsPtFsigFixedQuadr->Draw("P SAME");
fbVsPtFFFsigFixedQuadr->Draw("P SAME");

    

TLegend *legend1=new TLegend(0.17,0.72,0.42,0.88);
  legend1->SetBorderSize(0); legend1->SetFillColor(0); legend1->SetTextFont(42);
  legend1->SetFillStyle(0); legend1->SetMargin(0.25);
  legend1->SetEntrySeparation(0.15);
  legend1->SetTextSize(0.045);
  legend1->SetHeader("Fsig fixed in the likelihood fit");
  legend1->AddEntry(fbVsPtFsigFixed,"FF+FS","p");
  legend1->AddEntry(fbVsPtFFFsigFixed,"FF","p"); 
  legend1->AddEntry(fbVsPtFsigFixedQuadr,"FF+FS","p");
  legend1->AddEntry(fbVsPtFFFsigFixedQuadr,"FF","p"); 
  legend1->Draw();

resFb->SaveAs(Form("%s/fbResults.png",path.Data()));

    TCanvas *resFbME = new TCanvas("fbvsPtME","",1200,600);
    resFbME->Divide(2,1);
    resFbME->cd(1);
    fbVsPtME->GetYaxis()->SetRangeUser(0.,0.5);
    fbVsPtME->Draw("AP");
    fbVsPtFFME->Draw("P SAME");
    fbVsPtQuadrME->Draw("P SAME");
    fbVsPtFFQuadrME->Draw("P SAME");
      legend->Draw();
    resFbME->cd(2);
    fbVsPtFsigFixedME->GetYaxis()->SetRangeUser(0.,0.5);
    fbVsPtFsigFixedME->Draw("AP");
    fbVsPtFFFsigFixedME->Draw("P SAME");
    fbVsPtFsigFixedQuadrME->Draw("P SAME");
    fbVsPtFFFsigFixedQuadrME->Draw("P SAME");
    legend1->Draw();
    
    resFbME->SaveAs(Form("%s/fbResultsME.png",path.Data()));
    
TCanvas *resFsig = new TCanvas("fSigVsPt","fSigVsPt",1200,600); 
resFsig->Divide(2,1);
resFsig->cd(1);
fSigVsPt->SetTitle(" ");
fSigVsPt->GetXaxis()->SetTitle("p_{T}(GeV/c)");
fSigVsPt->GetYaxis()->SetTitle("f_{Sig}");
fSigVsPt->GetYaxis()->SetRangeUser(0.,0.2);
fSigVsPt->Draw("AP");
fSigVsPt1D->Draw("P SAME");

resFsig->cd(2);
fSigVsPt1DFF->GetYaxis()->SetRangeUser(0.,0.2);
fSigVsPt1DFF->Draw("AP");
fSigVsPtFF->Draw("P SAME");

TLegend *legend2=new TLegend(0.17,0.72,0.42,0.88);
  legend2->SetBorderSize(0); legend2->SetFillColor(0); legend2->SetTextFont(42);
  legend2->SetFillStyle(0); legend2->SetMargin(0.25);
  legend2->SetEntrySeparation(0.15);
  legend2->SetTextSize(0.045);
  legend2->SetHeader("FF candidates");
  legend2->AddEntry(fSigVsPtFF,"Fsig from likelihood fit","p");
  legend2->AddEntry(fSigVsPt1DFF,"Fsig from 1D fit","p"); 
  legend2->Draw();
resFsig->cd(1);
TLegend *legend3=new TLegend(0.17,0.72,0.42,0.88);
  legend3->SetBorderSize(0); legend3->SetFillColor(0); legend3->SetTextFont(42);
  legend3->SetFillStyle(0); legend3->SetMargin(0.25);
  legend3->SetEntrySeparation(0.15);
  legend3->SetTextSize(0.045);
  legend3->SetHeader("FF+FS candidates");
  legend3->AddEntry(fSigVsPtFF,"Fsig from likelihood fit","p");
  legend3->AddEntry(fSigVsPt1DFF,"Fsig from 1D fit","p"); 
  legend3->Draw();

resFsig->SaveAs(Form("%s/FsigResults.png",path.Data()));
 ///
 TCanvas *chi2Like = new TCanvas("chi2","chi2",1200,1200);
 chi2Like->Divide(2,2);
 chi2Like->cd(1);
 Chi2XVsPt->GetYaxis()->SetRangeUser(0.,5.);
 Chi2XVsPt->Draw("AP");
 Chi2XVsPtFF->Draw("P SAME");
 chi2Like->cd(2);
 Chi2MVsPt->GetYaxis()->SetRangeUser(0.,5.);
 Chi2MVsPt->Draw("AP");
 Chi2MVsPtFF->Draw("P SAME");

 chi2Like->cd(1);
  TLegend *legend4=new TLegend(0.17,0.72,0.42,0.88);
  legend4->SetBorderSize(0); legend4->SetFillColor(0); legend4->SetTextFont(42);
  legend4->SetFillStyle(0); legend4->SetMargin(0.25);
  legend4->SetEntrySeparation(0.15);
  legend4->SetTextSize(0.045);
  legend4->SetHeader("Projections on X");
  legend4->AddEntry(Chi2XVsPt,"FF+FS","p");
  legend4->AddEntry(Chi2XVsPtFF,"FF","p"); 
  legend4->Draw();

 chi2Like->cd(2);
  TLegend *legend5=new TLegend(0.17,0.72,0.42,0.88);
  legend5->SetBorderSize(0); legend5->SetFillColor(0); legend5->SetTextFont(42);
  legend5->SetFillStyle(0); legend5->SetMargin(0.25);
  legend5->SetEntrySeparation(0.15);
  legend5->SetTextSize(0.045);
  legend5->SetHeader("Projections on M");
  legend5->AddEntry(Chi2MVsPt,"FF+FS","p");
  legend5->AddEntry(Chi2MVsPtFF,"FF","p"); 
  legend5->Draw();

chi2Like->cd(3);
 Chi2XVsPtFsigFixed->GetYaxis()->SetRangeUser(0.,5.);
 Chi2XVsPtFsigFixed->Draw("AP");
 Chi2XVsPtFFFsigFixed->Draw("P SAME");


TLegend *legend6=new TLegend(0.17,0.72,0.42,0.88);
  legend6->SetBorderSize(0); legend6->SetFillColor(0); legend6->SetTextFont(42);
  legend6->SetFillStyle(0); legend6->SetMargin(0.25);
  legend6->SetEntrySeparation(0.15);
  legend6->SetTextSize(0.045);
  legend6->SetHeader("Projections on X (Fsig fixed)");
  legend6->AddEntry(Chi2XVsPtFsigFixed,"FF+FS","p");
  legend6->AddEntry(Chi2XVsPtFFFsigFixed,"FF","p");
  legend6->Draw();

 
  chi2Like->cd(4);
 Chi2MVsPtFsigFixed->GetYaxis()->SetRangeUser(0.,5.);
 Chi2MVsPtFsigFixed->Draw("AP");
 Chi2MVsPtFFFsigFixed->Draw("P SAME");


TLegend *legend7=new TLegend(0.17,0.72,0.42,0.88);
  legend7->SetBorderSize(0); legend7->SetFillColor(0); legend7->SetTextFont(42);
  legend7->SetFillStyle(0); legend7->SetMargin(0.25);
  legend7->SetEntrySeparation(0.15);
  legend7->SetTextSize(0.045);
  legend7->SetHeader("Projections on M (Fsig fixed)");
  legend7->AddEntry(Chi2MVsPtFsigFixed,"FF+FS","p");
  legend7->AddEntry(Chi2MVsPtFFFsigFixed,"FF","p");
  legend7->Draw();

chi2Like->SaveAs(Form("%s/Chi2LikeFitResults.png",path.Data()));

TCanvas *cFbcomp = new TCanvas();
    fbVsPtComp->GetYaxis()->SetRangeUser(0.,0.5);
fbVsPtComp->Draw("AP");
fbVsPtCompFsigFixed->Draw("P SAME");
    fbVsPtCompME->Draw("P SAME");
    fbVsPtCompFsigFixedME->Draw("P SAME");

    TLegend *legend8=new TLegend(0.17,0.72,0.42,0.88);
      legend8->SetBorderSize(0); legend8->SetFillColor(0); legend8->SetTextFont(42);
      legend8->SetFillStyle(0); legend8->SetMargin(0.25);
      legend8->SetEntrySeparation(0.15);
      legend8->SetTextSize(0.045);
    legend8->SetHeader("Invariant mass fit vs mixed event");
    legend8->AddEntry(fbVsPtComp,"bkg mass from fit","pe");
    legend8->AddEntry(fbVsPtCompFsigFixed,"bkg mass from fit, F_{sig} fixed","pe");
    legend8->AddEntry(fbVsPtCompME,"bkg mass from ME","pe");
    legend8->AddEntry(fbVsPtCompFsigFixedME,"bkg mass from ME, F_{sig} fixed","pe");
    legend8->Draw();
    
    TFile fout(Form("%s/fbResultsComb.root",massfilename.Data()),"RECREATE");
fbVsPtComp->Write();
fbVsPtCompME->Write();
fbVsPtCompLS->Write();
fbVsPtCompFsigFixed->Write();
fbVsPtCompFsigFixedME->Write();
fbVsPtCompFsigFixedLS->Write();
///
Chi2XVsPtCombME->Write();    
Chi2XVsPtQuadrCombME->Write();    
Chi2MVsPtCombME->Write();    
Chi2MVsPtQuadrCombME->Write();    
fout.Close();

TCanvas *cFbcompQuadratic = new TCanvas();
fbVsPtCompQuadr->GetYaxis()->SetRangeUser(0.,0.5);
fbVsPtCompQuadr->Draw("AP");
fbVsPtCompQuadrFsigFixed->Draw("P SAME");
    fbVsPtCompQuadrME->Draw("P SAME");
    fbVsPtCompQuadrFsigFixedME->Draw("P SAME");
    
    TLegend *legend9=new TLegend(0.17,0.72,0.42,0.88);
      legend9->SetBorderSize(0); legend9->SetFillColor(0); legend9->SetTextFont(42);
      legend9->SetFillStyle(0); legend9->SetMargin(0.25);
      legend9->SetEntrySeparation(0.15);
      legend9->SetTextSize(0.045);
    legend9->SetHeader("Invariant mass fit vs mixed event quadratic");
    legend9->AddEntry(fbVsPtCompQuadr,"bkg mass from fit","pe");
    legend9->AddEntry(fbVsPtCompQuadrFsigFixed,"bkg mass from fit, F_{sig} fixed","pe");
    legend9->AddEntry(fbVsPtCompQuadrME,"bkg mass from ME","pe");
    legend9->AddEntry(fbVsPtCompQuadrFsigFixedME,"bkg mass from ME, F_{sig} fixed","pe");
    legend9->Draw();
    
    TFile fout2(Form("%s/fbResultsCombQuadrWeights.root",massfilename.Data()),"RECREATE");
    fbVsPtCompQuadr->Write();
fbVsPtCompQuadrME->Write();
fbVsPtCompQuadrLS->Write();
fbVsPtCompQuadrFsigFixed->Write();
fbVsPtCompQuadrFsigFixedME->Write();
fbVsPtCompQuadrFsigFixedLS->Write();
fout2.Close();
    

/// create pdf
if(pdf){
gSystem->Exec(Form("cp showplots.tex %s",massfilename.Data()));
gSystem->cd(Form("%s",massfilename.Data()));
    gSystem->Exec(Form("pdflatex \"\\def\\mLOne{%1.1f}\\def\\mLTwo{%1.1f}\\def\\mLThree{%1.1f}\\def\\mHOne{%1.1f}\\def\\mHTwo{%1.1f}\\def\\mHThree{%1.1f}\\include{showplots}\"",massBins[0],massBins[1],massBins[2],massBins[3],massBins[4],massBins[5]));
        gSystem->Exec(Form("pdflatex \"\\def\\mLOne{%1.1f}\\def\\mLTwo{%1.1f}\\def\\mLThree{%1.1f}\\def\\mHOne{%1.1f}\\def\\mHTwo{%1.1f}\\def\\mHThree{%1.1f}\\include{showplots}\"",massBins[0],massBins[1],massBins[2],massBins[3],massBins[4],massBins[5]));
gSystem->Exec("mv showplots.pdf summaryPlots.pdf");
gSystem->Exec("rm -rf showplots.*");
    

/*
/// create pdf
if(pdf){
gSystem->Exec(Form("cp showplots.tex %s",massfilename.Data()));
gSystem->cd(Form("%s",massfilename.Data()));
gSystem->Exec(Form("pdflatex \"\\def\\mLOne{%1.1f}\\def\\mLTwo{%1.1f}\\def\\mHOne{%1.1f}\\def\\mHTwo{%1.1f}\\include{showplots}\"",massBins[0],massBins[1],massBins[2],massBins[3]));
gSystem->Exec(Form("pdflatex \"\\def\\mLOne{%1.1f}\\def\\mLTwo{%1.1f}\\def\\mHOne{%1.1f}\\def\\mHTwo{%1.1f}\\include{showplots}\"",massBins[0],massBins[1],massBins[2],massBins[3]));
gSystem->Exec("mv showplots.pdf summaryPlots.pdf");
gSystem->Exec("rm -rf showplots.*"); */
}

return;
}


////
TGraphErrors *BuildGraphFb(Double_t massBins[], TString type, Bool_t isQuadr, Bool_t fSigFixed, Bool_t mixedEvents, Int_t color, Int_t style){
TString massfilename="mass_";
TString massfilename1="";
for(int ij=0;ij<nmass; ij++) {
        if(ij != nmass-1) { massfilename.Append(Form("%1.2f_",massBins[ij])); massfilename1.Append(Form("m%d%1.2f_",ij+1,massBins[ij]));  }
        else {massfilename.Append(Form("%1.2f",massBins[ij])); massfilename1.Append(Form("m%d%1.2f",ij+1,massBins[ij]));}
}

Double_t fb, fberr, fsig, fsigerr, chi2x, chi2m;
TGraphErrors *fbVsPt = new TGraphErrors(nptbins+1);
// pt integrated
TString name = Form("%s/fbPt%1.2f_%1.2f_%s_type%s",massfilename.Data(),ptbins[0],ptbins[nptbins],massfilename1.Data(),type.Data());
if(isQuadr) name.Append("_quadratic");
if(fSigFixed) name.Append("_FsigFixed");
if(mixedEvents) name.Append("_ME");
ifstream filenameFb(Form("%s.txt",name.Data()));
filenameFb >> fb >> fberr >> fsig >> fsigerr >> chi2x >> chi2m;
fbVsPt->SetPoint(0,0.75,fb);
fbVsPt->SetPointError(0,0.75,fberr);
    filenameFb.close();
/// 
for(int ipt=0; ipt<nptbins; ipt++){
name = Form("%s/fbPt%1.2f_%1.2f_%s_type%s",massfilename.Data(),ptbins[ipt],ptbins[ipt+1],massfilename1.Data(),type.Data());
if(isQuadr) name.Append("_quadratic");
if(fSigFixed) name.Append("_FsigFixed");
if(mixedEvents) name.Append("_ME");
ifstream filenameFbPt(Form("%s.txt",name.Data()));
filenameFbPt >> fb >> fberr >> fsig >> fsigerr >> chi2x >> chi2m;
fbVsPt->SetPoint(ipt+1,(ptbins[ipt]+ptbins[ipt+1])/2.,fb);
fbVsPt->SetPointError(ipt+1,(ptbins[ipt+1]-ptbins[ipt])/2.,fberr);
    filenameFbPt.close();
}

fbVsPt->SetTitle(" ");
fbVsPt->GetXaxis()->SetTitle("p_{T}(GeV/c)");
fbVsPt->GetYaxis()->SetTitle("f_{B}");
fbVsPt->SetLineColor(color);
fbVsPt->SetMarkerColor(color);
fbVsPt->SetMarkerStyle(style);

return fbVsPt;
}

TGraphErrors *BuildGraphFSig(Double_t massBins[], TString type, Bool_t isQuadr, Bool_t fSigFixed, Int_t color, Int_t style){
TString massfilename="mass_";
TString massfilename1="";
for(int ij=0;ij<nmass; ij++) {
        if(ij != nmass-1) { massfilename.Append(Form("%1.2f_",massBins[ij])); massfilename1.Append(Form("m%d%1.2f_",ij+1,massBins[ij]));  }
        else {massfilename.Append(Form("%1.2f",massBins[ij])); massfilename1.Append(Form("m%d%1.2f",ij+1,massBins[ij]));}
}

Double_t fb, fberr, fsig, fsigerr, chi2x, chi2m;
TGraphErrors *fbVsPt = new TGraphErrors(nptbins+1);
// pt integrated
TString name = Form("%s/fbPt%1.2f_%1.2f_%s_type%s",massfilename.Data(),ptbins[0],ptbins[nptbins],massfilename1.Data(),type.Data());
if(isQuadr) name.Append("_quadratic");
if(fSigFixed) name.Append("_FsigFixed");
ifstream filenameFb(Form("%s.txt",name.Data()));
filenameFb >> fb >> fberr >> fsig >> fsigerr >> chi2x >> chi2m;
fbVsPt->SetPoint(0,0.75,fsig);
fbVsPt->SetPointError(0,0.75,fsigerr);
    filenameFb.close();
/// 
for(int ipt=0; ipt<nptbins; ipt++){
name = Form("%s/fbPt%1.2f_%1.2f_%s_type%s",massfilename.Data(),ptbins[ipt],ptbins[ipt+1],massfilename1.Data(),type.Data());
if(isQuadr) name.Append("_quadratic");
if(fSigFixed) name.Append("_FsigFixed");
ifstream filenameFbPt(Form("%s.txt",name.Data()));
filenameFbPt >> fb >> fberr >> fsig >> fsigerr >> chi2x >> chi2m;
fbVsPt->SetPoint(ipt+1,(ptbins[ipt]+ptbins[ipt+1])/2.,fsig);
fbVsPt->SetPointError(ipt+1,(ptbins[ipt+1]-ptbins[ipt])/2.,fsigerr);
    filenameFbPt.close();
}
fbVsPt->SetTitle(" ");
fbVsPt->GetXaxis()->SetTitle("p_{T}(GeV/c)");
fbVsPt->GetYaxis()->SetTitle("f_{Sig}");
fbVsPt->SetLineColor(color);
fbVsPt->SetMarkerColor(color);
fbVsPt->SetMarkerStyle(style);

return fbVsPt;
}


TGraphErrors *BuildGraphChi2X(Double_t massBins[], TString type, Bool_t isQuadr, Bool_t fSigFixed, Int_t color, Int_t style){
TString massfilename="mass_";
TString massfilename1="";
for(int ij=0;ij<nmass; ij++) {
        if(ij != nmass-1) { massfilename.Append(Form("%1.2f_",massBins[ij])); massfilename1.Append(Form("m%d%1.2f_",ij+1,massBins[ij]));  }
        else {massfilename.Append(Form("%1.2f",massBins[ij])); massfilename1.Append(Form("m%d%1.2f",ij+1,massBins[ij]));}
}

Double_t fb, fberr, fsig, fsigerr, chi2x, chi2m;
TGraphErrors *fbVsPt = new TGraphErrors(nptbins+1);
// pt integrated
TString name = Form("%s/fbPt%1.2f_%1.2f_%s_type%s",massfilename.Data(),ptbins[0],ptbins[nptbins],massfilename1.Data(),type.Data());
if(isQuadr) name.Append("_quadratic");
if(fSigFixed) name.Append("_FsigFixed");
ifstream filenameFb(Form("%s.txt",name.Data()));
filenameFb >> fb >> fberr >> fsig >> fsigerr >> chi2x >> chi2m;
fbVsPt->SetPoint(0,0.75,chi2x);
fbVsPt->SetPointError(0,0.75,0.);
    filenameFb.close();
/// 
for(int ipt=0; ipt<nptbins; ipt++){
name = Form("%s/fbPt%1.2f_%1.2f_%s_type%s",massfilename.Data(),ptbins[ipt],ptbins[ipt+1],massfilename1.Data(),type.Data());
if(isQuadr) name.Append("_quadratic");
if(fSigFixed) name.Append("_FsigFixed");
ifstream filenameFbPt(Form("%s.txt",name.Data()));
filenameFbPt >> fb >> fberr >> fsig >> fsigerr >> chi2x >> chi2m;
fbVsPt->SetPoint(ipt+1,(ptbins[ipt]+ptbins[ipt+1])/2.,chi2x);
fbVsPt->SetPointError(ipt+1,(ptbins[ipt+1]-ptbins[ipt])/2.,0.);
    filenameFbPt.close();
}

fbVsPt->SetTitle(" ");
fbVsPt->GetXaxis()->SetTitle("p_{T}(GeV/c)");
fbVsPt->GetYaxis()->SetTitle("#Chi^{2}/NDF");
fbVsPt->SetLineColor(color);
fbVsPt->SetMarkerColor(color);
fbVsPt->SetMarkerStyle(style);


return fbVsPt;
}

TGraphErrors *BuildGraphChi2M(Double_t massBins[], TString type, Bool_t isQuadr, Bool_t fSigFixed, Int_t color, Int_t style){
TString massfilename="mass_";
TString massfilename1="";
for(int ij=0;ij<nmass; ij++) {
        if(ij != nmass-1) { massfilename.Append(Form("%1.2f_",massBins[ij])); massfilename1.Append(Form("m%d%1.2f_",ij+1,massBins[ij]));  }
        else {massfilename.Append(Form("%1.2f",massBins[ij])); massfilename1.Append(Form("m%d%1.2f",ij+1,massBins[ij]));}
}

Double_t fb, fberr, fsig, fsigerr, chi2x, chi2m;
TGraphErrors *fbVsPt = new TGraphErrors(nptbins+1);
// pt integrated
TString name = Form("%s/fbPt%1.2f_%1.2f_%s_type%s",massfilename.Data(),ptbins[0],ptbins[nptbins],massfilename1.Data(),type.Data());
if(isQuadr) name.Append("_quadratic");
if(fSigFixed) name.Append("_FsigFixed");
ifstream filenameFb(Form("%s.txt",name.Data()));
filenameFb >> fb >> fberr >> fsig >> fsigerr >> chi2x >> chi2m;
fbVsPt->SetPoint(0,0.75,chi2m);
fbVsPt->SetPointError(0,0.75,0.);
    filenameFb.close();
///
for(int ipt=0; ipt<nptbins; ipt++){
name = Form("%s/fbPt%1.2f_%1.2f_%s_type%s",massfilename.Data(),ptbins[ipt],ptbins[ipt+1],massfilename1.Data(),type.Data());
if(isQuadr) name.Append("_quadratic");
if(fSigFixed) name.Append("_FsigFixed");
ifstream filenameFbPt(Form("%s.txt",name.Data()));
filenameFbPt >> fb >> fberr >> fsig >> fsigerr >> chi2x >> chi2m;
fbVsPt->SetPoint(ipt+1,(ptbins[ipt]+ptbins[ipt+1])/2.,chi2m);
fbVsPt->SetPointError(ipt+1,(ptbins[ipt+1]-ptbins[ipt])/2.,0.);
    filenameFbPt.close();
}
fbVsPt->SetTitle(" ");
fbVsPt->GetXaxis()->SetTitle("p_{T}(GeV/c)");
fbVsPt->GetYaxis()->SetTitle("#Chi^{2}/NDF");
fbVsPt->SetLineColor(color);
fbVsPt->SetMarkerColor(color);
fbVsPt->SetMarkerStyle(style);


return fbVsPt;
}


TGraphErrors *BuildGraphFSig1D(Double_t massBins[], TString type, Int_t color, Int_t style){
TString massfilename="mass_";
TString massfilename1="";
for(int ij=0;ij<nmass; ij++) {
        if(ij != nmass-1) { massfilename.Append(Form("%1.2f_",massBins[ij])); massfilename1.Append(Form("m%d%1.2f_",ij+1,massBins[ij]));  }
        else {massfilename.Append(Form("%1.2f",massBins[ij])); massfilename1.Append(Form("m%d%1.2f",ij+1,massBins[ij]));}
}

Double_t fsig1D, fsigerr1D;
TGraphErrors *fbVsPt = new TGraphErrors(nptbins+1);
// pt integrated
TString name = Form("%s/fSIGPt%1.2f_%1.2f_%1.2f_%1.2f_type%s.txt",massfilename.Data(),ptbins[0],ptbins[nptbins],massBins[0],massBins[nmass-1],type.Data());
ifstream filenameFsig(name);
filenameFsig >> fsig1D >> fsigerr1D;
fbVsPt->SetPoint(0,0.75,fsig1D);
fbVsPt->SetPointError(0,0.75,0.);
    filenameFsig.close();
/// 
for(int ipt=0; ipt<nptbins; ipt++){
name = Form("%s/fSIGPt%1.2f_%1.2f_%1.2f_%1.2f_type%s.txt",massfilename.Data(),ptbins[ipt],ptbins[ipt+1],massBins[0],massBins[nmass-1],type.Data());
ifstream filenameFsig(Form("%s",name.Data()));
filenameFsig >> fsig1D >> fsigerr1D;
fbVsPt->SetPoint(ipt+1,(ptbins[ipt]+ptbins[ipt+1])/2.,fsig1D);
fbVsPt->SetPointError(ipt+1,(ptbins[ipt+1]-ptbins[ipt])/2.,0.);
    filenameFsig.close();
}
fbVsPt->SetTitle(" ");
fbVsPt->GetXaxis()->SetTitle("p_{T}(GeV/c)");
fbVsPt->GetYaxis()->SetTitle("f_{Sig}");
fbVsPt->SetLineColor(color);
fbVsPt->SetMarkerColor(color);
fbVsPt->SetMarkerStyle(style);

return fbVsPt;
}

TGraphErrors *BuildGraphFbCompTypes(Double_t massBins[], TString types[], Bool_t isQuadr, Bool_t fSigFixed, Int_t bkgType,TString nameGraph, Int_t color, Int_t style){

TString massfilename="mass_";
TString massfilename1="";
for(int ij=0;ij<nmass; ij++) {
        if(ij != nmass-1) { massfilename.Append(Form("%1.2f_",massBins[ij])); massfilename1.Append(Form("m%d%1.2f_",ij+1,massBins[ij]));  }
        else {massfilename.Append(Form("%1.2f",massBins[ij])); massfilename1.Append(Form("m%d%1.2f",ij+1,massBins[ij]));}
}

Double_t fb, fberr, fsig, fsigerr, chi2x, chi2m;
TGraphErrors *fbVsPt = new TGraphErrors(nptbins+1);
// pt integrated
TString name = Form("%s/fbPt%1.2f_%1.2f_%s_type%s",massfilename.Data(),ptbins[0],ptbins[nptbins],massfilename1.Data(),types[0].Data());
if(isQuadr) name.Append("_quadratic");
if(fSigFixed) name.Append("_FsigFixed");
if(bkgType==1) name.Append("_ME");
else if(bkgType==2) name.Append("_LS");
ifstream filenameFb(Form("%s.txt",name.Data()));
filenameFb >> fb >> fberr >> fsig >> fsigerr >> chi2x >> chi2m;
fbVsPt->SetPoint(0,0.75,fb);
fbVsPt->SetPointError(0,0.75,fberr);
    filenameFb.close();
///
for(int ipt=0; ipt<nptbins; ipt++){
name = Form("%s/fbPt%1.2f_%1.2f_%s_type%s",massfilename.Data(),ptbins[ipt],ptbins[ipt+1],massfilename1.Data(),types[ipt+1].Data());
if(isQuadr) name.Append("_quadratic");
if(fSigFixed) name.Append("_FsigFixed");
    if(bkgType == 1) name.Append("_ME");
    else if(bkgType == 2) name.Append("_LS");
ifstream filenameFbPt(Form("%s.txt",name.Data()));
filenameFbPt >> fb >> fberr >> fsig >> fsigerr >> chi2x >> chi2m;
fbVsPt->SetPoint(ipt+1,(ptbins[ipt]+ptbins[ipt+1])/2.,fb);
fbVsPt->SetPointError(ipt+1,(ptbins[ipt+1]-ptbins[ipt])/2.,fberr);
    filenameFbPt.close();
}

fbVsPt->SetTitle(" ");
fbVsPt->SetName(nameGraph);
fbVsPt->GetXaxis()->SetTitle("p_{T}(GeV/c)");
fbVsPt->GetYaxis()->SetTitle("f_{B}");
fbVsPt->SetLineColor(color);
fbVsPt->SetMarkerColor(color);
fbVsPt->SetMarkerStyle(style);

return fbVsPt;

}

TGraphErrors *BuildGraphChi2XProjCompTypes(Double_t massBins[], TString types[], Bool_t isQuadr, Bool_t fSigFixed, Bool_t mixedEvents,TString nameGraph, Int_t color, Int_t style){

TString massfilename="mass_";
TString massfilename1="";
for(int ij=0;ij<nmass; ij++) {
        if(ij != nmass-1) { massfilename.Append(Form("%1.2f_",massBins[ij])); massfilename1.Append(Form("m%d%1.2f_",ij+1,massBins[ij]));  }
        else {massfilename.Append(Form("%1.2f",massBins[ij])); massfilename1.Append(Form("m%d%1.2f",ij+1,massBins[ij]));}
}

Double_t fb, fberr, fsig, fsigerr, chi2x, chi2m;
TGraphErrors *fbVsPt = new TGraphErrors(nptbins+1);
// pt integrated
TString name = Form("%s/fbPt%1.2f_%1.2f_%s_type%s",massfilename.Data(),ptbins[0],ptbins[nptbins],massfilename1.Data(),types[0].Data());
if(isQuadr) name.Append("_quadratic");
if(fSigFixed) name.Append("_FsigFixed");
if(mixedEvents) name.Append("_ME");
ifstream filenameFb(Form("%s.txt",name.Data()));
filenameFb >> fb >> fberr >> fsig >> fsigerr >> chi2x >> chi2m;
fbVsPt->SetPoint(0,0.75,chi2x);
fbVsPt->SetPointError(0,0.75,0.);
    filenameFb.close();
///
for(int ipt=0; ipt<nptbins; ipt++){
name = Form("%s/fbPt%1.2f_%1.2f_%s_type%s",massfilename.Data(),ptbins[ipt],ptbins[ipt+1],massfilename1.Data(),types[ipt+1].Data());
if(isQuadr) name.Append("_quadratic");
if(fSigFixed) name.Append("_FsigFixed");
    if(mixedEvents) name.Append("_ME");
ifstream filenameFbPt(Form("%s.txt",name.Data()));
filenameFbPt >> fb >> fberr >> fsig >> fsigerr >> chi2x >> chi2m;
fbVsPt->SetPoint(ipt+1,(ptbins[ipt]+ptbins[ipt+1])/2.,chi2x);
fbVsPt->SetPointError(ipt+1,(ptbins[ipt+1]-ptbins[ipt])/2.,0.);
    filenameFbPt.close();
}

fbVsPt->SetTitle(" ");
fbVsPt->SetName(nameGraph);
fbVsPt->GetXaxis()->SetTitle("p_{T}(GeV/c)");
fbVsPt->GetYaxis()->SetTitle("Chi^{2} (x-projection)");
fbVsPt->SetLineColor(color);
fbVsPt->SetMarkerColor(color);
fbVsPt->SetMarkerStyle(style);

return fbVsPt;

}

TGraphErrors *BuildGraphChi2MProjCompTypes(Double_t massBins[], TString types[], Bool_t isQuadr, Bool_t fSigFixed, Bool_t mixedEvents,TString nameGraph, Int_t color, Int_t style){

TString massfilename="mass_";
TString massfilename1="";
for(int ij=0;ij<nmass; ij++) {
        if(ij != nmass-1) { massfilename.Append(Form("%1.2f_",massBins[ij])); massfilename1.Append(Form("m%d%1.2f_",ij+1,massBins[ij]));  }
        else {massfilename.Append(Form("%1.2f",massBins[ij])); massfilename1.Append(Form("m%d%1.2f",ij+1,massBins[ij]));}
}

Double_t fb, fberr, fsig, fsigerr, chi2x, chi2m;
TGraphErrors *fbVsPt = new TGraphErrors(nptbins+1);
// pt integrated
TString name = Form("%s/fbPt%1.2f_%1.2f_%s_type%s",massfilename.Data(),ptbins[0],ptbins[nptbins],massfilename1.Data(),types[0].Data());
if(isQuadr) name.Append("_quadratic");
if(fSigFixed) name.Append("_FsigFixed");
if(mixedEvents) name.Append("_ME");
ifstream filenameFb(Form("%s.txt",name.Data()));
filenameFb >> fb >> fberr >> fsig >> fsigerr >> chi2x >> chi2m;
fbVsPt->SetPoint(0,0.75,chi2m);
fbVsPt->SetPointError(0,0.75,0.);
    filenameFb.close();
///
for(int ipt=0; ipt<nptbins; ipt++){
name = Form("%s/fbPt%1.2f_%1.2f_%s_type%s",massfilename.Data(),ptbins[ipt],ptbins[ipt+1],massfilename1.Data(),types[ipt+1].Data());
if(isQuadr) name.Append("_quadratic");
if(fSigFixed) name.Append("_FsigFixed");
    if(mixedEvents) name.Append("_ME");
ifstream filenameFbPt(Form("%s.txt",name.Data()));
filenameFbPt >> fb >> fberr >> fsig >> fsigerr >> chi2x >> chi2m;
fbVsPt->SetPoint(ipt+1,(ptbins[ipt]+ptbins[ipt+1])/2.,chi2m);
fbVsPt->SetPointError(ipt+1,(ptbins[ipt+1]-ptbins[ipt])/2.,0.);
    filenameFbPt.close();
}

fbVsPt->SetTitle(" ");
fbVsPt->SetName(nameGraph);
fbVsPt->GetXaxis()->SetTitle("p_{T}(GeV/c)");
fbVsPt->GetYaxis()->SetTitle("Chi^{2} (m-projection)");
fbVsPt->SetLineColor(color);
fbVsPt->SetMarkerColor(color);
fbVsPt->SetMarkerStyle(style);

return fbVsPt;

}

