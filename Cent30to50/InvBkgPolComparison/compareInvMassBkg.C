Int_t nptbins = 3;
Double_t ptbins[] = {1.5,3.0,5.0,10.0};
TString types[] = {"FF","FF","FF_FS","FF_FS"};
Int_t nmass = 4;
Double_t mbins[9][4] = {
        {2.50, 2.80, 3.20, 3.60},
        {2.50, 2.75, 3.16, 3.60},
        {2.50, 2.85, 3.25, 3.60},
        {2.50, 2.75, 3.25, 3.60},
        {2.50, 2.85, 3.16, 3.60}
        };
TGraphErrors *BuildGraphFbCompTypes(Int_t polOrder,Double_t massBins[], TString types[], Bool_t isQuadr, Bool_t fSigFixed, Bool_t mixedEvents,TString nameGraph, Int_t color, Int_t style);
void readFbResults(Int_t polOrder,Double_t massBins[], Bool_t pdf);
void plotVsMassEdges(TString nameGraph,Bool_t isQuadr);
void plotResults(Int_t polOrder,Int_t num, Bool_t createPdf);
void plotVsAllMassEdges();
void plotVsMassHypothesis(Int_t num);
void plotVsAllMassHypothesis();
void plotAllResults();
void computeRMSandAverages();


void computeRMSandAveragesForMassEdges(){
TString nameGraphs[] = {"fbVsPtComp","fbVsPtCompME"};
TString nameGraphs2[] = {"fbVsPtCombinedQuadrWeights","fbVsPtCombinedQuadrWeightsME"};
//TString nameGraphs[] = {"fbVsPtComp","fbVsPtCompME","fbVsPtCompFsigFixed","fbVsPtCompFsigFixedME"};
//TString nameGraphs[] = {"fbVsPtComp","fbVsPtCompME","fbVsPtCombinedQuadrWeights","fbVsPtCombinedQuadrWeightsME"};
TString titlePtString[] = {"1.5-10","1.5-3","3-5","5-10"};
TH1F *fbVal[nptbins+1]; 
    TGraphErrors *grAverage = 0x0;
    TGraphErrors *gRelativeUncertainties = 0x0;
        TString path = Form("mass_%1.2f_%1.2f_%1.2f_%1.2f",mbins[0][0],mbins[0][1],mbins[0][2],mbins[0][3]);
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
TCanvas *cfbVal = new TCanvas();
cfbVal->Divide(2,2);
for(int ipt=0; ipt<nptbins+1; ipt++){
cfbVal->cd(ipt+1);
fbVal[ipt]->GetXaxis()->SetTitle("f_{B}");
    fbVal[ipt]->GetYaxis()->SetTitle("Entries");
fbVal[ipt]->DrawCopy();
grAverage->GetPoint(ipt,xx,yy);
grAverage->SetPoint(ipt,xx, fbVal[ipt]->GetMean());
grAverage->SetPointError(ipt,grAverage->GetErrorX(ipt), maxValue[ipt]-minValue[ipt]);
gRelativeUncertainties->SetPoint(ipt,xx,1);
gRelativeUncertainties->SetPointError(ipt,grAverage->GetErrorX(ipt),(maxValue[ipt]-minValue[ipt])/fbVal[ipt]->GetMean());
}
cfbVal->Draw();
TCanvas *cAvRMS = new TCanvas();
cAvRMS->cd();
gRelativeUncertainties->GetYaxis()->SetTitle("Rel. uncertainty (%)");
gRelativeUncertainties->Draw();
cAvRMS->Draw();

TCanvas *cfbVsPtUnc = new TCanvas();
grAverage->GetYaxis()->SetTitle("Average +- RMS");
grAverage->Draw();
cAvRMS->SaveAs("ComparisonGraphs/RelativUncertaintiesOnInvMassBkg.pdf");

}

void computeRMSandAverages(){
    TString methodName="fbVsPtCompME";
    TString methodName2="fbVsPtCombinedQuadrWeightsME";
    TH1F *fbVal[nptbins+1]; 
    TString titlePtString[] = {"1.5-10","1.5-3","3-5","5-10"};
    

    
    TGraphErrors *grAverage = 0x0;
    TGraphErrors *gRelativeUncertainties = 0x0;
    printf("hei0\n\n");

for(int ipt=0; ipt<nptbins+1; ipt++) fbVal[ipt] = new TH1F(Form("Pt bin:%s",titlePtString[ipt].Data()), "", 100, 0.,1.);
    Double_t xx, yy;
for(int imass=0;imass<5;imass++){

    TString path = Form("mass_%1.2f_%1.2f_%1.2f_%1.2f",mbins[imass][0],mbins[imass][1],mbins[imass][2],mbins[imass][3]);
         TFile *fin = new TFile(Form("%s/fbResultsComb.root",path.Data()),"READ");
         TGraphErrors *gr = (TGraphErrors*)fin->Get(methodName);
         if(imass==0) grAverage = (TGraphErrors*)gr->Clone("grAverage");
         if(imass==0) gRelativeUncertainties = (TGraphErrors*)gr->Clone("gRelativeUncertainties");
         TFile *finquadratic = new TFile(Form("%s/fbResultsCombQuadrWeights.root",path.Data()),"READ");
         TGraphErrors *grquadratic = (TGraphErrors*)finquadratic->Get(methodName2);
for(int ipt=0; ipt<nptbins+1;ipt++){
gr->GetPoint(ipt, xx, yy);
fbVal[ipt]->Fill(yy);
grquadratic->GetPoint(ipt, xx, yy);
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
gRelativeUncertainties->GetYaxis()->SetTitle("Rel. uncertainty (%)");
gRelativeUncertainties->Draw();
cAvRMS->Draw();

TCanvas *cfbVsPtUnc = new TCanvas();
grAverage->Draw();

cAvRMS->SaveAs("ComparisonGraphs/RelativUncertaintiesOnXBkg.pdf");
cfbVal->SaveAs("ComparisonGraphs/FbValuesXBkg.pdf");
cfbVsPtUnc->SaveAs("ComparisonGraphs/AverageFbWithUncertaintyXBkg.pdf");

}

void plotVsAllMassHypothesis(){
 
    for(int i=0;i<5;i++) plotVsMassHypothesis(i);
    return;
        
        
}
void plotAllResults(){
 for(int i=0; i<3; i++) plotResults(i+3,0,kFALSE);
 return;
}

void plotResults(Int_t polOrder,Int_t num, Bool_t createPdf=kFALSE){

Double_t mbins[9][4] = {
        {2.50, 2.80, 3.20, 3.60},
        {2.50, 2.75, 3.16, 3.60},
        {2.50, 2.85, 3.25, 3.60},
        {2.50, 2.75, 3.25, 3.60},
        {2.50, 2.85, 3.16, 3.60}
        };
	readFbResults(polOrder,mbins[num],createPdf);
/////
}

void plotVsAllMassEdges(){
 
    
    TString nameGraphs[] = {"fbVsPtComp","fbVsPtCompME","fbVsPtCompFsigFixed","fbVsPtCompFsigFixedME"};
    TString nameGraphs2[] = {"fbVsPtCombinedQuadrWeights","fbVsPtCombinedQuadrWeightsFsigFixed","fbVsPtCombinedQuadrWeightsFsigFixedME","fbVsPtCombinedQuadrWeightsME"};
    for(int i=0; i<4; i++){
     plotVsMassEdges(nameGraphs[i],kFALSE);   
     plotVsMassEdges(nameGraphs2[i],kTRUE);  
    };
}

void plotVsMassHypothesis(Int_t num){
Double_t mbins[9][4] = {
        {2.50, 2.80, 3.20, 3.60},
        {2.50, 2.75, 3.16, 3.60},
        {2.50, 2.85, 3.25, 3.60},
        {2.50, 2.75, 3.25, 3.60},
        {2.50, 2.85, 3.16, 3.60}
        };
        TCanvas *canv = new TCanvas();
        gSystem->Exec(Form("mkdir -p ComparisonGraphs"));
            TLegend *legend=new TLegend(0.37,0.10,0.62,0.30);
  	legend->SetBorderSize(0); legend->SetFillColor(0); legend->SetTextFont(42);
  	legend->SetFillStyle(0); legend->SetMargin(0.25);
  	legend->SetEntrySeparation(0.15);
  	legend->SetTextSize(0.045);
        TString path = Form("mass_%1.2f_%1.2f_%1.2f_%1.2f",mbins[num][0],mbins[num][1],mbins[num][2],mbins[num][3]);
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

void plotVsMassEdges(TString nameGraph, Bool_t isQuadr){
  
	Double_t mbins[9][4] = {
        {2.50, 2.80, 3.20, 3.60},
        {2.50, 2.75, 3.16, 3.60},
        {2.50, 2.85, 3.25, 3.60},
        {2.50, 2.75, 3.25, 3.60},
        {2.50, 2.85, 3.16, 3.60}
        };
        TString quadrString = "";
        if(isQuadr) quadrString.Append("QuadrWeights");
gSystem->Exec(Form("mkdir -p ComparisonGraphs"));
        TLegend *legend=new TLegend(0.37,0.10,0.62,0.30);
  	legend->SetBorderSize(0); legend->SetFillColor(0); legend->SetTextFont(42);
  	legend->SetFillStyle(0); legend->SetMargin(0.25);
  	legend->SetEntrySeparation(0.15);
  	legend->SetTextSize(0.045);
    legend->SetHeader(Form("%s",nameGraph.Data()));

	Int_t colors[] = {1,2,4,6,7};
	Int_t styles[] = {20,24,21,32,28};
        TCanvas *fbout = new TCanvas();
	for(int im=0; im<5; im++){
        TString path = "mass_";
        for(int i=0; i<nmass-1; i++) path.Append(Form("%1.2f_",mbins[im][i]));
        path.Append(Form("%1.2f",mbins[im][nmass-1]));
        
        TFile *fin = new TFile(Form("%s/fbResultsComb%s.root",path.Data(),quadrString.Data()),"READ");
        fin->ls();
        TGraphErrors *gr = (TGraphErrors*)fin->Get(nameGraph);
	gr->SetLineColor(colors[im]);
	gr->SetMarkerColor(colors[im]);
	gr->SetMarkerStyle(styles[im]);
    gr->GetYaxis()->SetRangeUser(0,0.3);
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

void readFbResults(Int_t polOrder,Double_t massBins[], Bool_t pdf){

TString path = "mass_";
for(int i=0; i<nmass-1; i++) path.Append(Form("%1.2f_",massBins[i]));
path.Append(Form("%1.2f",massBins[nmass-1]));


//
TGraphErrors *fbVsPtComp = BuildGraphFbCompTypes(polOrder,massBins,types,kFALSE,kFALSE,kFALSE,"fbVsPtComp",1,20);
TGraphErrors *fbVsPtCompFsigFixed = BuildGraphFbCompTypes(polOrder,massBins,types,kFALSE,kTRUE,kFALSE,"fbVsPtCompFsigFixed",2,24);
TGraphErrors *fbVsPtCompME = BuildGraphFbCompTypes(polOrder,massBins,types,kFALSE,kFALSE,kTRUE,"fbVsPtCompME",4,21);
TGraphErrors *fbVsPtCompFsigFixedME = BuildGraphFbCompTypes(polOrder,massBins,types,kFALSE,kTRUE,kTRUE,"fbVsPtCompFsigFixedME",6,5);
TGraphErrors *fbVsPtCompQuadr = BuildGraphFbCompTypes(polOrder,massBins,types,kTRUE,kFALSE,kFALSE,"fbVsPtCombinedQuadrWeights",1,20);
TGraphErrors *fbVsPtCompQuadrFsigFixed = BuildGraphFbCompTypes(polOrder,massBins,types,kTRUE,kTRUE,kFALSE,"fbVsPtCombinedQuadrWeightsFsigFixed",2,24);
TGraphErrors *fbVsPtCompQuadrFsigFixedME = BuildGraphFbCompTypes(polOrder,massBins,types,kTRUE,kTRUE,kTRUE,"fbVsPtCombinedQuadrWeightsFsigFixedME",4,21);
TGraphErrors *fbVsPtCompQuadrME = BuildGraphFbCompTypes(polOrder,massBins,types,kTRUE,kFALSE,kTRUE,"fbVsPtCombinedQuadrWeightsME",6,5);



TString massfilename="mass_";
TString massfilename1="";
for(int ij=0;ij<nmass; ij++) { 
	if(ij != nmass-1) { massfilename.Append(Form("%1.2f_",massBins[ij])); massfilename1.Append(Form("m%d%1.2f_",ij+1,massBins[ij]));  }
	else {massfilename.Append(Form("%1.2f",massBins[ij])); massfilename1.Append(Form("m%d%1.2f",ij+1,massBins[ij]));}
}
TString polOrderName = Form("pol%d",polOrder);



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
    
    TFile fout(Form("%s/fbResultsComb.root",polOrderName.Data()),"RECREATE");
fbVsPtComp->Write();
fbVsPtCompME->Write();
fbVsPtCompFsigFixed->Write();
fbVsPtCompFsigFixedME->Write();
    
fout.Close();
cFbcomp->SaveAs(Form("%s/fbResults.png",polOrderName.Data()));

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
    
    TFile fout2(Form("%s/fbResultsCombQuadrWeights.root",polOrderName.Data()),"RECREATE");
    fbVsPtCompQuadr->Write();
fbVsPtCompQuadrME->Write();
fbVsPtCompQuadrFsigFixed->Write();
fbVsPtCompQuadrFsigFixedME->Write();
fout2.Close();
    cFbcompQuadratic->SaveAs(Form("%s/fbResultsQuadratic.png",polOrderName.Data()));

/// create pdf
if(pdf){
gSystem->Exec(Form("cp showplots.tex %s",massfilename.Data()));
gSystem->cd(Form("%s",massfilename.Data()));
gSystem->Exec(Form("pdflatex \"\\def\\mLOne{%1.1f}\\def\\mLTwo{%1.1f}\\def\\mHOne{%1.1f}\\def\\mHTwo{%1.1f}\\include{showplots}\"",massBins[0],massBins[1],massBins[2],massBins[3]));
gSystem->Exec(Form("pdflatex \"\\def\\mLOne{%1.1f}\\def\\mLTwo{%1.1f}\\def\\mHOne{%1.1f}\\def\\mHTwo{%1.1f}\\include{showplots}\"",massBins[0],massBins[1],massBins[2],massBins[3]));
gSystem->Exec("mv showplots.pdf summaryPlots.pdf");
gSystem->Exec("rm -rf showplots.*");
}

return;
}


TGraphErrors *BuildGraphFbCompTypes(Int_t polOrder,Double_t massBins[],TString types[], Bool_t isQuadr, Bool_t fSigFixed, Bool_t mixedEvents,TString nameGraph, Int_t color, Int_t style){

    TString massfilename="mass_";
TString massfilename1="";
for(int ij=0;ij<nmass; ij++) {
        if(ij != nmass-1) { massfilename.Append(Form("%1.2f_",massBins[ij])); massfilename1.Append(Form("m%d%1.2f_",ij+1,massBins[ij]));  }
        else {massfilename.Append(Form("%1.2f",massBins[ij])); massfilename1.Append(Form("m%d%1.2f",ij+1,massBins[ij]));}
}
TString polOrderName = Form("pol%d",polOrder);

Double_t fb, fberr, fsig, fsigerr, chi2x, chi2m;
TGraphErrors *fbVsPt = new TGraphErrors(nptbins+1);
// pt integrated
TString name = Form("%s/fbPt%1.2f_%1.2f_%s_type%s",polOrderName.Data(),ptbins[0],ptbins[nptbins],massfilename1.Data(),types[0].Data());
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
name = Form("%s/fbPt%1.2f_%1.2f_%s_type%s",polOrderName.Data(),ptbins[ipt],ptbins[ipt+1],massfilename1.Data(),types[ipt+1].Data());
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
fbVsPt->SetName(nameGraph);
fbVsPt->GetXaxis()->SetTitle("p_{T}(GeV/c)");
fbVsPt->GetYaxis()->SetTitle("f_{B}");
fbVsPt->SetLineColor(color);
fbVsPt->SetMarkerColor(color);
fbVsPt->SetMarkerStyle(style);

return fbVsPt;

}

