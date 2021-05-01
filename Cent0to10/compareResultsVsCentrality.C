void compareResultsVsCentrality(){

TFile *fc10to30 = new TFile("resultsVsCentrality/resultsFb10to30.root","READ");
TGraphErrors *fbStat10to30 = (TGraphErrors*)fc10to30->Get("fBstat");

TFile *fc30to50 = new TFile("resultsVsCentrality/resultsFb30to50.root","READ");
TGraphErrors *fbStat30to50 = (TGraphErrors*)fc30to50->Get("fBstat");

fbStat10to30->SetMarkerStyle(20);
fbStat10to30->SetMarkerColor(2);
fbStat10to30->SetLineColor(2);

fbStat30to50->SetMarkerStyle(24);
fbStat30to50->SetMarkerColor(4);
fbStat30to50->SetLineColor(4);


fbStat10to30->SetTitle(" ");
fbStat10to30->Draw("AP");
fbStat30to50->Draw("P SAME");

	TLegend *legend=new TLegend(0.17,0.70,0.29,0.88);
        legend->SetBorderSize(0); legend->SetFillColor(0); legend->SetTextFont(42);
        legend->SetFillStyle(0); legend->SetMargin(0.25);
        legend->SetEntrySeparation(0.15);
        legend->SetTextSize(0.045);
	legend->AddEntry(fbStat10to30,"10-30%","p");
	legend->AddEntry(fbStat30to50,"30-50%","p");
        legend->Draw();

return;
}
