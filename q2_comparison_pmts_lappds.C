{
TFile* _file0 = TFile::Open("out/q2comparison/30ps_proper_split/pmts_only/q2comparison.root");
TFile* _file1 = TFile::Open("out/q2comparison/30ps_proper_split/pmts_and_lappds/q2comparison.root");
TTree* tp = (TTree*)_file0->Get("q2tree");
TTree* tpl = (TTree*)_file1->Get("q2tree");
//tp->Draw("(abs(Q2Err)/TrueQ2)>>pq2e(50,0,3)");
//tpl->Draw("(abs(Q2Err)/TrueQ2)>>plq2e(50,0,3)");
tp->Draw("(abs(Q2Err))>>pq2e(50,0,1.2)");
tpl->Draw("(abs(Q2Err))>>plq2e(50,0,1.2)");
pq2e->SetLineColor(kBlue);
pq2e->SetTitle("128 PMTs");
plq2e->SetLineColor(kRed);
plq2e->SetTitle("5 LAPPDs + 128 PMTs");
plq2e->GetXaxis()->SetTitle("|#Delta Q^{2}| [(GeV/c)^{2}]");
plq2e->GetYaxis()->SetTitle("Percentage");
CenterTitles(plq2e);
plq2e->DrawNormalized("",100);
pq2e->DrawNormalized("same",100);
c1->BuildLegend();
Simulation();

TH1D* pcum = new TH1D("pcum;|#Delta Q^{2}|;Percentage","128 PMTs",50,0,2); // max 3 for relative error
TH1D* plcum = new TH1D("plcum","5 LAPPDs + 128 PMTs",50,0,3);
double totalentriesp=pq2e->GetEntries();
double totalentriespl=plq2e->GetEntries();
double runningsump=0, runningsumpl=0;
for(int bini=0; bini<52; bini++){
  runningsump+=pq2e->GetBinContent(bini);
  runningsumpl+=plq2e->GetBinContent(bini);
  pcum->SetBinContent(bini,(runningsump/totalentriesp)*100.);
  plcum->SetBinContent(bini,(runningsumpl/totalentriespl)*100.);
}
pcum->SetLineColor(kBlue);
pcum->GetXaxis()->SetTitle("|#Delta Q^{2}| [(GeV/c)^{2}]");
pcum->GetYaxis()->SetTitle("Percentage");
pcum->GetYaxis()->SetRangeUser(0,105);
CenterTitles(pcum);
plcum->SetLineColor(kRed);
pcum->Draw();
plcum->Draw("same");
c1->BuildLegend();
Simulation();

}



