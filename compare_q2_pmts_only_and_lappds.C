{
cout<<"beginning"<<endl;
const char* pwd = gSystem->pwd();
TFile* fileall = TFile::Open(TString::Format("%s/out/q2comparison/30ps_split/pmts_and_lappds/q2comparison.root",pwd));
TTree* t = (TTree*)fileall->Get("q2tree");
//TH1D* recoq2histpmtslappds=(TH1D*)fileall->Get("recoQ2hist");
t->Draw("Q2Err>>recoq2histpmtslappds(100,0,1)");
double totentriesnormalise=recoq2histpmtslappds->GetEntries();
double axismin=recoq2histpmtslappds->GetXaxis()->GetXmin();
double axismax=recoq2histpmtslappds->GetXaxis()->GetXmax();
int numbins = recoq2histpmtslappds->GetNbinsX();
TH1F cumulativerecoq2pmtslappds("cumulativerecoq2pmtslappds","Cumulative Error in Reconstructed Q^{2};Momentum Transfer [(MeV/c)^{2}];Cumulative Fraction of Events",numbins,axismin,axismax);
double runningsum=0;
cout<<"running"<<endl;
for(int i=0; i<numbins+1;i++){ // XXX NOTE WE HAVE TO GO OVER BY 1 BIN!
  runningsum+=recoq2histpmtslappds->GetBinContent(i);
  //cumulativerecoq2pmtslappds.SetBinContent(i,(runningsum/totentriesnormalise));
  cumulativerecoq2pmtslappds.SetBinContent(i,recoq2histpmtslappds->GetBinContent(i));
}

cout<<"file 2"<<endl;
TFile* filepmtsonly = TFile::Open(TString::Format("%s/out/q2comparison/30ps_split/pmts_only/q2comparison.root",pwd));
TTree* t2 = (TTree*)filepmtsonly->Get("q2tree");
//TH1D* recoq2histpmtsonly=(TH1D*)filepmtsonly->Get("recoQ2hist");
t2->Draw("Q2Err>>recoq2histpmtsonly(100,0,1)");
double totentriesnormalise2=recoq2histpmtsonly->GetEntries();
//double axismin=recoq2histpmtsonly->GetXaxis()->GetXmin();
//double axismax=recoq2histpmtsonly->GetXaxis()->GetXmax();
//int numbins = recoq2histpmtsonly->GetNbinsX();
TH1F cumulativerecoq2pmtsonly("cumulativerecoq2pmtsonly","Cumulative Error in Reconstructed Q^{2};Momentum Transfer [(MeV/c)^{2}];Cumulative Fraction of Events",numbins,axismin,axismax);
//double runningsum=0;
runningsum=0;
cout<<"running"<<endl;
for(int i=0; i<numbins+1;i++){
  runningsum+=recoq2histpmtsonly->GetBinContent(i);
  //cumulativerecoq2pmtsonly.SetBinContent(i,(runningsum/totentriesnormalise2));
  cumulativerecoq2pmtsonly.SetBinContent(i,recoq2histpmtsonly->GetBinContent(i));
}

TCanvas mycanv;
cumulativerecoq2pmtslappds.SetTitle("PMTs + LAPPDs");
cumulativerecoq2pmtslappds.SetLineColor(kRed);
cumulativerecoq2pmtsonly.SetLineColor(kBlue);
cumulativerecoq2pmtsonly.SetTitle("PMTs only");
cumulativerecoq2pmtslappds.Draw(); //DrawNormalized("",1);
cumulativerecoq2pmtsonly.Draw("same"); //DrawNormalized("same",1);
TLegend* aleg = mycanv.BuildLegend();

}
