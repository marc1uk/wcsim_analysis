{
TFile *_file0 = TFile::Open("/annie/app/users/moflaher/wcsim/root_work/out/q2comparison/30ps_proper_split/trueQEvertexinfo_ext2.root");
t = (TTree*)_file0->Get("vertextreenocuts");
std::vector<string> digitdets;
vector<string>* digitdetsp=&digitdets;
t->SetBranchAddress("DigitWhichDet",&digitdetsp);
std::vector<int> npmtdigits;
std::vector<int> nlappddigits;
for(int entryi=0; entryi<t->GetEntries(); entryi++){
  if((entryi%1000)==0) cout<<"entryi="<<entryi<<endl;
  t->GetEntry(entryi);
  int npmtdigitsthisevt=0;
  int nlappddigitsthisevt=0;
  for(int digiti=0; digiti<digitdets.size(); digiti++){
    string digitidet = digitdets.at(digiti);
    if(digitidet=="lappd_v0") nlappddigitsthisevt++;
    else npmtdigitsthisevt++;
  }
  npmtdigits.push_back(npmtdigitsthisevt);
  nlappddigits.push_back(nlappddigitsthisevt);
}
int maxhitsinevt= *max_element(nlappddigitshist.begin(), nlappddigitshist.end());
// actually this is far too big - LAPPD hits have tail up to 20k!! PMTs max ~200 (~nPMTS)! 
// to plot on same axis,up to ~2000 is sufficient.
TH1F* npmtdigitshist = new TH1F("npmtdigitshist","Num PMT Digits in an Event",200,0,maxhitsinevt/10);
for(auto acount : npmtdigits) npmtdigitshist->Fill(acount);
TH1F* nlappddigitshist = new TH1F("nlappddigitshist","Num LAPPD Digits in an Event",200,0,maxhitsinevt/10);
for(auto acount : nlappddigits) nlappddigitshist->Fill(acount);
npmtdigitshist->SetLineColor(kRed);
nlappddigitshist->SetLineColor(kBlue);
TCanvas* c1 = new TCanvas("c1", "c1",13,50,1174,727);
npmtdigitshist->Draw();
nlappddigitshist->Draw("same");
c1->BuildLegend();
}
