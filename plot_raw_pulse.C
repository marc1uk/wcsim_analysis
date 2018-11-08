//{
//TFile *_file0 = TFile::Open("/home/marc/anniegpvm/rawcopy/RAWDataR897S0p0.root");
//TTree* t  = (TTree*)_file0->Get("PMTData");
//uint16_t* z = new uint16_t[160000];
//TGraph* g = new TGraph(160000);
//t->SetBranchAddress("Data",z);
//t->GetEntry(12);
//for(int i=0; i<160000; i++){g->SetPoint(i,(float)i,(float)z[i]);}
//g->Draw();
//}
{
//TFile *_file0 = TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/RAWDataR4000S0p0.root");
TFile *_file0 = TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/RAWDataR0S0p0.root");
TTree* t  = (TTree*)_file0->Get("PMTData");
uint16_t* z = new uint16_t[108000];
TGraph* g = new TGraph(27000);
t->SetBranchAddress("Data",z);
t->GetEntry(12);
for(int i=0; i<27000; i++){g->SetPoint(i,(float)i,(float)z[i]);}
g->Draw();
}

// simulatons: total trigger time is 1350ns (950+400) / 2ns -> 675 samples/ minibuffer / channel
//-> 675*40 = 27,000 samples / channel = 108,000 points per full buffer.

std::vector<UShort_t> dats(108000,0);
PMTData->SetBranchAddress("Data",dats.data());
int maxentries=100;
for(int ei=0; ei<std::min(int(maxentries),int(PMTData->GetEntries())); ei++){
  PMTData->GetEntry(ei);
  TGraph* g = new TGraph();
  int si=0;
  for(auto aval : dats){
    g->SetPoint(si,si,aval);
    si++;
  }
  g->Draw("alp");
  gPad->WaitPrimitive();
  delete g;
}
