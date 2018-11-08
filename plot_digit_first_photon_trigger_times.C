{
//TFile* f= TFile::Open("/pnfs/annie/persistent/users/moflaher/wcsim/wcsim_0.2000.root");
//TFile* f=TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/out/ANNIEtest_10MeV_e-_Uni_Iso_annie_f717ff7d765d00225801751243aad1027d7924df_0.root");
//TFile* f=TFile::Open("/home/marc/LinuxSystemFiles/Bonsai/validation/Mahdi_Bonsaifiles/e5d1000PMT3.root");
TFile* f=TFile::Open("/annie/app/users/moflaher/wcsim/ANNIEp2v5_beamouts/wcsim_ANNIEp2v5_beam_0.root");
TTree* t = (TTree*)f->Get("wcsimT");
WCSimRootEvent* e=0;
TBranch* b=0;
t->SetBranchAddress("wcsimrootevent", &e, &b);
b->GetEntry(0);
WCSimRootTrigger* r=e->GetTrigger(0);
WCSimRootEventHeader* h=r->GetHeader();
h->GetDate();

long long int MAX_EVENTS_TO_PLOT=10000;
//int TRIGGER_OFFSET=0;
int TRIGGER_OFFSET=950; // historic offset if not disabled

TH1D hfirsttriggertime=TH1D("hfirsttriggertime","Time of event trigger",100,6660,6720);
TH1D hdigittimes=TH1D("hdigittimes","Times of digits within the first trigger",200,TRIGGER_OFFSET-50,TRIGGER_OFFSET+150);
TH1D hfirstphotonindigittimes=TH1D("hfirstphotonindigittimes","Times of the first photon within each Digit",60,0,60);

for(int i=0; i<std::min(t->GetEntries(),MAX_EVENTS_TO_PLOT); i++){
  for(int j=0; j<e->GetNumberOfEvents(); j++){
    b->GetEntry(i);
    r=e->GetTrigger(j);
    h=r->GetHeader();
    hfirsttriggertime.Fill(h->GetDate());  // note trigger time
    int ndigits=r->GetNcherenkovdigihits();
    cout<<"event "<<i<<" has "<<ndigits<<" digits ";
    float firsttime;
    double firsthittime;
    if(ndigits){
      WCSimRootCherenkovDigiHit* hit=(WCSimRootCherenkovDigiHit*)(r->GetCherenkovDigiHits()->At(0));
      firsttime=hit->GetT();
      if(j==0) hdigittimes.Fill(firsttime); // note first digit in first trigger
      std::vector<int> truephotonindices = hit->GetPhotonIds();
      int firstphotonsid = truephotonindices.at(0);
      WCSimRootCherenkovHitTime* thehittimeobject= (WCSimRootCherenkovHitTime*)r->GetCherenkovHitTimes()->At (firstphotonsid);
      firsthittime = thehittimeobject->GetTruetime();
      if(j==0) hfirstphotonindigittimes.Fill(firsthittime); // note first photon of first digit in first trigger
      cout<<"with first digit time "<<firsttime<<" and first photon time "<<firsthittime<<" ";
      cout<<"and date "<<h->GetDate()<<" correspnding to a first digit time of "<<(firsttime-950+(h->GetDate()))<<endl;
    }
  }
}

TCanvas c1("c1","canv1",600,600);
hfirsttriggertime.Draw();
TCanvas c2("c2","canv2",600,600);
hdigittimes.Draw();
TCanvas c3("c3","canv3",600,600);
hfirstphotonindigittimes.Draw();
}
