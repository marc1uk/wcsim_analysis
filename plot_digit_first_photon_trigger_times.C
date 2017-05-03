TFile* f= TFile::Open("/pnfs/annie/persistent/users/moflaher/wcsim/wcsim_0.2000.root")
TTree* t = (TTree*)f->Get("wcsimT")
WCSimRootEvent* e=0;
TBranch* b=0;
t->SetBranchAddress("wcsimrootevent", &e, &b);
b->GetEntry(0)
WCSimRootTrigger* r=e->GetTrigger(0)
WCSimRootEventHeader* h=r->GetHeader()
h->GetDate()

for(int i=0; i<100; i++){
  for(int j=0; j<e->GetNumberOfEvents(); j++){
    b->GetEntry(i);
    r=e->GetTrigger(j);
    h=r->GetHeader();
    int ndigits=r->GetNcherenkovdigihits();
    cout<<"event "<<i<<" has "<<ndigits<<" digits ";
    float firsttime;
    double firsthittime;
    if(ndigits){
      WCSimRootCherenkovDigiHit* hit=(WCSimRootCherenkovDigiHit*)(r->GetCherenkovDigiHits()->At(0));
      firsttime=hit->GetT();
      std::vector<int> truephotonindices = hit->GetPhotonIds();
      int firstphotonsid = truephotonindices.at(0);
      WCSimRootCherenkovHitTime* thehittimeobject= (WCSimRootCherenkovHitTime*)r->GetCherenkovHitTimes()->At (firstphotonsid);
      firsthittime = thehittimeobject->GetTruetime();
      cout<<"with first digit time "<<firsttime<<" and first photon time "<<firsthittime<<" ";
      cout<<"and date "<<h->GetDate()<<" correspnding to a first digit time of "<<(firsttime-950+(h->GetDate()))<<endl;
    }
  }
}


