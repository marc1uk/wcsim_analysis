{
TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_1k_photon_bomb_0.root");
TTree* t = (TTree*)f->Get("wcsimT");
WCSimRootEvent* e=0;
TBranch* b=0;
t->SetBranchAddress("wcsimrootevent", &e, &b);
b->GetEntry(0);
WCSimRootTrigger* r=e->GetTrigger(0);
WCSimRootEventHeader* h=r->GetHeader();
h->GetDate();

for(int i=0; i<min(1,(int)b->GetEntries()); i++){
  for(int j=0; j<min(1,(int)e->GetNumberOfEvents()); j++){
    b->GetEntry(i);
    r=e->GetTrigger(j);
    h=r->GetHeader();
    int numcherenkovhitttimes=r->GetNcherenkovhittimes();
    cout<<"event "<<i<<", trigger "<<j<<" has "<<numcherenkovhitttimes<<" ChereknovHitTime entries"<<endl;
    int ndigits=r->GetNcherenkovdigihits();
    for(int digiti=0; digiti<min(10,(int)ndigits); digiti++){
      WCSimRootCherenkovDigiHit* digihit=(WCSimRootCherenkovDigiHit*)(r->GetCherenkovDigiHits()->At(digiti));
      std::vector<int> truephotonindices = digihit->GetPhotonIds();
      cout<<"digit "<<digiti<<" has "<<truephotonindices.size()<<" true photons"<<endl;
      int numphotonsnotfrome=0;
      for(int truephoton=0; truephoton<min(5,(int)truephotonindices.size()); truephoton++){
        int thephotonsid = truephotonindices.at(truephoton);
        cout<<"     digit "<<digiti<<", photon "<<truephoton<<" has CherenkovHitTimeIndex "<<thephotonsid<<endl;
        WCSimRootCherenkovHitTime *thehittimeobject = (WCSimRootCherenkovHitTime*)r->GetCherenkovHitTimes()->At(thephotonsid);
        Int_t thephotonsparenttrackid = thehittimeobject->GetParentID();
        cout<<"          photon parentid is "<<thephotonsparenttrackid<<endl;
//        if(thephotonsparenttrackid!=1) {
//          cout<<"      digit "<<digiti<<", photon "<<truephoton<<" had parentid "<<thephotonsparenttrackid<<endl;
//          numphotonsnotfrome++;
//        }
//      if(numphotonsnotfrome!=0){  //this digit had contribution from OTHERS!
//        cout<<"################ "<<numphotonsnotfrome<<" PHOTONs NOT FROM PRIMARY e ############"<<endl;
//      }
      }
    }
  }
}

}

