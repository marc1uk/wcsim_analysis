{
//TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in_ratpac_comps/wcsim_1k_photon_bomb_0.root");
//TFile* f= TFile::Open("/pnfs/annie/persistent/users/moflaher/wcsim_lappd_24-09-17_BNB_Water_10k_22-05-17/wcsim_0.0.0.root");
//TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_photon_bomb_next_to_wall_0.root");
//TFile* f=TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/out/ANNIEtest_10MeV_e-_Uni_Iso_annie_f717ff7d765d00225801751243aad1027d7924df_0.root");
//TFile* f=TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in/MRD_muon_sample/ANNIEtest_MRD_muon_sample_1a0f480.root");
//TFile* f=TFile::Open("/annie/app/users/moflaher/wcsim/build/wcsim_photon_bomb_0.root");
//TFile* f=TFile::Open("/pnfs/annie/persistent/users/moflaher/wcsim/lappd/tankonly/wcsim_lappd_tankonly_24-09-17_BNB_Water_10k_22-05-17/wcsim_0.0.0.root");
//TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_0.root");
TFile* f=TFile::Open("/annie/app/users/moflaher/wcsim/ANNIEp2v5_beamouts/wcsim_ANNIEp2v5_beam_0.root");
TTree* t = (TTree*)f->Get("wcsimT");
WCSimRootEvent* e=0;
TBranch* b=0;
t->SetBranchAddress("wcsimrootevent", &e, &b);
b->GetEntry(0);
WCSimRootTrigger* r=e->GetTrigger(0);
WCSimRootEventHeader* h=r->GetHeader();
h->GetDate();

int maxentriestoprint=20;
int maxtriggerstoprint=3;
int maxprimariestoprint=1; // always 1
int maxtrackstoprint=10;
int maxdigitstoprint=0;
int maxphotonsperdigittoprint=5;
int maxphotonstoprint=0;
//neutrino vertices
cout<<"This run has "<<b->GetEntries()<<" entries"<<endl;
for(int i=0; i<min(maxentriestoprint,(int)b->GetEntries()); i++){
  //cout<<"NEW EVENT"<<endl<<"======="<<endl;
  b->GetEntry(i);
  //cout<<"This entry had "<<e->GetNumberOfEvents()<<" triggers"<<endl;
  //cout<<"This entry had "<< ( (e->HasSubEvents()) ? "1 or more" : "no" ) <<" delayed triggers"<<endl;
  //cout<<"This entry had "<<e->GetNumberOfSubEvents()<<" delayed triggers"<<endl;
  for(int j=0; j<min(maxtriggerstoprint,(int)e->GetNumberOfEvents()); j++){
    //cout<<"NEXT TRIGGER"<<endl;
    r=e->GetTrigger(j);
    h=r->GetHeader();
    //cout<<" >>> Trigger time was : "<<h->GetDate()<<endl; // TRIGGER TIME OF 0 MEANS NO NDIGITS TRIGGER
    int nprimaries = r->GetNvtxs();	// always 1
    //cout<<"This trigger had "<<nprimaries<<" vertex(es)"<<endl;
    for(int primaryi=0; primaryi<min(maxprimariestoprint,nprimaries); primaryi++){
      Float_t True_Vertex_X=r->GetVtx(0);
      Float_t True_Vertex_Y=r->GetVtx(1);
      Float_t True_Vertex_Z=r->GetVtx(2);
      int numtracks= r->GetNtrack();
      int numphotons = r->GetCherenkovHits()->GetEntries();
      int ndigits=r->GetCherenkovDigiHits()->GetEntries();
      cout<<"event "<<i<<" trigger "<<j<<" primary "<<primaryi<<" at ("
          <<True_Vertex_X<<", "<<True_Vertex_Y<<", "<<True_Vertex_Z<<") had:"<<endl
          <<"Num Tracks: "<<numtracks<<endl
          <<"Num Photons: " <<numphotons<<endl
          <<"Num Digits: "<<ndigits<<endl;
      // scan through the truth tracks, find the primary muon and pull vertex info from it
      int primarytrack=0;
      int primarycount=0;
      for(int track=0; track<min(maxtrackstoprint,numtracks); track++){
        WCSimRootTrack* nextrack = (WCSimRootTrack*)r->GetTracks()->At(track);
        //Int_t     GetIpnu()             pdg
        //Int_t     GetFlag()             -1: neutrino primary, -2: neutrino target, 0: other
        //Float_t   GetM()                mass
        //Float_t   GetP()                momentum magnitude
        //Float_t   GetE()                energy (inc rest mass^2)
        //Int_t     GetStartvol()         starting volume: 10 is tank, 20 is facc, 30 is mrd
        //Int_t     GetStopvol()          stopping volume: but these may not be set.
        //Float_t   GetDir(Int_t i=0)     momentum unit vector
        //Float_t   GetPdir(Int_t i=0)    momentum vector
        //Float_t   GetStop(Int_t i=0)    stopping vertex x,y,z for i=0-2, in cm
        //Float_t   GetStart(Int_t i=0)   starting vertex x,y,z for i=0-2, in cm
        //Int_t     GetParenttype()       parent pdg, 0 for primary.
        //Float_t   GetTime()             trj->GetGlobalTime(); stopping(?) time of particle
        //Int_t     GetId()               wcsim trackid
        Float_t Primary_Vertex_X, Primary_Vertex_Y, Primary_Vertex_Z;
        if(nextrack->GetFlag()==-1){ // primary neutrino. Start is BNB, so look at stop.)
          Primary_Vertex_X=nextrack->GetStop(0);
          Primary_Vertex_Y=nextrack->GetStop(1);
          Primary_Vertex_Z=nextrack->GetStop(2);
        } else {
          Primary_Vertex_X=nextrack->GetStart(0);
          Primary_Vertex_Y=nextrack->GetStart(1);
          Primary_Vertex_Z=nextrack->GetStart(2);
        }
        cout<<"  Track "<<track<<"{"
            <<"    Flag: "<<nextrack->GetFlag()
            <<"    PDG: "<<nextrack->GetIpnu()
            <<"    ID: "<<nextrack->GetId()
            <<"    ParentPDG: "<<nextrack->GetParenttype()
            <<"    Vertex: "<<Primary_Vertex_X<<", "<<Primary_Vertex_Y<<", "<<Primary_Vertex_Z<<")"
            <<"    }"<<endl;
//        if(nextrack->GetParenttype()==0&&nextrack->GetFlag()!=-1){
//          primarycount++;
//          cout<<"event "<<i<<" trigger "<<j<<" primary "<<primarytrack<<", PDG:"<<nextrack->GetIpnu()
//              <<", trackID "<<nextrack->GetId()<<", Flag:"<<nextrack->GetFlag()<<", at ("
//              <<Primary_Vertex_X<<", "<<Primary_Vertex_Y<<", "<<Primary_Vertex_Z<<")"<<endl;
//          ++primarytrack;
//        }  // conditional on tracks representing primary particles
      }  // loop over tracks
      
      double lastdigittime=0;
      // loop over digits
      cout<<"loop over digits:"<<endl;
      for(int digiti=0; digiti<min(maxdigitstoprint,(int)ndigits); digiti++){
        WCSimRootCherenkovDigiHit* digihit=(WCSimRootCherenkovDigiHit*)(r->GetCherenkovDigiHits()->At(digiti));
        std::vector<int> truephotonindices = digihit->GetPhotonIds();
        cout<<"  digit "<<digiti<<" at time "<<digihit->GetT()<<"ns has "<<truephotonindices.size()<<" true photons"<<endl;
        for(int photoni=0; photoni<min(maxphotonsperdigittoprint,(int)truephotonindices.size()); photoni++){
          int thephotonsid = truephotonindices.at(photoni);
          WCSimRootCherenkovHitTime *thehittimeobject = 
            (WCSimRootCherenkovHitTime*)r->GetCherenkovHitTimes()->At(thephotonsid);
          Int_t thephotonsparenttrackid = thehittimeobject->GetParentID();
          cout<<"    digit "<<digiti<<", photon "<<photoni<<" has parent "<<thephotonsparenttrackid<<endl;
        }
        
//        if(digihit->GetT() < lastdigittime){ cerr<<"DIGIT TIMES ARE NOT SORTED!!!!#############"<<endl; }
//        lastdigittime=digihit->GetT();
      }
      
      // loop over true hits
      cout<<"loop over photons:"<<endl;
      for(int photoni=0; photoni<min(maxphotonstoprint,(int)numphotons); photoni++){
        WCSimRootCherenkovHitTime *thehittimeobject = 
          (WCSimRootCherenkovHitTime*)r->GetCherenkovHitTimes()->At(photoni);
        Int_t thephotonsparenttrackid = thehittimeobject->GetParentID();
        cout<<"    photon "<<photoni<<" has parent "<<thephotonsparenttrackid<<endl;
      }
      
    }  // loop over primaries (only ever 1 neutrino vertex)
  } // loop over subevents
} // loop over events

}
