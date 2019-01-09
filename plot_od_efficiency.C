{
//TFile* f=TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/out/wcsim_odtest.root");
TFile* f=TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_0.root");
TTree* t = (TTree*)f->Get("wcsimT");
WCSimRootEvent* e=0;
TBranch* b=0;
t->SetBranchAddress("wcsimrootevent", &e, &b);
b->GetEntry(0);
WCSimRootTrigger* r=e->GetTrigger(0);
WCSimRootEventHeader* h=r->GetHeader();
h->GetDate();

TTree* geotree = (TTree*)f->Get("wcsimGeoT");
if(geotree==0){ cerr<<"NO GEOMETRY"<<endl; assert(false); }
WCSimRootGeom* geo = nullptr;
if (geotree->GetEntries() == 0) { cerr<<"geotree has no entries!"<<endl; exit(9); }
geotree->SetBranchAddress("wcsimrootgeom", &geo);
geotree->GetEntry(0);
if(geo==nullptr) { cerr<<"retrieved geo is still null!"<<endl; exit(9); }

TH2D numhitsvsemu("numhitsvsemu","Num OD Hits vs Muon Energy",10,0,10,100,0,2000);
TH2D odchargevsemu("odchargevsemu","Total OD Charge vs Muon Energy",100,0,100,100,0,2000);

int bignum = std::numeric_limits<int>::max();
int maxentriestoprint=bignum;
int maxtriggerstoprint=bignum;
int maxprimariestoprint=1; // always 1
int maxtrackstoprint=bignum;
int maxdigitstoprint=bignum;
int maxphotonsperdigittoprint=0;
int maxphotonstoprint=0;
//neutrino vertices
cout<<"This run has "<<b->GetEntries()<<" entries"<<endl;
for(int i=0; i<min(maxentriestoprint,(int)b->GetEntries()); i++){
  //cout<<"NEW EVENT"<<endl<<"======="<<endl;
  b->GetEntry(i);
  //cout<<"This entry had "<<e->GetNumberOfEvents()<<" triggers"<<endl;
  //cout<<"This entry had "<< ( (e->HasSubEvents()) ? "1 or more" : "no" ) <<" delayed triggers"<<endl;
  //cout<<"This entry had "<<e->GetNumberOfSubEvents()<<" delayed triggers"<<endl;
  int primmuidx=-1;
  double maxprime=-1;
  for(int j=0; j<min(maxtriggerstoprint,(int)e->GetNumberOfEvents()); j++){
    //cout<<"NEXT TRIGGER"<<endl;
    r=e->GetTrigger(j);
    h=r->GetHeader();
    //cout<<" >>> Trigger time was : "<<h->GetDate()<<endl; // TRIGGER TIME OF 0 MEANS NO NDIGITS TRIGGER
    int nprimaries = r->GetNvtxs();  // always 1
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
      for(int tracki=0; tracki<min(maxtrackstoprint,numtracks); tracki++){
        WCSimRootTrack* nextrack = (WCSimRootTrack*)r->GetTracks()->At(tracki);
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
        if(nextrack->GetFlag()!=0){ // primary neutrino. Start is BNB, so look at stop.)
//          std::cout<<"Primary "<<nextrack->GetE()<<"MeV PDG("<<nextrack->GetIpnu()
//                   <<") with flag "<<nextrack->GetFlag()<<" starting at ("
//                   <<nextrack->GetStart(0)<<", "
//                   <<nextrack->GetStart(1)<<", "
//                   <<nextrack->GetStart(2)<<") and travelling to ("
//                   <<nextrack->GetStop(0)<<", "
//                   <<nextrack->GetStop(1)<<", "
//                   <<nextrack->GetStop(2)<<")"<<std::endl;
        } else if(nextrack->GetIpnu()==13){
          if(nextrack->GetE()>maxprime){
            primmuidx = tracki;
          }
        }
      }  // loop over tracks
      if(primmuidx<0){
        std::cerr<<"Could not find any primary muon!!"<<std::endl;
        continue;
      }
      
      WCSimRootTrack* nextrack = (WCSimRootTrack*)r->GetTracks()->At(primmuidx);
      std::cout<<"Primary "<<nextrack->GetE()<<"MeV Muon selected starts at ("
               <<nextrack->GetStart(0)<<", "
               <<nextrack->GetStart(1)<<", "
               <<nextrack->GetStart(2)<<") and travelling to ("
               <<nextrack->GetStop(0)<<", "
               <<nextrack->GetStop(1)<<", "
               <<nextrack->GetStop(2)<<")"<<std::endl;
      double muE = nextrack->GetE();
      
      // loop over digits
      int numodhits=0;
      double totdigicharge=0;
      //cout<<"loop over digits:"<<endl;
      for(int digiti=0; digiti<min(maxdigitstoprint,(int)ndigits); digiti++){
        WCSimRootCherenkovDigiHit* digihit=(WCSimRootCherenkovDigiHit*)(r->GetCherenkovDigiHits()->At(digiti));
        int tubeid = digihit->GetTubeId();
        int tubetypeindex = geo->GetTubeIndex(tubeid);
        std::string PMT_type = geo->GetWCPMTNameAt(tubetypeindex);
        //std::cout<<"geo PMT_type was "<<PMT_type;
        WCSimRootPMT thepmt;
        try{
          thepmt = geo->GetPMT(tubeid);
        } catch (...){
          std::cerr<<"Could not geo->GetPMT("<<tubeid<<")"<<std::endl;
        }
        // alt way, check they agree
        std::string PMT_type_2 = thepmt.GetName();
        if(PMT_type!=PMT_type_2){
          std::cerr<<"PMT type name doesn't match! From PMT: "<<PMT_type_2<<", from geo: "<<PMT_type<<std::endl;
        }
        //std::cout<<", pmt PMT_type_2 was "<<PMT_type_2<<std::endl;
        
        if(PMT_type=="EMI9954KB"){
          numodhits++;
          totdigicharge += digihit->GetQ();
          
          std::vector<int> truephotonindices = digihit->GetPhotonIds();
          //cout<<"  digit "<<digiti<<" at time "<<digihit->GetT()<<"ns has "
          //    <<truephotonindices.size()<<" true photons"<<endl;
          for(int photoni=0; photoni<min(maxphotonsperdigittoprint,(int)truephotonindices.size()); photoni++){
            int thephotonsid = truephotonindices.at(photoni);
            WCSimRootCherenkovHitTime *thehittimeobject = 
              (WCSimRootCherenkovHitTime*)r->GetCherenkovHitTimes()->At(thephotonsid);
            Int_t thephotonsparenttrackid = thehittimeobject->GetParentID();
            //cout<<"    digit "<<digiti<<", photon "<<photoni<<" has parent "<<thephotonsparenttrackid<<endl;
          }
        }
      }
      
      numhitsvsemu.Fill(numodhits,muE);
      odchargevsemu.Fill(totdigicharge,muE);
      
      // loop over true hits
      //cout<<"loop over photons:"<<endl;
      for(int photoni=0; photoni<min(maxphotonstoprint,(int)numphotons); photoni++){
        WCSimRootCherenkovHitTime *thehittimeobject = 
          (WCSimRootCherenkovHitTime*)r->GetCherenkovHitTimes()->At(photoni);
        Int_t thephotonsparenttrackid = thehittimeobject->GetParentID();
        cout<<"    photon "<<photoni<<" has parent "<<thephotonsparenttrackid<<endl;
      }
      
    }  // loop over primaries (only ever 1 neutrino vertex)
  } // loop over subevents
} // loop over events

TCanvas c1;
numhitsvsemu.Draw("colz");
TCanvas c2;
odchargevsemu.Draw("colz");

}
