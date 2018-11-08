{
//TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_1k_thermal_neutron_test_scat_0.root");
//TFile* f= TFile::Open("/annie/app/users/moflaher/wcsim/build/wcsim_photon_bomb_0.root");
TFile* f = TFile::Open("/pnfs/annie/persistent/users/moflaher/wcsim/lappd/tankonly/wcsim_lappd_tankonly_24-09-17_BNB_Water_10k_22-05-17/wcsim_0.0.0.root");
//TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_0.root");

// Trigger times and CherenkovHitTimes (true photon times) are both in absolute event ns, 
// only digit times are relative, to the start of the trigger and with the historic offset
// true photon times are recorded independently of digits or triggers.

TTree* t = (TTree*)f->Get("wcsimT");
WCSimRootEvent* ev;
TBranch* evb =0;
t->SetBranchAddress("wcsimrootevent",&ev,&evb);
WCSimRootTrigger* tr;

TTree* geotree = (TTree*)f->Get("wcsimGeoT");
if(geotree==0){ cerr<<"NO GEOMETRY"<<endl; assert(false); }
WCSimRootGeom* geo = 0;
geotree->SetBranchAddress("wcsimrootgeom", &geo);
if (geotree->GetEntries() == 0) { cerr<<"geotree has no entries!"<<endl; exit(9); }
geotree->GetEntry(0);
int numpmts = geo->GetWCNumPMT();

gStyle->SetOptTitle(1);

TH1D phhist=TH1D("phhist","TrueTime of Photons;TrueTime (ns);Num Events",100,0,30);
double mintruetime=0;
double maxtruetime=0;

int maxentriestoprint=10;
//int maxtriggerstoprint=1; TODO
Int_t maxdigitstoprint=10;
Int_t maxphotstoprint=10;

Int_t treeentries = t->GetEntries();
for(int eventi=0; eventi<std::min(maxentriestoprint,treeentries); eventi++){
  evb->GetEntry(eventi);
  tr=ev->GetTrigger(0);
  cout<<"Trigger time "<<tr->GetHeader()->GetDate()<<endl;
  Int_t numphots = tr->GetCherenkovHits()->GetEntries();
  Int_t numphottimes = tr->GetCherenkovHitTimes()->GetEntriesFast();
  Int_t numdigs = tr->GetCherenkovDigiHits()->GetEntries();
  
  std::map<int,int> uniquepmtsmap{};
  for(int digiti=0; digiti<std::min(maxdigitstoprint,numdigs); digiti++){
    WCSimRootCherenkovDigiHit* digihit=(WCSimRootCherenkovDigiHit*)(tr->GetCherenkovDigiHits()->At(digiti));
    //WCSimRootChernkovDigiHit has methods GetTubeId(), GetT(), GetQ()
    cout<<"  Digit Time "<<digihit->GetT()<<endl;
    std::vector<int> truephotonindices = digihit->GetPhotonIds();
    
    int photoni=0;
    for( auto hittimeindex : truephotonindices){
      WCSimRootCherenkovHitTime* hittime = (WCSimRootCherenkovHitTime*)tr->GetCherenkovHitTimes()->At(hittimeindex);
      // WCSimRootCherenkovHitTime has methods GetTruetime() and GetParentID();
      double truetime = hittime->GetTruetime();
      cout<<"    Photon Time "<<truetime<<endl;
      if(truetime<mintruetime) mintruetime=truetime;
      if(truetime>maxtruetime) maxtruetime=truetime;
      phhist.Fill(truetime);
      photoni++;
      if(photoni>maxphotstoprint) break;
    }
  }
  
//  for(Int_t photoni=0; photoni<std::min(maxphotstoprint,numphots); photoni++){
//    WCSimRootCherenkovHit* hit = (WCSimRootCherenkovHit*)tr->GetCherenkovHits()->At(photoni);
//    //WCSimRootCherenkovHit has methods GetTubeId(), GetTotalPe(int). only really need int=0.
//    int tubeid = hit->GetTubeID();
//    if(hitsperpmt.count(tubeid)) hitsperpmt.at(tubeid)++;
//    else hitsperpmt.emplace(tubeid,1);
//    //int hittimeindex = hit->GetTotalPe(0);
//  }
  
  ev->ReInitialize();
}

cout<<"mintruetime="<<mintruetime<<endl;
cout<<"maxtruetime="<<maxtruetime<<endl;

TCanvas c1;
phhist.Draw();
}
