{
//TFile* f=TFile::Open("/pnfs/annie/persistent/users/moflaher/wcsim/lappd/tankonly/wcsim_lappd_tankonly_24-09-17_BNB_Water_10k_22-05-17/wcsim_0.0.0.root");
//TFile* fl=TFile::Open("/pnfs/annie/persistent/users/moflaher/wcsim/lappd/tankonly/wcsim_lappd_tankonly_24-09-17_BNB_Water_10k_22-05-17/wcsim_lappd_0.0.0.root");
TFile* f=TFile::Open("/pnfs/annie/persistent/users/moflaher/wcsim/lappd/tankonly/wcsim_lappd_tankonly_03-05-17_rhatcher/wcsim_0.1000.root");
TFile* fl=TFile::Open("/pnfs/annie/persistent/users/moflaher/wcsim/lappd/tankonly/wcsim_lappd_tankonly_03-05-17_rhatcher/wcsim_lappd_0.1000.root");


int maxentriestoprint=200000000;
int maxtriggerstoprint=200000000;
int maxprimariestoprint=1; // always 1
int maxtrackstoprint=200000000;
int maxdigitstoprint=200000000;
int maxphotonsperdigittoprint=200000000;
int maxphotonstoprint=200000000;

const int TRIGGER_OFFSET=950;
const int LAPPDHITSMAX=1000;

// for /pnfs/annie/persistent/users/moflaher/wcsim/lappd/tankonly/wcsim_lappd_tankonly_24-09-17_BNB_Water_10k_22-05-17/wcsim_0.0.0.root
//TH1D hdigittimes=TH1D("hdigittimes","Times of digits within the first trigger",150,-50,100);
//TH1D hlappdtimes=TH1D("hlappdtimes","Times of lappd hits within the first trigger",150,-50,100);
//TH1D hmuontimes=TH1D("hmuontimes","Times of primary muons",150,-50,100);

// for /pnfs/annie/persistent/users/moflaher/wcsim/lappd/tankonly/wcsim_lappd_tankonly_03-05-17_rhatcher/wcsim_0.1000.root
TH1D hdigittimes=TH1D("hdigittimes","Times of digits within the first trigger",150,6660,6720);
TH1D hlappdtimes=TH1D("hlappdtimes","Times of lappd hits within the first trigger",150,6660,6720);
TH1D hmuontimes=TH1D("hmuontimes","Times of primary muons",150,6660,6720);

gStyle->SetOptTitle(1);
gStyle->SetOptStat(1);

TTree* t = (TTree*)f->Get("wcsimT");
WCSimRootEvent* e=0;
TBranch* b=0;
t->SetBranchAddress("wcsimrootevent", &e, &b);
b->GetEntry(0);
int numWCSimEntries = b->GetEntries();
WCSimRootTrigger* r=e->GetTrigger(0);
WCSimRootEventHeader* h=r->GetHeader();

TTree* lappdtree= (TTree*)fl->Get("LAPPDTree");
int numlappdshitthisevt;
std::vector<double> lappdhittimes;
std::vector<double>* lappdhittimesp = &lappdhittimes;
double lappd_numphots[LAPPDHITSMAX];
lappdtree->SetBranchStatus("*",0);
lappdtree->SetBranchStatus("lappd_numhits",1);
lappdtree->SetBranchStatus("lappdhit_stripcoort",1);
lappdtree->SetBranchStatus("lappdhit_edep",1);
lappdtree->SetBranchAddress("lappd_numhits", &numlappdshitthisevt);
lappdtree->SetBranchAddress("lappdhit_stripcoort", &lappdhittimesp);
lappdtree->SetBranchAddress("lappdhit_edep", &lappd_numphots);

int numWCSimLAPPDEntries = lappdtree->GetEntries();
assert(numWCSimLAPPDEntries==numWCSimEntries);

cout<<"This run has "<<numWCSimEntries<<" entries"<<endl;
for(int i=0; i<min(maxentriestoprint,numWCSimEntries); i++){
  //cout<<"NEW EVENT"<<endl<<"======="<<endl;
  b->GetEntry(i);
  //cout<<"This entry had "<<e->GetNumberOfEvents()<<" triggers"<<endl;
  //cout<<"This entry had "<< ( (e->HasSubEvents()) ? "1 or more" : "no" ) <<" delayed triggers"<<endl;
  //cout<<"This entry had "<<e->GetNumberOfSubEvents()<<" delayed triggers"<<endl;
  
  for(int j=0; j<min(maxtriggerstoprint,(int)e->GetNumberOfEvents()); j++){
    //cout<<"NEXT TRIGGER"<<endl;
    r=e->GetTrigger(j);
    h=r->GetHeader();
    int triggertime = h->GetDate();
    int nprimaries = r->GetNvtxs();  // always 1 - number of neutrino vertices
    //cout<<"This trigger had "<<nprimaries<<" vertex(es)"<<endl;
    
    //for(int primaryi=0; primaryi<min(maxprimariestoprint,nprimaries); primaryi++){
      Float_t True_Vertex_X=r->GetVtx(0);
      Float_t True_Vertex_Y=r->GetVtx(1);
      Float_t True_Vertex_Z=r->GetVtx(2);
      
      int numtracks= r->GetNtrack();
      int ndigits=r->GetCherenkovDigiHits()->GetEntries();
//      int numphotons = r->GetCherenkovHits()->GetEntries();
//      cout<<"event "<<i<<" trigger "<<j<<" primary "<<primaryi<<" at ("
//          <<True_Vertex_X<<", "<<True_Vertex_Y<<", "<<True_Vertex_Z<<") had:"<<endl
//          <<"Num Tracks: "<<numtracks<<endl
//          <<"Num Photons: " <<numphotons<<endl
//          <<"Num Digits: "<<ndigits<<endl;
      
      // scan through the truth tracks, find the primary muon and pull vertex info from it
      int primarymuonindex=-1;
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
        //Float_t   GetTime()             trj->GetGlobalTime(); creation time of particle
        //Int_t     GetId()               wcsim trackid
        // Implemented only in later versions of WCSim files:
        //Float_t   GetEndE()             energy on stopping of particle tracking
        //Float_t   GetEndP()             momentum on stopping of particle tracking
        //Float_t   GetPdirEnd(Int_t i=0) direction vector on stop tracking
        //Float_t   GetStopTime()         global stopping time
        //Int_t     GetId()               wcsim trackid
        
        Float_t Primary_Vertex_X, Primary_Vertex_Y, Primary_Vertex_Z, Primary_Vertex_T;
        if(nextrack->GetParenttype()==0&&nextrack->GetFlag()==0){  // extract primary particles
          if(nextrack->GetIpnu()==13){                             // in particular primary muons.
            primarymuonindex=track;
            break;
          }  // conditional on being a muon
        }    // conditional on being a primary particle
      }  // loop over tracks
      
      if(primarymuonindex<0){ continue; } // there was no primary muon in this event, skip it.
      WCSimRootTrack* murack = (WCSimRootTrack*)r->GetTracks()->At(primarymuonindex);
      hmuontimes.Fill(murack->GetTime());
      std::cout<<"Track:"<<murack->GetTime();
      
      // Assuming we have a primary muon, get the PMT hit times.
      //cout<<"loop over digits:"<<endl;
      for(int digiti=0; digiti<min(maxdigitstoprint,(int)ndigits); digiti++){
        WCSimRootCherenkovDigiHit* digihit=(WCSimRootCherenkovDigiHit*)(r->GetCherenkovDigiHits()->At(digiti));
        hdigittimes.Fill(digihit->GetT()-TRIGGER_OFFSET+triggertime);
        if(digiti==0){ std::cout<<", digit: "<<(digihit->GetT()-TRIGGER_OFFSET+triggertime); }
        //cout<<"  digit "<<digiti<<" at time "<<digihit->GetT()<<"ns"<<endl;
      }
      
    //}  // loop over primaries (only ever 1 neutrino vertex)
  } // loop over subevents
  
  // now handle LAPPD entry
  lappdtree->GetEntry(i);
  //cout<<"Looping over "<<numlappdshitthisevt<<" LAPPD digits for this event"<<endl;
  int runningcount=0;
  for(int lappdi=0; lappdi<numlappdshitthisevt; lappdi++){
    // loop over LAPPDs that had at least one hit
    int numhitsthislappd=lappd_numphots[lappdi];
    //std::cout<<"lappd "<<lappdi<<" has "<<numhitsthislappd<<" hits"<<std::endl;
    int lastrunningcount=runningcount;
    
    // loop over all the hits on this lappd
    for(;runningcount<(lastrunningcount+numhitsthislappd); runningcount++){
      double lappdhittime  = lappdhittimes.at(runningcount);
      hlappdtimes.Fill(lappdhittime);
      if((lappdi==0)&&(runningcount==lastrunningcount)){ std::cout<<", lappd: "<<(lappdhittime); }
    } // end of loop over hits on this LAPPD
    //cout<<"Done looping over hits on this LAPPD"<<endl;
  } // end of loop over LAPPDs with a hit
  std::cout<<std::endl;
  
} // loop over events

TCanvas c1;
hdigittimes.Draw();
c1.SaveAs("DigitTimesHist.png");
TCanvas c2;
hlappdtimes.Draw();
c2.SaveAs("LAPPDTimesHist.png");
TCanvas c3;
hmuontimes.Draw();
c3.SaveAs("MuonTimesHist.png");

}
