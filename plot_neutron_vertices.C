{
TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_0.root");
TTree* t = (TTree*)f->Get("wcsimT");
WCSimRootEvent* e=0;
TBranch* b=0;
t->SetBranchAddress("wcsimrootevent", &e, &b);
b->GetEntry(0);
WCSimRootTrigger* r=e->GetTrigger(0);
for(int i=0; i<b->GetEntries(); i++){
for(int j=0; j<e->GetNumberOfEvents(); j++){
b->GetEntry(i);
r=e->GetTrigger(j);
int numtracks= r->GetNtrack();
// scan through the truth tracks, find neutrons and get their vertices
for(int track=0; track<numtracks; track++){
  WCSimRootTrack* nextrack = (WCSimRootTrack*)r->GetTracks()->At(track);
//          Int_t     GetIpnu()             pdg
//          Int_t     GetFlag()             -1: neutrino primary, -2: neutrino target, 0: other
//          Float_t   GetM()                mass
//          Float_t   GetP()                momentum magnitude
//          Float_t   GetE()                energy (inc rest mass^2)
//          Int_t     GetStartvol()         starting volume: 10 is tank, 20 is facc, 30 is mrd
//          Int_t     GetStopvol()          stopping volume: but these may not be set.
//          Float_t   GetDir(Int_t i=0)     momentum unit vector
//          Float_t   GetPdir(Int_t i=0)    momentum vector
//          Float_t   GetStop(Int_t i=0)    stopping vertex x,y,z for i=0-2, in cm
//          Float_t   GetStart(Int_t i=0)   starting vertex x,y,z for i=0-2, in cm
//          Int_t     GetParenttype()       parent pdg, 0 for primary.
//          Float_t   GetTime()             trj->GetGlobalTime(); stopping(?) time of particle
//          Int_t     GetId()               wcsim trackid
  if(nextrack->GetIpnu()!=2112) continue;
    Float_t Start_Vertex_X=nextrack->GetStart(0);
    Float_t Start_Vertex_Y=nextrack->GetStart(1);
    Float_t Start_Vertex_Z=nextrack->GetStart(2);
    Float_t Stop_Vertex_X=nextrack->GetStop(0);
    Float_t Stop_Vertex_Y=nextrack->GetStop(1);
    Float_t Stop_Vertex_Z=nextrack->GetStop(2);
    
    cout<<"event "<<i<<" trigger "<<j<<" track "<<track<<" started at ("<<Start_Vertex_X<<", "
        <<Start_Vertex_Y<<", "<<Start_Vertex_Z<<") and ended at ("<<Stop_Vertex_X<<", "<<Stop_Vertex_Y
        <<", "<<Stop_Vertex_Z<<")"<<endl;
    /* BEAM ONLY -  NEEDS FIXING 
    if(nextrack->GetFlag()==-2){ 
      track_scaling_factor=100.;
    } else {
      track_scaling_factor=0.1;
    }
    */
}  // loop over tracks
} // loop over subevents
} // loop over events
}
