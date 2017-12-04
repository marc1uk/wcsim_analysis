{
//TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_0.root");
//TTree* t = (TTree*)f->Get("wcsimT");
//TFile* f= TFile::Open("/pnfs/annie/persistent/users/moflaher/wcsim_particlegun_03-05-17/wcsim_0.root");
TChain* t= new TChain("wcsimT");
//t->Add("/pnfs/annie/persistent/users/moflaher/wcsim_particlegun_03-05-17/wcsim_[0-9].root");
t->Add("/pnfs/annie/persistent/users/moflaher/wcsim_tankonly_03-05-17_rhatcher/wcsim_100[0-9].root");
cout<<"loaded "<<t->GetEntries()<<" chain entries"<<endl;
t->LoadTree(0);
WCSimRootEvent* e=0;
TBranch* b=0;
TH3F* ahist=new TH3F("ahist","MuonStop",100,-500,500,100,-500,500,100,00,1000);
int treenum=-1;
for(int i=0; i<t->GetEntries(); i++){
int localentry = t->LoadTree(i);
int nextreenum=t->GetTreeNumber();
if(nextreenum!=treenum){
t->SetBranchAddress("wcsimrootevent", &e, &b);
treenum=nextreenum;
}
b->GetEntry(localentry);
for(int j=0; j<e->GetNumberOfEvents(); j++){
WCSimRootTrigger* r=e->GetTrigger(j);
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
  if((nextrack->GetIpnu()!=/*2112*/13)||(nextrack->GetFlag()!=0)/*||(nextrack->GetParenttype()!=0)*/) continue;
    Float_t Start_Vertex_X=nextrack->GetStart(0);
    Float_t Start_Vertex_Y=nextrack->GetStart(1);
    Float_t Start_Vertex_Z=nextrack->GetStart(2);
    Float_t Stop_Vertex_X=nextrack->GetStop(0);
    Float_t Stop_Vertex_Y=nextrack->GetStop(1);
    Float_t Stop_Vertex_Z=nextrack->GetStop(2);
    ahist->Fill(Stop_Vertex_X,Stop_Vertex_Y,Stop_Vertex_Z);
    
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
TCanvas c1;
ahist->Draw();
}
