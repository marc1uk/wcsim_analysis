{
//TFile* f= TFile::Open("/home/marc/anniegpvm/wcsim/root_work/in/energy_comparison/ANNIEtest_5MeV_mu-_Uni_Iso_0.root");
//TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in/ANNIEtest_200MeV_mu-_Uni_Iso_0.root");
TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/ANNIEtest_MRD_muon_sample_0.root");
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
// scan through the truth tracks, find primary muons and print their direction vectors
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
  if(nextrack->GetParenttype()==0&&nextrack->GetIpnu()==13){
    cout<<"DIR=("<<nextrack->GetDir(0)<<","<<nextrack->GetDir(1)<<","<<nextrack->GetDir(2)<<")"<<endl;
    cout<<"E="<<nextrack->GetE()<<endl;
    cout<<"Track Length="<<TMath::Sqrt(TMath::Power(nextrack->GetStop(0)-nextrack->GetStart(0),2)+
    					TMath::Power(nextrack->GetStop(1)-nextrack->GetStart(1),2)+
    					TMath::Power(nextrack->GetStop(2)-nextrack->GetStart(2),2))<<endl;
  }
}  // loop over tracks
} // loop over subevents
} // loop over events
}
