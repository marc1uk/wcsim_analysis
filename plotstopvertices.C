{
//TFile* f= TFile::Open("~/LinuxSystemFiles/WCSim/gitver/root_work/in/lappd_10k_bonsai_lin_mu/ANNIEtest_10MeV_Lin_Iso_0.root");
//TFile* f= TFile::Open("/pnfs/annie/persistent/users/moflaher/wcsim_tankonly_03-05-17_BNB_World_10k_29-06-17/wcsim_0.4293.root");
TFile* f= TFile::Open("/pnfs/annie/persistent/users/moflaher/wcsim_lappd_24-09-17_BNB_Water_10k_22-05-17/wcsim_0.818.3.root");
TTree* t = (TTree*)f->Get("wcsimT");
WCSimRootEvent* e=0;
TBranch* b=0;
t->SetBranchAddress("wcsimrootevent", &e, &b);
b->GetEntry(0);
WCSimRootTrigger* r=e->GetTrigger(0);
const Float_t tank_start = 15.70;
const Float_t tank_radius = 152.4;
const Float_t tank_halfheight = 198.;
const Float_t tank_yoffset = -14.46;
//primary vertices
// world extent in WCSim is +-600cm in all directions!
TH1F *hPrimaryVtx_X = new TH1F("Event Primary VTX_X", "Event Primary VTX_X", 200, -650, 650);
TH1F *hPrimaryVtx_Y = new TH1F("Event Primary VTX_Y", "Event Primary VTX_Y", 200, -650, 650);
TH1F *hPrimaryVtx_Z = new TH1F("Event Primary VTX_Z", "Event Primary VTX_Z", 200, -650, 650);
for(int i=0; i<b->GetEntries(); i++){
  for(int j=0; j</*e->GetNumberOfEvents()*/1; j++){
    b->GetEntry(i);
    r=e->GetTrigger(j);
    int numtracks= r->GetNtrack();
    int primarytrack=0;
    int primarycount=0;
    for(int track=0; track<numtracks; track++){
      WCSimRootTrack* nextrack = (WCSimRootTrack*)r->GetTracks()->At(track);
      if(nextrack->GetFlag()!=0) continue;
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
      Float_t Primary_Vertex_X = nextrack->GetStop(0);
      Float_t Primary_Vertex_Y = nextrack->GetStop(1);
      Float_t Primary_Vertex_Z = nextrack->GetStop(2);
      hPrimaryVtx_X->Fill(Primary_Vertex_X);
      hPrimaryVtx_Y->Fill(Primary_Vertex_Y);
      hPrimaryVtx_Z->Fill(Primary_Vertex_Z);
    }  // loop over tracks
  } // loop over subevents
} // loop over events
std::string drawoption="";
TCanvas c1;
c1.cd();
hPrimaryVtx_X->SetLineColor(kBlue);
hPrimaryVtx_X->Draw(drawoption.c_str());
TCanvas c2;
c2.cd();
hPrimaryVtx_Y->SetLineColor(kBlue);
hPrimaryVtx_Y->Draw(drawoption.c_str());
TCanvas c3;
c3.cd();
hPrimaryVtx_Z->SetLineColor(kBlue);
hPrimaryVtx_Z->Draw(drawoption.c_str());
}

