{
//TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in/wcsim_beam_test_08-04-17.root")
//TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in/wcsim_beam_test_w_genie_11-04-17.root");
//TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in/wcsim_tankonly_12-04-17_sample.root");
//TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in/wcsim_tankonly_13-04-17.root");
//TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in/wcsim_10MeV_iso_e_wDN_13-04-17.root");
TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_0.root");
TTree* t = (TTree*)f->Get("wcsimT");
WCSimRootEvent* e=0;
TBranch* b=0;
t->SetBranchAddress("wcsimrootevent", &e, &b);
b->GetEntry(0);
WCSimRootTrigger* r=e->GetTrigger(0);
WCSimRootEventHeader* h=r->GetHeader();
h->GetDate();
const Float_t tank_start = 15.70;
const Float_t tank_radius = 152.4;
const Float_t tank_halfheight = 198.;
const Float_t tank_yoffset = -14.46;
//neutrino vertices
TH1F *hTrueVtx_X = new TH1F("Event True VTX_X", "Event True VTX_X", 200, -200, 200);
TH1F *hTrueVtx_Y = new TH1F("Event True VTX_Y", "Event True VTX_Y", 200, -300, 300);
TH1F *hTrueVtx_Z = new TH1F("Event True VTX_Z", "Event True VTX_Z", 200, 0, 350);
//primary vertices
TH1F *hPrimaryVtx_X = new TH1F("Event Primary VTX_X", "Event Primary VTX_X", 200, -200, 200);
TH1F *hPrimaryVtx_Y = new TH1F("Event Primary VTX_Y", "Event Primary VTX_Y", 200, -300, 300);
TH1F *hPrimaryVtx_Z = new TH1F("Event Primary VTX_Z", "Event Primary VTX_Z", 200, 0, 350);
for(int i=0; i<b->GetEntries(); i++){
for(int j=0; j</*e->GetNumberOfEvents()*/1; j++){
b->GetEntry(i);
r=e->GetTrigger(j);
h=r->GetHeader();
int nprimaries = r->GetNvtxs();	// always 1
//for(int primaryi=0; primaryi<nprimaries; primaryi++){
Float_t vertex_scaling_factor=1.;	// 100.
Float_t True_Vertex_X=r->GetVtx(0)*vertex_scaling_factor;
Float_t True_Vertex_Y=r->GetVtx(1)*vertex_scaling_factor;
Float_t True_Vertex_Z=r->GetVtx(2)*vertex_scaling_factor;
cout<<"event "<<i<<" trigger "<<j<</*" primary "<<primaryi<<*/" at ("<<True_Vertex_X<<", "<<True_Vertex_Y<<", "<<True_Vertex_Z<<")"<<endl;
//if(std::to_string(True_Vertex_X).substr(0,5)=="-99.9" || TMath::Abs(True_Vertex_X)<0.00000001){ // we don't have neutrino vtx info, need to load tracks:
//cout<<"True_Vertex_X is "<<True_Vertex_X<<", replacing"<<endl;
int numtracks= r->GetNtrack();
// scan through the truth tracks, find the primary muon and pull vertex info from it
int primarytrack=0;
// first count the number of tracks and weight the primaries by this number.
int primarycount=0;
for(int track=0; track<numtracks; track++){
  WCSimRootTrack* nextrack = (WCSimRootTrack*)r->GetTracks()->At(track);
  if(nextrack->GetParenttype()==0&&nextrack->GetFlag()!=-1){ primarycount++; }
}
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
  if(nextrack->GetParenttype()==0){
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
    cout<<"event "<<i<<" trigger "<<j<<" primary "<<primarytrack<<", PDG:"<<nextrack->GetIpnu()<<", trackID "<<nextrack->GetId()<<", Flag:"<<nextrack->GetFlag()<<", at ("<<Primary_Vertex_X<<", "<<Primary_Vertex_Y<<", "<<Primary_Vertex_Z<<")"<<endl;
    Float_t track_scaling_factor=1.;
    /* BEAM ONLY -  NEEDS FIXING 
    if(nextrack->GetFlag()==-2){ 
      track_scaling_factor=100.;
    } else {
      track_scaling_factor=0.1;
    }
    */
    Primary_Vertex_X= Primary_Vertex_X*track_scaling_factor;
    Primary_Vertex_Y= Primary_Vertex_Y*track_scaling_factor;
    Primary_Vertex_Z= Primary_Vertex_Z*track_scaling_factor;
    // ^^^^^^^^**************************!!!
    hPrimaryVtx_X->Fill(Primary_Vertex_X, (1./primarycount));
    hPrimaryVtx_Y->Fill(Primary_Vertex_Y, (1./primarycount)); //+tank_halfheight
    hPrimaryVtx_Z->Fill(Primary_Vertex_Z, (1./primarycount)); //-tank_start-tank_radius
    ++primarytrack;
  }  // conditional on tracks representing primary particles
}  // loop over tracks
//} else {
  //cout<<"True_Vertex_X is "<<True_Vertex_X<<" !=-99.9, NOT replacing"<<endl;
    hTrueVtx_X->Fill(True_Vertex_X);
    hTrueVtx_Y->Fill(True_Vertex_Y); //+tank_halfheight
    hTrueVtx_Z->Fill(True_Vertex_Z); //-tank_start-tank_radius
//}  // conditional on genie info not loaded, use tracks instead. - now plot both. 
//}  // loop over primaries (only ever 1 neutrino vertex)
} // loop over subevents
} // loop over events
TCanvas c1;
c1.cd();
hPrimaryVtx_X->SetLineColor(kBlue);
hTrueVtx_X->SetLineColor(kRed);
string drawoption="HIST SAME";
hTrueVtx_X->Draw();
hPrimaryVtx_X->Draw(drawoption.c_str());
TCanvas c2;
c2.cd();
hPrimaryVtx_Y->SetLineColor(kBlue);
hTrueVtx_Y->SetLineColor(kRed);
hTrueVtx_Y->Draw();
hPrimaryVtx_Y->Draw(drawoption.c_str());
TCanvas c3;
c3.cd();
hPrimaryVtx_Z->SetLineColor(kBlue);
hTrueVtx_Z->SetLineColor(kRed);
hTrueVtx_Z->Draw();
hPrimaryVtx_Z->Draw(drawoption.c_str());
}

