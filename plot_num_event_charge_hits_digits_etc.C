{
//TFile* f=TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in/10k_1000_photons_bomb_423nm_23-04-17.root");
TFile* f=TFile::Open("/pnfs/annie/persistent/users/moflaher/wcsim_tankonly_03-05-17/wcsim_0.1000.root");
TTree* t = (TTree*)f->Get("wcsimT");
TTree* t2= (TTree*)f->Get("wcsimGeoT");
WCSimRootGeom* geo =0;
t2->SetBranchAddress("wcsimrootgeom",&geo);
t2->GetEntry(0);
int numpmts=geo->GetWCNumPMT();
cout<<"the geometry has "<<numpmts<<" PMTs"<<endl;
WCSimRootEvent* ev;
TBranch* evb =0;
t->SetBranchAddress("wcsimrootevent",&ev,&evb);
WCSimRootTrigger* tr;
TH1D hnumtubeshit=TH1D("hnumtubeshit","Number of PMTs with >0 Digits",300,0,300);
TH1D hnumphots=TH1D("hnumphots","Number Of Photon Hits per Event",100,0,100);
TH1D hnumdigits=TH1D("hnumdigits","Number Of Digits per Event",100,0,100);
TH1D hnumdigitizedphots=TH1D("hnumdigitizedphots","Number of Digitized Photons per Event",100,0,100);
TH1D hnumdigitizednoise=TH1D("hnumdigitizednoise","Number of Digitized Noise Hits per Event",5,0,5);
TH1D hnumphotsperdigit=TH1D("hnumphotsperdigit","Number Of Photons per Digit",10,0,10);
TH1D htotalchargeperdigit=TH1D("htotalchargeperdigit","Amount of Charge per Digit",100,0,30);
TH1D htotalchargeinevent=TH1D("htotalchargeinevent","Total Charge of All Digits in Event",100,0,100);
//TH1D hnumpesperevent=TH1D("hnumpesperevent","Number Of PEs per Event",100,0,200);
for(int i=0; i<100/*t->GetEntries()*/; i++){
evb->GetEntry(i);
  tr=ev->GetTrigger(0);
  Int_t numphots=tr->GetNcherenkovhittimes();
  hnumphots.Fill(numphots);
  Int_t numdigs = tr->GetNcherenkovdigihits();
  hnumdigits.Fill(numdigs);
  double totalchargeinevent=tr->GetSumQ();
  htotalchargeinevent.Fill(totalchargeinevent);
  int numtubeshitthisevent=tr->GetNumTubesHit();
  hnumtubeshit.Fill(numtubeshitthisevent);
  int totaldigitizedphotonsthisevent=0;
  int numnoisephots=0;
  for(int digii=0; digii<numdigs; digii++){
    WCSimRootCherenkovDigiHit* digihit=(WCSimRootCherenkovDigiHit*)(tr->GetCherenkovDigiHits()->At(digii));
    double thedigitcharge = digihit->GetQ();
    htotalchargeperdigit.Fill(thedigitcharge);
    std::vector<int> truephotonindices = digihit->GetPhotonIds();
    int numphotsinthisdigit = truephotonindices.size();
    hnumphotsperdigit.Fill(numphotsinthisdigit);
    totaldigitizedphotonsthisevent+=numphotsinthisdigit;
    for(int aphotohit=0; aphotohit<numphotsinthisdigit; aphotohit++){
      WCSimRootCherenkovHitTime* thehittime = (WCSimRootCherenkovHitTime*)(tr->GetCherenkovHitTimes()->At(aphotohit));
      int theparent=thehittime->GetParentID();
      //cout<<"hit "<<aphotohit<<" had parent "<<theparent<<endl;
      if(theparent==-1) numnoisephots++;
    }
    hnumdigitizednoise.Fill(numnoisephots);
  }
  hnumdigitizedphots.Fill(totaldigitizedphotonsthisevent);
}
TCanvas c1("c1","canv1",600,600);
hnumtubeshit.Draw();
TCanvas c2("c2","canv2",600,600);
hnumphots.Draw();
TCanvas c3("c3","canv3",600,600);
hnumdigits.Draw();
TCanvas c4("c4","canv4",600,600);
hnumdigitizedphots.Draw();
TCanvas c5("c5","canv5",600,600);
hnumphotsperdigit.Draw();
TCanvas c6("c6","canv6",600,600);
htotalchargeperdigit.Draw();
TCanvas c7("c7","canv7",600,600);
htotalchargeinevent.Draw();
//TCanvas c8("c8","canv8",600,600);
//hnumpesperevent.Draw();
TCanvas c9("c9","canv9",600,600);
hnumdigitizednoise.Draw();
}
