{
//TFile* f=TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/ANNIEtest_10MeV_Uni_Iso_0.root");
//TFile* f=TFile::Open("/home/marc/LinuxSystemFiles/Bonsai/createlike_nickwp/analysis_wcsim/ANNIEtest_10MeV_Uni_Iso_0.root");
//TFile* f=TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in/wcsim_tankonly_13-04-17.root");
//TFile* f=TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in/wcsim_10MeV_iso_e_wDN_13-04-17.root");
//TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_1k_photon_bomb_0.root");
//TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_1k_thermal_neutron_test_0.root");
//TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_1k_5MeV_electron_test_0.root");
//TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_1k_photon_bomb_test_0.root");
//TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_1k_thermal_neutron_test_scat_0.root");
//TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/Bonsai/validation/Mahdi_Bonsaifiles/e5d1000PMT3.root");
//TFile* f= TFile::Open("/annie/app/users/moflaher/wcsim/build/wcsim_photon_bomb_0.root");
//TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_0_rindex15.root");
TFile* f= TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_0.root");

#define FILE_VERSION 0

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
TH1D phhist=TH1D("phhist","Number Of Photon Hits per Event;Num Photon Hits;Events",100,0,1000);
TH1D dhhist=TH1D("dhhist","Number Of Digits per Event;Number of Digits;Events",100,0,300);
TH1D qhist = TH1D("qhist","Total Charge per Event;Total Charge;Events",100,0,1000);
TH1D pdhist=TH1D("pdhist","Number Of Photons per Digit;Photons per Digit;Digits",20,0,10);
TH1D pehist=TH1D("pehist","Charge per Digit;Charge;Digits",100,0,10);
TH1D tchist=TH1D("tchist","Num PMTs Hit per Event;Num PMTs Hit;Events",250,0,250);
TH1D tdhist=TH1D("tdhist","Num PMTs with a Digit per Event;Num PMTs With a Digit;Events",250,0,250);
TH1D hpphist=TH1D("hpphist","Num Photon Hits per PMT;Photons per PMT;Events",20,0,10);
TH1D dphist=TH1D("dphist","Num Direct Photon Hits per Event;Num Direct Photons;Events",100,0,100);
TH1D dshist=TH1D("dshist","Num Scattered Photon Hits per Event;Num Scattered Photons;Events",100,0,100);
for(int eventi=0; eventi<t->GetEntries(); eventi++){
  evb->GetEntry(eventi);
  tr=ev->GetTrigger(0);
  //Int_t numphots=tr->GetNcherenkovhittimes();
  Int_t numphots = tr->GetCherenkovHits()->GetEntries();
  Int_t numphottimes = tr->GetCherenkovHitTimes()->GetEntriesFast();
  //Int_t numdigs = tr->GetNcherenkovdigihits();
  Int_t numdigs = tr->GetCherenkovDigiHits()->GetEntries();
  if(numphots!=numphottimes){
    //cerr<<"NUM PHOTON TIMES ("<<numphottimes<<") != NUM PHOTON HITS ("<<numphots<<")?!"<<endl;
    //cerr<<"c.f. Num digits: "<<numdigs<<endl;
    // WCRawPMTSignalCollection, containing all true photon hits in the entire event, plus noise,
    // is saved in the WCSimRootEvent::CherenkovHitTimes collection. This represents every photon on
    // every PMT, signal and noise, before any rejection.
    // the WCRawPMTSignalCollection is also used to fill the WCSimRootEvent::CherenkovHits collection.
    // This represents the compression of the CherenkovHitTimes down to one entry per unique PMT hit.
    // The 'WCSimRootTrigger->NumTubesHit therefore represents the number of unique PMTs with at least one
    // photon, real or noise, before any rejection.
    // CherenkovHitTimes store a flattened array of photon hit times and parent ids, 
    // CherenkovHits stores a compressed array of tube ids and number of photons on that tube.
    // CherenkovDigiHits stores an array of all digits (multiple photon hits combined via digitizer integration)
    // that pass a series of rejections - digitizer deadtime, pe threshold, and trigger window.
    //
    // this means:
    // hits collection does not store all hits in hittimes collection: numphottimes > numphots
    // ^ because there are as many hittimes entries as photons on that unique pmt
    //
    // There are also sometimes fewer digits than unique PMTs with a 'CherenkovHit'. 
    // ^ because not all PMTs with a CherenkovHit (one or more photons) generate a digit.
    // Either the Hit may generate a digit, but the digit is not stored as it's outside any trigger window
    // Or the Hit was rejected, because it fell in the digitizer deadtime, or it was rejected by pe threshold.
  }
  
  //int numpmtshit = tr->GetNumTubesHit(); << THIS VARIABLE IS NOT WHAT IT SAYS IT IS.
  std::map<int,int> uniquepmtsmap{};
  int numdirectphots=0, numscatteredphots=0;
  for(int digiti=0; digiti<numdigs; digiti++){
    WCSimRootCherenkovDigiHit* digihit=(WCSimRootCherenkovDigiHit*)(tr->GetCherenkovDigiHits()->At(digiti));
    //WCSimRootChernkovDigiHit has methods GetTubeId(), GetT(), GetQ()
    pehist.Fill(digihit->GetQ());
    int thetubeid = digihit->GetTubeId();
    if(not uniquepmtsmap.count(thetubeid)) uniquepmtsmap.emplace(thetubeid,1);
    std::vector<int> truephotonindices = digihit->GetPhotonIds();
    pdhist.Fill(truephotonindices.size());
    
#if FILE_VERSION>2
    std::cout<<"scanning for scattered photons"<<std::endl;
    int photoni=0;
    for( auto hittimeindex : truephotonindices){
      WCSimRootCherenkovHitTime* hittime = 
        (WCSimRootCherenkovHitTime*)tr->GetCherenkovHitTimes()->At(hittimeindex);
      // WCSimRootCherenkovHitTime has methods GetTruetime() and GetParentID();
      // and when added, GetNumScatterings and GetScatterings();
      if(hittime->GetNumScatterings()==0){ numdirectphots++; }
      else {
        numscatteredphots++;
        if(eventi<10){
          std::map<std::string,int> scatteringevents = hittime->GetScatterings();
          cout<<"event "<<eventi<<" digit "<<digiti<<" photon "<<photoni<<" had "
              <<hittime->GetNumScatterings()<<" scatterings: ";
          for(auto ascatteringtype : scatteringevents){
            cout<<"{"<<ascatteringtype.first<<":"<<ascatteringtype.second<<"} ";
          }
          cout<<endl;
        }
      }
      photoni++;
    }
#else
  if(eventi==0&&digiti==0){
    std::cout<<"Skipping plots of scattered & direct photons as FILE_VERSION<2"<<std::endl;
  }
#endif
  }
  
  std::map<int,int> hitsperpmt{};
  for(Int_t photoni=0; photoni<numphots; photoni++){
    WCSimRootCherenkovHit* hit = (WCSimRootCherenkovHit*)tr->GetCherenkovHits()->At(photoni);
    //WCSimRootCherenkovHit has methods GetTubeId(), GetTotalPe(int). only really need int=0.
    int tubeid = hit->GetTubeID();
    if(hitsperpmt.count(tubeid)) hitsperpmt.at(tubeid)++;
    else hitsperpmt.emplace(tubeid,1);
    //int hittimeindex = hit->GetTotalPe(0);
  }
  std::vector<int> paddedhitsperpmt(numpmts,0);
  for(int tubeid=0; tubeid<numpmts; tubeid++){
    if(hitsperpmt.count(tubeid)) paddedhitsperpmt.at(tubeid)= hitsperpmt.at(tubeid);
  }
  
  //cerr<<"num unique pmts with Hits: "<<hitsperpmt.size()<<endl;
  //cout<<"event "<<i<<" had "<<numphots<<" photons and "<<numdigs<<" digits"<<endl;
  phhist.Fill(numphottimes);
  dhhist.Fill(numdigs);
  qhist.Fill(tr->GetSumQ()); // checks out ok.
  tchist.Fill(hitsperpmt.size()); // same as numphottimes
  tdhist.Fill(uniquepmtsmap.size());
  dphist.Fill(numdirectphots);
  dshist.Fill(numscatteredphots);
  ev->ReInitialize();
}
TCanvas c1;
phhist.Draw();
TCanvas c2;
dhhist.Draw();
TCanvas c3;
qhist.Draw();
TCanvas c4;
pdhist.Draw();
TCanvas c5;
pehist.Draw();
TCanvas c6;
tchist.Draw();
TCanvas c7;
tdhist.Draw();
#if FILE_VERSION>2
TCanvas c8;
dphist.Draw();
TCanvas c9;
dshist.Draw();
#endif
}
