// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THStack.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TPolyMarker3D.h"
#include "TPointSet3D.h"   // use TPointSet3D over TPolyMarker3D for ogl viewer (needed for gdml overlay)
#include "TGeoManager.h"
#include "TGLViewer.h"

// std headers
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <limits>

// WCSim headers
#include "../wcsim/include/WCSimRootEvent.hh"
#include "../wcsim/include/WCSimRootGeom.hh"
#include "../wcsim/include/WCSimPmtInfo.hh"
#include "../wcsim/include/WCSimLAPPDInfo.hh"
#include "../wcsim/include/WCSimEnumerations.hh"
#include "../wcsim/include/WCSimRootLinkDef.hh"
#include "../wcsim/include/WCSimRootOptions.hh"

// Other
#include "ColourWheel.hh"

#define FILE_VER 5

// setup static members for ColourWheel
const std::vector<EColor> ColourWheel::colours{kBlack, kBlue, (EColor)TColor::GetColorDark(kGreen), kRed, kViolet, kOrange, kMagenta,(EColor)(kAzure+2),(EColor)(kOrange+4),(EColor)(kViolet-6),(EColor)(kTeal-6)};
const std::vector<std::string> ColourWheel::colournames{"kBlack", "kBlue", "kGreen", "kRed", "kViolet", "kOrange", "kMagenta","kAzure","kOrange","kViolet","kTeal"};
void FillTankMapHist(WCSimRootGeom* geo, int tubeID, std::map<std::string, TH2D*> &maphistos, double weight=1);
void MakePMTmap(WCSimRootGeom* geo, std::map<int, std::pair<int,int> > &topcappositionmap, std::map<int, std::pair<int,int> > &bottomcappositionmap, std::map<int, std::pair<int,int> > &wallpositionmap, std::pair<int,int> &nbinswall, std::pair<int,int> &nbinstopcap, std::pair<int,int> &nbinsbottomcap);
#include "makepmtmaps_standalone.cxx"     // definition of this function
// needed for drawing tank 2D map histograms
std::map<int, std::pair<int,int> > topcappositionmap;
std::map<int, std::pair<int,int> > bottomcappositionmap;
std::map<int, std::pair<int,int> > wallpositionmap;
std::pair<int,int> nbinswall;
std::pair<int,int> nbinstopcap;
std::pair<int,int> nbinsbottomcap;

const double SPEED_OF_LIGHT=29.9792458; // cm/ns - match time and position units of PMT
const double REF_INDEX_WATER=1.42535;   // maximum Rindex gives maximum travel time

int print_event_summaries(TString filepathin=""){
  
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(1);
  ColourWheel thecolourwheel = ColourWheel();
  
  std::string filepath="";
  if(filepathin==""){
    filepath = "/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_0.root";
    //filepath = "/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_0.isomus.root";
    //filepath = "/annie/app/users/moflaher/wcsim/build/wcsim_0.root";
    //filepath = "/annie/app/users/moflaher/wcsim/build/wcsim_0.nosmearisoes.root";
    //filepath = "/annie/app/users/moflaher/wcsim/build/wcsim_0.nosmearphotonbomb.root";
    //filepath = "/annie/app/users/moflaher/wcsim/build/wcsim_0.beames.root";
    //filepath = "/annie/app/users/moflaher/wcsim/build/wcsim_0.isomus.root";
    //filepath = "/annie/app/users/moflaher/wcsim/build/wcsim_0.photonbomb.root";
    //filepath = "/pnfs/annie/persistent/users/moflaher/wcsim/multipmt/tankonly/wcsim_3-12-18_ANNIEp2v6_BNB_Water_10k_22-05-17/wcsim_0.0.0.root";
    //filepath = "/pnfs/annie/persistent/users/moflaher/wcsim/multipmt/tankonly/wcsim_3-12-18_ANNIEp2v6_BNB_Water_10k_22-05-17/wcsim_0.1006.8.root";
    //filepath = "/pnfs/annie/persistent/users/moflaher/wcsim/multipmt/tankonly/wcsim_3-12-18_ANNIEp2v6_BNB_Water_10k_22-05-17/wcsim_0.0.3.root";
    //filepath = "/pnfs/annie/persistent/users/moflaher/wcsim/multipmt/tankonly/wcsim_3-12-18_ANNIEp2v6_BNB_Water_10k_22-05-17/wcsim_0.*.root";
    //filepath="/pnfs/annie/persistent/users/moflaher/wcsim/lappd/tankonly/wcsim_lappd_tankonly_03-05-17_BNB_World_10k_29-06-17/wcsim_0.4000.root";
    //filepath="/pnfs/annie/persistent/users/moflaher/wcsim/lappd/tankonly/wcsim_lappd_tankonly_24-09-17_BNB_Water_10k_22-05-17/wcsim_0.0.0.root";
    //filepath="/pnfs/annie/persistent/users/moflaher/wcsim/multipmt/tankonly/wcsim_multipmt_tankonly_17-06-17_rhatcher/00362859-23ba-4f54-8b36-bb5b0a11d841-wcsim_0.2347.root";
    //filepath="/annie/app/users/moflaher/wcsim/mutestresultsanniep2v6/wcsim_0.root";
  }
  
  //////////////////////////////
  // Get Pointers from Trees
  //////////////////////////////
  
  //TFile* f=TFile::Open(filepath.c_str());
  //TTree* t = (TTree*)f->Get("wcsimT");
  TChain* t = new TChain("wcsimT");
  t->Add(filepath.c_str());
  WCSimRootEvent* e=new WCSimRootEvent();
  TBranch* b=0;
  t->SetBranchAddress("wcsimrootevent", &e, &b);
  //WCSimRootTrigger* r;
  //WCSimRootEventHeader* h;
  
  WCSimRootGeom* geo=nullptr;
  t->GetEntry(0);
  TFile* f = t->GetCurrentFile();
  TTree* gtree = (TTree*)f->Get("wcsimGeoT");
  gtree->SetBranchAddress("wcsimrootgeom",&geo);
  gtree->GetEntry(0);
  MakePMTmap(geo, topcappositionmap, bottomcappositionmap, wallpositionmap, nbinswall, nbinstopcap, nbinsbottomcap);
  
  //////////////////////////////
  // Verbosity Limiters
  //////////////////////////////
  
  int maxentriestoprint=100000;
  int maxtriggerstoprint=1;
  int maxprimariestoprint=1; // should be always 1.
  int maxtrackstoprint=0;
  int maxdigitstoprint=0;
  int maxphotonsperdigittoprint=0;
  int maxphotonstoprint=0;
  
  //////////////////////////////
  // Histogram Declarations
  //////////////////////////////
  
  bool absolute_digit_times=true; // whether to plot digit times with trigger time subtracted
  int digit_historic_offset=0;  // subtract from digit times
  TH1D* digittimes = new TH1D("digittimes","Digit Times",200,-20,50);
  TH1D* digittimes2 = new TH1D("digittimes2","Digit Times2",200,-5,100000); // extended trigger digits
  TH1D* photontimes = new TH1D("photontimes","Digit Times",100,-5,50); //-10,10
  TH1D* hcharge = new TH1D("hcharge","Charge",100,0,30);
  TH2D* hqvst = new TH2D("hqvst","charge vs digit time",60,-10,50,100,0,300);
  TH1D* htrigt = new TH1D("htrigt","delayed trigger times",200,0,100000);
  TH1D* hmuontimes = new TH1D("hmuontimes","primary muon times",200,0,100);
  TH1D* hphotindigits = new TH1D("hphotindigits","Times of photon hits in digit",500,-5,30);
  TH1D* hphotrelfirstphot = new TH1D("hphotrelfirstphot","Photon Times Rel to First Photon in Digit",300,-5,30);
  TH1D* hdigitrelfirstphot = new TH1D("hdigitrelfirstphot","Digit Times Rel to First Photon in Digit",300,-15,35);
  TH2D* htvstdiff = new TH2D("htvstdiff","digit time vs diff",100,-5,50,100,-20,20);
  TH2D* hqvstdiff = new TH2D("hqvstdiff","digit Q vs t diff",100,-5,50,100,-20,20);
  TH2D* hnphotvstdiff = new TH2D("hnphotvstdiff","Num Phots vs digit delay",100,-15,35,50,0,50);
  TH2D* hddelayvstdiff = new TH2D("hddelayvstdiff","digit delay vs t diff",100,-10,10,100,-5,10);
  TH1D* htimeresid = new TH1D("htimeresid","Digit Time residual",300,-500,2500);
  TH2D* htimeresidvsnsteps = new TH2D("htimeresidvsnsteps","time residual vs # Steps",100,-5,30,12,0,22);
  TH2D* htimeresidvsq = new TH2D("htimeresidvsq","time residual vs charge",100,-10,10,100,-1,10);
  TH1D* hnumsteps = new TH1D("hnumsteps","Number of Tracking Steps",100,-1,50);
  TH1D* hsmearing = new TH1D("hsmearing","Time Smearing",300,-30,30);
//  std::vector<TPointSet3D> timeresidsets;  // how to set the point colours?
//  TGraph2D axishack = TGraph2D();          // can use z coordinate as colour with colz, but 2D only
  gStyle->SetMarkerStyle(7); // critical, we won't be able to change the sizes/styles after the fact. 
  gStyle->SetPalette(84);    // colour gradient, something more sensible than rainbow - viridis
  // must call this BEFORE making tree
  TTree* pointsetree = new TTree("pointsettree","Time Residuals"); // easiest (only?) way to make coloured 3D plot
  double pointsetx, pointsety, pointsetz, pointsett;
  pointsetree->Branch("pointsetx",&pointsetx);
  pointsetree->Branch("pointsety",&pointsety);
  pointsetree->Branch("pointsetz",&pointsetz);
  pointsetree->Branch("pointsett",&pointsett);
  std::map<std::string,TH1D*> timediffvspmttype;
  std::map<int,TH1D*> timediffvsphotnum;
  std::map<int,TH1D*> timediffvsdigitnum;
  
  // Wall maps
  TH2D* histowall = new TH2D("map_wall", "Wall Distribution",nbinswall.first+2,-1,nbinswall.first+1,nbinswall.second+2,-1,nbinswall.second+1);
  TH2D* histotop = new TH2D("map_top","Top Cap Distribution",nbinstopcap.first+2,-1,nbinstopcap.first+1,nbinstopcap.second+2,-1,nbinstopcap.second+1);
  TH2D* histobottom = new TH2D("map_bottom","Bottom Cap Distribution",nbinsbottomcap.first+2,-1,nbinsbottomcap.first+1,nbinsbottomcap.second+2,-1,nbinsbottomcap.second+1);
  std::map<std::string, TH2D*> maphistos;
  maphistos.emplace("histowall",histowall);
  maphistos.emplace("histotop",histotop);
  maphistos.emplace("histobottom",histobottom);
  
  //////////////////////////////
  // BEGIN MAIN LOOP
  //////////////////////////////
  
  auto numentries = t->GetEntries();
  cout<<"This run has "<<numentries<<" entries"<<endl;
  for(int i=0; i<min(maxentriestoprint,static_cast<int>(numentries)); i++){
    cout<<"EVENT "<<i<<endl<<"======="<<endl;
    t->GetEntry(i);
    cout<<"This entry had "<<e->GetNumberOfEvents()<<" triggers"<<endl;
    cout<<"This entry had "<< ( (e->HasSubEvents()) ? "1 or more" : "no" ) <<" delayed triggers"<<endl;
    cout<<"This entry had "<<e->GetNumberOfSubEvents()<<" delayed triggers"<<endl;
    // XXX XXX XXX cherenkov hits for ALL triggers are stored in trigger 0! XXX XXX XXX 
    WCSimRootTrigger* firsttrig=e->GetTrigger(0);
    for(int j=0; j<min(maxtriggerstoprint,static_cast<int>(e->GetNumberOfEvents())); j++){
      cout<<"  TRIGGER "<<j<<endl;
      WCSimRootTrigger* r=e->GetTrigger(j);
      WCSimRootEventHeader* h=r->GetHeader();
      if(h->GetDate()>0) break;  // if no light in prompt window but a delayed Ndigits trigger....
      cout<<"  >>> Trigger time was : "<<h->GetDate()<<endl; // TRIGGER TIME OF 0 MEANS NO NDIGITS TRIGGER
      // only for old file versions: prompt trigger means time is 0
      if(j>0){ htrigt->Fill(h->GetDate()); }
      
      //////////////////////////////
      // Primary Vertices Loop
      //////////////////////////////
      // primary vertices are ONLY valid for first trigger. Should be always 1
      Float_t True_Vertex_X, True_Vertex_Y, True_Vertex_Z;
      if(j==0){
        int nprimaries=1;
        nprimaries = r->GetNvtxs();
        cout<<"  This trigger had "<<nprimaries<<" vertex(es)"<<endl;
        for(int primaryi=0; primaryi<min(maxprimariestoprint,max(1,nprimaries)); primaryi++){
          True_Vertex_X=r->GetVtxs(primaryi,0);
          True_Vertex_Y=r->GetVtxs(primaryi,1);
          True_Vertex_Z=r->GetVtxs(primaryi,2);
          cout<<"  Primary vertex "<<primaryi<<" was at ("
              <<True_Vertex_X<<", "<<True_Vertex_Y<<", "<<True_Vertex_Z<<")"<<endl;
        }
      }
      
      //////////////////////////////
      // Track Loop
      //////////////////////////////
      int numtracks= r->GetNtrack();
      int numphotons = firsttrig->GetCherenkovHits()->GetEntries();
      int ndigits=r->GetCherenkovDigiHits()->GetEntries();
      cout<<"   event "<<i<<" trigger "<<j<<" had:"<<endl
          <<"   Num Tracks: "<<numtracks<<endl
          <<"   Num Photons: " <<numphotons<<endl
          <<"   Num Digits: "<<ndigits<<endl;
      // scan through the truth tracks, find the primary muon and pull vertex info from it
      int primarytrack=0;
      int primarycount=0;
      for(int track=0; track<min(maxtrackstoprint,numtracks); track++){
        WCSimRootTrack* nextrack = (WCSimRootTrack*)r->GetTracks()->At(track);
        //Int_t     GetIpnu()             pdg
        //Int_t     GetFlag()             -1: neutrino primary, -2: neutrino target, 0: other
        //Float_t   GetM()                mass
        //Float_t   GetP()                momentum magnitude
        //Float_t   GetE()                energy (inc rest mass^2)
        //Float_t   GetEndE()             energy on stopping of particle tracking
        //Float_t   GetEndP()             momentum on stopping of particle tracking
        //Int_t     GetStartvol()         starting volume: 10 is tank, 20 is facc, 30 is mrd
        //Int_t     GetStopvol()          stopping volume: but these may not be set.
        //Float_t   GetDir(Int_t i=0)     momentum unit vector
        //Float_t   GetPdir(Int_t i=0)    momentum vector
        //Float_t   GetPdirEnd(Int_t i=0) direction vector on stop tracking
        //Float_t   GetStop(Int_t i=0)    stopping vertex x,y,z for i=0-2, in cm
        //Float_t   GetStart(Int_t i=0)   starting vertex x,y,z for i=0-2, in cm
        //Int_t     GetParenttype()       parent pdg, 0 for primary.
        //Float_t   GetTime()             trj->GetGlobalTime(); starting time of particle
        //Float_t   GetStopTime()
        //Int_t     GetId()               wcsim trackid
        
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
        cout<<"    Track "<<track<<"{"
            <<  " Flag: "<<nextrack->GetFlag()
            <<" | PDG: "<<nextrack->GetIpnu()
            <<" | ID: "<<nextrack->GetId()
            <<" | ParentPDG: "<<nextrack->GetParenttype()
            <<" | Vtx: ("<<nextrack->GetTime()<<", "<<Primary_Vertex_X<<", "<<Primary_Vertex_Y<<", "<<Primary_Vertex_Z<<")"
#if FILE_VER>5
            <<" | sProc: "<<nextrack->GetStartProcess()
            <<" | eProc: "<<nextrack->GetEndProcess()
#endif
            <<" }"<<endl;
//        if(nextrack->GetParenttype()==0&&nextrack->GetFlag()!=-1){
//          primarycount++;
//          cout<<"event "<<i<<" trigger "<<j<<" primary "<<primarytrack<<", PDG:"<<nextrack->GetIpnu()
//              <<", trackID "<<nextrack->GetId()<<", Flag:"<<nextrack->GetFlag()<<", at ("
//              <<Primary_Vertex_X<<", "<<Primary_Vertex_Y<<", "<<Primary_Vertex_Z<<")"<<endl;
//          ++primarytrack;
//        }  // conditional on tracks representing primary particles
      }  // loop over tracks
      
//        /////////////// muon time histogram
//        int primarymuonindex=-1;
//        for(int tracki=0; tracki<numtracks; tracki++){
//          WCSimRootTrack* nextrack = (WCSimRootTrack*)r->GetTracks()->At(tracki);
//          if(nextrack->GetFlag()==0&&nextrack->GetParenttype()==0&&nextrack->GetIpnu()==13){  // primary muon
//            primarymuonindex=tracki;
//            break;
//          }    // conditional on being a primary muon
//        }
//        if(primarymuonindex<0){ continue; } // there was no primary muon in this event, skip it.
//        WCSimRootTrack* mutrack = (WCSimRootTrack*)r->GetTracks()->At(primarymuonindex);
//        cout<<"PRIMARY MUON AT "<<mutrack->GetTime()<<" TO "<<mutrack->GetStopTime()<<endl;
//        hmuontimes->Fill(mutrack->GetTime());
//        /////////////// end histogram
      
      //////////////////////////////
      // DIGIT LOOP
      //////////////////////////////
      double lastdigittime=0;
      // loop over digits
      cout<<"    loop over digits:"<<endl;
      for(int digiti=0; digiti<min(maxdigitstoprint,(int)ndigits); digiti++){
        WCSimRootCherenkovDigiHit* digihit=(WCSimRootCherenkovDigiHit*)(r->GetCherenkovDigiHits()->At(digiti));
        std::vector<int> truephotonindices = digihit->GetPhotonIds();
        cout<<"      digit "<<digiti<<" at time "<<digihit->GetT()<<"ns has "<<truephotonindices.size()<<" true photons"<<endl;
        for(int photoni=0; photoni<min(maxphotonsperdigittoprint,(int)truephotonindices.size()); photoni++){
          int thephotonsid = truephotonindices.at(photoni);
          WCSimRootCherenkovHitTime *thehittimeobject = 
            (WCSimRootCherenkovHitTime*)firsttrig->GetCherenkovHitTimes()->At(thephotonsid);
          Int_t thephotonsparenttrackid = thehittimeobject->GetParentID();
          //cout<<"    digit "<<digiti<<", photon "<<photoni<<" has parent "<<thephotonsparenttrackid<<endl;
         cout<<"        digit "<<digiti<<", photon "<<photoni<<" has truetime "<<thehittimeobject->GetTruetime()<<endl;
        }
      }
      
      //////////////////////////////
      // DIGIT HISTOGRAM LOOP
      //////////////////////////////
      // first loop to find the earliest digit on this PMT
      int earliestdigit=-1;
      float earliestdigittime = std::numeric_limits<float>::max();
      for(int digiti=0; digiti<ndigits; digiti++){
        WCSimRootCherenkovDigiHit* digihit=(WCSimRootCherenkovDigiHit*)(r->GetCherenkovDigiHits()->At(digiti));
        float adigittime = digihit->GetT() - digit_historic_offset + h->GetDate();
        if(adigittime<earliestdigittime){ earliestdigittime=adigittime; earliestdigit=digiti; }
      }
      
      for(int digiti=0; digiti<ndigits; digiti++){
      //for(int digiti=earliestdigit; (earliestdigit>0&&digiti<(earliestdigit+1)); digiti++){
        //if(j==0) break;
        //if(digiti>4) break;
        WCSimRootCherenkovDigiHit* digihit=(WCSimRootCherenkovDigiHit*)(r->GetCherenkovDigiHits()->At(digiti));
        float adigittime;
        if(absolute_digit_times){ adigittime = digihit->GetT() - digit_historic_offset + h->GetDate(); }
        else{ adigittime = digihit->GetT() - digit_historic_offset; }
//          if(digihit->GetT()==0.f){
//            cout<<"digit "<<digiti<<" at "<<digihit->GetT()<<"/"<<(digihit->GetT() + h->GetDate())
//                <<" rel/abs, is digit "<<" in trigger "<<j<<"/"<<e->GetNumberOfEvents()
//                <<" has charge "<<digihit->GetQ()<<" from "<<digihit->GetPhotonIds().size()<<" photons,"
//                <<" and occurs on tube "<<digihit->GetTubeId()<<endl;
//          }
        //if(not (abs(adigittime)<0.5f&&digihit->GetQ()>0.&&digihit->GetQ()<2.5))
        digittimes->Fill(adigittime);
        digittimes2->Fill(adigittime);
        //if(not (abs(adigittime)<0.5f&&digihit->GetQ()>0.&&digihit->GetQ()<2.5))
        hqvst->Fill(adigittime,digihit->GetQ());
        hcharge->Fill(digihit->GetQ());
        
        if((j>0) && (std::floor(adigittime/8.)>digiti)){
          std::cerr<<"#### DIGIT "<<digiti<<" AT TIME "<<adigittime
                   <<", WINDOW "<<std::floor(adigittime/8.)<<std::endl;
        }
        
        /////////////////// (hit time - digit time) histo, split for pmt type
        WCSimRootPMT thepmt = geo->GetPMT(digihit->GetTubeId());
        std::string pmt_type = thepmt.GetName();
        TH1D* tdiffforthistype = nullptr;
        if(timediffvspmttype.count(pmt_type)!=0){
          tdiffforthistype = timediffvspmttype.at(pmt_type);
        } else {
          std::string hname = pmt_type + "_tdiff";
          tdiffforthistype =  new TH1D(hname.c_str(),hname.c_str(),200, -10,10);
          timediffvspmttype.emplace(pmt_type,tdiffforthistype);
        }
        /////////////////////
        
        // distance to pmt
        double distx = thepmt.GetPosition(0) /*- geo->GetWCOffset(0)*/ - True_Vertex_X; // includes offset
        double disty = thepmt.GetPosition(1) /*- geo->GetWCOffset(1)*/ - True_Vertex_Y;
        double distz = thepmt.GetPosition(2) /*- geo->GetWCOffset(2)*/ - True_Vertex_Z;
        double distancefromvertex = sqrt(pow(distx,2.)+pow(disty,2.)+pow(distz,2.));
        double htres = adigittime-((distancefromvertex*REF_INDEX_WATER)/SPEED_OF_LIGHT);
        htimeresid->Fill(htres);
        
        auto thephotonids = digihit->GetPhotonIds();
        int ndigitphots = thephotonids.size();
        
        // quick loop to find and note the time of the earliest photon in the digit
        double firstphotontime=999999;
//        std::cout<<"digit on tube "<<digihit->GetTubeId()<<" at time "<<adigittime<<" with photons {";
        for(auto aphotonid : thephotonids){
          WCSimRootCherenkovHitTime *thehittimeobject = 
            (WCSimRootCherenkovHitTime*)firsttrig->GetCherenkovHitTimes()->At(aphotonid);
          Int_t thephotonsparenttrackid = thehittimeobject->GetParentID();
          if(thephotonsparenttrackid==-1){ continue; } // NOISE HIT
          double ahittime = thehittimeobject->GetTruetime();
          //double ahittime = thehittimeobject->GetSmeartime();
          hsmearing->Fill(thehittimeobject->GetSmeartime() - thehittimeobject->GetTruetime());
          if(ahittime<firstphotontime) firstphotontime=ahittime;
//          std::cout<<ahittime<<", ";
        }
//        std::cout<<"}"<<std::endl;
        hdigitrelfirstphot->Fill(adigittime-firstphotontime);
        hnphotvstdiff->Fill(adigittime-firstphotontime,ndigitphots);
        
        // fill photon times only for photons that were in a digit
        //if(adigittime==0){
          int loopi=0;
          bool good_digit=true;
          for(auto aphotonid : thephotonids){
            //if(thephotonids.size()>1) break;
            //if(loopi>0) break;        // only the first photon in the digit
            //if(loopi==0){ continue; } // only delayed photons
            WCSimRootCherenkovHitTime *thehittimeobject = 
              (WCSimRootCherenkovHitTime*)firsttrig->GetCherenkovHitTimes()->At(aphotonid);
            Int_t thephotonsparenttrackid = thehittimeobject->GetParentID();
            if(thephotonsparenttrackid==-1){ good_digit=false; continue; } // NOISE HIT
            double ahittime = thehittimeobject->GetTruetime();  // XXX XXX XXX
            //double ahittime = thehittimeobject->GetSmeartime();
            hddelayvstdiff->Fill(adigittime-firstphotontime,ahittime-adigittime,1./ndigitphots);
            tdiffforthistype->Fill(adigittime-ahittime);
            //hphotindigits->Fill(ahittime-adigittime);
            //htvstdiff->Fill(adigittime,ahittime-adigittime);
            hqvstdiff->Fill(digihit->GetQ(),ahittime-adigittime);
//            auto thescatters = thehittimeobject->GetScatterings();
//            for(auto&& aprocess : thescatters){
//              std::cerr<<aprocess.first<<":"<<aprocess.second<<std::endl;
//            }
            if( (thehittimeobject->GetScatterings().size()==2) &&
                (thehittimeobject->GetScatterings().begin()->first=="FresnelReflection_WCBarrel_WCPMT") &&
                (thehittimeobject->GetScatterings().begin()->second==1) &&
                (thehittimeobject->GetScatterings().rbegin()->first.find("FresnelReflection_WCPMT_ANNIEp2v6-glassFaceWCPMT")!=std::string::npos) &&
                (thehittimeobject->GetScatterings().rbegin()->second==1)
              ){ // nothing but transportation, detection seems to involve reflection
              // calculate time residual
              hphotindigits->Fill(ahittime-adigittime);
              hphotrelfirstphot->Fill(ahittime-firstphotontime);
              htvstdiff->Fill(adigittime,ahittime-adigittime);
              //hqvstdiff->Fill(digihit->GetQ(),ahittime-adigittime);
              //hddelayvstdiff->Fill(adigittime-firstphotontime,ahittime-adigittime,1./ndigitphots);
              float trackendx, trackendy, trackendz;
              bool trackfound=false;
              for(int track=0; track<numtracks; track++){
                WCSimRootTrack* nextrack = (WCSimRootTrack*)r->GetTracks()->At(track);
                if(thephotonsparenttrackid==nextrack->GetId()){
                  trackfound=true;
                  trackendx=nextrack->GetStop(0);
                  trackendy=nextrack->GetStop(1);
                  trackendz=nextrack->GetStop(2);
//                  if((abs(nextrack->GetStart(0)- geo->GetWCOffset(0))>0.01) ||
//                     (abs(nextrack->GetStart(1)- geo->GetWCOffset(1))>0.01) ||
//                     (abs(nextrack->GetStart(2)- geo->GetWCOffset(2))>0.01)){
//                      std::cerr<<"TRACK WITH NON-ORIGIN START AT ("
//                               <<(nextrack->GetStart(0)- geo->GetWCOffset(0))<<", "
//                               <<(nextrack->GetStart(1)- geo->GetWCOffset(1))<<", "
//                               <<(nextrack->GetStart(2)- geo->GetWCOffset(2))<<")"<<std::endl;
//                  }
                  double distxphot = trackendx - True_Vertex_X; // includes offset
                  double distyphot = trackendy - True_Vertex_Y;
                  double distzphot = trackendz - True_Vertex_Z;
                  double distancefromvertexphot = sqrt(pow(distxphot,2.)+pow(distyphot,2.)+pow(distzphot,2.));
                  double traveltime= (distancefromvertexphot*REF_INDEX_WATER)/SPEED_OF_LIGHT;
                  double htresphot = ahittime-nextrack->GetTime()-traveltime;
                  //htimeresid->Fill(htresphot);
                  int nsteps = thehittimeobject->GetNumScatterings();
                  hnumsteps->Fill(nsteps);
                  if(not (nsteps==2&&htresphot<1)) htimeresidvsnsteps->Fill(htresphot,nsteps);
                  //if(htresphot<-5||htresphot>30){
                  //  std::cerr<<"VERY BAD HIT TIME RES: "<<htresphot<<std::endl;
                  //} else {
                    pointsetx=trackendx;
                    pointsety=trackendy;
                    pointsetz=trackendz;
                    pointsett=htresphot;
                    pointsetree->Fill();
                  //}
                  if(htresphot>5){
                    std::cerr<<"Bad hit time residual! Photon begins at ("
                               <<(nextrack->GetStart(0)- geo->GetWCOffset(0))<<", "
                               <<(nextrack->GetStart(1)- geo->GetWCOffset(1))<<", "
                               <<(nextrack->GetStart(2)- geo->GetWCOffset(2))<<")"
                               <<" and ends at ("
                               <<(nextrack->GetStop(0)- geo->GetWCOffset(0))<<", "<<", "
                               <<(nextrack->GetStop(1)- geo->GetWCOffset(1))<<", "<<", "
                               <<(nextrack->GetStop(2)- geo->GetWCOffset(2))<<", "<<")"
                               <<" on PMT "<<digihit->GetTubeId()<<" at ("
                               <<(thepmt.GetPosition(0) - geo->GetWCOffset(0))<<", "
                               <<(thepmt.GetPosition(1) - geo->GetWCOffset(1))<<", "
                               <<(thepmt.GetPosition(2) - geo->GetWCOffset(2))<<")"
                               <<", expected time is "<<traveltime+h->GetDate()
                               <<", hit time is "<<ahittime<<"ns"
                               <<", digit time is "<<digihit->GetT()
                               <<", trigger time is "<<h->GetDate()
                               <<", parent is "<<thephotonsparenttrackid
                               <<" which starts at "<<nextrack->GetTime()
                               <<std::endl;
                  }
                  break;
                }
              }
//              if(trackfound==false){
//                std::cerr<<"COULD NOT FIND TRACK FOR PHOTON "<<thehittimeobject->GetParentID()<<"!!!"<<std::endl;
//              }
            } else {
              good_digit=false;
            }
//            if(thehittimeobject->GetScatterings().size()>0){
//              std::cout<<"photon with non-transportation process(es): ";
//              std::map<std::string,int> procmap = thehittimeobject->GetScatterings();
//              for(auto&& aproc : procmap) std::cout<<(aproc.first)<<":"<<aproc.second<<", ";
//              std::cout<<endl;
//            } else {
//              std::cout<<"photon with no processes!"<<std::endl;
//            }
            
            TH1D* tdiffforthisphoti = nullptr;
            if(timediffvsphotnum.count(loopi)!=0){
              tdiffforthisphoti = timediffvsphotnum.at(loopi);
            } else {
              std::string hname = "tdiff_phot_"+std::to_string(loopi);
              tdiffforthisphoti =  new TH1D(hname.c_str(),hname.c_str(),200, -10,10);
              timediffvsphotnum.emplace(loopi,tdiffforthisphoti);
            }
            tdiffforthisphoti->Fill(ahittime-adigittime);
            
            TH1D* tdiffforthisdigiti = nullptr;
            if(timediffvsdigitnum.count(digiti)!=0){
              tdiffforthisdigiti = timediffvsdigitnum.at(digiti);
            } else {
              std::string hname = "tdiff_digit_"+std::to_string(digiti);
              tdiffforthisdigiti =  new TH1D(hname.c_str(),hname.c_str(),200, -10,10);
              timediffvsdigitnum.emplace(digiti,tdiffforthisdigiti);
            }
            tdiffforthisdigiti->Fill(ahittime-adigittime);
            
//            if((ahittime-adigittime)==0.f){
//              cerr<<"photon "<<loopi<<" in digit "<<digiti<<" has time "<<ahittime<<", equal to digit time at "
//                  <<digihit->GetT()<<"/"<<(digihit->GetT() + h->GetDate())
//                  <<" rel/abs, is in trigger "<<j<<"/"<<e->GetNumberOfEvents()
//                  <<" has charge "<<digihit->GetQ()<<" from "<<digihit->GetPhotonIds().size()<<" photons,"
//                  <<" and occurs on tube "<<digihit->GetTubeId()<<endl;
//            }
            
            ++loopi;
          }
          
          if(good_digit){
            //htimeresid->Fill(htres);
            htimeresidvsq->Fill(htres,digihit->GetQ());
          }
        //}
      }
      /////////////// end digit histogram
      
      //////////////////////////////
      // PHOTON LOOP
      //////////////////////////////
      // loop over true hits
      cout<<"    loop over photons:"<<endl;
      for(int photoni=0; photoni<min(maxphotonstoprint,(int)numphotons); photoni++){
        WCSimRootCherenkovHitTime *thehittimeobject = 
          (WCSimRootCherenkovHitTime*)firsttrig->GetCherenkovHitTimes()->At(photoni);
        Int_t thephotonsparenttrackid = thehittimeobject->GetParentID();
        cout<<"      photon "<<photoni<<" has parent "<<thephotonsparenttrackid<<endl;
      }
      
//        /////////////// photon histogram
        for(int photoni=0; photoni<numphotons; photoni++){
          WCSimRootCherenkovHitTime *thehittimeobject = 
            (WCSimRootCherenkovHitTime*)firsttrig->GetCherenkovHitTimes()->At(photoni);
//          hphotindigits->Fill(thehittimeobject->GetTruetime());
            auto thescatters = thehittimeobject->GetScatterings();
            for(auto&& aprocess : thescatters){
              std::cerr<<aprocess.first<<":"<<aprocess.second<<std::endl;
            }
        }
//        /////////////// end histogram
      
    } // loop over subevents
    e->ReInitialize();
  } // loop over events
  
  /////////////////////////
  // CONFIGURE HISTOGRAMS
  /////////////////////////
  THStack* tdiffbyphotstack = new THStack("tdiffbyphot","tdiffbyphot");
  for(auto histpair : timediffvsphotnum){
    TH1D* ahist = histpair.second;
    ahist->SetLineColor(thecolourwheel.GetNextColour());
    tdiffbyphotstack->Add(ahist);
  }
  THStack* tdiffbydigitstack = new THStack("tdiffbydigitstack","tdiffbydigitstack");
  for(auto histpair : timediffvsdigitnum){
    TH1D* ahist = histpair.second;
    ahist->SetLineColor(thecolourwheel.GetNextColour());
    tdiffbydigitstack->Add(ahist);
  }
  
  /////////////////////////
  // DRAW HISTOGRAMS
  /////////////////////////
  TCanvas* c1 = new TCanvas("c1");
  //digittimes->Draw();
  //htvstdiff->Draw("colz");
  //hqvstdiff->Draw("colz");
  hnphotvstdiff->Draw("colz");
  //hphotindigits->Draw();
  //hdigitrelfirstphot->Draw();
  //digittimes2->Draw();
  //hcharge->Draw();
  //hqvst->Draw("colz");
  //hmuontimes->Draw();
  //photontimes->Draw();
  //htimeresid->Draw();
  //hnumsteps->Draw();
  c1->Modified();
  c1->Update();
  
  TCanvas* c2 = new TCanvas("c2");
  //htrigt->Draw();
  //photontimes->Draw();
  //hphotindigits->Draw();
  //hphotrelfirstphot->Draw();
  //hdigitrelfirstphot->Draw();
  hsmearing->Draw();
  //hddelayvstdiff->Draw("colz");
  //hnphotvstdiff->Draw("colz");
  //htimeresidvsq->Draw("colz");
  //pointsetree->Draw("pointsetx:pointsety:pointsetz:pointsett","","colz");
  //htimeresidvsnsteps->Draw("colz");
  //hqvst->Draw("colz");
  //htvstdiff->Draw("colz");
  //tdiffbyphotstack->Draw();
  //tdiffbydigitstack->Draw();
//  int pmttypei=0;
//  for(auto histpair : timediffvspmttype){
//    TH1D* ahist = histpair.second;
//    ahist->SetLineColor(thecolourwheel.GetNextColour());
//    if(pmttypei==0) ahist->Draw();
//    else ahist->Draw("same");
//    ++pmttypei;
//  }
  c2->Modified();
  c2->Update();
  
//  TCanvas* cwall = new TCanvas("cwall");
//  histowall->Draw("colz");
//  TCanvas* ctopcap = new TCanvas("ctopcap");
//  histotop->Draw("colz");
//  TCanvas* cbottomcap = new TCanvas("cbottomcap");
//  histobottom->Draw("colz");
  
//  do{
//    gSystem->ProcessEvents();
//    std::this_thread::sleep_for(std::chrono::seconds(1));
//  } while (gROOT->FindObject("c1")!=nullptr);
  
  return 1;
}


void FillTankMapHist(WCSimRootGeom* geo, int tubeID, std::map<std::string, TH2D*> &maphistos, double weight=1){
  //Fill a bin on a 2D map of PMTs 
  WCSimRootPMT pmt = geo->GetPMT(tubeID);
  // WCSimRootPMT has members GetTubeNo(), GetCylLoc(), GetPosition(j), GetOrientation(j)
  // GetCylLoc(): 0=top cap, 2=bottom cap, 1=wall, 4=mrd, 5=veto, 3=obselete outer veto (shouldnt come up)
  // GetPosition(j), j=0..2: Returns x,y,z coordinates of the center of the sphere that forms the PMT.
  // GetOrientation(j), j=0..2: Returns the x,y,z components of a vector of the direction the PMT faces.
  // GetPMT(j) Returns a pmt object - NOT a pointer to a PMT object.
  //cout<<"Filling histogram for cylloc "<<pmt.GetCylLoc()<<" for tubeID "<<tubeID<<endl;
  
  double pmtx = pmt.GetPosition(0) - geo->GetWCOffset(0);
  double pmty = pmt.GetPosition(1) - geo->GetWCOffset(1);
  double pmtz = pmt.GetPosition(2) - geo->GetWCOffset(2);
  
  switch(pmt.GetCylLoc()){
    case 0: {
      if(topcappositionmap.count(tubeID)){
        std::pair<int,int> thebins = topcappositionmap.at(tubeID);
        TH2D* histotop=maphistos.at("histotop");
        if(histotop) histotop->Fill(thebins.first, thebins.second, weight);
      } else {cout<<"bad pmt: ID "<<tubeID<<" in CylLoc "<<pmt.GetCylLoc()<<endl;}
      break;
    }
    case 1: {
      if(wallpositionmap.count(tubeID)){
        std::pair<int,int> thebins = wallpositionmap.at(tubeID);
        TH2D* histowall=maphistos.at("histowall");
        if(histowall) histowall->Fill(thebins.first+0.5, thebins.second, weight);
      } else {cout<<"bad pmt: ID "<<tubeID<<" in CylLoc "<<pmt.GetCylLoc()<<endl;}
      break;
    }
    case 2: {
      if(bottomcappositionmap.count(tubeID)){
        std::pair<int,int> thebins = bottomcappositionmap.at(tubeID);
        TH2D* histobottom=maphistos.at("histobottom");
        if(histobottom) histobottom->Fill(thebins.first, thebins.second, weight);
      } else {cout<<"bad pmt: ID "<<tubeID<<" in CylLoc "<<pmt.GetCylLoc()<<endl;}
      break;
    }
    case 4: {
//        std::pair<int,int> thebins = mrdpositionmap.at(tubeID);
//        mrdhist->Fill(thebins.first, thebins.second, weight);
      break;
    }
    case 5: {
//        std::pair<int,int> thebins = faccpositionmap.at(tubeID);
//        facchist->Fill(thebins.first, thebins.second, weight);
      break;
    }
    default: {
      //cout<<"PMT "<<tubeID<<" has unknown location "<<pmt.GetCylLoc()<<"!"<<endl; 
      break;
    }
  }
}

void ClearMapHistos(std::map<std::string,TH2D*> maphistos){
  for(std::map<std::string,TH2D*>::iterator it= maphistos.begin(); it!=maphistos.end(); it++){
    it->second->Reset();
  }
}
