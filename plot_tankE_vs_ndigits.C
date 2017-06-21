/* vim:set noexpandtab tabstop=4 wrap */
/*TODO: to modify for new source files:
0. Re-generate and replace in this directory libWCSimRoot.so, .rootmap, .pcm files.
1. Re-enable #include RootOptions.hh
2. Disable timeArrayOffset lines.
*/
#ifndef VERBOSE
//#define VERBOSE
#endif
#include "TROOT.h"
#include "TSystem.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TApplication.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TLorentzVector.h"
#ifdef __MAKECINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif
#include "TLegend.h"
#include "TText.h"
#include "TColor.h"
#include <regex>
#include "TStyle.h"
#include <exception>	// for stdexcept
#include <vector>
#include <map>
#include <string>
#include <algorithm>	// remove and remove_if
#include <iostream>
#include <iomanip>
#include <fstream> 		//std::ofstream
#include <stdlib.h> 	/* atoi */
#include <valarray>
#ifdef __MAKECINT__
#pragma link C++ class std::map<std::string,bool>+;	// <<<< REQUIRED TO SAVE THIS KIND TO TREE
#endif
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
#include <time.h>         // clock_t, clock, CLOCKS_PER_SEC

// we need to #include all the WCSim headers.
#include "../wcsim/include/WCSimRootEvent.hh"
#include "../wcsim/include/WCSimRootGeom.hh"
#include "../wcsim/include/WCSimPmtInfo.hh"
#include "../wcsim/include/WCSimLAPPDInfo.hh"
#include "../wcsim/include/WCSimEnumerations.hh"
#include "../wcsim/include/WCSimRootLinkDef.hh"
#include "../wcsim/include/WCSimRootOptions.hh"

void MakePMTmap(WCSimRootGeom* geo, std::map<int, std::pair<int,int> > &topcappositionmap, std::map<int, std::pair<int,int> > &bottomcappositionmap, std::map<int, std::pair<int,int> > &wallpositionmap);
#include "makepmtmaps_standalone.cxx"     // definition of this function

void ColourPlotStyle();
void FillTankMapHist(WCSimRootGeom* geo, int tubeID, bool incone, std::map<std::string, TH2D*> &maphistos, double weight);
void ClearMapHistos(std::map<std::string,TH2D*> maphistos);	// clear the histograms

const Float_t MRD_width = (305./2.);      // half width of steel in cm
const Float_t MRD_height = (274./2.);     // half height of steel in cm
const Float_t MRD_layer2 = 290.755;       // position in wcsim coords of second scint layer in cm
const Float_t MRD_start = 325.5;          // position in wcsim coord of MRD front face in cm
const Float_t MRD_depth = 139.09;         // total depth of the MRD in cm
/* output from WCSim:
########## MRD front face: 325.5                     ##########
########## MRD total Z length: 139.09                ##########
########## MRD scintillator layer 0  (H) at z=336.08 ##########
########## MRD scintillator layer 1  (V) at z=348.19 ##########
########## MRD scintillator layer 2  (H) at z=360.30 ##########
########## MRD scintillator layer 3  (V) at z=372.41 ##########
########## MRD scintillator layer 4  (H) at z=384.52 ##########
########## MRD scintillator layer 5  (V) at z=396.63 ##########
########## MRD scintillator layer 6  (H) at z=408.74 ##########
########## MRD scintillator layer 7  (V) at z=420.85 ##########
########## MRD scintillator layer 8  (H) at z=432.96 ##########
########## MRD scintillator layer 9  (V) at z=445.07 ##########
########## MRD scintillator layer 10 (H) at z=457.18 ##########
*/
std::vector<double> mrdscintlayers{336.080, 348.190, 360.300, 372.410, 384.520, 396.630, 408.740, 420.850, 432.960, 445.070, 457.180 };

const Int_t numtankpmts=128+2*(26); // 26 pmts and lappds on each cap
const Int_t nummrdpmts=307;
const Int_t numvetopmts=26;
const Int_t caparraysize=8;         // pmts on the cap form an nxn grid where caparraysize=n
const Int_t pmtsperring=16;         // pmts around each ring of the main walls
const Int_t numpmtrings=8;          // num rings around the main walls

const Float_t tank_start = 15.70;          // front face of the tank in cm
const Float_t tank_radius = 152.4;         // tank radius in cm
const Float_t tank_halfheight = 198.;      // tank half height in cm
const Float_t tank_yoffset = -14.46;        // tank y offset in cm
/* from WCSimDetectorConfigs.cc
tankouterRadius= 1.524*m;
tankzoffset = 15.70*cm;
*/
const Float_t fidcutradius=tank_radius*0.8;			// fiducial volume is slightly lesss than the full tank radius
const Float_t fidcuty=50.;						// a meter total fiducial volume in the centre y
const Float_t fidcutz=0;	// fidcuial volume is before the tank centre.

// needed for drawing tank 2D map histograms
std::map<int, std::pair<int,int> > topcappositionmap;
std::map<int, std::pair<int,int> > bottomcappositionmap;
std::map<int, std::pair<int,int> > wallpositionmap;

const char* wcsimlibrarypath="/home/marc/LinuxSystemFiles/WCSim/gitver/wcsim/libWCSimRoot.so";
const char* outpath="/home/marc/LinuxSystemFiles/WCSim/gitver/root_work";
const char* inpath="/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in";

const Bool_t printneutrinoevent=false;

void plot_tankE_vs_ndigits(){
	ColourPlotStyle();
	
	// load WCSim library for reading WCSim files
	cout<<"loading "<<wcsimlibrarypath<<endl;
	gSystem->Load(wcsimlibrarypath);
	
	// ==============================================================================================
	// ==============================================================================================
	// Declare input paths, files, trees...
	TFile* wcsimfile=0;
	TString wcsimfilepath;
	TTree* wcsimT=0;
	Int_t numwcsimentries=0;
	
	// TChain for wcsim files
	TChain* c =  new TChain("wcsimT");
	TString chainpattern = TString::Format("%s/mu_tank_depos/*.root",inpath);
	cout<<"loading TChain entries from "<<chainpattern<<endl;
	c->Add(chainpattern);
	Int_t numents = c->GetEntries();
	cout<<"loaded "<<numents<<" entries in the chain"<<endl;
	if(numents<1){ return; }
	
	// wcsimT
	WCSimRootEvent* b=0, *m=0, *v=0;
	TBranch* bp=0, *mp=0, *vp=0;
	WCSimRootTrigger* atrigt=0, *atrigm=0, *atrigv=0;
	
	// geoT
	WCSimRootGeom* geo = 0; 
	
	// ==============================================================================================
	// ==============================================================================================
	// output file
	TFile* fileout = new TFile("TankEnergyVsCharge.root","RECREATE");
	fileout->cd();
	TTree* treeout = new TTree("treeout","Tank Event Properties");
	
	// store the file and event num info for further lookup later if needed
	double wcsimeventnum;
	TBranch* bWCSimEventNum = treeout->Branch("WCSimEventNum",&wcsimeventnum);
	// now information about the event
	TVector3 mustartvtx(0,0,0);
	TBranch* bMuonStartVtx = treeout->Branch("MuonStartVertex",&mustartvtx);
	TVector3 mustopvtx(0,0,0);
	TBranch* bMuonStopVtx = treeout->Branch("MuonStopVertex",&mustopvtx);
	double muonangle=0.;
	TBranch* bMuonAngle = treeout->Branch("MuonAngle",&muonangle);
	double mustartE=0.;
	TBranch* bMuonStartE = treeout->Branch("MuonStartEnergy",&mustartE);
	double mustartKE=0.;
	TBranch* bMuonStartKE = treeout->Branch("MuonStartKineticEnergy",&mustartKE);
	double mustopE=0.;
	TBranch* bMuonStopE = treeout->Branch("MuonStopEnergy",&mustopE);
	double mustopKE=0.;
	TBranch* bMuonStopKE = treeout->Branch("MuonStopKineticEnergy",&mustopKE);
	double mutracklengthintank=0.;
	TBranch* bMuonTrackLengthInTank = treeout->Branch("MuonTrackLengthInTank",&mutracklengthintank);
	int numTankDigits=0;
	TBranch* bNumTankDigits = treeout->Branch("TotalTankDigits",&numTankDigits);
	double totaltankcharge=0.;
	TBranch* bTotalTankCharge = treeout->Branch("TotalTankCharge",&totaltankcharge);
	int numtankdigitsfrommuon=0;
	TBranch* bNumTankDigitsFromMu = treeout->Branch("TankDigitsFromMuon",&numtankdigitsfrommuon);
	double tankchargefrommuon=0.;
	TBranch* bTankChargeFromMuon = treeout->Branch("TankChargeFromMuon",&tankchargefrommuon);
	std::vector<int> tanktubeshitbymu;
	TBranch* bTankTubesHitByMuon = treeout->Branch("TankTubesHitByMu",&tanktubeshitbymu);
	double fractionalchargeincone=0.;
	TBranch* bFractionOfMuonChargeInCone = 
		treeout->Branch("FractionOfMuonChargeInCone",&fractionalchargeincone);
	double upstreamcharge=0.;
	TBranch* bTotalUpstreamCharge = treeout->Branch("TotalUpstreamCharge",&upstreamcharge);
	double downstreamcharge=0.;
	TBranch* bTotalDownstreamCharge = treeout->Branch("TotalDownstreamCharge",&downstreamcharge);
	double topcapcharge=0.;
	TBranch* bTopCapCharge = treeout->Branch("TotalTopCapCharge",&topcapcharge);
	double bottomcapcharge=0.;
	TBranch* bBottomCapCharge = treeout->Branch("TotalBottomCapCharge",&bottomcapcharge);
	
	std::vector<TBranch*> thebranches{ bMuonStartVtx, bMuonStopVtx, bMuonStartE, bMuonTrackLengthInTank, bTankChargeFromMuon, bFractionOfMuonChargeInCone};
	int someit=0;
	bool haszombies=false;
	for(auto abranch : thebranches){
		if(abranch==0){ cout<<"branch "<<someit<<" is a zombie"<<endl; haszombies=true; }
		someit++;
	}
	assert(!haszombies&&"output file branches have zombies");
	
	// ===================================================================================================
	// ===================================================================================================
	// histograms file
	gROOT->cd();
	TFile* histofileout = new TFile("TankEvsQhistos.root","RECREATE");
	
	// separate one for wcsim, from digit integration
	TH1D* muedepositionswcsim = new TH1D("muedepositionswcsim","Distribution of Muon Energy Depositions In Tank ;Energy (PMT Q);Num Events",100,0.,1);
	TH1D* muedepositionsfidcut = new TH1D("muedepositionsfidcut","Distribution of Muon Energy Depositions In Tank (Fiducial);Energy (PMT Q);Num Events",100,0.,1);
	
	TH1D* fsltruetracklength = new TH1D("fsltruetracklength", "Distribution of True Track Lengths", 100, 0., 1500.);
	TH1D* fsltruetracklengthintank = new TH1D("fsltruetracklengthintank", "Distribution of True Track Lengths In Tank", 100, 0., 1500.);
	
	// test the hypothesis that track length in water can be estimated from total light in tank
	TH2D* tracklengthvsmuonlight = new TH2D("tracklengthvsmuonlight", "Muon Track Length vs Total Light from Muon", 100, 0., 1500., 100, 0., 50.);
	
	// Just to test the inside/outside cherenkov cone algorithm
	TH2D* inconehistowall = new TH2D("chargemap_incone_wall", "Charge Distribution Inside Cherenkov Cone (Wall)", pmtsperring+2,-1,pmtsperring+1,numpmtrings+2,-1,numpmtrings+1);
	TH2D* inconehistotop = new TH2D("chargemap_incone_top","Charge Distribution Inside Cherenkov Cone (Top Cap)",caparraysize+2,-1,caparraysize+1,caparraysize+2,-1,caparraysize+1);
	TH2D* inconehistobottom = new TH2D("chargemap_incone_bottom","Charge Distribution Inside Cherenkov Cone (Bottom Cap)",caparraysize+2,-1,caparraysize+1,caparraysize+2,-1,caparraysize+1);
	
	TH2D* outconehistowall = new TH2D("chargemap_outcone_wall", "Charge Distribution Outside Cherenkov Cone (Wall)", pmtsperring+2,-1,pmtsperring+1,numpmtrings+2,-1,numpmtrings+1);
	TH2D* outconehistotop = new TH2D("chargemap_outcone_top","Charge Distribution Outside Cherenkov Cone (Top Cap)",caparraysize+2,-1,caparraysize+1,caparraysize+2,-1,caparraysize+1);
	TH2D* outconehistobottom = new TH2D("chargemap_outcone_bottom","Charge Distribution Outside Cherenkov Cone (Bottom Cap)",caparraysize+2,-1,caparraysize+1,caparraysize+2,-1,caparraysize+1);
	
	// ===================================================================================================
	// ===================================================================================================
	// done declaring file: move to loading and processing
	gROOT->cd();
	cout<<"loading first wcsimT tree from "<<chainpattern<<" tchain"<<endl;
	c->LoadTree(0);
	Int_t treeNumber = -1;
	wcsimT = c->GetTree();
	Int_t thistreesentries = wcsimT->GetEntries();
	cout<<thistreesentries<<" entries in the first tree"<<endl;
	
	/*
	1. Load next g4dirt entry
	2. Check if genie primary, and volume is in tank - if not, continue
	3. If so, load genie entry.
	4. Check if interaction is QE - if not, continue
	5. If so, load wcsim detector response. 
	6. Load tracks, look for primary mu track through the MRD, record interaction details. 
	*/
	
	cout<<"looping over tchain entries"<<endl;
	//numents=10;
	for(Int_t inputEntry=0; inputEntry<numents; inputEntry++){
		/* 	1. Load next wcsimT entry */ 
		//==================================================================================================
		//==================================================================================================
#ifdef VERBOSE
		cout<<"loading entry "<<inputEntry<<endl;
#endif
		Long64_t localEntry = c->LoadTree(inputEntry);
		if( localEntry<0){ cout<<"end of tchain"<<endl; break; }
		Int_t nextTreeNumber = c->GetTreeNumber();
		if(treeNumber!=nextTreeNumber){
			cout<<"new tree: "<<nextTreeNumber<<endl;
			// this means we've switched file - need to load the new branch addresses (and potentially other trees)
			// first pull out the new file name
			wcsimT = c->GetTree();
			wcsimfile = wcsimT->GetCurrentFile();
			Int_t thistreesentries = wcsimT->GetEntries();
			cout<<"wcsimT has "<<thistreesentries<<" entries in this file"<<endl;
			
			// load the geometry tree and grab the geometry if we haven't already
			if(geo==0){
				TTree* geotree = (TTree*)wcsimfile->Get("wcsimGeoT");
				if(geotree==0){ cout<<"NO GEOMETRY IN FIRST FILE?"<<endl; assert(false); }
				geotree->SetBranchAddress("wcsimrootgeom", &geo);
				if (geotree->GetEntries() == 0) { cout<<"geotree has no entries!"<<endl; exit(9); }
				geotree->GetEntry(0);
				// pmtidtocopynum = makecopynummap(geo); // turns out this is a 1:1 mapping after all!
				//MakePMTmap(geo, topcappositionmap, bottomcappositionmap, wallpositionmap);
			}
			
			// wcsim trigger classes
			wcsimT->SetBranchAddress("wcsimrootevent",&b, &bp);
			wcsimT->SetBranchAddress("wcsimrootevent_mrd",&m, &mp);
			wcsimT->SetBranchAddress("wcsimrootevent_facc",&v, &vp);
			bp->SetAutoDelete(kTRUE);
			mp->SetAutoDelete(kTRUE);
			vp->SetAutoDelete(kTRUE);
			if(bp==0||mp==0||vp==0){ cout<<"branches are zombies!"<<endl; }
			
			treeNumber=nextTreeNumber;
		}
		
		//====================================================================================================
		//====================================================================================================
		// read only first subtrigger; delayed decay detector response is not interesting for primary FSL tracks
#ifdef VERBOSE
		cout<<"getting wcsim entry "<<inputEntry<<endl;
#endif
		wcsimT->GetEntry(inputEntry);
		wcsimeventnum=inputEntry;
		atrigt = b->GetTrigger(0);
		atrigm = m->GetTrigger(0);
		atrigv = v->GetTrigger(0);
		
		Int_t numtracks = atrigt->GetNtrack();
#ifdef VERBOSE
		cout<<"wcsim event had "<<numtracks<<" truth tracks"<<endl;
#endif
		
		for(int track=0; track<numtracks; track++){
			WCSimRootTrack* nextrack = (WCSimRootTrack*)atrigt->GetTracks()->At(track);
			cout<<"⇒";
			/* a WCSimRootTrack has methods: 
			Int_t     GetIpnu()             pdg
			Int_t     GetFlag()             -1: neutrino primary, -2: neutrino target, 0: other
			Float_t   GetM()                mass
			Float_t   GetP()                momentum magnitude
			Float_t   GetE()                energy (inc rest mass^2)
			Int_t     GetStartvol()         starting volume: 10 is tank, 20 is facc, 30 is mrd
			Int_t     GetStopvol()          stopping volume: but these may not be set.
			Float_t   GetDir(Int_t i=0)     momentum unit vector
			Float_t   GetPdir(Int_t i=0)    momentum vector
			Float_t   GetStop(Int_t i=0)    stopping vertex x,y,z for i=0-2, in cm
			Float_t   GetStart(Int_t i=0)   starting vertex x,y,z for i=0-2, in cm
			Int_t     GetParenttype()       parent pdg, 0 for primary.
			Float_t   GetTime()             trj->GetGlobalTime(); stopping(?) time of particle
			Int_t     GetId()               wcsim trackid
			*/
			
			/* 6. Load tracks, look for primary mu track through the MRD, record interaction details. */
			//===============================================================================================
			//===============================================================================================
			
			// ----------------------------------------------------------------------------------------------
			// Check if it's a primary muon starting in the tank
			// ----------------------------------------------------------------------------------------------
			// is it a (anti)muon?
			Int_t primarypdg = nextrack->GetIpnu();
#ifdef VERBOSE
			cout<<"primarypdg is "<<primarypdg<<endl;
#endif
			if(TMath::Abs(primarypdg)!=13) continue;       // not a muon
			
			// for now we use truth information
			// is it a primary?
			Int_t primaryparentpdg = nextrack->GetParenttype();
			if(primaryparentpdg!=0) continue;
			
			// does it start in the tank?
			Int_t primarystartvol = nextrack->GetStartvol();
//			// temporary override as these weren't correctly set in wcsim <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//			//TODO REMOVE THIS once 1M sample is regenerated
//			if(nextrack->GetStart(2)<tank_start){
//				primarystartvol = 20;						// start depth is facc or hall
//			} else if(nextrack->GetStart(2)>(tank_start+(2.*tank_radius))){
//				primarystartvol = 30;						// start depth is mrd or hall
//			} else {
//				primarystartvol = 10;						// start depth is tank
//			}
			
#ifdef VERBOSE
			cout<<"primarystartvol is "<<primarystartvol<<endl;
#endif
			if(primarystartvol!=10) continue;				// start volume is not the tank
			
			// does it stop in the mrd, or pass completely through the MRD? 
			Int_t primarystopvol = nextrack->GetStopvol();
//			// temporary override as these weren't correctly set in wcsim <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//			//TODO REMOVE THIS once 1M generated. 
//			// BUT need to think about how to check for 'range-out' mrd events. Maybe this is fine?
//			if(nextrack->GetStop(2)<tank_start){
//				primarystopvol = 20;						// start depth is facc or hall
//			} else if(nextrack->GetStop(2)>(tank_start+(2.*tank_radius))){
//				primarystopvol = 30;						// start depth is mrd or hall
//			} else {
//				primarystopvol = 10;						// start depth is tank
//			}
#ifdef VERBOSE
			cout<<"primarystopvol is "<<primarystopvol<<endl;
#endif
			
			TLorentzVector primarystartvertex(  nextrack->GetStart(0),
												nextrack->GetStart(1),
												nextrack->GetStart(2),
												0.);
			TLorentzVector primarystopvertex(   nextrack->GetStop(0),
												nextrack->GetStop(1),
												nextrack->GetStop(2),
												nextrack->GetTime());
			
			// for this analysis we want tracks that stopped by exiting the tank:
			double endpointx = primarystopvertex.X();
			double endpointy = primarystopvertex.Y()-tank_yoffset;
			double endpointz = primarystopvertex.Z()-tank_start-tank_radius;
			double endpointradius = sqrt(pow(endpointx,2) + pow(endpointz,2));
			// if radius isn't ~tank_radius OR endpoint y isn't ~tank_top/bottom, skip this entry
			if(!((abs(endpointradius-tank_radius)<3.)||(abs(abs(endpointy)-tank_halfheight)<15.))) continue;
			
			
			mustartvtx=primarystartvertex.Vect();
			mustopvtx=primarystopvertex.Vect();
			mustartE=nextrack->GetE();
			mustopE=nextrack->GetEndE();
			
			Float_t oppx = primarystopvertex.X() - primarystartvertex.X();
			Float_t adj = primarystopvertex.Z() - primarystartvertex.Z();
			Float_t avgtrackanglex = TMath::ATan(oppx/adj);
			Float_t oppy = primarystopvertex.Y() - primarystartvertex.Y();
			Float_t avgtrackangley = TMath::ATan(oppy/adj);
			
			TVector3 differencevector  = (primarystopvertex.Vect()-primarystartvertex.Vect());
			fsltruetracklength->Fill(differencevector.Mag());
			
			// ----------------------------------------------------------------------------------------------
			// retrieve any remaining information and record the event
			// ----------------------------------------------------------------------------------------------
			Float_t mu_rest_mass_E = 105.658*1000.;
			mustartKE = (mustartE-mu_rest_mass_E);  // starting energy (GeV) (p^2+m^2)
			mustopKE = (mustopE-mu_rest_mass_E);
			Float_t primarymomentummag = nextrack->GetP();                    // starting momentum
			TVector3 primarymomentumdir(nextrack->GetPdir(0),nextrack->GetPdir(1),nextrack->GetPdir(2));
			Float_t starttrackanglex = TMath::ATan(primarymomentumdir.X()/primarymomentumdir.Z());
			Float_t starttrackangley = TMath::ATan(primarymomentumdir.Y()/primarymomentumdir.Y());
			Float_t starttrackangle = TMath::Max(starttrackanglex,starttrackangley);
			
			Int_t primarytrackid = nextrack->GetId();
			
#ifdef VERBOSE
			cout<<"found a suitable primary track"<<endl;
#endif
			
			// ----------------------------------------------------------------------------------------------
			// calculate the track length in water FIXME: sometimes x sign is wrong, sometimes all are way out
			// ----------------------------------------------------------------------------------------------
			// to calculate track length _in water_ find distance from start vertex to the point
			// where it intercepts the tank. if this length > total track length, use total track length
			// otherwise use this length. 
			
			// first find out the z value where the tank would leave the radius of the tank
			// in the sample we're giving it, they were killed immmediately after exiting tank 
			// but do both methods for comparison
			Double_t xatziszero = 
			(mustartvtx.X() - (mustartvtx.Z()-tank_start-tank_radius)*TMath::Tan(avgtrackanglex));
			Double_t firstterm = -TMath::Tan(avgtrackanglex)*xatziszero;
			Double_t thirdterm = 1+TMath::Power(TMath::Tan(avgtrackanglex),2.);
			Double_t secondterm = (TMath::Power(tank_radius,2.)*thirdterm) - TMath::Power(xatziszero,2.);
			Double_t solution1 = (firstterm + TMath::Sqrt(secondterm))/thirdterm;
			Double_t solution2 = (firstterm - TMath::Sqrt(secondterm))/thirdterm;
			Double_t tankendpointz;
			if(primarystopvertex.Z() > mustartvtx.Z()){
				tankendpointz = solution1;	//forward going track
			} else {
				tankendpointz = solution2;	// backward going track
			}
			Double_t tankendpointx = TMath::Sqrt(TMath::Power(tank_radius,2)-TMath::Power(tankendpointz,2));
			// correct for tank z offset (do after tankendpointx, before tankendpointy)
			tankendpointz += tank_start+tank_radius;
			// now check if the particle would have exited through one of the caps before reaching this radius
			Double_t tankendpointy = 
			mustartvtx.Y() + (tankendpointz-mustartvtx.Z())*TMath::Tan(avgtrackangley);
			
			if(TMath::Abs(tankendpointy-tank_yoffset)>(tank_halfheight)){
				// this trajectory exits through the cap. Need to recalculate x, z exiting points...!
				if(primarystopvertex.Y()>mustartvtx.Y()){
					tankendpointy = tank_halfheight+tank_yoffset;	// by definition of leaving through cap
				} else {
					tankendpointy = -tank_halfheight+tank_yoffset;
				}
				tankendpointz = 
				mustartvtx.Z() + (tankendpointy-mustartvtx.Y())/TMath::Tan(avgtrackangley);
				tankendpointx = 
				mustartvtx.X() + (tankendpointz-mustartvtx.Z())*TMath::Tan(avgtrackanglex);
			} else {
				// this trajectory exited the tank by a side point; existing value is valid
			}
			
			Double_t maxtanktracklength = 
				TMath::Sqrt(TMath::Power(tank_radius*2.,2.)+TMath::Power(tank_halfheight*2.,2.));
			
			// we're now able to determine muon track length in the tank:
			mutracklengthintank = TMath::Sqrt(
				TMath::Power((tankendpointx-mustartvtx.X()),2)+
				TMath::Power((tankendpointy-mustartvtx.Y()),2)+
				TMath::Power((tankendpointz-mustartvtx.Z()),2) );
			
//			cout<<"*******************************************************"<<endl;
//			cout<<"Manual method: ("<<mutracklengthintank<<", end point ("
//				<<tankendpointx<<", "<<tankendpointy<<", "<<tankendpointz<<")"<<endl;
			
			mutracklengthintank = TMath::Sqrt(
				TMath::Power((mustopvtx.X()-mustartvtx.X()),2)+
				TMath::Power((mustopvtx.Y()-mustartvtx.Y()),2)+
				TMath::Power((mustopvtx.Z()-mustartvtx.Z()),2) );
//			cout<<"True method:   ("<<mutracklengthintank<<", end point ("
//				<<mustopvtx.X()<<", "<<mustopvtx.Y()<<", "<<mustopvtx.Z()<<")"<<endl;
//			cout<<"*******************************************************"<<endl;
			
#ifdef MUTRACKDEBUG
			cout<<"muon tank track length: ("<<(tankendpointx-primarystartvertex.X())
				<<", "<<(tankendpointy-primarystartvertex.Y())<<", "
				<<(tankendpointz-primarystartvertex.Z())<<") = "<<mutracklengthintank<<"cm total"<<endl;
			cout<<"muon tank exit point: ("<<tankendpointx<<", "<<tankendpointy<<", "<<tankendpointz<<") ";
			cout<<"muon start point : ("<<primarystartvertex.X()<<", "<<primarystartvertex.Y()
				<<", "<<primarystartvertex.Z()<<")"<<endl;
#endif
			if(mutracklengthintank > maxtanktracklength){
				cout<<"Track length is impossibly long!"<<endl;
				assert(false);
			}
			if(TMath::IsNaN(mutracklengthintank)){
				cout<<"NaN RESULT FROM MU TRACK LENGTH IN TANK?!"<<endl;
				assert(false);
			}
			fsltruetracklengthintank->Fill(mutracklengthintank);
			
			// ----------------------------------------------------------------------------------------------
			// digit analysis
			// ----------------------------------------------------------------------------------------------
#ifdef VERBOSE
			cout<<"Analysing tank digits"<<endl;
#endif
			Int_t numdigitsthisevent = atrigt->GetCherenkovDigiHits()->GetEntries();
			//cout<<"this event has "<<numdigitsthisevent<<" digits"<<endl;
			Double_t chargeinsidecherenkovcone=0;
			Double_t chargeoutsidecherenkovcone=0;
			tankchargefrommuon=0.;
			numTankDigits=numdigitsthisevent;
			numtankdigitsfrommuon=0;
			tanktubeshitbymu.clear();
			upstreamcharge=0.;
			downstreamcharge=0.;
			topcapcharge=0.;
			bottomcapcharge=0.;
			for(Int_t i=0; i<numdigitsthisevent; i++){
				// retrieve the digit information
				/////////////////////////////////
				WCSimRootCherenkovDigiHit* thedigihit = 
					(WCSimRootCherenkovDigiHit*)atrigt->GetCherenkovDigiHits()->At(i);
				int digitstubeid = thedigihit->GetTubeId()-1;
				double digitsq = thedigihit->GetQ();
				double digitst = thedigihit->GetT(); // this is time within the trigger window + 950ns
				std::vector<int> truephotonindices = thedigihit->GetPhotonIds();
				WCSimRootPMT pmt = geo->GetPMT(digitstubeid);
				double digitsx = pmt.GetPosition(0);
				double digitsy = pmt.GetPosition(1);
				double digitsz = pmt.GetPosition(2);
				int thepmtsloc = pmt.GetCylLoc();
				switch (thepmtsloc){
					case 0: topcapcharge+=digitsq; break;
					case 2: bottomcapcharge+=digitsq; break;
					case 1: ((digitsz-tank_start-tank_radius)<0) ? upstreamcharge+=digitsq : downstreamcharge+=digitsq; break;
				}
				
				// Make a plot of just the light from this muon, and the light from everything else
				///////////////////////////////////////////////////////////////////////////////////
				TVector3 digitvector(digitsx-primarystartvertex.X(), digitsy-primarystartvertex.Y(), 
					digitsz-primarystartvertex.Z());
				Double_t dotproduct = (digitvector.Unit()).Dot(differencevector.Unit());
				Double_t digitmuonangle = TMath::ACos(dotproduct);  // returns [0, Pi]
				if(digitmuonangle<((42.0*TMath::Pi())/180.0)){
					chargeinsidecherenkovcone+=digitsq;
					//FillTankMapHist(geo, digitstubeid, true, maphistos, digitsq);
					//hitpositions.Fill(digitsx, digitsz, digitsy);
				} else {
					chargeoutsidecherenkovcone+=digitsq;
					//FillTankMapHist(geo, digitstubeid, false, maphistos, digitsq);
				}
				
				// Calculate total charge from muon
				////////////////////////////////////
				// scan through the parents ID's for the photons contributing to this digit
				// and see if any of them are this muon - scale digit charge by this fraction
				/*int numphotonsfromthismuon=0;
				for(int truephoton=0; truephoton<truephotonindices.size(); truephoton++){
					int thephotonsid = truephotonindices.at(truephoton);
					WCSimRootCherenkovHitTime *thehittimeobject = 
						(WCSimRootCherenkovHitTime*)atrigt->GetCherenkovHitTimes()->At(thephotonsid);
					Int_t thephotonsparenttrackid = thehittimeobject->GetParentID();
					if(thephotonsparenttrackid==nextrack->GetId()) numphotonsfromthismuon++;
				}
				// to estimate the charge from the muon we should scale each digit's total charge
				// by the fraction of hits in the digit which come from the muon
				if(numphotonsfromthismuon!=0){   //this digit had contribution from the muon!
					//cout<<"muon digit "<<i<<" had charge "<< digitsq << " and " << truephotonindices.size()
					//	<<" true photons, of which "<<numphotonsfromthismuon <<" were from this muon"<<endl
					tankchargefrommuon += 
						digitsq * ((double)numphotonsfromthismuon/(double)truephotonindices.size());
					numtankdigitsfrommuon++;
					if(std::find(tanktubeshitbymu.begin(), tanktubeshitbymu.end(), digitstubeid)==tanktubeshitbymu.end())
						tanktubeshitbymu.push_back(digitstubeid);
				} */
				// Replace: everything in this sample comes from the muon - it's muon particle gun MC.
				tankchargefrommuon+=digitsq;
				numtankdigitsfrommuon++;
				if(std::find(tanktubeshitbymu.begin(), tanktubeshitbymu.end(), digitstubeid) == 
						tanktubeshitbymu.end()) tanktubeshitbymu.push_back(digitstubeid);
				
			}  // end loop over digits
			
			// debug: check the in and out of cone charge maps to check digits are being assigned correctly
			/*
			if(hitpositions.GetEntries() > 10){
				TCanvas c1("c1");
				c1.cd();
				TH1* histowall=(TH1*)maphistos.at("inconehistowall");
				histowall->Draw("colz");
				TCanvas c4("c4");
				c4.cd();
				histowall=maphistos.at("outconehistowall");
				histowall->Draw("colz");
				TCanvas c5("c5");
				c5.cd();
				hitpositions.SetMarkerStyle(20);
				hitpositions.Draw();
				gPad->WaitPrimitive();
				//std::this_thread::sleep_for (std::chrono::seconds(15));  // wait so we can look at histos
			}
			*/
			
			muedepositionswcsim->Fill(tankchargefrommuon);
			totaltankcharge=atrigt->GetSumQ();
			
			if((chargeinsidecherenkovcone+chargeoutsidecherenkovcone)!=0){
				fractionalchargeincone=
				chargeinsidecherenkovcone/(chargeinsidecherenkovcone+chargeoutsidecherenkovcone);
			} else {
				fractionalchargeincone=0;
			}
			
			// end of digit analysis
			////////////////////////////////////////////////
			
			treeout->Fill();
			break;
			
		}  // end of loop over tracks
		
	}  // end of loop over events
	
	//======================================================================================================
	//======================================================================================================
	
	gROOT->cd();
	
//	 draw other canvases
#ifdef WCSIMDEBUG
	TCanvas* debug1 = new TCanvas("debug1","debug1");
	debug1->cd();
	TLegend* aleg = new TLegend(0.5,0.78,0.7,0.88);
	aleg->AddEntry(fsltruetracklength,"Total track length","l");
	aleg->AddEntry(fsltruetracklengthintank,"Tank track length","l");
	fsltruetracklengthintank->SetLineColor(kRed);
	fsltruetracklengthintank->Draw();
	fsltruetracklength->SetLineColor(kBlue);
	fsltruetracklength->Draw("same");
	aleg->Draw();
	debug1->SaveAs("muon_track_lengths.png");
#endif
	if(fsltruetracklength) delete fsltruetracklength; fsltruetracklength=0;
	if(fsltruetracklengthintank) delete fsltruetracklengthintank; fsltruetracklengthintank=0;
	
#ifdef DRAW_HISTOS
	debug1->Clear();
	aleg->Clear();
	aleg->AddEntry(muedepositionswcsim,"Mu Tank E Dep (all)","l");
	muedepositionswcsim->SetLineColor(kBlue);
	muedepositionswcsim->Draw();
	aleg->Draw();
	debug1->SaveAs("muon_tank_energy_depositions.png");
	if(muedepositionswcsim) delete muedepositionswcsim; muedepositionswcsim=0;
	
	if(debug1) delete debug1; debug1=0;
	if(aleg) delete aleg; aleg=0;
	
#endif  // DRAW_HISTOS
	
	// cleanup
	// =======
	cout<<"cleanup"<<endl;
	if(c) c->ResetBranchAddresses();
	//cout<<"deleting wcsim chain"<<endl;
	if(c) delete c; c=0;					// ??
	
	// Save all the remaining histograms
	fileout->cd();
	std::vector<TH1*> otherhistos{
	(TH1*)muedepositionswcsim,
	(TH1*)muedepositionsfidcut,
	(TH1*)fsltruetracklength,
	(TH1*)fsltruetracklengthintank,
	(TH1*)tracklengthvsmuonlight,
	(TH1*)inconehistowall,
	(TH1*)inconehistotop,
	(TH1*)inconehistobottom,
	(TH1*)outconehistowall,
	(TH1*)outconehistotop,
	(TH1*)outconehistobottom
	};
	for(auto thehist : otherhistos){
		if(thehist) { thehist->Write(); delete thehist; }
	}
	if(histofileout) histofileout->Close(); delete histofileout; histofileout=0;
	
	// write and close file of event information
	treeout->SetEntries(bWCSimEventNum->GetEntries());
	treeout->Write();
	fileout->Close();
	delete fileout;
	fileout=0;
}

void ColourPlotStyle(){
	const Int_t NRGBs = 5;
	const Int_t NCont = 255;

	Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
	Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
	Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
	Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
	gStyle->SetNumberContours(NCont);
}

void FillTankMapHist(WCSimRootGeom* geo, int tubeID, bool incone, std::map<std::string, TH2D*> &maphistos, double weight=1){
	//Fill a bin on a 2D map of PMTs 
	WCSimRootPMT pmt = geo->GetPMT(tubeID);
	// WCSimRootPMT has members GetTubeNo(), GetCylLoc(), GetPosition(j), GetOrientation(j)
	// GetCylLoc(): 0=top cap, 2=bottom cap, 1=wall, 4=mrd, 5=veto, 3=obselete outer veto (shouldnt come up)
	// GetPosition(j), j=0..2: Returns x,y,z coordinates of the center of the sphere that forms the PMT.
	// GetOrientation(j), j=0..2: Returns the x,y,z components of a vector of the direction the PMT faces.
	// GetPMT(j) Returns a pmt object - NOT a pointer to a PMT object.
	//cout<<"Filling histogram for cylloc "<<pmt.GetCylLoc()<<" for tubeID "<<tubeID<<endl;
	switch(pmt.GetCylLoc()){
		case 0: {
			if(topcappositionmap.count(tubeID)){
				std::pair<int,int> thebins = topcappositionmap.at(tubeID);
				TH2D* histotop;
				if(incone){
					histotop=maphistos.at("inconehistotop");
				} else {
					histotop=maphistos.at("outconehistotop");
				}
				if(histotop) histotop->Fill(thebins.first, thebins.second, weight);
			} else {cout<<"bad pmt: ID "<<tubeID<<" in CylLoc "<<pmt.GetCylLoc()<<endl;}
			break;
		}
		case 1: {
			if(wallpositionmap.count(tubeID)){
				std::pair<int,int> thebins = wallpositionmap.at(tubeID);
				TH2D* histowall;
				if(incone) {
					histowall=maphistos.at("inconehistowall");
				} else {
					histowall=maphistos.at("outconehistowall");
				}
				if(histowall) histowall->Fill(thebins.first+0.5, thebins.second, weight);
			} else {cout<<"bad pmt: ID "<<tubeID<<" in CylLoc "<<pmt.GetCylLoc()<<endl;}
			break;
		}
		case 2: {
			if(bottomcappositionmap.count(tubeID)){
				std::pair<int,int> thebins = bottomcappositionmap.at(tubeID);
				TH2D* histobottom;
				if(incone) {
					histobottom=maphistos.at("inconehistobottom");
				} else {
					histobottom=maphistos.at("outconehistobottom");
				}
				if(histobottom) histobottom->Fill(thebins.first, thebins.second, weight);
			} else {cout<<"bad pmt: ID "<<tubeID<<" in CylLoc "<<pmt.GetCylLoc()<<endl;}
			break;
		}
		case 4: {
//				std::pair<int,int> thebins = mrdpositionmap.at(tubeID);
//				mrdhist->Fill(thebins.first, thebins.second, weight);
			break;
		}
		case 5: {
//				std::pair<int,int> thebins = faccpositionmap.at(tubeID);
//				facchist->Fill(thebins.first, thebins.second, weight);
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
