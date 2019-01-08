/* vim:set noexpandtab tabstop=4 wrap */
/*TODO: to modify for new source files:
0. Re-generate and replace in this directory libWCSimRoot.so, .rootmap, .pcm files.
1. Re-enable #include RootOptions.hh
*/
#ifndef VERBOSE
#define VERBOSE
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
#include "TStyle.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include <regex>
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

#include "MRDspecs.hh"
#include "MRDSubEvent_ReconstructionClasses.hh"

void set_plot_style();

const Int_t numtankpmts=128+2*(26);             // 26 pmts and lappds on each cap
const Int_t numvetopmts=26;
const Int_t caparraysize=8;                     // pmts on the cap form an nxn grid where caparraysize=n
const Int_t pmtsperring=16;                     // pmts around each ring of the main walls
const Int_t numpmtrings=8;                      // num rings around the main walls

const Float_t fidcutradius=tank_radius*0.8;     // fiducial volume is slightly lesss than the full tank radius
const Float_t fidcuty=50.;                      // a meter total fiducial volume in the centre y
const Float_t fidcutz=0.;                       // fidcuial volume is before the tank centre.

const char* wcsimpath="/annie/app/users/moflaher/wcsim/root_work/in/MRD_muon_sample";
const char* wcsimfilestring="wcsim_MRD_muon_sample_*.root";
const char* wcsimlibrarypath="/annie/app/users/moflaher/wcsim/wcsim/libWCSimRoot.so";
const char* outpath="/annie/app/users/moflaher/wcsim/root_work";

std::map<int,int> pmtidtocopynum;
std::map<int,int> makecopynummap(WCSimRootGeom* geo);

// for distance to pmt, efficiency by paddle/layer etc. This function looks up details by PMT number
bool ispaddlehoriz(int copyNo, int &panelnum, int &paddlenum);

void truthtracks(){
	
	gStyle->SetOptStat(0);
	
	// load WCSim library for reading WCSim files
	cout<<"loading "<<wcsimlibrarypath<<endl;
	gSystem->Load(wcsimlibrarypath);
	
	// ==============================================================================================
	// ==============================================================================================
	// Declare input paths, files, trees...
	TFile* wcsimfile=0;
	TTree* wcsimT=0;
	
	// wcsimT
	WCSimRootEvent* b=0, *m=0, *v=0;
	TBranch* bp=0, *mp=0, *vp=0;
	WCSimRootTrigger* atrigt=0, *atrigm=0, *atrigv=0;
	
	// geoT
	WCSimRootGeom* geo = 0; 
	
	// ===================================================================================================
	// ===================================================================================================
	// loading and processing
	gROOT->cd();

	// TChain for wcsim files - this will be the main driver of the loop - all it's events will be processed.
	TChain* c =  new TChain("wcsimT");
	TString chainpattern = TString::Format("%s/%s",wcsimpath,wcsimfilestring);
	cout<<"loading TChain entries from "<<chainpattern<<endl;
	c->Add(chainpattern);
	Int_t numents = c->GetEntries();
	cout<<"loaded "<<numents<<" entries in the chain"<<endl;
	if(numents<1){ return; }
	Int_t treeNumber = -1;
	
	TFile* outfile = new TFile("MRDpenetration2.root","RECREATE");
	outfile->cd();
	
	/* Declare histograms */
	// ===================================================================================================
	// ===================================================================================================
	// rate of energy loss = Edep / length. Only have energy at start and ends, so can only do avg for
	// stopped mus. avg rate of loss for axial penetration, or total penetration
	TH2D* penetrationcmvsenergyall = new TH2D("penetrationcmvsenergyall","Muon MRD Energy vs Penetration (All Muons)",100,0,300,100,0,3000);
	TH2D* penetrationhiststopall = new TH2D("penetrationhiststopall","Muon MRD Energy vs Penetration (Stopping Muons)",100,0,300,100,0,3000);
	// split into angles
	TH2D* penetrationhiststop1 = new TH2D("penetrationhiststop1","Muon MRD Energy vs Penetration (Stopping Muons with Angle < 1/5*(Pi/2))",100,0,150,100,0,1500);
	TH2D* penetrationhiststop2 = new TH2D("penetrationhiststop2","Muon MRD Energy vs Penetration (Stopping Muons with Angle < 2/5*(Pi/2))",100,0,150,100,0,1500);
	TH2D* penetrationhiststop3 = new TH2D("penetrationhiststop3","Muon MRD Energy vs Penetration (Stopping Muons with Angle < 3/5*(Pi/2))",100,0,150,100,0,1500);
	TH2D* penetrationhiststop4 = new TH2D("penetrationhiststop4","Muon MRD Energy vs Penetration (Stopping Muons with Angle < 3.5/5*(Pi/2))",100,0,150,100,0,1500);
	TH2D* penetrationhiststop5 = new TH2D("penetrationhiststop5","Muon MRD Energy vs Penetration (Stopping Muons with Angle < 4/5*(Pi/2))",100,0,150,100,0,1500);
	TH2D* penetrationhiststop6 = new TH2D("penetrationhiststop6","Muon MRD Energy vs Penetration (Stopping Muons with Angle < 4.5/5*(Pi/2))",100,0,150,100,0,1500);
	TH2D* penetrationhiststop7 = new TH2D("penetrationhiststop7","Muon MRD Energy vs Penetration (Stopping Muons with Angle < 5/5*(Pi/2))",100,0,150,100,0,1500);
	// make a vector of the angles and their errors
	// TODO: OK, for this to actually be sensible, the angular bins need to be the bins which will
	// actually be generated by the discrete paddle widths. 
	std::vector<double> theangles{0.5/5.*(TMath::Pi()/2.),1.5/5.*(TMath::Pi()/2.),2.5/5.*(TMath::Pi()/2.),3.25/5.*(TMath::Pi()/2.),3.75/5.*(TMath::Pi()/2.),4.25/5.*(TMath::Pi()/2.),4.75/5.*(TMath::Pi()/2.)};
	std::vector<double> theangleerrors{0.5/5.*(TMath::Pi()/2.),0.5/5.*(TMath::Pi()/2.),0.5/5.*(TMath::Pi()/2.),0.25/5.*(TMath::Pi()/2.),0.25/5.*(TMath::Pi()/2.),0.25/5.*(TMath::Pi()/2.),0.25/5.*(TMath::Pi()/2.)};
	TH2D* penetrationhistpeneall = new TH2D("penetrationhistpeneall","Muon MRD Energy vs Penetration (Fully Penetrating Muons)",100,0,300,100,0,3000);
	TH2D* penetrationhistotherall = new TH2D("penetrationhistotherall","Muon MRD Energy vs Penetration (Side Exiting Muons)",100,0,300,100,0,3000);
	
	// total detected charge vs muon E
	TH2D* chargevsedephist = new TH2D("chargevsedephist","Total Charge vs Muon Energy",100,0,1500,100,0,250);
	// chargevsedephist->Fill(numphots, muE);
	
	// faction of total charge deposited vs layer
	TH2D* fracofchargedephist = new TH2D("fracofchargedephist","Fraction of Charge Deposition vs Layer",12,0,12,10,0.,1.);
	// fracofchargedephist->Fill(dQ,layer);
	
	// charge per unit energy vs distance from pmt
	TH2D* chargevsdisthist = new TH2D("chargevsdisthist","Charge vs Distance from PMT",100,0,500,100,0,2);
	// chargevsdisthist->Fill(PMTQ/E, distfrompmt);
	
	TH1D* startposx = new TH1D("startposx","Start Position X",40,-20,20);
	TH1D* startposy = new TH1D("startposy","Start Position Y",40,tank_yoffset-20,tank_yoffset+20);
	TH1D* startposz = new TH1D("startposz","Start Position Z",20,MRD_start-10,MRD_start+10);
	TH1D* startdirx = new TH1D("startdirx","Start dirition X",100,-1.1,1.1);
	TH1D* startdiry = new TH1D("startdiry","Start dirition Y",100,-1.1,1.1);
	TH1D* startdirz = new TH1D("startdirz","Start dirition Z",100,-0.1,1.1);
	TH1D* startangleh = new TH1D("startangleh","Start Angle Horiz", 100, -TMath::Pi()*0.6, TMath::Pi()*0.6);
	TH1D* startanglev = new TH1D("startanglev","Start Angle Vert", 100, -TMath::Pi()*0.6, TMath::Pi()*0.6);
	
	TH1D* muonanglehist = new TH1D("muonanglehist","Muon Angles", 100, 0., TMath::Pi());
	
	TH1I* tubeshithist = new TH1I("tubeshithist","MRD PMT Tube IDs Hit",501,0,500);
	TH1I* layershithist = new TH1I("layershithist","MRD Layers Hit",12,0,11);
	
	// TODO: efficiency won't work with the Muon sample distribution, because they were generated starting from a small box 
	// in front centre of the MRD. Need to do efficiency measurement either with beam events, with maybe an isotropic
	// sample all around the MRD. 
	// detection efficiency vs x position (for each paddle)
	// MRD should have 90% efficiency > MRD thesis p26
	// Use the same range for all horiz / vertical paddles so that they can easily be drawn on top of each other
	// to mimick sciboone plot, have a graph per paddle. That's a lot of graphs!
	/*
	std::vector<std::vector<TH1D*>> efficiencygraphs;  // outer vector for layer, inner vector for paddle
	for(int pmtit=0; pmtit<nummrdpmts; pmtit++){
		// we'll need the layer and paddle nums
		int layernum=-1, paddlenum=-1;
		bool ishpaddle = ispaddlehoriz(pmtit, layernum, paddlenum);
		//std::pair<Double_t,Double_t> xextents = mrdcluster::paddle_extentsx.at(pmtit);
		//std::pair<Double_t,Double_t> yextents = mrdcluster::paddle_extentsy.at(pmtit);
		TH1D* histpoint = new TH1D(TString::Format("scinteff%d",pmtit),TString::Format("Efficiency Vs Position for Paddle %d",pmtit),150,-MRD_steel_width,MRD_steel_width);
		efficiencygraphs.at(layernum).at(paddlenum)=histpoint;
	}
	*/
	std::vector<int> numpaddlepassthrus(nummrdpmts,0);
	std::vector<int> numdetectedpassthrus(nummrdpmts,0);
	
	// TODO
	// run MRD sim through track reconstruction and see how effectively MRD track reconstruction is working
	// muon angle reconstruction with errors vs true
	
	/* Loop over WCSim file */
	// ===================================================================================================
	// ===================================================================================================
	cout<<"looping over wcsimT entries"<<endl;
//	numents=10;
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
			wcsimfilestring=wcsimfile->GetName();
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
		
		wcsimT->GetEntry(inputEntry);
		atrigt = b->GetTrigger(0);
		atrigm = m->GetTrigger(0);
		
		Int_t numtracks = atrigt->GetNtrack();
#ifdef VERBOSE
		cout<<"wcsim event had "<<numtracks<<" truth tracks"<<endl;
#endif
		
		/* 2. Load tracks, look for primary mu track through the MRD, record interaction details. */
		//===============================================================================================
		//===============================================================================================
		Int_t nummuontracksintank=0;
		bool muonentersMRD=false;
		bool muonstopsinMRD=false;
		bool muonrangesoutMRD=false;
		double mrdpenetrationcm=0.;
		int mrdpenetrationlayers=0;
		double mutracklengthinMRD=0.;
		// now scan through the truth tracks, find the primary muon and save the wcsim info from it
		TVector3 primarystartvertex(0,0,0);
		TVector3 primarystopvertex(0,0,0);
		Float_t avgtrackanglex, avgtrackangley;
		Double_t muonangle;
		
		// find the primary muon track:
		int maxmutrackindex=-1;
		int mutrackid=-1;
		Double_t muonenergy=0;        // we'll keep updating the info so it only records the highest energy muon
		int numprimarymus=0;
		for(int track=0; track<numtracks; track++){
			WCSimRootTrack* nextrack = (WCSimRootTrack*)atrigt->GetTracks()->At(track);
//#ifdef VERBOSE
//			cout<<"track "<<track<<" has pdg "<<nextrack->GetIpnu()<<" and ";
//			(nextrack->GetParenttype()) ? cout<<"is " : cout<<"is not " ;
//			cout<<"a primary"<<endl;
//#endif
			if(nextrack->GetFlag()!=0) continue;                // neutrino and neutrino target, not a track
			if(TMath::Abs(nextrack->GetIpnu())!=13) continue;   // not a muon
			if(nextrack->GetParenttype()!=0) continue;
			numprimarymus++;
			
//#ifdef VERBOSE
//			cout<<"primary muon found with energy "<<nextrack->GetE()<<"MeV, "<<endl;
//			cout<<"primarystartvol is "<<nextrack->GetStartvol()<<", startpoint is ("
//			    <<nextrack->GetStart(0)<<","<<nextrack->GetStart(1)<<","<<nextrack->GetStart(2)<<")"<<endl;
//			cout<<"primarystopvol is "<<nextrack->GetStopvol()<<", endpoint is ("
//			    <<nextrack->GetStop(0)<<","<<nextrack->GetStop(1)<<","<<nextrack->GetStop(2)<<")"<<endl;
//#endif
			if(!(muonenergy>nextrack->GetE())){
				muonenergy=nextrack->GetE();
				maxmutrackindex=track;
				mutrackid=nextrack->GetId();
			} else {
				continue;
			}
		}
		
		if(maxmutrackindex<0) continue;
#ifdef VERBOSE
			cout<<"found "<<numprimarymus<<" primary muons, the highest energy being track "<<maxmutrackindex
			    <<" which had an energy "<<muonenergy<<"MeV"<<endl;
#endif
		//for(int track=0; track<numtracks; track++){
		for(int track=maxmutrackindex; track==maxmutrackindex; track++){  // just one iteration with the primary mu.
			WCSimRootTrack* nextrack = (WCSimRootTrack*)atrigt->GetTracks()->At(track);
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
			
			// ----------------------------------------------------------------------------------------------
			// Check if it's a primary muon starting in the tank
			// ----------------------------------------------------------------------------------------------
			// is it a (anti)muon?
			Int_t primarypdg = nextrack->GetIpnu();
			if(TMath::Abs(primarypdg)!=13) continue;       // not a muon
			
			// is it a primary?
			if(nextrack->GetParenttype()!=0) continue;
			
#ifdef VERBOSE
			cout<<"primarystartvol is "<<nextrack->GetStartvol()<<", startpoint is ("
			    <<nextrack->GetStart(0)<<","<<nextrack->GetStart(1)<<","<<nextrack->GetStart(2)<<")"<<endl;
			cout<<"initial direction is ("<<nextrack->GetDir(0)<<","<<nextrack->GetDir(1)
			    <<","<<nextrack->GetDir(2)<<")"<<endl;
#endif
			
			// does it stop in the mrd, or pass completely through the MRD? 
			Int_t primarystopvol = nextrack->GetStopvol();
			// temporary override as these weren't correctly set in wcsim <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			//TODO REMOVE THIS once 1M generated. 
			// BUT need to think about how to check for 'range-out' mrd events. Maybe this is fine?
			if(nextrack->GetStop(2)<tank_start){
				primarystopvol = 20;                        // stop depth is facc or hall
			} else if(nextrack->GetStop(2)>(tank_start+(2.*tank_radius))){
				primarystopvol = 30;                        // stop depth is mrd or hall
			} else {
				primarystopvol = 10;                        // stop depth is tank
			}
#ifdef VERBOSE
			cout<<"primarystopvol is "<<primarystopvol<<", endpoint is ("
			    <<nextrack->GetStop(0)<<","<<nextrack->GetStop(1)<<","<<nextrack->GetStop(2)<<")"<<endl;
#endif
			
			primarystartvertex = TVector3(  nextrack->GetStart(0),
												nextrack->GetStart(1),
												nextrack->GetStart(2));
			primarystopvertex = TVector3(   nextrack->GetStop(0),
												nextrack->GetStop(1),
												nextrack->GetStop(2));
			
			Float_t oppx = primarystopvertex.X() - primarystartvertex.X();
			Float_t adj = primarystopvertex.Z() - primarystartvertex.Z();
			avgtrackanglex = TMath::ATan(oppx/adj);
			Float_t oppy = primarystopvertex.Y() - primarystartvertex.Y();
			avgtrackangley = TMath::ATan(oppy/adj);
			TVector3 differencevector  = primarystopvertex - primarystartvertex;
			TVector3 azaxisvector(0,0,1);
			muonangle=differencevector.Angle(azaxisvector);
#ifdef VERBOSE
			cout<<"muon angle is "<<((180./TMath::Pi())*muonangle)<<" degrees"<<endl;
#endif
			
			startposx->Fill(nextrack->GetStart(0));
			startposy->Fill(nextrack->GetStart(1));
			startposz->Fill(nextrack->GetStart(2));
			startdirx->Fill(nextrack->GetDir(0));
			startdiry->Fill(nextrack->GetDir(1));
			startdirz->Fill(nextrack->GetDir(2));
			muonanglehist->Fill(muonangle);
			
			TVector3 horizplane(nextrack->GetDir(0),0,nextrack->GetDir(2));
			double horizang=horizplane.Angle(azaxisvector);
			(nextrack->GetDir(0) < 0) ? horizang*=-1. : horizang*=1. ;
			TVector3 vertplane(0,nextrack->GetDir(1),nextrack->GetDir(2));
			double vertang=vertplane.Angle(azaxisvector);
			(nextrack->GetDir(1) < 0) ? vertang*=-1. : vertang*=1. ;
			startangleh->Fill(horizang);
			startanglev->Fill(vertang);
			
			// ----------------------------------------------------------------------------------------------
			// calculate muon MRD penetration 
			// ----------------------------------------------------------------------------------------------
			
			if(primarystopvertex.Z()<MRD_start){
				// the track stops before reaching the MRD
				muonentersMRD=false;
				muonstopsinMRD=false;
				muonrangesoutMRD=false;
				mrdpenetrationcm=0.;
				mrdpenetrationlayers=0;
				mutracklengthinMRD=0.;
			} else {
				// the track travels a depth past the MRD start, but may do so outside its angular acceptance
				Float_t trackZlengthbeforeMRDXexit;
				if(TMath::Abs(primarystopvertex.X())<MRD_steel_width){
					trackZlengthbeforeMRDXexit=primarystopvertex.Z()-primarystartvertex.Z();
				} else {
					double trackXdistanceinMRD=0.;
					if(primarystopvertex.X()>0){
						trackXdistanceinMRD=MRD_steel_width-primarystartvertex.X();
					} else {
						trackXdistanceinMRD=-MRD_steel_width-primarystartvertex.X();
					}
					trackZlengthbeforeMRDXexit= trackXdistanceinMRD/TMath::Tan(avgtrackanglex);
				}
				Float_t trackZlengthbeforeMRDYexit;
				if(TMath::Abs(primarystopvertex.Y())<MRD_steel_height){
					trackZlengthbeforeMRDYexit=primarystopvertex.Z()-primarystartvertex.Z();
				} else {
					double trackYdistanceinMRD=0.;
					if(primarystopvertex.Y()>0){
						trackYdistanceinMRD=MRD_steel_height-primarystartvertex.Y();
					} else {
						trackYdistanceinMRD=-MRD_steel_height-primarystartvertex.Y();
					}
					trackZlengthbeforeMRDYexit= trackYdistanceinMRD/TMath::Tan(avgtrackangley);
				}
				Float_t trackZlengthbeforeMRDexit=
					TMath::Min(trackZlengthbeforeMRDXexit,trackZlengthbeforeMRDYexit);
				double mrdexitpoint=primarystartvertex.Z()+trackZlengthbeforeMRDexit;
				if(mrdexitpoint<MRD_start){
					// track has such a steep angle it exits MRD acceptance before reaching its start
					muonentersMRD=false;
					muonstopsinMRD=false;
					muonrangesoutMRD=false;
					mrdpenetrationcm=0.;
					mrdpenetrationlayers=0;
					mutracklengthinMRD=0.;
				} else {
					muonentersMRD=true;
					mrdpenetrationcm=trackZlengthbeforeMRDexit-(MRD_start-primarystartvertex.Z());
					mrdpenetrationlayers=0;
					for(auto layerzval : mrdscintlayers){
						if(mrdexitpoint<layerzval) break;
						mrdpenetrationlayers++;
					}
					if(mrdexitpoint>mrdscintlayers.back()){
						muonrangesoutMRD=true;
						muonstopsinMRD=false;
					} else {
						muonrangesoutMRD=false;
						if((TMath::Abs(primarystopvertex.X())<MRD_steel_width)
							&&(TMath::Abs(primarystopvertex.Y())<MRD_steel_height)) muonstopsinMRD=true;
					}
					
					// calculate total distance travelled in MRD
					////////////////////////////////////////////
					double muXdistanceinMRD, muYdistanceinMRD;
					double muXexitpoint, muYexitpoint;
					double muXentrypoint, muYentrypoint;
					muXentrypoint = primarystartvertex.X() + 
						(MRD_start-primarystartvertex.Z())*TMath::Tan(avgtrackanglex);
					muYentrypoint = primarystartvertex.Y() + 
						(MRD_start-primarystartvertex.Z())*TMath::Tan(avgtrackangley);
					if(trackZlengthbeforeMRDexit==trackZlengthbeforeMRDXexit){
						// track exits the MRD through one of the sides
						muXexitpoint= (primarystopvertex.X()>0) ? MRD_steel_width : -MRD_steel_width;
						muYexitpoint=
							primarystartvertex.Y()+(trackZlengthbeforeMRDexit*TMath::Tan(avgtrackangley));
					} else {
						// track exits through the top or bottom of the MRD
						(primarystopvertex.Y()>0) ? muYexitpoint=MRD_steel_height : muYexitpoint=-MRD_steel_height;
						muXexitpoint=
							primarystartvertex.X()+(trackZlengthbeforeMRDexit*TMath::Tan(avgtrackanglex));
					}
					muXdistanceinMRD=muXexitpoint-muXentrypoint;
					muYdistanceinMRD=muYexitpoint-muYentrypoint;
					mutracklengthinMRD=
						TMath::Sqrt(TMath::Power(muXdistanceinMRD,2)+TMath::Power(muYdistanceinMRD,2)
						+TMath::Power(mrdpenetrationcm,2));
				}
			}
#ifdef VERBOSE
			if(muonentersMRD){
				cout<<"muon enters MRD, ";
				if(muonstopsinMRD){
					cout<<"and stops in MRD";
				} else if (muonrangesoutMRD){
					cout<<"and ranges out MRD";
				} else {
					cout<<"and leaves by the sides of the MRD";
				}
				cout<<endl;
			}
			cout<<"penetration was "<<mrdpenetrationcm<<"cm or "<<mrdpenetrationlayers<<" layers"<<endl;
			cout<<"total track length in MRD was "<<mutracklengthinMRD<<"cm"<<endl;
#endif
			
		}  // end of loop over tracks - by now we've got the info on the primary mu
		
		
		// mrd digit analysis
		// ----------------------------------------------------------------------------------------------
		int numMRDdigits = atrigm->GetCherenkovDigiHits()->GetEntries();
		double totMRDcharge = atrigm->GetSumQ();
		int numMRDdigitsfrommu=0;
		double MRDchargefrommuunscaled=0.;
		double MRDchargefrommu=0.;
		std::vector<int> mrdtubeshitbymu;
		std::map<int,double> chargeoneachpmtfrommuon;
		for(Int_t i=0; i<numMRDdigits; i++){
			// retrieve the digit information
			/////////////////////////////////
			WCSimRootCherenkovDigiHit* thedigihit = 
				(WCSimRootCherenkovDigiHit*)atrigm->GetCherenkovDigiHits()->At(i);
			int digitstubeid = thedigihit->GetTubeId()-1;	//XXX index is TubeId()-1
			int layernum=-1, paddlenum=-1;
			bool ishpaddle = ispaddlehoriz(digitstubeid, layernum, paddlenum);
			double digitsq = thedigihit->GetQ();
			// Calculate total charge from muon
			////////////////////////////////////
			std::vector<int> truephotonindices = thedigihit->GetPhotonIds();
			int numphotonsfromthismuon=0;
			for(int truephoton=0; truephoton<truephotonindices.size(); truephoton++){
				int thephotonsid = truephotonindices.at(truephoton);
				WCSimRootCherenkovHitTime *thehittimeobject = 
					(WCSimRootCherenkovHitTime*)atrigm->GetCherenkovHitTimes()->At(thephotonsid);
				Int_t thephotonsparenttrackid = thehittimeobject->GetParentID();
				if(thephotonsparenttrackid==mutrackid) numphotonsfromthismuon++;
			}
			// scale each digit charge by the fraction of hits in the digit which come from the muon
			if(numphotonsfromthismuon!=0){
				tubeshithist->Fill(digitstubeid);
				layershithist->Fill(layernum);
				double chargeinthisdigitfrommu = 
					digitsq * ((double)numphotonsfromthismuon/(double)truephotonindices.size());
				MRDchargefrommu+= chargeinthisdigitfrommu;
				MRDchargefrommuunscaled+=digitsq;
				numMRDdigitsfrommu++;
				if(std::find(mrdtubeshitbymu.begin(), mrdtubeshitbymu.end(), digitstubeid) 
					== mrdtubeshitbymu.end()){
					mrdtubeshitbymu.push_back(digitstubeid);
					chargeoneachpmtfrommuon.emplace(digitstubeid, chargeinthisdigitfrommu);
				} else {
					chargeoneachpmtfrommuon.at(digitstubeid) = 
						chargeoneachpmtfrommuon.at(digitstubeid)+chargeinthisdigitfrommu;
				}
			}
			
		}  // end loop over mrd digits
#ifdef VERBOSE
		cout<<"muon left "<<MRDchargefrommu<<" units of charge in the MRD ("<<MRDchargefrommuunscaled
		    <<" units without scaling) over "<<numMRDdigitsfrommu
		    <<" digits in "<<mrdtubeshitbymu.size()<<" tubes"<<endl;
		cout<<"There were "<<numMRDdigits<<" digits in this event, with "<<numMRDdigitsfrommu<<" digits from primary muon"<<endl;
#endif
		
		// end of mrd digit analysis
		////////////////////////////////////////////////
		
		/* 2. Fill the histograms. */
		//===============================================================================================
		//===============================================================================================
		penetrationcmvsenergyall->Fill(mrdpenetrationcm,muonenergy);
		if(muonstopsinMRD){ penetrationhiststopall->Fill(mrdpenetrationcm,muonenergy);
			// split into angular bins
		     if(muonangle<((TMath::Pi()/2.)*(1./5.))){ penetrationhiststop1->Fill(mrdpenetrationcm,muonenergy); }
		else if(muonangle<((TMath::Pi()/2.)*(2./5.))){ penetrationhiststop2->Fill(mrdpenetrationcm,muonenergy); }
		else if(muonangle<((TMath::Pi()/2.)*(3./5.))){ penetrationhiststop3->Fill(mrdpenetrationcm,muonenergy); }
		else if(muonangle<((TMath::Pi()/2.)*(3.5/5.))){penetrationhiststop4->Fill(mrdpenetrationcm,muonenergy); }
		else if(muonangle<((TMath::Pi()/2.)*(4./5.))){ penetrationhiststop5->Fill(mrdpenetrationcm,muonenergy); }
		else if(muonangle<((TMath::Pi()/2.)*(4.5/5.))){penetrationhiststop6->Fill(mrdpenetrationcm,muonenergy); }
		else if(muonangle<((TMath::Pi()/2.)*(5./5.))){ penetrationhiststop7->Fill(mrdpenetrationcm,muonenergy); }
		}
		else if(muonrangesoutMRD) penetrationhistpeneall->Fill(mrdpenetrationcm,muonenergy);
		else if(muonentersMRD) penetrationhistotherall->Fill(mrdpenetrationcm,muonenergy);
		
		if(muonstopsinMRD) chargevsedephist->Fill(muonenergy, MRDchargefrommu);
		
		std::vector<int> numdigitslayers(numpanels,0);
		std::vector<double> layercharges(numpanels,0);
		for(int pmtit=0; pmtit<nummrdpmts; pmtit++){
			// If we have hits, we need distance to PMT to calculate paddle extinction
			// Whether we have hits or not, we need to know distance to PMT to know
			// if the paddle was passed through or not. We use this to determine detection
			// efficiency
			
			/* calculate the distance of the track to this PMT; project track to appropriate z */
			// first get PMT position
			WCSimRootPMT thepmt = geo->GetMRDPMT(pmtit);
			double thexval = thepmt.GetPosition(0);
			double theyval = thepmt.GetPosition(1);
			double thezval = thepmt.GetPosition(2);
			if(thezval>primarystopvertex.Z()&&(chargeoneachpmtfrommuon.count(pmtit)!=0)){
				//cerr<<"MRD Hit tube further downstream than end of track!?"<<endl; 
				//cerr<<"Hit Tube is at z="<<thezval<<", track end Z is "<<primarystopvertex.Z()<<endl; 
				//assert(false);
			}
			// project the track forward to the Z position of the PMT to get X and Y at that depth
			double projectedxposition = primarystartvertex.X() + 
				(thezval-primarystartvertex.Z())*TMath::Tan(avgtrackanglex);
			double projectedyposition = primarystartvertex.Y() + 
				(thezval-primarystartvertex.Z())*TMath::Tan(avgtrackangley);
			double thedistancefrommu=-1.;
			
			// Decide if we need to use the X or Y distance based on the PMT ID - i.e.
			// this looks up paddle orientation, layer and paddle number
			int layernum=-1, paddlenum=-1;
			bool ishpaddle = ispaddlehoriz(pmtit, layernum, paddlenum);
			if(layernum<0) cerr<<"COULD NOT FIND PMT "<<pmtit<<endl;
			//int paddleisvertical = mrdcluster::paddle_orientations.at(pmtit);  // alternative using mrdcluster statics
			//int layernum = mrdcluster::paddle_layers.at(pmtit);                // alternative using mrdcluster statics
			//                                                                   // doesn't store paddle number in layer
			//                                                                   // needed for efficiency plots.
			// is it a horizontal or vertical PMT?
			if(ishpaddle){
				thedistancefrommu = TMath::Abs(thexval-projectedxposition);
			} else {
				thedistancefrommu = TMath::Abs(theyval-projectedyposition);
			}
			// _should_ this paddle have seen a hit? Use this to establish efficiency of detection
			// need to compare the (X,Y) position of the track at that depth to the x and y extents of
			// the paddle corresponding to this PMT
			std::pair<Double_t,Double_t> xextents = mrdcluster::paddle_extentsx.at(pmtit);
			std::pair<Double_t,Double_t> yextents = mrdcluster::paddle_extentsy.at(pmtit);
			if( projectedxposition>xextents.first && projectedxposition<xextents.second
				&& projectedyposition>yextents.first && projectedyposition<yextents.second){
				// projected track position was within this paddle's bounds! It should have struck it!
				numpaddlepassthrus.at(pmtit)++;
				if(chargeoneachpmtfrommuon.count(pmtit)!=0) numdetectedpassthrus.at(pmtit)++;
			} else {
				if(chargeoneachpmtfrommuon.count(pmtit)!=0){
					cerr<<"HITS ON PMT EVEN THOUGH TRACK DID NOT PASS THROUGH THIS PMT!!"<<endl;
					cerr<<"projected position is ("<<projectedxposition<<", "<<projectedyposition<<")"<<endl
					    <<"paddle extent is ("<<xextents.first<<"->"<<xextents.second<<", "
					    <<yextents.first<<"->"<<yextents.second<<")"<<endl;
				}
			}
			
			if(chargeoneachpmtfrommuon.count(pmtit)!=0){
				// this PMT has some charge from muon hits
				double thechargefrommu = chargeoneachpmtfrommuon.at(pmtit);
				chargevsdisthist->Fill(thedistancefrommu,thechargefrommu);
				layercharges.at(layernum)+=thechargefrommu;
				numdigitslayers.at(layernum)++;
			}
		}
		for(int layerit=0; layerit<numpanels; layerit++){
			//cout<<"event "<<inputEntry<<" had "<<numdigitslayers.at(layerit)<<" digits on layer "<<layerit<<endl;
			if(muonstopsinMRD&&numMRDdigitsfrommu!=0){
				fracofchargedephist->Fill(layerit,(numdigitslayers.at(layerit)/numMRDdigitsfrommu));
			}
		}
		
	}  // end of loop over events
	
	
	// this has to go OUTSIDE the event loop.
	cout<<"doing efficiency calculation"<<endl;
	std::vector<double> paddleefficiency(nummrdpmts,0.);
	for(int pmtit=0; pmtit<nummrdpmts; pmtit++){
		if(numpaddlepassthrus.at(pmtit)>0) {
			paddleefficiency.at(pmtit) = (numdetectedpassthrus.at(pmtit)/numpaddlepassthrus.at(pmtit));
		} else {
			// this paddle had no tracks go through it, we can't determine efficiency.
			paddleefficiency.at(pmtit)=-1.;
		}
		
		int layernum=-1, paddlenum=-1;
		bool ishpaddle = ispaddlehoriz(pmtit, layernum, paddlenum);
		// ok, to actually do this we need an efficiency measurement PER POSITION PER PADDLE.
		// no way we have enough stats for that. We would also need to modify the numpaddlepassthrus and numdetectedpassthrus
		// above to be a vector to replace the int, with the vector representing number of (detected) hits 
		// within a given position bin. 
//		TH1D* thepaddleeffhist = efficiencygraphs.at(layer).at(paddle);
//		thepaddleeffhist.SetBinContent(thepaddleeffhist.FindBin(projectedxposition),paddleefficiency.at(pmtit));
	}
	
	
	/* Draw the histograms */
	//======================================================================================================
	//======================================================================================================
	TCanvas* canv1 =  new TCanvas("c1","c1");
	penetrationcmvsenergyall->Draw();
	TCanvas* canv2 = new TCanvas("c2","c2");
	penetrationhiststopall->Draw();
	TCanvas* canv3 = new TCanvas("c3","c3");
	penetrationhistpeneall->Draw();
	TCanvas* canv4 = new TCanvas("c4","c4");
	chargevsedephist->Draw(); 
	TCanvas* canv5 = new TCanvas("c5","c5");
	fracofchargedephist->Draw();
	TCanvas* canv6 = new TCanvas("c6","c6");
	chargevsdisthist->Draw();
	
	// write histograms to file
	outfile->cd();
	penetrationcmvsenergyall->Write();
	// the penetration vs energy plot including all angles
	penetrationhiststopall->Write();
	// the penetration vs energy plots for various angular ranges
	penetrationhiststop1->SetLineColor(kRed+1);
	penetrationhiststop2->SetLineColor(kMagenta);
	penetrationhiststop3->SetLineColor(kViolet-3);
	penetrationhiststop4->SetLineColor(kAzure-3);
	penetrationhiststop5->SetLineColor(kBlue+1);
	penetrationhiststop6->SetLineColor(8);
	penetrationhiststop7->SetLineColor(kOrange-2);
	penetrationhiststop1->Write();
	penetrationhiststop2->Write();
	penetrationhiststop3->Write();
	penetrationhiststop4->Write();
	penetrationhiststop5->Write();
	penetrationhiststop6->Write();
	penetrationhiststop7->Write();
	
	TCanvas* canv7 = new TCanvas("c7","c7");
	penetrationhiststop1->SetTitle("MRD Energy vs Muon Penetration for Varying Incident Angle");
	penetrationhiststop1->Draw("box");
	penetrationhiststop2->Draw("box, same");
	penetrationhiststop3->Draw("box, same");
	penetrationhiststop4->Draw("box, same");
	penetrationhiststop5->Draw("box, same");
	penetrationhiststop6->Draw("box, same");
	penetrationhiststop7->Draw("box, same");
	TLegend *leg = new TLegend(0.6798749,0.1278459,0.8925965,0.4220665,NULL,"brNDC");
	auto e1=leg->AddEntry(penetrationhiststop1,"0 - (1/5)Pi/2","l");
	auto e2=leg->AddEntry(penetrationhiststop2,"(1/5)Pi/2 - (2/5)Pi/2","l");
	auto e3=leg->AddEntry(penetrationhiststop3,"(2/5)Pi/2 - (3/5)Pi/2","l");
	auto e4=leg->AddEntry(penetrationhiststop4,"(3/5)Pi/2 - (3.5/5)Pi/2","l");
	auto e5=leg->AddEntry(penetrationhiststop5,"(3.5/5)Pi/2 - (4/5)Pi/2","l");
	auto e6=leg->AddEntry(penetrationhiststop6,"(4/5)Pi/2 - (4.5/5)Pi/2","l");
	auto e7=leg->AddEntry(penetrationhiststop7,"(4.5/5)Pi/2 - (5/5)Pi/2","l");
	leg->Draw();
	canv7->SaveAs("penetration2.png");
	
	tubeshithist->Write();
	
	// now measure the gradient and width of each histogram to get a (Penetration/Unit Energy) vs Angle with appropriate errors
	TF1* fitfunc;
	std::vector<double> penetrationconstvsangle;
	std::vector<double> penetrationgradvsangle;
	std::vector<double> penetrationerrorvsangle;
	// theangles and theangleerrors define the bin centers and widths in x direction
	Double_t* fitpars = new Double_t[3];
#ifdef VERBOSE
	cout<<"doing fit"<<endl;
	penetrationhiststop1->Fit("pol1");
	cout<<"getting fit function"<<endl;
#else
	penetrationhiststop1->Fit("pol1", "Q");
#endif
	fitfunc=(TF1*)penetrationhiststop1->GetFunction("pol1");
#ifdef VERBOSE
	cout<<"fit function is "<<fitfunc<<", getting fit params"<<endl;
#endif
	assert(fitfunc);
	fitfunc->GetParameters(fitpars);
#ifdef VERBOSE
	cout<<"fitpars is at "<<fitpars<<", putting into vectors"<<endl;
#endif
	penetrationconstvsangle.push_back(fitpars[0]);
	penetrationgradvsangle.push_back(fitpars[1]);
	penetrationerrorvsangle.push_back(fitfunc->GetParError(1));
	penetrationhiststop2->Fit("pol1", "Q");
	fitfunc=(TF1*)penetrationhiststop2->GetFunction("pol1");
	assert(fitfunc);
	fitfunc->GetParameters(fitpars);
	penetrationconstvsangle.push_back(fitpars[0]);
	penetrationgradvsangle.push_back(fitpars[1]);
	penetrationerrorvsangle.push_back(fitfunc->GetParError(1));
	penetrationhiststop3->Fit("pol1", "Q");
	fitfunc=(TF1*)penetrationhiststop3->GetFunction("pol1");
	fitfunc->GetParameters(fitpars);
	penetrationconstvsangle.push_back(fitpars[0]);
	penetrationgradvsangle.push_back(fitpars[1]);
	penetrationerrorvsangle.push_back(fitfunc->GetParError(1));
	penetrationhiststop4->Fit("pol1", "Q");
	fitfunc=(TF1*)penetrationhiststop4->GetFunction("pol1");
	fitfunc->GetParameters(fitpars);
	penetrationconstvsangle.push_back(fitpars[0]);
	penetrationgradvsangle.push_back(fitpars[1]);
	penetrationerrorvsangle.push_back(fitfunc->GetParError(1));
	penetrationhiststop5->Fit("pol1", "Q");
	fitfunc=(TF1*)penetrationhiststop5->GetFunction("pol1");
	fitfunc->GetParameters(fitpars);
	penetrationconstvsangle.push_back(fitpars[0]);
	penetrationgradvsangle.push_back(fitpars[1]);
	penetrationerrorvsangle.push_back(fitfunc->GetParError(1));
	penetrationhiststop6->Fit("pol1", "Q");
	fitfunc=(TF1*)penetrationhiststop6->GetFunction("pol1");
	fitfunc->GetParameters(fitpars);
	penetrationconstvsangle.push_back(fitpars[0]);
	penetrationgradvsangle.push_back(fitpars[1]);
	penetrationerrorvsangle.push_back(fitfunc->GetParError(1));
	penetrationhiststop7->Fit("pol1", "Q");
	fitfunc=(TF1*)penetrationhiststop7->GetFunction("pol1");
	fitfunc->GetParameters(fitpars);
	penetrationconstvsangle.push_back(fitpars[0]);
	penetrationgradvsangle.push_back(fitpars[1]);
	penetrationerrorvsangle.push_back(fitfunc->GetParError(1));
	
	// now make the graph of the results:
	TGraphErrors* penetrationgraph = new TGraphErrors(theangles.size(), &theangles[0], &penetrationgradvsangle[0], &theangleerrors[0], &penetrationerrorvsangle[0]);
	penetrationgraph->GetYaxis()->SetTitle("Rate of Energy Loss [MeV/cm]");
	penetrationgraph->GetXaxis()->SetTitle("Angle of Incidence [rads]");
	penetrationgraph->Write(); // Draw with option "PA" to draw the error bars and axes but no line
	TF1* f2 = new TF1("f2","[0]+[1]*x+[2]*x*x",0,1.6);
	f2->SetParameters(5,0.1,0.02);
	//f2->SetParameters(0.22,-0.1,-0.02); // need to set suitable initial parameters for fit to work
	gStyle->SetOptFit(0111); //use this to show the fit statisticson the graph.
	TFitResultPtr fitresult = penetrationgraph->Fit("f2"); // above found by trial & error in TBrowser
	fitresult->Write();
	
	TMatrixDSym cov = fitresult->GetCovarianceMatrix();  // covariance matrix for the TGraphErrors fit
	Double_t chi2   = fitresult->Chi2();                 // fit chi2
	Double_t par0   = fitresult->Value(0);               // fit offset
	Double_t err0   = fitresult->ParError(0);            // fit offset error
	Double_t par1   = fitresult->Value(1);               // fit linear component
	Double_t err1   = fitresult->ParError(1);            // fit linear component error
	Double_t par2   = fitresult->Value(2);               // fit square component
	Double_t err3   = fitresult->ParError(2);            // fit square component error
	fitresult->Print("V");                               // verbose print including covariance matrix
	
//	EXT PARAMETER                                   STEP         FIRST   
//	NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
//	p0           2.01418e-01   1.34709e-02   2.66739e-06  -5.47476e-02   
//	p1          -4.03848e-02   4.57988e-02   2.88394e-06  -9.53750e-02   
//	p2          -5.79169e-02   2.97912e-02   2.33052e-06  -1.36306e-01   
	
	
	TCanvas* canv8 = new TCanvas("c8","c8");
	chargevsedephist->Draw("colz");
	canv8->SaveAs("chargevsedep2.png");
	
	penetrationhistpeneall->Write();
	penetrationhistotherall->Write();
	chargevsedephist->Write();
	fracofchargedephist->Write();
	chargevsdisthist->Write();
	
	startposx->Write();
	startposy->Write();
	startposz->Write();
	startdirx->Write();
	startdiry->Write();
	startdirz->Write();
	startangleh->Write();
	startanglev->Write();
	muonanglehist->Write();
	
	/* Done, cleanup */
	//======================================================================================================
	//======================================================================================================
	gROOT->cd();
	// cleanup
	// =======
	cout<<"cleanup"<<endl;
	
	cout<<"resetting wcsimT branches"<<endl;
	if(wcsimT) wcsimT->ResetBranchAddresses();
	cout<<"closing wcsimfile"<<endl;
	if(wcsimfile) wcsimfile->Close();
	// should clean up wcsimT
	
	cout<<"closing output file"<<endl;
	outfile->Close();
	delete outfile;
	
}

std::map<int,int> makecopynummap(WCSimRootGeom* geo){
	// before we can look up the extent,layer and paddle num etc of a MRD paddle by it's WCSimRootGeom tube ID
	// we need to map the tube ID (arbitrarily the order in which ConstructGeometryTables found it) with the copyNo
	// (the order in which it was created). We should be able to do this by mapping their positions.
	std::map<int,int> tubeidtocopynum;
	double pmtradius = geo->GetMRDPMTRadius()/10.;  // returns in MM 
//	cout<<"MRD PMT radius is "<<pmtradius<<endl;
	for(int pmtit=0; pmtit<nummrdpmts; pmtit++){
		WCSimRootPMT thepmt = geo->GetMRDPMT(pmtit);
		double thexval = thepmt.GetPosition(0);
		double theyval = thepmt.GetPosition(1);
		double thezval = thepmt.GetPosition(2);
		// we have the position, but not the orientation... 
//		cout<<"MRD PMT "<<pmtit<<" has position ("<<thexval<<", "<<theyval<<", "<<thezval<<")"<<endl;
		
		double copyxval, copyyval, copyzval;
		bool foundcopy=false;
		int copynum;
		double totalerror;
		for(copynum=0; copynum<nummrdpmts; copynum++){
			copyzval = mrdcluster::paddle_originz.at(copynum)/10.;
			if(TMath::Abs(copyzval-thezval)<6.){   // we're in the correct layer - layer separation is 12cm.
				std::pair<Double_t,Double_t> xextents = mrdcluster::paddle_extentsx.at(copynum);
				std::pair<Double_t,Double_t> yextents = mrdcluster::paddle_extentsy.at(copynum);
				int copyorient = mrdcluster::paddle_orientations.at(copynum);
				if(copyorient==0) {  // 0 is a horizontal paddle
					copyyval = mrdcluster::paddle_originy.at(copynum)/10.;
					double copyxorig = mrdcluster::paddle_originx.at(copynum)/10.; // this is the PADDLE ORIGIN. 
					if(copyxorig>0) copyxval = (xextents.second)/10. + 41.868;      // RHS (largest x) for RH paddle (x>0)
					else copyxval = (xextents.first)/10. - 41.868;                 // the 41 accounts for end of taper+lightguide
				} else {             // paddle is vertical
					copyxval = mrdcluster::paddle_originx.at(copynum)/10.;
					double copyyorig = mrdcluster::paddle_originy.at(copynum)/10.;
					if(copyyorig>0) copyyval = (yextents.second)/10. + 41.868;
					else copyyval = (yextents.first)/10. - 41.868;
				}
				//cout<<"copynum "<<copynum<<" has position ("<<copyxval<<", "<<copyyval<<", "<<copyzval<<")"<<endl;
				// position matching is not _exact_, but generally within 1cm, Use a tolerance of the radius, for generosity
				totalerror = TMath::Sqrt(TMath::Power(copyyval-theyval,2)+TMath::Power(copyxval-thexval,2));
				if(totalerror < pmtradius){
					foundcopy=true;
					break;
				}
			}
		}
		
//		cout<<"matched copynum "<<copynum<<" has position ("<<copyxval<<", "<<copyyval<<", "<<copyzval<<")"<<endl;
//		cout<<"The error is ("<<copyxval-thexval<<", "<<copyyval-theyval<<", "<<copyzval-thezval<<") = "<<totalerror<<endl;
		if(!foundcopy){ cerr<<"Did not find copynum for MRD Tube ID "+to_string(pmtit)<<endl; assert(false); }
		tubeidtocopynum.emplace(pmtit,copynum);
		//cout<<pmtit<<" : "<<copynum<<endl;
	}
	
	return tubeidtocopynum;
}

bool ispaddlehoriz(int copyNo, int &panelnum, int &paddlenum){
	// which pair of panels
	int panelpairnum = floor(copyNo/(numpaddlesperpanelv+numpaddlesperpanelh));
	// copy num within a pair of panels
	int panelnumrem = copyNo - panelpairnum*(numpaddlesperpanelv+numpaddlesperpanelh);
	//int panelnum; passed by ref and set as a return value
	//int paddlenum; passed by ref and set as return value
	bool isvpaddle=false, ishpaddle=false;
	// first layer is a horizontal layer (leading horizontal layer removed)
	if(panelnumrem>(numpaddlesperpanelh-1)){
		panelnum = (panelpairnum*2) +1;
		isvpaddle = true;
		paddlenum = panelnumrem-numpaddlesperpanelh;
	} else {
		panelnum = (panelpairnum*2);
		ishpaddle = true;
		paddlenum = panelnumrem;
	}
	// scint paddles 0,1 are a vertical pair; X offset is the same for every pair
	int pairnum = floor(paddlenum/2);
	return ishpaddle;
}

void set_plot_style()
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}
