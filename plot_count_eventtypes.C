/* vim:set noexpandtab tabstop=4 wrap */
/*TODO: to modify for new source files:
0. Re-generate and replace in this directory libWCSimRoot.so, .rootmap, .pcm files.
1. Re-enable #include RootOptions.hh
*/
#ifndef VERBOSE
//#define VERBOSE
#endif
#ifndef WCSIMDEBUG
//#define WCSIMDEBUG
#endif
#ifndef MUTRACKDEBUG
//#define MUTRACKDEBUG 1
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

// genie headers
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepUtils.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Interaction/Interaction.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGLibrary.h"
#include "Conventions/Constants.h"
//#include "Numerical/GSFunc.h" 
// ensure in the evironment we have set
//ROOT_INCLUDE_PATH=${ROOT_INCLUDE_PATH}:${GENIE}/../include/GENIE
//ROOT_LIBRARY_PATH=${ROOT_LIBRARY_PATH}:${GENIE}/../lib

#include "genieinfo_struct.cxx"           // definition of a struct to hold genie info
// function to fill the into
void GetGenieEntryInfo(genie::EventRecord* gevtRec, genie::Interaction* genieint, GenieInfo& thegenieinfo);

const Float_t MRD_width = (305./2.);      // half width of steel in cm
const Float_t MRD_height = (274./2.);     // half height of steel in cm
const Float_t MRD_layer2 = 290.755;       // position in wcsim coords of second scint layer in cm
const Float_t MRD_start = 325.5;          // position in wcsim coord of MRD front face in cm
const Float_t MRD_depth = 139.09;         // total depth of the MRD in cm
/* output from WCSim:
########## MRD front face: 325.5                      ##########
########## MRD total Z length: 139.09                 ##########
########## MRD scintillator layer 0  (H) at z=266.535 ##########
########## MRD scintillator layer 1  (V) at z=278.645 ##########
########## MRD scintillator layer 2  (H) at z=290.755 ##########
########## MRD scintillator layer 3  (V) at z=302.865 ##########  // TODO: these layers are all too
########## MRD scintillator layer 4  (H) at z=314.975 ##########  // small by MRD_depth/2.
########## MRD scintillator layer 5  (V) at z=327.085 ##########  // now fixed in WCSim, but these nums
########## MRD scintillator layer 6  (H) at z=339.195 ##########  // need regenerating.
########## MRD scintillator layer 7  (V) at z=351.305 ##########  // Once done, update mrdscintlayers
########## MRD scintillator layer 8  (H) at z=363.415 ##########
########## MRD scintillator layer 9  (V) at z=375.525 ##########
########## MRD scintillator layer 10 (H) at z=387.635 ########## */
std::vector<double> mrdscintlayers{336.080, 348.190, 360.300, 372.410, 384.520, 396.630, 408.740, 420.850, 432.960, 445.070, 457.180 }; // this now includes the MRD_depth/2. offset, hence difference from above

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

const char* dirtpath="/pnfs/annie/persistent/users/moflaher/g4dirt";
const char* geniepath="/pnfs/annie/persistent/users/rhatcher/genie";
//const char* wcsimpath="/pnfs/annie/persistent/users/moflaher/wcsim";   // first set of 1M, various issues
const char* wcsimpath="/pnfs/annie/persistent/users/moflaher/wcsim_tankonly_03-05-17";
const char* wcsimlibrarypath="/annie/app/users/moflaher/wcsim/wcsim/libWCSimRoot.so";
const char* outpath="/annie/app/users/moflaher/wcsim/root_work";

const Bool_t printneutrinoevent=false;

void truthtracks(){
	
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
	
	TFile* dirtfile=0;
	std::string dirtfilename;
	TTree* tankflux=0;
	TTree* tankmeta=0;
	
	TFile* geniefile=0;
	TString geniefilepath;
	TTree* gtree=0;
	Int_t numgenietentries=0;
	
	// TChain for dirt files - this will be the main driver of the loop - all it's events will be processed.
	TChain* c =  new TChain("tankflux");
	TString chainpattern = TString::Format("%s/annie_tank_flux.*.root",dirtpath);	//1000->*
	cout<<"loading TChain entries from "<<chainpattern<<endl;
	c->Add(chainpattern);
	Int_t numents = c->GetEntries();
	cout<<"loaded "<<numents<<" entries in the chain"<<endl;
	if(numents<1){ return; }
	
	// tankflux
	Int_t genieentry=0;
	TBranch* genieentrybranch=0;
	Int_t ntankbranchval=0;
	TBranch *nTankBranch=0;
	Int_t* nuprimarybranchval=0;
	TBranch* nuprimaryBranch=0;
	Char_t vertexmaterial[100];
	TBranch* vertexmaterialbranch=0;
	
	// tankmeta
	TBranch* geniefilenamebranch=0;
	Char_t geniefilename[100];
	TBranch* potsbranch=0;
	Double_t pots;
	Double_t totalpots=0;       // count of POTs in all processed files
	
	// gtree
	cout<<"creating new genie::NtpMCEventRecord"<<endl;
	genie::NtpMCEventRecord* genierecordval = new genie::NtpMCEventRecord;
	TBranch* genierecordBranch=0;
	
	// wcsimT
	WCSimRootEvent* b=0, *m=0, *v=0;
	TBranch* bp=0, *mp=0, *vp=0;
	WCSimRootTrigger* atrigt=0, *atrigm=0, *atrigv=0;
	
	// geoT
	WCSimRootGeom* geo = 0; 
	
	// information from genie:
	GenieInfo thegenieinfo;
	
	// ===================================================================================================
	// ===================================================================================================
	// done declaring file: move to loading and processing
	gROOT->cd();
	cout<<"loading first tankflux tree from "<<chainpattern<<" tchain"<<endl;
	c->LoadTree(0);
	Int_t treeNumber = -1;
	tankflux = c->GetTree();
	Int_t thistreesentries = tankflux->GetEntries();
	cout<<thistreesentries<<" entries in the first tree"<<endl;
	
	/*
	1. Load next g4dirt entry
	2. Check if genie primary, and volume is in tank - if not, continue
	3. If so, load genie entry.
	4. Check if interaction is QE - if not, continue
	5. If so, load wcsim detector response. 
	6. Load tracks, look for primary mu track through the MRD, record interaction details. 
	*/
	
	std::map<std::string,bool> eventtypes;
	Double_t numfiducialevents=0.;
	Double_t numfiducialNCevents=0.;
	Double_t numfiducialCCevents=0.;
	Double_t numfiducialCCQEevents=0.;
	Double_t numfiducialeventsmrdentering=0.;
	Double_t numfiducialNCeventsmrdentering=0.;
	Double_t numfiducialCCeventsmrdentering=0.;
	Double_t numfiducialCCQEeventsmrdentering=0.;
	Double_t numfiducialeventsmrdstopping=0.;
	Double_t numfiducialNCeventsmrdstopping=0.;
	Double_t numfiducialCCeventsmrdstopping=0.;
	Double_t numfiducialCCQEeventsmrdstopping=0.;
	Double_t numfiducialeventsmrdpenetrating=0.;
	Double_t numfiducialNCeventsmrdpenetrating=0.;
	Double_t numfiducialCCeventsmrdpenetrating=0.;
	Double_t numfiducialCCQEeventsmrdpenetrating=0.;
	
	
	cout<<"looping over tchain entries"<<endl;
//	numents=1000;
	Int_t wcsimTentry;
	// since WCSim only propagated tank events, there is no longer a 1:1 mapping between event numbers
	// in dirt files and WCSim files. As long as the selection criterion for dirt events is the same here
	// as in WCSim's PrimaryGeneratorAction, we can select the dirt files, and then just pull the next 
	// wcsim entry
	for(Int_t inputEntry=0; inputEntry<numents; inputEntry++){
		/* 	1. Load next g4dirt entry */ 
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
			// this means we've switched file - need to load the new meta tree and genie tree.
			// first pull out the new file name
			tankflux = c->GetTree();
			dirtfile = tankflux->GetCurrentFile();
			dirtfilename=dirtfile->GetName();
			thistreesentries = tankflux->GetEntries();
			cout<<"tankflux has "<<thistreesentries<<" entries in this file"<<endl;
			
			// retrieve genie filename from the meta tree, and open the corresponding genie file
			tankmeta = (TTree*)dirtfile->Get("tankmeta");
			tankmeta->SetBranchAddress("inputFluxName", geniefilename, &geniefilenamebranch);
			geniefilenamebranch->GetEntry(0);
			geniefilepath = TString::Format("%s/%s",geniepath,geniefilename);
			cout<<"corresponding genie file is "<<geniefilepath<<", loading this file"<<endl;
			if(geniefile) geniefile->Close(); geniefile=0;
			geniefile = TFile::Open(geniefilepath);
			if(!geniefile){
				cerr<<"genie file "<<geniefilepath<<" doesn't exist!"<<endl; 
				inputEntry += thistreesentries;	// skip the loop iterator forward by all the entries in this file
				continue; 
			}
			gtree = (TTree*)geniefile->Get("gtree");
			if(!gtree){cerr<<"gtree doesn't exist!"<<endl; break; }
			numgenietentries = gtree->GetEntries();
			cout<<"gtree has "<<numgenietentries<<" entries in this file"<<endl;
			if(numgenietentries<1){cerr<<"gtree has no entries!"<<endl; break; }
			
			// use regexp to pull out the file number needed for identifying the corresponding wcsim file
			std::match_results<string::const_iterator> submatches;
			// filename is of the form "annie_tank_flux.####.root"
			// #### is input file num. Need this to match against genie/wcsim file names
			std::regex theexpression (".*/[^0-9]+\\.([0-9]+)\\.root");
			cout<<"matching regex for filename "<<dirtfilename<<endl;
			std::regex_match (dirtfilename, submatches, theexpression);
			std::string submatch = (std::string)submatches[0];	// match 0 is 'whole match' or smthg
			if(submatch==""){ cerr<<"unrecognised input file pattern: "<<dirtfilename<<endl; return; }
			submatch = (std::string)submatches[1];
			cout<<"extracted submatch is "<<submatch<<endl;
			int filenum = atoi(submatch.c_str());
			
			// use filenum to open the corresponding wcsim file
			wcsimfilepath = TString::Format("%s/wcsim_0.%d.root",wcsimpath,filenum);
			cout<<"corresponding wcsim file is "<<wcsimfilepath<<endl;
			if(wcsimfile) wcsimfile->Close(); wcsimfile=0;
			wcsimfile = TFile::Open(wcsimfilepath);
			if(!wcsimfile){
				cerr<<"wcsimfile "<<wcsimfilepath<<" doesn't exist!"<<endl; 
				inputEntry += thistreesentries;	// skip iterator forward by all the entries in this file
				continue; 
			}
			// load the geometry tree and grab the geometry if we haven't already
			if(geo==0){
				TTree* geotree = (TTree*)wcsimfile->Get("wcsimGeoT");
				if(geotree==0){ cerr<<"NO GEOMETRY IN FIRST FILE?"<<endl; assert(false); }
				geotree->SetBranchAddress("wcsimrootgeom", &geo);
				if (geotree->GetEntries() == 0) { cerr<<"geotree has no entries!"<<endl; exit(9); }
				geotree->GetEntry(0);
			}
			// load the next set of wcsim event info
			wcsimT = (TTree*)wcsimfile->Get("wcsimT");
			if(!wcsimT){cerr<<"wcsimT doesn't exist!"<<endl; break; }
			numwcsimentries = wcsimT->GetEntries();
			cout<<"wcsimT has "<<numwcsimentries<<" entries in this file"<<endl;
			if(numwcsimentries<1){cerr<<"wcsimT has no entries!"<<endl; break; }
			wcsimTentry=-1;
			
			/* Set the branch addresses for the new trees */
			// tankflux:
			// genie file entry number for each entry, to get the genie intx info
			c->SetBranchAddress("entry",&genieentry,&genieentrybranch);
			// number of primaries (so we can create appropriately sized array)
			c->SetBranchAddress("ntank",&ntankbranchval, &nTankBranch);
			// material of vertex - identify as 'TankWater' to pull only primaries in the tank
			c->SetBranchAddress("vtxmat",&vertexmaterial, &vertexmaterialbranch);
			// array of whether particle is a genie primary
			nuprimaryBranch=c->GetBranch("primary");
			// tankmeta:
			// POTs in this genie file, so we can normalise wcsim event frequency
			tankmeta->SetBranchAddress("inputTotalPOTs", &pots, &potsbranch);
			potsbranch->GetEntry(0);
			Double_t lasttotpots=totalpots;
			if(TMath::IsNaN(pots)==0){ totalpots += pots; }
			if(totalpots<lasttotpots){ 
				cerr<<"ERROR! TOTAL POTS CAME DOWN FROM "<<lasttotpots<<" TO "<<totalpots<<endl; return; 
			} else {
				cout<<"adding "<<pots<<" POTs to the running total, making "<<totalpots<<" POTs so far"<<endl;
			}
			
			// wcsimT:
			// wcsim trigger classes
			wcsimT->SetBranchAddress("wcsimrootevent",&b, &bp);
			wcsimT->SetBranchAddress("wcsimrootevent_mrd",&m, &mp);
			wcsimT->SetBranchAddress("wcsimrootevent_facc",&v, &vp);
			bp->SetAutoDelete(kTRUE);
			mp->SetAutoDelete(kTRUE);
			vp->SetAutoDelete(kTRUE);
			if(bp==0||mp==0||vp==0){ cout<<"branches are zombies!"<<endl; }
			
			// gtree:
			gtree->SetBranchAddress("gmcrec",&genierecordval,&genierecordBranch);
			
			treeNumber=nextTreeNumber;
		}
		
		/* 2. Check if genie primary, and volume is in tank - if not, continue */
		//====================================================================================================
		//====================================================================================================
#ifdef VERBOSE
		cout<<"processing inputEntry "<<inputEntry<<", localEntry "<<localEntry
		    <<"/"<<thistreesentries<<" in tree "<<treeNumber<<endl;
#endif
		nTankBranch->GetEntry(localEntry);
		vertexmaterialbranch->GetEntry(localEntry);
		if(strcmp(vertexmaterial,"TankWater")!=0){ cout<<"neutrino vtx not in tank"<<endl; continue; }
		if(nuprimarybranchval){delete[] nuprimarybranchval;}
		nuprimarybranchval = new Int_t[ntankbranchval];
		nuprimaryBranch->SetAddress(nuprimarybranchval);
		nuprimaryBranch->GetEntry(localEntry);
		
		Bool_t primariesinthisentry=false;
		for(int i=0;i<ntankbranchval;i++){
			if(nuprimarybranchval[i]==1){ primariesinthisentry=true; break; }
		}
		if(!primariesinthisentry){ cout<<"wcsim primaries not genie primaries"<<endl; continue; }
		// These selection criteria are the WCSim PrimaryGeneratorAction ones. Any event that passes here
		// will have created a WCSimT entry:
		wcsimTentry++;   // do this now in case we introduce any 'continue' statements later
		
		/* 3. If so, load genie entry. */
		//====================================================================================================
		//====================================================================================================
#ifdef VERBOSE
		cout<<"getting genie info"<<endl;
#endif
		if(localEntry>(numgenietentries-1)){ cerr<<"can't load localEntry "<<localEntry
								 <<" from "<<geniefilepath<<" gtree: not enough entries!"<<endl; continue; }
		genieentrybranch->GetEntry(localEntry);
		genierecordBranch->GetEntry(genieentry);
		genie::EventRecord* gevtRec = genierecordval->event;
		genie::Interaction* genieint = gevtRec->Summary();
		
		GetGenieEntryInfo(gevtRec, genieint, thegenieinfo);  // fill thegenieinfo struct with all the genie info
		
		eventtypes=thegenieinfo.eventtypes;
		
		/* 4.5 Do fiducial volume cut: */
		//====================================================================================================
		//====================================================================================================
		Bool_t isinfiducialvol=false;
		if( (TMath::Sqrt(TMath::Power(thegenieinfo.genie_x, 2) 
			+ TMath::Power(thegenieinfo.genie_z-tank_start-tank_radius,2)) < fidcutradius) && 
			(TMath::Abs(thegenieinfo.genie_y-tank_yoffset) < fidcuty) && 
			((thegenieinfo.genie_z-tank_start-tank_radius) < fidcutz) ){
			isinfiducialvol=true;
		}
		if(!isinfiducialvol) continue;
		
		/*5. primary vertex in the tank fiducial volume: load wcsim detector response. */
		//====================================================================================================
		//====================================================================================================
		// read only first subtrigger; delayed decay detector response is not interesting for primary FSL tracks
#ifdef VERBOSE
		cout<<"getting wcsim entry "<<wcsimTentry<<endl;
#endif
		if(wcsimTentry>(numwcsimentries-1)){ cout<<"can't load wcsimT entry "<<wcsimTentry
				<<" from "<<wcsimfilepath<<" wcsimT - not enough entries!"<<endl; continue; }
		wcsimT->GetEntry(wcsimTentry);
		atrigt = b->GetTrigger(0);
		
		Int_t numtracks = atrigt->GetNtrack();
#ifdef VERBOSE
		cout<<"wcsim event had "<<numtracks<<" truth tracks"<<endl;
#endif
		
		Int_t nummuontracksintank=0;
		Double_t muonenergy=0;        // we'll keep updating the info so it only records the highest energy muon
		bool muonentersMRD=false;
		bool muonstopsinMRD=false;
		bool muonrangesoutMRD=false;
		double mrdpenetrationcm=0.;
		int mrdpenetrationlayers=0;
		double mutracklengthinMRD=0.;
		
//		/////////////////////////////////////////////////////
//		int numprimarymuons=0;
//		int primarymutrackid=-1;
//		Double_t maxmuonenergy=0.;
//		Int_t mutrackindex=-1;
//		for(int track=0; track<numtracks; track++){
//			WCSimRootTrack* nextrack = (WCSimRootTrack*)atrigt->GetTracks()->At(track);
//			Int_t primarypdg = nextrack->GetIpnu();
//			if(nextrack->GetFlag()!=0) continue;
//			if(primarypdg!=13) continue;                // not a muon
//			if(nextrack->GetE()>maxmuonenergy){
//				maxmuonenergy=nextrack->GetE();
//				mutrackindex=track;
//			} else {
//				//continue;
//			}
//			if(nextrack->GetParenttype()!=0) continue;  // not a primary
//			numprimarymuons++;
//			primarymutrackid=track;
//		}
//		//assert((numprimarymuons<2)&&"MORE THAN ONE PRIMARY MUON");
//		 // this does not trigger - can assume only one primary muon.
//		//assert((mutrackindex==primarymutrackid)&&"HIGHEST ENERGY MUON!=PRIMARY MUON");
//		// this does trigger - we cannot assume the primary muon is the highest E muon
//		///////////////////////////////////////////////
		
		// now scan through the truth tracks, find the primary muon and save the wcsim info from it
		for(int track=0; track<numtracks; track++){
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
			
			// is it a primary?
			//if(nextrack->GetParenttype()!=0) continue;     // not a primary
			
			// does it start in the tank?
			Int_t primarystartvol = nextrack->GetStartvol();
			// temporary override as these weren't correctly set in wcsim <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			//TODO REMOVE THIS once 1M sample is regenerated
			if(nextrack->GetStart(2)<tank_start){
				primarystartvol = 20;                       // start depth is facc or hall
			} else if(nextrack->GetStart(2)>(tank_start+(2.*tank_radius))){
				primarystartvol = 30;                       // start depth is mrd or hall
			} else {
				primarystartvol = 10;                       // start depth is tank
			}
			
#ifdef VERBOSE
			cout<<"primarystartvol is "<<primarystartvol<<endl;
#endif
			if(primarystartvol!=10) continue;               // start volume is not the tank
			nummuontracksintank++;
			// skip if we already evaluated a muon of higher energy, otherwise process this muon track
			//if(nextrack->GetE()>muonenergy){ muonenergy=nextrack->GetE(); } else { continue; }
			
			// does it stop in the mrd, or pass completely through the MRD? 
			Int_t primarystopvol = nextrack->GetStopvol();
			// temporary override as these weren't correctly set in wcsim <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			//TODO REMOVE THIS once 1M generated. 
			// BUT need to think about how to check for 'range-out' mrd events. Maybe this is fine?
			if(nextrack->GetStop(2)<tank_start){
				primarystopvol = 20;                        // start depth is facc or hall
			} else if(nextrack->GetStop(2)>(tank_start+(2.*tank_radius))){
				primarystopvol = 30;                        // start depth is mrd or hall
			} else {
				primarystopvol = 10;                        // start depth is tank
			}
#ifdef VERBOSE
			cout<<"primarystopvol is "<<primarystopvol<<endl;
#endif
			
			TLorentzVector primarystartvertex(  nextrack->GetStart(0),
												nextrack->GetStart(1),
												nextrack->GetStart(2),
												thegenieinfo.genie_t);
			TLorentzVector primarystopvertex(   nextrack->GetStop(0),
												nextrack->GetStop(1),
												nextrack->GetStop(2),
												nextrack->GetTime());
			
			Float_t oppx = primarystopvertex.X() - primarystartvertex.X();
			Float_t adj = primarystopvertex.Z() - primarystartvertex.Z();
			Float_t avgtrackanglex = TMath::ATan(oppx/adj);
			Float_t oppy = primarystopvertex.Y() - primarystartvertex.Y();
			Float_t avgtrackangley = TMath::ATan(oppy/adj);
			
			TVector3 differencevector  = (primarystopvertex.Vect()-primarystartvertex.Vect());
			
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
				if(TMath::Abs(primarystopvertex.X())<MRD_width){
					trackZlengthbeforeMRDXexit=primarystopvertex.Z()-primarystartvertex.Z();
				} else {
					double trackXdistanceinMRD=0.;
					if(primarystopvertex.X()>0){
						trackXdistanceinMRD=MRD_width-primarystartvertex.X();
					} else {
						trackXdistanceinMRD=-MRD_width-primarystartvertex.X();
					}
					trackZlengthbeforeMRDXexit= trackXdistanceinMRD/TMath::Tan(avgtrackanglex);
				}
				Float_t trackZlengthbeforeMRDYexit;
				if(TMath::Abs(primarystopvertex.Y())<MRD_height){
					trackZlengthbeforeMRDYexit=primarystopvertex.Z()-primarystartvertex.Z();
				} else {
					double trackYdistanceinMRD=0.;
					if(primarystopvertex.Y()>0){
						trackYdistanceinMRD=MRD_height-primarystartvertex.Y();
					} else {
						trackYdistanceinMRD=-MRD_height-primarystartvertex.Y();
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
						if((TMath::Abs(primarystopvertex.X())<MRD_width)
							&&(TMath::Abs(primarystopvertex.Y())<MRD_height)) muonstopsinMRD=true;
					}
				}
			}
			
		//}  // end of loop over tracks << moved after additions: we'll record ALL muons, 
		//   not just highest energy
		
		// ==============================================================================================
		// check if we found any suitable muon tracks
		// ==============================================================================================
		
#ifdef VERBOSE
			cout<<"This event had "<<nummuontracksintank<<" muon tracks. Using the highest energy track."<<endl;
#endif
		if(nummuontracksintank>0){
			if(muonentersMRD) numfiducialeventsmrdentering++;
			if(eventtypes.at("IsWeakNC")&&muonentersMRD) numfiducialNCeventsmrdentering++;
			if(eventtypes.at("IsWeakCC")&&muonentersMRD) numfiducialCCeventsmrdentering++;
			if(eventtypes.at("IsWeakCC")&&eventtypes.at("IsQuasiElastic")&&muonentersMRD) numfiducialCCQEeventsmrdentering++;
			if(muonstopsinMRD) numfiducialeventsmrdstopping++;
			if(eventtypes.at("IsWeakNC")&&muonstopsinMRD) numfiducialNCeventsmrdstopping++;
			if(eventtypes.at("IsWeakCC")&&muonstopsinMRD) numfiducialCCeventsmrdstopping++;
			if(eventtypes.at("IsWeakCC")&&eventtypes.at("IsQuasiElastic")&&muonstopsinMRD) numfiducialCCQEeventsmrdstopping++;
			if(muonrangesoutMRD) numfiducialeventsmrdpenetrating++;
			if(eventtypes.at("IsWeakNC")&&muonrangesoutMRD) numfiducialNCeventsmrdpenetrating++;
			if(eventtypes.at("IsWeakCC")&&muonrangesoutMRD) numfiducialCCeventsmrdpenetrating++;
			if(eventtypes.at("IsWeakCC")&&eventtypes.at("IsQuasiElastic")&&muonrangesoutMRD) numfiducialCCQEeventsmrdpenetrating++;
		}
		} // moved end of loop over tracks here.
		numfiducialevents++;
		if(eventtypes.at("IsWeakNC")) numfiducialNCevents++;
		if(eventtypes.at("IsWeakCC")) numfiducialCCevents++;
		if(eventtypes.at("IsWeakCC")&&eventtypes.at("IsQuasiElastic")) numfiducialCCQEevents++;
		
	}  // end of loop over events
	
	//======================================================================================================
	//======================================================================================================
	
	Double_t numbeamspills = totalpots/(4.0 * TMath::Power(10.,12.));
	Double_t numbeamspillsperday = (24.*60.*60.*1000.)/133.3333;	// 24 hours in ms / 133.33 ms between spills
	Double_t numdays = numbeamspills/numbeamspillsperday;
	cout<<"Results based on "<<totalpots<<" POTs, or "<<numbeamspills<<" beam spills, or "<<numdays<<" days of data"<<endl;
	
	cout<<"There were "<<numfiducialevents<<" total fiducial events, "<<numfiducialNCevents<<" NC and "
	    <<numfiducialCCevents<<" CC, with "<<numfiducialCCQEevents<<" CCQE events."<<endl;
	cout<<"In percentages that's "<<((numfiducialNCevents/numfiducialevents)*100.)<<"% NC, "
	    <<((numfiducialCCevents/numfiducialevents)*100.)<<"% CC, "
	    <<" and "<<((numfiducialCCQEevents/numfiducialevents)*100.)<<"% CCQE, "<<endl;
	cout<<"Entering the MRD there were "<<numfiducialeventsmrdentering<<" events in total ("
	    <<((numfiducialeventsmrdentering/numfiducialCCevents)*100.)<<"% of fiducial CC events), of which "
	    <<numfiducialNCeventsmrdentering<<"("<<((numfiducialNCeventsmrdentering/numfiducialeventsmrdentering)*100.)<<"%) were NC, "
	    <<numfiducialCCeventsmrdentering<<"("<<((numfiducialCCeventsmrdentering/numfiducialeventsmrdentering)*100.)<<"%) were CC, "
	    <<"and "<<numfiducialCCQEeventsmrdentering<<"("<<((numfiducialCCQEeventsmrdentering/numfiducialeventsmrdentering)*100.)
	    <<"%) were CCQE"<<endl;
	
	cout<<"Stopping in the MRD there were "<<numfiducialeventsmrdstopping<<" events in total ("
	    <<((numfiducialeventsmrdstopping/numfiducialCCevents)*100.)<<"% of fiducial CC events), of which "
	    <<numfiducialNCeventsmrdstopping<<"("<<((numfiducialNCeventsmrdstopping/numfiducialeventsmrdstopping)*100.)<<"%) were NC, "
	    <<numfiducialCCeventsmrdstopping<<"("<<((numfiducialCCeventsmrdstopping/numfiducialeventsmrdstopping)*100.)<<"%) were CC, "
	    <<"and "<<numfiducialCCQEeventsmrdstopping<<"("<<((numfiducialCCQEeventsmrdstopping/numfiducialeventsmrdstopping)*100.)
	    <<"%) were CCQE"<<endl;
	
	cout<<"Fully penetrating the MRD there were "<<numfiducialeventsmrdpenetrating<<" events in total ("
	    <<((numfiducialeventsmrdpenetrating/numfiducialCCevents)*100.)<<"% of fiducial CC events), of which "
	    <<numfiducialNCeventsmrdpenetrating<<"("<<((numfiducialNCeventsmrdpenetrating/numfiducialeventsmrdpenetrating)*100.)<<"%) were NC, "
	    <<numfiducialCCeventsmrdpenetrating<<"("<<((numfiducialCCeventsmrdpenetrating/numfiducialeventsmrdpenetrating)*100.)<<"%) were CC, "
	    <<"and "<<numfiducialCCQEeventsmrdpenetrating<<"("<<((numfiducialCCQEeventsmrdpenetrating/numfiducialeventsmrdpenetrating)*100.)
	    <<"%) were CCQE"<<endl;
	
	double numfiducialCCeventsmrdsideexit = numfiducialCCeventsmrdentering-numfiducialCCeventsmrdstopping-numfiducialCCeventsmrdpenetrating;
	double numfiducialCCQEeventsmrdsideexit = numfiducialCCQEeventsmrdentering-numfiducialCCQEeventsmrdstopping-numfiducialCCQEeventsmrdpenetrating;
	
	cout<<"All results are fiducial:"<<endl;
	
	cout<<"There were "<<numfiducialCCevents<<" total CC events."<<endl;
	cout<<"There were "<<numfiducialCCeventsmrdentering<<" CC events that entered the MRD ("
	    <<((numfiducialCCeventsmrdentering/numfiducialCCevents)*100.)<<"% of all CC events)."<<endl
	    <<numfiducialCCeventsmrdstopping<<" CC events stopped in the MRD ("
	    <<((numfiducialCCeventsmrdstopping/numfiducialCCeventsmrdentering)*100.)<<"% of CC entering events, "
	    <<((numfiducialCCeventsmrdstopping/numfiducialCCevents)*100.)<<"% of all CC events)."<<endl
	    <<numfiducialCCeventsmrdpenetrating<<" CC events fully penetrated the MRD ("
	    <<((numfiducialCCeventsmrdpenetrating/numfiducialCCeventsmrdentering)*100.)<<"% of CC entering events, "
	    <<((numfiducialCCeventsmrdpenetrating/numfiducialCCevents)*100.)<<"% of all CC events)."<<endl
	    <<"That leaves "<<numfiducialCCeventsmrdsideexit<<" CC events that exited the sides of the MRD ("
	    <<((numfiducialCCeventsmrdsideexit/numfiducialCCeventsmrdentering)*100.)<<"% of CC entering events, "
	    <<((numfiducialCCeventsmrdsideexit/numfiducialCCevents)*100.)<<"% of all CC events)."<<endl;
	    
	cout<<"There were "<<numfiducialCCQEevents<<" total CCQE events."<<endl;
	cout<<"There were "<<numfiducialCCQEeventsmrdentering<<" CCQE events that entered the MRD ("
	    <<((numfiducialCCQEeventsmrdentering/numfiducialCCQEevents)*100.)<<"% of all CCQE events)."<<endl
	    <<numfiducialCCQEeventsmrdstopping<<" CCQE events stopped in the MRD ("
	    <<((numfiducialCCQEeventsmrdstopping/numfiducialCCQEeventsmrdentering)*100.)
	    <<"% of CCQE entering events, "
	    <<((numfiducialCCQEeventsmrdstopping/numfiducialCCQEevents)*100.)<<"% of all CCQE events)."<<endl
	    <<numfiducialCCQEeventsmrdpenetrating<<" CCQE events fully penetrated the MRD ("
	    <<((numfiducialCCQEeventsmrdpenetrating/numfiducialCCQEeventsmrdentering)*100.)
	    <<"% of CCQE entering events, "
	    <<((numfiducialCCQEeventsmrdpenetrating/numfiducialCCQEevents)*100.)<<"% of all CCQE events)."<<endl
	    <<"That leaves "<<numfiducialCCQEeventsmrdsideexit<<" CCQE events that exited the sides of the MRD ("
	    <<((numfiducialCCQEeventsmrdsideexit/numfiducialCCQEeventsmrdentering)*100.)
	    <<"% of CCQE entering events, "
	    <<((numfiducialCCQEeventsmrdsideexit/numfiducialCCQEevents)*100.)<<"% of all CCQE events)."<<endl;
	
	gROOT->cd();
	// cleanup
	// =======
	cout<<"cleanup"<<endl;
	if(c) c->ResetBranchAddresses();
	//cout<<"resetting tankflux branches"<<endl;
	if(tankflux) tankflux->ResetBranchAddresses();
	//cout<<"resetting tankmeta branches"<<endl;
	if(tankmeta) tankmeta->ResetBranchAddresses();
	cout<<"closing dirtfile"<<endl;
	if(dirtfile) dirtfile->Close(); // do we need to do this with a TChain?
	//cout<<"deleting dirtfile"<<endl;
	if(c) delete c; c=0;					// ??
	
	//cout<<"resetting gtree branches"<<endl;
	if(gtree) gtree->ResetBranchAddresses();
	cout<<"closing geniefile"<<endl;
	if(geniefile) geniefile->Close();
	// should clean up gtree
	//cout<<"deleting genierecordval"<<endl;
	if(genierecordval) delete genierecordval; genierecordval=0;
	
	//cout<<"resetting wcsimT branches"<<endl;
	if(wcsimT) wcsimT->ResetBranchAddresses();
	cout<<"closing wcsimfile"<<endl;
	if(wcsimfile) wcsimfile->Close();
	// should clean up wcsimT
	//cout<<"deleting nuprimarybranchval array"<<endl;
	if(nuprimarybranchval) delete[] nuprimarybranchval; nuprimarybranchval=0; // ? branch array
	
}

void GetGenieEntryInfo(genie::EventRecord* gevtRec, genie::Interaction* genieint, GenieInfo &thegenieinfo){
	// process information:
	/*TString*/ thegenieinfo.procinfostring = genieint->ProcInfo().AsString();
	/*TString*/ thegenieinfo.scatteringtypestring = genieint->ProcInfo().ScatteringTypeAsString();
	/*TString*/ thegenieinfo.interactiontypestring = genieint->ProcInfo().InteractionTypeAsString();
	thegenieinfo.eventtypes.at("IsQuasiElastic") = genieint->ProcInfo().IsQuasiElastic();
	thegenieinfo.eventtypes.at("IsResonant") = genieint->ProcInfo().IsResonant();
	thegenieinfo.eventtypes.at("IsDeepInelastic") = genieint->ProcInfo().IsDeepInelastic();
	thegenieinfo.eventtypes.at("IsCoherent") = genieint->ProcInfo().IsCoherent();
	thegenieinfo.eventtypes.at("IsDiffractive") = genieint->ProcInfo().IsDiffractive();
	thegenieinfo.eventtypes.at("IsInverseMuDecay") = genieint->ProcInfo().IsInverseMuDecay();
	thegenieinfo.eventtypes.at("IsIMDAnnihilation") = genieint->ProcInfo().IsIMDAnnihilation();
	thegenieinfo.eventtypes.at("IsSingleKaon") = genieint->ProcInfo().IsSingleKaon();
	thegenieinfo.eventtypes.at("IsNuElectronElastic") = genieint->ProcInfo().IsNuElectronElastic();
	thegenieinfo.eventtypes.at("IsEM") = genieint->ProcInfo().IsEM();
	thegenieinfo.eventtypes.at("IsWeakCC") = genieint->ProcInfo().IsWeakCC();
	thegenieinfo.eventtypes.at("IsWeakNC") = genieint->ProcInfo().IsWeakNC();
	thegenieinfo.eventtypes.at("IsMEC") = genieint->ProcInfo().IsMEC();
	/*Int_t*/ thegenieinfo.neutinteractioncode = genie::utils::ghep::NeutReactionCode(gevtRec);
	/*Int_t*/ thegenieinfo.nuanceinteractioncode  = genie::utils::ghep::NuanceReactionCode(gevtRec);
	/*TLorentzVector**/ thegenieinfo.genieVtx = gevtRec->Vertex();
	/*Double_t*/ thegenieinfo.genie_x = thegenieinfo.genieVtx->X() * 100.;         // same info as nuvtx in g4dirt file
	/*Double_t*/ thegenieinfo.genie_y = thegenieinfo.genieVtx->Y() * 100.;         // GENIE uses meters
	/*Double_t*/ thegenieinfo.genie_z = thegenieinfo.genieVtx->Z() * 100.;         // GENIE uses meters
	/*Double_t*/ thegenieinfo.genie_t = thegenieinfo.genieVtx->T() * 1000000000;   // GENIE uses seconds
	
	// neutrino information:
	/*Double_t*/ thegenieinfo.probeenergy = genieint->InitState().ProbeE(genie::kRfLab);  // GeV
	/*Int_t*/ thegenieinfo.probepdg = genieint->InitState().Probe()->PdgCode();
	/*TString*/ thegenieinfo.probepartname = genieint->InitState().Probe()->GetName();
	/*TLorentzVector**/ thegenieinfo.probemomentum = gevtRec->Probe()->P4();
	/*TVector3*/ thegenieinfo.probethreemomentum = thegenieinfo.probemomentum->Vect();
	/*TVector3*/ thegenieinfo.probemomentumdir = thegenieinfo.probethreemomentum.Unit();
	/*Double_t*/ thegenieinfo.probeanglex = 
		TMath::ATan(thegenieinfo.probethreemomentum.X()/thegenieinfo.probethreemomentum.Z());
	/*Double_t*/ thegenieinfo.probeangley = 
		TMath::ATan(thegenieinfo.probethreemomentum.Y()/thegenieinfo.probethreemomentum.Z());
	/*Double_t*/ thegenieinfo.probeangle = TMath::Max(thegenieinfo.probeanglex,thegenieinfo.probeangley);
	// n.b.  genieint->InitState().Probe != gevtRec->Probe()
	
	// target nucleon:
	/*genie::GHepParticle**/ thegenieinfo.targetnucleon = gevtRec->HitNucleon();
	/*int*/ thegenieinfo.targetnucleonpdg = genieint->InitState().Tgt().HitNucPdg();
	/*TString*/ thegenieinfo.targetnucleonname="";
	if ( genie::pdg::IsNeutronOrProton(thegenieinfo.targetnucleonpdg) ) {
		TParticlePDG * p = genie::PDGLibrary::Instance()->Find(thegenieinfo.targetnucleonpdg);
		thegenieinfo.targetnucleonname = p->GetName();
	} else {
		thegenieinfo.targetnucleonname = std::to_string(thegenieinfo.targetnucleonpdg);
	}
	/*TVector3*/ thegenieinfo.targetnucleonthreemomentum=TVector3(0.,0.,0.);
	/*Double_t*/ thegenieinfo.targetnucleonenergy=0.;
	if(thegenieinfo.targetnucleon){
		TLorentzVector* targetnucleonmomentum = thegenieinfo.targetnucleon->P4();
		thegenieinfo.targetnucleonthreemomentum = targetnucleonmomentum->Vect();
		thegenieinfo.targetnucleonenergy = targetnucleonmomentum->Energy(); //GeV
	}
	
	// target nucleus:
	/*Int_t*/ thegenieinfo.targetnucleuspdg = genieint->InitState().Tgt().Pdg();
	/*TParticlePDG**/ thegenieinfo.targetnucleus = 
		genie::PDGLibrary::Instance()->Find(thegenieinfo.targetnucleuspdg);
	/*TString*/ thegenieinfo.targetnucleusname = "unknown";
	if(thegenieinfo.targetnucleus){ thegenieinfo.targetnucleusname = thegenieinfo.targetnucleus->GetName(); }
	/*Int_t*/ thegenieinfo.targetnucleusZ = genieint->InitState().Tgt().Z();
	/*Int_t*/ thegenieinfo.targetnucleusA = genieint->InitState().Tgt().A();
	
	// remnant nucleus:
	int remnucpos = gevtRec->RemnantNucleusPosition(); 
	/*TString*/ thegenieinfo.remnantnucleusname="n/a";
	/*Double_t*/ thegenieinfo.remnantnucleusenergy=-1.;
	if(remnucpos>-1){
		thegenieinfo.remnantnucleusname = gevtRec->Particle(remnucpos)->Name();
		thegenieinfo.remnantnucleusenergy = gevtRec->Particle(remnucpos)->Energy(); //GeV
	}
	
	// final state lepton:
	int fsleppos = gevtRec->FinalStatePrimaryLeptonPosition();
	/*TString*/ thegenieinfo.fsleptonname="n/a";
	/*Double_t*/ thegenieinfo.fsleptonenergy=-1.;
	if(fsleppos>-1){
		thegenieinfo.fsleptonname = gevtRec->Particle(fsleppos)->Name();
		thegenieinfo.fsleptonenergy = gevtRec->Particle(fsleppos)->Energy();
	}
	
	// other remnants: TODO: this information is NOT being correctly read in
	/*Int_t*/ thegenieinfo.numfsprotons = genieint->ExclTag().NProtons();
	/*Int_t*/ thegenieinfo.numfsneutrons = genieint->ExclTag().NNeutrons();
	/*Int_t*/ thegenieinfo.numfspi0 = genieint->ExclTag().NPi0();
	/*Int_t*/ thegenieinfo.numfspiplus = genieint->ExclTag().NPiPlus();
	/*Int_t*/ thegenieinfo.numfspiminus = genieint->ExclTag().NPiMinus();
	
	// kinematic information
	Double_t NucleonM  = genie::constants::kNucleonMass; 
	// Calculate kinematic variables "as an experimentalist would measure them; 
	// neglecting fermi momentum and off-shellness of bound nucleons"
	/*TLorentzVector**/ thegenieinfo.k1 = gevtRec->Probe()->P4();
	/*TLorentzVector**/ thegenieinfo.k2 = gevtRec->FinalStatePrimaryLepton()->P4();
	/*Double_t*/ thegenieinfo.costhfsl = TMath::Cos( thegenieinfo.k2->Vect().Angle(thegenieinfo.k1->Vect()) );
	/*Double_t*/ thegenieinfo.fslanglegenie = thegenieinfo.k2->Vect().Angle(thegenieinfo.k1->Vect());
	// q=k1-k2, 4-p transfer
	/*TLorentzVector*/ thegenieinfo.q  = (*(thegenieinfo.k1))-(*(thegenieinfo.k2));
//	/*Double_t*/ thegenieinfo.Q2 = genieint->Kine().Q2();    // not set in our GENIE files!
	// momemtum transfer
	/*Double_t*/ thegenieinfo.Q2 = -1 * thegenieinfo.q.M2();
	// E transfer to the nucleus
	/*Double_t*/ thegenieinfo.Etransf  = (thegenieinfo.targetnucleon) ? thegenieinfo.q.Energy() : -1;
	// Bjorken x
	/*Double_t*/ thegenieinfo.x  = 
		(thegenieinfo.targetnucleon) ? 0.5*thegenieinfo.Q2/(NucleonM*thegenieinfo.Etransf) : -1;
	// Inelasticity, y = q*P1/k1*P1
	/*Double_t*/ thegenieinfo.y  = 
		(thegenieinfo.targetnucleon) ? thegenieinfo.Etransf/thegenieinfo.k1->Energy() : -1;
	// Hadronic Invariant mass ^ 2
	/*Double_t*/ thegenieinfo.W2 = 
	(thegenieinfo.targetnucleon) ? (NucleonM*NucleonM + 2*NucleonM*thegenieinfo.Etransf - thegenieinfo.Q2) : -1;
	
	if(printneutrinoevent){
		cout<<"This was a "<< thegenieinfo.procinfostring <<" interaction of a "
			<<thegenieinfo.probeenergy<<"GeV " << thegenieinfo.probepartname << " on a "; 
		
		if( thegenieinfo.targetnucleonpdg==2212 || thegenieinfo.targetnucleonpdg==2122 ){
			cout<<thegenieinfo.targetnucleonname<<" in a ";
		} else {
			cout<<"PDG-Code " << thegenieinfo.targetnucleonpdg<<" in a ";
		}
		
		if( thegenieinfo.targetnucleusname!="unknown"){ cout<<thegenieinfo.targetnucleusname<<" nucleus, "; }
		else { cout<<"Z=["<<thegenieinfo.targetnucleusZ<<","<<thegenieinfo.targetnucleusA<<"] nucleus, "; }
		
		if(remnucpos>-1){
			cout<<"producing a "<<thegenieinfo.remnantnucleusenergy<<"GeV "<<thegenieinfo.remnantnucleusname;
		} else { cout<<"with no remnant nucleus"; }  // DIS on 16O produces no remnant nucleus?!
		
		if(fsleppos>-1){
			cout<<" and a "<<thegenieinfo.fsleptonenergy<<"GeV "<<thegenieinfo.fsleptonname<<endl;
		} else{ cout<<" and no final state leptons"<<endl; }
		
		cout<<endl<<"Q^2 was "<<thegenieinfo.Q2<<"(GeV/c)^2, with final state lepton"
			<<" ejected at Cos(θ)="<<thegenieinfo.costhfsl<<endl;
		cout<<"Additional final state particles included "<<endl;
		cout<< "   N(p) = "       << thegenieinfo.numfsprotons
			<< "   N(n) = "       << thegenieinfo.numfsneutrons
			<< endl
			<< "   N(pi^0) = "    << thegenieinfo.numfspi0
			<< "   N(pi^+) = "    << thegenieinfo.numfspiplus
			<< "   N(pi^-) = "    << thegenieinfo.numfspiminus
			<<endl;
	}
}

