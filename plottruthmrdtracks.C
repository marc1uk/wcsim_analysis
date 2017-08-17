/* vim:set noexpandtab tabstop=4 wrap */
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// XXX XXX Version of WCSim used to generate the file XXX XXX
// XXX XXX      THIS MUST BE SET BEFORE CALLING       XXX XXX
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define FILE_VERSION 2
/*
Version 1:
wcsim_wdirt_07-02-17, 200 PMTs + 200 LAPPDs, including dirt intx.

Version 2:
wcsim_tankonly_03-05-17, 200 PMTs + 200 LAPPDs, tank only, with bug fixes, bad lappd pulse timing resoln

Version 3:
wcsim_tankonly_17-06-17, 120 PMTs of 3 different types (LUX, Watchboy, LBNE, 8inHQE), no LAPPDs.
*/

#ifndef VERBOSE
//#define VERBOSE
#endif
#ifndef WCSIMDEBUG
//#define WCSIMDEBUG
#endif
#ifndef MUTRACKDEBUG
//#define MUTRACKDEBUG
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
#pragma link C++ class std::vector<TVector3>+;
#endif
#include "TLegend.h"
#include "TText.h"
#include "TColor.h"
#include "TRandom3.h"
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
#if FILE_VERSION>2
#include "../wcsim/include/WCSimRootOptions.hh"
#endif

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

void MakePMTmap(WCSimRootGeom* geo, std::map<int, std::pair<int,int> > &topcappositionmap, std::map<int, std::pair<int,int> > &bottomcappositionmap, std::map<int, std::pair<int,int> > &wallpositionmap);
#include "makepmtmaps_standalone.cxx"     // definition of this function

void ColourPlotStyle();
void FillTankMapHist(WCSimRootGeom* geo, int tubeID, bool incone, std::map<std::string, TH2D*> &maphistos, double weight);
void ClearMapHistos(std::map<std::string,TH2D*> maphistos);	// clear the histograms

#include "genieinfo_struct.cxx"           // definition of a struct to hold genie info
// function to fill the into
void GetGenieEntryInfo(genie::EventRecord* gevtRec, genie::Interaction* genieint, GenieInfo& thegenieinfo);

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

float* Getlappdqpe();
static void SKIDigitizerThreshold(double& pe,int& iflag);
TRandom3* mrand = new TRandom3();         // needed to generate charge and other things

// not currently used but these should be stored in and retrieved from geo.
const Int_t numtankpmts=128+2*(26); // 26 pmts and lappds on each cap
const Int_t nummrdpmts=307;
const Int_t numvetopmts=26;
// these are used for making the pmt map.
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

#if FILE_VERSION==1
const char* wcsimpath="/pnfs/annie/persistent/users/moflaher/wcsim_wdirt_07-02-17";
#elif FILE_VERSION==2
const char* wcsimpath="/pnfs/annie/persistent/users/moflaher/wcsim_tankonly_03-05-17";
//const char* wcsimpath="/annie/app/users/moflaher/wcsim/build";
#elif FILE_VERSION==3
const char* wcsimpath="/pnfs/annie/persistent/users/moflaher/wcsim_tankonly_17-06-17";
#endif

const char* dirtpath="/pnfs/annie/persistent/users/moflaher/g4dirt";
const char* geniepath="/pnfs/annie/persistent/users/rhatcher/genie";
//const char* wcsimpath="/pnfs/annie/persistent/users/moflaher/wcsim";  // first 1M sample, various issues
//const char* wcsimpath="/annie/app/users/moflaher/wcsim/build";
const char* wcsimlibrarypath="/annie/app/users/moflaher/wcsim/wcsim/libWCSimRoot.so";
const char* outpath="/annie/app/users/moflaher/wcsim/root_work";
//const char* outpath="/annie/app/users/moflaher/wcsim/root_work/temp";

const Bool_t printneutrinoevent=false;

void truthtracks(){
	ColourPlotStyle();
	
	// load WCSim library for reading WCSim files
	cout<<"loading "<<wcsimlibrarypath<<endl;
	gSystem->Load(wcsimlibrarypath);
	// load genie for reading genie files - done in a separate macro, before calling ACLiC on this file
//	TString script_dir = gSystem->Getenv("GENIE");
//	script_dir += "/src/scripts/gcint/";
//	TString curr_dir = gSystem->pwd();
//	gSystem->cd(script_dir.Data());
//	gROOT->ProcessLine(".x loadincs.C");
//	gROOT->ProcessLine(".x loadlibs.C");
//	gSystem->cd(curr_dir.Data());
	
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
	
	TFile* lappdfile=0;
	TString lappdfilepath;
	TTree* lappdtree=0;
	Int_t numlappdentries=0;
	
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
	//Int_t* pdgbranchval=0;    // retrieve this info from wcsimT
	//TBranch* pdgBranch=0;
	
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
	
	// lappdtree
	// some of these are c-style arrays of all hits in the event, some are vectors of same.
	// for c-style arrays, avoid re-allocations of memory by just setting one array big enough for all events
	int lappd_evtnum;
	int lappd_numhitsthisevt;
	int LAPPDHITSMAX=1000;
	int lappd_hittile[LAPPDHITSMAX];            // "lappdhit_objnum"
	double lappd_hittilesposx[LAPPDHITSMAX];    // "lappdhit_x"
	double lappd_hittilesposy[LAPPDHITSMAX];    // could get this from geo, since we have the tile ID
	double lappd_hittilesposz[LAPPDHITSMAX];    // 
	double lappd_numphots[LAPPDHITSMAX];        // "lappdhit_edep" - # digits on lappd
	std::vector<double> lappd_hitcharge;        // "lappdhit_totalpes_perlappd2" - # hits on lappd
	std::vector<double> lappd_hitpeposx;        // position of the hit in the tile's ref frame [mm]
	std::vector<double> lappd_hitpeposy;        // "strip_coorx"
	std::vector<double> lappd_hitpeposz;
	std::vector<double> lappd_hittruetime;      // "strip_coort"
	// need pointers to set the branch addresses
	std::vector<double>* lappd_hitchargep=&lappd_hitcharge;
	std::vector<double>* lappd_hitpeposxp=&lappd_hitpeposx;
	std::vector<double>* lappd_hitpeposyp=&lappd_hitpeposy;
	std::vector<double>* lappd_hitpeposzp=&lappd_hitpeposz;
	std::vector<double>* lappd_hittruetimep=&lappd_hittruetime;
#if FILE_VERSION>3
	std::vector<double> lappd_hitglobalposx;    // hit position in global coordinates
	std::vector<double> lappd_hitglobalposy;
	std::vector<double> lappd_hitglobalposz;
	std::vector<double>* lappd_hitglobalposxp=&lappd_hitglobalposx;
	std::vector<double>* lappd_hitglobalposyp=&lappd_hitglobalposy;
	std::vector<double>* lappd_hitglobalposzp=&lappd_hitglobalposz;
#endif
	
	// not retrieved:
	//LAPPDtree->Branch("lappdhit_totalpes_perevt", &lappdhit_totalpes_perevt, "lappdhit_totalpes_perevt/I");
	//LAPPDtree->Branch("lappdhit_process",lappdhit_process,"lappdhit_process[lappd_numhits]/I");
	//LAPPDtree->Branch("lappdhit_particleID",lappdhit_particleID,"lappdhit_particleID[lappd_numhits]/I");
	//LAPPDtree->Branch("lappdhit_trackID",lappdhit_trackID,"lappdhit_trackID[lappd_numhits]/I");
	//LAPPDtree->Branch("lappdhit_stripnum", &lappdhit_stripnum); // hit strip
	//LAPPDtree->Branch("lappdhit_truetime2",&lappdhit_truetime2); // duplicate of strip_coort
	//LAPPDtree->Branch("lappdhit_smeartime2", &lappdhit_smeartime2); // smeared time. bad in current version.
	//LAPPDtree->Branch("lappdhit_primaryParentID2",&lappdhit_primaryParentID2);
	//LAPPDtree->Branch("lappdhit_NoOfneighstripsHit", &lappdhit_NoOfneighstripsHit); // size of below vectors
	//LAPPDtree->Branch("lappdhit_neighstripnum", &lappdhit_neighstripnum); // strip num
	//LAPPDtree->Branch("lappdhit_neighstrippeak", &lappdhit_neighstrippeak); // pulse amplitude
	//LAPPDtree->Branch("lappdhit_neighstrip_time", &lappdhit_neighstrip_time); // smeared hit time (0 of below)
	//LAPPDtree->Branch("lappdhit_neighstrip_lefttime", &lappdhit_neighstrip_lefttime); // rel. pulse ETA @ LHS
	//LAPPDtree->Branch("lappdhit_neighstrip_righttime", &lappdhit_neighstrip_righttime); // "    "   "   @ RHS
	
	// geoT
	WCSimRootGeom* geo = 0;
	int numpmts;
	int numlappds;
	std::vector<std::string> PMTNames;
	std::vector<Int_t> NumPMTsByType;
	
	// information from genie:
	GenieInfo thegenieinfo;
	
	// ==============================================================================================
	// ==============================================================================================
	// output file
	TFile* fileout = new TFile(TString::Format("%s/EventDistributions.root",outpath),"RECREATE");
	fileout->cd();
	TTree* treeout = new TTree("treeout","Tank Event Properties");
	
	// store the file and event num info for further lookup later if needed
	TString geniefilestring;
	TBranch* bGenieFileString = treeout->Branch("GenieFile",&geniefilestring);
	int genieeventnum=-1;
	TBranch* bGenieEventNum = treeout->Branch("GenieEventNum",&genieeventnum);
	TString dirtfilestring;
	TBranch* bDirtFileString = treeout->Branch("DirtFile",&dirtfilestring);
	int dirteventnum=-1;
	TBranch* bDirtEventNum = treeout->Branch("DirtEventNum",&dirteventnum);
	TString wcsimfilestring;
	TBranch* bWCSimFileString = treeout->Branch("WCSimFile",&wcsimfilestring);
	TString lappdfilestring;
	TBranch* bLAPPDFileString = treeout->Branch("LAPPDFile",&lappdfilestring);
	int wcsimeventnum=-1;
	TBranch* bWCSimEventNum = treeout->Branch("WCSimEventNum",&wcsimeventnum);
	// now information about the event
	std::map<std::string,bool> eventtypes;
	std::map<std::string,bool>* eventtypesp=&eventtypes;
	TBranch* bEventType = treeout->Branch("TypesMap",&eventtypesp);
	// using a map is a bad idea! It means you can't use the event type in tree->Draw calls!!! 
	// plus you need to load a CollectionProxy by putting the following code:
	//		#include <map>
	//		#ifdef __MAKECINT__
	//		#pragma link C++ class std::map<std::string,bool>+;
	//		#endif
	// into a file and calling '.L thefile.C+' before you can even read that tree with Tree/Branch->GetEntry!
	// so use the genie method, and save a bunch of bools
	bool IsQuasiElastic=false;
	TBranch* bIsQuasiElastic = treeout->Branch("IsQuasiElastic",&IsQuasiElastic);
	bool IsResonant=false;
	TBranch* bIsResonant = treeout->Branch("IsResonant",&IsResonant);
	bool IsDeepInelastic=false;
	TBranch* bIsDeepInelastic = treeout->Branch("IsDeepInelastic",&IsDeepInelastic);
	bool IsCoherent=false;
	TBranch* bIsCoherent = treeout->Branch("IsCoherent",&IsCoherent);
	bool IsDiffractive=false;
	TBranch* bIsDiffractive = treeout->Branch("IsDiffractive",&IsDiffractive);
	bool IsInverseMuDecay=false;
	TBranch* bIsInverseMuDecay = treeout->Branch("IsInverseMuDecay",&IsInverseMuDecay);
	bool IsIMDAnnihilation=false;
	TBranch* bIsIMDAnnihilation = treeout->Branch("IsIMDAnnihilation",&IsIMDAnnihilation);
	bool IsSingleKaon=false;
	TBranch* bIsSingleKaon = treeout->Branch("IsSingleKaon",&IsSingleKaon);
	bool IsNuElectronElastic=false;
	TBranch* bIsNuElectronElastic = treeout->Branch("IsNuElectronElastic",&IsNuElectronElastic);
	bool IsEM=false;
	TBranch* bIsEM = treeout->Branch("IsEM",&IsEM);
	bool IsWeakCC=false;
	TBranch* bIsWeakCC = treeout->Branch("IsWeakCC",&IsWeakCC);
	bool IsWeakNC=false;
	TBranch* bIsWeakNC = treeout->Branch("IsWeakNC",&IsWeakNC);
	bool IsMEC=false;
	TBranch* bIsMEC = treeout->Branch("IsMEC",&IsMEC);
	// ok, moving on
	bool isintank=false;
	TBranch* bInTank = treeout->Branch("NuVtxInTank",&isintank);
	bool isinfiducialvol=false;
	TBranch* bInFidVol = treeout->Branch("NuVtxInFidVol",&isinfiducialvol);
	double eventq2=0.;
	TBranch* bEventQ2 = treeout->Branch("EventQ2",&eventq2);
	double eventEnu=0.;
	TBranch* bEventEnu = treeout->Branch("NeutrinoEnergy",&eventEnu);
	int neutrinopdg=0;
	TBranch* bNeutrinoPdg = treeout->Branch("NeutrinoPDG",&neutrinopdg);
	bool hasmuon=false;
	TBranch* bEventHasMuon = treeout->Branch("EventHasMuon",&hasmuon);
	TVector3 mustartvtx(0,0,0);
	TBranch* bMuonStartVtx = treeout->Branch("MuonStartVertex",&mustartvtx);
	TVector3 mustopvtx(0,0,0);
	TBranch* bMuonStopVtx = treeout->Branch("MuonStopVertex",&mustopvtx);
	double muonangle=0.;
	TBranch* bMuonAngle = treeout->Branch("MuonAngle",&muonangle);
	double mustartE=0.;
	TBranch* bMuonStartE = treeout->Branch("MuonStartEnergy",&mustartE);
//	double muendE=0.;
//	TBranch* bMuonEndE = treeout->Branch("MuonEndEnergy",&muendE);  // information not saved in WCSim
	bool muonentersMRD=false;
	TBranch* bMuonEntersMRD = treeout->Branch("MuonEntersMRD",&muonentersMRD);
	bool muonstopsinMRD=false;
	TBranch* bMuonStopsInMRD = treeout->Branch("MuonStopsInMRD",&muonstopsinMRD);
	bool muonrangesoutMRD=false;
	TBranch* bMuonRangesOutMRD = treeout->Branch("MuonRangesOutMRD",&muonrangesoutMRD);
	double mrdpenetrationcm=0.;
	TBranch* bMuonMrdPenetrationInCm = treeout->Branch("MuonMrdPenetrationInCm",&mrdpenetrationcm);
	int mrdpenetrationlayers=0;
	TBranch* bMuonMrdPenetrationLayers = treeout->Branch("MuonMrdPenetrationLayers",&mrdpenetrationlayers);
	double mutracklengthinMRD=0.;
	TBranch* bMuonTrackLengthInMRD = treeout->Branch("MuonTrackLengthInMRD",&mutracklengthinMRD);
	double mutracklengthintank=0.;
	TBranch* bMuonTrackLengthInTank = treeout->Branch("MuonTrackLengthInTank",&mutracklengthintank);
	// TODO: add LAPPD hit info
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
	// save presence of other particles: how many events have additional pions, etc
	// store number of charged particle tracks in the event! and their energy, direction... 
	int numpizerotracks=-1;
	TBranch* bNumPiZeroTracks = treeout->Branch("NumPiZeroTracks",&numpizerotracks);
	int numpiplustracks=-1;
	TBranch* bNumPiPlusTracks = treeout->Branch("NumPiPlusTracks",&numpiplustracks);
	int numpiminustracks=-1;
	TBranch* bNumPiMinusTracks = treeout->Branch("NumPiMinusTracks",&numpiminustracks);
	int nummutracks=-1;
	TBranch* bNumMuTracks = treeout->Branch("NumMuTracks",&nummutracks);
	int numgammatracks=-1;
	TBranch* bNumGammaTracks = treeout->Branch("NumGammaTracks",&numgammatracks);
	int numneutrontracks=-1;
	TBranch* bNumNeutronTracks = treeout->Branch("NumNeutronTracks",&numneutrontracks);
	int numprotontracks=-1;
	TBranch* bNumProtonTracks = treeout->Branch("NumProtonTracks",&numprotontracks);
	// information from the MRD PMTs
	int numMRDdigits=0;
	TBranch* bTotalMrdDigits = treeout->Branch("TotalMrdDigits",&numMRDdigits);
	double totMRDcharge=0.;
	TBranch* bTotalMrdCharge = treeout->Branch("TotalMrdCharge",&totMRDcharge);
	int numMRDdigitsfrommu=0;
	TBranch* bMrdDigitsFromMuon = treeout->Branch("MrdDigitsFromMuon",&numMRDdigitsfrommu);
	double MRDchargefrommu=0.;
	TBranch* bMrdChargeFromMuon = treeout->Branch("MrdChargeFromMuon",&MRDchargefrommu);
	std::vector<int> mrdtubeshitbymu;
	TBranch* bMrdTubesHitByMuon = treeout->Branch("MrdTubesHitByMuon",&mrdtubeshitbymu);
	// store information about the other primary tracks in the event.
	// this may be useful for extracting events that have something else going on, and what
	std::vector<int> trackpdg;
	std::vector<int>* trackpdgp = &trackpdg;
	TBranch* bTrackPDG = treeout->Branch("TrackPDG",&trackpdgp);
	std::vector<TVector3> trackstartpos;
	std::vector<TVector3>* trackstartposp = &trackstartpos;
	TBranch* bTrackStartPos = treeout->Branch("TrackStartPos",&trackstartposp);
	std::vector<TVector3> trackstoppos;
	std::vector<TVector3>* trackstopposp = &trackstoppos;
	TBranch* bTrackStopPos = treeout->Branch("TrackStopPos",&trackstopposp);
	std::vector<TVector3> trackstartmom;
	std::vector<TVector3>* trackstartmomp = &trackstartmom;
	TBranch* bTrackStartMom = treeout->Branch("TrackStartMom",&trackstartmomp);
	std::vector<double> trackstartE;
	std::vector<double>* trackstartEp = &trackstartE;
	TBranch* bTrackStartE = treeout->Branch("TrackStartE",&trackstartEp);
	std::vector<int> trackstartvol;
	std::vector<int>* trackstartvolp = &trackstartvol;
	TBranch* bTrackStartVol = treeout->Branch("TrackStartVol",&trackstartvolp);
	std::vector<int> trackstopvol;
	std::vector<int>* trackstopvolp = &trackstopvol;
	TBranch* bTrackStopVol = treeout->Branch("TrackStopVol",&trackstopvolp);
	std::vector<int> trackparenttype;
	std::vector<int>* trackparenttypep = &trackparenttype;
	TBranch* bTrackParentType = treeout->Branch("TrackParentType",&trackparenttypep);
	std::vector<double> trackstoptime;
	std::vector<double>* trackstoptimep = &trackstoptime;
	TBranch* bTrackStopTime = treeout->Branch("TrackStopTime",&trackstoptimep);
	
	// pmt masking
	
	std::vector<TBranch*> thebranches{ bInTank, bEventType, bEventQ2, bEventEnu, bNeutrinoPdg, bInFidVol, bEventHasMuon, bMuonStartVtx, bMuonStopVtx, bMuonStartE, bMuonTrackLengthInTank, bMuonMrdPenetrationInCm, bMuonMrdPenetrationLayers, bMuonEntersMRD, bMuonStopsInMRD, bMuonRangesOutMRD, bMuonTrackLengthInMRD, bTankChargeFromMuon, bFractionOfMuonChargeInCone};
	int someit=0;
	bool haszombies=false;
	for(auto abranch : thebranches){
		if(abranch==0){ cout<<"branch "<<someit<<" is a zombie"<<endl; haszombies=true; }
		someit++;
	}
	assert(!haszombies&&"output file branches have zombies");
	
	// ==============================================================================================
	// ==============================================================================================
	// file for outputting true vertices and digits for tank reconstruction efforts
	gROOT->cd();
	TFile* flateventfileout = new TFile(TString::Format("%s/trueQEvertexinfo.root",outpath), "RECREATE");
	flateventfileout->cd();
	TLorentzVector filemuonstartvertex(0.,0.,0.,0.);
	TLorentzVector filemuonstopvertex(0.,0.,0.,0.);
	TVector3 filemuondirectionvector(0.,0.,0.);
	std::vector<ROOT::Math::XYZTVector>  filedigitvertices;
	std::vector<ROOT::Math::XYZTVector>* filedigitverticesp = &filedigitvertices;
//	std::vector<TLorentzVector> filedigitvertices;
//	std::vector<TLorentzVector>* filedigitverticesp = &filedigitvertices;
	std::vector<Double_t> filedigitQs;
	std::vector<Double_t>* filedigitQsp = &filedigitQs;
	std::vector<std::string> filedigitsensortypes;
	std::vector<std::string>* filedigitsensortypesp=&filedigitsensortypes;
	std::vector<Double_t> filedigittsmears;
	std::vector<Double_t>* filedigittsmearsp=&filedigittsmears;
	TTree* vertextreenocuts = new TTree("vertextreenocuts","All True Tank QE Events");
	TBranch* MuonStartBranch = vertextreenocuts->Branch("MuonStartVertex",&filemuonstartvertex);
	TBranch* MuonStopBranch = vertextreenocuts->Branch("MuonStopVertex", &filemuonstopvertex);
	TBranch* MuonDirectionBranch = vertextreenocuts->Branch("MuonDirection", &filemuondirectionvector);
	TBranch* DigitVertexBranch = vertextreenocuts->Branch("DigitVertices", "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >", &filedigitverticesp);
//	TBranch* DigitVertexBranch = vertextreenocuts->Branch("DigitVertices", "std::vector<TLorentzVector>", &filedigitverticesp);
	TBranch* DigitChargeBranch = vertextreenocuts->Branch("DigitCharges", &filedigitQsp);
	TBranch* DigitDetTypeBranch = vertextreenocuts->Branch("DigitWhichDet",&filedigitsensortypesp);
	TBranch* DigitSmearBranch = vertextreenocuts->Branch("DigitTimeSmear",&filedigittsmearsp);
#ifdef LAPPD_DEBUG
	std::vector<double> intileposx;
	std::vector<double> intileposy;
	std::vector<double> poserrx;
	std::vector<double> poserry;
	std::vector<double> poserrz;
	std::vector<int> tileorient;
	std::vector<int> octagonside;
	std::vector<double>* intileposxp=&intileposx;
	std::vector<double>* intileposyp=&intileposy;
#if FILE_VERSION>3
	std::vector<double>* poserrxp=&poserrx;
	std::vector<double>* poserryp=&poserry;
	std::vector<double>* poserrzp=&poserrz;
#endif
	std::vector<int>* tileorientp=&tileorient;
	std::vector<int>* octagonsidep=&octagonside;
	TBranch* LAPPD_intileposx = vertextreenocuts->Branch("LAPPD_intileposx",&intileposxp);
	TBranch* LAPPD_intileposy = vertextreenocuts->Branch("LAPPD_intileposy",&intileposyp);
#if FILE_VERSION>3
	TBranch* LAPPD_poserrx = vertextreenocuts->Branch("LAPPD_poserrx",&poserrxp);
	TBranch* LAPPD_poserry = vertextreenocuts->Branch("LAPPD_poserry",&poserryp);
	TBranch* LAPPD_poserrz = vertextreenocuts->Branch("LAPPD_poserrz",&poserrzp);
#endif
	TBranch* LAPPD_tileorient = vertextreenocuts->Branch("LAPPD_tileorient",&tileorientp);
	TBranch* LAPPD_octagonside = vertextreenocuts->Branch("LAPPD_octagonside",&octagonsidep);
#endif
	
	if(MuonStartBranch==0||MuonStopBranch==0||MuonDirectionBranch==0||DigitVertexBranch==0||DigitChargeBranch==0||DigitDetTypeBranch==0||DigitSmearBranch==0){ 
		cerr<<"branches are zombies argh!"<<endl; 
		cout<<"MuonStartBranch="<<MuonStartBranch<<endl
			<<"MuonStopBranch="<<MuonStopBranch<<endl
			<<"MuonDirectionBranch="<<MuonDirectionBranch<<endl
			<<"DigitVertexBranch="<<DigitVertexBranch<<endl
			<<"DigitChargeBranch="<<DigitChargeBranch<<endl
			<<"DigitDetTypeBranch="<<DigitDetTypeBranch<<endl
			<<"DigitSmearBranch="<<DigitSmearBranch<<endl;
		assert(false&&"branches are zombies argh!");
	}
	
	// Jingbo tree of fiducial neutrino events
	TTree* vertextreefiducialcut = new TTree("vertextreefiducialcut","True Tank QE Events in Fiducial Volume");
	TBranch* MuonStartBranchFid = vertextreefiducialcut->Branch("MuonStartVertex",&filemuonstartvertex);
	TBranch* MuonStopBranchFid = vertextreefiducialcut->Branch("MuonStopVertex", &filemuonstopvertex);
	TBranch* MuonDirectionBranchFid = vertextreefiducialcut->Branch("MuonDirection", &filemuondirectionvector);
	TBranch* DigitVertexBranchFid = vertextreefiducialcut->Branch("DigitVertices", "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >", &filedigitverticesp);
//	TBranch* DigitVertexBranchFid = vertextreefiducialcut->Branch("DigitVertices", "std::vector<TLorentzVector>", &filedigitverticesp);
	TBranch* DigitChargeBranchFid = vertextreefiducialcut->Branch("DigitCharges", &filedigitQsp);
	TBranch* DigitDetTypeBranchFid = vertextreefiducialcut->Branch("DigitWhichDet",&filedigitsensortypesp);
	TBranch* DigitSmearBranchFid = vertextreefiducialcut->Branch("DigitTimeSmear",&filedigittsmearsp);
	
	TTree* vertextreefiducialmrd = new TTree("vertextreefiducialmrd","True Tank QE Events in Fiducial Volume With Muon in MRD");
	TBranch* MuonStartBranchFidMRD = vertextreefiducialmrd->Branch("MuonStartVertex",&filemuonstartvertex);
	TBranch* MuonStopBranchFidMRD = vertextreefiducialmrd->Branch("MuonStopVertex", &filemuonstopvertex);
	TBranch* MuonDirectionBranchFidMRD = vertextreefiducialmrd->Branch("MuonDirection", &filemuondirectionvector);
	TBranch* DigitVertexBranchFidMRD = vertextreefiducialmrd->Branch("DigitVertices", "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >", &filedigitverticesp);
//	TBranch* DigitVertexBranchFidMRD = vertextreefiducialmrd->Branch("DigitVertices", "std::vector<TLorentzVector>", &filedigitverticesp);
	TBranch* DigitChargeBranchFidMRD = vertextreefiducialmrd->Branch("DigitCharges", &filedigitQsp);
	TBranch* DigitDetTypeBranchFidMRD = vertextreefiducialmrd->Branch("DigitWhichDet",&filedigitsensortypesp);
	TBranch* DigitSmearBranchFidMRD = vertextreefiducialmrd->Branch("DigitTimeSmear",&filedigittsmearsp);
	
	// some counters
	Double_t numneutrinoeventsintank=0.;
	Double_t numneutrinoeventsintankmrdstopped=0.;
	Double_t numCCQEneutrinoeventsintank=0.;
	Double_t numCCQEneutrinoeventsinfidvol=0.;
	Double_t numCCQEneutrinoeventsinfidvolmrd=0.;
	Double_t numCCQEneutrinoeventsinfidvolmrdstopped=0.;
	Double_t numCCQEneutrinoeventsmrdstopped=0;
	Double_t numNCneutrinoeventsintank=0;
	Double_t numNCneutrinoeventsintankmrd=0;
	Double_t numCCneutrinoeventsintank=0;
	Double_t numCCneutrinoeventsintankmrd=0;
	Double_t nummuontracksintank=0.;
	Double_t nummuontracksinfidvol=0.;
	Double_t nummuontracksintankpassedcut=0.;
	Double_t nummuontracksinmrd=0.;
	Double_t nummuontracksinfidvolmrd=0.;
	
	// ===================================================================================================
	// ===================================================================================================
	// histograms file
	gROOT->cd();
	TFile* histofileout = new TFile(TString::Format("%s/TruthHistos.root",outpath),"RECREATE");
	
	// genie histograms
	TH1D* incidentneutrinoenergiesall = new TH1D("incidentneutrinoenergiesall","Distribution of Probe Neutrino Energies;Energy (GeV);Num Events",100,0,0.);
	TH1D* incidentneutrinoenergiesaccepted = new TH1D("incidentneutrinoenergiesaccepted","Distribution of Accepted Probe Neutrino Energies;Energy (GeV);Num Events",100,0,0.);
	TH1D* fslanglesall = new TH1D("fslanglesall","Distribution of Final State Lepton Angles;Angle (rads);Num Events",100,0.,TMath::Pi());
	TH1D* fslanglesaccepted = new TH1D("fslanglesaccepted","Distribution of Accepted Final State Lepton Angles;Angle (rads);Num Events",100,0.,TMath::Pi());
	TH1D* fslenergiesall = new TH1D("fslenergiesall","Distribution of Final State Lepton Energies;Energy (GeV);Num Events",100,0.,3.);
	TH1D* fslenergiesaccepted = new TH1D("fslenergiesaccepted","Distribution of Accepted Final State Lepton Energies;Energy (GeV);Num Events",100,0.,3.);
	TH1D* eventq2all = new TH1D("eventq2all","Distribution of Event Q^2 Values;Q^2 (GeV/c)^2;Num Events",100,0.,3.);
	TH1D* eventq2accepted = new TH1D("eventq2accepted","Distribution of Accepted Event Q^2 Values;Q^2 (GeV/c)^2;Num Events",100,0.,3.);
	
	TH3D* neutrinovertex = new TH3D("neutrinovertex","Distribution of Neutrino Vertices in the tank",100,-tank_radius,tank_radius,100,tank_start,tank_start+(tank_radius*2),100,-tank_halfheight,tank_halfheight);
	TH3D* neutrinovertexQE = new TH3D("neutrinovertexQE","Distribution of QE Neutrino Vertices in the tank",100,-tank_radius,tank_radius,100,tank_start,tank_start+(tank_radius*2),100,-tank_halfheight,tank_halfheight);
	TH3D* neutrinovertexQEaccepted = new TH3D("neutrinovertexQEaccepted","Distribution of QE Neutrino Vertices with an Accepted MRD Track",100,-tank_radius,tank_radius,100,tank_start,tank_start+(tank_radius*2),100,-tank_halfheight,tank_halfheight);
//	TH3D* neutronstopvertex = new TH3D("neutronstopvertex","Distribution of Primary Neutron Stopping Vertices in the tank",100,-tank_radius,tank_radius,100,tank_start,tank_start+(tank_radius*2),100,-tank_halfheight,tank_halfheight);
//	TH3D* neutronstopvertexaccepted = new TH3D("neutronstopvertexaccepted","Distribution of Primary Neutron Stopping Vertices With an Accepted MRD Track",100,-tank_radius,tank_radius,100,tank_start,tank_start+(tank_radius*2),100,-tank_halfheight,tank_halfheight);
	
	// WCSim histograms
	TH1D* incidentneutrinoenergiesacceptedwcsim = new TH1D("incidentneutrinoenergiesacceptedwcsim","Distribution of Probe Neutrino Energies;Energy (GeV);Num Events",100,0,0.);
	TH1D* fslanglesacceptedwcsim = new TH1D("fslanglesacceptedwcsim","Distribution of Accepted Final State Lepton Angles;Angle (rads);Num Events",100,0.,TMath::Pi());
	TH1D* fslenergiesacceptedwcsim = new TH1D("fslenergiesacceptedwcsim","Distribution of Accepted Final State Lepton Energies;Energy (GeV);Num Events",100,0.,3.);
	TH1D* eventq2acceptedwcsim = new TH1D("eventq2acceptedwcsim","Distribution of Accepted Event Q^2 Values;Q^2 (GeV/c)^2;Num Events",100,0.,3.);
	
	// separate one for wcsim, from digit integration
	TH1D* muedepositionswcsim = new TH1D("muedepositionswcsim","Distribution of Muon Energy Depositions In Tank ;Energy (PMT Q);Num Events",100,0.,1);
	TH1D* muedepositionsfidcut = new TH1D("muedepositionsfidcut","Distribution of Muon Energy Depositions In Tank (Fiducial);Energy (PMT Q);Num Events",100,0.,1);
	
	// Fiducial cut versions
	// genie ones
	TH1D* incidentneutrinoenergiesallfidcut = new TH1D("incidentneutrinoenergiesallfidcut","Distribution of Probe Neutrino Energies Fiducial;Energy (GeV);Num Events",100,0,0.);
	TH1D* incidentneutrinoenergiesacceptedfidcut = new TH1D("incidentneutrinoenergiesacceptedfidcut","Distribution of Accepted Probe Neutrino Energies;Energy (GeV);Num Events",100,0,0.);
	TH1D* fslanglesallfidcut = new TH1D("fslanglesallfidcut","Distribution of Accepted Final State Lepton Angles;Angle (rads);Num Events",100,0.,TMath::Pi());
	TH1D* fslanglesacceptedfidcut = new TH1D("fslanglesacceptedfidcut","Distribution of Accepted Final State Lepton Angles;Angle (rads);Num Events",100,0.,TMath::Pi());
	TH1D* fslenergiesallfidcut = new TH1D("fslenergiesallfidcut","Distribution of Accepted Final State Lepton Energies;Energy (GeV);Num Events",100,0.,3.);
	TH1D* fslenergiesacceptedfidcut = new TH1D("fslenergiesacceptedfidcut","Distribution of Accepted Final State Lepton Energies;Energy (GeV);Num Events",100,0.,3.);
	TH1D* eventq2allfidcut = new TH1D("eventq2allfidcut","Distribution of Accepted Event Q^2 Values;Q^2 (GeV/c)^2;Num Events",100,0.,3.);
	TH1D* eventq2acceptedfidcut = new TH1D("eventq2acceptedfidcut","Distribution of Accepted Event Q^2 Values;Q^2 (GeV/c)^2;Num Events",100,0.,3.);
	
	// with reconstructed values
	TH1D* incidentneutrinoenergiesacceptedwcsimfidcut = new TH1D("incidentneutrinoenergiesacceptedwcsimfidcut","Distribution of Accepted Probe Neutrino Energies;Energy (GeV);Num Events",100,0,0.);
	TH1D* fslanglesacceptedwcsimfidcut = new TH1D("fslanglesacceptedwcsimfidcut","Distribution of Accepted Final State Lepton Angles;Angle (rads);Num Events",100,0.,TMath::Pi());
	TH1D* fslenergiesacceptedwcsimfidcut = new TH1D("fslenergiesacceptedwcsimfidcut","Distribution of Accepted Final State Lepton Energies;Energy (GeV);Num Events",100,0.,3.);
	TH1D* eventq2acceptedwcsimfidcut = new TH1D("eventq2acceptedwcsimfidcut","Distribution of Accepted Event Q^2 Values;Q^2 (GeV/c)^2;Num Events",100,0.,3.);
	
	
	// yet to work due to all parents being 0.
	TH1D* muedepositionsacceptedwcsim = new TH1D("muedepositionsacceptedwcsim","Distribution of Muon Energy Depositions In Tank, with MRD Selection;Energy (PMT Q);Num Events",100,0.,100);
	
//	// debugging:
	TH1D* fsltruetracklength = new TH1D("fsltruetracklength", "Distribution of True Track Lengths", 100, 0., 1500.);
	TH1D* fsltruetracklengthintank = new TH1D("fsltruetracklengthintank", "Distribution of True Track Lengths In Tank", 100, 0., 1500.);
#ifdef WCSIMDEBUG
	TH3D* tankstartvertices = new TH3D("tankstartvertices", "Distribution of Tank Starting Vertices", 100, -150.,150., 100, -150., 150., 100, -100., 800.);
	TH3D* vetostartvertices = new TH3D("vetostartvertices", "Distribution of Veto Starting Vertices", 100, -150.,150., 100, -150., 150., 100, -100., 800.);
	TH3D* mrdstartvertices = new TH3D("mrdstartvertices", "Distribution of MRD Starting Vertices", 100, -150.,150., 100, -150., 150., 100, -100., 800.);
	TH3D* tankstopvertices = new TH3D("tankstopvertices", "Distribution of Tank Stopping Vertices", 100, -150.,150., 100, -150., 150., 100, -100., 800.);
	TH3D* vetostopvertices = new TH3D("vetostopvertices", "Distribution of Veto Stopping Vertices", 100, -150.,150., 100, -150., 150., 100, -100., 800.);
	TH3D* mrdstopvertices = new TH3D("mrdstopvertices", "Distribution of MRD Stopping Vertices", 100, -150.,150., 100, -150., 150., 100, -100., 800.);
#endif
	
	// test the hypothesis that track length in water can be estimated from total light in tank
	TH2D* tracklengthvsmuonlight = new TH2D("tracklengthvsmuonlight", "Muon Track Length vs Total Light from Muon", 100, 0., 1500., 100, 0., 50.);
	
	// EFFECTS OF PIONS IN FINAL STATE
	// ===============================
	// record map of hits on the wall with both charge and time of the hits
	TH3D* chargemap_nopions = new TH3D("chargemap_nopions", "Charge Distribution for CC0pi events", pmtsperring+2,-1,pmtsperring+1,numpmtrings+2,-1,numpmtrings+1, 100, 0., 1400.);
	
	// Just to test the inside/outside cherenkov cone algorithm
	TH2D* inconehistowall = new TH2D("chargemap_incone_wall", "Charge Distribution Inside Cherenkov Cone (Wall)", pmtsperring+2,-1,pmtsperring+1,numpmtrings+2,-1,numpmtrings+1);
	TH2D* inconehistotop = new TH2D("chargemap_incone_top","Charge Distribution Inside Cherenkov Cone (Top Cap)",caparraysize+2,-1,caparraysize+1,caparraysize+2,-1,caparraysize+1);
	TH2D* inconehistobottom = new TH2D("chargemap_incone_bottom","Charge Distribution Inside Cherenkov Cone (Bottom Cap)",caparraysize+2,-1,caparraysize+1,caparraysize+2,-1,caparraysize+1);
	
	TH2D* outconehistowall = new TH2D("chargemap_outcone_wall", "Charge Distribution Outside Cherenkov Cone (Wall)", pmtsperring+2,-1,pmtsperring+1,numpmtrings+2,-1,numpmtrings+1);
	TH2D* outconehistotop = new TH2D("chargemap_outcone_top","Charge Distribution Outside Cherenkov Cone (Top Cap)",caparraysize+2,-1,caparraysize+1,caparraysize+2,-1,caparraysize+1);
	TH2D* outconehistobottom = new TH2D("chargemap_outcone_bottom","Charge Distribution Outside Cherenkov Cone (Bottom Cap)",caparraysize+2,-1,caparraysize+1,caparraysize+2,-1,caparraysize+1);
	
	std::map<std::string, TH2D*> maphistos;
	maphistos.emplace("inconehistowall",inconehistowall);
	maphistos.emplace("inconehistotop",inconehistotop);
	maphistos.emplace("inconehistobottom",inconehistobottom);
	maphistos.emplace("outconehistowall",outconehistowall);
	maphistos.emplace("outconehistotop",outconehistotop);
	maphistos.emplace("outconehistobottom",outconehistobottom);
	
	
	// CHECKING DIGIT CHARGES AND TIMES
	// ================================
	TH1D* digitsqpmthist = new TH1D("digitsqpmthist","Digit Charges for PMTs",100,0.,30.);
	TH1D* digitsqlappdhist = new TH1D("digitsqlappdhist","Digit Charges for LAPPDs",100,0.,30.);
	TH1D* digitstpmthist = new TH1D("digitstpmthist","Digit Times for PMTs",1000,800.,3000.);
	TH1D* digitstlappdhist = new TH1D("digitstlappdhist","Digit Times for LAPPDs",1000,800.,3000.);
	TH1D* digittsmearpmthist = new TH1D("digittsmearpmthist","Digit T Smearings for PMTs",100,0.,5.);
	TH1D* digittsmearlappdhist = new TH1D("digittsmearlappdhist","Digit T Smearings for LAPPDs",100,0.,0.1);
	TH2D* pmttimesmearvsqhist = new TH2D("pmttimesmearvsqhist","PMT Q vs T Smearing",1000,0.,3.,100,0.,100.);
	TH2D* lappdtimesmearvsqhist = new TH2D("lappdtimesmearvsqhist","LAPPD Q vs T Smearing",1000,0.,0.1,100,0.,100.);
	
	
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
	
	cout<<"looping over tchain entries"<<endl;
//	numents=10000;
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
				cerr<<"this genie file doesn't exist!"<<endl; 
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
				cerr<<"this wcsimfile doesn't exist!"<<endl; 
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
				MakePMTmap(geo, topcappositionmap, bottomcappositionmap, wallpositionmap);
				numpmts = geo->GetWCNumPMT();
				numlappds = geo->GetWCNumLAPPD();
#if FILE_VERSION>2
				//TODO: save these
				PMTNames = geo->GetPMTNames();
				NumPMTsByType = geo->GetPmtCounts();
#else
				PMTNames=std::vector<std::string>{"PMT8inch"};
				NumPMTsByType = std::vector<int>{numpmts};
#endif
			}
			// load the next set of wcsim event info
			wcsimT = (TTree*)wcsimfile->Get("wcsimT");
			if(!wcsimT){cerr<<"wcsimT doesn't exist!"<<endl; break; }
			numwcsimentries = wcsimT->GetEntries();
			cout<<"wcsimT has "<<numwcsimentries<<" entries in this file"<<endl;
			if(numwcsimentries<1){cerr<<"wcsimT has no entries!"<<endl; break; }
			wcsimTentry=-1;
			
			// use the filenum to open the corresponding lappd file
			lappdfilepath = TString::Format("%s/wcsim_lappd_0.%d.root",wcsimpath,filenum);
			cout<<"corresponding lappd file is "<<lappdfilepath<<endl;
			if(lappdfile) lappdfile->Close(); lappdfile=0;
			lappdfile = TFile::Open(lappdfilepath);
			if(!lappdfile){
				cerr<<"this lappdfile doesn't exist!"<<endl; 
				inputEntry += thistreesentries;	// skip iterator forward by all the entries in this file
				continue; 
			}
			//load the next set of lappd event info
			lappdtree = (TTree*)lappdfile->Get("LAPPDTree");
			if(!lappdtree){cerr<<"lappdtree doesn't exist!"<<endl; break; }
			numlappdentries = lappdtree->GetEntries();
			cout<<"lappdtree has "<<numlappdentries<<" entries in this file"<<endl;
			if(numlappdentries!=numwcsimentries){
				cerr<<"NUM LAPPD TREE ENTRIES != NUM WCSIMT ENTRIES!!!"<<endl;
				break;
			}
			if(numlappdentries<1){cerr<<"lappdtree has no entries!"<<endl; break;}
			
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
			if(bp==0||mp==0||vp==0){ cerr<<"branches are zombies!"<<endl; break; }
			
			// lappdtree:
			int branchesok=0;
			branchesok =lappdtree->SetBranchAddress("lappdevt",       &lappd_evtnum);
			if(branchesok<0) cerr<<"lappdevt branch error "<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappd_numhits",  &lappd_numhitsthisevt);
			if(branchesok<0) cerr<<"lappd_numhits="<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappdhit_stripcoorx",    &lappd_hitpeposxp);
			if(branchesok<0) cerr<<"lappd_stripcoorx="<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappdhit_stripcoory",    &lappd_hitpeposyp);
			if(branchesok<0) cerr<<"lappd_stripcoory="<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappdhit_stripcoorz",    &lappd_hitpeposzp);
			if(branchesok<0) cerr<<"lappd_stripcoorz="<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappdhit_stripcoort",    &lappd_hittruetimep);
			if(branchesok<0) cerr<<"lappd_stripcoort="<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappdhit_objnum",        &lappd_hittile);
			if(branchesok<0) cerr<<"lappd_objnum="<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappdhit_x",             &lappd_hittilesposx);
			if(branchesok<0) cerr<<"lappd_z="<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappdhit_y",             &lappd_hittilesposy);
			if(branchesok<0) cerr<<"lappd_y="<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappdhit_z",             &lappd_hittilesposz);
			if(branchesok<0) cerr<<"lappd_z="<<branchesok<<endl;
#if FILE_VERSION>3
			branchesok =lappdtree->SetBranchAddress("lappdhit_globalcoorx",   &lappd_hitglobalposxp);
			if(branchesok<0) cerr<<"lappd_hitglobalx="<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappdhit_globalcoory",   &lappd_hitglobalposyp);
			if(branchesok<0) cerr<<"lappd_hitglobaly="<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappdhit_globalcoorz",   &lappd_hitglobalposzp);
			if(branchesok<0) cerr<<"lappd_hitglobalz="<<branchesok<<endl;
#endif
			branchesok =lappdtree->SetBranchAddress("lappdhit_edep",          &lappd_numphots);
			if(branchesok<0) cerr<<"lappd_edep="<<branchesok<<endl;
			// lappdhit_edep is an c-style array of doubles of #true photon hits on each LAPPD in an event
			branchesok =lappdtree->SetBranchAddress("lappdhit_totalpes_perlappd2", &lappd_hitchargep);
			// lappdhit_totalpes_perlappd2 is the same thing in a vector, retrieved from the digits not SD hits
			// FIXME should store the charges
			if(branchesok<0) cerr<<"lappd_hitcharge="<<branchesok<<endl;
			
			// gtree:
			gtree->SetBranchAddress("gmcrec",&genierecordval,&genierecordBranch);
			
			treeNumber=nextTreeNumber;
		}
		
		geniefilestring=geniefilepath;
		dirtfilestring=TString(dirtfilename);
		wcsimfilestring=wcsimfilepath;
		lappdfilestring=lappdfilepath;
		dirteventnum=localEntry;
		
		/* 2. Check if genie primary, and volume is in tank - if not, continue */
		//====================================================================================================
		//====================================================================================================
#ifdef VERBOSE
		cout<<"processing inputEntry "<<inputEntry<<", localEntry "<<localEntry
		    <<"/"<<thistreesentries<<" in tree "<<treeNumber<<endl;
#endif
		nTankBranch->GetEntry(localEntry);
		vertexmaterialbranch->GetEntry(localEntry);
		if(strcmp(vertexmaterial,"TankWater")!=0){ /*cout<<"neutrino vtx not in tank"<<endl;*/ continue; }
		if(nuprimarybranchval){delete[] nuprimarybranchval;}
		nuprimarybranchval = new Int_t[ntankbranchval];
		nuprimaryBranch->SetAddress(nuprimarybranchval);
		nuprimaryBranch->GetEntry(localEntry);
		
		Bool_t primariesinthisentry=false;
		for(int i=0;i<ntankbranchval;i++){
			if(nuprimarybranchval[i]==1){ primariesinthisentry=true; break; }
		}
		if(!primariesinthisentry){ cout<<"wcsim primaries not genie primaries"<<endl; continue; }	// dirt recorded particles weren't genie primaries
		// These selection criteria are the WCSim PrimaryGeneratorAction ones. Any event that passes here
		// will have created a WCSimT entry:
		wcsimTentry++;   // do this now in case we introduce any 'continue' statements later
		numneutrinoeventsintank++;
		isintank=true;
		
		/* 3. If so, load genie entry. */
		//====================================================================================================
		//====================================================================================================
#ifdef VERBOSE
		cout<<"getting genie info"<<endl;
#endif
		if(localEntry>(numgenietentries-1)){ cout<<"can't load localEntry "<<localEntry
								 <<" from "<<geniefilepath<<" gtree: not enough entries!"<<endl; continue; }
		genieentrybranch->GetEntry(localEntry);
		genierecordBranch->GetEntry(genieentry);
		genie::EventRecord* gevtRec = genierecordval->event;
		genie::Interaction* genieint = gevtRec->Summary();
		
		//cout<<"scraping event info"<<endl;
		GetGenieEntryInfo(gevtRec, genieint, thegenieinfo);  // fill thegenieinfo struct with all the genie info
		//cout<<"done scraping info"<<endl;
		
		genieeventnum=genieentry;
		eventtypes=thegenieinfo.eventtypes;
		// then all the bools
		IsQuasiElastic=thegenieinfo.eventtypes.at("IsQuasiElastic");
		IsResonant=thegenieinfo.eventtypes.at("IsResonant");
		IsDeepInelastic=thegenieinfo.eventtypes.at("IsDeepInelastic");
		IsCoherent=thegenieinfo.eventtypes.at("IsCoherent");
		IsDiffractive=thegenieinfo.eventtypes.at("IsDiffractive");
		IsInverseMuDecay=thegenieinfo.eventtypes.at("IsInverseMuDecay");
		IsIMDAnnihilation=thegenieinfo.eventtypes.at("IsIMDAnnihilation");
		IsSingleKaon=thegenieinfo.eventtypes.at("IsSingleKaon");
		IsNuElectronElastic=thegenieinfo.eventtypes.at("IsNuElectronElastic");
		IsEM=thegenieinfo.eventtypes.at("IsEM");
		IsWeakCC=thegenieinfo.eventtypes.at("IsWeakCC");
		IsWeakNC=thegenieinfo.eventtypes.at("IsWeakNC");
		IsMEC=thegenieinfo.eventtypes.at("IsMEC");
		// ok, and the other info
		eventq2=thegenieinfo.Q2;
		eventEnu=thegenieinfo.probeenergy;
		neutrinopdg=thegenieinfo.probepdg;
		muonangle=thegenieinfo.fslanglegenie;
		if(eventtypes.at("IsWeakCC")) numCCneutrinoeventsintank++;
		else if(eventtypes.at("IsWeakNC")) numNCneutrinoeventsintank++;
		if(eventtypes.at("IsWeakCC")&&eventtypes.at("IsQuasiElastic")){ numCCQEneutrinoeventsintank++; }
		
		/* 4. Check if interaction is QE - if not, continue */
		//====================================================================================================
		//====================================================================================================
		//if(!(eventtypes.at("IsWeakCC")&&eventtypes.at("IsQuasiElastic"))) continue; 
		// disable for now: we'll check how many NC and CC-Other events have muons<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		
#ifdef VERBOSE
		cout<<"filling incident neutrino histograms"<<endl;
#endif
		incidentneutrinoenergiesall->Fill(thegenieinfo.probeenergy);
		//incidentneutrinoanglesall->Fill(thegenieinfo.probeangle);
		fslanglesall->Fill(thegenieinfo.fslanglegenie);
		fslenergiesall->Fill(thegenieinfo.fsleptonenergy);
		eventq2all->Fill(thegenieinfo.Q2);
		neutrinovertex->Fill(thegenieinfo.genie_x, thegenieinfo.genie_y, thegenieinfo.genie_z);
		if(eventtypes.at("IsWeakCC") && eventtypes.at("IsQuasiElastic"))
			neutrinovertexQE->Fill(thegenieinfo.genie_x, thegenieinfo.genie_y, thegenieinfo.genie_z);
		
		/* 4.5 Do fiducial volume cut: */
		//====================================================================================================
		//====================================================================================================
		isinfiducialvol=false;
		if( (TMath::Sqrt(TMath::Power(thegenieinfo.genie_x, 2) 
			+ TMath::Power(thegenieinfo.genie_z-tank_start-tank_radius,2)) < fidcutradius) && 
			(TMath::Abs(thegenieinfo.genie_y-tank_yoffset) < fidcuty) && 
			((thegenieinfo.genie_z-tank_start-tank_radius) < fidcutz) ){
			isinfiducialvol=true;
			incidentneutrinoenergiesallfidcut->Fill(thegenieinfo.probeenergy);
			fslanglesallfidcut->Fill(thegenieinfo.fslanglegenie);
			fslenergiesallfidcut->Fill(thegenieinfo.fsleptonenergy);
			eventq2allfidcut->Fill(thegenieinfo.Q2);
			numCCQEneutrinoeventsinfidvol++;
		}
		
		/*5. primary neutrino true QE vertex in the tank: load wcsim detector response. */
		//====================================================================================================
		//====================================================================================================
		// read only first subtrigger; delayed decay detector response is not interesting for primary FSL tracks
#ifdef VERBOSE
		cout<<"getting wcsim entry "<<wcsimTentry<<endl;
#endif
		if(wcsimTentry>(numwcsimentries-1)){ cout<<"can't load wcsimT entry "<<wcsimTentry
				<<" from "<<wcsimfilepath<<" wcsimT - not enough entries!"<<endl; break; }
		wcsimT->GetEntry(wcsimTentry);
		atrigt = b->GetTrigger(0);
		atrigm = m->GetTrigger(0);
		atrigv = v->GetTrigger(0);
		
		wcsimeventnum=wcsimTentry;
		int eventnumcheck=atrigt->GetHeader()->GetEvtNum();
		
		// Get LAPPD entries as well:
		// first we need to allocate necessary dynamic arrays
#ifdef VERBOSE
		cout<<"getting lappd entry"<<wcsimTentry<<endl;
#endif
		if(wcsimTentry>(numlappdentries-1)){ cout<<"can't load lappdtree entry "<<wcsimTentry
				<<" from "<<lappdfilepath<<" lappdtree - not enough entries!"<<endl; continue; }
		lappdtree->GetEntry(wcsimTentry);
		if(lappd_evtnum!=eventnumcheck){
			cerr<<"mismatch between lappd_evtnum="<<lappd_evtnum<<", and wcsim header event num, "
			<<"eventnumcheck="<<eventnumcheck<<". For ref: wcsimTentry="<<wcsimTentry<<endl;
		}
		
		// process WCSim truth tracks
		Int_t numtracks = atrigt->GetNtrack();
#ifdef VERBOSE
		cout<<"wcsim event had "<<numtracks<<" truth tracks"<<endl;
#endif
		
		std::vector<Float_t> neutrinoenergiesvector;
		std::vector<Float_t> primaryenergiesvector;
		std::vector<Double_t> scatteringanglesvector;
		std::vector<Int_t> acceptedtrackids;
		std::vector<Double_t> q2vector;
		std::vector<Double_t> muonenergydepositions;
		
		/* yet TODO - scan for neutron captures here? */  // 
		
		numpizerotracks=0;
		numpiplustracks=0;
		numpiminustracks=0;
		nummutracks=0;
		numgammatracks=0;
		numneutrontracks=0;
		numprotontracks=0;
		
		// now scan through the truth tracks, find the primary muon and save the wcsim info from it
		// UPDATED; to pull just the highest energy primary muon track, not all muons
		// (there are on average ~1.4 muons per event!) let's just do a quick scan first and pull
		// out only the highest energy track.
		Double_t maxmuonenergy=0.;
		Int_t mutrackindex=-1;
		trackpdg.clear();
		trackstartpos.clear();
		trackstoppos.clear();
		trackstartmom.clear();
		trackstartE.clear();
		trackstartvol.clear();
		trackstopvol.clear();
		trackparenttype.clear();
		trackstoptime.clear();
		for(int track=0; track<numtracks; track++){
			WCSimRootTrack* nextrack = (WCSimRootTrack*)atrigt->GetTracks()->At(track);
			Int_t primarypdg = nextrack->GetIpnu();
			switch (primarypdg){
				case 111: numpizerotracks++; break;
				case 211: numpiplustracks++; break;
				case -211: numpiminustracks++; break;
				case 13: nummutracks++; break;
				case 22: numgammatracks++; break;
				case 2112: numneutrontracks++; break;
				case 2212: numprotontracks++; break;
			}
			// record all tracks so we can tell when analysing later what else happened in this event
			trackpdg.push_back(nextrack->GetIpnu());
			TVector3 startpos(nextrack->GetStart(0),nextrack->GetStart(1),nextrack->GetStart(2));
			trackstartpos.push_back(startpos);
			TVector3 stoppos(nextrack->GetStop(0),nextrack->GetStop(1),nextrack->GetStop(2));
			trackstoppos.push_back(stoppos);
			TVector3 startmom(nextrack->GetPdir(0),nextrack->GetPdir(1),nextrack->GetPdir(2));
			trackstartmom.push_back(startmom);
			trackstartE.push_back(nextrack->GetE());
			trackstartvol.push_back(nextrack->GetStartvol());
			trackstopvol.push_back(nextrack->GetStopvol());
			trackparenttype.push_back(nextrack->GetParenttype());
			trackstoptime.push_back(nextrack->GetTime());
			
			if(TMath::Abs(primarypdg)!=13) continue;                // not a muon
			if(nextrack->GetParenttype()!=0) continue;              // not a primary
			if(nextrack->GetStartvol()!=10) continue;               // track doesn't start in tank
			Int_t primarystartvol;
			if(nextrack->GetStart(2)<tank_start){
				primarystartvol = 20;                               // start depth is facc or hall
			} else if(nextrack->GetStart(2)>(tank_start+(2.*tank_radius))){
				primarystartvol = 30;                               // start depth is mrd or hall
			} else {
				primarystartvol = 10;                               // start depth is tank
			}
			if(primarystartvol!=10){
				cerr<<"START VOLUME IS 10 BUT DEPTH ISN'T CORRECT FOR TANK!?"<<endl;
				cerr<<"start depth is "<<nextrack->GetStart(2)<<endl;
				cerr<<"tank starts at "<<(tank_start)<<", ends at "<<(tank_start+(2.*tank_radius))<<endl;
				//assert(false);
				continue;
			}
			nummuontracksintank++;
			if(nextrack->GetE()>maxmuonenergy){
				maxmuonenergy=nextrack->GetE();
				mutrackindex=track;
			} else {
				continue;
			}
		}
		if(mutrackindex<0) continue;                                  // there was no primary muon
		//for(int track=0; track<numtracks; track++){                 // disable loop over tracks
		for(int track=mutrackindex; track==mutrackindex; track++){    // just one iteration on highest E mu
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
			Float_t   GetTime()             trj->GetGlobalTime(); starting time of particle
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
			switch (primarypdg){
				case 111: numpizerotracks++; break;
				case 211: numpiplustracks++; break;
				case -211: numpiminustracks++; break;
				case 13: nummutracks++; break;
				case 22: numgammatracks++; break;
				case 2112: numneutrontracks++; break;
				case 2212: numprotontracks++; break;
			}
			if(TMath::Abs(primarypdg)!=13) continue;       // not a muon
			
			// for now we use truth information
			// is it a primary?
			Int_t primaryparentpdg = nextrack->GetParenttype();
			if(primaryparentpdg!=0) continue;
			
			// does it start in the tank?
			Int_t primarystartvol;
#if FILE_VERSION>1
			primarystartvol = nextrack->GetStartvol();
#else
			if(nextrack->GetStart(2)<tank_start){
				primarystartvol = 20;						// start depth is facc or hall
			} else if(nextrack->GetStart(2)>(tank_start+(2.*tank_radius))){
				primarystartvol = 30;						// start depth is mrd or hall
			} else {
				primarystartvol = 10;						// start depth is tank
			}
#endif
			
#ifdef VERBOSE
			cout<<"primarystartvol is "<<primarystartvol<<endl;
#endif
			if(primarystartvol!=10) continue;				// start volume is not the tank
			nummuontracksintank++;
			
			// does it stop in the mrd, or pass completely through the MRD? 
			Int_t primarystopvol;
#if FILE_VERSION>1
			primarystopvol = nextrack->GetStopvol();
#else
			// Do we need to think about 'range-out' mrd events. Maybe this is preferable?
			if(nextrack->GetStop(2)<tank_start){
				primarystopvol = 20;						// start depth is facc or hall
			} else if(nextrack->GetStop(2)>(tank_start+(2.*tank_radius))){
				primarystopvol = 30;						// start depth is mrd or hall
			} else {
				primarystopvol = 10;						// start depth is tank
			}
#endif
			
#ifdef VERBOSE
			cout<<"primarystopvol is "<<primarystopvol<<endl;
#endif
			
			TLorentzVector primarystartvertex(  nextrack->GetStart(0),
												nextrack->GetStart(1),
												nextrack->GetStart(2),
												nextrack->GetTime());
#if FILE_VERSION>2
			TLorentzVector primarystopvertex(   nextrack->GetStop(0),
												nextrack->GetStop(1),
												nextrack->GetStop(2),
												nextrack->GetStopTime());
#else
			TLorentzVector primarystopvertex(   nextrack->GetStop(0),
												nextrack->GetStop(1),
												nextrack->GetStop(2),
												-1); // not stored prior to this
#endif

#ifdef WCSIMDEBUG
			switch (primarystartvol){
			case 10:
				tankstartvertices->Fill(primarystartvertex.X(), primarystartvertex.Y(), primarystartvertex.Z());
				break;
			case 20:
				vetostartvertices->Fill(primarystartvertex.X(), primarystartvertex.Y(), primarystartvertex.Z());
				break;
			case 30:
				mrdstartvertices->Fill(primarystartvertex.X(), primarystartvertex.Y(), primarystartvertex.Z());
				break;
			}
			switch (primarystopvol){
			case 10:
				tankstopvertices->Fill(primarystopvertex.X(), primarystopvertex.Y(), primarystopvertex.Z());
				break;
			case 20:
				vetostopvertices->Fill(primarystopvertex.X(), primarystopvertex.Y(), primarystopvertex.Z());
				break;
			case 30:
				mrdstopvertices->Fill(primarystopvertex.X(), primarystopvertex.Y(), primarystopvertex.Z());
				break;
			}
#endif
			
			hasmuon=true;
			mustartvtx=primarystartvertex.Vect();
			mustopvtx=primarystopvertex.Vect();
			mustartE=nextrack->GetE();
			
			Float_t oppx = primarystopvertex.X() - primarystartvertex.X();
			Float_t adj = primarystopvertex.Z() - primarystartvertex.Z();
			Float_t avgtrackanglex = TMath::ATan(oppx/adj);
			Float_t oppy = primarystopvertex.Y() - primarystartvertex.Y();
			Float_t avgtrackangley = TMath::ATan(oppy/adj);
			
			TVector3 differencevector  = (primarystopvertex.Vect()-primarystartvertex.Vect());
			fsltruetracklength->Fill(differencevector.Mag());
			
			// continue if stopping volume is in either the tank or veto, or track is backward going.
			//if( primarystopvol==10 || primarystopvol==20 || primarystopvertex.Z() < primarystartvertex.Z() )
			//	continue;  
			// disabled for now, we'll record all muons regardless of penetration <<<<<<<<<<<<<<<<<<<<<<<<<<<<
			
			// ----------------------------------------------------------------------------------------------
			// calculate muon MRD penetration 
			// ----------------------------------------------------------------------------------------------
			
			// Check if the track penetrates enough layers of the MRD
			/////////////////////////////////////////////////////////
			// The mrd is as wide as the tank. We can ensure a track enters the mrd by projecting
			// the track forward from the start vertex at the angle between start and stop vertices,
			// and requiring:
			// 1) at z=MRD_start, x is between +/-(MRD_width/2);
			// 2) z endpoint is > MRD_start
			// For at least 2 layers of penetration, as above but with requirement on x @ z=MRD_layer2
			// For range-out, as above with requirement on x @ z=MRD_start+MRD_depth && 
			//   z endpoint is > MRD_start+MRD_depth
			
			Float_t xatmrd = primarystartvertex.X() + (MRD_layer2-primarystartvertex.Z())*TMath::Tan(avgtrackanglex);
			Float_t yatmrd = primarystartvertex.Y() + (MRD_layer2-primarystartvertex.Z())*TMath::Tan(avgtrackangley);
			
#ifdef WCSIMDEBUG
			cout<<"primary start vertex: ("<<primarystartvertex.X()<<", "<<primarystartvertex.Y()
				<<","<<primarystartvertex.Z()<<")"<<endl;
			cout<<"primary stop vertex: ("<<primarystopvertex.X()<<", "<<primarystopvertex.Y()
				<<","<<primarystopvertex.Z()<<")"<<endl;
			cout<<"oppx = "<<oppx<<endl;
			cout<<"adj = "<<adj<<endl;
			cout<<"angle = "<<avgtrackanglex<<endl;
			cout<<"tan(angle) = "<<TMath::Tan(avgtrackanglex)<<endl;
			cout<<"projected x at z="<<MRD_layer2<<" is "<<xatmrd<<endl;
			cout<<"xatmrd="<<xatmrd<<", MRD_width="<<MRD_width<<endl;
			cout<<"yatmrd="<<yatmrd<<", MRD_height="<<MRD_height<<endl;
#endif
			//if((TMath::Abs(xatmrd)>MRD_width)||(TMath::Abs(yatmrd)>MRD_height)) 
			//	continue;	// track does not meet MRD penetration requirement
			// if we got to here, we have a muon that stops in, or penetrates at least 2 MRD layers!
			
			
			// Alternatively, calculate the MRD penetration & place the cuts afterwards
			////////////////////////////////////////////////////////////////////////////
			
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
						muXexitpoint= (primarystopvertex.X()>0) ? MRD_width : -MRD_width;
						muYexitpoint=
							primarystartvertex.Y()+(trackZlengthbeforeMRDexit*TMath::Tan(avgtrackangley));
					} else {
						// track exits through the top or bottom of the MRD
						(primarystopvertex.Y()>0) ? muYexitpoint=MRD_height : muYexitpoint=-MRD_height;
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
			
			// ----------------------------------------------------------------------------------------------
			// retrieve any remaining information and record the event
			// ----------------------------------------------------------------------------------------------
			Float_t mu_rest_mass_E = 105.658;
			Float_t primaryenergy = (nextrack->GetE()-mu_rest_mass_E)/1000.;  // starting energy (GeV) (p^2+m^2)
			Float_t primarymomentummag = nextrack->GetP();                    // starting momentum
			TVector3 primarymomentumdir(nextrack->GetPdir(0),nextrack->GetPdir(1),nextrack->GetPdir(2));
			Float_t starttrackanglex = TMath::ATan(primarymomentumdir.X()/primarymomentumdir.Z());
			Float_t starttrackangley = TMath::ATan(primarymomentumdir.Y()/primarymomentumdir.Y());
			Float_t starttrackangle = TMath::Max(starttrackanglex,starttrackangley);
			//Double_t scatteringangle = thegenieinfo.k1->Angle(primarymomentumdir);
			Double_t scatteringangle = thegenieinfo.k1->Angle(differencevector);
			// scatteringangle for primaries by definition should be the same as fslanglegenie
			Double_t neutrinoenergyguess = thegenieinfo.probeenergy;	//TODO: how do we estimate this?
			TVector3 momtrans = (primarymomentummag*primarymomentumdir)-TVector3(0.,0.,neutrinoenergyguess);
			Double_t calculatedq2 = -1 * momtrans.Mag2();
			
			// Add this primary information to a vector. We _should_ only find one suitable primary
			// and can break once it is found, but... let's just see. Keep wcsim trackID in case.
			Int_t primarytrackid = nextrack->GetId();
			
#ifdef VERBOSE
			cout<<"found a suitable primary track"<<endl;
#endif
			neutrinoenergiesvector.push_back(neutrinoenergyguess);
			scatteringanglesvector.push_back(scatteringangle);
			primaryenergiesvector.push_back(primaryenergy);
			acceptedtrackids.push_back(primarytrackid);
			q2vector.push_back(calculatedq2);
			nummuontracksinmrd++;
			
			// ----------------------------------------------------------------------------------------------
			// calculate the track length in water
			// ----------------------------------------------------------------------------------------------
			// to calculate track length _in water_ find distance from start vertex to the point
			// where it intercepts the tank. if this length > total track length, use total track length
			// otherwise use this length. 
			
			// first find out the z value where the tank would leave the radius of the tank
#ifdef MUTRACKDEBUG
			cout<<"z0 = "<<thegenieinfo.genie_z-tank_start-tank_radius<<", x0 = "<<thegenieinfo.genie_x<<endl;
#endif
			Double_t xatziszero = 
			(thegenieinfo.genie_x - (thegenieinfo.genie_z-tank_start-tank_radius)*TMath::Tan(avgtrackanglex));
			Double_t firstterm = -TMath::Tan(avgtrackanglex)*xatziszero;
			Double_t thirdterm = 1+TMath::Power(TMath::Tan(avgtrackanglex),2.);
			Double_t secondterm = (TMath::Power(tank_radius,2.)*thirdterm) - TMath::Power(xatziszero,2.);
			Double_t solution1 = (firstterm + TMath::Sqrt(secondterm))/thirdterm;
			Double_t solution2 = (firstterm - TMath::Sqrt(secondterm))/thirdterm;
			Double_t tankendpointz;
			if(primarystopvertex.Z() > primarystartvertex.Z()){
				tankendpointz = solution1;	//forward going track
			} else {
				tankendpointz = solution2;	// backward going track
			}
			// correct for tank z offset
			tankendpointz += tank_start+tank_radius;
			Double_t tankendpointx = thegenieinfo.genie_x + (tankendpointz-thegenieinfo.genie_z)*TMath::Tan(avgtrackanglex);
			// now check if the particle would have exited through one of the caps before reaching this radius
			Double_t tankendpointy = 
			thegenieinfo.genie_y + (tankendpointz-thegenieinfo.genie_z)*TMath::Tan(avgtrackangley);
#ifdef MUTRACKDEBUG
			cout<<"tank start: "<<tank_start<<endl;
			cout<<"tank end: "<<(tank_start+2*tank_radius)<<endl;
			cout<<"tank radius: "<<tank_radius<<endl;
			cout<<"start dir =("<<primarymomentumdir.X()<<", "<<primarymomentumdir.Y()
				<<", "<<primarymomentumdir.Z()<<")"<<endl;
			cout<<"avgtrackanglex="<<avgtrackanglex<<endl;
			cout<<"avgtrackangley="<<avgtrackangley<<endl;
			cout<<"xatziszero="<<xatziszero<<endl;
			cout<<"firstterm="<<firstterm<<endl;
			cout<<"thirdterm="<<thirdterm<<endl;
			cout<<"secondterm="<<secondterm<<endl;
			cout<<"solution1="<<solution1<<endl;
			cout<<"solution2="<<solution2<<endl<<endl;
			
			cout<<"values before cap exit check"<<endl;
			cout<<"tankendpointz="<<tankendpointz<<endl;
			cout<<"tankendpointx="<<tankendpointx<<endl;
			cout<<"tankendpointy="<<tankendpointy<<endl;
#endif
			
			if(TMath::Abs(tankendpointy-tank_yoffset)>(tank_halfheight)){
				// this trajectory exits through the cap. Need to recalculate x, z exiting points...!
				if(primarystopvertex.Y()>primarystartvertex.Y()){
					tankendpointy = tank_halfheight+tank_yoffset;	// by definition of leaving through cap
				} else {
					tankendpointy = -tank_halfheight+tank_yoffset;
				}
				tankendpointz = 
				thegenieinfo.genie_z + (tankendpointy-thegenieinfo.genie_y)/TMath::Tan(avgtrackangley);
				tankendpointx = 
				thegenieinfo.genie_x + (tankendpointz-thegenieinfo.genie_z)*TMath::Tan(avgtrackanglex);
			} else {
				// this trajectory exited the tank by a side point; existing value is valid
			}
			Double_t maxtanktracklength = 
				TMath::Sqrt(TMath::Power(tank_radius*2.,2.)+TMath::Power(tank_halfheight*2.,2.));
#ifdef MUTRACKDEBUG
			cout<<"values after cap exit check"<<endl;
			cout<<"tankendpointz="<<tankendpointz<<endl;
			cout<<"tankendpointx="<<tankendpointx<<endl;
			cout<<"tankendpointy="<<tankendpointy<<endl;
			cout<<"max tank track length is "<<maxtanktracklength<<endl;
#endif
			
			// we're now able to determine muon track length in the tank:
			mutracklengthintank = TMath::Sqrt(
				TMath::Power((tankendpointx-thegenieinfo.genie_x),2)+
				TMath::Power((tankendpointy-thegenieinfo.genie_y),2)+
				TMath::Power((tankendpointz-thegenieinfo.genie_z),2) );
#ifdef MUTRACKDEBUG
			cout<<"muon tank track length: ("<<(tankendpointx-thegenieinfo.genie_x)
				<<", "<<(tankendpointy-thegenieinfo.genie_y)<<", "
				<<(tankendpointz-thegenieinfo.genie_z)<<") = "<<mutracklengthintank<<"cm total"<<endl;
			cout<<"muon tank exit point: ("<<tankendpointx<<", "<<tankendpointy<<", "<<tankendpointz<<") ";
			cout<<"muon start point : ("<<thegenieinfo.genie_x<<", "<<thegenieinfo.genie_y
				<<", "<<thegenieinfo.genie_z<<")"<<endl;
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
			if(mutracklengthintank>50){ nummuontracksintankpassedcut++; }
			
			
			// ----------------------------------------------------------------------------------------------
			// digit analysis
			// ----------------------------------------------------------------------------------------------
#ifdef VERBOSE
			cout<<"Analysing tank digits"<<endl;
#endif
			Int_t numdigitsthisevent = atrigt->GetCherenkovDigiHits()->GetEntries();
			//cout<<"this event has "<<numdigitsthisevent<<" digits"<<endl;
			filedigitvertices.clear();
			filedigitQs.clear();
			filedigittsmears.clear();
			filedigitsensortypes.clear();
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
			//TH3D hitpositions = TH3D("hitpositions","tiel",100,-tank_radius,tank_radius, 100,tank_start,tank_start+2*tank_radius,100,-tank_halfheight+tank_yoffset,tank_halfheight+tank_yoffset);
			//ClearMapHistos(maphistos);
			for(Int_t i=0; i<numdigitsthisevent; i++){
				// retrieve the digit information
				/////////////////////////////////
				WCSimRootCherenkovDigiHit* thedigihit = 
					(WCSimRootCherenkovDigiHit*)atrigt->GetCherenkovDigiHits()->At(i);
				int digitstubeid = thedigihit->GetTubeId()-1;
				double digitsq = thedigihit->GetQ();
				double digitst = thedigihit->GetT(); // this is time within the trigger window + 950ns
				WCSimRootEventHeader* trigheader=atrigt->GetHeader();
				double triggertime=trigheader->GetDate();
				double absolutedigitst=digitst-950.+triggertime;
#if FILE_VERSION<2
				// to be able to match digits to HitTimes to get true times or parents 
				// we need to correct the indices returned by GetPhotonIds with the offset of the
				// CherenkovHitTime entry matching the first CherenkovHit entry on the same PMT...!
				int timeArrayOffset=-1;
				int ncherenkovhits=atrigt->GetNcherenkovhits();
				for(int ipmt = 0; ipmt < ncherenkovhits; ipmt++) {
					WCSimRootCherenkovHit* hitobject = 
						(WCSimRootCherenkovHit*)atrigt->GetCherenkovHits()->At(ipmt);
					int tubeNumber = hitobject->GetTubeID()-1;
					if(tubeNumber==digitstubeid) {
						timeArrayOffset = hitobject->GetTotalPe(0);
						break;
					}
				}
#endif
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
				
				// if we want photon times; do this
				////////////////////////////////////
				/*
				std::vector<double> thetruetimes;
				double firstphotontime = 0.;
				for(int truephoton=0; truephoton<truephotonindices.size(); truephoton++){
					int thephotonsid = truephotonindices.at(truephoton);
#if FILE_VERSION<2
					thephotonsid+=timeArrayOffset;  // this is where we need to add in the offset!
#endif
					WCSimRootCherenkovHitTime *thehittimeobject = 
						(WCSimRootCherenkovHitTime*)atrigt->GetCherenkovHitTimes()->At(thephotonsid);
					double ahittime = thehittimeobject->GetTruetime();
					thetruetimes.push_back(ahittime);
				}
				firstphotontime = thetruetimes.at(0);
				*/
				
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
				// and see if any of them are this muon
				//cout<<endl<<"calculating muon energy deposition, mu track id "<<nextrack->GetId()<<endl;
				int numphotonsfromthismuon=0;
				//cout<<"   digit "<<i<<" has "<<truephotonindices.size()<<" true photons"<<endl;
				for(int truephoton=0; truephoton<truephotonindices.size(); truephoton++){
					int thephotonsid = truephotonindices.at(truephoton);
#if FILE_VERSION<2
					thephotonsid+=timeArrayOffset;
#endif
					WCSimRootCherenkovHitTime *thehittimeobject = 
						(WCSimRootCherenkovHitTime*)atrigt->GetCherenkovHitTimes()->At(thephotonsid);
					Int_t thephotonsparenttrackid = thehittimeobject->GetParentID();
					//cout<<"      HitTimeIndex "<<thephotonsid<<", HitTimeObject "<<thehittimeobject
					//	<<", HitTimeParent "<<thephotonsparenttrackid<<endl;
					if(thephotonsparenttrackid==nextrack->GetId()) {
						numphotonsfromthismuon++;
						//cout<<"################ FOUND A DIGIT ############"<<endl;
						//cout<<"the digits Q is "<<thedigihit->GetQ()<<endl;
						//return;
					}
				}
				// to estimate the charge from the muon we should scale each digit's total charge
				// by the fraction of hits in the digit which come from the muon
				if(numphotonsfromthismuon!=0){   //this digit had contribution from the muon!
					//cout<<"muon digit "<<i<<" had charge "<< digitsq << " and " << truephotonindices.size()
					//	<<" true photons, of which "<<numphotonsfromthismuon <<" were from this muon"<<endl
					tankchargefrommuon += 
						digitsq * ((double)numphotonsfromthismuon/(double)truephotonindices.size());
					numtankdigitsfrommuon++;
					if(std::find(tanktubeshitbymu.begin(), tanktubeshitbymu.end(),
						 digitstubeid)==tanktubeshitbymu.end())
						tanktubeshitbymu.push_back(digitstubeid);
				}
				
				// add the digit info for simplified file format
				/////////////////////////////////////////////////
				//TLorentzVector* adigitvector = new TLorentzVector(digitsx,digitsy,digitsz,absolutedigitst);
				ROOT::Math::XYZTVector adigitvector = ROOT::Math::XYZTVector(digitsx,digitsy,digitsz,absolutedigitst);
				filedigitvertices.push_back(adigitvector);
//				TLorentzVector adigitvector = TLorentzVector(digitsx,digitsy,digitsz,absolutedigitst);
//				filedigitvertices.push_back(adigitvector);
				filedigitQs.push_back(digitsq);
				std::string digitspst;
#if FILE_VERSION<3
				digitspst = "PMT8inch";
#else
				int tubetypeindex = geo->GetTubeIndex(digitstubeid);
				digitspst = geo->GetWCPMTNameAt(tubetypeindex);
#endif
				filedigitsensortypes.push_back(digitspst);
				double timingConstant;
				double timingResConstant;
				double timingResMinimum;
				double timingResPower;
				// WCSimPMTObject.cc:255
				// float timingResolution = timingResConstant + sqrt(timingConstant/Q);
				// if (timingResolution < 0.58) timingResolution=0.58;
				// float Smearing_factor = G4RandGauss::shoot(0.0,timingResolution);
				// --> HitTimeSmearing = returned as Smearing_factor
				// WCSimWCPMT.cc:173
				// time_PMT = time_true + PMT->HitTimeSmearing(Q);
				if(digitspst=="PMT8inch"){                    // SK 8inch
						timingConstant=1.890;
						timingResConstant=0.33;
						timingResMinimum=0.58;
				}
				else if(digitspst=="R7081"){                  // 10" LUX/Watchboy
						timingConstant = 1.890;
						timingResConstant=0.33;
						timingResMinimum=0.58;                // TTS: 3.2ns;
				}
				else if(digitspst=="D784KFLB"){               // 11" HQE LBNE
						timingConstant = 1.890;
						timingResConstant=0.33;
						timingResMinimum=0.58;                // TTS: 1.98ns;
				}
				else if(digitspst=="R5912HQE"){               // 8"  HQE new
						timingConstant = 1.890;
						timingResConstant=0.33;
						timingResMinimum=0.58;                // TTS: 2.4ns;
				}
				else if(digitspst=="FlatFacedPMT2inch"){      // Veto / MRD PMT ver 1
						timingConstant = 1.890;
						timingResConstant=0.33;
						timingResMinimum=0.58;                // TTS: 3.0ns;
				}
				else if(digitspst=="lappd_v1"){               // LAPPD
						timingConstant = 0.04;                // not yet implemented
						timingResConstant=0.001;
						timingResMinimum=0.05;                // TTS: 50ps;
						timingResPower=-0.7;
				}
				double digitqforsmearing = (digitsq > 0.5) ? digitsq : 0.5;
				double filedigittsmeared;
				if(digitspst.find("lappd")== std::string::npos){
					filedigittsmeared = timingResConstant + sqrt(timingConstant/digitqforsmearing);
					if (filedigittsmeared < timingResMinimum) filedigittsmeared=timingResMinimum;
					// XXX saturates at 0.5 and 2.27 << due to digitqs minimum
				} else {
					digitqforsmearing = (digitsq > 0.5) ? digitsq : 0.5;
					// based on a roughly similar form that gives ~60ps for Q~30 and ~5ps for Q~1
					filedigittsmeared = // XXX if changing ensure it is mirrored below in lappd digit sec
						timingResConstant + timingConstant*pow(digitqforsmearing,timingResPower);
					if (filedigittsmeared < timingResMinimum) filedigittsmeared=timingResMinimum;
				}
				filedigittsmears.push_back(filedigittsmeared);
				//float Smearing_factor = G4RandGauss::shoot(0.0,filedigittsmeared);
				digitsqpmthist->Fill(digitsq);
				digitstpmthist->Fill(absolutedigitst);
				digittsmearpmthist->Fill(filedigittsmeared);
				pmttimesmearvsqhist->Fill(filedigittsmeared,digitsq);

#ifdef LAPPD_DEBUG
				intileposx.push_back(0);
				intileposy.push_back(0);
#ifdef FILE_VERSION>3
				poserrx.push_back(0);
				poserry.push_back(0);
				poserrz.push_back(0);
#endif
				tileorient.push_back(0);
				octagonside.push_back(0);
#endif
				
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
			muonenergydepositions.push_back(tankchargefrommuon);
			totaltankcharge=atrigt->GetSumQ();
			
			if((chargeinsidecherenkovcone+chargeoutsidecherenkovcone)!=0){
				fractionalchargeincone=
				chargeinsidecherenkovcone/(chargeinsidecherenkovcone+chargeoutsidecherenkovcone);
			} else {
				fractionalchargeincone=0;
			}
			
			// end of digit analysis
			////////////////////////////////////////////////
			
			// ----------------------------------------------------------------------------------------------
			// lappd digit analysis
			// ----------------------------------------------------------------------------------------------
			// FIXME: the vectors store all digits in the event regardless of which trigger they are in. 
			// we should scan through and only retrieve / fill jingbo file vectors with events within
			// the trigger window.
			
			int runningcount=0;
			for(int lappdi=0; lappdi<lappd_numhitsthisevt; lappdi++){
				// loop over LAPPDs that had at least one hit
				int LAPPDID = lappd_hittile[lappdi];
				double tileposx = lappd_hittilesposx[lappdi];  // position of LAPPD in global coords
				double tileposy = lappd_hittilesposy[lappdi];
				double tileposz = lappd_hittilesposz[lappdi];
				WCSimRootPMT pmt = geo->GetLAPPD(LAPPDID-1);
				double pmtx = pmt.GetPosition(0);              // verified these are equivalent ^
				double pmty = pmt.GetPosition(1);
				double pmtz = pmt.GetPosition(2);
				int thepmtsloc = pmt.GetCylLoc();
				double Rinnerstruct=270.9/2.; // cm octagonal inner structure radius = 106.64"
				double Rthresh=Rinnerstruct*pow(2.,-0.5);
				bool printthis=false;
				double tileangle;
				int theoctagonside;
				switch (thepmtsloc){
					case 0: break; // top cap
					case 2: break; // bottom cap
					case 1: // wall
						// we need to account for the angle of the LAPPD within the tank
						// determine the angle based on it's position
						double octangle1=TMath::Pi()*(3./8.);
						double octangle2=TMath::Pi()*(1./8.);
						pmtz+=-tank_start-tank_radius;
							 if(pmtx<-Rthresh&&pmtz<0)         {tileangle=-octangle1; theoctagonside=0;}
						else if(-Rthresh<pmtx&&pmtx<0&&pmtz<0) {tileangle=-octangle2; theoctagonside=1;}
						else if(0<pmtx&&pmtx<Rthresh&&pmtz<0)  {tileangle=octangle2;  theoctagonside=2;}
						else if(Rthresh<pmtx&&pmtz<0)          {tileangle=octangle1;  theoctagonside=3;}
						else if(pmtx<-Rthresh&&pmtz>0)         {tileangle=octangle1;  theoctagonside=4;}
						else if(-Rthresh<pmtx&&pmtx<0&&pmtz>0) {tileangle=octangle2;  theoctagonside=5;}
						else if(0<pmtx&&pmtx<Rthresh&&pmtz>0)  {tileangle=-octangle2; theoctagonside=6;}
						else if(Rthresh<pmtx&&pmtz>0)          {tileangle=-octangle1; theoctagonside=7;}
						break;
				}
				int numhitsthislappd=lappd_numphots[lappdi];
				int lastrunningcount=runningcount;
				// loop over all the hits on this lappd
				for(;runningcount<(lastrunningcount+numhitsthislappd); runningcount++){
					double peposx = lappd_hitpeposx.at(runningcount);      // position of hit on tile in LAPPD
					double peposy = lappd_hitpeposy.at(runningcount);
					//double peposz   = lappd_hitpeposz.at(runningcount);  // no apparent meaning
					double digitsx, digitsy, digitsz;
					switch (thepmtsloc){
						case 0: // top cap
							digitsx  = tileposx - peposx;  // perfect
							digitsy  = tileposy;           // mostly error +9.2? but some down to -9.2?
							digitsz  = tileposz - peposy;  // perfect
							break;
						case 2: // bottom cap
							digitsx  = tileposx + peposx;  // perfect
							digitsy  = tileposy;           // mostly error -9.2? but some up to +9.2?
							digitsz  = tileposz - peposy;  // perfect
							break;
						case 1: // wall
							digitsy  = tileposy + peposx; // perfect
							if(theoctagonside<4){
									digitsx  = tileposx - peposy*TMath::Cos(tileangle);
									digitsz  = tileposz - peposy*TMath::Sin(tileangle);
							} else {
									digitsx  = tileposx + peposy*TMath::Cos(tileangle);
									digitsz  = tileposz + peposy*TMath::Sin(tileangle);
							}
					}
#ifdef LAPPD_DEBUG
					static int printcount1=0;
					static int printcount2=0;
					static int printcount3=0;
					static std::vector<int>printedtiles(30,-1);
					if( printthis &&
					   ((printcount1<10&&thepmtsloc==1) ||
						(printcount2<10&&thepmtsloc==0) ||
						(printcount3<10&&thepmtsloc==2)) &&
					   (find(printedtiles.begin(),printedtiles.end(),LAPPDID)==printedtiles.end())){
						switch (thepmtsloc){
							case 1: printcount1++; cout<<"wall "; break;
							case 0: printcount2++; cout<<"top cap "; break;
							case 2: printcount3++; cout<<"bottom cap "; break;
						}
						cout<<"lappd position is: ("
							<<tileposx<<", "<<tileposy<<", "<<tileposz<<")"<<endl;
						cout<<"pepos is: ("<<peposx<<", "<<peposy<<")"<<endl;
						cout<<"derived position is: ("
							<<digitsx<<", "<<digitsy<<", "<<digitsz<<")"<<endl
#ifdef FILE_VERSION>3
							<<"true global position is: ("<<lappd_hitglobalposx.at(runningcount)
							<<", "<<lappd_hitglobalposy.at(runningcount)
							<<", "<<lappd_hitglobalposz.at(runningcount)<<")"<<endl;
#endif
					}
#endif
					double digitst  = lappd_hittruetime.at(runningcount);
					WCSimRootEventHeader* trigheader=atrigt->GetHeader();
					double triggertime=trigheader->GetDate();
#if FILE_VERSION<10 // dont know what file version this will be fixed in 
					// n.b. all lappd digits are stored regardless of being in or outside trigger window.
					// so this time may be way out, this digit may even be part of a different trigger.
					// some sort of loop over trigger times to sort it here....
					// (we could check if abs time > trigger time for next trig...)
#endif
					double absolutedigitst=digitst+950.-triggertime;
					ROOT::Math::XYZTVector adigitvector =  // convert mm to cm 
						ROOT::Math::XYZTVector(digitsx/10.,digitsy/10.,digitsz/10.,absolutedigitst);
					
#ifdef LAPPD_DEBUG
					intileposx.push_back(peposx);
					intileposy.push_back(peposy);
#ifdef FILE_VERSION>3
					poserrx.push_back(digitsx-lappd_hitglobalposx.at(runningcount));
					poserry.push_back(digitsy-lappd_hitglobalposy.at(runningcount));
					poserrz.push_back(digitsz-lappd_hitglobalposz.at(runningcount));
#endif
					tileorient.push_back(thepmtsloc);
					octagonside.push_back(theoctagonside);
#endif
#if FILE_VERSION>10 // don't know when this will be fixed
					double digitsq = lappd_hitcharge.at(runningcount);
#else
					// lappd hit charge is not stored - only a count of pes within the whole event is stored.
					// for now we'll calculate the charge here using the formula WCSim would have used:
					double random = mrand->Rndm();
					double random2 = mrand->Rndm(); 
					float *qpe0 = Getlappdqpe();
					int irand;
					for(irand = 0; irand < 501; irand++){ if (random <= *(qpe0+irand)) break; }
					if(irand==500) random = mrand->Rndm();
					double digitsq = (double(irand-50) + random2)/22.83;
					// FIXME: digit charges for LAPPDs have a big peak at <1 charge, because there's
					// no digitizer integration and rejection step. Don't know if we can transfer this... 
//					int iflag;
//					SKIDigitizerThreshold(digitsq,iflag);
//					if(iflag!=0) continue; // this digit gets rejected by digitizer ... if it were a PMT
//					// TODO need to move all vector push_back's after this if we want to 'continue' 
//					// in order to keep all vectors the same size
//					digitsq*0.985;         // efficiency in WCSimWCDigitizer
#endif
					filedigitvertices.push_back(adigitvector);
					filedigitQs.push_back(digitsq);
					filedigitsensortypes.push_back("lappd_v0"); // bad timing resolution used for pulses
					double timingConstant = 0.04;
					double timingResConstant=0.001;
					double timingResMinimum=0.005;              // TTS: 50ps;
					double timingResPower=-0.7;
					double digitqforsmearing = (digitsq > 0.5) ? digitsq : 0.5;
					// based on a roughly similar form that gives ~60ps for Q~1 and ~5ps for Q~30
					double filedigittsmeared = // XXX if changing this ensure it is mirrored above.
						timingResConstant + timingConstant*pow(digitqforsmearing,timingResPower);
					// saturates @ 0.065 with majority of entries.
					if (filedigittsmeared < timingResMinimum) filedigittsmeared=timingResMinimum;
					filedigittsmears.push_back(filedigittsmeared);
					digitsqlappdhist->Fill(digitsq);
					digitstlappdhist->Fill(absolutedigitst);
					digittsmearlappdhist->Fill(filedigittsmeared);
					lappdtimesmearvsqhist->Fill(filedigittsmeared,digitsq);
				}
			}
			
			// end of lappd digit analysis
			////////////////////////////////////////////////
			
			// Fill simplified file for reconstruction dev
			///////////////////////////////////////////////
			filemuonstartvertex = primarystartvertex;
			filemuonstopvertex = primarystopvertex;
			filemuondirectionvector = differencevector.Unit();
			vertextreenocuts->Fill();
			if(isinfiducialvol){
					vertextreefiducialcut->Fill();
					nummuontracksinfidvol++;
			}
			
			// ----------------------------------------------------------------------------------------------
			// mrd digit analysis
			// ----------------------------------------------------------------------------------------------
			numMRDdigits = atrigm->GetCherenkovDigiHits()->GetEntries();
			totMRDcharge = atrigm->GetSumQ();
			numMRDdigitsfrommu=0;
			MRDchargefrommu=0.;
			mrdtubeshitbymu.clear();
			for(Int_t i=0; i<numMRDdigits; i++){
				// retrieve the digit information
				/////////////////////////////////
				WCSimRootCherenkovDigiHit* thedigihit = 
					(WCSimRootCherenkovDigiHit*)atrigm->GetCherenkovDigiHits()->At(i);
				int digitstubeid = thedigihit->GetTubeId()-1;
				double digitsq = thedigihit->GetQ();
#if FILE_VERSION<2
				// we need to correct the indices returned by GetPhotonIds with the PMT offset
				// to be able to retrieve digit parents and extract digits from primary mu
				int timeArrayOffset=-1;
				int ncherenkovhits=atrigm->GetNcherenkovhits();
				for(int ipmt = 0; ipmt < ncherenkovhits; ipmt++) {
					WCSimRootCherenkovHit* hitobject = 
						(WCSimRootCherenkovHit*)atrigm->GetCherenkovHits()->At(ipmt);
					int tubeNumber = hitobject->GetTubeID()-1;
					if(tubeNumber==digitstubeid) {
						timeArrayOffset = hitobject->GetTotalPe(0);
						break;
					}
				}
#endif
				// Calculate total charge from muon
				////////////////////////////////////
				std::vector<int> truephotonindices = thedigihit->GetPhotonIds();
				int numphotonsfromthismuon=0;
				for(int truephoton=0; truephoton<truephotonindices.size(); truephoton++){
					int thephotonsid = truephotonindices.at(truephoton);
#if FILE_VERSION<2
					thephotonsid+=timeArrayOffset;
#endif
					WCSimRootCherenkovHitTime *thehittimeobject = 
						(WCSimRootCherenkovHitTime*)atrigm->GetCherenkovHitTimes()->At(thephotonsid);
					Int_t thephotonsparenttrackid = thehittimeobject->GetParentID();
					if(thephotonsparenttrackid==nextrack->GetId()) {
						numphotonsfromthismuon++;
					}
				}
				// scale each digit charge by the fraction of hits in the digit which come from the muon
				if(numphotonsfromthismuon!=0){
					MRDchargefrommu+= 
						digitsq * ((double)numphotonsfromthismuon/(double)truephotonindices.size());
					numMRDdigitsfrommu++;
					if(std::find(mrdtubeshitbymu.begin(), mrdtubeshitbymu.end(), digitstubeid)==mrdtubeshitbymu.end())
						mrdtubeshitbymu.push_back(digitstubeid);
				}
				
			}  // end loop over mrd digits
			
			// end of mrd digit analysis
			////////////////////////////////////////////////
			
//			for(std::map<std::string,bool>::iterator mapit=eventtypes.begin(); mapit!=eventtypes.end(); mapit++)
//				cout<<"eventtypes."<<mapit->first<<" = "<<mapit->second<<endl;
//			cout<<"isintank="<<isintank<<endl;
//			cout<<"isinfiducialvol="<<isinfiducialvol<<endl;
//			cout<<"eventq2="<<eventq2<<endl;
//			cout<<"eventEnu="<<eventEnu<<endl;
//			cout<<"neutrinopdg="<<neutrinopdg<<endl;
//			cout<<"hasmuon="<<hasmuon<<endl;
//			cout<<"mustartvtx=("<<mustartvtx.X()<<","<<mustartvtx.Y()<<","<<mustartvtx.Z()<<")"<<endl;
//			cout<<"mustopvtx=("<<mustopvtx.X()<<","<<mustopvtx.Y()<<","<<mustopvtx.Z()<<")"<<endl;
//			cout<<"mustartE="<<mustartE<<endl;
//			cout<<"muonentersMRD="<<muonentersMRD<<endl;
//			cout<<"muonstopsinMRD="<<muonstopsinMRD<<endl;
//			cout<<"muonrangesoutMRD="<<muonrangesoutMRD<<endl;
//			cout<<"mrdpenetrationcm="<<mrdpenetrationcm<<endl;
//			cout<<"mrdpenetrationlayers="<<mrdpenetrationlayers<<endl;
//			cout<<"mutracklengthinMRD="<<mutracklengthinMRD<<endl;
//			cout<<"mutracklengthintank="<<mutracklengthintank<<endl;
//			cout<<"tankchargefrommuon="<<tankchargefrommuon<<endl;
//			cout<<"fractionalchargeincone="<<fractionalchargeincone<<endl;
			
			bGenieFileString->Fill();
			bGenieEventNum->Fill();
			bDirtFileString->Fill();
			bDirtEventNum->Fill();
			bWCSimFileString->Fill();
			bLAPPDFileString->Fill();
			bWCSimEventNum->Fill();
			bInTank->Fill();
			bEventType->Fill();
			bIsQuasiElastic->Fill();
			bIsResonant->Fill();
			bIsDeepInelastic->Fill();
			bIsCoherent->Fill();
			bIsDiffractive->Fill();
			bIsInverseMuDecay->Fill();
			bIsIMDAnnihilation->Fill();
			bIsSingleKaon->Fill();
			bIsNuElectronElastic->Fill();
			bIsEM->Fill();
			bIsWeakCC->Fill();
			bIsWeakNC->Fill();
			bIsMEC->Fill();
			bEventQ2->Fill();
			bEventEnu->Fill();
			bNeutrinoPdg->Fill();
			bInFidVol->Fill();
			bEventHasMuon->Fill();
			bMuonStartVtx->Fill();
			bMuonStopVtx->Fill();
			bMuonAngle->Fill();
			bMuonStartE->Fill();
			bMuonTrackLengthInTank->Fill();
			bMuonMrdPenetrationInCm->Fill();
			bMuonMrdPenetrationLayers->Fill();
			bMuonEntersMRD->Fill();
			bMuonStopsInMRD->Fill();
			bMuonRangesOutMRD->Fill();
			bMuonTrackLengthInMRD->Fill();
			bNumTankDigits->Fill();
			bTankChargeFromMuon->Fill();
			bNumTankDigitsFromMu->Fill();
			bTankTubesHitByMuon->Fill();
			bFractionOfMuonChargeInCone->Fill();
			bTotalTankCharge->Fill();
			bTotalUpstreamCharge->Fill();
			bTotalDownstreamCharge->Fill();
			bTopCapCharge->Fill();
			bBottomCapCharge->Fill();
			bNumPiZeroTracks->Fill();
			bNumPiPlusTracks->Fill();
			bNumPiMinusTracks->Fill();
			bNumGammaTracks->Fill();
			bNumNeutronTracks->Fill();
			bNumProtonTracks->Fill();
			bNumMuTracks->Fill();
			bTotalMrdDigits->Fill();
			bTotalMrdCharge->Fill();
			bMrdDigitsFromMuon->Fill();
			bMrdChargeFromMuon->Fill();
			bMrdTubesHitByMuon->Fill();
			bTrackPDG->Fill();
			bTrackStartPos->Fill();
			bTrackStopPos->Fill();
			bTrackStartMom->Fill();
			bTrackStartE->Fill();
			bTrackStartVol->Fill();
			bTrackStopVol->Fill();
			bTrackParentType->Fill();
			bTrackStopTime->Fill();
			
			// if we've got as far as this (i.e. haven't triggered a 'continue' above)
			// then we've had a muon track. Let's just consider one muon for now
			// otherwise we need to add _ALL_ this stuff to vectors and scan through that. 
			// break; // or maybe not... mayb it'll work just as is. 
			
		}  // end of loop over tracks
		
		// ==============================================================================================
		// check if we found any suitable muon tracks
		// ==============================================================================================
		
		if(scatteringanglesvector.size()>1) {
			cerr<<"Mutliple accepted final state leptons!"<<endl;
			for(int i=0; i<scatteringanglesvector.size(); i++){
				cout<<"wcsim lepton angle "<<i<<" = "<<scatteringanglesvector.at(i)<<endl;
				cout<<"wcsim lepton energy "<<i<<" = "<<primaryenergiesvector.at(i)<<endl;
				cout<<"wcsim lepton trackid "<<i<<" = "<<acceptedtrackids.at(i)<<endl;
			}
			cout<<"compare to:"<<endl;
			cout<<"genie lepton angle = "<<thegenieinfo.fslanglegenie<<endl;
			cout<<"genie lepton energy = "<<thegenieinfo.fsleptonenergy<<endl;
			
			std::vector<Int_t>::iterator minit = std::min_element(acceptedtrackids.begin(), acceptedtrackids.end());		
			Int_t indextouse = std::distance(acceptedtrackids.begin(),minit);
			incidentneutrinoenergiesacceptedwcsim->Fill(neutrinoenergiesvector.at(indextouse));
			fslanglesacceptedwcsim->Fill(scatteringanglesvector.at(indextouse));
			fslenergiesacceptedwcsim->Fill(primaryenergiesvector.at(indextouse));
			eventq2acceptedwcsim->Fill(q2vector.at(indextouse));
			muedepositionsacceptedwcsim->Fill(muonenergydepositions.at(indextouse));
			
			// fiducial cut versions
			if(isinfiducialvol){
				incidentneutrinoenergiesacceptedwcsimfidcut->Fill(neutrinoenergiesvector.at(indextouse));
				fslanglesacceptedwcsimfidcut->Fill(scatteringanglesvector.at(indextouse));
				fslenergiesacceptedwcsimfidcut->Fill(primaryenergiesvector.at(indextouse));
				eventq2acceptedwcsimfidcut->Fill(q2vector.at(indextouse));
			}
		} else if (scatteringanglesvector.size()==1) { // just one matched wcsim track. 
			incidentneutrinoenergiesacceptedwcsim->Fill(neutrinoenergiesvector.at(0));
			fslanglesacceptedwcsim->Fill(scatteringanglesvector.at(0));
			fslenergiesacceptedwcsim->Fill(primaryenergiesvector.at(0));
			eventq2acceptedwcsim->Fill(q2vector.at(0));
			muedepositionsacceptedwcsim->Fill(muonenergydepositions.at(0));
			// fiducial cut versions
			if(isinfiducialvol){
				// fiducial cut versions
				incidentneutrinoenergiesacceptedwcsimfidcut->Fill(neutrinoenergiesvector.at(0));
				fslanglesacceptedwcsimfidcut->Fill(scatteringanglesvector.at(0));
				fslenergiesacceptedwcsimfidcut->Fill(primaryenergiesvector.at(0));
				eventq2acceptedwcsimfidcut->Fill(q2vector.at(0));
			}
#ifdef VERBOSE
			cout<<"wcsim lepton angle = "<<scatteringanglesvector.at(0)<<endl;
			cout<<"wcsim lepton energy = "<<primaryenergiesvector.at(0)<<endl;
			cout<<"wcsim lepton trackid = "<<acceptedtrackids.at(0)<<endl;
			cout<<"compare to:"<<endl;
			cout<<"genie lepton angle = "<<thegenieinfo.fslanglegenie<<endl;
			cout<<"genie lepton energy = "<<thegenieinfo.fsleptonenergy<<endl;
#endif
		} else {
#ifdef VERBOSE
			cout<<"no accepted wcsim track"<<endl;
#endif
		}
		if(scatteringanglesvector.size()!=0){
#ifdef VERBOSE
			cout<<"MATCHED A TRACK YEEEEAH"<<endl; /*return;*/
#endif
			// by getting this far we've found a primary muon that penetrated the MRD.
			// fill genie 'accepted' histograms.
			incidentneutrinoenergiesaccepted->Fill(thegenieinfo.probeenergy);
			//incidentneutrinoanglesaccepted->Fill(thegenieinfo.probeangle);
			fslanglesaccepted->Fill(thegenieinfo.fslanglegenie);
			fslenergiesaccepted->Fill(thegenieinfo.fsleptonenergy);
			eventq2accepted->Fill(thegenieinfo.Q2);
			neutrinovertexQEaccepted->Fill(thegenieinfo.genie_x, thegenieinfo.genie_y, thegenieinfo.genie_z);
			// if it passes the fiducial cut, we've also accepted it, fill.
			if(isinfiducialvol){
				// genie values
				incidentneutrinoenergiesacceptedfidcut->Fill(thegenieinfo.probeenergy);
				fslanglesacceptedfidcut->Fill(thegenieinfo.fslanglegenie);
				fslenergiesacceptedfidcut->Fill(thegenieinfo.fsleptonenergy);
				eventq2acceptedfidcut->Fill(thegenieinfo.Q2);
				// fill flat tree for reconstruction dev
				vertextreefiducialmrd->Fill();
				numCCQEneutrinoeventsinfidvolmrd++;
				nummuontracksinfidvolmrd+=scatteringanglesvector.size();
			}
		}
		
		flateventfileout->Write("",TObject::kOverwrite);
	}  // end of loop over events
	
	//======================================================================================================
	//======================================================================================================
	
	cout<<"generating scaled histograms"<<endl;
	Double_t numbeamspills = totalpots/(4.0 * TMath::Power(10.,12.));
	Double_t numbeamspillsperday = (24.*60.*60.*1000.)/133.3333;	// 24 hours in ms / 133.33 ms between spills
	Double_t numdays = numbeamspills/numbeamspillsperday;
	cout<<"Results based on "<<totalpots<<" POTs, or "<<numbeamspills<<" beam spills, or "<<numdays<<" days of data"<<endl;
	cout<<"There were "<<numneutrinoeventsintank<<" neutrino interactions in the tank, of which "<<numCCQEneutrinoeventsintank<<" were true CCQE events."<<endl;
	cout<<"Of those, "<<numCCQEneutrinoeventsinfidvol<<" were within the fiducial volume."<<endl;
	cout<<"Of those in turn, "<<numCCQEneutrinoeventsinfidvolmrd<<" produced an accepted MRD muon"<<endl;
	cout<<"There were "<<nummuontracksintank<<" muons in the tank, of which "
		<<nummuontracksinfidvol<<" were from (CCQE?) events in the fiducial volume."<<endl;
	/* cout<<"There were "<<nummuontracksinmrd<<" muons from tank events which passed through 3 MRD layers"
		<<" of which "<<nummuontracksinfidvolmrd<<" originated from vertices in the fiducial volume"<<endl;*/
	// ^^^ this is currently disabled.
	
	
//	numneutrinoeventsintank
//	numCCQEneutrinoeventsintank
//	numCCQEneutrinoeventsinfidvol
//	numCCQEneutrinoeventsinfidvolmrd
//	nummuontracksintank
//	nummuontracksinfidvol
	
	// TODO: neutrinovertex, neutrinovertexQE, neutrinovertexQEaccepted TH3D plots: save projections
	
	gROOT->cd();
	// create scaled histograms with bin contents scaled to the number of POTs in input files:
	TH1D* incidentneutrinoenergiesallscaled = new TH1D(*incidentneutrinoenergiesall);
	TH1D* incidentneutrinoenergiesacceptedscaled = new TH1D(*incidentneutrinoenergiesaccepted);
	TH1D* incidentneutrinoenergiesacceptedwcsimscaled = new TH1D(*incidentneutrinoenergiesacceptedwcsim);
	TH1D* fslanglesallscaled = new TH1D(*fslanglesall);
	TH1D* fslanglesacceptedscaled = new TH1D(*fslanglesaccepted);
	TH1D* fslanglesacceptedwcsimscaled = new TH1D(*fslanglesacceptedwcsim);
	TH1D* fslenergiesallscaled = new TH1D(*fslenergiesall);
	TH1D* fslenergiesacceptedscaled = new TH1D(*fslenergiesaccepted);
	TH1D* fslenergiesacceptedwcsimscaled = new TH1D(*fslenergiesacceptedwcsim);
	TH1D* eventq2allscaled = new TH1D(*eventq2all);
	TH1D* eventq2acceptedscaled = new TH1D(*eventq2accepted);
	TH1D* eventq2acceptedwcsimscaled = new TH1D(*eventq2acceptedwcsim);
	
	TH1D* incidentneutrinoenergiesallfidcutscaled = new TH1D(*incidentneutrinoenergiesallfidcut);
	TH1D* incidentneutrinoenergiesacceptedfidcutscaled = new TH1D(*incidentneutrinoenergiesacceptedfidcut);
	TH1D* incidentneutrinoenergiesacceptedwcsimfidcutscaled = new TH1D(*incidentneutrinoenergiesacceptedwcsimfidcut);
	TH1D* fslanglesallfidcutscaled = new TH1D(*fslanglesallfidcut);
	TH1D* fslanglesacceptedfidcutscaled = new TH1D(*fslanglesacceptedfidcut);
	TH1D* fslanglesacceptedwcsimfidcutscaled = new TH1D(*fslanglesacceptedwcsimfidcut);
	TH1D* fslenergiesallfidcutscaled = new TH1D(*fslenergiesallfidcut);
	TH1D* fslenergiesacceptedfidcutscaled = new TH1D(*fslenergiesacceptedfidcut);
	TH1D* fslenergiesacceptedwcsimfidcutscaled = new TH1D(*fslenergiesacceptedwcsimfidcut);
	TH1D* eventq2allfidcutscaled = new TH1D(*eventq2allfidcut);
	TH1D* eventq2acceptedfidcutscaled = new TH1D(*eventq2acceptedfidcut);
	TH1D* eventq2acceptedwcsimfidcutscaled = new TH1D(*eventq2acceptedwcsimfidcut);
	
	//TODO: save histos to root file so they can be scaled arbitrarily without regenerating!
	
	TH1D* placeholder=0;
	
	std::vector<TH1D*> scaledhistopointers {
		incidentneutrinoenergiesallscaled,
		incidentneutrinoenergiesacceptedscaled,
		incidentneutrinoenergiesacceptedwcsimscaled,
		fslanglesallscaled,
		fslanglesacceptedscaled,
		fslanglesacceptedwcsimscaled,
		fslenergiesallscaled,
		fslenergiesacceptedscaled,
		fslenergiesacceptedwcsimscaled,
		eventq2allscaled,
		eventq2acceptedscaled,
		eventq2acceptedwcsimscaled,
		
		incidentneutrinoenergiesallfidcutscaled,
		incidentneutrinoenergiesacceptedfidcutscaled,
		incidentneutrinoenergiesacceptedwcsimfidcutscaled,
		fslanglesallfidcutscaled,
		fslanglesacceptedfidcutscaled,
		fslanglesacceptedwcsimfidcutscaled,
		fslenergiesallfidcutscaled,
		fslenergiesacceptedfidcutscaled,
		fslenergiesacceptedwcsimfidcutscaled,
		eventq2allfidcutscaled,
		eventq2acceptedfidcutscaled,
		eventq2acceptedwcsimfidcutscaled,
		
		incidentneutrinoenergiesallscaled,
		incidentneutrinoenergiesallfidcutscaled,
		placeholder,
		fslanglesallscaled,
		fslanglesallfidcutscaled,
		placeholder,
		fslenergiesallscaled,
		fslenergiesallfidcutscaled,
		placeholder,
		eventq2allscaled,
		eventq2allfidcutscaled,
		placeholder,
		
		incidentneutrinoenergiesallscaled,
		incidentneutrinoenergiesacceptedfidcutscaled,
		placeholder,
		fslanglesallscaled,
		fslanglesacceptedfidcutscaled,
		placeholder,
		fslenergiesallscaled,
		fslenergiesacceptedfidcutscaled,
		placeholder,
		eventq2allscaled, 
		eventq2acceptedfidcutscaled,
		placeholder
	};
	
	std::vector<TH1D*> histopointers {
		incidentneutrinoenergiesall,
		incidentneutrinoenergiesaccepted,
		incidentneutrinoenergiesacceptedwcsim,
		fslanglesall,
		fslanglesaccepted,
		fslanglesacceptedwcsim,
		fslenergiesall,
		fslenergiesaccepted,
		fslenergiesacceptedwcsim,
		eventq2all,
		eventq2accepted,
		eventq2acceptedwcsim,
		
		incidentneutrinoenergiesallfidcut,
		incidentneutrinoenergiesacceptedfidcut,
		incidentneutrinoenergiesacceptedwcsimfidcut,
		fslanglesallfidcut,
		fslanglesacceptedfidcut,
		fslanglesacceptedwcsimfidcut,
		fslenergiesallfidcut,
		fslenergiesacceptedfidcut,
		fslenergiesacceptedwcsimfidcut,
		eventq2allfidcut,
		eventq2acceptedfidcut,
		eventq2acceptedwcsimfidcut,
		
		incidentneutrinoenergiesall,
		incidentneutrinoenergiesallfidcut,
		placeholder,
		fslanglesall,
		fslanglesallfidcut,
		placeholder,
		fslenergiesall,
		fslenergiesallfidcut,
		placeholder,
		eventq2all,
		eventq2allfidcut,
		placeholder,
		
		incidentneutrinoenergiesall,
		incidentneutrinoenergiesacceptedfidcut,
		placeholder,
		fslanglesall,
		fslanglesacceptedfidcut,
		placeholder,
		fslenergiesall,
		fslenergiesacceptedfidcut,
		placeholder,
		eventq2all,
		eventq2acceptedfidcut,
		placeholder
	};
	
	std::vector<std::string> legendstrings {
		// first 4 graphs, 3 plots each, compare effect of MRD cut on all true tank QE events
		"all incident",
		"MRD cut (truth)",
		"MRD cut (reco)",
		"all incident",
		"MRD cut (truth)",
		"MRD cut (reco)",
		"all incident",
		"MRD cut (truth)",
		"MRD cut (reco)",
		"all incident",
		"MRD cut (truth)",
		"MRD cut (reco)",
		// next 4 graphs, 3 plots each, compare effect of MRD cut on all true QE events within fiducial volume
		"fiducial incident",
		"fiducial w/ MRD cut (truth)",
		"fiducial w/ MRD cut (reco)",
		"fiducial incident",
		"fiducial w/ MRD cut (truth)",
		"fiducial w/ MRD cut (reco)",
		"fiducial incident",
		"fiducial w/ MRD cut (truth)",
		"fiducial w/ MRD cut (reco)",
		"fiducial incident",
		"fiducial w/ MRD cut (truth)",
		"fiducial w/ MRD cut (reco)",
		// next 4 graphs, 2 plots each, compare effect of fiducial cut on true tank QE all events
		"all incident", 
		"fiducial incident",
		"all incident", 
		"fiducial incident",
		"all incident", 
		"fiducial incident",
		"all incident", 
		"fiducial incident",
		// final 4 graphs, 2 plots each, compare effect of combined fiducial+MRD cuts on all true tank QE events
		"all incident", 
		"fiducial w/ MRD cut (truth)",
		"all incident", 
		"fiducial w/ MRD cut (truth)",
		"all incident", 
		"fiducial w/ MRD cut (truth)",
		"all incident", 
		"fiducial w/ MRD cut (truth)"
	};
	
	Double_t norm = 1.;
	for(int i=0; i<scaledhistopointers.size(); i++){
		TH1D* temp = scaledhistopointers.at(i);
		TH1D* temp2 = histopointers.at(i);
		if(temp==0||temp2==0) continue;
		TString thetitle = temp2->GetTitle();
		cout<<thetitle<<" has "<<temp2->GetEntries()<<" entries"<<endl;
//		for(int j=1; j<(temp->GetNbinsX()-2); j++){
//			Double_t thecontent = temp2->GetBinContent(j);
//			if(thecontent!=0){ thecontent /= temp->GetEntries(); /*numdays*/}
//			temp->SetBinContent(j,thecontent);
//		}
		temp->Scale(norm/temp->Integral());
		temp->GetYaxis()->SetRangeUser(0.,0.1);
		TString tempname = temp2->GetName();
		temp->SetName(TString::Format("%s_scaled",tempname.Data()));
	}
	
	// draw the original histograms
	// ============================
//	Double_t win_width=500;
//	Double_t win_height=400;
//	Double_t win_scale=0.7;
//	Int_t n_wide=2;
//	Int_t n_high=3;
//	TCanvas* c1 = new TCanvas("c1","c1",win_width*n_wide*win_scale,win_height*n_high*win_scale);
	cout<<"drawing histograms"<<endl;
	std::vector<TCanvas*> canvaspointers;
	std::vector<TLegend*> legendpointers;
	Int_t legendindex=0;
#ifdef DRAW_HISTOS
	for(int i=0; i<histopointers.size(); i++){
		TLegend* leg;
		TCanvas* canv;
		TH1D* hist = histopointers.at(i);
		if(hist==0) continue;
		TH1D* upnext=0;
		if((i+1)<histopointers.size()){
			upnext = histopointers.at(i+1);
		}
		if(i%3==0){
			canv = new TCanvas(TString::Format("c%d",i),TString::Format("c%d",i));
			canvaspointers.push_back(canv);
			canv->cd();
			leg = new TLegend(0.5,0.78,0.7,0.88);
			legendpointers.push_back(leg);
			leg->AddEntry(hist,legendstrings.at(legendindex).c_str(),"l");
			legendindex++;
			hist->SetLineColor(kRed);
			hist->Draw();
		} else if(i%3==1){
			hist->SetLineColor(kBlue);
			leg->AddEntry(hist,legendstrings.at(legendindex).c_str(),"l");
			legendindex++;
			hist->Draw("same");
			if(upnext==0){
				leg->Draw();
			}
		} else {
			hist->SetLineColor(kViolet);
			hist->SetLineStyle(2);
			leg->AddEntry(hist,legendstrings.at(legendindex).c_str(),"l");
			legendindex++;
			hist->Draw("same");
			leg->Draw();
		}
	}
#endif
	
	// draw the scaled histograms
	// ==========================
	std::vector<TCanvas*> scaledcanvaspointers;
	std::vector<TLegend*> scaledlegendpointers;
	legendindex=0;
#ifdef DRAW_HISTOS
	for(int i=0; i<scaledhistopointers.size(); i++){
		TLegend* leg;
		TCanvas* canv;
		TText* numentstitle;
		TH1D* hist = scaledhistopointers.at(i);
		if(hist==0) continue;
		TH1D* upnext=0;
		if((i+1)<scaledhistopointers.size()){
			upnext = scaledhistopointers.at(i+1);
		}
		if(i%3==0){
			canv = new TCanvas(TString::Format("cc%d",i),TString::Format("cc%d",i));
			scaledcanvaspointers.push_back(canv);
			canv->cd();
			leg = new TLegend(0.5,0.78,0.7,0.88);
			scaledlegendpointers.push_back(leg);
			leg->AddEntry(hist,legendstrings.at(legendindex).c_str(),"l");
			legendindex++;
			hist->SetLineColor(kRed);
			hist->Draw();
		} else if(i%3==1){
			hist->SetLineColor(kBlue);
			leg->AddEntry(hist,legendstrings.at(legendindex).c_str(),"l");
			legendindex++;
			hist->Draw("same");
			if(upnext==0){
				leg->Draw();
			}
		} else {
			hist->SetLineColor(kViolet);
			hist->SetLineStyle(2);
			leg->AddEntry(hist,legendstrings.at(legendindex).c_str(),"l");
			legendindex++;
			hist->Draw("same");
			leg->Draw();
		}
	}
	
//	 draw other canvases
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
	aleg->AddEntry(muedepositionsacceptedwcsim,"Mu Tank E Dep (accepted)","l");
	muedepositionswcsim->SetLineColor(kBlue);
	muedepositionswcsim->Draw();
	muedepositionsacceptedwcsim->SetLineColor(kRed);
	muedepositionsacceptedwcsim->Draw("same");
	aleg->Draw();
	debug1->SaveAs("muon_tank_energy_depositions.png");
	if(muedepositionswcsim) delete muedepositionswcsim; muedepositionswcsim=0;
	if(muedepositionsacceptedwcsim) delete muedepositionsacceptedwcsim; muedepositionsacceptedwcsim=0;
	
	if(debug1) delete debug1; debug1=0;
	if(aleg) delete aleg; aleg=0;
	cout<<"total of "<<nummuontracksintank<<" muon tracks in the tank, of which "
		<<nummuontracksintankpassedcut<<" had a track length > 50cm in the tank"<<endl;
	
	
#ifdef WCSIMDEBUG
	TCanvas* debug2 = new TCanvas("debug2","debug2");
	debug2->cd();
	tankstartvertices->SetMarkerColor(kRed);
	tankstartvertices->SetMarkerStyle(7);
	tankstartvertices->SetMarkerSize(1);
	tankstartvertices->Draw();
	vetostartvertices->SetMarkerColor(kBlue);
	vetostartvertices->SetMarkerStyle(7);
	vetostartvertices->SetMarkerSize(1);
	vetostartvertices->Draw("same");
	mrdstartvertices->SetMarkerColor(kViolet);
	mrdstartvertices->SetMarkerStyle(7);
	mrdstartvertices->SetMarkerSize(1);
	mrdstartvertices->Draw("same");
	debug2->SaveAs("start vertex distributions.png");
	delete tankstartvertices;
	delete vetostartvertices;
	delete mrdstartvertices;
	delete debug2;
	
	TCanvas* debug3 = new TCanvas("debug3","debug3");
	debug3->cd();
	tankstopvertices->SetMarkerColor(kRed);
	tankstopvertices->SetMarkerStyle(7);
	tankstopvertices->SetMarkerSize(1);
	tankstopvertices->Draw();
	vetostopvertices->SetMarkerColor(kBlue);
	vetostopvertices->SetMarkerStyle(7);
	vetostopvertices->SetMarkerSize(1);
	vetostopvertices->Draw("same");
	mrdstopvertices->SetMarkerColor(kViolet);
	mrdstopvertices->SetMarkerStyle(7);
	mrdstopvertices->SetMarkerSize(1);
	mrdstopvertices->Draw("same");
	debug3->SaveAs("stop vertex distributions.png");
	delete tankstopvertices;
	delete vetostopvertices;
	delete mrdstopvertices;
	delete debug3;
#endif  // WCSIMDEBUG
#endif  // DRAW_HISTOS

	if(flateventfileout){
		flateventfileout->cd();
		digitsqpmthist->Write();
		digitsqlappdhist->Write();
		digitstpmthist->Write();
		digitstlappdhist->Write();
		digittsmearpmthist->Write();
		digittsmearlappdhist->Write();
		lappdtimesmearvsqhist->Write();
		pmttimesmearvsqhist->Write();
	}
	
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
	
	// save, then delete canvases and legends
	cout<<"saving histograms to "<<outpath<<endl;
	//cout<<"deleting canvases and legends"<<endl;
	for(int i=0; i<canvaspointers.size(); i++){
		TCanvas* canv = canvaspointers.at(i);
		if(canv) canv->SaveAs(TString::Format("%s/histograms_%d.png",outpath,i));
		if(canv) delete canv; canv=0;
		TLegend* leg = legendpointers.at(i);
		if(leg) delete leg; leg=0;
	}
	
	// save, then delete scaled canvases and legends
	//cout<<"deleting scaled canvases and legends"<<endl;
	for(int i=0; i<scaledcanvaspointers.size(); i++){
		TCanvas* canv = scaledcanvaspointers.at(i);
		if(canv) canv->SaveAs(TString::Format("%s/scaled_histograms_%d.png",outpath,i));
		if(canv) delete canv; canv=0;
		TLegend* leg = scaledlegendpointers.at(i);
		if(leg) delete leg; leg=0;
	}
	
	// delete histograms - do this after the canvas version, which saves them before deleting them
	//cout<<"deleting histograms"<<endl;
	histofileout->cd();
	for(int i=0; i<histopointers.size(); i++){
		TH1D* temp = histopointers.at(i);
		if(temp) temp->Write();
		// some histograms are in the vector twice, so only delete it if it doesn't appear again
		if(std::count((histopointers.begin()+i), histopointers.end(), histopointers.at(i))==0) delete temp;
	}
	
	// delete scaled histograms
	//cout<<"deleting scaled histograms"<<endl;
	for(int i=0; i<scaledhistopointers.size(); i++){
		TH1D* temp = scaledhistopointers.at(i);
		if(temp) temp->Write(); 
		// some histograms are in the vector twice, so only delete it if it doesn't appear again
		if(std::count((scaledhistopointers.begin()+i), scaledhistopointers.end(), scaledhistopointers.at(i))==0) delete temp;
	}
	cout<<"end"<<endl;
	
	// Save all the remaining histograms
	std::vector<TH1*> otherhistos{
	(TH1*)neutrinovertex,
	(TH1*)neutrinovertexQE,
	(TH1*)neutrinovertexQEaccepted,
	(TH1*)muedepositionswcsim,
	(TH1*)muedepositionsfidcut,
	(TH1*)muedepositionsacceptedwcsim,
	(TH1*)fsltruetracklength,
	(TH1*)fsltruetracklengthintank,
#ifdef WCSIMDEBUG
	(TH1*)tankstartvertices,
	(TH1*)vetostartvertices,
	(TH1*)mrdstartvertices,
	(TH1*)tankstopvertices,
	(TH1*)vetostopvertices,
	(TH1*)mrdstopvertices,
#endif
	(TH1*)tracklengthvsmuonlight,
	(TH1*)chargemap_nopions,
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
	
	// delete flat file output
	//cout<<"deleting flat file output "<<flateventfileout<<endl;
	if(flateventfileout){
		flateventfileout->Close(); delete flateventfileout;
	}
	
	// write and close file of event information
	fileout->cd();
	treeout->SetEntries(bInTank->GetEntries());
	treeout->Write();
	fileout->Close();
	delete fileout;
	fileout=0;
}


/////////////////////////////
//    http://doxygen.genie-mc.org
//    if (genieint->ProcInfo().IsQuasiElastic() && genieint->ProcInfo().IsWeakCC()) { ... }
//    double Ethr = genieint->PhaseSpace().Threshold();
//    GHepParticle* probeneutrino = gevtRec->Probe();
//    GHepParticle* finalstatelepton = gevtRec->FinalStatePrimaryLepton();
//    GHepParticle* targetnucleon = gevtRec->HitNucleon();
//    PDG1000 nuclear codes: 10L[ZZZ][AAA]I
//    TString nuname = genieint->InitState().Probe()->GetName();
//    Int_t nupdg = gevtRec->Probe()->Pdg();
//    TLorentzVector* nuvtxval = gevtRec->Vertex();
//    nuvtx volumeName ?
//    nuvtx material ?
//    target energy?
/////////////////////////////
//    $GENIE/src/stdapp/gNtpConv.cxx

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
	/*
	getting the neut reaction code results in the printing of a bunch of surplus info, e.g:
1501283211 NOTICE GHepUtils : [n] <GHepUtils.cxx::NeutReactionCode (106)> : Current event is RES or DIS with W<2
1501283211 NOTICE GHepUtils : [n] <GHepUtils.cxx::NeutReactionCode (153)> : Num of primary particles: 
 p = 1, n = 0, pi+ = 0, pi- = 1, pi0 = 0, eta = 0, K+ = 0, K- = 0, K0 = 0, Lambda's = 0, gamma's = 0
	if we could redirect and capture this (rather than printing it to stdout) it might actually be useful,
	as extracting number of other hadrons doesn't work! but for now, just turn it off to reduce verbosity.
	*/
	/*Int_t*/ thegenieinfo.neutinteractioncode = -1; //genie::utils::ghep::NeutReactionCode(gevtRec);
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

float* Getlappdqpe(){
  static float qpe0[501]= {
    // 1
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
    0.000000, 0.000129, 0.000754, 0.004060, 0.028471,
    // 2
    0.068449, 0.115679, 0.164646, 0.203466, 0.235631,
    0.262351, 0.282064, 0.303341, 0.320618, 0.338317,
    0.357825, 0.371980, 0.385820, 0.398838, 0.413595,
    0.428590, 0.444387, 0.461685, 0.482383, 0.502369,
    0.520779, 0.540011, 0.559293, 0.579354, 0.599337,
    0.619580, 0.639859, 0.659807, 0.679810, 0.699620,
    0.718792, 0.737382, 0.755309, 0.772042, 0.788232,
    0.803316, 0.817861, 0.831148, 0.844339, 0.855532,
    0.866693, 0.876604, 0.886067, 0.894473, 0.902150,
    0.909515, 0.915983, 0.922050, 0.927418, 0.932492,
    // 3
    0.936951, 0.940941, 0.944660, 0.948004, 0.951090,
    0.953833, 0.956576, 0.958886, 0.961134, 0.963116,
    0.964930, 0.966562, 0.968008, 0.969424, 0.970687,
    0.971783, 0.972867, 0.973903, 0.974906, 0.975784,
    0.976632, 0.977438, 0.978190, 0.978891, 0.979543,
    0.980124, 0.980666, 0.981255, 0.981770, 0.982227,
    0.982701, 0.983146, 0.983566, 0.983975, 0.984357,
    0.984713, 0.985094, 0.985404, 0.985739, 0.986049,
    0.986339, 0.986630, 0.986922, 0.987176, 0.987431,
    0.987655, 0.987922, 0.988173, 0.988414, 0.988639,
    // 4
    0.988856, 0.989065, 0.989273, 0.989475, 0.989662,
    0.989828, 0.990007, 0.990172, 0.990327, 0.990497,
    0.990645, 0.990797, 0.990981, 0.991135, 0.991272,
    0.991413, 0.991550, 0.991673, 0.991805, 0.991928,
    0.992063, 0.992173, 0.992296, 0.992406, 0.992514,
    0.992632, 0.992733, 0.992837, 0.992954, 0.993046,
    0.993148, 0.993246, 0.993354, 0.993458, 0.993549,
    0.993656, 0.993744, 0.993836, 0.993936, 0.994033,
    0.994134, 0.994222, 0.994307, 0.994413, 0.994495,
    0.994572, 0.994659, 0.994739, 0.994816, 0.994886,
    // 5
    0.994970, 0.995032, 0.995110, 0.995178, 0.995250,
    0.995321, 0.995383, 0.995464, 0.995532, 0.995609,
    0.995674, 0.995750, 0.995821, 0.995889, 0.995952,
    0.996010, 0.996071, 0.996153, 0.996218, 0.996283,
    0.996335, 0.996384, 0.996431, 0.996484, 0.996537,
    0.996597, 0.996655, 0.996701, 0.996745, 0.996802,
    0.996860, 0.996917, 0.996962, 0.997014, 0.997079,
    0.997114, 0.997165, 0.997204, 0.997250, 0.997295,
    0.997335, 0.997379, 0.997418, 0.997454, 0.997488,
    0.997530, 0.997573, 0.997606, 0.997648, 0.997685,
    // 6
    0.997725, 0.997762, 0.997795, 0.997835, 0.997866,
    0.997898, 0.997941, 0.997966, 0.997997, 0.998039,
    0.998065, 0.998104, 0.998128, 0.998153, 0.998179,
    0.998205, 0.998223, 0.998254, 0.998293, 0.998319,
    0.998346, 0.998374, 0.998397, 0.998414, 0.998432,
    0.998456, 0.998482, 0.998511, 0.998532, 0.998553,
    0.998571, 0.998594, 0.998614, 0.998638, 0.998669,
    0.998693, 0.998715, 0.998743, 0.998762, 0.998793,
    0.998812, 0.998834, 0.998857, 0.998872, 0.998888,
    0.998904, 0.998926, 0.998946, 0.998963, 0.998983,
    // 7
    0.999007, 0.999027, 0.999044, 0.999064, 0.999079,
    0.999096, 0.999120, 0.999133, 0.999152, 0.999160,
    0.999174, 0.999188, 0.999206, 0.999221, 0.999234,
    0.999248, 0.999263, 0.999276, 0.999286, 0.999300,
    0.999313, 0.999321, 0.999331, 0.999347, 0.999356,
    0.999369, 0.999381, 0.999394, 0.999402, 0.999415,
    0.999427, 0.999433, 0.999446, 0.999458, 0.999472,
    0.999484, 0.999499, 0.999513, 0.999522, 0.999532,
    0.999540, 0.999550, 0.999559, 0.999567, 0.999574,
    0.999588, 0.999599, 0.999613, 0.999618, 0.999627,
    // 8
    0.999635, 0.999639, 0.999652, 0.999662, 0.999667,
    0.999671, 0.999678, 0.999682, 0.999688, 0.999693,
    0.999698, 0.999701, 0.999706, 0.999711, 0.999718,
    0.999722, 0.999727, 0.999732, 0.999737, 0.999740,
    0.999746, 0.999750, 0.999754, 0.999763, 0.999766,
    0.999769, 0.999774, 0.999780, 0.999784, 0.999788,
    0.999796, 0.999803, 0.999807, 0.999809, 0.999815,
    0.999820, 0.999827, 0.999830, 0.999833, 0.999833,
    0.999836, 0.999839, 0.999842, 0.999845, 0.999850,
    0.999853, 0.999857, 0.999860, 0.999865, 0.999870,
    // 9
    0.999873, 0.999877, 0.999880, 0.999882, 0.999883,
    0.999886, 0.999888, 0.999889, 0.999895, 0.999896,
    0.999897, 0.999901, 0.999902, 0.999905, 0.999907,
    0.999907, 0.999909, 0.999911, 0.999911, 0.999912,
    0.999913, 0.999914, 0.999917, 0.999919, 0.999921,
    0.999923, 0.999927, 0.999929, 0.999931, 0.999933,
    0.999936, 0.999942, 0.999942, 0.999944, 0.999947,
    0.999947, 0.999948, 0.999949, 0.999952, 0.999955,
    0.999957, 0.999957, 0.999961, 0.999962, 0.999963,
    0.999963, 0.999963, 0.999964, 0.999965, 0.999965,
    // 10
    0.999965, 0.999965, 0.999966, 0.999968, 0.999969,
    0.999971, 0.999972, 0.999972, 0.999973, 0.999975,
    0.999975, 0.999975, 0.999975, 0.999975, 0.999975,
    0.999975, 0.999979, 0.999979, 0.999980, 0.999982,
    0.999983, 0.999985, 0.999986, 0.999987, 0.999987,
    0.999988, 0.999989, 0.999989, 0.999989, 0.999989,
    0.999990, 0.999990, 0.999992, 0.999993, 0.999994,
    0.999994, 0.999994, 0.999994, 0.999994, 0.999995,
    0.999995, 0.999995, 0.999996, 0.999996, 0.999996,
    0.999996, 0.999998, 0.999999, 1.000000, 1.000000,
    // Dummy element for noticing if the loop reached the end of the array
    0.0 
  };
  return qpe0;
}

void SKIDigitizerThreshold(double& pe,int& iflag){
	double x = pe+0.1; iflag=0;
	double thr; double RDUMMY,err;
	if ( x<1.1) {
	  thr = std::min(1.0,
			 -0.06374+x*(3.748+x*(-63.23+x*(452.0+x*(-1449.0+x*(2513.0
									+x*(-2529.+x*(1472.0+x*(-452.2+x*(51.34+x*2.370))))))))));
	} else {
	  thr = 1.0;
	}
	RDUMMY = mrand->Rndm();
	if (thr < RDUMMY) {
	  pe = 0.0;
	  iflag = 1;
	}
	else {
	  err = mrand->Gaus(0.0,0.03);
	  /////      call rngaus(0.0, 0.03, err);
	  pe = pe+err;
	}
}
