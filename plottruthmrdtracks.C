/* vim:set noexpandtab tabstop=4 wrap */
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// XXX XXX FILE_VERSION of wcsim used to generate the file   XXX XXX
// XXX XXX          and Git source file COMMIT               XXX XXX
// XXX XXX          MUST BE SET BEFORE CALLING               XXX XXX
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define FILE_VERSION 3
/*
Version 1:
wcsim_wdirt_07-02-17, 200 PMTs + 200 LAPPDs, including dirt intx.

Version 2:
wcsim_tankonly_03-05-17, 200 PMTs + 200 LAPPDs, tank only, with bug fixes, bad lappd pulse timing resoln

Version 3:
wcsim_tankonly_17-06-17, 120 PMTs of 3 different types (LUX, Watchboy, LBNE, 8inHQE), no LAPPDs, endpoint energy and mom, no vector of names from WCSimRootGeom

Version 4:
wcsim ..., 200 PMTs + 200 LAPPDs, global position of LAPPD hits, vector of PMT names from WCSimRootGeom
*/

#ifndef USE_GRID
#define USE_GRID
#endif

// For wcsim processing of vincent's genie files, each genie file of 10k events had to be split up into 10 WCSim files.
// Since we have 10 wcsim files per dirt/genie file, we'll need to use a chain to load the wcsim files.
#ifndef SPLITWCSIMFILES
//#define SPLITWCSIMFILES
#endif

#ifndef NOLAPPDS
#define NOLAPPDS // no LAPPD files
#endif

#ifndef TANKONLY
#define TANKONLY // can only be disabled if it was disabled for WCSim as well
#endif

#ifndef VERBOSE
//#define VERBOSE
#endif

#ifndef TIME_EVENTS
//#define TIME_EVENTS
#endif

#ifndef WCSIMDEBUG
//#define WCSIMDEBUG
#endif

#ifndef MUTRACKDEBUG
//#define MUTRACKDEBUG
#endif

#ifndef MUTRACKLENGTHDEBUG
//#define MUTRACKLENGTHDEBUG
#endif

#ifndef PARTICLEGUNEVENTS
//#define PARTICLEGUNEVENTS // currently setup to match files of format wcsim_####.root and wcsim_lappd_####.root
//#define NOGENIE           // this needs to be changed in several places! search for "_####.root"
#endif

#ifndef NOGENIE
//#define NOGENIE
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
#include "TF1.h"
#include "TStopwatch.h"
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
#endif // FILE_VERSION>2

// genie headers
#ifndef NOGENIE
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
#endif //defined NOGENIE

void MakePMTmap(WCSimRootGeom* geo, std::map<int, std::pair<int,int> > &topcappositionmap, std::map<int, std::pair<int,int> > &bottomcappositionmap, std::map<int, std::pair<int,int> > &wallpositionmap);
#include "makepmtmaps_standalone.cxx"     // definition of this function

void ColourPlotStyle();
void FillTankMapHist(WCSimRootGeom* geo, int tubeID, int incone, std::map<std::string, TH2D*> &maphistos, double weight);
void ClearMapHistos(std::map<std::string,TH2D*> maphistos);	// clear the histograms

// helper functions to find the MRD intersection points
bool CheckLineBox( TVector3 L1, TVector3 L2, TVector3 B1, TVector3 B2, TVector3 &Hit, TVector3 &Hit2, bool &error);
int inline InBox( TVector3 Hit, TVector3 B1, TVector3 B2, const int Axis);
int inline GetIntersection( float fDst1, float fDst2, TVector3 P1, TVector3 P2, TVector3 &Hit);
// helper function for finding tank intersection points
bool CheckTankIntercepts( TVector3 primarystartvertex, TVector3 primarystopvertex, double avgtrackgradx, double avgtrackgrady, int trackstartvol, int trackstopvol, TVector3 &Hit, TVector3 &Hit2);
// unused helper functions
bool CheckLineCircle( TVector3 trackstart, TVector3 trackend, TVector3 tankorigin, double tankradius, TVector3 &Hit, TVector3 &Hit2, bool &error);
bool solveQuadratic(const double &a, const double &b, const double &c, double &x0, double &x1);

void PrintVector(TVector3 avec);

std::string PdgToString(int code);
std::map<int,std::string>* GeneratePdgMap();
static std::map<int,std::string> pdgcodetoname;

#ifndef NOGENIE
#include "genieinfo_struct.cpp"           // definition of a struct to hold genie info
// function to fill the into
void GetGenieEntryInfo(genie::EventRecord* gevtRec, genie::Interaction* genieint, GenieInfo& thegenieinfo, bool printneutrinoevent);
void PrintNeutrinoEvent(GenieInfo &thegenieinfo);
#endif // !defined NOGENIE

const Float_t MRD_width = (305./2.);      // half width of steel in cm
const Float_t MRD_height = (274./2.);     // half height of steel in cm
const Float_t MRD_layer2 = 290.755;       // position in wcsim coords of second scint layer in cm
const Float_t MRD_start = 325.5;          // position in wcsim coord of MRD front face in cm
const Float_t MRD_depth = 139.09;         // total depth of the MRD in cm
const Float_t MRD_end = MRD_start + MRD_depth;
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
Float_t mu_rest_mass_E = 105.658;

float* Getlappdqpe();
static void SKIDigitizerThreshold(double& pe,int& iflag);
TRandom3* mrand = new TRandom3();             // needed to generate charge and other things

// not currently used but these should be stored in and retrieved from geo.
Int_t numtankpmts=128+2*(26);                 // 26 pmts and lappds on each cap
Int_t nummrdpmts=307;
Int_t numvetopmts=26;
// these are used for making the pmt map.
const Int_t caparraysize=8;                   // pmts on the cap form an nxn grid where caparraysize=n
const Int_t pmtsperring=16;                   // pmts around each ring of the main walls
const Int_t numpmtrings=8;                    // num rings around the main walls

const Float_t tank_start = 15.70;             // front face of the tank in cm
const Float_t tank_radius = 152.4;            // tank radius in cm
const Float_t tank_halfheight = 198.;         // tank half height in cm
const Float_t tank_yoffset = -14.46;          // tank y offset in cm
/* from WCSimDetectorConfigs.cc
tankouterRadius= 1.524*m;
tankzoffset = 15.70*cm;
*/
const Float_t fidcutradius=tank_radius*0.8;   // fiducial volume is slightly lesss than the full tank radius
const Float_t fidcuty=50.;                    // a meter total fiducial volume in the centre y
const Float_t fidcutz=0;                      // fidcuial volume is before the tank centre.

// needed for drawing tank 2D map histograms
std::map<int, std::pair<int,int> > topcappositionmap;
std::map<int, std::pair<int,int> > bottomcappositionmap;
std::map<int, std::pair<int,int> > wallpositionmap;

#ifdef VERBOSE
const Bool_t printneutrinoevent=true;
#else 
const Bool_t printneutrinoevent=false;
#endif

void truthtracks(const char* wcsimpathin="", const char* dirtpathin="", const char* geniepathin="", const char* outpathin=""){
	
#if FILE_VERSION==0
	//const char awcsimpath[]="/pnfs/annie/persistent/users/moflaher/wcsim";  // first 1M sample, various issues
	const double triggeroffset=950;
	// is the implementation of this the same as VERSION_1? not sure...
#elif FILE_VERSION==1
	//const char awcsimpath[]="/pnfs/annie/persistent/users/moflaher/wcsim/wcsim_wdirt_07-02-17_rhatcher";
	const char awcsimpath[]="/pnfs/annie/persistent/users/moflaher/wcsim/lappd/wdirt/wcsim_lappd_wdirt_07-02-17_rhatcher";
	const double triggeroffset=950;
#elif FILE_VERSION==2 // lappd branch commit: 744169e224cc9d8573ea269132104c84b459466c
	//const char awcsimpath[]="/pnfs/annie/persistent/users/moflaher/wcsim/wcsim_tankonly_03-05-17_rhatcher";
	const char awcsimpath[]="/pnfs/annie/persistent/users/moflaher/wcsim/lappd/wcsim_lappd_tankonly_03-05-17_rhatcher";
	const double triggeroffset=950;
	
	//const char awcsimpath[]="/pnfs/annie/persistent/users/moflaher/wcsim/wcsim_tankonly_03-05-17_BNB_World_10k_29-06-17";
	//const char awcsimpath[]="/pnfs/annie/persistent/users/moflaher/wcsim/lappd/tankonly/wcsim_lappd_tankonly_03-05-17_BNB_World_10k_29-06-17";
#elif FILE_VERSION==3 // add multiple PMT types
	//const char awcsimpath[]="/pnfs/annie/persistent/users/moflaher/wcsim/wcsim_tankonly_17-06-17_rhatcher";
	const char awcsimpath[]="/pnfs/annie/persistent/users/moflaher/wcsim/multipmt/tankonly/wcsim_multipmt_tankonly_17-06-17_rhatcher";
	const double triggeroffset=950;
#elif FILE_VERSION==4 // add global position of LAPPD hits
	const char awcsimpath[]="/pnfs/annie/persistent/users/moflaher/wcsim/lappd/tankonly/wcsim_lappd_tankonly_24-09-17_BNB_Water_10k_22-05-17";
	const double triggeroffset=950;
#elif FILE_VERSION==5 
	// yet to be run
	const char awcsimpath="";
	const double triggeroffset=0; // set to 0 in fc66fe34ae9b49a72428e9f8d51de5ade93c2534, 26/09/17
	// fix PMT names in WCSimRootGeom... 
#else // if FILE_VERSION >5
	const char awcsimpath="";
	const double triggeroffset=0; // set to 0 in fc66fe34ae9b49a72428e9f8d51de5ade93c2534, 26/09/17
#endif // FILE_VERSION switch
	// generic paths not tied to a specific file version
	//const char awcsimpath[]="/annie/app/users/moflaher/wcsim/build";
	//std::string wcsimpathstlstring = std::string(gSystem->pwd())+"/in/muongun_beamsim";
	////const char awcsimpath[]=wcsimpathstlstring.c_str(); << doesn't work to create a char[].
	//char awcsimpath[150]; snprintf(&awcsimpath[0],150,"%s",wcsimpathstlstring.c_str());
	
	const char adirtpath[]="/pnfs/annie/persistent/users/moflaher/g4dirt_rhatcher";
	const char ageniepath[]="/pnfs/annie/persistent/users/rhatcher/genie";
	//const char adirtpath[]="/pnfs/annie/persistent/users/moflaher/g4dirt_vincentsgenie/world";
	//const char ageniepath[]="/pnfs/annie/persistent/users/vfischer/genie/BNB_World_10k_29-06-17";
	
	//std::string outpathstlstring=std::string(gSystem->pwd())+"/out/muongun_beamsim";
	std::string outpathstlstring=std::string(gSystem->pwd())+"/temp";
	char aoutpath[150]; snprintf(&aoutpath[0],150,"%s",outpathstlstring.c_str());
	//const char aoutpath[]=outpathstlstring.c_str(); << doesn't work to create a char[]. 
	//const char aoutpath[]=gSystem->pwd();
	//const char aoutpath[]="/annie/app/users/moflaher/wcsim/root_work";
	//const char aoutpath[]="/annie/app/users/moflaher/wcsim/root_work/temp";
	
	// ACTUAL PATHS SET HERE //
	// ===================== //
#ifndef USE_GRID
	const char* wcsimpath = (strcmp(wcsimpathin,"")!=0) ? wcsimpathin : &awcsimpath[0];
	const char* dirtpath = (strcmp(dirtpathin,"")!=0) ? dirtpathin : &adirtpath[0];
	const char* geniepath = (strcmp(geniepathin,"")!=0) ? geniepathin : &ageniepath[0];
	const char* outpath = (strcmp(outpathin,"")!=0) ? outpathin : &aoutpath[0];
#else // if defined USE_GRID
	// FIXME cannot mount root files from pnfs directly. Always use PWD
//	const char* wcsimpath = gSystem->Getenv("WCSIMDIR");
//	const char* dirtpath = gSystem->Getenv("DIRTDIR");
//	const char* geniepath = gSystem->Getenv("GENIEDIR");
	const char* PWD = gSystem->pwd();
	const char* wcsimpath = PWD;
	const char* dirtpath = PWD;
	const char* geniepath = PWD;
	const char* outpath = gSystem->Getenv("LOCALOUTDIR");
#endif // defined USE_GRID
	if(wcsimpath==nullptr||dirtpath==nullptr||geniepath==nullptr||outpath==nullptr||
		strcmp(wcsimpath,"")==0||strcmp(dirtpath,"")==0||strcmp(geniepath,"")==0||strcmp(outpath,"")==0){
		cerr<<"Path not set! wcsimpath="<<wcsimpath<<", dirtpath="<<dirtpath<<", geniepath="<<geniepath<<", outpath="<<outpath
			<<"FILE_VERSION="<<FILE_VERSION<<endl;
		assert(false);
	}
	
	ColourPlotStyle();
	
	// load WCSim library for reading WCSim files
	//const char* wcsimlibrarypath="/annie/app/users/moflaher/wcsim/wcsim/libWCSimRoot.so";
	std::string wcsimlibrarypath=std::string(gSystem->pwd())+"/../wcsim/libWCSimRoot.so";
	cout<<"loading "<<wcsimlibrarypath<<endl;
	gSystem->Load(wcsimlibrarypath.c_str());
#ifndef NOGENIE
	// load genie for reading genie files - done in a separate macro, before calling ACLiC on this file
//	TString script_dir = gSystem->Getenv("GENIE");
//	script_dir += "/src/scripts/gcint/";
//	TString curr_dir = gSystem->pwd();
//	gSystem->cd(script_dir.Data());
//	gROOT->ProcessLine(".x loadincs.C");
//	gROOT->ProcessLine(".x loadlibs.C");
//	gSystem->cd(curr_dir.Data());
#endif // !defined NOGENIE
	
	// ==============================================================================================
	// Declare input paths, files, trees...
	// ==============================================================================================
	TFile* wcsimfile=0;
	TString wcsimfilepath="n/a";
	TTree* wcsimT=0;
	TChain* wcsimchain=0;
	Int_t numwcsimentries=0;
	
	TFile* dirtfile=0;
	std::string dirtfilename="n/a";
	TTree* tankflux=0;
	TTree* tankmeta=0;
	
	TFile* lappdfile=0;
	TString lappdfilepath="n/a";
	TTree* lappdtree=0;
	TChain* lappdchain=0;
	Int_t numlappdentries=0;
	
	TFile* geniefile=0;
	TString geniefilepath="n/a";
	TTree* gtree=0;
	Int_t numgenietentries=0;
	
#ifndef PARTICLEGUNEVENTS
	// TChain for dirt files - this will be the main driver of the loop - all it's events will be processed.
	TChain* c =  new TChain("tankflux");
#ifndef USE_GRID
	TString chainpattern = TString::Format("%s/annie_tank_flux.*.root",dirtpath);	//1000->*
#else // if defined USE_GRID
// looks like you cannot mount root files from pnfs directly. 
//	const char* fnumchars = gSystem->Getenv("FILENUM");
//	if(fnumchars==nullptr||strcmp(fnumchars,"")==0)
//		{cerr<<"FILENUM environmental variable not defined!!"<<endl; assert(false);}
//	TString chainpattern = TString::Format("%s/annie_tank_flux.%s.root",dirtpath,fnumchars);
	TString chainpattern = TString::Format("%s/annie_tank_flux.*.root",PWD);
#endif // defined USE_GRID
	cout<<"loading TChain entries from "<<chainpattern<<endl;
	c->Add(chainpattern);
	Int_t numents = c->GetEntries();
	cout<<"loaded "<<numents<<" entries in the chain"<<endl;
	if(numents<1){ return; }
#else // if defined PARTICLEGUNEVENTS
	// TChain for wcsim files - this will be the main driver of the loop - all it's events will be processed.
	TChain* c =  new TChain("wcsimT");
	//TString chainpattern = TString::Format("%s/wcsim_0.*.root",wcsimpath);        // for files "wcsim_0.####.root"
	//TString chainpattern = TString::Format("%s/wcsim_*([0-9]).root",wcsimpath);   // matches below in bash but not root
	TString chainpattern = TString::Format("%s/wcsim_[0-9]+.root",wcsimpath);       // for files "wcsim_####.root"
	cout<<"loading TChain entries from "<<chainpattern<<endl;
	c->Add(chainpattern);
	Int_t numents = c->GetEntries();
	cout<<"loaded "<<numents<<" entries in the chain"<<endl;
	if(numents<1){ return; }
#endif // defined PARTICLEGUNEVENTS
	
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
	
#ifndef NOGENIE
	// gtree
	cout<<"creating new genie::NtpMCEventRecord"<<endl;
	genie::NtpMCEventRecord* genierecordval = new genie::NtpMCEventRecord;
	TBranch* genierecordBranch=0;
#endif // !defined NOGENIE
	
	// wcsimT
	WCSimRootEvent* b=0, *m=0, *v=0;
	TBranch* bp=0, *mp=0, *vp=0;
	WCSimRootTrigger* atrigt=0, *atrigm=0, *atrigv=0;
	
	// lappdtree
	// some of these are c-style arrays of all hits in the event, some are vectors of same.
	// for c-style arrays, avoid re-allocations of memory by just setting one array big enough for all events
	int lappd_evtnum;
	int lappd_numtileshitthisevt;
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
#endif // FILE_VERSION>3
	
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
	
#ifndef NOGENIE
	// information from genie:
	GenieInfo thegenieinfo;
#endif // !defined NOGENIE
	
	// ==============================================================================================
	// EventDistributions output file
	// ==============================================================================================
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
	int wcsimtriggernum=-1;
	TBranch* bWCSimTriggerNum = treeout->Branch("WCSimTriggerNum",&wcsimtriggernum);
	// now information about the event
	TLorentzVector interactionvertex(0.,0.,0.,0.);
	TBranch* bInteractionVtx = treeout->Branch("NeutrinoVertex",&interactionvertex);
	std::map<std::string,bool> eventtypes;
	std::map<std::string,bool>* eventtypesp=&eventtypes;
	TBranch* bEventType = treeout->Branch("TypesMap",&eventtypesp);
	// using a map is a bad idea! It means you can't use the event type in tree->Draw calls!!!
	// NOR CAN YOU EVEN DO tree->Show() or tree/Branch->GetEntry() calls. 
	// To be able to do these, you first need to load a CollectionProxy by putting the following code
	// into a file and calling '.L thefile.C+' 
	//		#include <map>
	//		#ifdef __MAKECINT__
	//		#pragma link C++ class std::map<std::string,bool>+;
	//		#endif
	// Instead, use the genie method, and save a bunch of bools, and a string
	std::string interactiontypestring;
	TBranch* bIntxTypeString = treeout->Branch("InteractionTypeString",&interactiontypestring);
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
	// ok, moving on - remainder of neutrino interaction information
	bool isintank=false;
	TBranch* bInTank = treeout->Branch("NuVtxInTank",&isintank);
	bool isinfiducialvol=false;
	TBranch* bInFidVol = treeout->Branch("NuVtxInFidVol",&isinfiducialvol);
	double eventq2=-1;
	TBranch* bEventQ2 = treeout->Branch("EventQ2",&eventq2);
	double eventEnu=-1;
	TBranch* bEventEnu = treeout->Branch("NeutrinoEnergy",&eventEnu);
	int neutrinopdg=-1;
	TBranch* bNeutrinoPdg = treeout->Branch("NeutrinoPDG",&neutrinopdg);
	double geniefsle=-1;
	TBranch* bGenieFslE = treeout->Branch("GenieFslE",&geniefsle);
	double geniefslangle=-1;
	TBranch* bGenieFslAng = treeout->Branch("GenieFslAngle",&geniefslangle);
	// summary information about the event in the detector
	// TODO: add LAPPD hit info
	int numTankDigits=-1;
	TBranch* bNumTankDigits = treeout->Branch("TotalTankDigits",&numTankDigits);
	double totaltankcharge=-1;
	TBranch* bTotalTankCharge = treeout->Branch("TotalTankCharge",&totaltankcharge);
	double upstreamcharge=-1;
	TBranch* bTotalUpstreamCharge = treeout->Branch("TotalUpstreamCharge",&upstreamcharge);
	double downstreamcharge=-1;
	TBranch* bTotalDownstreamCharge = treeout->Branch("TotalDownstreamCharge",&downstreamcharge);
	double topcapcharge=-1;
	TBranch* bTopCapCharge = treeout->Branch("TotalTopCapCharge",&topcapcharge);
	double bottomcapcharge=-1;
	TBranch* bBottomCapCharge = treeout->Branch("TotalBottomCapCharge",&bottomcapcharge);
	int numMRDdigits=-1;
	TBranch* bTotalMrdDigits = treeout->Branch("TotalMrdDigits",&numMRDdigits);
	double totMRDcharge=-1;
	TBranch* bTotalMrdCharge = treeout->Branch("TotalMrdCharge",&totMRDcharge);
	int numFACCdigits=-1;
	TBranch* bTotalVetoDigits = treeout->Branch("TotalVetoDigits",&numFACCdigits);
	double totFACCcharge=-1;
	TBranch* bTotalVetoCharge = treeout->Branch("TotalVetoCharge",&totFACCcharge);
	// store information about the other primary tracks in the event.
	// this may be useful for extracting events that have something else going on, and what
	// count what other tracks were present
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
	// also record the number of neutron captures
	int numneutroncaptures=-1;
	TBranch* bNumNeutronCaptures = treeout->Branch("NumNeutronCaptures",&numneutroncaptures);
	// pointer to the primary muon track, if present
	int primarymuonindex=-1;
	TBranch* bPrimaryMuonIndex = treeout->Branch("PrimaryMuonIndex",&primarymuonindex);
	int primaryneutrinoindex=-1;
	TBranch* bPrimaryNeutrinoIndex = treeout->Branch("PrimaryNeutrinoIndex",&primaryneutrinoindex);
	// pointer to the highest energy muon... in case not, or in case it's different?
	int highestenergymuonindex=-1;
	TBranch* bHEMuon = treeout->Branch("HighestEnergyMuIndex",&highestenergymuonindex);
	// vectors of information about the tracks
	std::map<std::string,int> particlestosave;
	// THESE STRINGS MUST MATCH THOSE IN PdgToString
	particlestosave.emplace("Pion-",0);
	particlestosave.emplace("Pion+",0);
	particlestosave.emplace("Pion0",0);
	particlestosave.emplace("Muon-",0);
	particlestosave.emplace("Muon+",0);
	particlestosave.emplace("Gamma",0);
	particlestosave.emplace("Neutron",0);
	particlestosave.emplace("Proton",0);
	std::map<std::string,int>* particlestosavep = &particlestosave;
	TBranch* bParticlesToSave = treeout->Branch("SavedParticleTypes",&particlestosavep); // constant entry
	std::vector<int> trackpdg;
	std::vector<int>* trackpdgp = &trackpdg;
	TBranch* bTrackPDG = treeout->Branch("TrackPDG",&trackpdgp);
	std::vector<TVector3> trackstartpos;
	std::vector<TVector3>* trackstartposp = &trackstartpos;
	TBranch* bTrackStartPos = treeout->Branch("TrackStartPos",&trackstartposp);
	std::vector<int> trackstartvol;
	std::vector<int>* trackstartvolp = &trackstartvol;
	TBranch* bTrackStartVol = treeout->Branch("TrackStartVol",&trackstartvolp);
	std::vector<double> trackstarttime;
	std::vector<double>* trackstarttimep = &trackstarttime;
	TBranch* bTrackStartTime = treeout->Branch("TrackStartTime",&trackstarttimep);
	std::vector<double> trackstartE;
	std::vector<double>* trackstartEp = &trackstartE;
	TBranch* bTrackStartE = treeout->Branch("TrackStartE",&trackstartEp);
	std::vector<TVector3> trackstartmom;
	std::vector<TVector3>* trackstartmomp = &trackstartmom;
	TBranch* bTrackStartMom = treeout->Branch("TrackStartMom",&trackstartmomp);
	std::vector<TVector3> trackstoppos;
	std::vector<TVector3>* trackstopposp = &trackstoppos;
	TBranch* bTrackStopPos = treeout->Branch("TrackStopPos",&trackstopposp);
	std::vector<int> trackstopvol;
	std::vector<int>* trackstopvolp = &trackstopvol;
	TBranch* bTrackStopVol = treeout->Branch("TrackStopVol",&trackstopvolp);
	std::vector<double> trackstoptime;
	std::vector<double>* trackstoptimep = &trackstoptime;
	TBranch* bTrackStopTime = treeout->Branch("TrackStopTime",&trackstoptimep);
	std::vector<double> trackstopE;
	std::vector<double>* trackstopEp=&trackstopE;
	TBranch* bTrackStopE = treeout->Branch("TrackStopE",&trackstopEp);
	std::vector<TVector3> trackstopmom;
	std::vector<TVector3>* trackstopmomp = &trackstopmom;
	TBranch* bTrackStopMom = treeout->Branch("TrackStopMom",&trackstopmomp);
	std::vector<double> trackangle; // redundant... but maybe useful for MRD E
	std::vector<double>* trackanglep = &trackangle;
	TBranch* bTrackAngle = treeout->Branch("TrackAvgAngle",&trackanglep);
	std::vector<int> trackparenttype;
	std::vector<int>* trackparenttypep = &trackparenttype;
	TBranch* bTrackParentType = treeout->Branch("TrackParentType",&trackparenttypep);
	// detector response to this particle
	std::vector<int> numtankdigitsfromtrack;
	std::vector<int>* numtankdigitsfromtrackp = &numtankdigitsfromtrack;
	TBranch* bTrackNumTankDigits = treeout->Branch("TankDigitsFromTrack",&numtankdigitsfromtrackp);
	std::vector<double> tankchargefromtrack;
	std::vector<double>* tankchargefromtrackp = &tankchargefromtrack;
	TBranch* bTrackTankCharge = treeout->Branch("TankChargeFromTrack",&tankchargefromtrackp);
	std::vector<double> fractionaltrackchargeincone;
	std::vector<double>* fractionaltrackchargeinconep = &fractionaltrackchargeincone;
	TBranch* bTrackFractionOfChargeInCone = treeout->Branch("FractionOfTrackChargeInCone",&fractionaltrackchargeinconep);
	std::vector<double> upstreamchargefromtrack;
	std::vector<double>* upstreamchargefromtrackp = &upstreamchargefromtrack;
	TBranch* bTrackUpstreamCharge = treeout->Branch("UpstreamChargeFromTrack",&upstreamchargefromtrackp);
	std::vector<double> downstreamchargefromtrack;
	std::vector<double>* downstreamchargefromtrackp = &downstreamchargefromtrack;
	TBranch* bTrackDownstreamCharge = treeout->Branch("DownstreamChargeFromTrack",&downstreamchargefromtrackp);
	std::vector<double> topcapchargefromtrack;
	std::vector<double>* topcapchargefromtrackp = &topcapchargefromtrack;
	TBranch* bTrackTopCapCharge = treeout->Branch("TopCapChargeFromTrack",&topcapchargefromtrackp);
	std::vector<double> bottomcapchargefromtrack;
	std::vector<double>* bottomcapchargefromtrackp = &bottomcapchargefromtrack;
	TBranch* bTrackBottomCapCharge = treeout->Branch("BottomCapChargeFromTrack",&bottomcapchargefromtrackp);
	std::vector<std::vector<int>> tanktubeshitbytrack;
	std::vector<std::vector<int>>* tanktubeshitbytrackp = &tanktubeshitbytrack;
	TBranch* bTrackTankTubesHitByTrack = treeout->Branch("TankTubesHitByTrack",&tanktubeshitbytrackp);
	std::vector<double> tracklengthintank;
	std::vector<double>* tracklengthintankp = &tracklengthintank;
	TBranch* bTrackTrackLengthInTank = treeout->Branch("TrackLengthInTank",&tracklengthintankp);
	std::vector<bool> trackentersMRD;
	std::vector<bool>* trackentersMRDp = &trackentersMRD;
	TBranch* bTracktrackEntersMRD = treeout->Branch("TrackEntersMRD",&trackentersMRDp);
	std::vector<bool> trackstopsinMRD;
	std::vector<bool>* trackstopsinMRDp = &trackstopsinMRD;
	TBranch* bTracktrackStopsInMRD = treeout->Branch("TrackStopsInMRD",&trackstopsinMRDp);
	std::vector<bool> trackpenetratesMRD;
	std::vector<bool>* trackpenetratesMRDp = &trackpenetratesMRD;
	TBranch* bTracktrackRangesOutMRD = treeout->Branch("TrackRangesOutMRD",&trackpenetratesMRDp);
	std::vector<double> trackmrdpenetrationcm;
	std::vector<double>* trackmrdpenetrationcmp = &trackmrdpenetrationcm;
	TBranch* bTracktrackMrdPenetrationInCm = treeout->Branch("TrackMrdPenetrationInCm",&trackmrdpenetrationcmp);
	std::vector<int> trackmrdpenetrationlayers;
	std::vector<int>* trackmrdpenetrationlayersp = &trackmrdpenetrationlayers;
	TBranch* bTracktrackMrdPenetrationLayers = treeout->Branch("TrackMrdPenetrationLayers",&trackmrdpenetrationlayersp);
	std::vector<double> tracklengthinMRD;
	std::vector<double>* tracklengthinMRDp = &tracklengthinMRD;
	TBranch* bTrackTrackLengthInMRD = treeout->Branch("TrackLengthInMRD",&tracklengthinMRDp);
	std::vector<int> numMRDdigitsfromtrack;
	std::vector<int>* numMRDdigitsfromtrackp = &numMRDdigitsfromtrack;
	TBranch* bTrackMrdDigitsFromTrack = treeout->Branch("MrdDigitsFromTrack",&numMRDdigitsfromtrackp);
	std::vector<double> MRDchargefromtrack;
	std::vector<double>* MRDchargefromtrackp = &MRDchargefromtrack;
	TBranch* bTrackMrdChargeFromTrack = treeout->Branch("MrdChargeFromTrack",&MRDchargefromtrackp);
	std::vector<std::vector<int>> mrdtubeshitbytrack;
	std::vector<std::vector<int>>* mrdtubeshitbytrackp = &mrdtubeshitbytrack;
	TBranch* bTrackMrdTubesHitByMuon = treeout->Branch("MrdTubesHitByTrack",&mrdtubeshitbytrackp);
	std::vector<int> numFACCdigitsfromtrack;
	std::vector<int>* numFACCdigitsfromtrackp = &numFACCdigitsfromtrack;
	TBranch* bTrackVetoDigitsFromTrack = treeout->Branch("VetoDigitsFromTrack",&numFACCdigitsfromtrackp);
	std::vector<double> FACCchargefromtrack;
	std::vector<double>* FACCchargefromtrackp = &FACCchargefromtrack;
	TBranch* bTrackVetoChargeFromTrack = treeout->Branch("VetoChargeFromTrack",&FACCchargefromtrackp);
	std::vector<std::vector<int>> facctubeshitbytrack;
	std::vector<std::vector<int>>* facctubeshitbytrackp = &facctubeshitbytrack;
	TBranch* bTrackVetoTubesHitByMuon = treeout->Branch("VetoTubesHitByTrack",&facctubeshitbytrackp);
	
//	std::vector<TBranch*> thebranches{ bInTank, bEventType, bEventQ2, bEventEnu, bNeutrinoPdg, bInFidVol, bMuonStartVtx, bMuonStopVtx, bMuonStartE, bMuonTrackLengthInTank, bMuonMrdPenetrationInCm, bMuonMrdPenetrationLayers, bMuonEntersMRD, bMuonStopsInMRD, bMuonRangesOutMRD, bMuonTrackLengthInMRD, bTankChargeFromMuon, bFractionOfMuonChargeInCone};
//	int someit=0;
//	bool haszombies=false;
//	for(auto abranch : thebranches){
//		if(abranch==0){ cout<<"branch "<<someit<<" is a zombie"<<endl; haszombies=true; }
//		someit++;
//	}
//	assert(!haszombies&&"output file branches have zombies");
	
	// ==============================================================================================
	// flat file for true vertices and digits for tank reconstruction efforts
	// ==============================================================================================
	gROOT->cd();
	TFile* flateventfileout = new TFile(TString::Format("%s/trueQEvertexinfo.root",outpath), "RECREATE");
	flateventfileout->cd();
	double fileeventnum=-1;
	double filetriggernum=-1;
	double fileneutrinoE=-1;
	string fileinteractiontypestring="";
	int fileneutcode=-1;
	double filemomtrans=-1;
	double filemuonenergy=-1;
	double filemuonangle=-1;
	double filepathlengthtotal=-1;
	double filepathlengthinwater=-1;
	double filepathlengthinmrd=-1;
	double fileenergylossinmrd=-1;
	TLorentzVector filemuonstartvertex(0.,0.,0.,0.);
	TLorentzVector filemuonstopvertex(0.,0.,0.,0.);
	TVector3 filemuondirectionvector(0.,0.,0.);
	std::vector<ROOT::Math::XYZTVector>  filedigitvertices;
	std::vector<ROOT::Math::XYZTVector>* filedigitverticesp = &filedigitvertices;
	std::vector<Double_t> filedigitQs;
	std::vector<Double_t>* filedigitQsp = &filedigitQs;
	std::vector<std::string> filedigitsensortypes;
	std::vector<std::string>* filedigitsensortypesp=&filedigitsensortypes;
	std::vector<Double_t> filedigittsmears;
	std::vector<Double_t>* filedigittsmearsp=&filedigittsmears;
	std::vector<int> filedigitPMTIDs;
	std::vector<int>* filedigitPMTIDsp=&filedigitPMTIDs;
	TTree* vertextreenocuts = new TTree("vertextreenocuts","All True Tank QE Events");
	TBranch* FileEventNumBranch = vertextreenocuts->Branch("EventNum",&fileeventnum);
	TBranch* FileTriggerNumBranch = vertextreenocuts->Branch("SubTriggerNum",&filetriggernum);
	TBranch* NeutrinoEnergyBranch = vertextreenocuts->Branch("NeutrinoEnergy",&fileneutrinoE);
	TBranch* InteractionTypeBranch = vertextreenocuts->Branch("InteractionType",&fileinteractiontypestring);
	TBranch* NeutCodeBranch = vertextreenocuts->Branch("NeutCode",&fileneutcode);
	TBranch* MomTransBranch = vertextreenocuts->Branch("MomentumTransfer",&filemomtrans);
	TBranch* MuonEnergyBranch = vertextreenocuts->Branch("MuonEnergy",&filemuonenergy);
	TBranch* MuonAngleBranch = vertextreenocuts->Branch("MuonAngle",&filemuonangle);
	TBranch* TotalTrackLengthBranch = vertextreenocuts->Branch("TotalTrackLength",&filepathlengthtotal);
	TBranch* TrackLengthInWaterBranch = vertextreenocuts->Branch("TrackLengthInWater",&filepathlengthinwater);
	TBranch* TrackLengthInMrdBranch = vertextreenocuts->Branch("TrackLengthInMrd",&filepathlengthinmrd);
#ifdef MUTRACKLENGTHDEBUG
	double filemuxtracklength, filemuytracklength, filemuztracklength;
	TBranch* TrackLengthInMrdXBranch = vertextreenocuts->Branch("TrackLengthInMrdX",&filemuxtracklength);
	TBranch* TrackLengthInMrdYBranch = vertextreenocuts->Branch("TrackLengthInMrdY",&filemuytracklength);
	TBranch* TrackLengthInMrdZBranch = vertextreenocuts->Branch("TrackLengthInMrdZ",&filemuztracklength);
#endif // defined MUTRACKLENGTHDEBUG
	TBranch* EnergyLossInMrdBranch =  vertextreenocuts->Branch("EnergyLossInMrd",&fileenergylossinmrd);
	TBranch* MuonStartBranch = vertextreenocuts->Branch("MuonStartVertex",&filemuonstartvertex);
	TBranch* MuonStopBranch = vertextreenocuts->Branch("MuonStopVertex", &filemuonstopvertex);
	TBranch* MuonDirectionBranch = vertextreenocuts->Branch("MuonDirection", &filemuondirectionvector);
	TBranch* DigitVertexBranch = vertextreenocuts->Branch("DigitVertices", "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >", &filedigitverticesp);
	TBranch* DigitChargeBranch = vertextreenocuts->Branch("DigitCharges", &filedigitQsp);
	TBranch* DigitDetTypeBranch = vertextreenocuts->Branch("DigitWhichDet",&filedigitsensortypesp);
	TBranch* DigitSmearBranch = vertextreenocuts->Branch("DigitTimeSmear",&filedigittsmearsp);
	TBranch* DigitPmtIdBranch = vertextreenocuts->Branch("DigitPmtId",&filedigitPMTIDsp);
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
#endif // FILE_VERSION>3
	std::vector<int>* tileorientp=&tileorient;
	std::vector<int>* octagonsidep=&octagonside;
	TBranch* LAPPD_intileposx = vertextreenocuts->Branch("LAPPD_intileposx",&intileposxp);
	TBranch* LAPPD_intileposy = vertextreenocuts->Branch("LAPPD_intileposy",&intileposyp);
#if FILE_VERSION>3
	TBranch* LAPPD_poserrx = vertextreenocuts->Branch("LAPPD_poserrx",&poserrxp);
	TBranch* LAPPD_poserry = vertextreenocuts->Branch("LAPPD_poserry",&poserryp);
	TBranch* LAPPD_poserrz = vertextreenocuts->Branch("LAPPD_poserrz",&poserrzp);
#endif // FILE_VERSION>3
	TBranch* LAPPD_tileorient = vertextreenocuts->Branch("LAPPD_tileorient",&tileorientp);
	TBranch* LAPPD_octagonside = vertextreenocuts->Branch("LAPPD_octagonside",&octagonsidep);
#endif // defined LAPPD_DEBUG
	
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
	TBranch* FileEventNumBranchFid = vertextreefiducialcut->Branch("EventNum",&fileeventnum);
	TBranch* FileTriggerNumBranchFid = vertextreefiducialcut->Branch("SubTriggerNum",&filetriggernum);
	TBranch* NeutrinoEnergyBranchFid = vertextreefiducialcut->Branch("NeutrinoEnergy",&fileneutrinoE);
	TBranch* InteractionTypeBranchFid = vertextreefiducialcut->Branch("InteractionType",&fileinteractiontypestring);
	TBranch* NeutCodeBranchFid = vertextreefiducialcut->Branch("NeutCode",&fileneutcode);
	TBranch* MomTransBranchFid = vertextreefiducialcut->Branch("MomentumTransfer",&filemomtrans);
	TBranch* MuonEnergyBranchFid = vertextreefiducialcut->Branch("MuonEnergy",&filemuonenergy);
	TBranch* MuonAngleBranchFid = vertextreefiducialcut->Branch("MuonAngle",&filemuonangle);
	TBranch* TotalTrackLengthBranchFid = vertextreefiducialcut->Branch("TotalTrackLength",&filepathlengthtotal);
	TBranch* TrackLengthInWaterBranchFid = vertextreefiducialcut->Branch("TrackLengthInWater",&filepathlengthinwater);
	TBranch* TrackLengthInMrdBranchFid = vertextreefiducialcut->Branch("TrackLengthInMrd",&filepathlengthinmrd);
	TBranch* EnergyLossInMrdBranchFid =  vertextreefiducialcut->Branch("EnergyLossInMrd",&fileenergylossinmrd);
	TBranch* MuonStartBranchFid = vertextreefiducialcut->Branch("MuonStartVertex",&filemuonstartvertex);
	TBranch* MuonStopBranchFid = vertextreefiducialcut->Branch("MuonStopVertex", &filemuonstopvertex);
	TBranch* MuonDirectionBranchFid = vertextreefiducialcut->Branch("MuonDirection", &filemuondirectionvector);
	TBranch* DigitVertexBranchFid = vertextreefiducialcut->Branch("DigitVertices", "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >", &filedigitverticesp);
	TBranch* DigitChargeBranchFid = vertextreefiducialcut->Branch("DigitCharges", &filedigitQsp);
	TBranch* DigitDetTypeBranchFid = vertextreefiducialcut->Branch("DigitWhichDet",&filedigitsensortypesp);
	TBranch* DigitSmearBranchFid = vertextreefiducialcut->Branch("DigitTimeSmear",&filedigittsmearsp);
	TBranch* DigitPmtIdBranchFid = vertextreefiducialcut->Branch("DigitPmtId",&filedigitPMTIDsp);
	
	TTree* vertextreefiducialmrd = new TTree("vertextreefiducialmrd","True Tank QE Events in Fiducial Volume With Muon in MRD");
	TBranch* FileEventNumBranchFidMRD = vertextreefiducialmrd->Branch("EventNum",&fileeventnum);
	TBranch* FileTriggerNumBranchFidMRD = vertextreefiducialmrd->Branch("SubTriggerNum",&filetriggernum);
	TBranch* NeutrinoEnergyBranchFidMRD = vertextreefiducialmrd->Branch("NeutrinoEnergy",&fileneutrinoE);
	TBranch* InteractionTypeBranchFidMRD = vertextreefiducialmrd->Branch("InteractionType",&fileinteractiontypestring);
	TBranch* NeutCodeBranchFidMRD = vertextreefiducialmrd->Branch("NeutCode",&fileneutcode);
	TBranch* MomTransBranchFidMRD = vertextreefiducialmrd->Branch("MomentumTransfer",&filemomtrans);
	TBranch* MuonEnergyBranchFidMRD = vertextreefiducialmrd->Branch("MuonEnergy",&filemuonenergy);
	TBranch* MuonAngleBranchFidMRD = vertextreefiducialmrd->Branch("MuonAngle",&filemuonangle);
	TBranch* TotalTrackLengthBranchFidMRD = vertextreefiducialmrd->Branch("TotalTrackLength",&filepathlengthtotal);
	TBranch* TrackLengthInWaterBranchFidMRD = vertextreefiducialmrd->Branch("TrackLengthInWater",&filepathlengthinwater);
	TBranch* TrackLengthInMrdBranchFidMRD = vertextreefiducialmrd->Branch("TrackLengthInMrd",&filepathlengthinmrd);
	TBranch* EnergyLossInMrdBranchFidMRD =  vertextreefiducialmrd->Branch("EnergyLossInMrd",&fileenergylossinmrd);
	TBranch* MuonStartBranchFidMRD = vertextreefiducialmrd->Branch("MuonStartVertex",&filemuonstartvertex);
	TBranch* MuonStopBranchFidMRD = vertextreefiducialmrd->Branch("MuonStopVertex", &filemuonstopvertex);
	TBranch* MuonDirectionBranchFidMRD = vertextreefiducialmrd->Branch("MuonDirection", &filemuondirectionvector);
	TBranch* DigitVertexBranchFidMRD = vertextreefiducialmrd->Branch("DigitVertices", "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >", &filedigitverticesp);
	TBranch* DigitChargeBranchFidMRD = vertextreefiducialmrd->Branch("DigitCharges", &filedigitQsp);
	TBranch* DigitDetTypeBranchFidMRD = vertextreefiducialmrd->Branch("DigitWhichDet",&filedigitsensortypesp);
	TBranch* DigitSmearBranchFidMRD = vertextreefiducialmrd->Branch("DigitTimeSmear",&filedigittsmearsp);
	TBranch* DigitPmtIdBranchFidMRD = vertextreefiducialmrd->Branch("DigitPmtId",&filedigitPMTIDsp);
	
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
	// histograms file
	// ===================================================================================================
	gROOT->cd();
	TFile* histofileout = new TFile(TString::Format("%s/TruthHistos.root",outpath),"RECREATE");
	
	// genie histograms
	TH1D* incidentneutrinoenergiesall = new TH1D("incidentneutrinoenergiesall","Distribution of Probe Neutrino Energies;Energy (GeV);Num Events",100,0,0.);
	TH1D* incidentneutrinoenergiesaccepted = new TH1D("incidentneutrinoenergiesaccepted","Distribution of Accepted Probe Neutrino Energies;Energy (GeV);Num Events",100,0,0.);
	TH1D* incidentneutrinoanglesall = new TH1D("incidentneutrinoanglesall","Distribution of Probe Neutrino Angles from z;Angle (rads);Num Events",100,-TMath::Pi(),TMath::Pi());
	TH1D* incidentneutrinoanglesaccepted = new TH1D("incidentneutrinoanglesaccepted","Distribution of Accepted Probe Neutrino Angles from z;Angle (rads);Num Events",100,-TMath::Pi(),TMath::Pi());
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
	
	// WCSim ones
	TH1D* incidentneutrinoenergiesacceptedwcsimfidcut = new TH1D("incidentneutrinoenergiesacceptedwcsimfidcut","Distribution of Accepted Probe Neutrino Energies;Energy (GeV);Num Events",100,0,0.);
	TH1D* fslanglesacceptedwcsimfidcut = new TH1D("fslanglesacceptedwcsimfidcut","Distribution of Accepted Final State Lepton Angles;Angle (rads);Num Events",100,0.,TMath::Pi());
	TH1D* fslenergiesacceptedwcsimfidcut = new TH1D("fslenergiesacceptedwcsimfidcut","Distribution of Accepted Final State Lepton Energies;Energy (GeV);Num Events",100,0.,3.);
	TH1D* eventq2acceptedwcsimfidcut = new TH1D("eventq2acceptedwcsimfidcut","Distribution of Accepted Event Q^2 Values;Q^2 (GeV/c)^2;Num Events",100,0.,3.);
	
	// should be fixed in later file versions...?
	TH1D* muedepositionsacceptedwcsim = new TH1D("muedepositionsacceptedwcsim","Distribution of Muon Energy Depositions In Tank, with MRD Selection;Energy (PMT Q);Num Events",100,0.,100);
	
//	// debugging
	TH1D* fsltruetracklength = new TH1D("fsltruetracklength", "Distribution of True Track Lengths", 100, 0., 1500.);
	TH1D* fsltruetracklengthintank = new TH1D("fsltruetracklengthintank", "Distribution of True Track Lengths In Tank", 100, 0., 1500.);
#ifdef WCSIMDEBUG
	TH3D* tankstartvertices = new TH3D("tankstartvertices", "Distribution of Tank Starting Vertices", 100, -150.,150., 100, -150., 150., 100, -100., 800.);
	TH3D* vetostartvertices = new TH3D("vetostartvertices", "Distribution of Veto Starting Vertices", 100, -150.,150., 100, -150., 150., 100, -100., 800.);
	TH3D* mrdstartvertices = new TH3D("mrdstartvertices", "Distribution of MRD Starting Vertices", 100, -150.,150., 100, -150., 150., 100, -100., 800.);
	TH3D* tankstopvertices = new TH3D("tankstopvertices", "Distribution of Tank Stopping Vertices", 100, -150.,150., 100, -150., 150., 100, -100., 800.);
	TH3D* vetostopvertices = new TH3D("vetostopvertices", "Distribution of Veto Stopping Vertices", 100, -150.,150., 100, -150., 150., 100, -100., 800.);
	TH3D* mrdstopvertices = new TH3D("mrdstopvertices", "Distribution of MRD Stopping Vertices", 100, -150.,150., 100, -150., 150., 100, -100., 800.);
#endif // defined WCSIMDEBUG
	
	// test the hypothesis that track length in water can be estimated from total light in tank
	TH2D* tracklengthvsmuonlight = new TH2D("tracklengthvsmuonlight", "Muon Track Length vs Total Light from Muon", 100, 0., 1500., 100, 0., 50.);
	
	// Effects Of Pions In Final State - map of hits on the wall with both charge and time of the hits
	TH3D* chargemap_nopions = new TH3D("chargemap_nopions", "Charge Distribution for CC0pi events", pmtsperring+2,-1,pmtsperring+1,numpmtrings+2,-1,numpmtrings+1, 100, 0., 1400.);
	
	// Just to test the inside/outside cherenkov cone algorithm
	TH2D* inconehistowall = new TH2D("chargemap_incone_wall", "Charge Distribution Inside Cherenkov Cone (Wall)", pmtsperring+2,-1,pmtsperring+1,numpmtrings+2,-1,numpmtrings+1);
	TH2D* inconehistotop = new TH2D("chargemap_incone_top","Charge Distribution Inside Cherenkov Cone (Top Cap)",caparraysize+2,-1,caparraysize+1,caparraysize+2,-1,caparraysize+1);
	TH2D* inconehistobottom = new TH2D("chargemap_incone_bottom","Charge Distribution Inside Cherenkov Cone (Bottom Cap)",caparraysize+2,-1,caparraysize+1,caparraysize+2,-1,caparraysize+1);
	
	TH2D* outconehistowall = new TH2D("chargemap_outcone_wall", "Charge Distribution Outside Cherenkov Cone (Wall)", pmtsperring+2,-1,pmtsperring+1,numpmtrings+2,-1,numpmtrings+1);
	TH2D* outconehistotop = new TH2D("chargemap_outcone_top","Charge Distribution Outside Cherenkov Cone (Top Cap)",caparraysize+2,-1,caparraysize+1,caparraysize+2,-1,caparraysize+1);
	TH2D* outconehistobottom = new TH2D("chargemap_outcone_bottom","Charge Distribution Outside Cherenkov Cone (Bottom Cap)",caparraysize+2,-1,caparraysize+1,caparraysize+2,-1,caparraysize+1);

	TH2D* tothistowall = new TH2D("chargemap_total_wall", "Total Charge Distribution (Wall)", pmtsperring+2,-1,pmtsperring+1,numpmtrings+2,-1,numpmtrings+1);
	TH2D* tothistotop = new TH2D("chargemap_total_top","Total Charge Distribution (Top Cap)",caparraysize+2,-1,caparraysize+1,caparraysize+2,-1,caparraysize+1);
	TH2D* tothistobottom = new TH2D("chargemap_total_bottom","Total Charge Distribution (Bottom Cap)",caparraysize+2,-1,caparraysize+1,caparraysize+2,-1,caparraysize+1);
	
	// wall map histos
	std::map<std::string, TH2D*> maphistos;
	maphistos.emplace("inconehistowall",inconehistowall);
	maphistos.emplace("inconehistotop",inconehistotop);
	maphistos.emplace("inconehistobottom",inconehistobottom);
	maphistos.emplace("outconehistowall",outconehistowall);
	maphistos.emplace("outconehistotop",outconehistotop);
	maphistos.emplace("outconehistobottom",outconehistobottom);
	maphistos.emplace("tothistowall",tothistowall);
	maphistos.emplace("tothistotop",tothistotop);
	maphistos.emplace("tothistobottom",tothistobottom);
	
	
	// Checking Digit Charges And Times
	TH1D* digitsqpmthist = new TH1D("digitsqpmthist","Digit Charges for PMTs",100,0.,30.);
	TH1D* digitsqlappdhist = new TH1D("digitsqlappdhist","Digit Charges for LAPPDs",100,0.,30.);
	TH1D* digitstpmthist = new TH1D("digitstpmthist","Digit Times for PMTs",1000,800.,3000.);
	TH1D* digitstlappdhist = new TH1D("digitstlappdhist","Digit Times for LAPPDs",1000,800.,3000.);
	TH1D* digittsmearpmthist = new TH1D("digittsmearpmthist","Digit T Smearings for PMTs",100,0.,5.);
	TH1D* digittsmearlappdhist = new TH1D("digittsmearlappdhist","Digit T Smearings for LAPPDs",100,0.,0.1);
	TH2D* pmttimesmearvsqhist = new TH2D("pmttimesmearvsqhist","PMT Q vs T Smearing",1000,0.,3.,100,0.,100.);
	TH2D* lappdtimesmearvsqhist = new TH2D("lappdtimesmearvsqhist","LAPPD Q vs T Smearing",1000,0.,0.1,100,0.,100.);
	
	
	// ===================================================================================================
	// Loop over dirt entries
	// ===================================================================================================
	gROOT->cd();
	cout<<"loading first tree from "<<chainpattern<<" tchain"<<endl;
	c->LoadTree(0);
	Int_t treeNumber = -1;
	Int_t thistreesentries;
	Int_t wcsimtreeNumber = -1;
	Int_t numentswcsimchain;
	
	/*
	1. Load next g4dirt entry
	2. Check if genie primary, and volume is in tank - if not, continue
	3. If so, load genie entry.
	4. Check if interaction is QE - if not, continue
	5. If so, load wcsim detector response. 
	6. Load tracks, look for primary mu track through the MRD, record interaction details. 
	*/
	
#ifdef TIME_EVENTS
	TStopwatch* timer = new TStopwatch();
#endif // defined TIME_EVENTS
	
	cout<<"looping over tchain entries"<<endl;
	numents=100;
	Int_t wcsimTentry;
	// since WCSim only propagated tank events, there is no longer a 1:1 mapping between event numbers
	// in dirt files and WCSim files. As long as the selection criterion for dirt events is the same here
	// as in WCSim's PrimaryGeneratorAction, we can select the dirt files, and then just pull the next 
	// wcsim entry
	Int_t inputEntry=0;
#ifdef SPLITWCSIMFILES
	const char* inputoffsetchars = gSystem->Getenv("INFILE_OFFSET");
	if(inputoffsetchars==nullptr||strcmp(inputoffsetchars,"")==0){
		cerr<<"INPUTOFFSET NOT DEFINED!!"<<endl; assert(false);
	}
	inputEntry = atoi(inputoffsetchars);
	// we could set the num entries to process based on SPLITFACTOR, but this shouldn't be necessary
	// as the loop will exit normally once we reach the end of the wcsimfile.
#endif // defined USE_GRID
	for(; inputEntry<numents; inputEntry++){
#ifdef TIME_EVENTS
		timer->Start();
#endif // defined TIME_EVENTS
		//==================================================================================================
		// Load next g4dirt entry
		//==================================================================================================
#ifdef VERBOSE
		cout<<"loading entry "<<inputEntry<<endl;
#endif // VERBOSE
		Long64_t localEntry = c->LoadTree(inputEntry);
		if( localEntry<0){ cout<<"end of tchain"<<endl; break; }
		Int_t nextTreeNumber = c->GetTreeNumber();
		if(treeNumber!=nextTreeNumber){
			cout<<"new tree: "<<nextTreeNumber<<endl;
			//==================================================================================================
			// Open All Files And Load Trees
			//==================================================================================================
#ifndef PARTICLEGUNEVENTS
			// this means we've switched file - need to load the new meta tree and genie tree.
			// first pull out the new file name
			tankflux = c->GetTree();
			dirtfile = tankflux->GetCurrentFile();
			dirtfilename=dirtfile->GetName();
			thistreesentries = tankflux->GetEntries();
			cout<<"tankflux has "<<thistreesentries<<" entries in this file"<<endl;
			
#ifndef NOGENIE
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
#endif // !defined NOGENIE
			
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
#ifndef SPLITWCSIMFILES
			wcsimfilepath = TString::Format("%s/wcsim_0.%d.root",wcsimpath,filenum);
			cout<<"corresponding wcsim file is "<<wcsimfilepath<<endl;
			if(wcsimfile) wcsimfile->Close(); wcsimfile=0;
			wcsimfile = TFile::Open(wcsimfilepath);
			if(!wcsimfile){
				cerr<<"this wcsimfile doesn't exist!"<<endl; 
				inputEntry += thistreesentries;	// skip iterator forward by all the entries in this file
				continue; 
			}
			// load the next set of wcsim event info
			wcsimT = (TTree*)wcsimfile->Get("wcsimT");
			if(!wcsimT){cerr<<"wcsimT doesn't exist!"<<endl; break; }
#else // if defined SPLITWCSIMFILES
			wcsimfilepath = TString::Format("%s/wcsim_0.%d.*.root",wcsimpath,filenum);
			cout<<"corresponding wcsim files are "<<wcsimfilepath<<endl;
			if(wcsimchain) { wcsimchain->ResetBranchAddresses(); delete wcsimchain; }
			wcsimchain =  new TChain("wcsimT");
			wcsimchain->Add(wcsimfilepath);
			numentswcsimchain = wcsimchain->GetEntries();
			cout<<"loaded "<<numentswcsimchain<<" entries in wcsimchain"<<endl;
			if(numentswcsimchain<1){ 
				cerr<<"these wcsim files don't exist, or have no entries!"<<endl;
				inputEntry += thistreesentries;
				continue;
			}
			// need to load the first tree to get the geometry if we haven't already
			wcsimchain->LoadTree(0);
			wcsimtreeNumber=wcsimchain->GetTreeNumber();
			wcsimT = wcsimchain->GetTree();
			wcsimfile = wcsimT->GetCurrentFile();
			wcsimfilepath=wcsimfile->GetName();
#endif // defined SPLITWCSIMFILES
			numwcsimentries = wcsimT->GetEntries();
			cout<<"wcsimT has "<<numwcsimentries<<" entries in this file"<<endl;
			if(numwcsimentries<1){cerr<<"wcsimT has no entries!"<<endl; break; }
			wcsimTentry=-1;
#else // if defined PARTICLEGUNEVENTS
			wcsimT = c->GetTree();
			wcsimfile = wcsimT->GetCurrentFile();
			wcsimfilestring=wcsimfile->GetName();
			thistreesentries = wcsimT->GetEntries();
			numwcsimentries=thistreesentries;
			cout<<"wcsimT has "<<thistreesentries<<" entries in this file"<<endl;
			
			// use regexp to pull out the file number needed for identifying the corresponding lappd file
			std::match_results<string::const_iterator> submatches;
			// filename is of the form "wcsim_0.####.root"
			// #### is input file num. Need this to match against lappd file names
			//std::regex theexpression (".*/[^0-9_]+_0\\.([0-9]+)\\.root");  // matches "wcsim_0.####.root"
			std::regex theexpression (".*/[^0-9_]+_([0-9]+)\\.root");        // matches "wcsim_####.root"
			cout<<"matching regex for filename "<<wcsimfilestring<<endl;
			std::string wcsimfilename(wcsimfilestring.Data());
			std::regex_match (wcsimfilename, submatches, theexpression);
			std::string submatch = (std::string)submatches[0];	// match 0 is 'whole match' or smthg
			if(submatch==""){ cerr<<"unrecognised input file pattern: "<<dirtfilename<<endl; return; }
			submatch = (std::string)submatches[1];
			cout<<"extracted submatch is "<<submatch<<endl;
			int filenum = atoi(submatch.c_str());
#endif // defined PARTICLEGUNEVENTS
			
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
#if FILE_VERSION>3
				//TODO: save these
				PMTNames = geo->GetPMTNames();
				NumPMTsByType = geo->GetWCNumPmts();
				// for files with fPmtTypeName, we can use that. 
				// for files before b534c8135093494cce6c421419a00e09df2478eb, 
				// GetWCPMTNameAt(index) method of WCSimRootGeom was broken, so this info may not be available
				
#else // if FILE_VERSION<=3 (for LAPPD branch files, all PMTs are 8")
				PMTNames=std::vector<std::string>{"PMT8inch"};
				NumPMTsByType = std::vector<int>{numpmts};
#endif // FILE_VERSION<=3
			}
			
#ifndef NOLAPPDS
			// use the filenum to open the corresponding lappd file
#ifndef SPLITWCSIMFILES
#ifndef PARTICLEGUNEVENTS
			lappdfilepath = TString::Format("%s/wcsim_lappd_0.%d.root",wcsimpath,filenum);
#else // if defined PARTICLEGUNEVENTS
			lappdfilepath = TString::Format("%s/wcsim_lappd_%d.root",wcsimpath,filenum);       // for files "wcsim_lappd_####.root"
#endif // defined PARTICLEGUNEVENTS
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
#else  // if defined SPLITWCSIMFILES
			lappdfilepath = TString::Format("%s/wcsim_lappd_0.%d.*.root",wcsimpath,filenum);
			cout<<"corresponding lappd files are "<<lappdfilepath<<endl;
			if(lappdchain) { lappdchain->ResetBranchAddresses(); delete lappdchain; }
			lappdchain =  new TChain("LAPPDTree");
			// for files "wcsim_lappd_0.YYYY.###.root"
			lappdchain->Add(lappdfilepath);
			Int_t numentslappdchain = lappdchain->GetEntries();
			cout<<"loaded "<<numentslappdchain<<" entries in lappdchain"<<endl;
			if(numentslappdchain<1){ 
				cerr<<"these lappd files don't exist, or have no entries!"<<endl;
				inputEntry += thistreesentries;
				continue;
			}
			// need to load the first tree to compare number of events
			lappdchain->LoadTree(0);
			lappdtree = lappdchain->GetTree();
			lappdfile = lappdtree->GetCurrentFile();
			lappdfilepath=lappdfile->GetName();
#endif // defined SPLITWCSIMFILES
			numlappdentries = lappdtree->GetEntries();
			cout<<"lappdtree has "<<numlappdentries<<" entries in this file"<<endl;
			if(numlappdentries!=numwcsimentries){
				cerr<<"NUM LAPPD TREE ENTRIES != NUM WCSIMT ENTRIES!!!"<<endl;
				break;
			}
			if(numlappdentries<1){cerr<<"lappdtree has no entries!"<<endl; break;}
#else
			lappdfilepath="";
#endif // ndefined NOLAPPDS
			
#ifndef PARTICLEGUNEVENTS
			//==================================================================================================
			// Set Tree Branch Addresses
			//==================================================================================================
			// tankflux:
			// genie file entry number for each entry, to get the genie intx info
			c->SetBranchAddress("entry",&genieentry,&genieentrybranch);
			// number of primaries (so we can create appropriately sized array)
			c->SetBranchAddress("ntank",&ntankbranchval, &nTankBranch);
			// material of vertex - identify as 'TankWater' to pull only primaries in the tank
			c->SetBranchAddress("vtxmat",&vertexmaterial, &vertexmaterialbranch);
			// array of whether particle is a genie primary
			nuprimaryBranch=c->GetBranch("primary");
#ifndef NOGENIE
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
			
			// gtree:
			gtree->SetBranchAddress("gmcrec",&genierecordval,&genierecordBranch);
#endif // !defined NOGENIE
#endif // !defined PARTICLEGUNEVENTS
			
			// wcsimT:
			// wcsim trigger classes
			wcsimT->SetBranchAddress("wcsimrootevent",&b, &bp);
			wcsimT->SetBranchAddress("wcsimrootevent_mrd",&m, &mp);
			wcsimT->SetBranchAddress("wcsimrootevent_facc",&v, &vp);
			bp->SetAutoDelete(kTRUE);
			mp->SetAutoDelete(kTRUE);
			vp->SetAutoDelete(kTRUE);
			if(bp==0||mp==0||vp==0){ cerr<<"branches are zombies!"<<endl; break; }
			
#ifndef NOLAPPDS
			// lappdtree:
			int branchesok=0;
			branchesok =lappdtree->SetBranchAddress("lappdevt",       &lappd_evtnum);
			if(branchesok<0) cerr<<"lappdevt branch error "<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappd_numhits",  &lappd_numtileshitthisevt);
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
#endif // FILE_VERSION>3
			branchesok =lappdtree->SetBranchAddress("lappdhit_edep",          &lappd_numphots);
			if(branchesok<0) cerr<<"lappd_edep="<<branchesok<<endl;
			// lappdhit_edep is an c-style array of doubles of #true photon hits on each LAPPD in an event
			branchesok =lappdtree->SetBranchAddress("lappdhit_totalpes_perlappd2", &lappd_hitchargep);
			// lappdhit_totalpes_perlappd2 is the same as lappdhit_edep in a vector. It's retrieved from the "digits" rather than SD hits, but only after each "digit" has it's time smeared and randpe generated. This is before digitizer integration, so still photon hits, but is LAPPD event-wide (ie not the rand num pe per hit, just total number of photon hits on this lappd). For PMTs (including MRD etc) this may also include dark noise hits, but for LAPPDs this isn't implemented.
			// FIXME should store the charges
			if(branchesok<0) cerr<<"lappd_hitcharge="<<branchesok<<endl;
#endif // ndef NOLAPPDS
			
			treeNumber=nextTreeNumber;
		}
		// Done Loading Files, Trees and Setting Branch Addresses
		
		geniefilestring=geniefilepath;
		dirtfilestring=TString(dirtfilename);
		wcsimfilestring=wcsimfilepath;
		lappdfilestring=lappdfilepath;
		
		//====================================================================================================
		// Get G4dirt event info
		//====================================================================================================
		
#ifdef VERBOSE
		cout<<"processing inputEntry "<<inputEntry<<", localEntry "<<localEntry
		    <<"/"<<thistreesentries<<" in tree "<<treeNumber<<endl;
#endif // VERBOSE
		
#ifndef PARTICLEGUNEVENTS // skip this whole dirt + genie event section if using WCSim particle gun events
		
		dirteventnum=localEntry;
		nTankBranch->GetEntry(localEntry);
		vertexmaterialbranch->GetEntry(localEntry);
		if(strcmp(vertexmaterial,"TankWater")!=0){
#ifdef VERBOSE
			cout<<"neutrino vtx not in tank"<<endl;
#endif // VERBOSE
#ifdef TANKONLY
			continue;
#endif // def TANKONLY
		} else {
			numneutrinoeventsintank++;
			isintank=true;
		}
		
		if(nuprimarybranchval){delete[] nuprimarybranchval;}
		nuprimarybranchval = new Int_t[ntankbranchval];
		nuprimaryBranch->SetAddress(nuprimarybranchval);
		nuprimaryBranch->GetEntry(localEntry);
		
		Bool_t primariesinthisentry=false;
		for(int i=0;i<ntankbranchval;i++){
			if(nuprimarybranchval[i]==1){ primariesinthisentry=true; break; }
		}
		if(!primariesinthisentry){
			cout<<"dirt particles were not genie primaries"<<endl;
#ifdef TANKONLY
			continue; // we shouldn't hit this, we should've already continued above.
#endif // def TANKONLY
		}
		
		// These selection criteria are the WCSim PrimaryGeneratorAction ones
		// (i.e TANKONLY: vtx material is tankwater, dirt particles are genie primaries == same as saying intx vertex is in tank)
		// Any event that passes here will have created a WCSimT entry
		// do this now to keep file processing in sync, in case we introduce any 'continue' statements later
		wcsimTentry++;
#ifdef VERBOSE
		//cout<<"INCREMENTED WCSIMTENTRYNUM TO "<<wcsimTentry<<endl;
#endif // VERBOSE
		
		//====================================================================================================
		// Get Genie event info
		//====================================================================================================
#ifndef NOGENIE
#ifdef VERBOSE
		cout<<"getting genie info"<<endl;
#endif // VERBOSE
		if(localEntry>(numgenietentries-1)){ cerr<<"can't load localEntry "<<localEntry
								 <<" from "<<geniefilepath<<" gtree: not enough entries!"<<endl; continue; }
		genieentrybranch->GetEntry(localEntry);
		genierecordBranch->GetEntry(genieentry);
		genie::EventRecord* gevtRec = genierecordval->event;
		genie::Interaction* genieint = gevtRec->Summary();
		
		GetGenieEntryInfo(gevtRec, genieint, thegenieinfo, printneutrinoevent);  // fill thegenieinfo struct with all the genie info
		if(printneutrinoevent) PrintNeutrinoEvent(thegenieinfo);
		
		genieeventnum=genieentry;
		interactionvertex=*(thegenieinfo.IntxVtx);
		eventtypes=thegenieinfo.eventtypes;
		interactiontypestring=thegenieinfo.procinfostring;
		// thegenieinfo.procinfostring gives format "<DIS - Weak[CC]>" for which symbols might not be ideal. Strip them.
		interactiontypestring=interactiontypestring.substr(1,interactiontypestring.length() - 2);
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
		geniefslangle=thegenieinfo.fslangle;
		geniefsle=thegenieinfo.fsleptonenergy;
		if(isintank&&eventtypes.at("IsWeakCC")) numCCneutrinoeventsintank++;
		else if(isintank&&eventtypes.at("IsWeakNC")) numNCneutrinoeventsintank++;
		if(isintank&&eventtypes.at("IsWeakCC")&&eventtypes.at("IsQuasiElastic")){ numCCQEneutrinoeventsintank++; }
		
#ifdef VERBOSE
		cout<<"filling incident neutrino histograms"<<endl;
#endif // VERBOSE
		incidentneutrinoenergiesall->Fill(thegenieinfo.probeenergy);
		incidentneutrinoanglesall->Fill(thegenieinfo.probeangle);
		fslanglesall->Fill(thegenieinfo.fslangle);
		fslenergiesall->Fill(thegenieinfo.fsleptonenergy);
		eventq2all->Fill(thegenieinfo.Q2);
		neutrinovertex->Fill(thegenieinfo.Intx_x, thegenieinfo.Intx_y, thegenieinfo.Intx_z);
		if(eventtypes.at("IsWeakCC") && eventtypes.at("IsQuasiElastic"))
			neutrinovertexQE->Fill(thegenieinfo.Intx_x, thegenieinfo.Intx_y, thegenieinfo.Intx_z);
		
		isinfiducialvol=false;
		if( (TMath::Sqrt(TMath::Power(thegenieinfo.Intx_x, 2) 
			+ TMath::Power(thegenieinfo.Intx_z-tank_start-tank_radius,2)) < fidcutradius) && 
			(TMath::Abs(thegenieinfo.Intx_y-tank_yoffset) < fidcuty) && 
			((thegenieinfo.Intx_z-tank_start-tank_radius) < fidcutz) ){
			isinfiducialvol=true;
			incidentneutrinoenergiesallfidcut->Fill(thegenieinfo.probeenergy);
			fslanglesallfidcut->Fill(thegenieinfo.fslangle);
			fslenergiesallfidcut->Fill(thegenieinfo.fsleptonenergy);
			eventq2allfidcut->Fill(thegenieinfo.Q2);
			if(eventtypes.at("IsWeakCC") && eventtypes.at("IsQuasiElastic")) numCCQEneutrinoeventsinfidvol++;
		}
#endif // !defined NOGENIE
		
#else  // if defined PARTICLEGUNEVENTS
		wcsimTentry=localEntry;
#endif // defined PARTICLEGUNEVENTS - end of skipping whole dirt and genie sections for wcsim particle gun
		
		//====================================================================================================
		// Set WCSim Branch Addresses - for split wcsim files need to do this more often than dirt / genie files
		//====================================================================================================
#ifdef VERBOSE
		cout<<"getting wcsim entry "<<wcsimTentry<<endl;
#endif // VERBOSE
		
#ifdef SPLITWCSIMFILES
		Long64_t localwcsimTentry = wcsimchain->LoadTree(wcsimTentry);
#ifndef NOLAPPDS
		lappdchain->LoadTree(wcsimTentry);
#endif // ndef NOLAPPDS
		if( localwcsimTentry<0){
			cerr<<"end of wcsimchain!! should have reloaded before this"<<endl;
			inputEntry += (numentswcsimchain-thistreesentries);
			continue;
		}
		Int_t nextwcsimTreeNumber = wcsimchain->GetTreeNumber();
		if(wcsimtreeNumber!=nextwcsimTreeNumber){
			cout<<"new wcsim tree: "<<nextwcsimTreeNumber<<endl;
			// LOAD NEW TREES
			wcsimT = wcsimchain->GetTree();
			wcsimfile = wcsimT->GetCurrentFile();
			wcsimfilepath = wcsimfile->GetName();
			numwcsimentries = wcsimT->GetEntries();
			cout<<"wcsimT has "<<numwcsimentries<<" entries in this file"<<endl;
			if(numwcsimentries<1){cerr<<"wcsimT has no entries!"<<endl; break; }
			
#ifndef NOLAPPDS
			lappdtree = lappdchain->GetTree();
			lappdfile = lappdtree->GetCurrentFile();
			lappdfilepath = lappdfile->GetName();
			numlappdentries = lappdtree->GetEntries();
			cout<<"lappdtree has "<<numlappdentries<<" entries in this file"<<endl;
			if(numlappdentries<1){cerr<<"lappdtree has no entries!"<<endl; break; }
#else
			lappdfilepath = "";
#endif // ndef NOLAPPDS
			
			// SET NEW BRANCH ADDRESSES
			// wcsimT:
			wcsimT->SetBranchAddress("wcsimrootevent",&b, &bp);
			wcsimT->SetBranchAddress("wcsimrootevent_mrd",&m, &mp);
			wcsimT->SetBranchAddress("wcsimrootevent_facc",&v, &vp);
			bp->SetAutoDelete(kTRUE);
			mp->SetAutoDelete(kTRUE);
			vp->SetAutoDelete(kTRUE);
			if(bp==0||mp==0||vp==0){ cerr<<"branches are zombies!"<<endl; break; }
			
#ifndef NOLAPPDS
			// lappdtree:
			int branchesok=0;
			branchesok =lappdtree->SetBranchAddress("lappdevt",                    &lappd_evtnum);
			if(branchesok<0) cerr<<"lappdevt branch error "<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappd_numhits",               &lappd_numtileshitthisevt);
			if(branchesok<0) cerr<<"lappd_numhits="<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappdhit_stripcoorx",         &lappd_hitpeposxp);
			if(branchesok<0) cerr<<"lappd_stripcoorx="<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappdhit_stripcoory",         &lappd_hitpeposyp);
			if(branchesok<0) cerr<<"lappd_stripcoory="<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappdhit_stripcoorz",         &lappd_hitpeposzp);
			if(branchesok<0) cerr<<"lappd_stripcoorz="<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappdhit_stripcoort",         &lappd_hittruetimep);
			if(branchesok<0) cerr<<"lappd_stripcoort="<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappdhit_objnum",             &lappd_hittile);
			if(branchesok<0) cerr<<"lappd_objnum="<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappdhit_x",                  &lappd_hittilesposx);
			if(branchesok<0) cerr<<"lappd_z="<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappdhit_y",                  &lappd_hittilesposy);
			if(branchesok<0) cerr<<"lappd_y="<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappdhit_z",                  &lappd_hittilesposz);
			if(branchesok<0) cerr<<"lappd_z="<<branchesok<<endl;
#if FILE_VERSION>3
			branchesok =lappdtree->SetBranchAddress("lappdhit_globalcoorx",        &lappd_hitglobalposxp);
			if(branchesok<0) cerr<<"lappd_hitglobalx="<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappdhit_globalcoory",        &lappd_hitglobalposyp);
			if(branchesok<0) cerr<<"lappd_hitglobaly="<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappdhit_globalcoorz",        &lappd_hitglobalposzp);
			if(branchesok<0) cerr<<"lappd_hitglobalz="<<branchesok<<endl;
#endif // FILE_VERSION>3
			branchesok =lappdtree->SetBranchAddress("lappdhit_edep",               &lappd_numphots);
			if(branchesok<0) cerr<<"lappd_edep="<<branchesok<<endl;
			branchesok =lappdtree->SetBranchAddress("lappdhit_totalpes_perlappd2", &lappd_hitchargep);
			if(branchesok<0) cerr<<"lappd_hitcharge="<<branchesok<<endl;
#endif // ndef NOLAPPDS
			
			wcsimtreeNumber=nextwcsimTreeNumber;
			// when loading a new file, we need to re-synchronize dirt/genie event num = wcsim event num
			// we should always have wcsimTentry<localEntry, as wcsimTentry only skips some outside tank
			// localEntries. We now skip wcsimTentry back up to localEntry.
			wcsimTentry=localEntry;
			localwcsimTentry = wcsimchain->LoadTree(wcsimTentry); // this shouldn't result in a new tree.
		}
		wcsimeventnum=localwcsimTentry;
		wcsimT->GetEntry(localwcsimTentry);
#else  // if !defined SPLITWCSIMFILES
		wcsimeventnum=wcsimTentry;
		if(wcsimTentry>(numwcsimentries-1)){ cout<<"can't load wcsimT entry "<<wcsimTentry
				<<" from "<<wcsimfilepath<<" wcsimT - not enough entries!"<<endl; break; }
		wcsimT->GetEntry(wcsimTentry);
#endif // !defined SPLITWCSIMFILES
		
#ifndef NOLAPPDS
#ifdef VERBOSE
		cout<<"getting lappd entry "<<wcsimTentry<<endl;
#endif // VERBOSE
#ifndef SPLITWCSIMFILES
		if(wcsimTentry>(numlappdentries-1)){ cerr<<"can't load lappdtree entry "<<wcsimTentry
				<<" from "<<lappdfilepath<<" lappdtree - not enough entries!"<<endl; continue; }
		lappdtree->GetEntry(wcsimTentry);
#else // if defined SPLITWCSIMFILES
		lappdtree->GetEntry(localwcsimTentry);
#endif // defined SPLITWCSIMFILES
		int eventnumcheck=b->GetTrigger(0)->GetHeader()->GetEvtNum();
		if(lappd_evtnum!=eventnumcheck){
			cerr<<"mismatch between lappd_evtnum="<<lappd_evtnum<<", and wcsim header event num, "
			<<"eventnumcheck="<<eventnumcheck<<". For ref: wcsimTentry="<<wcsimTentry<<endl;
		}
#endif // ndef NOLAPPDS
		// End Set Branch Addresses for split wcsim files
		
		//====================================================================================================
		// Get WCSim trigger and Truth track info
		//====================================================================================================
		//
		for(int wcsimtriggernum=0; wcsimtriggernum<(b->GetNumberOfEvents()); wcsimtriggernum++){
#ifdef VERBOSE
			cout<<"Analysing WCSim event "<<wcsimTentry<<", trigger "<<wcsimtriggernum<<endl;
#endif
			atrigt = b->GetTrigger(wcsimtriggernum);
			(wcsimtreeNumber < m->GetNumberOfEvents() ) ? atrigm = m->GetTrigger(wcsimtriggernum) : atrigm = nullptr;
			(wcsimtreeNumber < v->GetNumberOfEvents() ) ? atrigv = v->GetTrigger(wcsimtriggernum) : atrigv = nullptr;
			
			// check consistency of genie vtx and wcsim vertex, as sanity check
			if( (wcsimtriggernum==0) && 
				(    ((thegenieinfo.Intx_x-atrigt->GetVtx(0))>0.1)  // float comparisons...
				  || ((thegenieinfo.Intx_y-atrigt->GetVtx(1))>0.1)
				  || ((thegenieinfo.Intx_z-atrigt->GetVtx(2))>0.1)
				)
			){
				cerr<<"WCSim vertex is not the same as genie vertex! Entries may be misaligned!"<<endl
					<<"Genie Vtx: ("<<thegenieinfo.Intx_x<<", "<<thegenieinfo.Intx_y<<", "<<thegenieinfo.Intx_z<<"), "<<endl
					<<"WCSim Vtx: ("<<atrigt->GetVtx(0)<<", "<<atrigt->GetVtx(1)<<", "<<atrigt->GetVtx(2)<<")"<<endl;
					assert(false);
			}
			
			// process WCSim truth tracks
			Int_t numtracks = atrigt->GetNtrack();
#ifdef VERBOSE
			cout<<"wcsim event had "<<numtracks<<" truth tracks"<<endl;
#endif // VERBOSE
			
			numpizerotracks=0;
			numpiplustracks=0;
			numpiminustracks=0;
			nummutracks=0;
			numgammatracks=0;
			numneutrontracks=0;
			numprotontracks=0;
			
			trackpdg.clear();
			trackstartpos.clear();
			trackstartvol.clear();
			trackstarttime.clear();
			trackstartE.clear();
			trackstartmom.clear();
			trackstoppos.clear();
			trackstopvol.clear();
			trackstoptime.clear();
			trackstopE.clear();
			trackstopmom.clear();
			trackparenttype.clear();
			numtankdigitsfromtrack.clear();
			tankchargefromtrack.clear();
			fractionaltrackchargeincone.clear();
			upstreamchargefromtrack.clear();
			downstreamchargefromtrack.clear();
			topcapchargefromtrack.clear();
			bottomcapchargefromtrack.clear();
			trackangle.clear();
			trackentersMRD.clear();
			trackstopsinMRD.clear();
			trackpenetratesMRD.clear();
			trackmrdpenetrationcm.clear();
			trackmrdpenetrationlayers.clear();
			tracklengthinMRD.clear();
			tracklengthintank.clear();
			numMRDdigitsfromtrack.clear();
			MRDchargefromtrack.clear();
			
			std::vector<int> trackindicestoscan;
			Double_t maxmuonenergy=0.;
			Double_t maxprimarymuonenergy=0.;
			primaryneutrinoindex=-1;
			primarymuonindex=-1;
			highestenergymuonindex=-1;
			
			std::map<int,int> eventparticles;
			
			// now scan through the truth tracks, find the primary muon and save the wcsim info from it
			// UPDATED; to pull just the highest energy *primary* muon track, not all muons
			// (there are on average ~1.4 muons per event!) let's just do a quick scan first and pull
			// out only the highest energy track.
			// UPDATED 2; record all information about all interesting particles. 
			for(int track=0; track<numtracks; track++){
				WCSimRootTrack* nextrack = (WCSimRootTrack*)atrigt->GetTracks()->At(track);
				/* a WCSimRootTrack has methods: 
				Int_t     GetIpnu()             pdg
				Int_t     GetFlag()             -1: neutrino primary, -2: neutrino target, 0: other
				Float_t   GetM()                mass
				Float_t   GetP()                momentum magnitude
				Float_t   GetE()                energy (inc rest mass^2)
				Float_t   GetEndE()             energy on stopping of particle tracking
				Float_t   GetEndP()             momentum on stopping of particle tracking
				Int_t     GetStartvol()         starting volume: 10 is tank, 20 is facc, 30 is mrd
				Int_t     GetStopvol()          stopping volume: but these may not be set.
				Float_t   GetDir(Int_t i=0)     momentum unit vector
				Float_t   GetPdir(Int_t i=0)    momentum vector
				Float_t   GetPdirEnd(Int_t i=0) direction vector on stop tracking
				Float_t   GetStop(Int_t i=0)    stopping vertex x,y,z for i=0-2, in cm
				Float_t   GetStart(Int_t i=0)   starting vertex x,y,z for i=0-2, in cm
				Int_t     GetParenttype()       parent pdg, 0 for primary.
				Float_t   GetTime()             trj->GetGlobalTime(); starting time of particle
				Float_t   GetStopTime()
				Int_t     GetId()               wcsim trackid
				GetFlag=-1; neutrino or... first(?) primary?. This should be the neutrino. Only stop vertex is stored. 
				GetFlag=-2; primary target. This may or may not be set? ???? what info stored?
				GetFlag= 0; any other track. ALL PRIMARY TRACKS ALSO GET STORED HERE. TO avoid duplication, only count these tracks!!!
				*/
				if(nextrack->GetFlag()==-1) primaryneutrinoindex=track; // always 0?
				if(nextrack->GetFlag()!=0) continue; // flags -1 and -2 are neutrino and target, or other special. don't double count.
				if((abs(nextrack->GetIpnu())==11)&&(nextrack->GetStopvol()!=10)) continue; // ignore electrons not in the tank
				// Too many electrons!! (>200!) But we do want to keep Michel electrons!
				Int_t primarypdg = nextrack->GetIpnu();
				if(eventparticles.count(primarypdg)==0) eventparticles.emplace(primarypdg,1); else eventparticles.at(primarypdg)++;
				switch (primarypdg){
					case 111: numpizerotracks++; break;
					case 211: numpiplustracks++; break;
					case -211: numpiminustracks++; break;
					case 13: nummutracks++; break;
					case 22: numgammatracks++; break;
					case 2112: numneutrontracks++; break;
					case 2212: numprotontracks++; break;
				}
				
				// check if this track is of interest to us:
				if(particlestosave.count(PdgToString(nextrack->GetIpnu()))==0) continue; // SET SELECTION CRITERIA HERE
				
				// old versions had broken start / stop volumes
				Int_t primarystartvol;
				Int_t primarystopvol;
#if FILE_VERSION<2   // fixed on 2nd Mar 17 in 595163fd20592621161415eee1824d4b3c9744f2
				if(nextrack->GetStart(2)<tank_start){
					primarystartvol = 20;                               // start depth is facc or hall
				} else if(nextrack->GetStart(2)>(tank_start+(2.*tank_radius))){
					primarystartvol = 30;                               // start depth is mrd or hall
//				} else if((TMath::Sqrt(TMath::Power(nextrack->GetStart(0),2)+
//						     TMath::Power(nextrack->GetStart(2)-tank_start-tank_radius,2))<tank_radius) &&
//						  (TMath::Abs(nextrack->GetStart(1))<tank_halfheight)){
//				primarystartvol = 10;                               // start is tank // why was this disabled?
				} else {
//					primarystartvol = 30;                               // start depth is about tank, but outside tank region - hall
					primarystartvol = 10;                               // start depth is tank
				}
				if(primarystartvol!=10){
					cerr<<"START VOLUME IS 10 BUT DEPTH ISN'T CORRECT FOR TANK!?"<<endl;
					cerr<<"start depth is "<<nextrack->GetStart(2)<<endl;
					cerr<<"tank starts at "<<(tank_start)<<", ends at "<<(tank_start+(2.*tank_radius))<<endl;
					//assert(false);
					continue;
				}
				
				if(nextrack->GetStop(2)<tank_start){
					primarystopvol = 20;						// start depth is facc or hall
				} else if(nextrack->GetStop(2)>(tank_start+(2.*tank_radius))){
					primarystopvol = 30;						// start depth is mrd or hall
				} else {
					primarystopvol = 10;						// start depth is tank
				}
#else
				primarystartvol = nextrack->GetStartvol();
				primarystopvol = nextrack->GetStopvol();
#endif
				// old checks
				//if(TMath::Abs(primarypdg)!=13) continue;                  // not a muon
				//if(nextrack->GetParenttype()!=0) continue;                // not a primary
				//if(nextrack->GetStartvol()!=10) continue;                 // track doesn't start in tank
				
				// mark for saving further data later
				trackindicestoscan.push_back(track);
				// record track info we can access here
				trackparenttype.push_back(nextrack->GetParenttype());
				trackpdg.push_back(nextrack->GetIpnu());
				TVector3 startpos(nextrack->GetStart(0),nextrack->GetStart(1),nextrack->GetStart(2));
				trackstartpos.push_back(startpos);
				trackstartvol.push_back(primarystartvol);
				trackstarttime.push_back(nextrack->GetTime());
				trackstartE.push_back(nextrack->GetE());
				TVector3 startmom(nextrack->GetPdir(0),nextrack->GetPdir(1),nextrack->GetPdir(2));
				trackstartmom.push_back(startmom);
				
				// sanity check that primary particles vertices match dirt vertices in case of tank
				if((trackparenttype.back()==0)&&(strcmp(vertexmaterial,"TankWater")==0)&&(primarystartvol!=10)){
					cerr<<"g4dirt vertexmaterial is TankWater but WCSim primary particle has startvol not tank!?!"<<endl;
				}
				
				TVector3 stoppos(nextrack->GetStop(0),nextrack->GetStop(1),nextrack->GetStop(2));
				trackstoppos.push_back(stoppos);
				trackstopvol.push_back(primarystopvol);
#if FILE_VERSION>2
				trackstoptime.push_back(nextrack->GetStopTime());
#else
				trackstoptime.push_back(-1.);
#endif // FILE_VERSION>2
				trackstopE.push_back(nextrack->GetEndE());
				TVector3 stopmom(nextrack->GetPdirEnd(0),nextrack->GetPdirEnd(1),nextrack->GetPdirEnd(2));
				trackstopmom.push_back(stopmom);
				
				if((TMath::Abs(primarypdg)==13)&&(primarystartvol==10)) nummuontracksintank++; // debug counter
				if((TMath::Abs(primarypdg)==13)&&(nextrack->GetE()>maxmuonenergy)){ highestenergymuonindex=trackindicestoscan.size()-1;
					cout<<" setting highestenergymuonindex="<<track<<endl;}
				if((TMath::Abs(primarypdg)==13)&&(nextrack->GetParenttype()==0)&&(nextrack->GetE()>maxprimarymuonenergy)) primarymuonindex=trackindicestoscan.size()-1;
#ifdef VERBOSE
				//cout<<"   noting track "<<track<<"/"<<numtracks<<" of type "<<PdgToString(primarypdg)<<endl;
#endif // VERBOSE
			} // end intial scan over wcsim truth tracks
#ifdef VERBOSE
			cout<<"entry had "<<endl;
			for(auto apair : eventparticles) cout<<"   "<<apair.second<<" "<<PdgToString(apair.first)<<"s"<<endl;
#endif
			
			int filemuonindex=-1; // which one to save to the flat file, if any
			if((primarymuonindex!=-1) || (highestenergymuonindex!=-1)){
				if(primarymuonindex!=-1) filemuonindex=primarymuonindex;
				else filemuonindex=highestenergymuonindex;
			} else {
#ifdef VERBOSE
				cout<<"no primary muon track in this event!"<<endl;
#endif // VERBOSE
			}
			if(trackindicestoscan.size()==0 && strcmp(vertexmaterial,"TankWater")==0){
				cerr<<"found no tracks to note, even though vertex was in tank!?!"<<endl;
			}
			
			//====================================================================================================
			// Scan over relevant WCSim truth tracks for additional calculated info
			//====================================================================================================
			
			// TODO depreciated vectors for filling some histograms
			std::vector<Float_t> neutrinoenergiesvector;
			std::vector<Float_t> primaryenergiesvector;
			std::vector<Double_t> scatteringanglesvector;
			std::vector<Int_t> acceptedtrackids;
			std::vector<Double_t> q2vector;
			std::vector<Double_t> muonenergydepositions;
			
			cout<<"recording "<<trackindicestoscan.size()<<" tracks"<<endl;
			
			for(int track=0; track<trackindicestoscan.size(); track++){
				WCSimRootTrack* nextrack = (WCSimRootTrack*)atrigt->GetTracks()->At(trackindicestoscan.at(track));
				int atrackstartvol = trackstartvol.at(track);
				int atrackstopvol = trackstopvol.at(track);
#ifdef VERBOSE
				cout<<"track "<<track<<" is a ";
				if(trackparenttype.at(track)==0) cout<<"primary ";
				cout<<PdgToString(trackpdg.at(track));
				cout<<" starting in the ";
				int stv = trackstartvol.at(track);
				if(stv==10) cout<<"Tank"<<endl;
				else if(stv==20) cout<<"FACC"<<endl;
				else if(stv==30) cout<<"MRD"<<endl;
				else cout<<"volume "<<stv<<endl;
#endif // VERBOSE
				
				// world extent in WCSim is +-600cm in all directions!
				TVector3 primarystartvertex = trackstartpos.at(track);
				TVector3 primarystopvertex = trackstoppos.at(track);
				TVector3 differencevector = (primarystopvertex-primarystartvertex);
				fsltruetracklength->Fill(differencevector.Mag());
				
				double atrackangle = differencevector.Angle(TVector3(0,0,1));
				trackangle.push_back(atrackangle);
				
				// continue if stopping volume is in either the tank or veto, or track is backward going.
				//if( atrackstopvol==10 || atrackstopvol==20 || primarystopvertex.Z() < primarystartvertex.Z() )
				//	continue;  
				
#if defined NOGENIE || defined PARTICLEGUNEVENTS
				// sanity check for correct synchronization between dirt and wcsim files.
				if(  (nextrack->GetParenttype()==0)                     &&  // only expect vertices to align for primaries
					!(abs(primarystartvertex.X()-thegenieinfo.Intx_x)<1 &&
					  abs(primarystartvertex.Y()-thegenieinfo.Intx_y)<1 &&
					  abs(primarystartvertex.Z()-thegenieinfo.Intx_z)<1) ){
					cerr<<"GENIE VERTEX IS IN FIDUCIAL VOLUME BUT PRIMARY MUON VERTEX ISN'T?!"<<endl
						<<"Genie vertex: ("<<thegenieinfo.Intx_x<<", "<<thegenieinfo.Intx_y<<", "<<thegenieinfo.Intx_z<<")"<<endl
						<<"Muon vertex: ("<<primarystartvertex.X()<<", "<<primarystartvertex.Y()<<", "<<primarystartvertex.Z()<<")"<<endl
						<<"Trigger vertex: ("<<atrigt->GetVtx(0)<<", "<<atrigt->GetVtx(1)<<", "<<atrigt->GetVtx(2)<<")"<<endl
						<<"dirt file is "<<dirtfilestring<<", genie file is "<<geniefilestring<<", wcsim file is "<<wcsimfilepath<<endl
						<<"dirt event num is "<<dirteventnum<<", wcsim event num is "<<wcsimTentry<<endl;
#ifdef SPLITWCSIMFILES
					cerr<<", wcsim chain local entry is "<<localwcsimTentry<<endl;
#endif // def SPLITWCSIMFILES
						assert(false);
				}
				
				// use wcsim primary vertex to count fiducial events when not using genie
				isinfiducialvol=false;
				if( (TMath::Sqrt(TMath::Power(primarystartvertex.X(), 2) 
					+ TMath::Power(primarystartvertex.Z()-tank_start-tank_radius,2)) < fidcutradius) && 
					(TMath::Abs(primarystartvertex.Y()-tank_yoffset) < fidcuty) && 
					((primarystartvertex.Z()-tank_start-tank_radius) < fidcutz) ){
					isinfiducialvol=true;
					numCCQEneutrinoeventsinfidvol++; // can't strictly say it's CCQE without genie.
				}
#endif // defined NOGENIE || defined PARTICLEGUNEVENTS
				
#ifdef WCSIMDEBUG
				// debug histograms
				switch (atrackstartvol){
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
				switch (atrackstopvol){
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
#endif // defined WCSIMDEBUG
				
				//====================================================================================================
				// calculate MRD penetration 
				//====================================================================================================
				
				// The mrd is as wide as the tank. We can ensure a track enters the mrd by projecting
				// the track forward from the start vertex at the angle between start and stop vertices,
				// and requiring:
				// 1) at z=MRD_start, x is between +/-(MRD_width/2);
				// 2) z endpoint is > MRD_start
				// For at least 2 layers of penetration, as above but with requirement on x @ z=MRD_layer2
				// For range-out, as above with requirement on x @ z=MRD_start+MRD_depth && 
				//   z endpoint is > MRD_start+MRD_depth
				
				Float_t oppx = primarystopvertex.X() - primarystartvertex.X();
				Float_t adj = primarystopvertex.Z() - primarystartvertex.Z();
				Double_t avgtrackgradx = (adj!=0) ? (oppx/adj) : 1000000;
				Float_t avgtrackanglex = TMath::ATan(avgtrackgradx);
				Float_t oppy = primarystopvertex.Y() - primarystartvertex.Y();
				Double_t avgtrackgrady = (adj!=0) ? (oppy/adj) : 1000000;
				Float_t avgtrackangley = TMath::ATan(avgtrackgrady);
				
				Float_t xatmrd = primarystartvertex.X() + (MRD_layer2-primarystartvertex.Z())*avgtrackgradx;
				Float_t yatmrd = primarystartvertex.Y() + (MRD_layer2-primarystartvertex.Z())*avgtrackgrady;
				
#ifdef WCSIMDEBUG
				cout<<"primary start vertex: ("<<primarystartvertex.X()<<", "<<primarystartvertex.Y()
					<<","<<primarystartvertex.Z()<<")"<<endl;
				cout<<"primary stop vertex: ("<<primarystopvertex.X()<<", "<<primarystopvertex.Y()
					<<","<<primarystopvertex.Z()<<")"<<endl;
				cout<<"oppx = "<<oppx<<endl;
				cout<<"adj = "<<adj<<endl;
				cout<<"angle = "<<avgtrackanglex<<endl;
				cout<<"tan(angle) = "<<avgtrackgradx<<endl;
				cout<<"projected x at z="<<MRD_layer2<<" is "<<xatmrd<<endl;
				cout<<"xatmrd="<<xatmrd<<", MRD_width="<<MRD_width<<endl;
				cout<<"yatmrd="<<yatmrd<<", MRD_height="<<MRD_height<<endl;
#endif // defined WCSIMDEBUG
				
				// variables to be saved to EventDistributions
				TVector3 MRDentrypoint(0,0,0), MRDexitpoint(0,0,0), MuTrackInMRD(0,0,0);
				double atracklengthinmrd, mrdpenetrationcm;
				int mrdpenetrationlayers;
				bool atrackentersmrd, atrackstopsinmrd, atrackpenetratesmrd;
				// check for intercept and record entry point
				bool checkboxlinerror=false;
				
				//new version based on external function calls
				///////////////////////////////////////////////
				// bool CheckLineBox( TVector3 L1, TVector3 L2, TVector3 B1, TVector3 B2, 
				//					  TVector3 &Hit, TVector3 &Hit2, bool &error)
				// returns true if line (L1, L2) intersects with the box (B1, B2)
				// returns intersection with smaller Z in Hit
				// if 2 interceptions are found, returns interception with larger Z,
				// if 1 interception is found, returns L2 (stopping point).
				// error returns true if >2 intercepts are found, or other error.
				atrackentersmrd  =  CheckLineBox( primarystartvertex, primarystopvertex, 
												TVector3(-MRD_width,-MRD_height,MRD_start), 
												TVector3(MRD_width,MRD_height,MRD_end),
												MRDentrypoint, MRDexitpoint, checkboxlinerror );
				// sanity check: XXX DISABLE TO ALLOW TRACKS STARTING IN THE MRD
				//assert(MRDentrypoint!=primarystartvertex&&"track starts in MRD!?");
				// check if MRD stops in the MRD
				atrackstopsinmrd = ( abs(primarystopvertex.X())<MRD_width&&
								   abs(primarystopvertex.Y())<MRD_height&&
								   primarystopvertex.Z()>MRD_start&&
								   primarystopvertex.Z()<(MRD_start+MRD_depth) );
				if(atrackentersmrd){
					if(abs(nextrack->GetIpnu())==13) nummuontracksinmrd++;
					atrackpenetratesmrd = ((MRDentrypoint.Z()==MRD_start)&&(MRDexitpoint.Z()==MRD_end));
					MuTrackInMRD = (MRDexitpoint-MRDentrypoint);
					atracklengthinmrd = MuTrackInMRD.Mag();
					mrdpenetrationcm = MuTrackInMRD.Z();
					mrdpenetrationlayers=0;
					for(auto layerzval : mrdscintlayers){
						if(MRDexitpoint.Z()<layerzval) break;
						mrdpenetrationlayers++;
					}
				} else {
					atrackpenetratesmrd=false;
					mrdpenetrationcm=0.;
					mrdpenetrationlayers=0;
					atracklengthinmrd=0.;
				}
				
				double muxtracklength=MuTrackInMRD.X();
				double muytracklength=MuTrackInMRD.Y();
				double muztracklength=MuTrackInMRD.Z();
#ifdef MUTRACKLENGTHDEBUG
				if(track==filemuonindex){ // set these flat file variables when processing the corresponding track
					filemuxtracklength=muxtracklength;
					filemuytracklength=muytracklength;
					filemuztracklength=muztracklength;
				}
				cout<<"particle travels "<<atracklengthinmrd<<"cm before exiting MRD bounds"<<endl
					<<"particle travels "<<muxtracklength<<"cm before leaving X bounds"<<endl
					<<"particle travles "<<muytracklength<<"cm before leaving Y bounds"<<endl;
#endif // defined MUTRACKLENGTHDEBUG
				double maxtracklengthinMRD = TMath::Sqrt(
						TMath::Power(MRD_width*2,2) +
						TMath::Power(MRD_height*2,2) +
						TMath::Power(MRD_depth,2) );
				if((atracklengthinmrd>maxtracklengthinMRD)||
					(muxtracklength>((MRD_width*2)*1.01)) || (muytracklength>((MRD_height*2)*1.01))||
					(mrdpenetrationcm>(MRD_depth*1.01))   || checkboxlinerror ){ // stupid float inaccuracies...
					cerr<<"MRD track length is bad!"<<endl
						<<"Max track length is "<<maxtracklengthinMRD<<"cm"
						<<", this track is "<<atracklengthinmrd<<"cm, "
						<<"distances are ("<<muxtracklength<<", "
						<<muytracklength<<", "<<mrdpenetrationcm<<")"<<endl
						<<"compare to maximum extents ("<<(2*MRD_width)
						<<", "<<(2*MRD_height)<<", "<<MRD_depth<<")"<<endl
						<<"Track goes = ("<<MRDentrypoint.X()<<", "<<MRDentrypoint.Y()<<", "
						<<MRDentrypoint.Z()<<") -> ("<<MRDexitpoint.X()<<", "<<MRDexitpoint.Y()<<", "
						<<MRDexitpoint.Z()<<") and ";
						if(atrackstopsinmrd) cerr<<"stops in the MRD"<<endl;
						else if(atrackpenetratesmrd) cerr<<"penetrates the MRD"<<endl;
						else cerr<<"exits the side of the MRD"<<endl;
					cerr<<"MRD width is "<<MRD_width<<", MRD height is "<<MRD_height
						<<", MRD start is "<<MRD_start<<", MRD end is "<<(MRD_start+MRD_depth)<<endl
						<<"total path is from ("<<primarystartvertex.X()<<", "<<primarystartvertex.Y()
						<<", "<<primarystartvertex.Z()<<") -> ("<<primarystopvertex.X()<<", "
						<<primarystopvertex.Y()<<", "<<primarystopvertex.Z()<<")"<<endl
						<<"avgtrackangley="<<avgtrackangley<<", avgtrackanglex="<<avgtrackanglex<<endl
						<<"CheckLineBox error is "<<checkboxlinerror<<endl;
					
					assert(false);
				}
				
				trackentersMRD.push_back(atrackentersmrd);
				trackstopsinMRD.push_back(atrackstopsinmrd);
				trackpenetratesMRD.push_back(atrackpenetratesmrd);
				trackmrdpenetrationcm.push_back(mrdpenetrationcm);
				trackmrdpenetrationlayers.push_back(mrdpenetrationlayers);
				tracklengthinMRD.push_back(atracklengthinmrd);
				
				//====================================================================================================
				// calculate the track length in water
				//====================================================================================================
				// to calculate track length _in water_ find distance from start vertex to the point
				// where it intercepts the tank. if this length > total track length, use total track length
				// otherwise use this length. 
				
				// first check if the start and endpoints are in the tank - if so, there is no tank exit point
				// and track length in tank is total length. 
				double atracklengthintank;
				if((atrackstartvol==10)&&(atrackstopvol==10)){ // stop volume is in the tank
					atracklengthintank=differencevector.Mag();
				} else {
					// either start or endpoint may be outside the track. 
					// We need to find relevant intercepts for the relevant points, in the process checking
					// there is indeed tank interception
					TVector3 tankentryvtx, tankexitvtx;
					bool interceptstank = CheckTankIntercepts(primarystartvertex, primarystopvertex, avgtrackgradx, avgtrackgrady, atrackstartvol, atrackstopvol, tankexitvtx, tankentryvtx);
					if(atrackstartvol==10) tankentryvtx=primarystartvertex;
					if(atrackstopvol==10) tankexitvtx=primarystopvertex;
					
					// we're now able to determine muon track length in the tank:
					atracklengthintank = TMath::Sqrt(
						TMath::Power((tankexitvtx.X()-tankentryvtx.X()),2)+
						TMath::Power((tankexitvtx.Y()-tankentryvtx.Y()),2)+
						TMath::Power((tankexitvtx.Z()-tankentryvtx.Z()),2) );
					
					Double_t maxtanktracklength = 
					TMath::Sqrt(TMath::Power(tank_radius*2.,2.)+TMath::Power(tank_halfheight*2.,2.));
#ifdef MUTRACKDEBUG
					cout<<"max tank track length is "<<maxtanktracklength<<endl;
					cout<<"muon tank track length: ("<<(tankexitvtx.X()-tankentryvtx.X())
						<<", "<<(tankexitvtx.Y()-tankentryvtx.Y())<<", "
						<<(tankexitvtx.Z()-tankentryvtx.Z())<<") = "<<atracklengthintank<<"cm total"<<endl;
					cout<<"muon tank exit point: ("<<tankexitvtx.X()<<", "<<tankexitvtx.Y()<<", "<<tankexitvtx.Z()<<") ";
					cout<<"track start point : ("<<tankentryvtx.X()<<", "<<tankentryvtx.Y()
						<<", "<<tankentryvtx.Z()<<")"<<endl
						<<"track stop point : ("<<primarystopvertex.X()<<", "<<primarystopvertex.Y()
						<<", "<<primarystopvertex.Z()<<")"<<endl;
#endif // defined MUTRACKDEBUG
					if(atracklengthintank > maxtanktracklength){
						cout<<"Track length is impossibly long!"<<endl;
						assert(false);
					}
					if(atracklengthintank > differencevector.Mag()){
						cout<<"Track length in tank is greater than total track length"<<endl;
						assert(false);
					}
					if(TMath::IsNaN(atracklengthintank)){
						cout<<"NaN RESULT FROM MU TRACK LENGTH IN TANK?!"<<endl;
						assert(false);
					}
				
				}
				tracklengthintank.push_back(atracklengthintank);
				
				fsltruetracklengthintank->Fill(atracklengthintank);
				if(atracklengthintank>50){ nummuontracksintankpassedcut++; }
				
				//====================================================================================================
				// track level tank digit analysis
				//====================================================================================================
#ifdef VERBOSE
				cout<<"Analysing tank digits for track level info"<<endl;
#endif // VERBOSE
				int numdigitsthisevent=atrigt->GetCherenkovDigiHits()->GetEntries();
				
				// clear counters for this track
				int numtankdigitsfromthistrack=0;
				double tankchargefromthistrack=0;
				double afractionaltrackchargeincone=0;
				double aupstreamchargefromtrack=0;
				double adownstreamchargefromtrack=0;
				double atopcapchargefromtrack=0;
				double abottomcapchargefromtrack=0;
				
				std::vector<int> tanktubeshitbyatrack;
				double chargeinsidecherenkovcone=0;  // not stored
				double chargeoutsidecherenkovcone=0; // not stored
				// we store only the fractional charge in cone
				
				//TH3D hitpositions = TH3D("hitpositions","tiel",100,-tank_radius,tank_radius, 100,tank_start,tank_start+2*tank_radius,100,-tank_halfheight+tank_yoffset,tank_halfheight+tank_yoffset);
				//ClearMapHistos(maphistos);
				
				// Loop over digits
				// ================
				for(Int_t i=0; i<numdigitsthisevent; i++){
					WCSimRootCherenkovDigiHit* thedigihit = 
						(WCSimRootCherenkovDigiHit*)atrigt->GetCherenkovDigiHits()->At(i);
					if(thedigihit==nullptr) cerr<<"DIGIT IS NULL!"<<endl;
					int digitstubeid = thedigihit->GetTubeId()-1;
					
					// add the digit charge to the counters of upstream/downstream/cap charge
					double digitsq = thedigihit->GetQ();
					WCSimRootPMT pmt = geo->GetPMT(digitstubeid);
					double digitsx = pmt.GetPosition(0);
					double digitsy = pmt.GetPosition(1);
					double digitsz = pmt.GetPosition(2);
					int thepmtsloc = pmt.GetCylLoc();
					switch (thepmtsloc){
						case 0: atopcapchargefromtrack+=digitsq; break;
						case 2: abottomcapchargefromtrack+=digitsq; break;
						case 1: ((digitsz-tank_start-tank_radius)<0) ? aupstreamchargefromtrack+=digitsq : adownstreamchargefromtrack+=digitsq; break;
					}
					
					// measure the number of digits inside and outside this track's cherenkov cone
					///////////////////////////////////////////////////////////////////////////////////
					TVector3 digitvector(digitsx-primarystartvertex.X(), digitsy-primarystartvertex.Y(), 
						digitsz-primarystartvertex.Z());
					Double_t dotproduct = (digitvector.Unit()).Dot(differencevector.Unit());
					Double_t digitmuonangle = TMath::ACos(dotproduct);  // returns [0, Pi]
					if(digitmuonangle<((42.0*TMath::Pi())/180.0)){
						chargeinsidecherenkovcone+=digitsq;
						//FillTankMapHist(geo, digitstubeid, 1, maphistos, digitsq);
						//hitpositions.Fill(digitsx, digitsz, digitsy);
					} else {
						chargeoutsidecherenkovcone+=digitsq;
						//FillTankMapHist(geo, digitstubeid, 0, maphistos, digitsq);
					}
					
					// Calculate total charge from this track
					/////////////////////////////////////////
					// scan through the parents ID's for the photons contributing to this digit
					// and see if any of them are this track
					//cout<<endl<<"calculating track energy deposition, track id "<<nextrack->GetId()<<endl;
					
					// To match digits to their parent particles we need the corresponding CherenkovHitTimes
					//---------------------------------------
					int ncherenkovhits=b->GetTrigger(0)->GetCherenkovHits()->GetEntries(); //atrigt->GetNcherenkovhits();
					int nhittimes = b->GetTrigger(0)->GetCherenkovHitTimes()->GetEntries();
#if FILE_VERSION<2
					// The CherenkovHitTimes is a flattened array (over PMTs) of arrays (over photons)
					// The  PhotonIds available from a digit are the indices within the subarray for that PMT
					// we therefore need we need to correct these indices with the start of the pmt's subarray
					// This may be found by scanning the CherenkovHits array (over PMTs), in which the offset we need is stored as the 'GetTotalPe(0)' member
					int timeArrayOffset=-1;
					for(int ihit = 0; ihit < ncherenkovhits; ihit++) {
						WCSimRootCherenkovHit* hitobject = 
							(WCSimRootCherenkovHit*)b->GetTrigger(0)->GetCherenkovHits()->At(ihit);
						int tubeNumber = hitobject->GetTubeID()-1;
						if(tubeNumber==digitstubeid) {
							timeArrayOffset = hitobject->GetTotalPe(0);
							break;
						}
					}
#endif // FILE_VERSION<2
					
					std::vector<int> truephotonindices = thedigihit->GetPhotonIds();
					//cout<<"   digit "<<i<<" has "<<truephotonindices.size()<<" true photons"<<endl;
					int numphotonsfromthistrack=0;
					for(int truephoton=0; truephoton<truephotonindices.size(); truephoton++){
						int thephotonsid = truephotonindices.at(truephoton);
						//cout<<"the photons id = "<<thephotonsid<<endl;
#if FILE_VERSION<2
						thephotonsid+=timeArrayOffset;
#endif // FILE_VERSION<2
						//cout<<"getting cherenkov hit "<<thephotonsid<<"/"<<nhittimes<<endl;
						WCSimRootCherenkovHitTime *thehittimeobject = 
							(WCSimRootCherenkovHitTime*)(b->GetTrigger(0)->GetCherenkovHitTimes()->At(thephotonsid));
						//cout<<"thehittimeobject="<<thehittimeobject<<endl;
						if(thehittimeobject==nullptr) cerr<<"HITTIME IS NULL"<<endl;
						Int_t thephotonsparenttrackid = (thehittimeobject) ? thehittimeobject->GetParentID() : -1;
						//cout<<"      HitTimeIndex "<<thephotonsid<<", HitTimeObject "<<thehittimeobject
						//	<<", HitTimeParent "<<thephotonsparenttrackid<<endl;
						if(thephotonsparenttrackid==nextrack->GetId()) {
							numphotonsfromthistrack++;
							//cout<<"################ FOUND A DIGIT ############"<<endl;
							//cout<<"the digits Q is "<<thedigihit->GetQ()<<endl;
							//return;
						}
					}
					// to estimate the charge from the muon we scale each digit's total charge
					// by the fraction of hits in the digit which come from the track
					if(numphotonsfromthistrack!=0){   //this digit had some contribution from the track!
						numtankdigitsfromthistrack++;
						//cout<<"muon digit "<<i<<" had charge "<< digitsq << " and " << truephotonindices.size()
						//	<<" true photons, of which "<<numphotonsfromthistrack <<" were from this track"<<endl
						tankchargefromthistrack += 
							digitsq * ((double)numphotonsfromthistrack/(double)truephotonindices.size());
						if(std::find(tanktubeshitbyatrack.begin(), tanktubeshitbyatrack.end(),
							 digitstubeid)==tanktubeshitbyatrack.end())
							tanktubeshitbyatrack.push_back(digitstubeid);
					}
					numtankdigitsfromtrack.push_back(numtankdigitsfromthistrack);
					tankchargefromtrack.push_back(tankchargefromthistrack);
				
				}  // end loop over digits
				
				if((chargeinsidecherenkovcone+chargeoutsidecherenkovcone)!=0){
					afractionaltrackchargeincone=
					chargeinsidecherenkovcone/(chargeinsidecherenkovcone+chargeoutsidecherenkovcone);
				} else {
					afractionaltrackchargeincone=0;
				}
				fractionaltrackchargeincone.push_back(afractionaltrackchargeincone);
				upstreamchargefromtrack.push_back(aupstreamchargefromtrack);
				downstreamchargefromtrack.push_back(adownstreamchargefromtrack);
				topcapchargefromtrack.push_back(atopcapchargefromtrack);
				bottomcapchargefromtrack.push_back(abottomcapchargefromtrack);
				tanktubeshitbytrack.push_back(tanktubeshitbyatrack);
				
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
				
				muedepositionswcsim->Fill(tankchargefromthistrack);
				
				// end of track level tank digit analysis
				////////////////////////////////////////////////
				
/*
#ifndef NOLAPPDS
				//====================================================================================================
				// track level lappd digit analysis - nothing stored here
				//====================================================================================================
				// FIXME: the vectors store all digits in the event regardless of which trigger they are in. 
				// we should scan through and only retrieve / fill jingbo file vectors with events within
				// the trigger window.
				cout<<"Analysing LAPPD digits"<<endl;
				int runningcount=0;
				
				for(int lappdi=0; lappdi<lappd_numtileshitthisevt; lappdi++){
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
#if FILE_VERSION>3
								<<"true global position is: ("<<lappd_hitglobalposx.at(runningcount)
								<<", "<<lappd_hitglobalposy.at(runningcount)
								<<", "<<lappd_hitglobalposz.at(runningcount)<<")"<<endl;
#endif // FILE_VERSION>3
						}
#endif // defined LAPPD_DEBUG
						double digitst  = lappd_hittruetime.at(runningcount);
#if FILE_VERSION>10 // dont know what file version this will be fixed in (check, has been fixed?)
						// but currently lappd digits store absolute digit time, NOT the time within a trigger....
						WCSimRootEventHeader* trigheader=atrigt->GetHeader();
						double triggertime=trigheader->GetDate();
						double absolutedigitst=digitst+triggeroffset-triggertime;
						// ALSO all lappd digits are stored regardless of being in or outside trigger window.
						// so this time may be way out, this digit may even be part of a different trigger.
						// some sort of loop over trigger times to sort it here....
						// (we could check if abs time > trigger time for next trig...)
#else // if FILE_VERSION<=10
						double absolutedigitst=digitst;
#endif // FILE_VERSION<=10
						ROOT::Math::XYZTVector adigitvector =  // convert mm to cm 
							ROOT::Math::XYZTVector(digitsx/10.,digitsy/10.,digitsz/10.,absolutedigitst);
						
#ifdef LAPPD_DEBUG
						intileposx.push_back(peposx);
						intileposy.push_back(peposy);
#if FILE_VERSION>3
						poserrx.push_back(digitsx-lappd_hitglobalposx.at(runningcount));
						poserry.push_back(digitsy-lappd_hitglobalposy.at(runningcount));
						poserrz.push_back(digitsz-lappd_hitglobalposz.at(runningcount));
#endif // FILE_VERSION>3
						tileorient.push_back(thepmtsloc);
						octagonside.push_back(theoctagonside);
#endif // defined LAPPD_DEBUG
#if FILE_VERSION>10 // don't know when this will be fixed
						double digitsq = lappd_hitcharge.at(runningcount);
#else // if FILE_VERSION<=10
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
//						int iflag;
//						SKIDigitizerThreshold(digitsq,iflag);
//						if(iflag!=0) continue; // this digit gets rejected by digitizer ... if it were a PMT
//						// TODO need to move all vector push_back's after this if we want to 'continue' 
//						// in order to keep all vectors the same size
//						digitsq*0.985;         // efficiency in WCSimWCDigitizer
#endif // FILE_VERSION<=10
						filedigitPMTIDs.push_back(LAPPDID+numpmts);
						filedigitvertices.push_back(adigitvector);
						filedigitQs.push_back(digitsq);
						filedigitsensortypes.push_back("lappd_v0"); // bad timing resolution used for pulses
						double timingConstant = 0.12; // 0.04
						double timingResConstant=0.00; // 0.001
						double timingResMinimum=0.005;              // TTS: 50ps;
						double timingResPower=-0.9; // -0.7
						double digitqforsmearing = (digitsq > 0.5) ? digitsq : 0.5;
						// based on a roughly similar form that gives ~60ps for Q~1 and ~5ps for Q~30
						double filedigittsmeared = // XXX if changing this ensure it is mirrored above.
							timingResConstant + timingConstant*pow(digitqforsmearing,timingResPower);
						// saturates @ 0.065 with majority of entries.
						if (filedigittsmeared < timingResMinimum) filedigittsmeared=timingResMinimum;
						// increasing power makes dropoff flatter, but also raises tail.
						// increasing prescale pulls distribution closer to y axis, 
						// making dropoff less right-angled...?
						// use wolfram alpha: 
						// to view overview
						// plot(y=max(0.00 + 0.12*(x^-0.9),0.005)) from x=0 to 30, y=0 to 0.065
						// to view Q~1
						// plot(y=max(0.00 + 0.12*(x^-0.9),0.005)) from x=0 to 3, y=0 to 0.12
						// to view tail
						// plot(y=max(0.00 + 0.12*(x^-0.9),0.005)) from x=25 to 30, y=0 to 0.01
						filedigittsmears.push_back(filedigittsmeared);
						digitsqlappdhist->Fill(digitsq);
						digitstlappdhist->Fill(absolutedigitst);
						digittsmearlappdhist->Fill(filedigittsmeared);
						lappdtimesmearvsqhist->Fill(filedigittsmeared,digitsq);
					}
				}
				
				// end of track level lappd digit analysis
				////////////////////////////////////////////////
#endif // ndef NOLAPPDS
*/
				
				//====================================================================================================
				// track level mrd digit analysis
				//====================================================================================================
				int anumMRDdigitsfromtrack=0;
				double aMRDchargefromtrack=0.;
				std::vector<int> mrdtubeshitbyatrack;
				
				// Loop over MRD digits
				// ====================
#ifdef VERBOSE
				cout<<"Analysing MRD digits"<<endl;
#endif // VERBOSE
				(atrigm!=nullptr) ? numMRDdigits = atrigm->GetCherenkovDigiHits()->GetEntries() : 0;
				for(Int_t i=0; i<numMRDdigits; i++){
					WCSimRootCherenkovDigiHit* thedigihit = 
						(WCSimRootCherenkovDigiHit*)atrigm->GetCherenkovDigiHits()->At(i);
					int digitstubeid = thedigihit->GetTubeId()-1;
					
#if FILE_VERSION<2
					// The CherenkovHitTimes is a flattened array (over PMTs) of arrays (over photons)
					// The  PhotonIds available from a digit are the indices within the subarray for that PMT
					// we therefore need we need to correct these indices with the start of the pmt's subarray
					// This may be found by scanning the CherenkovHits array (over PMTs),
					// in which the offset we need is stored as the 'GetTotalPe(0)' member
					int timeArrayOffset=-1;
					int ncherenkovhits=atrigm->GetNcherenkovhits();
					for(int ipmt = 0; ipmt < ncherenkovhits; ipmt++) {
						WCSimRootCherenkovHit* hitobject = 
							(WCSimRootCherenkovHit*)m->GetTrigger(0)->GetCherenkovHits()->At(ipmt);
						int tubeNumber = hitobject->GetTubeID()-1;
						if(tubeNumber==digitstubeid) {
							timeArrayOffset = hitobject->GetTotalPe(0);
							break;
						}
					}
#endif // FILE_VERSION<2
				
					// Calculate total charge from muon
					////////////////////////////////////
					std::vector<int> truephotonindices = thedigihit->GetPhotonIds();
					int numphotonsfromthismuon=0;
					for(int truephoton=0; truephoton<truephotonindices.size(); truephoton++){
						int thephotonsid = truephotonindices.at(truephoton);
#if FILE_VERSION<2
						thephotonsid+=timeArrayOffset;
#endif // FILE_VERSION<2
						WCSimRootCherenkovHitTime *thehittimeobject = 
							(WCSimRootCherenkovHitTime*)m->GetTrigger(0)->GetCherenkovHitTimes()->At(thephotonsid);
						Int_t thephotonsparenttrackid = thehittimeobject->GetParentID();
						if(thephotonsparenttrackid==nextrack->GetId()) {
							numphotonsfromthismuon++;
						}
					}
					
					// scale each digit charge by the fraction of hits in the digit which come from the muon
					double digitsq = thedigihit->GetQ();
					if(numphotonsfromthismuon!=0){
						aMRDchargefromtrack+= 
							digitsq * ((double)numphotonsfromthismuon/(double)truephotonindices.size());
						anumMRDdigitsfromtrack++;
						if(std::find(mrdtubeshitbyatrack.begin(), mrdtubeshitbyatrack.end(), 
										digitstubeid)==mrdtubeshitbyatrack.end())
							mrdtubeshitbyatrack.push_back(digitstubeid);
					}
				
				}  // end loop over mrd digits
				
				numMRDdigitsfromtrack.push_back(anumMRDdigitsfromtrack);
				MRDchargefromtrack.push_back(aMRDchargefromtrack);
				mrdtubeshitbytrack.push_back(mrdtubeshitbyatrack);
				
				// end of track level mrd digit analysis
				////////////////////////////////////////////////
				
				//====================================================================================================
				// track level veto digit analysis
				//====================================================================================================
				int anumFACCdigitsfromtrack=0;
				double aFACCchargefromtrack=0.;
				std::vector<int> facctubeshitbyatrack;
				
				// Loop over FACC digits
				// ====================
#ifdef VERBOSE
				cout<<"Analysing FACC digits"<<endl;
#endif // VERBOSE
				(atrigv!=nullptr) ? numFACCdigits = atrigv->GetCherenkovDigiHits()->GetEntries() : 0;
				for(Int_t i=0; i<numFACCdigits; i++){
					WCSimRootCherenkovDigiHit* thedigihit = 
						(WCSimRootCherenkovDigiHit*)atrigv->GetCherenkovDigiHits()->At(i);
					int digitstubeid = thedigihit->GetTubeId()-1;
					
#if FILE_VERSION<2
					// The CherenkovHitTimes is a flattened array (over PMTs) of arrays (over photons)
					// The  PhotonIds available from a digit are the indices within the subarray for that PMT
					// we therefore need we need to correct these indices with the start of the pmt's subarray
					// This may be found by scanning the CherenkovHits array (over PMTs), 
					// in which the offset we need is stored as the 'GetTotalPe(0)' member
					int timeArrayOffset=-1;
					int ncherenkovhits=atrigv->GetNcherenkovhits();
					for(int ipmt = 0; ipmt < ncherenkovhits; ipmt++) {
						WCSimRootCherenkovHit* hitobject = 
							(WCSimRootCherenkovHit*)v->GetTrigger(0)->GetCherenkovHits()->At(ipmt);
						int tubeNumber = hitobject->GetTubeID()-1;
						if(tubeNumber==digitstubeid) {
							timeArrayOffset = hitobject->GetTotalPe(0);
							break;
						}
					}
#endif // FILE_VERSION<2
				
					// Calculate total charge from muon
					////////////////////////////////////
					std::vector<int> truephotonindices = thedigihit->GetPhotonIds();
					int numphotonsfromthismuon=0;
					for(int truephoton=0; truephoton<truephotonindices.size(); truephoton++){
						int thephotonsid = truephotonindices.at(truephoton);
#if FILE_VERSION<2
						thephotonsid+=timeArrayOffset;
#endif // FILE_VERSION<2
						WCSimRootCherenkovHitTime *thehittimeobject = 
							(WCSimRootCherenkovHitTime*)v->GetTrigger(0)->GetCherenkovHitTimes()->At(thephotonsid);
						if(thehittimeobject==nullptr) cerr<<"HITTIME IS NULL"<<endl;
						Int_t thephotonsparenttrackid = (thehittimeobject) ? thehittimeobject->GetParentID() : -1;
						if(thephotonsparenttrackid==nextrack->GetId()) {
							numphotonsfromthismuon++;
						}
					}
					
					// scale each digit charge by the fraction of hits in the digit which come from the muon
					double digitsq = thedigihit->GetQ();
					if(numphotonsfromthismuon!=0){
						aFACCchargefromtrack+= 
							digitsq * ((double)numphotonsfromthismuon/(double)truephotonindices.size());
						anumFACCdigitsfromtrack++;
						if(std::find(facctubeshitbyatrack.begin(), facctubeshitbyatrack.end(), 
										digitstubeid)==facctubeshitbyatrack.end())
							facctubeshitbyatrack.push_back(digitstubeid);
					}
				
				}  // end loop over facc digits
				
				numFACCdigitsfromtrack.push_back(anumFACCdigitsfromtrack);
				FACCchargefromtrack.push_back(aFACCchargefromtrack);
				facctubeshitbytrack.push_back(facctubeshitbyatrack);
				
				// end of track level veto digit analysis
				////////////////////////////////////////////////
				
				//====================================================================================================
				// Print output info (don't fill until after setting remainder of event level variables)
				//====================================================================================================
//				for(std::map<std::string,bool>::iterator mapit=eventtypes.begin(); mapit!=eventtypes.end(); mapit++)
//					cout<<"eventtypes."<<mapit->first<<" = "<<mapit->second<<endl;
//				cout<<"isintank="<<isintank<<endl;
//				cout<<"isinfiducialvol="<<isinfiducialvol<<endl;
//				cout<<"eventq2="<<eventq2<<endl;
//				cout<<"eventEnu="<<eventEnu<<endl;
//				cout<<"neutrinopdg="<<neutrinopdg<<endl;
//				cout<<"muonentersMRD="<<muonentersMRD<<endl;
//				cout<<"muonstopsinMRD="<<muonstopsinMRD<<endl;
//				cout<<"muonrangesoutMRD="<<muonrangesoutMRD<<endl;
//				cout<<"mrdpenetrationcm="<<mrdpenetrationcm<<endl;
//				cout<<"mrdpenetrationlayers="<<mrdpenetrationlayers<<endl;
//				cout<<"mutracklengthinMRD="<<mutracklengthinMRD<<endl;
//				cout<<"atracklengthintank="<<atracklengthintank<<endl;
//				cout<<"tankchargefrommuon="<<tankchargefrommuon<<endl;
//				cout<<"fractionalchargeincone="<<fractionalchargeincone<<endl;
				
				// ----------------------------------------------------------------------------------------------
				// put things into vectors? Mostly depreciated, but still used for histograms?
				// ----------------------------------------------------------------------------------------------
				// this was after all the continue statements to record all muons that passed cuts
#ifndef NOGENIE
				Double_t scatteringangle = thegenieinfo.fslangle;
				Double_t neutrinoenergy = thegenieinfo.probeenergy;
				Double_t calculatedq2 = thegenieinfo.Q2;
#else // if defined NOGENIE
				Double_t scatteringangle = differencevector.Angle(TVector3(0,0,1));  // reco only, assume neutrino || z
				Double_t neutrinoenergy = -1;                                        // don't estimate here
				Double_t calculatedq2 = -1;                                          // don't estimate here
#endif // defined NOGENIE
				
				neutrinoenergiesvector.push_back(neutrinoenergy);
				primaryenergiesvector.push_back(trackstartE.at(track)); // RELATIVISTIC - INCLUDES REST MASS ENERGY
				scatteringanglesvector.push_back(scatteringangle);
				acceptedtrackids.push_back(nextrack->GetId());
				q2vector.push_back(calculatedq2);
				muonenergydepositions.push_back(tankchargefromthistrack);
			
			} // end of loop over tracks to analyze
			
#ifdef VERBOSE
			cout<<"End of track level analysis"<<endl;
#endif // VERBOSE
		
			//====================================================================================================
			// event level tank digit analysis
			//====================================================================================================
#ifdef VERBOSE
			cout<<"Analysing tank digits for event level info"<<endl;
#endif // VERBOSE
			numTankDigits = atrigt->GetCherenkovDigiHits()->GetEntries();
			if(numTankDigits==0){
				cerr<<"Event had no tank digits!!"<<endl;
				cerr<<"Event had "<<b->GetTrigger(0)->GetCherenkovHitTimes()->GetEntries()<<" cherenkov hits"<<endl;
				cerr<<"Event was a "<<interactiontypestring<<" interaction of a "<<eventEnu<<"GeV "<<PdgToString(neutrinopdg)
					<<" in the ";
				TVector3 primarystartvertex2(atrigt->GetVtx(0),atrigt->GetVtx(1),atrigt->GetVtx(2)); // same as neutrino stop vtx
				WCSimRootTrack* neutrinotrack = (WCSimRootTrack*)atrigt->GetTracks()->At(primaryneutrinoindex);
				if(neutrinotrack->GetStopvol()==20) cerr<<"Veto/Hall";
				else if(neutrinotrack->GetStopvol()==10) cerr<<"Tank";
				else if(neutrinotrack->GetStopvol()==30) cerr<<"MRD/Hall";
				else cerr<<"volume "<<(neutrinotrack->GetStartvol());
				cerr<<" with vertex material "<<vertexmaterial
					<<" and tank-centred event vertex ("<<primarystartvertex2.X()<<", "<<primarystartvertex2.Y()-tank_yoffset
					<<", "<<(primarystartvertex2.Z()-tank_start-tank_radius)<<")"
					<<",  r="<<sqrt(pow(primarystartvertex2.X(),2.)+pow(primarystartvertex2.Z()-tank_start-tank_radius,2.))
					<<endl<<"c.f. tank_start="<<tank_start<<", tank_radius="<<tank_radius<<endl;
				cerr<<"Event primaries were: "<<endl;
				for(int tracki=0; tracki<trackpdg.size(); tracki++){
					if(trackparenttype.at(tracki)!=0) continue;
					cerr<<trackstartE.at(tracki)<<"MeV "<<PdgToString(trackpdg.at(tracki))<<" in ";
					if(trackstartvol.at(tracki)==10) cerr<<"Tank";
					else if(trackstartvol.at(tracki)==20) cerr<<"FACC";
					else if(trackstartvol.at(tracki)==30) cerr<<"MRD";
					else cerr<<"Hall";
					cerr<<endl;
				}
				//assert(false);
			}
			
			// EventDistributions digit info
			totaltankcharge=atrigt->GetSumQ();
			upstreamcharge=0.;
			downstreamcharge=0.;
			topcapcharge=0.;
			bottomcapcharge=0.;
			// clear flat file digit vectors
			filedigitvertices.clear();
			filedigitQs.clear();
			filedigittsmears.clear();
			filedigitsensortypes.clear();
			filedigitPMTIDs.clear();
#ifdef LAPPD_DEBUG
			intileposxp->clear();
			intileposyp->clear();
			poserrxp->clear();
			poserryp->clear();
			poserrzp->clear();
			tileorientp->clear();
			octagonsidep->clear();
#endif
			
			//TH3D hitpositions = TH3D("hitpositions","hitpositions",100,-tank_radius,tank_radius, 100,tank_start,tank_start+2*tank_radius,100,-tank_halfheight+tank_yoffset,tank_halfheight+tank_yoffset);
			//ClearMapHistos(maphistos);
			
			// Loop over digits
			// ================
#ifdef VERBOSE
			cout<<"looping over "<<numTankDigits<<" tank digits for this event"<<endl;
#endif
			for(Int_t i=0; i<numTankDigits; i++){
				WCSimRootCherenkovDigiHit* thedigihit = 
					(WCSimRootCherenkovDigiHit*)atrigt->GetCherenkovDigiHits()->At(i);
				int digitstubeid = thedigihit->GetTubeId()-1;
				
				// add the digit charge to the counters of upstream/downstream/cap charge
				double digitsq = thedigihit->GetQ();
				WCSimRootPMT pmt = geo->GetPMT(digitstubeid);
				int thepmtsloc = pmt.GetCylLoc();
				switch (thepmtsloc){
					case 0: topcapcharge+=digitsq; break;
					case 2: bottomcapcharge+=digitsq; break;
					case 1: ((pmt.GetPosition(2)-tank_start-tank_radius)<0) ? upstreamcharge+=digitsq : downstreamcharge+=digitsq; break;
				}
				
				// map of all digit hits
				//FillTankMapHist(geo, digitstubeid, 2, maphistos, digitsq);
				
//				// Get digit first photon time if desired
//				//---------------------------------------
//				int firstphotonindex = thedigihit->GetPhotonIds().front();
//	#if FILE_VERSION<2
//				// to be able to match digits to HitTimes to get true times or parents 
//				// we need to correct the indices returned by GetPhotonIds with the offset of the
//				// CherenkovHitTime entry matching the first CherenkovHit entry on the same PMT...!
//				int timeArrayOffset=-1;
//				int ncherenkovhits=atrigt->GetNcherenkovhits();
//				for(int ipmt = 0; ipmt < ncherenkovhits; ipmt++) {
//					WCSimRootCherenkovHit* hitobject = 
//						(WCSimRootCherenkovHit*)b->GetTrigger(0)->GetCherenkovHits()->At(ipmt);
//					int tubeNumber = hitobject->GetTubeID()-1;
//					if(tubeNumber==digitstubeid) {
//						timeArrayOffset = hitobject->GetTotalPe(0);
//						break;
//					}
//				}
//				firstphotonindex += timeArrayOffset;
//	#endif // FILE_VERSION<2
//				double firstphotontime = 
//					(WCSimRootCherenkovHitTime*)b->GetTrigger(0)->GetCherenkovHitTimes()->At(thephotonsid)->GetTruetime();
				
				// add the digit info for simplified file format
				/////////////////////////////////////////////////
				filedigitPMTIDs.push_back(digitstubeid);
				double digitsx = pmt.GetPosition(0);
				double digitsy = pmt.GetPosition(1);
				double digitsz = pmt.GetPosition(2);
				double digitst = thedigihit->GetT(); // time within trigger window + triggeroffset
				WCSimRootEventHeader* trigheader=atrigt->GetHeader();
				double triggertime=trigheader->GetDate();
				double absolutedigitst=digitst-triggeroffset+triggertime;
				ROOT::Math::XYZTVector adigitvector = ROOT::Math::XYZTVector(digitsx,digitsy,digitsz,absolutedigitst);
				//cout<<"adding digit to flat file vectors"<<endl;
				filedigitvertices.push_back(adigitvector);
				filedigitQs.push_back(digitsq);
				std::string digitspst;
#if FILE_VERSION<4
				digitspst = "PMT8inch";
#else // if FILE_VERSION>=4
				int tubetypeindex = geo->GetTubeIndex(digitstubeid+1); // add back to get to ids from 1
				digitspst = geo->GetWCPMTNameAt(tubetypeindex);
				if(digitspst=="INDEX_OUT_OF_BOUNDS"){
					cerr<<"PMT TYPE NAME REQUEST IS OUT OF BOUNDS!?"<<endl;
					cerr<<"TubeId = "<<digitstubeid<<", TypeTypeIndex = "<<tubetypeindex<<endl;
					int pmi;
					for(pmi=0; ; pmi++){
						if(geo->GetWCPMTRadiusAt(pmi)<0) break;
					}
					cerr<<"counted "<<pmi<<" type counts"
						<<", cf. PMT name vector has "<<geo->GetPMTNames().size()<<" names"<<endl;
						assert(false);
				}
#endif // FILE_VERSION>=4
				filedigitsensortypes.push_back(digitspst);
				
				// calculate time resolution applied based on charge
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
				//time_PMT = time_true + PMT->HitTimeSmearing(Q);
				digitsqpmthist->Fill(digitsq);
				digitstpmthist->Fill(absolutedigitst);
				digittsmearpmthist->Fill(filedigittsmeared);
				pmttimesmearvsqhist->Fill(filedigittsmeared,digitsq);
				
#ifdef LAPPD_DEBUG
				intileposx.push_back(0);
				intileposy.push_back(0);
#if FILE_VERSION>3
				poserrx.push_back(0);
				poserry.push_back(0);
				poserrz.push_back(0);
#endif // FILE_VERSION<=4
				tileorient.push_back(0);
				octagonside.push_back(0);
#endif // defined LAPPD_DEBUG
				
				
			}  // end loop over digits
#ifdef VERBOSE
			cout<<"end loop over tank PMT digits"<<endl;
#endif
			
			// debug: check the charge maps
			/*
			if(hitpositions.GetEntries() > 10){
				TCanvas c1("c1");
				c1.cd();
				TH1* histowall=(TH1*)maphistos.at("tothistowall");
				histowall->Draw("colz");
				TCanvas c4("c4");
				c4.cd();
				TH1* histotop=maphistos.at("tothistotop");
				histowall->Draw("colz");
				TCanvas c8("c8");
				c8.cd();
				TH1* histobottom=maphistos.at("tothistobottom");
				histobottom->Draw("colz");
				TCanvas c5("c5");
				c5.cd();
				hitpositions.SetMarkerStyle(20);
				hitpositions.Draw();
				gPad->WaitPrimitive();
				//std::this_thread::sleep_for (std::chrono::seconds(15));  // wait so we can look at histos
			}
			*/
			
			// end of tank digit analysis
			////////////////////////////////////////////////
			
#ifndef NOLAPPDS
			//====================================================================================================
			// event level lappd digit analysis
			//====================================================================================================
#ifdef VERBOSE
			cout<<"Analysing LAPPD digits for event level info"<<endl;
#endif
			// FIXME: the vectors store all digits in the event regardless of which trigger they are in. 
			// we should scan through and only retrieve / fill jingbo file vectors with events within
			// the trigger window.
#ifdef VERBOSE
			cout<<"Looping over "<<lappd_numtileshitthisevt<<" LAPPD digits for this event"<<endl;
#endif
			int runningcount=0;
			for(int lappdi=0; lappdi<lappd_numtileshitthisevt; lappdi++){
				// loop over LAPPDs that had at least one hit
				int LAPPDID = lappd_hittile[lappdi];
#ifdef VERBOSE
				cout<<"Getting info for LAPPDID="<<LAPPDID<<endl;
#endif
				double tileposx = lappd_hittilesposx[lappdi];  // position of LAPPD in global coords
				double tileposy = lappd_hittilesposy[lappdi];  // IN MM
				double tileposz = lappd_hittilesposz[lappdi];
				WCSimRootPMT pmt = geo->GetLAPPD(LAPPDID-1);
				double pmtx = pmt.GetPosition(0);              // verified these are equivalent ^, BUT IN CM
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
						// use pmtx (in cm) to compare to Rthresh (in cm)
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
#ifdef VERBOSE
				cout<<"Looping over the hits on this LAPPD"<<endl;
#endif
				for(;runningcount<(lastrunningcount+numhitsthislappd); runningcount++){
					double peposx = lappd_hitpeposx.at(runningcount);      // position of hit on tile in LAPPD
					double peposy = lappd_hitpeposy.at(runningcount);
					//double peposz   = lappd_hitpeposz.at(runningcount);  // no apparent meaning
					double digitsx, digitsy, digitsz;
					switch (thepmtsloc){
						case 0: // top cap
							// use tileposx (in mm) to add to peposx (in mm)
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
#if FILE_VERSION>3
							<<"true global position is: ("<<lappd_hitglobalposx.at(runningcount)
							<<", "<<lappd_hitglobalposy.at(runningcount)
							<<", "<<lappd_hitglobalposz.at(runningcount)<<")"<<endl;
#endif // FILE_VERSION>3
					}
#endif // defined LAPPD_DEBUG
					double digitst  = lappd_hittruetime.at(runningcount);
#if FILE_VERSION>10 // dont know what file version this will be fixed in (check, has been fixed?)
					// but currently lappd digits store absolute digit time, NOT the time within a trigger....
					WCSimRootEventHeader* trigheader=atrigt->GetHeader();
					double triggertime=trigheader->GetDate();
					double absolutedigitst=digitst+triggeroffset-triggertime;
					// ALSO all lappd digits are stored regardless of being in or outside trigger window.
					// so this time may be way out, this digit may even be part of a different trigger.
					// some sort of loop over trigger times to sort it here....
					// (we could check if abs time > trigger time for next trig...)
#else // if FILE_VERSION<=10
					double absolutedigitst=digitst;
#endif // FILE_VERSION<=10
					ROOT::Math::XYZTVector adigitvector =  // convert mm to cm 
						ROOT::Math::XYZTVector(digitsx/10.,digitsy/10.,digitsz/10.,absolutedigitst);
				
#ifdef LAPPD_DEBUG
					intileposx.push_back(peposx);
					intileposy.push_back(peposy);
#if FILE_VERSION>3
					poserrx.push_back(digitsx-lappd_hitglobalposx.at(runningcount));
					poserry.push_back(digitsy-lappd_hitglobalposy.at(runningcount));
					poserrz.push_back(digitsz-lappd_hitglobalposz.at(runningcount));
#endif // FILE_VERSION>3
					tileorient.push_back(thepmtsloc);
					octagonside.push_back(theoctagonside);
#endif // defined LAPPD_DEBUG
#if FILE_VERSION>10 // don't know when this will be fixed
					double digitsq = lappd_hitcharge.at(runningcount);
#else // if FILE_VERSION<=10
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
//						int iflag;
//						SKIDigitizerThreshold(digitsq,iflag);
//						if(iflag!=0) continue; // this digit gets rejected by digitizer ... if it were a PMT
//						// TODO need to move all vector push_back's after this if we want to 'continue' 
//						// in order to keep all vectors the same size
//						digitsq*0.985;         // efficiency in WCSimWCDigitizer
#endif // FILE_VERSION<=10
					filedigitPMTIDs.push_back(LAPPDID+numpmts);
					filedigitvertices.push_back(adigitvector);
					filedigitQs.push_back(digitsq);
					filedigitsensortypes.push_back("lappd_v0"); // bad timing resolution used for pulses
					double timingConstant = 0.12; // 0.04
					double timingResConstant=0.00; // 0.001
					double timingResMinimum=0.005;              // TTS: 50ps;
					double timingResPower=-0.9; // -0.7
					double digitqforsmearing = (digitsq > 0.5) ? digitsq : 0.5;
					// based on a roughly similar form that gives ~60ps for Q~1 and ~5ps for Q~30
					double filedigittsmeared = // XXX if changing this ensure it is mirrored above.
						timingResConstant + timingConstant*pow(digitqforsmearing,timingResPower);
					// saturates @ 0.065 with majority of entries.
					if (filedigittsmeared < timingResMinimum) filedigittsmeared=timingResMinimum;
					// increasing power makes dropoff flatter, but also raises tail.
					// increasing prescale pulls distribution closer to y axis, 
					// making dropoff less right-angled...?
					// use wolfram alpha: 
					// to view overview
					// plot(y=max(0.00 + 0.12*(x^-0.9),0.005)) from x=0 to 30, y=0 to 0.065
					// to view Q~1
					// plot(y=max(0.00 + 0.12*(x^-0.9),0.005)) from x=0 to 3, y=0 to 0.12
					// to view tail
					// plot(y=max(0.00 + 0.12*(x^-0.9),0.005)) from x=25 to 30, y=0 to 0.01
					filedigittsmears.push_back(filedigittsmeared);
					digitsqlappdhist->Fill(digitsq);
					digitstlappdhist->Fill(absolutedigitst);
					digittsmearlappdhist->Fill(filedigittsmeared);
					lappdtimesmearvsqhist->Fill(filedigittsmeared,digitsq);
				} // end of loop over hits on this LAPPD
#ifdef VERBOSE
				cout<<"Done looping over hits on this LAPPD"<<endl;
#endif
			} // end of loop over LAPPDs with a hit
			// end of event level lappd digit analysis
			////////////////////////////////////////////////
#endif // ndef NOLAPPDS
			
			//====================================================================================================
			// event level mrd digit analysis
			//====================================================================================================
#ifdef VERBOSE
			cout<<"Analysing MRD digits for event level info"<<endl;
#endif
			totMRDcharge = (atrigm) ? atrigm->GetSumQ() : 0;
			
			// end of event level mrd digit analysis
			////////////////////////////////////////////////
			
			//====================================================================================================
			// event level facc digit analysis
			//====================================================================================================
#ifdef VERBOSE
			cout<<"Analysing FACC digits for event level info"<<endl;
#endif
			totFACCcharge = (atrigv) ? atrigv->GetSumQ() : 0;
			
			// end of event level facc digit analysis
			////////////////////////////////////////////////
			
			//====================================================================================================
			// Fill flat file for reconstruction dev with only the primary muon
			//====================================================================================================
#ifdef VERBOSE
			cout<<"Filling flat file"<<endl;
#endif
			// flat file stores the main muon form CCQE events
			if(filemuonindex!=-1){
				fileeventnum=wcsimeventnum;
				filetriggernum=wcsimtriggernum;
#ifndef NOGENIE
				fileneutrinoE=thegenieinfo.probeenergy;
				fileinteractiontypestring=thegenieinfo.procinfostring;
				// thegenieinfo.procinfostring gives format "<DIS - Weak[CC]>" for which symbols might not be ideal. Strip them.
				fileinteractiontypestring=fileinteractiontypestring.substr(1,fileinteractiontypestring.length() - 2);
				fileneutcode=thegenieinfo.neutinteractioncode;
				filemomtrans=thegenieinfo.Q2;
				filemuonenergy= thegenieinfo.fsleptonenergy;
				filemuonangle= thegenieinfo.fslangle;
#endif // !defined NOGENIE
				
				filemuonstartvertex=TLorentzVector(trackstartpos.at(filemuonindex),trackstarttime.at(filemuonindex));
				filemuonstopvertex=TLorentzVector(trackstoppos.at(filemuonindex),trackstoptime.at(filemuonindex));;
				filemuondirectionvector= (trackstoppos.at(filemuonindex)-trackstartpos.at(filemuonindex)).Unit();
				filepathlengthtotal=(trackstoppos.at(filemuonindex)-trackstartpos.at(filemuonindex)).Mag();
				filepathlengthinwater=tracklengthintank.at(filemuonindex);
				filepathlengthinmrd=tracklengthinMRD.at(filemuonindex);
				// best fit funciton as of time of writing
				TF1 MRDenergyvspenetration=TF1("af","expo(0)+pol0(2)+([3]/([4]-x))",0,1.6);
				MRDenergyvspenetration.SetParameters(-3.62645, 3.75503, 2.68525, 3.59244, 1.66969);
				double dEdx=MRDenergyvspenetration.Eval(filemuonangle);
				fileenergylossinmrd=filepathlengthinmrd*dEdx;
				
				vertextreenocuts->Fill();
				if(isinfiducialvol){
						vertextreefiducialcut->Fill();
						nummuontracksinfidvol++;
				}
			}
			
			// ===================================================================================
			// fill histograms and increment track counters
			// fill vertextreefiducialmrd if applicable
			// ===================================================================================
#ifdef VERBOSE
			cout<<"Filling histograms"<<endl;
			cout<<"filemuonindex="<<filemuonindex<<endl;
#endif
			
			Int_t indextouse=filemuonindex;
//			if(scatteringanglesvector.size()>1){
//				cerr<<"Mutliple accepted final state leptons!"<<endl;
//				for(int i=0; i<scatteringanglesvector.size(); i++){
//					cout<<"wcsim lepton angle "<<i<<" = "<<scatteringanglesvector.at(i)<<endl;
//					cout<<"wcsim lepton energy "<<i<<" = "<<primaryenergiesvector.at(i)<<endl;
//					cout<<"wcsim lepton trackid "<<i<<" = "<<acceptedtrackids.at(i)<<endl;
//				}
//#ifndef NOGENIE
//				cout<<"compare to:"<<endl;
//				cout<<"genie lepton angle = "<<thegenieinfo.fslangle<<endl;
//				cout<<"genie lepton energy = "<<thegenieinfo.fsleptonenergy<<endl;
//#endif // !defined NOGENIE
//				std::vector<Int_t>::iterator minit = std::min_element(acceptedtrackids.begin(), acceptedtrackids.end());
//				indextouse = std::distance(acceptedtrackids.begin(),minit);
//			}
			if(scatteringanglesvector.size()>0 && filemuonindex>0){
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
			} else {
#ifdef VERBOSE
				cout<<"no accepted wcsim track"<<endl;
#endif // VERBOSE
			}
			if(scatteringanglesvector.size()>0 && filemuonindex>0){
#ifndef NOGENIE
				// fill genie 'accepted' histograms.
				incidentneutrinoenergiesaccepted->Fill(thegenieinfo.probeenergy);
				incidentneutrinoanglesaccepted->Fill(thegenieinfo.probeangle);
				fslanglesaccepted->Fill(thegenieinfo.fslangle);
				fslenergiesaccepted->Fill(thegenieinfo.fsleptonenergy);
				eventq2accepted->Fill(thegenieinfo.Q2);
				if(eventtypes.at("IsWeakCC") && eventtypes.at("IsQuasiElastic"))
					neutrinovertexQEaccepted->Fill(thegenieinfo.Intx_x, thegenieinfo.Intx_y, thegenieinfo.Intx_z);
				// if it passes the fiducial cut, we've also accepted it, fill.
#endif // !defined NOGENIE
				if(isinfiducialvol){
#ifndef NOGENIE
					// genie values
					incidentneutrinoenergiesacceptedfidcut->Fill(thegenieinfo.probeenergy);
					fslanglesacceptedfidcut->Fill(thegenieinfo.fslangle);
					fslenergiesacceptedfidcut->Fill(thegenieinfo.fsleptonenergy);
					eventq2acceptedfidcut->Fill(thegenieinfo.Q2);
					if( (isinfiducialvol&&trackentersMRD.at(filemuonindex)) &&
						(eventtypes.at("IsWeakCC")) && 
						(eventtypes.at("IsQuasiElastic")) ){
						numCCQEneutrinoeventsinfidvolmrd++;
					}
#endif // !defined NOGENIE
					// fill flat tree for reconstruction dev
					if(trackentersMRD.at(filemuonindex)) vertextreefiducialmrd->Fill();
					nummuontracksinfidvolmrd+=scatteringanglesvector.size();
				}
			}
			
#ifdef TIME_EVENTS
			timer->Stop();
			cout<<"This event took "<<timer->RealTime()<<"ms, or "<<timer->CpuTime()<<"ms CPU equivalent time"<<endl; 
#endif // defined TIME_EVENTS
			
#ifdef VERBOSE
			cout<<"DONE WITH EVENT: filling EventDistributions tree"<<endl;
#endif
			treeout->Fill();
			//fileout->cd();
			//treeout->Write("",TObject::kOverwrite);
			
			//flateventfileout->cd();
			//vertextreenocuts->Write("",TObject::kOverwrite);
			//vertextreefiducialcut->Write("",TObject::kOverwrite);
			//vertextreefiducialmrd->Write("",TObject::kOverwrite);
			
		} // end of loop over WCSim Triggers
	}
	// end of loop over events
#ifdef VERBOSE
			cout<<"End of loop over events, writing files"<<endl;
#endif
	
	//====================================================================================================
	// Write trees to files
	//====================================================================================================
	fileout->cd();
//	treeout->Write();
	treeout->Write("",TObject::kOverwrite);
	
	flateventfileout->cd();
//	vertextreenocuts->Write();
//	vertextreefiducialcut->Write();
//	vertextreefiducialmrd->Write();
	vertextreenocuts->Write("",TObject::kOverwrite);
	vertextreefiducialcut->Write("",TObject::kOverwrite);
	vertextreefiducialmrd->Write("",TObject::kOverwrite);
	
	//====================================================================================================
	// Make histograms
	//====================================================================================================
#ifdef VERBOSE
	cout<<"Filling summary histograms"<<endl;
#endif
	gROOT->cd();
	
	cout<<"generating scaled histograms"<<endl;
	Double_t numbeamspills = totalpots/(4.0 * TMath::Power(10.,12.));
	Double_t numbeamspillsperday = (24.*60.*60.*1000.)/133.3333;	// 24 hours in ms / 133.33 ms between spills
	Double_t numdays = numbeamspills/numbeamspillsperday;
	cout<<"Results based on "<<totalpots<<" POTs, or "<<numbeamspills<<" beam spills, or "<<numdays<<" days of data"<<endl;
	cout<<"There were "<<numneutrinoeventsintank<<" neutrino interactions in the tank, of which "<<numCCQEneutrinoeventsintank<<" were true CCQE events."<<endl;
	cout<<"Of those, "<<numCCQEneutrinoeventsinfidvol<<" were within the fiducial volume."<<endl;
	cout<<"Of those in turn, "<<numCCQEneutrinoeventsinfidvolmrd<<" produced an accepted MRD muon"<<endl;
	cout<<"There were "<<nummuontracksintank<<" muons in the tank, of which "
		<<nummuontracksinfidvol<<" were from events in the fiducial volume."<<endl;
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
#endif // defined DRAW_HISTOS
	
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
	
	// draw other canvases
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
	
	debug1->Clear();
	aleg->Clear();
	aleg->AddEntry(incidentneutrinoanglesall,"Incident Neutrino Angles","l");
	aleg->AddEntry(incidentneutrinoanglesaccepted,"Incident Neutrino Angles Accepted","l");
	incidentneutrinoanglesall->SetLineColor(kRed);
	incidentneutrinoanglesall->Draw();
	incidentneutrinoanglesaccepted->SetLineColor(kBlue);
	incidentneutrinoanglesaccepted->Draw("same");
	aleg->Draw();
	debug1->SaveAs("neutrino_angles.png");
	
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
	delete debug3;
#endif  // defined WCSIMDEBUG
#endif  // defined DRAW_HISTOS
	
	//====================================================================================================
	// Write histograms to file and cleanup
	//====================================================================================================
	
	cout<<"writing flat event file"<<endl;
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
	
	// save, then delete canvases and legends
	cout<<"saving histogram images to '"<<outpath<<"' then deleting canvases and legends"<<endl;
	for(int i=0; i<canvaspointers.size(); i++){
		TCanvas* canv = canvaspointers.at(i);
		if(canv) canv->SaveAs(TString::Format("%s/histograms_%d.png",outpath,i));
		if(canv) delete canv; canv=0;
		TLegend* leg = legendpointers.at(i);
		if(leg) delete leg; leg=0;
	}
	
	// save, then delete scaled canvases and legends
	cout<<"saving scaled histogram images then deleting canvases and legends"<<endl;
	for(int i=0; i<scaledcanvaspointers.size(); i++){
		TCanvas* canv = scaledcanvaspointers.at(i);
		if(canv) canv->SaveAs(TString::Format("%s/scaled_histograms_%d.png",outpath,i));
		if(canv) delete canv; canv=0;
		TLegend* leg = scaledlegendpointers.at(i);
		if(leg) delete leg; leg=0;
	}
	
	// XXX Histogram deletion: since histograms are duplicated and pointers copied around, make sure to keep this updated!!! XXX
	// As of 27-09-17 this is correct: all points are stored either in histopointers, scaledhistopointers, or otherhistos
	// and each of those lists are unique
	// delete histograms - do this after the canvas version, which saves them before deleting them
	cout<<"writing histograms to histo file then deleting"<<endl;
	histofileout->cd();
	for(int i=0; i<histopointers.size(); i++){
		TH1D* temp = histopointers.at(i);
		if(temp) temp->Write("",TObject::kOverwrite);
		// some histograms are in the vector twice, so only delete it if it doesn't appear again
		if(std::count((histopointers.begin()+i), histopointers.end(), histopointers.at(i))==0) delete temp;
	}
	
	// delete scaled histograms
	cout<<"writing scaled histograms to histo file then deleting"<<endl;
	for(int i=0; i<scaledhistopointers.size(); i++){
		TH1D* temp = scaledhistopointers.at(i);
		if(temp) temp->Write("",TObject::kOverwrite);
		// some histograms are in the vector twice, so only delete it if it doesn't appear again
		if(std::count((scaledhistopointers.begin()+i), scaledhistopointers.end(), scaledhistopointers.at(i))==0) delete temp;
	}
	
	cout<<"writing otherhistos to histo file then deleting"<<endl;
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
#endif // defined WCSIMDEBUG
	(TH1*)tracklengthvsmuonlight,
	(TH1*)chargemap_nopions,
	(TH1*)inconehistowall,
	(TH1*)inconehistotop,
	(TH1*)inconehistobottom,
	(TH1*)outconehistowall,
	(TH1*)outconehistotop,
	(TH1*)outconehistobottom,
	(TH1*)tothistowall,
	(TH1*)tothistotop,
	(TH1*)tothistobottom
	};
	for(auto thehist : otherhistos){
		if(thehist) { thehist->Write(); delete thehist; }
	}
	cout<<"end writing"<<endl;
	
	//====================================================================================================
	// cleanup
	//====================================================================================================
	cout<<"cleanup"<<endl;
#ifndef PARTICLEGUNEVENTS
	// dirtfile trees:
	cout<<"resetting tankflux branches"<<endl;
	if(tankflux) tankflux->ResetBranchAddresses();
	cout<<"resetting tankmeta branches"<<endl;
	if(tankmeta) tankmeta->ResetBranchAddresses();
	cout<<"closing dirtfile"<<endl;
	if(dirtfile) dirtfile->Close(); // do we need to do this with a TChain?
	
#ifndef NOGENIE
	// geniefile tree:
	cout<<"resetting gtree branches"<<endl;
	if(gtree) gtree->ResetBranchAddresses();
	cout<<"closing geniefile"<<endl;
	if(geniefile) geniefile->Close();
	
	// associated genie objects declared with 'new'
	cout<<"deleting genierecordval"<<endl;
	if(genierecordval) delete genierecordval; genierecordval=0;
	cout<<"deleting nuprimarybranchval array"<<endl;
	if(nuprimarybranchval) delete[] nuprimarybranchval; nuprimarybranchval=0; // ? branch array
#endif // !defined NOGENIE
#endif // !defined PARTICLEGUNEVENTS
	
	// wcsim tree:
	cout<<"resetting wcsimT branches"<<endl;
#ifndef SPLITWCSIMFILES
	if(wcsimT) wcsimT->ResetBranchAddresses();
	cout<<"closing wcsimfile"<<endl;
	if(wcsimfile) wcsimfile->Close();
#else // if defined SPLITWCSIMFILES
	wcsimchain->ResetBranchAddresses();
	delete wcsimchain;
#endif // defined SPLITWCSIMFILES
	
	// lappd tree:
	cout<<"resetting lappdtree branches"<<endl;
#ifndef SPLITWCSIMFILES
	if(lappdtree) lappdtree->ResetBranchAddresses();
	cout<<"closing lappdfile"<<endl;
	if(lappdfile) lappdfile->Close();
#else // if defined SPLITWCSIMFILES
	if(lappdchain){ lappdchain->ResetBranchAddresses(); delete lappdchain; }
#endif // defined SPLITWCSIMFILES
	
	// The chain... this could either be tankflux (normal) or wcsimT (particle gun).
	// I think having reset both, this does not need to be dealt with?
	cout<<"resettng chain branches and deleting chain"<<endl;
	if(c) c->ResetBranchAddresses();
	//if(c) delete c; c=0;   // .... does it need to be deleted? ... leave it for application cleanup.
	
	// delete flat file output
	cout<<"deleting flat file output "<<flateventfileout->GetName()<<endl;
	if(flateventfileout){
		flateventfileout->Close(); delete flateventfileout; flateventfileout=0;
	}
	
	// delete histogram file out
	cout<<"deleting histofileout"<<endl;
	if(histofileout){ 
		histofileout->Close(); delete histofileout; histofileout=0;
	}
	
	// write and close file of event information
	cout<<"deleting EventDistributions fileout"<<endl;
	if(fileout){
		fileout->Close(); delete fileout; fileout=0;
	}
	
	cout<<"end of application. Goodbye."<<endl;
	gApplication->Terminate(); // so that it ends the job on the grid.
}

//====================================================================================================
//====================================================================================================
//====================================================================================================


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

#ifndef NOGENIE
void GetGenieEntryInfo(genie::EventRecord* gevtRec, genie::Interaction* genieint, GenieInfo &thegenieinfo, bool printneutrinoevent=false){
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
	///*Int_t*/ thegenieinfo.neutinteractioncode = genie::utils::ghep::NeutReactionCode(gevtRec);
	/*Int_t*/ thegenieinfo.nuanceinteractioncode  = genie::utils::ghep::NuanceReactionCode(gevtRec);
	/*TLorentzVector**/ thegenieinfo.IntxVtx = gevtRec->Vertex();
	/*Double_t*/ thegenieinfo.Intx_x = thegenieinfo.IntxVtx->X() * 100.;         // same info as nuvtx in g4dirt file
	/*Double_t*/ thegenieinfo.Intx_y = thegenieinfo.IntxVtx->Y() * 100.;         // GENIE uses meters
	/*Double_t*/ thegenieinfo.Intx_z = thegenieinfo.IntxVtx->Z() * 100.;         // GENIE uses meters
	/*Double_t*/ thegenieinfo.Intx_t = thegenieinfo.IntxVtx->T() * 1000000000;   // GENIE uses seconds
	
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
	/*Double_t*/ thegenieinfo.fslangle = thegenieinfo.k2->Vect().Angle(thegenieinfo.k1->Vect());
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
	
//	if(printneutrinoevent){
//		cout<<"   This was a "<< thegenieinfo.procinfostring <<" (neut code "<<thegenieinfo.neutinteractioncode
//			<<") interaction of a "
//			<<thegenieinfo.probeenergy<<"GeV " << thegenieinfo.probepartname << " on a "; 
//		
//		if( thegenieinfo.targetnucleonpdg==2212 || thegenieinfo.targetnucleonpdg==2112 ){
//			cout<<thegenieinfo.targetnucleonname<<" in a ";
//		} else {
//			cout<<"PDG-Code " << thegenieinfo.targetnucleonpdg<<" in a ";
//		}
//		
//		if( thegenieinfo.targetnucleusname!="unknown"){ cout<<thegenieinfo.targetnucleusname<<" nucleus, "; }
//		else { cout<<"Z=["<<thegenieinfo.targetnucleusZ<<","<<thegenieinfo.targetnucleusA<<"] nucleus, "; }
//		
//		if(remnucpos>-1){
//			cout<<endl<<"   producing a "<<thegenieinfo.remnantnucleusenergy<<"GeV "<<thegenieinfo.remnantnucleusname;
//		} else {
//			cout<<endl<<"   with no remnant nucleus";
//		}  // DIS on 16O produces no remnant nucleus?!
//		
//		if(fsleppos>-1){
//			cout<<" and a "<<thegenieinfo.fsleptonenergy<<"GeV "<<thegenieinfo.fsleptonname<<" final state lepton"<<endl;
//		} else{ cout<<" and no final state leptons"<<endl; }
//		
//		cout<<"   Q^2 was "<<thegenieinfo.Q2<<"(GeV/c)^2";
//		if(fsleppos>-1){
//			cout<<" and final state lepton was ejected at Cos(θ)="<<thegenieinfo.costhfsl<<endl;
//		}
//		cout<<"   Additional final state particles included:"<<endl
//			<< "     N(p) = "       << thegenieinfo.numfsprotons
//			<< "     N(n) = "       << thegenieinfo.numfsneutrons
//			<< "     N(pi^0) = "    << thegenieinfo.numfspi0
//			<< "     N(pi^+) = "    << thegenieinfo.numfspiplus
//			<< "     N(pi^-) = "    << thegenieinfo.numfspiminus
//			<<endl;
//	}
}

void PrintNeutrinoEvent(GenieInfo &thegenieinfo){
	cout<<"   This was a "<< thegenieinfo.procinfostring <<" (neut code "<<thegenieinfo.neutinteractioncode
		<<") interaction of a "
		<<thegenieinfo.probeenergy<<"GeV " << thegenieinfo.probepartname << " on a "; 
	
	if( thegenieinfo.targetnucleonpdg==2212 || thegenieinfo.targetnucleonpdg==2112 ){
		cout<<thegenieinfo.targetnucleonname<<" in a ";
	} else {
		cout<<"PDG-Code " << thegenieinfo.targetnucleonpdg<<" in a ";
	}
	
	if( thegenieinfo.targetnucleusname!="unknown"){ cout<<thegenieinfo.targetnucleusname<<" nucleus, "; }
	else { cout<<"Z=["<<thegenieinfo.targetnucleusZ<<","<<thegenieinfo.targetnucleusA<<"] nucleus, "; }
	
	//if(remnucpos>-1){
	if(thegenieinfo.remnantnucleusname!="n/a"){
		cout<<endl<<"   producing a "<<thegenieinfo.remnantnucleusenergy<<"GeV "<<thegenieinfo.remnantnucleusname;
	} else {
		cout<<endl<<"   with no remnant nucleus";
	}  // DIS on 16O produces no remnant nucleus?!
	
	//if(fsleppos>-1){
	if(thegenieinfo.fsleptonname!="n/a"){
		cout<<" and a "<<thegenieinfo.fsleptonenergy<<"GeV "<<thegenieinfo.fsleptonname<<" final state lepton"<<endl;
	} else{ cout<<" and no final state leptons"<<endl; }
	
	cout<<"   Q^2 was "<<thegenieinfo.Q2<<"(GeV/c)^2";
	//if(fsleppos>-1){
	if(thegenieinfo.fsleptonname!="n/a"){
		cout<<" and final state lepton was ejected at Cos(θ)="<<thegenieinfo.costhfsl<<endl;
	}
	cout<<"   Additional final state particles included:"<<endl
		<< "     N(p) = "       << thegenieinfo.numfsprotons
		<< "     N(n) = "       << thegenieinfo.numfsneutrons
		<< "     N(pi^0) = "    << thegenieinfo.numfspi0
		<< "     N(pi^+) = "    << thegenieinfo.numfspiplus
		<< "     N(pi^-) = "    << thegenieinfo.numfspiminus
		<<endl;
}
#endif

void FillTankMapHist(WCSimRootGeom* geo, int tubeID, int incone, std::map<std::string, TH2D*> &maphistos, double weight=1){
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
				if(incone==1){
					histotop=maphistos.at("inconehistotop");
				} else if(incone==0) {
					histotop=maphistos.at("outconehistotop");
				} else {
					histotop=maphistos.at("tothistotop");
				}
				if(histotop) histotop->Fill(thebins.first, thebins.second, weight);
			} else {cout<<"bad pmt: ID "<<tubeID<<" in CylLoc "<<pmt.GetCylLoc()<<endl;}
			break;
		}
		case 1: {
			if(wallpositionmap.count(tubeID)){
				std::pair<int,int> thebins = wallpositionmap.at(tubeID);
				TH2D* histowall;
				if(incone==1) {
					histowall=maphistos.at("inconehistowall");
				} else if(incone==0){
					histowall=maphistos.at("outconehistowall");
				} else {
					histowall=maphistos.at("tothistowall");
				}
				if(histowall) histowall->Fill(thebins.first+0.5, thebins.second, weight);
			} else {cout<<"bad pmt: ID "<<tubeID<<" in CylLoc "<<pmt.GetCylLoc()<<endl;}
			break;
		}
		case 2: {
			if(bottomcappositionmap.count(tubeID)){
				std::pair<int,int> thebins = bottomcappositionmap.at(tubeID);
				TH2D* histobottom;
				if(incone==1){
					histobottom=maphistos.at("inconehistobottom");
				} else if(incone==0){
					histobottom=maphistos.at("outconehistobottom");
				} else {
					histobottom=maphistos.at("tothistobottom");
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

//============================================================================

// a test to see if a projected point in a plane is within a box in that plane
int inline InBox( TVector3 Hit, TVector3 B1, TVector3 B2, const int Axis) {
	if ( Axis==1 && Hit.Z() > B1.Z() && Hit.Z() < B2.Z() && Hit.Y() > B1.Y() && Hit.Y() < B2.Y()) return 1;
	if ( Axis==2 && Hit.Z() > B1.Z() && Hit.Z() < B2.Z() && Hit.X() > B1.X() && Hit.X() < B2.X()) return 1;
	if ( Axis==3 && Hit.X() > B1.X() && Hit.X() < B2.X() && Hit.Y() > B1.Y() && Hit.Y() < B2.Y()) return 1;
	return 0;
}

// projects the hitpoint by adding a scaled vector to the start point
int inline GetIntersection( float fDst1, float fDst2, TVector3 P1, TVector3 P2, TVector3 &Hit) {
	if ( (fDst1 * fDst2) >= 0.0f) return 0;
	if ( fDst1 == fDst2) return 0; 
	Hit = P1 + (P2-P1) * ( -fDst1/(fDst2-fDst1) );
	return 1;
}

// returns true if line (L1, L2) intersects with the box (B1, B2)
// returns intersection point in Hit
bool CheckLineBox( TVector3 L1, TVector3 L2, TVector3 B1, TVector3 B2, TVector3 &Hit, TVector3 &Hit2, bool &error){
	error=false;
	
#ifdef MUTRACKDEBUG
	cout<<"CheckLineBox called with start:"; PrintVector(L1); cout<<" and stop "; PrintVector(L2);
	cout<<", c.f. box corners are "; PrintVector(B1); cout<<" and "; PrintVector(B2); cout<<endl;
#endif
	// check if it misses the box entirely by being on one side of a plane over entire track
	if (L2.X() <= B1.X() && L1.X() <= B1.X()) return false;
	if (L2.X() >= B2.X() && L1.X() >= B2.X()) return false;
	if (L2.Y() <= B1.Y() && L1.Y() <= B1.Y()) return false;
	if (L2.Y() >= B2.Y() && L1.Y() >= B2.Y()) return false;
	if (L2.Z() <= B1.Z() && L1.Z() <= B1.Z()) return false;
	if (L2.Z() >= B2.Z() && L1.Z() >= B2.Z()) return false;
#ifdef MUTRACKDEBUG
	cout<<"passes initial tests as not missing box"<<endl;
#endif
	
	// check if it's inside the box to begin with (classed as an interception at start vtx)
	if (L1.X() > B1.X() && L1.X() < B2.X() &&
		L1.Y() > B1.Y() && L1.Y() < B2.Y() &&
		L1.Z() > B1.Z() && L1.Z() < B2.Z()){
#ifdef MUTRACKDEBUG
		cout<<"starts in the box"<<endl;
#endif
		Hit = L1;
	}
	// check if it's inside the box at endpoint (classed as an interception at end vtx)
	if (L2.X() > B1.X() && L2.X() < B2.X() &&
		L2.Y() > B1.Y() && L2.Y() < B2.Y() &&
		L2.Z() > B1.Z() && L2.Z() < B2.Z()){
		Hit2 = L2;
#ifdef MUTRACKDEBUG
		cout<<"ends in the box"<<endl;
#endif
	}
	if(Hit==L1&&Hit2==L2) return true;
	
	// check for an interception in X, Y then Z.
	//if ( (GetIntersection( L1.X()-B1.X(), L2.X()-B1.X(), L1, L2, Hit) && InBox( Hit, B1, B2, 1 ))
	//  || (GetIntersection( L1.Y()-B1.Y(), L2.Y()-B1.Y(), L1, L2, Hit) && InBox( Hit, B1, B2, 2 ))
	//  || (GetIntersection( L1.Z()-B1.Z(), L2.Z()-B1.Z(), L1, L2, Hit) && InBox( Hit, B1, B2, 3 ))
	//  || (GetIntersection( L1.X()-B2.X(), L2.X()-B2.X(), L1, L2, Hit) && InBox( Hit, B1, B2, 1 ))
	//  || (GetIntersection( L1.Y()-B2.Y(), L2.Y()-B2.Y(), L1, L2, Hit) && InBox( Hit, B1, B2, 2 ))
	//  || (GetIntersection( L1.Z()-B2.Z(), L2.Z()-B2.Z(), L1, L2, Hit) && InBox( Hit, B1, B2, 3 )))
	//	return true;
	
	// Above seems to assume there will only be one interception!!
	// e.g. if X has an interception, there are no checks for Z interception - if it enters
	// the front face and exits the side, only the side exit will be returned. 
	// Instead, note all interception points and return the first (and second if it exists)
	std::vector<TVector3> interceptions;
	bool anyinterception=false;
	bool thisinterception;
	
#ifdef MUTRACKDEBUG
	cout<<"finding interceptions"<<endl;
#endif
	thisinterception=
	GetIntersection( L1.X()-B1.X(), L2.X()-B1.X(), L1, L2, Hit) && InBox(Hit, B1, B2, 1);
	if(thisinterception){ interceptions.push_back(Hit); anyinterception=true; }
	thisinterception=
	GetIntersection( L1.Y()-B1.Y(), L2.Y()-B1.Y(), L1, L2, Hit) && InBox( Hit, B1, B2, 2 );
	if(thisinterception){ interceptions.push_back(Hit); anyinterception=true; }
	thisinterception=
	GetIntersection( L1.Z()-B1.Z(), L2.Z()-B1.Z(), L1, L2, Hit) && InBox( Hit, B1, B2, 3 );
	if(thisinterception){ interceptions.push_back(Hit); anyinterception=true; }
	thisinterception=
	GetIntersection( L1.X()-B2.X(), L2.X()-B2.X(), L1, L2, Hit) && InBox( Hit, B1, B2, 1 );
	if(thisinterception){ interceptions.push_back(Hit); anyinterception=true; }
	thisinterception=
	GetIntersection( L1.Y()-B2.Y(), L2.Y()-B2.Y(), L1, L2, Hit) && InBox( Hit, B1, B2, 2 );
	if(thisinterception){ interceptions.push_back(Hit); anyinterception=true; }
	thisinterception=
	GetIntersection( L1.Z()-B2.Z(), L2.Z()-B2.Z(), L1, L2, Hit) && InBox( Hit, B1, B2, 3 );
	if(thisinterception){ interceptions.push_back(Hit); anyinterception=true; }
#ifdef MUTRACKDEBUG
	cout<<"found "<<interceptions.size()<<" interceptions";
	if(interceptions.size()>0) { cout<<" at "; PrintVector(interceptions.at(0)); }
	if(interceptions.size()>1) { cout<<" and "; PrintVector(interceptions.at(1)); }
	cout<<endl;
#endif
	
	if(interceptions.size()>2){
		cerr<<"CheckLineBox found more than two intercepts?! They are at:"<<endl;
		for(auto&& avec : interceptions)
			cerr<<"("<<avec.X()<<", "<<avec.Y()<<", "<<avec.Z()<<")"<<endl;
		error=true;
		//assert(false); // leave for later so we can print debug info.
		return false;
	} else if(interceptions.size()==2){
		auto vec1 = interceptions.at(0);
		auto vec2 = interceptions.at(1);
		if(vec1.Z()<vec2.Z()){
			Hit=vec1;
			Hit2=vec2;
		} else {
			Hit=vec2;
			Hit2=vec1;
		}
		return true;
	} else if(interceptions.size()==1) {
		if(Hit==L1){ // track starts in mrd - found intercept is exit point
			Hit2=interceptions.at(0);
			return true;
		} else {     // track starts outside mrd - found intercept is entry point
			Hit=interceptions.at(0);
			return true;
		}
	} else {
		return false;
	}
}

//=========================================================================================
// A test to see whether a projected line within a plane intersects a circle within the plane

bool CheckTankIntercepts( TVector3 primarystartvertex, TVector3 primarystopvertex, double avgtrackgradx, double avgtrackgrady, int trackstartvol, int trackstopvol, TVector3 &Hit, TVector3 &Hit2){
#ifdef MUTRACKDEBUG
	cout<<"CheckLineCircle with primarystartvertex=("<<primarystartvertex.X()<<", "<<primarystartvertex.Y()
		<<", "<<primarystartvertex.Z()<<"), primarystopvertex=("<<primarystopvertex.X()<<", "
		<<primarystopvertex.Y()<<", "<<primarystopvertex.Z()<<")"<<endl;
#endif
	
	// simple checks to save time
	if( 
		( (abs(primarystartvertex.X()) > tank_radius) && (abs(primarystopvertex.X()) > tank_radius) ) ||
		( (abs(primarystartvertex.Y()-tank_yoffset) > tank_halfheight) && (abs(primarystopvertex.Y()-tank_yoffset) > tank_yoffset) ) ||
		( (primarystartvertex.Z()>(tank_start+2.*tank_radius)) && (primarystopvertex.Z()>(tank_start+2.*tank_radius)) )
	){ return false; }
	
	// first check for the track being in the z plane, as this will produce infinite gradients
	if(abs(primarystopvertex.Z()-primarystartvertex.Z())<0.1){
		// we have a simpler case: the tank is simply a box of height tankheight
		// and width = length of a chord at the given z
#ifdef MUTRACKDEBUG
		cout<<"track is in z plane, using alt method"<<endl;
#endif
		double base = abs(primarystartvertex.Z()-tank_start-tank_radius);
		double chordlen = sqrt(pow(tank_radius,2.)-pow(base,2.)); // half the chord length, strictly. 
#ifdef MUTRACKDEBUG
		cout<<"half width of tank in plane z="<<primarystartvertex.Z()<<" is "<<chordlen<<endl;
#endif
		
		// simple checks to save time
		if(base>tank_radius) return false; // track misses tank in z plane
		if( ((primarystartvertex.Y()-tank_yoffset)<-tank_halfheight) &&
		    ((primarystopvertex.Y() -tank_yoffset)<-tank_halfheight) ){
		    return false; // track is entirely below tank
		} else 
		if( ((primarystartvertex.Y()-tank_yoffset)>tank_halfheight) &&
		    ((primarystopvertex.Y() -tank_yoffset)>tank_halfheight) ){
		    return false; // track is entirely above tank
		} else 
		if( (primarystartvertex.X()<-chordlen) &&
		    (primarystopvertex.Z( )<-chordlen) ){
		    return false; // track is entirely to left of tank
		} else 
		if( (primarystartvertex.X()>chordlen) &&
		    (primarystopvertex.X() >chordlen) ){
		    return false; // track is entirely to right of tank
		}
		
#ifdef MUTRACKDEBUG
		cout<<"track intercepted tank somewhere"<<endl;
#endif
		// second possibility: track is parallel to x-axis: then d*/dz=inf AND dx/dy=inf
		if((primarystopvertex.Y()-primarystartvertex.Y())<0.1){
#ifdef MUTRACKDEBUG
			cout<<"track ran parallel to x axis"<<endl;
#endif
			// even simpler: entry and exit points are just the tank walls in that z plane
			bool isrightgoing = (primarystopvertex.X()>primarystartvertex.X());
			if(trackstartvol!=10){
				(isrightgoing) ? Hit2.SetX(-chordlen) : Hit2.SetX(chordlen);
				Hit2.SetY(primarystartvertex.Y()); Hit2.SetZ(primarystartvertex.Z());
			}
			if(trackstopvol !=10){
				(isrightgoing) ? Hit.SetX(chordlen) : Hit.SetX(-chordlen);
				Hit.SetY(primarystartvertex.Y()); Hit.SetZ(primarystartvertex.Z());
			}
			
		} else {
			
			// track intercepted the tank at least once. find where
			double trackdxdy = (primarystopvertex.X()-primarystartvertex.X())/(primarystopvertex.Y()-primarystartvertex.Y());
#ifdef MUTRACKDEBUG
			cout<<"trackdxdy="<<trackdxdy<<endl;
#endif
			// there are 4 possible intersection points (2x, 2y), of which only 1 or 2 may be relevant
			double xentry=-1., xexit=-1., yentry=-1., yexit=-1.;
			bool entryfound=false, exitfound=false;
		
			// entry point
			if(trackstartvol!=10){ // if startvol is not in tank we must have had an entry point
#ifdef MUTRACKDEBUG
					cout<<"track started outside tank, finding entry point"<<endl;
#endif
				// find y entry - track must start outside tank y bounds to have a cap entry
				if(abs(primarystartvertex.Y()-tank_yoffset)>tank_halfheight){
					bool isupgoing = (primarystopvertex.Y()>primarystartvertex.Y());
					(isupgoing) ? yentry=-tank_halfheight : yentry=tank_halfheight;
					double distancetoyentry =  yentry - primarystartvertex.Y(); // signed
					cout<<"distancetoyentry="<<distancetoyentry<<endl;
					double xatyentry = primarystartvertex.X()+(trackdxdy*distancetoyentry);
					if(abs(xatyentry)<chordlen){ Hit2.SetX(xatyentry); Hit2.SetY(yentry+tank_yoffset); entryfound=true; }
					// else track entered via a wall
#ifdef MUTRACKDEBUG
					cout<<"track started outside caps ";
					if(entryfound) cout<<"and entered at ("<<Hit2.X()<<", "<<(Hit2.Y()-tank_yoffset)<<")"<<endl;
					else cout<<"but entered via the wall, as xatyentry="<<xatyentry<<endl;
#endif
				} // else track did not start outside tank y bounds: no cap entry
			
				// find wall entry - track must start outside tank x bounds to have a wall entry
				if((!entryfound) && abs(primarystartvertex.X())>chordlen){
					bool isrightgoing = (primarystopvertex.X()>primarystartvertex.X());
					(isrightgoing) ? xentry=-chordlen : xentry=chordlen;
					double distancetoxentry =  xentry - (primarystartvertex.X()); // signed
					double yatxentry = primarystartvertex.Y()+(distancetoxentry/trackdxdy);
					if(abs(yatxentry)<tank_halfheight){ Hit2.SetX(xentry); Hit2.SetY(yatxentry); entryfound=true; }
					// else track entered via a cap
#ifdef MUTRACKDEBUG
					cout<<"track started outside walls ";
					if(entryfound) cout<<"and entered at ("<<Hit2.X()<<", "<<(Hit2.Y()-tank_yoffset)<<")"<<endl;
					else cout<<"but yatxentry="<<yatxentry<<endl;
#endif
				} // else track did not start outside tank x bounds: no wall entry
			
				if(!entryfound) cerr<<"could not find track entry point!?"<<endl;
				else Hit2.SetZ(primarystartvertex.Z());
#ifdef MUTRACKDEBUG
				if(entryfound) cout<<"setting entry Z to "<<(Hit2.Z()-tank_start-tank_radius)<<endl;
#endif
			}
		
			// repeat for exit point
			if(trackstopvol!=10){
#ifdef MUTRACKDEBUG
					cout<<"track ended outside tank, finding exit point"<<endl;
#endif
				// find y exit - track must stop outside tank y bounds to have a cap exit
				if(abs(primarystopvertex.Y()-tank_yoffset)>tank_halfheight){
					bool isupgoing = (primarystopvertex.Y()>primarystartvertex.Y());
					(isupgoing) ? yexit=tank_halfheight : yexit=-tank_halfheight;
					double distancetoyexit =  yexit - primarystartvertex.Y(); // signed
					cout<<"distancetoyexit="<<distancetoyexit<<endl;
					double xatyexit = primarystartvertex.X()+(trackdxdy*distancetoyexit);
					if(abs(xatyexit)<chordlen){ Hit.SetX(xatyexit); Hit.SetY(yexit+tank_yoffset); exitfound=true; }
					// else track entered via a wall
#ifdef MUTRACKDEBUG
					cout<<"track stopped outside caps ";
					if(exitfound) cout<<"and exited at ("<<Hit.X()<<", "<<(Hit.Y()-tank_yoffset)<<")"<<endl;
					else cout<<"but exited via the wall, as xatyexit="<<xatyexit<<endl;
#endif
				} // else track did not start outside tank y bounds: no cap exit
		
				// find wall exit - track must start outside tank x bounds to have a wall exit
				if((!exitfound) && abs(primarystopvertex.X())>chordlen){
					bool isrightgoing = (primarystopvertex.X()>primarystartvertex.X());
					(isrightgoing) ? xexit=chordlen : xexit=-chordlen;
					double distancetoxexit =  xexit - (primarystartvertex.X()); // signed
					double yatxexit = primarystartvertex.Y()+(distancetoxexit/trackdxdy);
					if(abs(yatxexit)<tank_halfheight){ Hit.SetX(xexit); Hit.SetY(yatxexit); exitfound=true; }
					// else track entered via a cap
#ifdef MUTRACKDEBUG
					cout<<"track stopped outside walls ";
					if(entryfound) cout<<"and exited at ("<<Hit.X()<<", "<<(Hit.Y()-tank_yoffset)<<")"<<endl;
					else cout<<"but yatxexit="<<yatxexit<<endl;
#endif
				} // else track did not start outside tank y bounds: no wall exit
			
				if(!exitfound) cerr<<"could not find track exit point!?"<<endl;
				else Hit.SetZ(primarystartvertex.Z());
#ifdef MUTRACKDEBUG
				if(exitfound) cout<<"setting exit Z to "<<(Hit.Z()-tank_start-tank_radius)<<endl;
#endif
			}
			if( ((trackstopvol!=10)&&(!exitfound)) || ((trackstartvol!=10)&&(!entryfound)) ) return false;
		
		}
		
	} else {  // track was not in z plane: use old method, based on dx/dz, dy/dz
		
		// parameterize the track by it's z position, and find the 
		// z value where the tank would leave the radius of the tank
		Double_t xattankcentre = primarystartvertex.X() - (primarystartvertex.Z()-tank_start-tank_radius)*avgtrackgradx;
		Double_t firstterm = -avgtrackgradx*xattankcentre;
		Double_t thirdterm = 1+TMath::Power(avgtrackgradx,2.);
		Double_t secondterm = (TMath::Power(tank_radius,2.)*thirdterm) - TMath::Power(xattankcentre,2.);
		if(secondterm<=0){
#ifdef MUTRACKDEBUG
			cout<<"Tank miss!"<<endl;
#endif
			return false; // Doens't hit the tank.
		}
		Double_t solution1 = (firstterm + TMath::Sqrt(secondterm))/thirdterm;
		Double_t solution2 = (firstterm - TMath::Sqrt(secondterm))/thirdterm;
		Double_t tankstartpointz, tankendpointz;
		if(primarystopvertex.Z() > primarystartvertex.Z()){
			tankendpointz = solution1;    //forward going track
			tankstartpointz = solution2;  // more upstream intercept is entry point
		} else {
			tankendpointz = solution2;    // backward going track
			tankstartpointz = solution1;  // more upstream intercept is exit point
		}
		// correct for tank z offset
		tankendpointz += tank_start+tank_radius;
		tankstartpointz += tank_start+tank_radius;
#ifdef MUTRACKDEBUG
		// for sanity check:
		cout<<"z'1 = "<<solution1<<", z1 = "<<(solution1+primarystartvertex.Z())
			<<", x1 = "<<(primarystartvertex.X()+(avgtrackgradx*(solution1+primarystartvertex.Z())))
			<<", r1 = "<<sqrt(pow(solution1+primarystartvertex.Z(),2.)+pow((primarystartvertex.X()+(avgtrackgradx*(solution1+primarystartvertex.Z()))),2.))<<endl;
		cout<<"z'2 = "<<solution2<<", z2 = "<<(solution2+primarystartvertex.Z())
			<<", x2 = "<<(primarystartvertex.X()+(avgtrackgradx*(solution2+primarystartvertex.Z())))
			<<", r2 = "<<sqrt(pow(solution2+primarystartvertex.Z(),2.)+pow((primarystartvertex.X()+(avgtrackgradx*(solution2+primarystartvertex.Z()))),2.))<<endl;
#endif
		
		// calculate x by projecting track along parameterization
		Double_t tankendpointx = primarystartvertex.X() + (tankendpointz-primarystartvertex.Z())*avgtrackgradx;
		Double_t tankstartpointx = primarystartvertex.X() + (tankstartpointz-primarystartvertex.Z())*avgtrackgradx;
		
		// and same for y
		Double_t tankendpointy = primarystartvertex.Y() + (tankendpointz-primarystartvertex.Z())*avgtrackgrady;
		Double_t tankstartpointy = primarystartvertex.Y() + (tankstartpointz-primarystartvertex.Z())*avgtrackgrady;
		
#ifdef MUTRACKDEBUG
		cout<<"tank start: "<<tank_start<<endl;
		cout<<"tank end: "<<(tank_start+2*tank_radius)<<endl;
		cout<<"tank radius: "<<tank_radius<<endl;
		TVector3 trackdir = (primarystopvertex - primarystartvertex).Unit();
		cout<<"start dir =("<<trackdir.X()<<", "<<trackdir.Y()
			<<", "<<trackdir.Z()<<")"<<endl;
		cout<<"avgtrackgradx="<<avgtrackgradx<<endl;
		cout<<"avgtrackgrady="<<avgtrackgrady<<endl;
		cout<<"xattankcentre="<<xattankcentre<<endl;
		cout<<"firstterm="<<firstterm<<endl;
		cout<<"thirdterm="<<thirdterm<<endl;
		cout<<"secondterm="<<secondterm<<endl;
		cout<<"tank intercept z solution1="<<solution1<<endl;
		cout<<"tank intercept z solution2="<<solution2<<endl<<endl;
		
		cout<<"values before cap exit check"<<endl;
		cout<<"tankstartpointz="<<tankstartpointz<<endl;
		cout<<"tankstartpointx="<<tankstartpointx<<endl;
		cout<<"tankstartpointy="<<tankstartpointy<<endl;
		
		cout<<"tankendpointz="<<tankendpointz<<endl;
		cout<<"tankendpointx="<<tankendpointx<<endl;
		cout<<"tankendpointy="<<tankendpointy<<endl;
		
		cout<<"TMath::Abs(primarystartvertex.Y()-tank_yoffset)="
			<<TMath::Abs(primarystartvertex.Y()-tank_yoffset)
			<<"(tank_halfheight)="<<(tank_halfheight)<<endl;
#endif // MUTRACKDEBUG
		
		// now check if the particle would have exited through one of the caps before reaching this point
		// if projected y value at the radial exit is outside the tank, and the track started inside the tank,
		// then the tank must have left prior to this point via one of the caps
		if( TMath::Abs(tankendpointy-tank_yoffset)>(tank_halfheight) && 
			TMath::Abs(primarystartvertex.Y()-tank_yoffset)<(tank_halfheight)){
			// ^ second condition due to slight (mm) differences in tank y height in WCSim.
			// this trajectory exits through the cap. Need to recalculate x, z exiting points...!
			if(primarystopvertex.Y()>primarystartvertex.Y()){
				tankendpointy = tank_halfheight+tank_yoffset;	// by definition of leaving through cap
			} else {
				tankendpointy = -tank_halfheight+tank_yoffset;
			}
			tankendpointz = 
			primarystartvertex.Z() + (tankendpointy-primarystartvertex.Y())/avgtrackgrady;
			tankendpointx = 
			primarystartvertex.X() + (tankendpointz-primarystartvertex.Z())*avgtrackgradx;
		} // else this trajectory exited the tank by a side point; existing value is valid
		
		// repeat for the tankstartpoint
		if( TMath::Abs(tankstartpointy-tank_yoffset)>(tank_halfheight) && 
			TMath::Abs(primarystartvertex.Y()-tank_yoffset)<(tank_halfheight)){
			// ^ second condition due to slight (mm) differences in tank y height in WCSim.
			// this trajectory exits through the cap. Need to recalculate x, z exiting points...!
			if(primarystopvertex.Y()>primarystartvertex.Y()){
				tankstartpointy = tank_halfheight+tank_yoffset;	// by definition of leaving through cap
			} else {
				tankstartpointy = -tank_halfheight+tank_yoffset;
			}
			tankstartpointz = 
			primarystartvertex.Z() + (tankstartpointy-primarystartvertex.Y())/avgtrackgrady;
			tankstartpointx = 
			primarystartvertex.X() + (tankstartpointz-primarystartvertex.Z())*avgtrackgradx;
		} // else this trajectory entered the tank by a side point; existing value is valid
		
#ifdef MUTRACKDEBUG
		cout<<"values after cap exit check"<<endl;
		cout<<"tankstartpointz="<<tankstartpointz<<endl;
		cout<<"tankstartpointx="<<tankstartpointx<<endl;
		cout<<"tankstartpointy="<<tankstartpointy<<endl;
		
		cout<<"tankendpointz="<<tankendpointz<<endl;
		cout<<"tankendpointx="<<tankendpointx<<endl;
		cout<<"tankendpointy="<<tankendpointy<<endl;
#endif // defined MUTRACKDEBUG
		
		// return only the appropriate vertices
		if(trackstopvol!=10) Hit=TVector3(tankendpointx,tankendpointy,tankendpointz);
		if(trackstartvol!=10) Hit2=TVector3(tankstartpointx,tankstartpointy,tankstartpointz);
	}
	
	return true;
}

//=========================================================================================

// unused code - untested
// A test to see whether a projected line within a plane intersects a circle within the plane
bool CheckLineCircle( TVector3 trackstart, TVector3 trackend, TVector3 tankorigin, double tankradius, TVector3 &Hit, TVector3 &Hit2, bool &error){
	error=false;
	trackstart.SetY(0);
	trackend.SetY(0);
	TVector3 trackunitdir = (trackend-trackstart).Unit();
	
	double t0, t1; // distance along the line segment to intersections if present
	
	// geometric solution
	TVector3 L = tankorigin - trackstart;
	double tca = L.Dot(trackunitdir);
	//if( tca < 0 ) return false;                  // track points away from tank centre
	// but don't break yet: if the track started in the tank, we may still have an exit point!
	double d2 = L.Mag2() - tca*tca;             // square of distance of closest approach to tank centre
	double tankradius2 = pow(tankradius,2.);
	if (d2 > tankradius2) return false;         // track does not intercept tank
	double thc = sqrt(tankradius2 - d2);        // distance along line segment to point of closest approach
	t0 = tca - thc;
	t1 = tca + thc;
	
//	// analytic solution
//	Vec3f L = orig - center;
//	double a = dir.dotProduct(dir);
//	double b = 2 * dir.dotProduct(L);
//	double c = L.dotProduct(L) - tankradius2;
//	if (!solveQuadratic(a, b, c, t0, t1)) return false;
//	if (t0 > t1) std::swap(t0, t1);
	
	
//	t0<0 && t1<0 - outside and going away
//	t0>0 && t1>0 - outside and entering            << starts outside (check it enters)
//	t0<0 && t1>0 - inside and leaving              << starts inside  (check it leaves)
//	t0>0 && t1<0 - not possible
	
	double tracklength = (trackend-trackstart).Mag();
	TVector3 dir = (trackend-trackstart).Unit();
	if ((t0<0) && (t1<0))           { return false; }     // track starts outside and points away from tank
	if (t0>tracklength)             { return false; }     // track starts outside and ends before entering
	if ((t0>0) && (t1>tracklength)) { return false; }     // track starts inside and stops before exiting
	if ((t0>0) && (t1<0)) { cerr<<"t0>0 & t1<0 - CheckLineCircle impossible error!"<<endl; return false; }
	
	// remaining cases have at least one intercept
	if (t0>0) Hit  = trackstart + t0*dir;         // track entry point
	if (t1>0) Hit2 = trackstart + t1*dir;         // track exit point
	
	return true;
}

bool solveQuadratic(const double &a, const double &b, const double &c, double &x0, double &x1){
	double discr = b*b - 4*a*c;
	if(discr < 0){ return false; }
	else if(discr==0){ x0=-0.5*(b/a); x1=x0; }
	else {
		double q = (b>0) ? -0.5*(b+sqrt(discr)) : -0.5*(b-sqrt(discr));
		x0 = q/a;
		x1 = c/q;
	}
	if (x0>x1) std::swap(x0, x1);
	
	return true;
}

//=========================================================================================

std::string PdgToString(int code){
	if(pdgcodetoname.size()==0) GeneratePdgMap();
	if(pdgcodetoname.count(code)!=0){
		return pdgcodetoname.at(code);
	} else {
		cerr<<"unknown pdg code "<<code<<endl;
		return std::to_string(code);
	}
}

std::map<int,std::string>* GeneratePdgMap(){
	if(pdgcodetoname.size()!=0) return &pdgcodetoname;
	pdgcodetoname.emplace(2212,"Proton");
	pdgcodetoname.emplace(-2212,"Anti Proton");
	pdgcodetoname.emplace(11,"Electron");
	pdgcodetoname.emplace(-11,"Positron");
	pdgcodetoname.emplace(12,"Electron Neutrino");
	pdgcodetoname.emplace(-12,"Anti Electron Neutrino");
	pdgcodetoname.emplace(22,"Gamma");
	pdgcodetoname.emplace(2112,"Neutron");
	pdgcodetoname.emplace(-2112,"Anti Neutron");
	pdgcodetoname.emplace(-13,"Muon+");
	pdgcodetoname.emplace(13,"Muon-");
	pdgcodetoname.emplace(130,"Kaonlong");
	pdgcodetoname.emplace(211,"Pion+");
	pdgcodetoname.emplace(-211,"Pion-");
	pdgcodetoname.emplace(321,"Kaon+");
	pdgcodetoname.emplace(-321,"Kaon-");
	pdgcodetoname.emplace(3122,"Lambda");
	pdgcodetoname.emplace(-3122,"Antilambda");
	pdgcodetoname.emplace(310,"Kaonshort");
	pdgcodetoname.emplace(3112,"Sigma-");
	pdgcodetoname.emplace(3222,"Sigma+");
	pdgcodetoname.emplace(3212,"Sigma0");
	pdgcodetoname.emplace(111,"Pion0");
	pdgcodetoname.emplace(311,"Kaon0");
	pdgcodetoname.emplace(-311,"Antikaon0");
	pdgcodetoname.emplace(14,"Muon Neutrino");
	pdgcodetoname.emplace(-14,"Anti Muon Neutrino");
	pdgcodetoname.emplace(-3222,"Anti Sigma-");
	pdgcodetoname.emplace(-3212,"Anti Sigma0");
	pdgcodetoname.emplace(-3112,"Anti Sigma+");
	pdgcodetoname.emplace(3322,"Xsi0");
	pdgcodetoname.emplace(-3322,"Anti Xsi0");
	pdgcodetoname.emplace(3312,"Xsi-");
	pdgcodetoname.emplace(-3312,"Xsi+");
	pdgcodetoname.emplace(3334,"Omega-");
	pdgcodetoname.emplace(-3334,"Omega+");
	pdgcodetoname.emplace(-15,"Tau+");
	pdgcodetoname.emplace(15,"Tau-");
	pdgcodetoname.emplace(100,"OpticalPhoton");
	pdgcodetoname.emplace(3328,"Alpha");
	pdgcodetoname.emplace(3329,"Deuteron");
	pdgcodetoname.emplace(3330,"Triton");
	pdgcodetoname.emplace(3351,"Li7");
	pdgcodetoname.emplace(3331,"C10");
	pdgcodetoname.emplace(3345,"B11");
	pdgcodetoname.emplace(3332,"C12");
	pdgcodetoname.emplace(3350,"C13");
	pdgcodetoname.emplace(3349,"N13");
	pdgcodetoname.emplace(3340,"N14");
	pdgcodetoname.emplace(3333,"N15");
	pdgcodetoname.emplace(3334,"N16");
	pdgcodetoname.emplace(3335,"O16");
	pdgcodetoname.emplace(3346,"Al27");
	pdgcodetoname.emplace(3341,"Fe54");
	pdgcodetoname.emplace(3348,"Mn54");
	pdgcodetoname.emplace(3342,"Mn55");
	pdgcodetoname.emplace(3352,"Mn56");
	pdgcodetoname.emplace(3343,"Fe56");
	pdgcodetoname.emplace(3344,"Fe57");
	pdgcodetoname.emplace(3347,"Fe58");
	pdgcodetoname.emplace(3353,"Eu154");
	pdgcodetoname.emplace(3336,"Gd158");
	pdgcodetoname.emplace(3337,"Gd156");
	pdgcodetoname.emplace(3338,"Gd157");
	pdgcodetoname.emplace(3339,"Gd155");
	return &pdgcodetoname;
}

void PrintVector(TVector3 avec){
	cout<<"("<<avec.X()<<", "<<avec.Y()<<", "<<avec.Z()<<")";
}
