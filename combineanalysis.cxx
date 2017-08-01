/* vim:set noexpandtab tabstop=4 wrap */
/*TODO: to modify for new source files:
0. Re-generate and replace in this directory libWCSimRoot.so, .rootmap, .pcm files.
1. Re-enable #include RootOptions.hh
2. Disable timeArrayOffset lines.
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

// the WCSim analysis headers
#include "wcsimanalysis.hh"

// bonsai class based on this output root file 
#include "BonsaiEventClass.C" // TODO

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

void FillTankMapHist(WCSimRootGeom* geo, int tubeID, bool incone, std::map<std::string, TH2D*> &maphistos, double weight);
void ClearMapHistos(std::map<std::string,TH2D*> maphistos);	// clear the histograms

#include "genieinfo_struct.cxx"           // definition of a struct to hold genie info
// function to fill the into
void GetGenieEntryInfo(genie::EventRecord* gevtRec, genie::Interaction* genieint, GenieInfo& thegenieinfo);

const Float_t fidcutradius=tank_radius*0.8;     // fiducial volume is slightly lesss than the full tank radius
const Float_t fidcuty=50.;                      // a meter total fiducial volume in the centre y
const Float_t fidcutz=0;                        // fidcuial volume is before the tank centre.
const Float_t mu_rest_mass_E = 105.658;         // in MeV

// needed for drawing tank 2D map histograms
std::map<int, std::pair<int,int> > topcappositionmap;
std::map<int, std::pair<int,int> > bottomcappositionmap;
std::map<int, std::pair<int,int> > wallpositionmap;

const char* dirtpath="/pnfs/annie/persistent/users/moflaher/g4dirt";
const char* geniepath="/pnfs/annie/persistent/users/rhatcher/genie";
//const char* wcsimpath="/pnfs/annie/persistent/users/moflaher/wcsim";  // first 1M sample, various issues
//const char* wcsimpath="/annie/app/users/moflaher/wcsim/build";
const char* wcsimpath="/pnfs/annie/persistent/users/moflaher/wcsim_tankonly_17-06-17";
const char* wcsimlibrarypath="/annie/app/users/moflaher/wcsim/wcsim/libWCSimRoot.so";
//const char* outpath="/annie/app/users/moflaher/wcsim/root_work";
const char* outpath="/pnfs/annie/persistent/users/moflaher/wcsim_tankonly_17-06-17_ana";
//const char* analysispath="/annie/app/users/moflaher/wcsim/root_work";
const char* analysispath="/pnfs/annie/persistent/users/moflaher/wcsim_tankonly_17-06-17_ana";
const char* analysislibrarypath="/annie/app/users/moflaher/wcsim/root_work/analysiscaller_cxx.so";
//const char* bonsaipath="/annie/app/users/moflaher/bonsai/Bonsai_v0/results";
const char* bonsaipath="/pnfs/annie/persistent/users/moflaher/wcsim_tankonly_17-06-17_ana";
//const char* bonsaiclasspath="/annie/app/users/moflaher/bonsai/Bonsai_v0/cBonsaiEvent.C"; // TODO

const Bool_t printneutrinoevent=false;

double CalculateNeutrinoEnergy(double recoMuonEnergy, double recoMuonAngle);
double CalculateEventQ2(double recoMuonEnergy, double recoNeutrinoEnergy, double recoMuonAngle);

void truthtracks(){
	//ColourPlotStyle();
	
	// load WCSim library for reading WCSim files
	cout<<"loading "<<analysislibrarypath<<endl;
	gSystem->Load(analysislibrarypath);
	cout<<"loading "<<wcsimlibrarypath<<endl;
	//gSystem->Load(wcsimlibrarypath);
	//cout<<"loading "<<bonsaiclasspath<<endl; TODO
	//std::string loadclassstring = "\".L " + bonsaiclasspath + "\"";
	//gROOT->ProcessLine(loadclassstring.c_str());
	
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
	
	TFile* bonsaifile=0;
	TString bonsaifilepath;
	TTree* bonsaitree=0;
	Int_t numbonsaientries=0;
	
	TFile* mrdfile=0;
	TString mrdfilepath;
	TTree* mrdtree=0;
	Int_t nummrdentries=0;
	
	TFile* vetofile=0;
	TString vetofilepath;
	TTree* vetotree=0;
	Int_t numvetoentries=0;
	
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
	genie::NtpMCEventRecord* genierecordval = new genie::NtpMCEventRecord;
	TBranch* genierecordBranch=0;
	
	// wcsimT
	WCSimRootEvent* b=0, *m=0, *v=0;
	TBranch* bp=0, *mp=0, *vp=0;
	WCSimRootTrigger* atrigt=0, *atrigm=0, *atrigv=0;
	
	// bonsaitree
	cBonsaiEvent* bonsaievent = new cBonsaiEvent();
	
	// mrdtree
	Int_t numMrdEvents;
	Int_t numMrdTracks;
	TClonesArray* mrdevents = new TClonesArray("cMRDSubEvent");
	
	// vetotree
	//Int_t numVetoEvents;
	//TClonesArray* vetoevents = new TClonesArray("cVetoEvent");
	
	// geoT
	WCSimRootGeom* geo = 0; 
	
	// information from genie:
	GenieInfo thegenieinfo;
	
	// ==============================================================================================
	// ==============================================================================================
	// output file
	TFile* fileout = new TFile("CombinedAnalysis.root","RECREATE");
	fileout->cd();
	TTree* treeout = new TTree("recotree","Combined Event Reconstruction");
	
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
	int wcsimeventnum=-1;
	TBranch* bWCSimEventNum = treeout->Branch("WCSimEventNum",&wcsimeventnum);
	int wcsimtriggernum=-1;
	TBranch* bWCSimTriggerNum = treeout->Branch("WCSimTriggerNum",&wcsimtriggernum);
	TString mrdfilestring;
	TBranch* bMrdFileString = treeout->Branch("MrdTrackFile",&mrdfilestring);
	TString vetofilestring;
	TBranch* bVetoFileString = treeout->Branch("VetoEventFile",&vetofilestring);
	TString bonsaifilestring;
	TBranch* bBonsaiFileStrong = treeout->Branch("BonsaiFile",&bonsaifilestring);
	// genie information
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
	// wcsim info
	TVector3 mustartvtx(0,0,0);
	TBranch* bMuonStartVtx = treeout->Branch("MuonStartVertex",&mustartvtx);
	double mustarttime;
	TBranch* bMuonStartT = treeout->Branch("MuonStartT",&mustarttime);
	TVector3 mustopvtx(0,0,0);
	TBranch* bMuonStopVtx = treeout->Branch("MuonStopVertex",&mustopvtx);
	double muonangle=0.;
	TBranch* bMuonAngle = treeout->Branch("MuonAngle",&muonangle);
	double mustartE=0.;
	TBranch* bMuonStartE = treeout->Branch("MuonStartEnergy",&mustartE);
	double muendE=0.;
	TBranch* bMuonEndE = treeout->Branch("MuonEndEnergy",&muendE);
	double mutracklengthintank=0.; // keep this so we can do cuts on it, under the assumption our reconstruction will improve
	TBranch* bMuonTrackLengthInTank = treeout->Branch("MuonTrackLengthInTank",&mutracklengthintank);
	int numTankDigits=0;
	TBranch* bNumTankDigits = treeout->Branch("TotalTankDigits",&numTankDigits);
	double totaltankcharge=0.;
	TBranch* bTotalTankCharge = treeout->Branch("TotalTankCharge",&totaltankcharge);
	int numpizerotracks=-1;
	TBranch* bNumPiZeroTracks = treeout->Branch("NumPiZeroTracks",&numpizerotracks);
	int numpiplustracks=-1;
	TBranch* bNumPiPlusTracks = treeout->Branch("NumPiPlusTracks",&numpiplustracks);
	int numpiminustracks=-1;
	TBranch* bNumPiMinusTracks = treeout->Branch("NumPiMinusTracks",&numpiminustracks);
	int nummutracks=-1;
	TBranch* bNumMuTracks = treeout->Branch("NumMuTracks",&nummutracks);
	int numneutrontracks=-1;
	TBranch* bNumNeutronTracks = treeout->Branch("NumNeutronTracks",&numneutrontracks);
	int numprotontracks=-1;
	TBranch* bNumProtonTracks = treeout->Branch("NumProtonTracks",&numprotontracks);
	// bonsai info
	bool hadBonsaiEvent;
	TBranch* bHadBonsaiEvent = treeout->Branch("HadBonsaiEvent",&hadBonsaiEvent);
	TVector3 bonsaiVertex;
	TBranch* bBonsaiVertex = treeout->Branch("BonsaiVertex",&bonsaiVertex);
	double bonsaiVertexT;
	TBranch* bBonsaiVertexT = treeout->Branch("BonsaiVertexT",&bonsaiVertexT);
	double bonsaiVertexError;
	TBranch* bBonsaiVertexError = treeout->Branch("BonsaiVertexError",&bonsaiVertexError);
	bool bonsaiVertexInTank;
	TBranch* bBonsaiVertexInTank = treeout->Branch("BonsaiVertexInTank",&bonsaiVertexInTank);
	TVector3 bonsaiDirection;
	TBranch* bBonsaiDirection = treeout->Branch("BonsaiDirection",&bonsaiDirection);
	Double_t bonsaiDirectionError;
	TBranch* bBonsaiDirectionError = treeout->Branch("BonsaiDirectionError",&bonsaiDirectionError);
	TVector3 bonsaiTankExit;
	TBranch* bBonsaiTankExit = treeout->Branch("BonsaiTankExit",&bonsaiTankExit);
	bool bonsaiInterceptsMrd;
	TBranch* bBonsaiInterceptsMrd = treeout->Branch("BonsaiInterceptsMrd",&bonsaiInterceptsMrd);
	TVector3 bonsaiMrdEntry;
	TBranch* bBonsaiMrdEntry = treeout->Branch("BonsaiMrdEntry",&bonsaiMrdEntry);
	double bonsaiTrackLengthInTank;
	TBranch* bBonsaiTrackLengthInTank = treeout->Branch("BonsaiTrackLengthInTank",&bonsaiTrackLengthInTank);
	double bonsaiEnergyLoss;
	TBranch* bBonsaiEnergyLoss = treeout->Branch("BonsaiEnergyLoss",&bonsaiEnergyLoss);
	double bonsaiEnergyLossError;
	TBranch* bBonsaiEnergyLossError = treeout->Branch("BonsaiEnergyLossError",&bonsaiEnergyLossError);
	// mrd info
	bool hadMrdEvent; // may have had a bunch of digits but no reconstructable track
	TBranch* bHadMrdEvent = treeout->Branch("HadMrdEvent",&hadMrdEvent);
	bool hadMrdTrack;
	TBranch* bHadMrdTrack = treeout->Branch("HadMrdTrack",&hadMrdTrack);
	TVector3 mrdEntryVertex;
	TBranch* bMrdEntryVertex = treeout->Branch("MrdEntryVertex",&mrdEntryVertex);
	double mrdEntryTime;
	TBranch* bMrdEntryTime = treeout->Branch("MrdEntryTime",&mrdEntryTime);
	TVector3 mrdStopVertex;   // could be exit or stopping point
	TBranch* bMrdStopVertex = treeout->Branch("MrdStopVertex",&mrdStopVertex);
	double mrdDirectionError;
	TBranch* bMrdDirectionError = treeout->Branch("MrdDirectionError",&mrdDirectionError);
	double mrdTrackLength;
	TBranch* bMrdTrackLength = treeout->Branch("MrdTrackLength",&mrdTrackLength);
	double mrdEnergyLoss;
	TBranch* bMrdEnergyLoss = treeout->Branch("MrdEnergyLoss",&mrdEnergyLoss);
	double mrdEnergyLossError;
	TBranch* bMrdEnergyLossError = treeout->Branch("MrdEnergyLossError",&mrdEnergyLossError);
	bool mrdPenetrated;
	TBranch* bMrdPenetrated = treeout->Branch("MrdPenetrated",&mrdPenetrated);
	bool mrdStopped;
	TBranch* bMrdStopped = treeout->Branch("MrdStopped",&mrdStopped);
	bool mrdSideExit;
	TBranch* bMrdSideExit = treeout->Branch("MrdSideExit",&mrdSideExit);
	// veto info
	bool hadVetoEvent;
	TBranch* bHadVetoEvent = treeout->Branch("HadVetoEvent",&hadVetoEvent);
	double vetoEventTime;
	TBranch* bVetoEventTime = treeout->Branch("VetoEventTime",&vetoEventTime);
	// combined reconstruction info
	bool hadRecoEvent;
	TBranch* bHadRecoEvent = treeout->Branch("HadRecoEvent",&hadRecoEvent);
	bool recoMuonStopped;
	TBranch* bRecoMuonStopped = treeout->Branch("RecoMuonStopped",&recoMuonStopped);
	double recoMuonEnergy;
	TBranch* bRecoMuonEnergy = treeout->Branch("RecoMuonEnergy",&recoMuonEnergy);
	double recoMuonEnergyError;
	TBranch* bRecoMuonEnergyError = treeout->Branch("RecoMuonEnergyError",&recoMuonEnergyError);
	double recoMuonAngle;
	TBranch* bRecoMuonAngle = treeout->Branch("RecoMuonAngle",&recoMuonAngle);
	double recoMuonAngleError;
	TBranch* bRecoMuonAngleError = treeout->Branch("RecoMuonAngleError",&recoMuonAngleError);
	double recoEventQ2;
	TBranch* bRecoEventQ2 = treeout->Branch("RecoEventQ2",&recoEventQ2);
	double recoEventQ2Error;
	TBranch* bRecoEventQ2Error = treeout->Branch("RecoEventQ2Error",&recoEventQ2Error);
	double recoNeutrinoEnergy;
	TBranch* bRecoNeutrinoEnergy = treeout->Branch("RecoNeutrinoEnergy",&recoNeutrinoEnergy);
	double recoNeutrinoEnergyError;
	TBranch* bRecoNeutrinoEnergyError = treeout->Branch("RecoNeutrinoEnergyError",&recoNeutrinoEnergyError);
	TVector3 recoVertex;
	TBranch* bRecoVertex = treeout->Branch("RecoVertex",&recoVertex); // combined bonsai+MRD fit
	TVector3 recoVertexError;
	TBranch* bRecoVertexError = treeout->Branch("RecoVertexError",&recoVertexError);
	// compatibility testing
	TVector3 bonsaiMrdTankVtxDiff;
	TBranch* bBonsaiMrdTankVtxDiff = treeout->Branch("BonsaiMrdTankVtxDiff",&bonsaiMrdTankVtxDiff);
	double bonsaiMrdAngDiff;
	TBranch* bBonsaiMrdAngDiff = treeout->Branch("BonsaiMrdAngDiff",&bonsaiMrdAngDiff);
	
/*
	std::vector<TBranch*> thebranches{ bInTank, bEventType, bEventQ2, bEventEnu, bNeutrinoPdg, bInFidVol, bEventHasMuon, bMuonStartVtx, bMuonStopVtx, bMuonStartE, bMuonTrackLengthInTank, bMuonMrdPenetrationInCm, bMuonMrdPenetrationLayers, bMuonEntersMRD, bMuonStopsInMRD, bMuonRangesOutMRD, bMuonTrackLengthInMRD, bTankChargeFromMuon, bFractionOfMuonChargeInCone};
	int someit=0;
	bool haszombies=false;
	for(auto abranch : thebranches){
		if(abranch==0){ cout<<"branch "<<someit<<" is a zombie"<<endl; haszombies=true; }
		someit++;
	}
	assert(!haszombies&&"output file branches have zombies");
*/
	
	// ======================================================================================================
	// ======================================================================================================
//	// Just to test the inside/outside cherenkov cone algorithm
//	TH2D* inconehistowall = new TH2D("chargemap_incone_wall", "Charge Distribution Inside Cherenkov Cone (Wall)", pmtsperring+2,-1,pmtsperring+1,numpmtrings+2,-1,numpmtrings+1);
//	TH2D* inconehistotop = new TH2D("chargemap_incone_top","Charge Distribution Inside Cherenkov Cone (Top Cap)",caparraysize+2,-1,caparraysize+1,caparraysize+2,-1,caparraysize+1);
//	TH2D* inconehistobottom = new TH2D("chargemap_incone_bottom","Charge Distribution Inside Cherenkov Cone (Bottom Cap)",caparraysize+2,-1,caparraysize+1,caparraysize+2,-1,caparraysize+1);
//	
//	TH2D* outconehistowall = new TH2D("chargemap_outcone_wall", "Charge Distribution Outside Cherenkov Cone (Wall)", pmtsperring+2,-1,pmtsperring+1,numpmtrings+2,-1,numpmtrings+1);
//	TH2D* outconehistotop = new TH2D("chargemap_outcone_top","Charge Distribution Outside Cherenkov Cone (Top Cap)",caparraysize+2,-1,caparraysize+1,caparraysize+2,-1,caparraysize+1);
//	TH2D* outconehistobottom = new TH2D("chargemap_outcone_bottom","Charge Distribution Outside Cherenkov Cone (Bottom Cap)",caparraysize+2,-1,caparraysize+1,caparraysize+2,-1,caparraysize+1);
//	
//	std::map<std::string, TH2D*> maphistos;
//	maphistos.emplace("inconehistowall",inconehistowall);
//	maphistos.emplace("inconehistotop",inconehistotop);
//	maphistos.emplace("inconehistobottom",inconehistobottom);
//	maphistos.emplace("outconehistowall",outconehistowall);
//	maphistos.emplace("outconehistotop",outconehistotop);
//	maphistos.emplace("outconehistobottom",outconehistobottom);
	
	// ======================================================================================================
	// ======================================================================================================
	
	// done declaring file: move to loading and processing
	gROOT->cd();
	cout<<"loading first tankflux tree from "<<chainpattern<<" tchain"<<endl;
	c->LoadTree(0);
	Int_t treeNumber = -1;
	tankflux = c->GetTree();
	Int_t thistreesentries = tankflux->GetEntries();
	cout<<thistreesentries<<" entries in the first tree"<<endl;
	
	/*
	1.  Load next g4dirt entry
	2.  Load associated genie genie
	    If interaction is in tank:
	4.  WCSimRootEvent:
		4a.i.  Load true tracks, look for primary mu track, record track details.
	5a. cBonsaiEvent:
		5a.i.  Load events, note times and tank exit vertices.
		5a.i.  Save bonsai event info.
	5b. cMRDTrack:
		5b.i.  Load subevents, check for tracks.
		5b.ii. Load tracks, find the highest energy track consistent with bonsai event or tank if none.
		5b.ii. Save track details.
		5b.iii.Combine MRD + bonsai information to improve muon constraints.
	5c. cVetoEvent info
		5c.i.  Load events, note times.
		5c.ii. Save ...?
	6.  Reconstruct & save neutrino event info.
	7.  Close & cleanup.
	*/
	
	cout<<"looping over tchain entries"<<endl;
//	numents=10000;
	int maxfilenum=2000;
	Int_t wcsimTentry;
	// since WCSim only propagated tank events, there is no longer a 1:1 mapping between event numbers
	// in dirt files and WCSim files. As long as the selection criterion for dirt events is the same here
	// as in WCSim's PrimaryGeneratorAction, we can select the dirt files, and then just pull the next 
	// wcsim entry
	for(Int_t inputEntry=0; inputEntry<numents; inputEntry++){
		//===================================================================================================
		/* 	1. Load next g4dirt entry */ 
		//===================================================================================================
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
				cout<<"this genie file doesn't exist!"<<endl; 
				inputEntry += thistreesentries;	// skip loop iterator forward by num entries in this file
				continue; 
			}
			gtree = (TTree*)geniefile->Get("gtree");
			if(!gtree){cout<<"gtree doesn't exist!"<<endl; break; }
			numgenietentries = gtree->GetEntries();
			cout<<"gtree has "<<numgenietentries<<" entries in this file"<<endl;
			if(numgenietentries<1){cout<<"gtree has no entries!"<<endl; break; }
			
			// use regexp to pull out the file number needed for identifying the corresponding wcsim file
			std::match_results<string::const_iterator> submatches;
			// filename is of the form "annie_tank_flux.####.root"
			// #### is input file num. Need this to match against genie/wcsim file names
			std::regex theexpression (".*/[^0-9]+\\.([0-9]+)\\.root");
			cout<<"matching regex for filename "<<dirtfilename<<endl;
			std::regex_match (dirtfilename, submatches, theexpression);
			std::string submatch = (std::string)submatches[0];	// match 0 is 'whole match' or smthg
			if(submatch==""){ cout<<"unrecognised input file pattern: "<<dirtfilename<<endl; return; }
			submatch = (std::string)submatches[1];
			cout<<"extracted submatch is "<<submatch<<endl;
			int filenum = atoi(submatch.c_str());
			if(filenum>maxfilenum) break;
			
			// use filenum to open the corresponding wcsim file
			wcsimfilepath = TString::Format("%s/wcsim_0.%d.root",wcsimpath,filenum);
			cout<<"corresponding wcsim file is "<<wcsimfilepath<<endl;
			if(wcsimfile) wcsimfile->Close(); wcsimfile=0;
			wcsimfile = TFile::Open(wcsimfilepath);
			if(!wcsimfile){
				cout<<"wcsimfile "<<wcsimfilepath<<" doesn't exist!"<<endl; 
				inputEntry += thistreesentries;	// skip iterator forward by all the entries in this file
				continue; 
			}
			// load the geometry tree and grab the geometry if we haven't already
			if(geo==0){
				TFile* f = TFile::Open("/pnfs/annie/persistent/users/moflaher/wcsim_wdirt_17-06-17/wcsim_0.1000.root");
				TTree* geotree = (TTree*)f->Get("wcsimGeoT"); // TODO temporary override
				//TTree* geotree = (TTree*)wcsimfile->Get("wcsimGeoT");
				if(geotree==0){ cout<<"NO GEOMETRY IN FIRST FILE?"<<endl; assert(false); }
				geotree->SetBranchAddress("wcsimrootgeom", &geo);
				if (geotree->GetEntries() == 0) { cout<<"geotree has no entries!"<<endl; exit(9); }
				geotree->GetEntry(0);
				MakePMTmap(geo, topcappositionmap, bottomcappositionmap, wallpositionmap);
			}
			// load the next set of wcsim event info
			wcsimT = (TTree*)wcsimfile->Get("wcsimT");
			if(!wcsimT){cout<<"wcsimT doesn't exist!"<<endl; break; }
			numwcsimentries = wcsimT->GetEntries();
			cout<<"wcsimT has "<<numwcsimentries<<" entries in this file"<<endl;
			if(numwcsimentries==0){cout<<"wcsimT has no entries!"<<endl; break; }
			wcsimTentry=-1;
			
			// load reconstructed MRD info file
			mrdfilepath = TString::Format("%s/mrdtrackfile.%d.root",analysispath,filenum);
			//cout<<"corresponding mrd track file is "<<mrdfilepath<<endl;
			mrdfile = TFile::Open(mrdfilepath);
			if(!mrdfile){
				cout<<"mrdfile "<<mrdfilepath<<" doesn't exist!"<<endl; 
				inputEntry += thistreesentries;	// skip iterator forward by all the entries in this file
				continue; 
			}
			// load the set of reconstructed mrdtracks
			mrdtree = (TTree*)mrdfile->Get("mrdtree");
			if(!mrdtree){cout<<"mrdtree doesn't exist!"<<endl; break; }
			nummrdentries = mrdtree->GetEntries();
			if(nummrdentries==0){cout<<"mrdtree has no entries!"<<endl; break; }
			
			// load Veto event file
			/*
			vetofilepath = TString::Format("%s/vetotrackfile.%d.root",analysispath,filenum);
			//cout<<"corresponding veto event file is "<<wcsimfilepath<<endl;
			vetofile = TFile::Open(vetofilepath);
			if(!vetofile){
				cout<<"vetofile "<<vetofilepath<<" doesn't exist!"<<endl; 
				inputEntry += thistreesentries;	// skip iterator forward by all the entries in this file
				continue; 
			}
			// load the set of reconstructed mrdtracks
			vetotree = (TTree*)vetofile->Get("vetotree");
			if(!vetotree{cout<<"vetotree doesn't exist!"<<endl; break; }
			numvetoentries = vetotree->GetEntries();
			if(numvetoentries==0){cout<<"vetotree has no entries!"<<endl; break; }
			*/
			
			// load bonsai info
			bonsaifilepath = TString::Format("%s/bonsaiout.%d.root",bonsaipath,filenum);
			//cout<<"corresponding bonsai file is "<<bonsaifilepath<<endl;
			bonsaifile = TFile::Open(bonsaifilepath);
			if(!bonsaifile){
				cout<<"bonsaifile "<<bonsaifilepath<<" doesn't exist!"<<endl; 
				inputEntry += thistreesentries;	// skip iterator forward by all the entries in this file
				continue; 
			}
			// load the set of reconstructed mrdtracks
			bonsaitree = (TTree*)bonsaifile->Get("bonsaitree");
			if(!bonsaitree){cout<<"bonsaitree doesn't exist!"<<endl; break; }
			numbonsaientries = bonsaitree->GetEntries();
			//if(numbonsaientries==0){cout<<"bonsaitree has no entries!"<<endl; break; }
			
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
			wcsimT->SetBranchAddress("wcsimrootevent",&b, &bp);
			wcsimT->SetBranchAddress("wcsimrootevent_mrd",&m, &mp);
			wcsimT->SetBranchAddress("wcsimrootevent_facc",&v, &vp);
			bp->SetAutoDelete(kTRUE);
			mp->SetAutoDelete(kTRUE);
			vp->SetAutoDelete(kTRUE);
			if(bp==0||mp==0||vp==0){ cout<<"branches are zombies!"<<endl; }
			
			// gtree:
			gtree->SetBranchAddress("gmcrec",&genierecordval,&genierecordBranch);
			
			// mrdtree:
			mrdtree->SetBranchAddress("nummrdsubeventsthisevent",&numMrdEvents);
			mrdtree->SetBranchAddress("nummrdtracksthisevent",&numMrdTracks);
			mrdtree->SetBranchAddress("subeventsinthisevent",&mrdevents);
			
			// vetotree: TODO TODO TODO
			//vetotree->SetBranchAddress("numvetoeventsthisevent",&numVetoEvents);
			//vetotree->SetBranchAddress("vetoeventsinthisevent",&vetoevents);
			
			// bonsaitree:
			
			bonsaievent->Init(bonsaitree);
			bonsaievent->DisableBranches(); // reduces the amount of reading necessary
			
			treeNumber=nextTreeNumber;
		} // end of load new tree
		
		geniefilestring=geniefilepath;
		dirtfilestring=TString(dirtfilename);
		wcsimfilestring=wcsimfilepath;
		mrdfilestring=mrdfilepath;
		vetofilestring=vetofilepath;
		bonsaifilestring=bonsaifilepath;
		dirteventnum=localEntry;
		
		//===================================================================================================
		/* 2. retrieve genie info to test if vertex was in tank. TODO DISABLE FOR FILES W/DIRT NEUTRINOS */
		//===================================================================================================
#ifdef VERBOSE
		cout<<"processing inputEntry "<<inputEntry<<", localEntry "<<localEntry
		    <<"/"<<thistreesentries<<" in tree "<<treeNumber<<endl;
#endif
		nTankBranch->GetEntry(localEntry);
		vertexmaterialbranch->GetEntry(localEntry);
		// it doesn't make sense to count how many non-tank events there were if we didn't simulate them
		// they only contribute background, but they don't contribute background because they weren't simulated!
		// TODO remove this when analysing sample including dirt interactions.
		if(strcmp(vertexmaterial,"TankWater")!=0){ /*cout<<"neutrino vtx not in tank"<<endl;*/ continue; }
		if(nuprimarybranchval){delete[] nuprimarybranchval;}
		nuprimarybranchval = new Int_t[ntankbranchval];
		nuprimaryBranch->SetAddress(nuprimarybranchval);
		nuprimaryBranch->GetEntry(localEntry);
		
		Bool_t primariesinthisentry=false;
		for(int i=0;i<ntankbranchval;i++){
			if(nuprimarybranchval[i]==1){ primariesinthisentry=true; break; }
		}
		if(!primariesinthisentry){ cerr<<"wcsim primaries not genie primaries"<<endl; continue; }
		// These selection criteria are the WCSim PrimaryGeneratorAction ones. Any event that passes here
		// will have created a WCSimT entry:
		wcsimTentry++;   // do this now in case we introduce any 'continue' statements later
		isintank=true;
		
		//===================================================================================================
		/* 3. load remaining genie info. */
		//===================================================================================================
#ifdef VERBOSE
		cout<<"getting genie info"<<endl;
#endif
		if(localEntry>(numgenietentries-1)){ cerr<<"can't load localEntry "<<localEntry
								 <<" from "<<geniefilepath<<" gtree: not enough entries!"<<endl; continue; }
		genieentrybranch->GetEntry(localEntry);
		genierecordBranch->GetEntry(genieentry);
		genie::EventRecord* gevtRec = genierecordval->event;
		genie::Interaction* genieint = gevtRec->Summary();
		
		// fill thegenieinfo struct with all the genie info
		GetGenieEntryInfo(gevtRec, genieint, thegenieinfo);
		
		// fill the branches from the struct members
		genieeventnum=genieentry;
		eventtypes=thegenieinfo.eventtypes;
		// all the bools
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
		// the other info
		eventq2=thegenieinfo.Q2;
		eventEnu=thegenieinfo.probeenergy;
		neutrinopdg=thegenieinfo.probepdg;
		muonangle=thegenieinfo.fslanglegenie;
		hasmuon=(thegenieinfo.fsleptonname=="mu-");
		mustartE=thegenieinfo.fsleptonenergy;
		
		isinfiducialvol=false;
		if( (TMath::Sqrt(TMath::Power(thegenieinfo.genie_x, 2) 
			+ TMath::Power(thegenieinfo.genie_z-tank_start-tank_radius,2)) < fidcutradius) && 
			(TMath::Abs(thegenieinfo.genie_y-tank_yoffset) < fidcuty) && 
			((thegenieinfo.genie_z-tank_start-tank_radius) < fidcutz) ){
			isinfiducialvol=true;
		}
		
		//===================================================================================================
		/* 4. load wcsim detector response. */
		//===================================================================================================
#ifdef VERBOSE
		cout<<"getting wcsim entry "<<wcsimTentry<<endl;
#endif
		if(wcsimTentry>(numwcsimentries-1)){ cerr<<"can't load wcsimT entry "<<wcsimTentry
				<<" from "<<wcsimfilepath<<" wcsimT - not enough entries!"<<endl; continue; }
		wcsimT->GetEntry(wcsimTentry);
		// read only first subtrigger; delayed events are not interesting for primary FSL tracks
		atrigt = b->GetTrigger(0);
		atrigm = m->GetTrigger(0);
		atrigv = v->GetTrigger(0);
		
		wcsimeventnum=wcsimTentry;
		
		Int_t numtracks = atrigt->GetNtrack();
#ifdef VERBOSE
		cout<<"wcsim event had "<<numtracks<<" truth tracks"<<endl;
#endif
		numpizerotracks=0;
		numpiplustracks=0;
		numpiminustracks=0;
		nummutracks=0;
		numneutrontracks=0;
		numprotontracks=0;
		
		// pull just the primary muon track - there should only be one.
		// (there are on average ~1.4 muons per event!)
		Int_t mutrackindex=-1;
		for(int track=0; track<numtracks; track++){
			WCSimRootTrack* nextrack = (WCSimRootTrack*)atrigt->GetTracks()->At(track);
			Int_t primarypdg = nextrack->GetIpnu();
			if(nextrack->GetFlag()!=0) continue;
			switch (primarypdg){
				case 111: numpizerotracks++; break;
				case 211: numpiplustracks++; break;
				case -211: numpiminustracks++; break;
				case 13: nummutracks++; break;
				case 2112: numneutrontracks++; break;
				case 2212: numprotontracks++; break;
			}
			if(TMath::Abs(primarypdg)!=13) continue;                // not a muon
			if(nextrack->GetParenttype()!=0) continue;              // not a primary
			mutrackindex=track;
			break;
		}
		if(mutrackindex<0){
		  // this should be the same as hasmuon==false set by genie
		  // record tank/mrd events even for DIS/CC1Pi
		  mustartvtx=TVector3(0,0,0);
		  mustarttime=0.;
		  mustopvtx=TVector3(0,0,0);
		  muendE=0;
		  mutracklengthintank=0;
		} else {
			
			//===============================================================================================
			/* 6. Retrieve additional primary muon details from WCSim event. */
			//===============================================================================================
			WCSimRootTrack* nextrack = (WCSimRootTrack*)atrigt->GetTracks()->At(mutrackindex);
			
			TLorentzVector primarystartvertex(  nextrack->GetStart(0),
												nextrack->GetStart(1),
												nextrack->GetStart(2),
												nextrack->GetTime());
			TLorentzVector primarystopvertex(   nextrack->GetStop(0),
												nextrack->GetStop(1),
												nextrack->GetStop(2),
												nextrack->GetStopTime());
			
			mustartvtx=primarystartvertex.Vect();
			mustarttime=primarystartvertex.T();
			mustopvtx=primarystopvertex.Vect();
			TVector3 differencevector  = (primarystopvertex.Vect()-primarystartvertex.Vect());
			muendE=nextrack->GetEndE();
			
			// ----------------------------------------------------------------------------------------------
			// calculate the track length in water - TODO remove once tank reco is good
			// ----------------------------------------------------------------------------------------------
			// to calculate track length in water find distance from start vertex to the point
			// where it intercepts the tank. if this length > total track length, use total track length
			// otherwise use this length. 
			
			// calculate track angles - assuming a straight track
			Float_t oppx = primarystopvertex.X() - primarystartvertex.X();
			Float_t adj = primarystopvertex.Z() - primarystartvertex.Z();
			Float_t avgtrackanglex = (oppx/adj);
			Float_t oppy = primarystopvertex.Y() - primarystartvertex.Y();
			Float_t avgtrackangley = (oppy/adj);
			
			// first find out the z value where the tank would leave the radius of the tank
#ifdef MUTRACKDEBUG
			cout<<"z0 = "<<genie_z-tank_start-tank_radius<<", x0 = "<<genie_x<<endl;
#endif
			Double_t xatziszero = 
			(thegenieinfo.genie_x - (thegenieinfo.genie_z-tank_start-tank_radius)*(avgtrackanglex));
			Double_t firstterm = -(avgtrackanglex)*xatziszero;
			Double_t thirdterm = 1+TMath::Power((avgtrackanglex),2.);
			Double_t secondterm = (TMath::Power(tank_radius,2.)*thirdterm) - TMath::Power(xatziszero,2.);
			Double_t solution1 = (firstterm + TMath::Sqrt(secondterm))/thirdterm;
			Double_t solution2 = (firstterm - TMath::Sqrt(secondterm))/thirdterm;
			Double_t tankendpointz;
			if(primarystopvertex.Z() > primarystartvertex.Z()){
				tankendpointz = solution1;	//forward going track
			} else {
				tankendpointz = solution2;	// backward going track
			}
			// correct for tank z offset (do after tankendpointx, before tankendpointy)
			tankendpointz += tank_start+tank_radius;
			Double_t tankendpointx = thegenieinfo.genie_x + (tankendpointz-thegenieinfo.genie_z)*(avgtrackanglex);
			// now check if the particle would have exited through one of the caps before reaching this radius
			Double_t tankendpointy = 
			thegenieinfo.genie_y + (tankendpointz-thegenieinfo.genie_z)*(avgtrackangley);
			
#ifdef MUTRACKDEBUG
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
				thegenieinfo.genie_z + (tankendpointy-thegenieinfo.genie_y)/(avgtrackangley);
				tankendpointx = 
				thegenieinfo.genie_x + (tankendpointz-thegenieinfo.genie_z)*(avgtrackanglex);
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
				cerr<<"Track length is impossibly long!"<<endl;
				assert(false);
			}
			if(TMath::IsNaN(mutracklengthintank)){
				cerr<<"NaN result from mu track length in tank!"<<endl;
				assert(false);
			}
			
			// ----------------------------------------------------------------------------------------------
			// digit analysis
			// ----------------------------------------------------------------------------------------------
#ifdef VERBOSE
			cout<<"Analysing tank digits"<<endl;
#endif
			numTankDigits = atrigt->GetCherenkovDigiHits()->GetEntries();
			totaltankcharge=atrigt->GetSumQ();
			
			/* just incase we want to draw events, leave this here.
			// it needs work to be enabled.
			ClearMapHistos(maphistos);
			for(Int_t i=0; i<numTankDigits; i++){
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
				
				WCSimRootPMT pmt = geo->GetPMT(digitstubeid);
				double digitsx = pmt.GetPosition(0);
				double digitsy = pmt.GetPosition(1);
				double digitsz = pmt.GetPosition(2);
				int thepmtsloc = pmt.GetCylLoc();
				switch (thepmtsloc){
					case 0: topcapcharge+=digitsq; break;
					case 2: bottomcapcharge+=digitsq; break;
					case 1: ((digitsz-tank_start-tank_radius)<0) ? upstreamcharge+=digitsq : downstreamcharge+=digitsq; break;
					FillTankMapHist(geo, digitstubeid, true, maphistos, digitsq);
				}
			}  // end loop over digits
			
			TH1* histowall=(TH1*)maphistos.at("inconehistowall");
			if(histowall.GetEntries() > 10){
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
			
		}  // end if if WCSim true tracks had a primary muon
		
		// ==================================================================================================
		// load bonsai event
		// ==================================================================================================
		
		bonsaievent->GetEntry(wcsimTentry);
		// check this is the correct event
		/*
		if(
		bool correctevent=false;
			//bonsaievent->wcsimfilestring==wcsimfilestring&& //<<TODO: fix in bonsai files.
			bonsaievent->EventId==wcsimeventnum&&
			bonsaievent->SubtriggerId==0 //<<TODO: fix in bonsai files: all 32163???
			)
		correctevent=true;
		if(!correctevent){
			cerr<<"bonsai event did not match wcsimT event!"<<endl
				<<"wcsimfile: "<<wcsimfilestring
				<<", bonsai's wcsimfile: "<<bonsaievent->wcsimfilestring
				<<", wcsimTentry: "<<wcsimTentry
				<<", bonsai's event entry: "<<bonsaievent->EventId
				<<", wcsimTentry subevent: "<<0
				<<", bonsai subentry = "<<bonsaievent->SubtriggerId
				<<endl;
				assert(false);
		}
		*/
		// for the file:
		hadBonsaiEvent=bonsaievent->Vertex_Found;
		if(hadBonsaiEvent){
			// bonsai only returns one vertex, so fill the info now
			// - no need to wait until we see if it's mrd matched
			bonsaiVertex               = bonsaievent->RecoVtx_Vec;
			bonsaiVertexError          = 80.; // TODO ?? how to calculate?
			bonsaiVertexT              = bonsaievent->RecoVtx_T;
			bonsaiVertexInTank         = bonsaievent->Vertex_In_Tank;
			bonsaiDirection            = bonsaievent->RecoDir_Vec;
			bonsaiDirectionError       = bonsaievent->RecoDir_Error;
			bonsaiTankExit             = bonsaievent->Tank_Exit_Point; // could be updated if mrd track found
			TVector3 diffvec           = bonsaiTankExit-bonsaiVertex;
			bonsaiTrackLengthInTank    = diffvec.Mag();
			bonsaiEnergyLoss           = bonsaievent->TankEnergyLoss; // [MeV]
			bonsaiEnergyLossError      = bonsaievent->TankEnergyLossError;
			bonsaiInterceptsMrd        = bonsaievent->Intercepts_MRD;
			bonsaiMrdEntry             = bonsaievent->Projected_MRDEntryPoint;
		} else {
			bonsaiVertex               = TVector3(0,0,0);
			bonsaiVertexError          = 0;
			bonsaiVertexT              = 0;
			bonsaiVertexInTank         = false;
			bonsaiDirection            = TVector3(0,0,0);
			bonsaiDirectionError       = 0;
			bonsaiTankExit             = TVector3(0,0,0);
			bonsaiTrackLengthInTank    = 0;
			bonsaiEnergyLoss           = 0;
			bonsaiEnergyLossError      = 0;
			bonsaiInterceptsMrd        = false;
			bonsaiMrdEntry             = TVector3(0,0,0);
		}
		
		// ==================================================================================================
		// load mrd event
		// ==================================================================================================
		
		TDatabasePDG db;
		Double_t muonmass = (db.GetParticle(13)->Mass())*1000.;      // converted to MeV
		// need this when calculating neutrino event info, which uses relativistic muon energy
		
		if(numMrdEvents>0) hadMrdEvent=true;
		// in case we find no suitable tracks, pre-fill with blanks
		mrdEntryVertex=TVector3(0,0,0);
		mrdEntryTime= 0;
		mrdStopVertex=TVector3(0,0,0);
		mrdDirectionError=0;
		mrdTrackLength=0;
		mrdEnergyLoss=0;
		mrdEnergyLossError=0;
		mrdPenetrated=false;
		mrdStopped=false;
		mrdSideExit=false;
		recoVertex=TVector3(0,0,0);
		recoVertexError=TVector3(0,0,0);
		bonsaiMrdAngDiff=0.;
		
		TVector3 mrdtankvertex;
		
		int matchedsubevent=-1;
		int matchedtrack=-1;
		
		if(numMrdTracks>0){ // quick skip
			hadMrdTrack=true;
			mrdevents->Clear();
			mrdtree->GetEntry(wcsimTentry);
			assert(numMrdEvents==mrdevents->GetEntriesFast()
				&&"Num MRD SubEvents in TClonesArray does not match claimed number!");
			
			// need to scan over all mrd tracks and find the "best" match.
			// this will be defined as the longest track, consistent with bonsai event if available,
			// or just with a tank interception otherwise.
			double maxtracklength=-1.;
			for(int subev=0; subev<numMrdEvents; subev++){
				// TClonesArray->At() returns a pointer
				cMRDSubEvent* asubevent =(cMRDSubEvent*)mrdevents->At(subev);
				std::vector<cMRDTrack> &mrdtracks = *(asubevent->GetTracks());
				
				// loop over tracks in the mrd subevent
				for(int tracki=0; tracki<mrdtracks.size(); tracki++){
					auto anmrdtrack = mrdtracks.at(tracki);
					
					// is this track longer than our current max?
					if(anmrdtrack.GetTrackLength()<maxtracklength) continue;
					
					// if we had a bonsai event, is it consistent with the bonsai vertex?
					if(hadBonsaiEvent){
						// Look for an overlap region... this gets tricky.
						
						// we know there will be an overlap region if the mrd track comes within 
						// bonsaiVertexError of the bonsaiVertex.
						double closestapp = anmrdtrack.GetClosestApproach(bonsaiVertex, 0);
						if(closestapp<bonsaiVertexError){
							// this is a simple case - best fit mrd track is consistent with bonsai region
							
							// next require the direction is compatible
							TVector3 mrdtrackdir = anmrdtrack.GetStopVertex() - anmrdtrack.GetStartVertex();
							mrdDirectionError = anmrdtrack.GetTrackAngleError();
							bonsaiMrdAngDiff = TMath::ACos(bonsaiDirection.Dot(mrdtrackdir.Unit()));
							if(bonsaiMrdAngDiff>(bonsaiDirectionError+mrdDirectionError)) continue;
							
							// this is a new longest consistent track!
							// the projected mrd vertex will be the point of closest approach
							mrdtankvertex = anmrdtrack.GetClosestPoint(bonsaiVertex, 0);
							
							// the global best fit point is between bonsai vertex and projected mrd vertex
							recoVertex = bonsaiVertex + 0.5*(bonsaiVertex-mrdtankvertex);
							
							// some rough estimates of error: TODO how to define?
							TVector3 vertexdiff = bonsaiVertex - mrdtankvertex;
							double vtxdiffmag = vertexdiff.Mag();
							double xmax, xmin, ymax, ymin;
							anmrdtrack.GetProjectionLimits(recoVertex.Z(),xmax, xmin, ymax, ymin);
							double mrdvtxerrmag = sqrt(pow(xmax-xmin,2)+pow(ymax-ymin,2));
							double errormag = 0.5*sqrt((vtxdiffmag+bonsaiVertexError)*
														(vtxdiffmag+mrdvtxerrmag));
							recoVertexError=TVector3(errormag,errormag,errormag); // XXX???
							
							// since this is new best match: retrieve/override details
							maxtracklength     = anmrdtrack.GetTrackLength();
							mrdEntryVertex     = anmrdtrack.GetMrdEntryPoint();
							mrdEntryTime       = anmrdtrack.GetStartTime();
							mrdStopVertex      = anmrdtrack.GetStopVertex();
							mrdTrackLength     = anmrdtrack.GetTrackLength();
							mrdEnergyLoss      = anmrdtrack.GetEnergyLoss();
							mrdEnergyLossError = anmrdtrack.GetEnergyLossError();
							mrdPenetrated      = anmrdtrack.GetIsPenetrating();
							mrdStopped         = anmrdtrack.GetIsStopped();
							mrdSideExit        = anmrdtrack.GetIsSideExit();
							
							matchedsubevent=subev;
							matchedtrack=tracki;
						
						} else {
							// even if best fit line does not cross, we can still have consistency, provided
							// the region spanned by mrd track error intersects the bonsai allowed region.
							// This region exists if either mrd track upper or lower boundary lines 
							// approaches within bonsaiVertexError of bonsaiVertex.
							
							// we have 4 corners that define the square based pryamid projection
							// from mrd. TODO: How do we find the closest point?
							// we can make an estimate by taking 2D case in the plane of bonsai vertex Z.
							// We can then find the closest point fairly easily from the projection corners:
							double xmax, xmin, ymax, ymin;
							anmrdtrack.GetProjectionLimits(recoVertex.Z(), xmax, xmin, ymax, ymin);
							double xclosest, yclosest;
							if(xmax>bonsaiVertex.X()&&xmin<bonsaiVertex.X()){
								xclosest=bonsaiVertex.X();
							} else {
								if(abs(bonsaiVertex.X()-xmin)>abs(bonsaiVertex.X()-xmax)){
									xclosest=xmax;
								} else {
									xclosest=xmin;
								}
							}
							if(ymax>bonsaiVertex.Y()&&ymin<bonsaiVertex.Y()){
								yclosest=bonsaiVertex.Y();
							} else {
								if(abs(bonsaiVertex.Y()-ymin)>abs(bonsaiVertex.Y()-ymax)){
									yclosest=ymax;
								} else {
									yclosest=ymin;
								}
							}
							// check the closest point within the mrd allowed region is within
							// the bonsai error radius
							double radiusofclosest = sqrt(pow(bonsaiVertex.X()-xclosest,2)
														 +pow(bonsaiVertex.Y()-yclosest,2));
							if(radiusofclosest>bonsaiVertexError) continue;
							mrdtankvertex =TVector3(xclosest,yclosest,bonsaiVertex.Z());
							
							// vertices are consistent: next require the direction is compatible
							TVector3 mrdtrackdir = anmrdtrack.GetStopVertex() - anmrdtrack.GetStartVertex();
							mrdDirectionError = anmrdtrack.GetTrackAngleError();
							bonsaiMrdAngDiff = TMath::ACos(bonsaiDirection.Dot(mrdtrackdir.Unit()));
							if(bonsaiMrdAngDiff>(bonsaiDirectionError+mrdDirectionError)) continue;
							
							// the best fit point will be somewhere within the overlap region. 
							// TODO: for now let's place it midway: halfway between closest point 
							// (representing mrd outer limit) and bonsaiVertexError (representing 
							// bonsai outer limit).
							double magofdistance = radiusofclosest + ((bonsaiVertexError - radiusofclosest)/2.);
							TVector3 dirofdistance = mrdtankvertex - bonsaiVertex;
							dirofdistance.SetMag(magofdistance);
							recoVertex = bonsaiVertex+dirofdistance;
							
							// TODO: error calculation... 
							double errormag = 0.5*sqrt((magofdistance+bonsaiVertexError)*
														(magofdistance/*+mrdvtxerrmag*/)); //XXX???
							recoVertexError = TVector3(errormag,errormag,errormag); // XXX???
							
							// since this is new best match: retrieve/override details
							maxtracklength     = anmrdtrack.GetTrackLength();
							mrdEntryVertex     = anmrdtrack.GetMrdEntryPoint();
							mrdEntryTime       = anmrdtrack.GetStartTime();
							mrdStopVertex      = anmrdtrack.GetStopVertex();
							mrdTrackLength     = anmrdtrack.GetTrackLength();
							mrdEnergyLoss      = anmrdtrack.GetEnergyLoss();
							mrdEnergyLossError = anmrdtrack.GetEnergyLossError();
							mrdPenetrated      = anmrdtrack.GetIsPenetrating();
							mrdStopped         = anmrdtrack.GetIsStopped();
							mrdSideExit        = anmrdtrack.GetIsSideExit();
							
							matchedsubevent=subev;
							matchedtrack=tracki;
						}
					} else {
						// no bonsai event
						
						if(anmrdtrack.GetInterceptsTank()){
							// this is a simple case - best fit mrd track is consistent with tank region
							// this is a new longest consistent track!
							
							// we find the entry and exit points of the tank:
							TVector3 entrypoint;
							TVector3 exitpoint;
							bool valid = anmrdtrack.CheckTankIntercept(&entrypoint, &exitpoint,0);
							recoVertex = entrypoint+0.5*(exitpoint-entrypoint);
							
							// TODO i don't even know.
							double xmax, xmin, ymax, ymin;
							anmrdtrack.GetProjectionLimits(recoVertex.Z(),xmax, xmin, ymax, ymin);
							double errormag = sqrt(pow(entrypoint.Z()-exitpoint.Z(),2)
												  +pow(xmax-xmin,2)+pow(ymax-ymin,2));
							recoVertexError=TVector3(errormag,errormag,errormag);
							
							// n.b. no directionality check required, but still save error
							mrdDirectionError = anmrdtrack.GetTrackAngleError();
							
							// since this is new best match: retrieve/override details
							maxtracklength     = anmrdtrack.GetTrackLength();
							mrdEntryVertex     = anmrdtrack.GetMrdEntryPoint();
							mrdEntryTime       = anmrdtrack.GetStartTime();
							mrdStopVertex      = anmrdtrack.GetStopVertex();
							mrdTrackLength     = anmrdtrack.GetTrackLength();
							mrdEnergyLoss      = anmrdtrack.GetEnergyLoss();
							mrdEnergyLossError = anmrdtrack.GetEnergyLossError();
							mrdPenetrated      = anmrdtrack.GetIsPenetrating();
							mrdStopped         = anmrdtrack.GetIsStopped();
							mrdSideExit        = anmrdtrack.GetIsSideExit();
							
							matchedsubevent=subev;
							matchedtrack=tracki;
						
						} else {
//							// TODO THIS IS TOO COMPLICATED RIGHT NOW. NO VERTEX FOUND. 
//							// even if best fit line does not cross, we can still have consistency, provided
//							// the region spanned by mrd track error intersects the tank.
//							
//							// we need to check all 4 corners: max and min of h projection
//							// each with max and min of v projection. Up to 2 of them may intersect the tank.
//							// TODO 
//							// TVector3 entrypoint;
//							// TVector3 exitpoint;
//							// bool corner1 = anmrdtrack.CheckTankIntercept(&entrypoint, &exitpoint,1);
//							
//								
//							// if neither are within tank_radius of tank origin, there's no overlap region
//							if((maxclosestapp>tank_radius)&&(minclosestapp>tank_radius)) continue;
//							
//							// n.b. no directionality check required, but still save error
//							mrdDirectionError = anmrdtrack.GetTrackAngleError();
//							
//							// since the best fit point is inconsistent, only one line can intersect
//							double closestapp = std::min(minclosestapp, maxclosestapp);
//							int maxmin = (closestapp==minclosestapp) ? 1 : -1;
//							
//							// TODO what we *should* do is find the smallest angle from best fit line that
//							// intercepts the tank, and find where it does so.
//							// since it's less complicated, just do as per best fit line, 
//							// but with whichever limit intercepted the tank.
//							TVector3 entrypoint;
//							TVector3 exitpoint;
//							bool valid = anmrdtrack.CheckTankIntercept(&entrypoint, &exitpoint,maxmin);
//							recoVertex = entrypoint+0.5*(exitpoint-entrypoint);
//							
//							TVector3 tankorigin=(0,0,tank_start+tank_radius);
//							double bestfitclosestapp = anmrdtrack.GetClosestApproach(tankorigin, 0);
//							recoVertexError=closestapp-tank_radius;
//							
//							// since this is new best match: retrieve/override details
//							maxtracklength     = anmrdtrack.GetTrackLength();
//							mrdEntryVertex     = anmrdtrack.GetMrdEntryPoint();
//							mrdEntryTime       = anmrdtrack.GetStartTime();
//							mrdStopVertex      = anmrdtrack.GetStopVertex();
//							mrdTrackLength     = anmrdtrack.GetTrackLength();
//							mrdEnergyLoss      = anmrdtrack.GetEnergyLoss();
//							mrdEnergyLossError = anmrdtrack.GetEnergyLossError();
//							mrdPenetrated      = anmrdtrack.GetIsPenetrating();
//							mrdStopped         = anmrdtrack.GetIsStopped();
//							mrdSideExit        = anmrdtrack.GetIsSideExit();
//							
//							matchedsubevent=subev;
//							matchedtrack=tracki;
						} // else: best fit does not intercept tank
					} // else : no bonsai event case
				} // loop over tracks in this mrd subevent
			} // loop over mrd subevents (note: these are not separate triggers)
			
		}
		
		if(hadBonsaiEvent&&(matchedtrack>0)){
			if(true){ // TODO store MRD track even if it doesn't match bonsai <<<<<<<<<<<<<<<<<<
				// combine the information from bonsaiVertex and recalculate improved fit for 
				// mrd values such as angle, energy loss etc. 
				cMRDSubEvent* asubevent =(cMRDSubEvent*)mrdevents->At(matchedsubevent);
				auto anmrdtrack = asubevent->GetTracks()->at(matchedtrack);
				anmrdtrack.AddTrackPoint(bonsaiVertex,
					TVector3(bonsaiVertexError,bonsaiVertexError,bonsaiVertexError));
				anmrdtrack.DoTGraphErrorsFit();
				recoMuonAngle = anmrdtrack.GetTrackAngle();
				recoMuonAngleError = anmrdtrack.GetTrackAngleError();
				anmrdtrack.CheckIfStopping();
				// stopping status may have change due to adjustment of track angle affecting endpoint
				 // and fiducial requirements
				recoMuonStopped = anmrdtrack.GetIsStopped();
				anmrdtrack.CalculateEnergyLoss();
				// XXX where do we start using relativistic energy? need to add in rest mass
				recoMuonEnergy = anmrdtrack.GetEnergyLoss()+bonsaiEnergyLoss + pow(muonmass,2.);
				recoMuonEnergyError = anmrdtrack.GetEnergyLossError()+bonsaiEnergyLossError;
			}
			
			// whether or not the events were consistent, if we had both, fill the difference details
			bonsaiMrdTankVtxDiff = bonsaiVertex-mrdtankvertex;
		}
		
		
		// ==================================================================================================
		// load veto event // TODO
		// ==================================================================================================
		
		/* veto events are empty until dirt events are simulated */
		/*
		vetoevent->Clear();
		vetotree->GetEntry(wcsimTentry);
		*/
		hadVetoEvent  = false;
		vetoEventTime = 0.;
		
		// ==================================================================================================
		// combine events to calculate neutrino reconstructed information
		// ==================================================================================================
		
		if(recoVertex!=TVector3(0,0,0)){
			
			hadRecoEvent            = true;
			//recoVertex            = set above
			//recoVertexError       = set above
			//recoMuonEnergy        = set above
			//recoMuonEnergyError   = set above
			//recoMuonAngle         = set above
			//recoMuonAngleError    = set above
			recoNeutrinoEnergy      = CalculateNeutrinoEnergy(recoMuonEnergy, recoMuonAngle);
			double muangleerrorsign = (recoMuonAngle>0.) ? 1. : -1.;
			double recoEnumax       = CalculateNeutrinoEnergy(recoMuonEnergy+recoMuonEnergyError,
											recoMuonAngle-(muangleerrorsign*recoMuonAngleError));
			double recoEnumin       = CalculateNeutrinoEnergy(recoMuonEnergy-recoMuonEnergyError,
											recoMuonAngle+(muangleerrorsign*recoMuonAngleError));
			recoNeutrinoEnergyError = recoEnumax-recoEnumin;
			recoEventQ2             = CalculateEventQ2(recoMuonEnergy, recoNeutrinoEnergy, recoMuonAngle);
			double evq2max          = CalculateEventQ2(recoMuonEnergy+recoMuonEnergyError, 
										recoEnumin, recoMuonAngle-(muangleerrorsign*recoMuonAngleError));
			double evq2min          = CalculateEventQ2(recoMuonEnergy-recoMuonEnergyError, 
										recoEnumax, recoMuonAngle+(muangleerrorsign*recoMuonAngleError));
			recoEventQ2Error        = evq2max-evq2min;
		} else {
			// no consistent event with mrd and tank track
			hadRecoEvent            = false;
			recoMuonEnergy          = 0;
			recoMuonEnergyError     = 0;
			recoMuonAngle           = 0;
			recoMuonAngleError      = 0;
			recoEventQ2             = 0;
			recoEventQ2Error        = 0;
			recoNeutrinoEnergy      = 0;
			recoNeutrinoEnergyError = 0;
			recoVertex              = TVector3(0,0,0);
			recoVertexError         = TVector3(0,0,0);
		}
		
		// ==================================================================================================
		// should have set all the information now - fill the branches 
		// ==================================================================================================
		treeout->Fill();
		
	}  // end of loop over events
	
	// ======================================================================================================
	// ======================================================================================================
	
//	Double_t numbeamspills = totalpots/(4.0 * TMath::Power(10.,12.));
//	Double_t numbeamspillsperday = (24.*60.*60.*1000.)/133.3333;	// 24 hours in ms / 133.33 ms between spills
//	Double_t numdays = numbeamspills/numbeamspillsperday;
//	cout<<"Results based on "<<totalpots<<" POTs, or "<<numbeamspills<<" beam spills, or "<<numdays<<" days of data"<<endl;
//	cout<<"There were "<<numneutrinoeventsintank<<" neutrino interactions in the tank, of which "<<numCCQEneutrinoeventsintank<<" were true CCQE events."<<endl;
//	cout<<"Of those, "<<numCCQEneutrinoeventsinfidvol<<" were within the fiducial volume."<<endl;
//	cout<<"Of those in turn, "<<numCCQEneutrinoeventsinfidvolmrd<<" produced an accepted MRD muon"<<endl;
//	cout<<"There were "<<nummuontracksintank<<" muons in the tank, of which "
//		<<nummuontracksinfidvol<<" were from (CCQE?) events in the fiducial volume."<<endl;
	
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
	
	//cout<<"resetting bonsaitree branches"<<endl;
	if(bonsaitree) bonsaitree->ResetBranchAddresses();
	cout<<"closing bonsai file"<<endl;
	if(bonsaifile) bonsaifile->Close();
	
	//cout<<"resetting mrdfile branches"<<endl;
	if(mrdtree) mrdtree->ResetBranchAddresses();
	cout<<"closing mrdtrack file"<<endl;
	if(mrdfile) mrdfile->Close();
	
	//cout<<"resetting veto event branches"<<endl;
	if(vetotree) vetotree->ResetBranchAddresses();
	cout<<"closing veto event file"<<endl;
	if(vetofile) vetofile->Close();
	
	// write and close file of event information
	cout<<"writing and closing output file"<<endl;
	fileout->cd();
	treeout->SetEntries(bInTank->GetEntries());
	treeout->Write();
	fileout->Close();
	delete fileout;
	fileout=0;
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

double CalculateNeutrinoEnergy(double recoMuonEnergy, double recoMuonAngle){
	TDatabasePDG db;
	Double_t neutronmass = (db.GetParticle(2112)->Mass())*1000.; // converted to MeV
	Double_t protonmass = (db.GetParticle(2212)->Mass())*1000.; // converted to MeV
	Double_t muonmass = (db.GetParticle(13)->Mass())*1000.;      // converted to MeV
	Double_t O16bindingEnergy = 7.9762086875; // MeV (per nucleon), from http://tinyurl.com/y8m9s4z6
	Double_t boundneutronmass = neutronmass-O16bindingEnergy;
	
	// calculate neutrino energy: quasi-elastic, fermi-gas model
	double numerator = (2.*boundneutronmass*recoMuonEnergy) -
				(pow(boundneutronmass,2.) + pow(muonmass,2.) - pow(protonmass,2.));
	double denominator = 2.*( boundneutronmass - recoMuonEnergy + 
				sqrt( ( pow(recoMuonEnergy,2.) - (pow(muonmass,2.)*TMath::Cos(recoMuonAngle)) ) ));
	double recoNeutrinoEnergy = numerator / denominator;
	return recoNeutrinoEnergy;
}

double CalculateEventQ2(double recoMuonEnergy, double recoNeutrinoEnergy, double recoMuonAngle){
	TDatabasePDG db;
	Double_t neutronmass = (db.GetParticle(2112)->Mass())*1000.; // converted to MeV
	Double_t protonmass = (db.GetParticle(2212)->Mass())*1000.; // converted to MeV
	Double_t muonmass = (db.GetParticle(13)->Mass())*1000.;      // converted to MeV
	Double_t O16bindingEnergy = 7.9762086875; // MeV (per nucleon), from http://tinyurl.com/y8m9s4z6
	Double_t boundneutronmass = neutronmass-O16bindingEnergy;
	
	double part1 = recoMuonEnergy - (sqrt(pow(recoMuonEnergy,2.)-pow(muonmass,2.))*TMath::Cos(recoMuonAngle));
	double eventq2 = -pow(muonmass,2.) + 2.*recoNeutrinoEnergy*part1;
	return eventq2;
}
