/* vim:set noexpandtab tabstop=4 wrap */

#ifndef VERBOSE
#define VERBOSE
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

// we need to #include all the WCSim headers.
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimPmtInfo.hh"
#include "WCSimLAPPDInfo.hh"
#include "WCSimEnumerations.hh"
#include "WCSimRootLinkDef.hh"

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
#include "makepmtmaps_standalone.cxx"	// definition of this function

void ColourPlotStyle();

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
########## MRD scintillator layer 3  (V) at z=302.865 ##########
########## MRD scintillator layer 4  (H) at z=314.975 ##########
########## MRD scintillator layer 5  (V) at z=327.085 ##########
########## MRD scintillator layer 6  (H) at z=339.195 ##########
########## MRD scintillator layer 7  (V) at z=351.305 ##########
########## MRD scintillator layer 8  (H) at z=363.415 ##########
########## MRD scintillator layer 9  (V) at z=375.525 ##########
########## MRD scintillator layer 10 (H) at z=387.635 ########## */

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
const Float_t fidcutz=tank_start+tank_radius;	// fidcuial volume is before the tank centre.

// needed for drawing tank 2D map histograms
std::map<int, std::pair<int,int> > topcappositionmap;
std::map<int, std::pair<int,int> > bottomcappositionmap;
std::map<int, std::pair<int,int> > wallpositionmap;

const char* dirtpath="/pnfs/annie/persistent/users/moflaher/g4dirt";
const char* geniepath="/pnfs/annie/persistent/users/rhatcher/genie";
const char* wcsimpath="/pnfs/annie/persistent/users/moflaher/wcsim";
const char* wcsimlibrarypath="/annie/app/users/moflaher/wcsim/wcsim/libWCSimRoot.so";
const char* outpath="/annie/app/users/moflaher/wcsim/root_work";

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
	
	// Declare paths, files, trees...
	TFile* wcsimfile=0;
	TString wcsimfilepath;
	TTree* wcsimT=0;
	Int_t numwcsimentries=0;
	
	TFile* dirtfile=0;
	TTree* tankflux=0;
	TTree* tankmeta=0;
	
	TFile* geniefile=0;
	TString geniefilepath;
	TTree* gtree=0;
	Int_t numgenietentries=0;
	
	// TChain for dirt files - this will be the main driver of the loop - all it's events will be processed.
	TChain* c =  new TChain("tankflux");
	TString chainpattern = TString::Format("%s/annie_tank_flux.*.root",dirtpath);
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
	//Int_t* pdgbranchval=0;	// retrieve this info from wcsimT
	//TBranch* pdgBranch=0;
	
	// tankmeta
	TBranch* geniefilenamebranch=0;
	Char_t geniefilename[100];
	TBranch* potsbranch=0;
	Double_t pots;
	Double_t totalpots=0;					// count of POTs in all processed files
	
	Double_t numneutrinoeventsintank=0.;
	Double_t numQEneutrinoeventsintank=0.;
	Double_t numQEneutrinoeventsinfidvol=0.;
	Double_t numQEneutrinoeventsinfidvolmrd=0.;
	Double_t nummuontracksintank=0.;
	Double_t nummuontracksinfidvol=0.;
	Double_t nummuontracksintankpassedcut=0.;
	Double_t nummuontracksinmrd=0.;
	Double_t nummuontracksinfidvolmrd=0.;
	
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
	//gROOT->cd();
	TFile* histofileout = new TFile("TruthHistos.root","RECREATE");
	
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
	
	// information from WCSim: 
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
	TH3D* chargemap_nopions = new TH2D("chargemap_nopions", "Charge Distribution for CC0pi events", pmtsperring+2,-1,pmtsperring+1,numpmtrings+2,-1,numpmtrings+1, 100, 0., 1400.);
	
	// Just to test the inside/outside cherenkov cone algorithm
	TH2D* chargemap_incone = new TH2D("chargemap_incone", "Charge Distribution Inside Cherenkov Cone", pmtsperring+2,-1,pmtsperring+1,numpmtrings+2,-1,numpmtrings+1);
	TH2D* chargemap_outcone = new TH2D("chargemap_outcone", "Charge Distribution Outside Cherenkov Cone", pmtsperring+2,-1,pmtsperring+1,numpmtrings+2,-1,numpmtrings+1);
	
	// create the file for outputting true vertices and digits for tank reconstruction efforts
	TFile* flateventfileout = new TFile("trueQEvertexinfo.root", "RECREATE");
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
	TTree* vertextreenocuts = new TTree("vertextreenocuts","All True Tank QE Events");
	TBranch* MuonStartBranch = vertextreenocuts->Branch("MuonStartVertex",&filemuonstartvertex);
	TBranch* MuonStopBranch = vertextreenocuts->Branch("MuonStopVertex", &filemuonstopvertex);
	TBranch* MuonDirectionBranch = vertextreenocuts->Branch("MuonDirection", &filemuondirectionvector);
	TBranch* DigitVertexBranch = vertextreenocuts->Branch("DigitVertices", "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >", &filedigitverticesp);
//	TBranch* DigitVertexBranch = vertextreenocuts->Branch("DigitVertices", "std::vector<TLorentzVector>", &filedigitverticesp);
	TBranch* DigitChargeBranch = vertextreenocuts->Branch("DigitCharges", &filedigitQsp);
	if(MuonStartBranch==0||MuonStopBranch==0||MuonDirectionBranch==0||DigitVertexBranch==0||DigitChargeBranch==0){ 
		cerr<<"branches are zombies argh!"<<endl; 
		cout<<"MuonStartBranch="<<MuonStartBranch<<endl
			<<"MuonStopBranch="<<MuonStopBranch<<endl
			<<"MuonDirectionBranch="<<MuonDirectionBranch<<endl
			<<"DigitVertexBranch="<<DigitVertexBranch<<endl
			<<"DigitChargeBranch="<<DigitChargeBranch<<endl;
		assert(false&&"branches are zombies argh!");
	}
	
	TTree* vertextreefiducialcut = new TTree("vertextreefiducialcut","True Tank QE Events in Fiducial Volume");
	TBranch* MuonStartBranchFid = vertextreefiducialcut->Branch("MuonStartVertex",&filemuonstartvertex);
	TBranch* MuonStopBranchFid = vertextreefiducialcut->Branch("MuonStopVertex", &filemuonstopvertex);
	TBranch* MuonDirectionBranchFid = vertextreefiducialcut->Branch("MuonDirection", &filemuondirectionvector);
	TBranch* DigitVertexBranchFid = vertextreefiducialcut->Branch("DigitVertices", "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >", &filedigitverticesp);
//	TBranch* DigitVertexBranchFid = vertextreefiducialcut->Branch("DigitVertices", "std::vector<TLorentzVector>", &filedigitverticesp);
	TBranch* DigitChargeBranchFid = vertextreefiducialcut->Branch("DigitCharges", &filedigitQsp);
	
	TTree* vertextreefiducialmrd = new TTree("vertextreefiducialmrd","True Tank QE Events in Fiducial Volume With Muon in MRD");
	TBranch* MuonStartBranchFidMRD = vertextreefiducialmrd->Branch("MuonStartVertex",&filemuonstartvertex);
	TBranch* MuonStopBranchFidMRD = vertextreefiducialmrd->Branch("MuonStopVertex", &filemuonstopvertex);
	TBranch* MuonDirectionBranchFidMRD = vertextreefiducialmrd->Branch("MuonDirection", &filemuondirectionvector);
	TBranch* DigitVertexBranchFidMRD = vertextreefiducialmrd->Branch("DigitVertices", "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >", &filedigitverticesp);
//	TBranch* DigitVertexBranchFidMRD = vertextreefiducialmrd->Branch("DigitVertices", "std::vector<TLorentzVector>", &filedigitverticesp);
	TBranch* DigitChargeBranchFidMRD = vertextreefiducialmrd->Branch("DigitCharges", &filedigitQsp);
	
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
	for(Int_t inputEntry=0; inputEntry<numents; inputEntry++){
		/* 	1. Load next g4dirt entry */ 
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
			std::string nextdirtfilename(dirtfile->GetName());
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
				inputEntry += thistreesentries;	// skip the loop iterator forward by all the entries in this file
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
			cout<<"matching regex for filename "<<nextdirtfilename<<endl;
			std::regex_match (nextdirtfilename, submatches, theexpression);
			std::string submatch = (std::string)submatches[0];	// match 0 is 'whole match' or smthg
			if(submatch==""){ cout<<"unrecognised input file pattern: "<<nextdirtfilename<<endl; return; }
			submatch = (std::string)submatches[1];
			cout<<"extracted submatch is "<<submatch<<endl;
			int filenum = atoi(submatch.c_str());
			
			// use filenum to open the corresponding wcsim file
			wcsimfilepath = TString::Format("%s/wcsim_0.%d.root",wcsimpath,filenum);
			cout<<"corresponding wcsim file is "<<wcsimfilepath<<endl;
			if(wcsimfile) wcsimfile->Close(); wcsimfile=0;
			wcsimfile = TFile::Open(wcsimfilepath);
			if(!wcsimfile){
				cout<<"this wcsimfile doesn't exist!"<<endl; 
				inputEntry += thistreesentries;	// skip the loop iterator forward by all the entries in this file
				continue; 
			}
			// load the geometry tree and grab the geometry if we haven't already
			if(geo==0){
				TTree* geotree = (TTree*)wcsimfile->Get("wcsimGeoT");
				if(geotree==0){ cout<<"NO GEOMETRY IN FIRST FILE?"<<endl; assert(false); }
				geotree->SetBranchAddress("wcsimrootgeom", &geo);
				if (geotree->GetEntries() == 0) { cout<<"geotree has no entries!"<<endl; exit(9); }
				geotree->GetEntry(0);
				MakePMTmap(WCSimRootGeom* geo, topcappositionmap, bottomcappositionmap, wallpositionmap);
			}
			// load the next set of wcsim event info
			wcsimT = (TTree*)wcsimfile->Get("wcsimT");
			if(!wcsimT){cout<<"wcsimT doesn't exist!"<<endl; break; }
			numwcsimentries = wcsimT->GetEntries();
			cout<<"wcsimT has "<<numwcsimentries<<" entries in this file"<<endl;
			if(numwcsimentries<1){cout<<"wcsimT has no entries!"<<endl; break; }
			
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
				cout<<"ERROR! TOTAL POTS CAME DOWN FROM "<<lasttotpots<<" TO "<<totalpots<<endl; return; 
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
#ifdef VERBOSE
		cout<<"processing inputEntry "<<inputEntry<<", localEntry "<<localEntry
		    <<"/"<<thistreesentries<<" in tree "<<treeNumber<<endl;
#endif
		
		/* 2. Check if genie primary, and volume is in tank - if not, continue */
		// cout<<"getting wcsim info"<<endl;
		nTankBranch->GetEntry(localEntry);
		vertexmaterialbranch->GetEntry(localEntry);
		if(strcmp(vertexmaterial,"TankWater")!=0){ continue; }	// neutrino intx wasn't in tank water
		if(nuprimarybranchval){delete[] nuprimarybranchval;}
		nuprimarybranchval = new Int_t[ntankbranchval];
		nuprimaryBranch->SetAddress(nuprimarybranchval);
		nuprimaryBranch->GetEntry(localEntry);
		
		Bool_t primariesinthisentry=false;
		for(int i=0;i<ntankbranchval;i++){
			if(nuprimarybranchval[i]==1){ primariesinthisentry=true; break; }
		}
		if(!primariesinthisentry){ continue; }	// dirt recorded particles weren't genie primaries
		
		/* 3. If so, load genie entry. */
#ifdef VERBOSE
		cout<<"getting genie info"<<endl;
#endif
		if(localEntry>numgenietentries){ cout<<"can't load localEntry "<<localEntry
								 <<"from "<<geniefilepath<<" gtree: not enough entries!"<<endl; continue; }
		genieentrybranch->GetEntry(localEntry);
		genierecordBranch->GetEntry(genieentry);
		genie::EventRecord* gevtRec = genierecordval->event;
		genie::Interaction* genieint = gevtRec->Summary();
		
		// process information:
		TString procinfostring = genieint->ProcInfo().AsString();
		TString scatteringtypestring = genieint->ProcInfo().ScatteringTypeAsString();
		TString interactiontypestring = genieint->ProcInfo().InteractionTypeAsString();
		Bool_t isQE = genieint->ProcInfo().IsQuasiElastic();
		Bool_t isWeakCC = genieint->ProcInfo().IsWeakCC();
		Int_t neutinteractioncode = genie::utils::ghep::NeutReactionCode(gevtRec);
		Int_t nuanceinteractioncode  = genie::utils::ghep::NuanceReactionCode(gevtRec);
		TLorentzVector* genieVtx = gevtRec->Vertex();
		Double_t genie_x = genieVtx->X() * 100.;         // same info as nuvtx in g4dirt file
		Double_t genie_y = genieVtx->Y() * 100.;         // GENIE uses meters
		Double_t genie_z = genieVtx->Z() * 100.;         // GENIE uses meters
		Double_t genie_t = genieVtx->T() * 1000000000;   // GENIE uses seconds for time
		
		// neutrino information:
		Double_t probeenergy = genieint->InitState().ProbeE(genie::kRfLab);	// GeV
		Int_t probepdg = genieint->InitState().Probe()->PdgCode();
		TString probepartname = genieint->InitState().Probe()->GetName();
		TLorentzVector* probemomentum = gevtRec->Probe()->P4();
		TVector3 probethreemomentum = probemomentum->Vect();
		TVector3 probemomentumdir = probethreemomentum.Unit();
		Double_t probeanglex = TMath::ATan(probethreemomentum.X()/probethreemomentum.Z());
		Double_t probeangley = TMath::ATan(probethreemomentum.Y()/probethreemomentum.Z());
		Double_t probeangle = TMath::Max(probeanglex,probeangley);
		// n.b.  genieint->InitState().Probe != gevtRec->Probe()
		
		// target nucleon:
		genie::GHepParticle* targetnucleon = gevtRec->HitNucleon();
		int targetnucleonpdg = genieint->InitState().Tgt().HitNucPdg();
		TString targetnucleonname;
		if ( genie::pdg::IsNeutronOrProton(targetnucleonpdg) ) {
			TParticlePDG * p = genie::PDGLibrary::Instance()->Find(targetnucleonpdg);
			targetnucleonname = p->GetName();
		} else {
			targetnucleonname = targetnucleonpdg;
		}
		TVector3 targetnucleonthreemomentum=TVector3(0.,0.,0.);
		Double_t targetnucleonenergy=0.;
		if(targetnucleon){
			TLorentzVector* targetnucleonmomentum = targetnucleon->P4();
			targetnucleonthreemomentum = targetnucleonmomentum->Vect();
			targetnucleonenergy = targetnucleonmomentum->Energy(); //GeV
		}
		
		// target nucleus:
		Int_t targetnucleuspdg = genieint->InitState().Tgt().Pdg();
		TParticlePDG* targetnucleus = genie::PDGLibrary::Instance()->Find( targetnucleuspdg );
		TString targetnucleusname = "unknown";
		if(targetnucleus){ targetnucleusname = targetnucleus->GetName(); }
		Int_t targetnucleusZ = genieint->InitState().Tgt().Z();
		Int_t targetnucleusA = genieint->InitState().Tgt().A();
		
		// remnant nucleus:
		int remnucpos = gevtRec->RemnantNucleusPosition(); 
		TString remnantnucleusname="n/a";
		Double_t remnantnucleusenergy=-1.;
		if(remnucpos>-1){
			remnantnucleusname = gevtRec->Particle(remnucpos)->Name();
			remnantnucleusenergy = gevtRec->Particle(remnucpos)->Energy(); //GeV
		}
		
		// final state lepton:
		int fsleppos = gevtRec->FinalStatePrimaryLeptonPosition();
		TString fsleptonname="n/a";
		Double_t fsleptonenergy=-1.;
		if(fsleppos>-1){
			fsleptonname = gevtRec->Particle(fsleppos)->Name();
			fsleptonenergy = gevtRec->Particle(fsleppos)->Energy();
		}
		
		// other remnants: TODO: this information is NOT being correctly read in
		Int_t numfsprotons = genieint->ExclTag().NProtons();
		Int_t numfsneutrons = genieint->ExclTag().NNeutrons();
		Int_t numfspi0 = genieint->ExclTag().NPi0();
		Int_t numfspiplus = genieint->ExclTag().NPiPlus();
		Int_t numfspiminus = genieint->ExclTag().NPiMinus();
		
		// kinematic information
		Double_t M  = genie::constants::kNucleonMass; 
		// Calculate kinematic variables "as an experimentalist would measure them; 
		// neglecting fermi momentum and off-shellness of bound nucleons"
		TLorentzVector& k1 = *(gevtRec->Probe()->P4());
		TLorentzVector& k2 = *(gevtRec->FinalStatePrimaryLepton()->P4());
		Double_t costhfsl = TMath::Cos( k2.Vect().Angle(k1.Vect()) );
		Double_t fslanglegenie = k2.Vect().Angle(k1.Vect());
		TLorentzVector q  = k1-k2;                              // q=k1-k2, 4-p transfer
		//Double_t neutrinoq2 = genieint->Kine().Q2();          // not set in our GENIE files!
		Double_t Q2 = -1 * q.M2();                              // momemtum transfer
		Double_t Ev  = (targetnucleon) ? q.Energy()       : -1; // v (E transfer to the nucleus)
		Double_t x  = (targetnucleon) ? 0.5*Q2/(M*Ev)     : -1; // Bjorken x
		Double_t y  = (targetnucleon) ? Ev/k1.Energy()    : -1; // Inelasticity, y = q*P1/k1*P1
		Double_t W2 = (targetnucleon) ? M*M + 2*M*Ev - Q2 : -1; // Hadronic Invariant mass ^ 2
		
		if(printneutrinoevent){
			cout<<"This was a "<< procinfostring <<" interaction of a "<<probeenergy<<"GeV " 
				<< probepartname << " on a "; 
			
			if( targetnucleonpdg==2212 || targetnucleonpdg==2122 ){ cout<<targetnucleonname<<" in a "; }
			else { cout<<"PDG-Code " << targetnucleonpdg<<" in a "; }
			
			if( targetnucleusname!="unknown"){ cout<<targetnucleusname<<" nucleus, "; }
			else { cout<<"Z=["<<targetnucleusZ<<","<<targetnucleusA<<"] nucleus, "; }
			
			if(remnucpos>-1){ cout<<"producing a "<<remnantnucleusenergy<<"GeV "<<remnantnucleusname; }
			else { cout<<"with no remnant nucleus"; }  // DIS on 16O produces no remnant nucleus?!
			
			if(fsleppos>-1){ cout<<" and a "<<fsleptonenergy<<"GeV "<<fsleptonname<<endl; }
			else{ cout<<" and no final state leptons"<<endl; }
			
			cout<<endl<<"Q^2 was "<<Q2<<"(GeV/c)^2, with final state lepton ejected at Cos(θ)="<<costhfsl<<endl;
			cout<<"Additional final state particles included "<<endl;
			cout<< "   N(p) = "       << numfsprotons
				<< "   N(n) = "       << numfsneutrons
				<< endl
				<< "   N(pi^0) = "    << numfspi0
				<< "   N(pi^+) = "    << numfspiplus
				<< "   N(pi^-) = "    << numfspiminus
				<<endl;
		}
		
		numneutrinoeventsintank++;
		/* 4. Check if interaction is QE - if not, continue */
		if (!(isQE && isWeakCC)){ continue; }
#ifdef VERBOSE
		cout<<"filling incident neutrino histograms"<<endl;
#endif
		incidentneutrinoenergiesall->Fill(probeenergy);
		//incidentneutrinoanglesall->Fill(probeangle);
		fslanglesall->Fill(fslanglegenie);
		fslenergiesall->Fill(fsleptonenergy);
		eventq2all->Fill(Q2);
		neutrinovertex->Fill(genie_x, genie_y, genie_z);
		neutrinovertexQE->Fill(genie_x, genie_y, genie_z);
		numQEneutrinoeventsintank++;
		
		/* 4.5 Do fiducial volume cut: */
		if( (TMath::Sqrt(TMath::Power(genie_x,2)+TMath::Power(genie_z,2)) < fidcutradius) && 
			(TMath::Abs(genie_y-tank_yoffset) < fidcuty) && 
			(TMath::Abs(genie_z-tank_start-tank_radius) < fidcutz) ){
			incidentneutrinoenergiesallfidcut->Fill(probeenergy);
			fslanglesallfidcut->Fill(fslanglegenie);
			fslenergiesallfidcut->Fill(fsleptonenergy);
			eventq2allfidcut->Fill(Q2);
			numQEneutrinoeventsinfidvol++;
		}
		
		/*5. primary neutrino true QE vertex in the tank: load wcsim detector response. */
		// read only first subtrigger; delayed decay detector response is not interesting for primary FSL tracks
#ifdef VERBOSE
		cout<<"getting wcsim entry"<<endl;
#endif
		if(localEntry>numwcsimentries){ cout<<"can't load localEntry "<<localEntry
				<<"from "<<wcsimfilepath<<" wcsimT - not enough entries!"<<endl; continue; }
		wcsimT->GetEntry(localEntry);
		atrigt = b->GetTrigger(0);
		atrigm = m->GetTrigger(0);
		atrigv = v->GetTrigger(0);
		
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
		
		/* yet TODO - scan for neutron captures here? */
		
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
			// for now we use truth information
			// is it a primary?
			Int_t primaryparentpdg = nextrack->GetParenttype();
			if(primaryparentpdg!=0) continue;
			
			// is it a (anti)muon?
			Int_t primarypdg = nextrack->GetIpnu();
#ifdef VERBOSE
			cout<<"primarypdg is "<<primarypdg<<endl;
#endif
			if(TMath::Abs(primarypdg)!=13) continue; 		// not a muon
			
			// does it start in the tank?
			Int_t primarystartvol = nextrack->GetStartvol();
			// temporary override as these weren't correctly set in wcsim #################
			if(nextrack->GetStart(2)<tank_start){
				primarystartvol = 20;						// start depth is facc or hall
			} else if(nextrack->GetStart(2)>(tank_start+(2.*tank_radius))){
				primarystartvol = 30;						// start depth is mrd or hall
			} else {
				primarystartvol = 10;						// start depth is tank
			}
			
#ifdef VERBOSE
			cout<<"primarystartvol is "<<primarystartvol<<endl;
#endif
			if(primarystartvol!=10) continue;				// start volume is not the tank
			
			// does it stop in the mrd, or pass completely through the MRD? 
			Int_t primarystopvol = nextrack->GetStopvol();
			// temporary override as these weren't correctly set in wcsim #################
			if(nextrack->GetStop(2)<tank_start){
				primarystopvol = 20;						// start depth is facc or hall
			} else if(nextrack->GetStop(2)>(tank_start+(2.*tank_radius))){
				primarystopvol = 30;						// start depth is mrd or hall
			} else {
				primarystopvol = 10;						// start depth is tank
			}
#ifdef VERBOSE
			cout<<"primarystopvol is "<<primarystopvol<<endl;
#endif
			
			TLorentzVector primarystartvertex(nextrack->GetStart(0),nextrack->GetStart(1),nextrack->GetStart(2), genie_t);
			TLorentzVector primarystopvertex(nextrack->GetStop(0), nextrack->GetStop(1), nextrack->GetStop(2), nextrack->GetTime());
			
			Float_t oppx = primarystopvertex.X() - primarystartvertex.X();
			Float_t adj = primarystopvertex.Z() - primarystartvertex.Z();
			Float_t avgtrackanglex = TMath::ATan(oppx/adj);
			Float_t oppy = primarystopvertex.Y() - primarystartvertex.Y();
			Float_t avgtrackangley = TMath::ATan(oppy/adj);
			
			TVector3 differencevector  = (primarystopvertex.Vect()-primarystartvertex.Vect());
			fsltruetracklength->Fill(differencevector.Mag());
			// to calculate track length _in water_ find distance from start vertex to the point
			// where it intercepts the tank. if this length > total track length, use total track length
			// otherwise use this length. 
			
			// first find out the z value where the tank would leave the radius of the tank
#ifdef MUTRACKDEBUG
			cout<<"z0 = "<<genie_z-tank_start-tank_radius<<", x0 = "<<genie_x<<endl;
#endif
			Double_t xatziszero = (genie_x - (genie_z-tank_start-tank_radius)*TMath::Tan(avgtrackanglex));
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
			Double_t tankendpointx = TMath::Sqrt(TMath::Power(tank_radius,2)-TMath::Power(tankendpointz,2));
			// correct for tank z offset (do after tankendpointx, before tankendpointy)
			tankendpointz += tank_start+tank_radius;
			// now check if the particle would have exited through one of the caps before reaching this radius
			Double_t tankendpointy = genie_y + (tankendpointz-genie_z)*TMath::Tan(avgtrackangley);
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
				tankendpointz = genie_z + (tankendpointy-genie_y)/TMath::Tan(avgtrackangley);
				tankendpointx = genie_x + (tankendpointz-genie_z)*TMath::Tan(avgtrackanglex);
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
			Double_t mutracklengthintank = TMath::Sqrt(
				TMath::Power((tankendpointx-genie_x),2)+
				TMath::Power((tankendpointy-genie_y),2)+
				TMath::Power((tankendpointz-genie_z),2) );
#ifdef MUTRACKDEBUG
			cout<<"muon tank track length: ("<<(tankendpointx-genie_x)<<", "<<(tankendpointy-genie_y)<<", "
				<<(tankendpointz-genie_z)<<") = "<<mutracklengthintank<<"cm total"<<endl;
			cout<<"muon tank exit point: ("<<tankendpointx<<", "<<tankendpointy<<", "<<tankendpointz<<") ";
			cout<<"muon start point : ("<<genie_x<<", "<<genie_y<<", "<<genie_z<<")"<<endl;
#endif
			if(mutracklengthintank > maxtanktracklength){
				cout<<"Track length is impossibly long!"<<endl;
				assert(false);
			}
			if(TMath::IsNaN(mutracklengthintank)){
				cout<<"NaN RESULT FROM MU TRACK LENGTH IN TANK?!"<<endl;
				assert(false);
			}
			nummuontracksintank++;
			fsltruetracklengthintank->Fill(mutracklengthintank);
			if(mutracklengthintank>50){ nummuontracksintankpassedcut++; }
			
			///////////////////////////////////////////////////
			
			// Now how about calculating the energy detected by tank PMTs from this track?
			// search for all digits, select those whose photons have this muon as their parent,
			// and make a list of those digits and their Q. Then, sum the Q's. 
			//cout<<endl<<"calculating muon energy deposition, mu track id "<<nextrack->GetId()<<endl;
			std::vector<int> numphotsinthedigits;
			std::vector<int> numphotsfromthemuoninthedigits;
			std::vector<double> chargesinthedigits;
			Int_t numdigitsthisevent = atrigt->GetCherenkovDigiHits()->GetEntries();
			//cout<<"this event has "<<numdigitsthisevent<<" digits"<<endl;

			filedigitvertices.clear();
			filedigitQs.clear();
			Double_t chargeinsidecherenkovcone=0;
			Double_t chargeoutsidecherenkovcone=0;
			for(Int_t i=0; i<numdigitsthisevent; i++){
				// retrieve the digit information
				// ==============================
				// 1. note all in a flat tree reconstruction dev
				// 2. add up charges in digits whose parent photons are the primary muon, for muon tank Edep
				WCSimRootCherenkovDigiHit* thedigihit = 
					(WCSimRootCherenkovDigiHit*)atrigt->GetCherenkovDigiHits()->At(i);
				double digitsq = thedigihit->GetQ();
				/////////////////
				// digit times not working; for now save one (or more) of the true hit times
				/////////////////
				//double digitst = thedigihit->GetT(); 
				double digitst = 0.;
				std::vector<int> truephotonindices = thedigihit->GetPhotonIds();
				std::vector<double> thetruetimes;
				for(int truephoton=0; truephoton</*truephotonindices.size()*/1; truephoton++){
					int thephotonsid = truephotonindices.at(truephoton);
					WCSimRootCherenkovHitTime *thehittimeobject = 
						(WCSimRootCherenkovHitTime*)atrigt->GetCherenkovHitTimes()->At(thephotonsid);
					double ahittime = thehittimeobject->GetTruetime();
					thetruetimes.push_back(ahittime);
				}
				std::valarray<double> truetimes2(&thetruetimes.front(), thetruetimes.size());
				digitst = (truetimes2.sum())/(truetimes2.size());
				/////////////////
				Int_t digitstubeid = thedigihit->GetTubeId();
				WCSimRootPMT pmt = geo->GetPMT(digitstubeid);
				double digitsx = pmt.GetPosition(0);
				double digitsy = pmt.GetPosition(1);
				double digitsz = pmt.GetPosition(2);
				//TLorentzVector* adigitvector = new TLorentzVector(digitsx,digitsy,digitsz,digitst);
				ROOT::Math::XYZTVector adigitvector = ROOT::Math::XYZTVector(digitsx,digitsy,digitsz,digitst);
				filedigitvertices.push_back(adigitvector);
//				TLorentzVector adigitvector = TLorentzVector(digitsx,digitsy,digitsz,digitst);
//				filedigitvertices.push_back(adigitvector);
				filedigitQs.push_back(digitsq);
//TODO: Scan through tubes/digits contribution and look for light from this muon trackID
// Can we make a plot of just the light from this muon, and the light from everything else?
// Moreover, add up all Q from digits on tubes inside & outside the muon ring
// for this we need to check if the digit lies within the cherenkov angle from the muon start vertex
// ***********************************************************************************************
				TVector3 digitvector(digitsx-genie_x, digitsy-genie_y, digitsz-genie_z);
				Double_t dotproduct = (digitvector.Unit()).Dot(differencevector.Unit());
				Double_t digitmuonangle = TMath::ACos(dotproduct);
				if(digitmuonangle<((42.0*TMath::Pi())/180.0)){
					chargeinsidecherenkovcone+=digitsq;
					FillTankMapHist(geo, digitstubeid, chargemap_incone); 
					
					//FillTankMapHist(WCSimRootGeom* geo, int tubeID, TH1* histowall=0, TH1* histotop=0, TH1* histobottom=0, double weight=1)
					
				} else { 
					chargeoutsidecherenkovcone+=digitsq;
					chargemap_outcone->Fill(
					
					FillTankMapHist(WCSimRootGeom* geo, int tubeID, TH1* histowall, TH1* histotop, TH1* histobottom, double weight=1)
				}
				
				// scan through the parents ID's for the photons contributing to this digit
				// and see if any of them are this muon
				std::vector<int> truephotonindices = thedigihit->GetPhotonIds();
				int numphotonsfromthismuon=0;
				cout<<"   digit "<<i<<" has "<<truephotonindices.size()<<" true photons"<<endl;
				for(int truephoton=0; truephoton<truephotonindices.size(); truephoton++){
					int thephotonsid = truephotonindices.at(truephoton);
					WCSimRootCherenkovHitTime *thehittimeobject = 
						(WCSimRootCherenkovHitTime*)atrigt->GetCherenkovHitTimes()->At(thephotonsid);
					Int_t thephotonsparenttrackid = thehittimeobject->GetParentID();
					cout<<"      HitTimeIndex "<<thephotonsid<<", HitTimeObject "<<thehittimeobject
						<<", HitTimeParent "<<thephotonsparenttrackid<<endl;
					if(thephotonsparenttrackid==nextrack->GetId()) {
						numphotonsfromthismuon++;
						cout<<"################ FOUND A DIGIT ############"<<endl;
						cout<<"the digits Q is "<<thedigihit->GetQ()<<endl;
						return;
					}
				}
				if(numphotonsfromthismuon!=0){	//this digit had contribution from the muon!
					numphotsinthedigits.push_back(truephotonindices.size());
					numphotsfromthemuoninthedigits.push_back(numphotonsfromthismuon);
					chargesinthedigits.push_back(thedigihit->GetQ());
				}
			}
			// Save vertex and digit information to a flat tree for reconstruction development
			filemuonstartvertex = primarystartvertex;
			filemuonstopvertex = primarystopvertex;
			filemuondirectionvector = differencevector.Unit();
		
			MuonStartBranch->Fill();
			MuonStopBranch->Fill();
			MuonDirectionBranch->Fill();
			DigitVertexBranch->Fill();
			DigitChargeBranch->Fill();
			if( (TMath::Sqrt(TMath::Power(genie_x,2)+TMath::Power(genie_z,2)) < fidcutradius) && 
				(TMath::Abs(genie_y-tank_yoffset) < fidcuty) && 
				(TMath::Abs(genie_z-tank_start-tank_radius) < fidcutz) ){
					MuonStartBranchFid->Fill();
					MuonStopBranchFid->Fill();
					MuonDirectionBranchFid->Fill();
					DigitVertexBranchFid->Fill();
					DigitChargeBranchFid->Fill();
					nummuontracksinfidvol++;
			}
			
			// Now sum up the energy by adding the charges in the digits, scaled by photonic contribution
			Double_t energydepositionofmuon=0.;
			for(int idigit=0; idigit<chargesinthedigits.size(); idigit++){
				cout<<"muon digit with charge "<<chargesinthedigits.at(idigit)
					<<" scaled by contribution factor "
					<<(numphotsinthedigits.at(idigit)/numphotsfromthemuoninthedigits.at(idigit))<<endl;
				energydepositionofmuon += chargesinthedigits.at(idigit) * 
					(numphotsinthedigits.at(idigit)/numphotsfromthemuoninthedigits.at(idigit));
			}
			muedepositionswcsim->Fill(energydepositionofmuon);
			muonenergydepositions.push_back(energydepositionofmuon);
			if(energydepositionofmuon!=0) return;
			
			////////////////////////////////////////////////
			
			// continue if stopping volume is in either the tank or veto, or track is backward going.
			if( primarystopvol==10 || primarystopvol==20 || primarystopvertex.Z() < primarystartvertex.Z() )
				continue;
			
			/* We next need to check it intercepts the MRD. The mrd is as wide as the tank. 
			We can ensure a track enters the mrd by projecting the track forward from the start vertex, 
			at the angle between start and stop vertices, and requiring:
			1) at z=MRD_start, x is between +/-(MRD_width/2);
			2) z endpoint is > MRD_start
			To require at least 2-layers of penetration, as above but with z=MRD_layer2 */
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
			
			if((TMath::Abs(xatmrd)>MRD_width)||(TMath::Abs(yatmrd)>MRD_height)) 
				continue;	// track does not meet MRD penetration requirement
			
			// if we got to here, we have a muon that stops in, or penetrates at least 2 MRD layers!
			// retrieve any remaining information and record the event
			Float_t mu_rest_mass_E = 105.658;
			Float_t primaryenergy = (nextrack->GetE()-mu_rest_mass_E)/1000.;	// starting energy (GeV) (p^2+m^2)
			Float_t primarymomentummag = nextrack->GetP();		// starting momentum
			TVector3 primarymomentumdir(nextrack->GetPdir(0),nextrack->GetPdir(1),nextrack->GetPdir(2));
			Float_t starttrackanglex = TMath::ATan(primarymomentumdir.X()/primarymomentumdir.Z());
			Float_t starttrackangley = TMath::ATan(primarymomentumdir.Y()/primarymomentumdir.Y());
			Float_t starttrackangle = TMath::Max(starttrackanglex,starttrackangley);
			//Double_t scatteringangle = k1.Angle(primarymomentumdir);
			Double_t scatteringangle = k1.Angle(differencevector);
			// scatteringangle for primaries by definition should be the same as fslanglegenie
			Double_t neutrinoenergyguess = probeenergy;	//TODO: how do we estimate this?
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
		}
		
		if(scatteringanglesvector.size()>1) {
			cerr<<"Mutliple accepted final state leptons!"<<endl;
			for(int i=0; i<scatteringanglesvector.size(); i++){
				cout<<"wcsim lepton angle "<<i<<" = "<<scatteringanglesvector.at(i)<<endl;
				cout<<"wcsim lepton energy "<<i<<" = "<<primaryenergiesvector.at(i)<<endl;
				cout<<"wcsim lepton trackid "<<i<<" = "<<acceptedtrackids.at(i)<<endl;
			}
			cout<<"compare to:"<<endl;
			cout<<"genie lepton angle = "<<fslanglegenie<<endl;
			cout<<"genie lepton energy = "<<fsleptonenergy<<endl;
			
			std::vector<Int_t>::iterator minit = std::min_element(acceptedtrackids.begin(), acceptedtrackids.end());		
			Int_t indextouse = std::distance(acceptedtrackids.begin(),minit);
			incidentneutrinoenergiesacceptedwcsim->Fill(neutrinoenergiesvector.at(indextouse));
			fslanglesacceptedwcsim->Fill(scatteringanglesvector.at(indextouse));
			fslenergiesacceptedwcsim->Fill(primaryenergiesvector.at(indextouse));
			eventq2acceptedwcsim->Fill(q2vector.at(indextouse));
			muedepositionsacceptedwcsim->Fill(muonenergydepositions.at(indextouse));
			
			// fiducial cut versions
			if( (TMath::Sqrt(TMath::Power(genie_x,2)+TMath::Power(genie_z,2)) < fidcutradius) && 
				(TMath::Abs(genie_y-tank_yoffset) < fidcuty) && 
				(TMath::Abs(genie_z-tank_start-tank_radius) < fidcutz) ){
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
			if( (TMath::Sqrt(TMath::Power(genie_x,2)+TMath::Power(genie_z,2)) < fidcutradius) && 
				(TMath::Abs(genie_y-tank_yoffset) < fidcuty) && 
				(TMath::Abs(genie_z-tank_start-tank_radius) < fidcutz) ){
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
			cout<<"genie lepton angle = "<<fslanglegenie<<endl;
			cout<<"genie lepton energy = "<<fsleptonenergy<<endl;
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
			incidentneutrinoenergiesaccepted->Fill(probeenergy);
			//incidentneutrinoanglesaccepted->Fill(probeangle);
			fslanglesaccepted->Fill(fslanglegenie);
			fslenergiesaccepted->Fill(fsleptonenergy);
			eventq2accepted->Fill(Q2);
			neutrinovertexQEaccepted->Fill(genie_x, genie_y, genie_z);
			// if it passes the fiducial cut, we've also accepted it, fill.
			if( (TMath::Sqrt(TMath::Power(genie_x,2)+TMath::Power(genie_z,2)) < fidcutradius) && 
				(TMath::Abs(genie_y-tank_yoffset) < fidcuty) && 
				(TMath::Abs(genie_z-tank_start-tank_radius) < fidcutz) ){
				// genie values
				incidentneutrinoenergiesacceptedfidcut->Fill(probeenergy);
				fslanglesacceptedfidcut->Fill(fslanglegenie);
				fslenergiesacceptedfidcut->Fill(fsleptonenergy);
				eventq2acceptedfidcut->Fill(Q2);
				// fill flat tree for reconstruction dev
				MuonStartBranchFidMRD->Fill();
				MuonStopBranchFidMRD->Fill();
				MuonDirectionBranchFidMRD->Fill();
				DigitVertexBranchFidMRD->Fill();
				DigitChargeBranchFidMRD->Fill();
				numQEneutrinoeventsinfidvolmrd++;
				nummuontracksinfidvolmrd+=scatteringanglesvector.size();
			}
		}
		
		vertextreenocuts->SetEntries(MuonStartBranch->GetEntries());
		vertextreefiducialcut->SetEntries(MuonStartBranchFid->GetEntries());
		vertextreefiducialmrd->SetEntries(MuonStartBranchFidMRD->GetEntries());
		flateventfileout->Write("",TObject::kOverwrite);
//		flateventfileout->cd();
//		vertextreenocuts->Write();
//		vertextreefiducialcut->Write();
//		vertextreefiducialmrd->Write();
	}
	
	cout<<"generating scaled histograms"<<endl;
	Double_t numbeamspills = totalpots/(4.0 * TMath::Power(10.,12.));
	Double_t numbeamspillsperday = (24.*60.*60.*1000.)/133.3333;	// 24 hours in ms / 133.33 ms between spills
	Double_t numdays = numbeamspills/numbeamspillsperday;
	cout<<"Results based on "<<totalpots<<" POTs, or "<<numbeamspills<<" beam spills, or "<<numdays<<" days of data"<<endl;
	cout<<"There were "<<numneutrinoeventsintank<<" neutrino interactions in the tank, of which "<<numQEneutrinoeventsintank<<" were true QE events."<<endl;
	cout<<"Of those, "<<numQEneutrinoeventsinfidvol<<"were within the fiducial volume."<<endl;
	cout<<"Of those in turn, "<<numQEneutrinoeventsinfidvolmrd<<"produced an accepted MRD muon"<<endl;
	cout<<"There were "<<nummuontracksintank<<" muons in the tank, of which "
		<<nummuontracksinfidvol<<" were from events in the fiducial volume."<<endl;
	cout<<"There were "<<nummuontracksinmrd<<" muons from tank events which passed through 3 MRD layers"
		<<" of which "<<nummuontracksinfidvolmrd<<" originated from vertices in the fiducial volume"<<endl;
	
	
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
	
	// draw the scaled histograms
	// ==========================
	std::vector<TCanvas*> scaledcanvaspointers;
	std::vector<TLegend*> scaledlegendpointers;
	legendindex=0;
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
	if(fsltruetracklength) delete fsltruetracklength; fsltruetracklength=0;
	if(fsltruetracklengthintank) delete fsltruetracklengthintank; fsltruetracklengthintank=0;
	
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
#endif
	
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
		canv->SaveAs(TString::Format("%s/histograms_%d.png",outpath,i));
		if(canv) delete canv; canv=0;
		TLegend* leg = legendpointers.at(i);
		if(leg) delete leg; leg=0;
	}
	
	// save, then delete scaled canvases and legends
	//cout<<"deleting scaled canvases and legends"<<endl;
	for(int i=0; i<scaledcanvaspointers.size(); i++){
		TCanvas* canv = scaledcanvaspointers.at(i);
		canv->SaveAs(TString::Format("%s/scaled_histograms_%d.png",outpath,i));
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
	if(histofileout) histofileout->Close(); delete histofileout; histofileout=0;
	
	// delete flat file output
	//cout<<"deleting flat file output"<<endl;
	if(flateventfileout) flateventfileout->Close(); delete flateventfileout;
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

void FillTankMapHist(WCSimRootGeom* geo, int tubeID, TH1* histowall=0, TH1* histotop=0, TH1* histobottom=0, double weight=1){
	//Fill a bin on a 2D map of PMTs 
	WCSimRootPMT pmt = geo->GetPMT(tubeID - 1);	// TUBE ID NEEDS -1 IN GEO FILE
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
				if(histotop) histotop->Fill(thebins.first, thebins.second, digihit->GetT());
			} else {cout<<"bad pmt: ID "<<tubeID<<" in CylLoc "<<pmt.GetCylLoc()<<endl;}
			break;
		}
		case 1: {
			if(wallpositionmap.count(tubeID)){
				std::pair<int,int> thebins = wallpositionmap.at(tubeID);
				if(histowall) histowall->Fill(thebins.first+0.5, thebins.second, digihit->GetT());
			} else {cout<<"bad pmt: ID "<<tubeID<<" in CylLoc "<<pmt.GetCylLoc()<<endl;}
			break;
		}
		case 2: {
			if(bottomcappositionmap.count(tubeID)){
				std::pair<int,int> thebins = bottomcappositionmap.at(tubeID);
				if(histobottom) histobottom->Fill(thebins.first, thebins.second, digihit->GetT());
			} else {cout<<"bad pmt: ID "<<tubeID<<" in CylLoc "<<pmt.GetCylLoc()<<endl;}
			break;
		}
		case 4: {
//				std::pair<int,int> thebins = mrdpositionmap.at(tubeID);
//				mrdhist->Fill(thebins.first, thebins.second);
			break;
		}
		case 5: {
//				std::pair<int,int> thebins = faccpositionmap.at(tubeID);
//				facchist->Fill(thebins.first, thebins.second);
			break;
		}
		default: {
			//cout<<"PMT "<<tubeID<<" has unknown location "<<pmt.GetCylLoc()<<"!"<<endl; 
			break;
		}
	}
}
