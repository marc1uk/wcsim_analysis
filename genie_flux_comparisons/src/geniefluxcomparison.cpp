/* vim:set noexpandtab tabstop=2 wrap */
//C++
#include <iostream>
#include <vector>
#include <string>
#include <algorithm> // std::find
#include <limits>
//ROOT
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TString.h>
#include <TVector3.h>
#include <TPie.h>
#include <TAxis.h>
//GENIE
#include <FluxDrivers/GSimpleNtpFlux.h>
#include <FluxDrivers/GNuMIFlux.h>
// Custom
#include "MRDspecs.hh"
#include "ColourWheel.hh"
#ifndef FLUX_STAGE
//#define FLUX_STAGE 0 // gsimple / bnb_annie_* files (rhatcher versions)
#define FLUX_STAGE 1   // gntp files      (zarko versions)
#endif

#if FLUX_STAGE==1
// for fluxstage=1 (gntp files) we can retrieve info about the neutrino interaction
#include "../../genieinfo_struct.cpp"     // definition of a struct to hold genie info
#include <GHEP/GHepUtils.h>               // neut reaction codes
#include <PDG/PDGLibrary.h>
#include <Ntuple/NtpMCEventRecord.h>
#include <Conventions/Constants.h>
// function to fill the info
void GetGenieEntryInfo(genie::EventRecord* gevtRec, genie::Interaction* genieint, GenieInfo& thegenieinfo, bool printneutrinoevent=false);
#endif

// type conversion functions:
static std::map<int,std::string> pdgcodetoname;
static std::map<int,std::string> decaymap;
static std::map<int,std::string> gnumicodetoname;
static std::map<int,std::string>* GenerateGnumiMap();
static std::map<int,std::string>* GeneratePdgMap();
static std::map<int,std::string>* GenerateDecayMap();
std::string GnumiToString(int code);
std::string PdgToString(int code);
std::string DecayTypeToString(int code);
std::string MediumToString(int code);
TPie* GeneratePieFromHisto(TH1F* histo, int verbose);
TPie* GeneratePieFromHisto(std::string histoname, int verbose); // unused

using std::cout;
using std::cerr;
using std::endl;

//  main(string filepattern, int verbose, string fileoutname, string fileoutdir)
int main(int argc, char* argv[]){
	
//#if FLUX_STAGE==1
//	// turn down genie verbosity
//	genie::GHepRecord::SetPrintLevel(-1); // doesn't work - disable neut code in GetGenieEntryInfo
//#endif
	
	std::string filepattern = (argc>1) ? std::string(argv[1]) : "";
	char * pEnd;
	int verbose = (argc>2) ? strtol(argv[2],&pEnd, 0) : 0;
	if(argc>2 && pEnd==nullptr){
		cerr<<"invalid parameters: "<<endl
				<<"main(string filepattern, int verbose, string fileoutname, string fileoutdir)"<<endl;
				return 1;
	}
	std::string fileoutname = (argc>3) ? std::string(argv[2]) : "";
	std::string dirout = (argc>4) ? std::string(argv[3]) : "";
	
	int fluxstage = FLUX_STAGE;
	// create output file
	if(dirout=="") dirout = gSystem->pwd();
	if(fileoutname=="") fileoutname = "bnb_flux_comparisons.root";
	std::string fileoutpath = dirout + "/" + fileoutname;
	
	cout<<"saving results to "<<fileoutpath<<endl;
	TFile* fileout = new TFile(fileoutpath.c_str(),"RECREATE");
	
	// common input/output variables:
	int parentpdg;
	std::string parenttypestring;
	int parentdecaymode; // some arbitrary number that maps to a decay mode string.
	std::string parentdecaystring; // descriptive string. Should we store a map of the translation?
	float parentdecayvtx_x, parentdecayvtx_y, parentdecayvtx_z;
	TVector3 parentdecayvtx;
	float parentdecaymom_x, parentdecaymom_y, parentdecaymom_z;
	TVector3 parentdecaymom;
	float parentprodmom_x, parentprodmom_y, parentprodmom_z;
	TVector3 parentprodmom;
	int parentprodmedium;     // they're all 0
	std::string parentprodmediumstring; // do we even have this mapping?
	int parentpdgattgtexit;
	std::string parenttypestringattgtexit;
	TVector3 parenttgtexitmom;
	float parenttgtexitmom_x, parenttgtexitmom_y, parenttgtexitmom_z;
	
	std::string currentfilestring;
	uint currentevtnum; // entry num in local tree
	int fluxver;    // 0 = old flux, 1 = new flux
	
	long long maxentriestoanalyse=std::numeric_limits<long long>::max(); //20000;
	
	// variables in output tree
	cout<<"creating treeout"<<endl;
	TTree* treeout = new TTree("nutree","neutrino and parent properties");
	treeout->SetDirectory(fileout);
	
	treeout->Branch("file",&currentfilestring);
	treeout->Branch("fluxver",&fluxver);
	treeout->Branch("evtnum",&currentevtnum);
	treeout->Branch("ParentPdg",&parentpdg);
	treeout->Branch("ParentTypeString",&parenttypestring);
	treeout->Branch("ParentDecayMode",&parentdecaymode);
	treeout->Branch("ParentDecayString",&parentdecaystring);
	treeout->Branch("ParentDecayVtx",&parentdecayvtx);
	treeout->Branch("ParentDecayVtx_X",&parentdecayvtx_x);
	treeout->Branch("ParentDecayVtx_Y",&parentdecayvtx_y);
	treeout->Branch("ParentDecayVtx_Z",&parentdecayvtx_z);
	treeout->Branch("ParentDecayMom",&parentdecaymom);
	treeout->Branch("ParentDecayMom_X",&parentdecaymom_x);
	treeout->Branch("ParentDecayMom_Y",&parentdecaymom_y);
	treeout->Branch("ParentDecayMom_Z",&parentdecaymom_z);
	treeout->Branch("ParentProdMom",&parentprodmom);
	treeout->Branch("ParentProdMom_X",&parentprodmom_x);
	treeout->Branch("ParentProdMom_Y",&parentprodmom_y);
	treeout->Branch("ParentProdMom_Z",&parentprodmom_z);
	//treeout->Branch("ParentProdMedium",&parentprodmedium);
	//treeout->Branch("ParentProdMediumString",&parentprodmediumstring);
	treeout->Branch("ParentPdgAtTgtExit",&parentpdgattgtexit);
	treeout->Branch("ParentTypeAtTgtExitString",&parenttypestringattgtexit);
	treeout->Branch("ParentTgtExitMom",&parenttgtexitmom);
	treeout->Branch("ParentTgtExitMom_X",&parenttgtexitmom_x);
	treeout->Branch("ParentTgtExitMom_Y",&parenttgtexitmom_y);
	treeout->Branch("ParentTgtExitMom_Z",&parenttgtexitmom_z);
	
#if FLUX_STAGE==1
	// store the neutrino info from gntp files
	// a load of variables to specify interaction type
	bool IsQuasiElastic=false;
	bool IsResonant=false;
	bool IsDeepInelastic=false;
	bool IsCoherent=false;
	bool IsDiffractive=false;
	bool IsInverseMuDecay=false;
	bool IsIMDAnnihilation=false;
	bool IsSingleKaon=false;
	bool IsNuElectronElastic=false;
	bool IsEM=false;
	bool IsWeakCC=false;
	bool IsWeakNC=false;
	bool IsMEC=false;
	std::string interactiontypestring="";
	int neutcode=-1;
	treeout->Branch("IsQuasiElastic",&IsQuasiElastic);
	treeout->Branch("IsResonant",&IsResonant);
	treeout->Branch("IsDeepInelastic",&IsDeepInelastic);
	treeout->Branch("IsCoherent",&IsCoherent);
	treeout->Branch("IsDiffractive",&IsDiffractive);
	treeout->Branch("IsInverseMuDecay",&IsInverseMuDecay);
	treeout->Branch("IsIMDAnnihilation",&IsIMDAnnihilation);
	treeout->Branch("IsSingleKaon",&IsSingleKaon);
	treeout->Branch("IsNuElectronElastic",&IsNuElectronElastic);
	treeout->Branch("IsEM",&IsEM);
	treeout->Branch("IsWeakCC",&IsWeakCC);
	treeout->Branch("IsWeakNC",&IsWeakNC);
	treeout->Branch("IsMEC",&IsMEC);
	treeout->Branch("InteractionTypeString",&interactiontypestring);
	treeout->Branch("NeutCode",&neutcode);
	// ok, moving on
	double nuIntxVtx_X; // cm
	double nuIntxVtx_Y; // cm
	double nuIntxVtx_Z; // cm
	double nuIntxVtx_T; // ns
	bool isintank=false;
	bool isinfiducialvol=false;
	double eventq2=-1;
	double eventEnu=-1;
	int neutrinopdg=-1;
	double muonenergy=-1;
	double muonangle=-1;
	std::string fsleptonname; // assumed to be muon, but we should confirm
	// these may not be properly copied... 
	int numfsprotons;
	int numfsneutrons;
	int numfspi0;
	int numfspiplus;
	int numfspiminus;
	treeout->Branch("NuIntxVtx_X",&nuIntxVtx_X);
	treeout->Branch("NuIntxVtx_Y",&nuIntxVtx_Y);
	treeout->Branch("NuIntxVtx_Z",&nuIntxVtx_Z);
	treeout->Branch("NuIntxVtx_T",&nuIntxVtx_T);
	treeout->Branch("NuVtxInTank",&isintank);
	treeout->Branch("NuVtxInFidVol",&isinfiducialvol);
	treeout->Branch("EventQ2",&eventq2);
	treeout->Branch("NeutrinoEnergy",&eventEnu);
	treeout->Branch("NeutrinoPDG",&neutrinopdg);
	treeout->Branch("MuonEnergy",&muonenergy);
	treeout->Branch("MuonAngle",&muonangle);
	treeout->Branch("FSLeptonName",&fsleptonname);
	treeout->Branch("NumFSProtons",&numfsprotons);
	treeout->Branch("NumFSNeutrons",&numfsneutrons);
	treeout->Branch("NumFSPi0",&numfspi0);
	treeout->Branch("NumFSPiPlus",&numfspiplus);
	treeout->Branch("NumFSPiMinus",&numfspiminus);
#endif // FLUX_STAGE==1
	
	TChain* oldflux;
	genie::flux::GNuMIFluxPassThroughInfo* gnumipassthruentry  = nullptr;
	genie::NtpMCEventRecord* genieintx = nullptr; // = new genie::NtpMCEventRecord;
	if(fluxstage==0){ // use stage 0 bnb_annie_*.root flux files
		oldflux = new TChain("h10");
		oldflux->Add("/pnfs/annie/persistent/flux/bnb/bnb_annie_*.root");
		cout<<"old flux has "<<oldflux->GetEntries()<<" entries"<<endl;
		cout<<"setting oldflux branch addresses"<<endl;
		oldflux->SetBranchAddress("ptype",&parentpdg);
		oldflux->SetBranchAddress("Ndecay",&parentdecaymode);
		oldflux->SetBranchAddress("Vx",&parentdecayvtx_x);
		oldflux->SetBranchAddress("Vy",&parentdecayvtx_y);
		oldflux->SetBranchAddress("Vz",&parentdecayvtx_z);
		oldflux->SetBranchAddress("pdpx",&parentdecaymom_x);
		oldflux->SetBranchAddress("pdpy",&parentdecaymom_y);
		oldflux->SetBranchAddress("pdpz",&parentdecaymom_z);
		oldflux->SetBranchAddress("ppdxdz",&parentprodmom_x);  // needs fixing in fluxver 0, fluxstage 0
		oldflux->SetBranchAddress("ppdydz",&parentprodmom_y);  // needs fixing in fluxver 0, fluxstage 0
		oldflux->SetBranchAddress("pppz",&parentprodmom_z);
		//oldflux->SetBranchAddress("ppmedium",&parentprodmedium);
		oldflux->SetBranchAddress("tptype",&parentpdgattgtexit);
		oldflux->SetBranchAddress("tpx",&parenttgtexitmom_x);
		oldflux->SetBranchAddress("tpy",&parenttgtexitmom_y);
		oldflux->SetBranchAddress("tpz",&parenttgtexitmom_z);
	} else {    // use stage 1 gntp.*.ghep.root genie files
		oldflux = new TChain("gtree");
		oldflux->Add("/pnfs/annie/persistent/users/vfischer/genie/BNB_Water_10k_22-05-17/gntp.10??.ghep.root"); // XXX replace with 100? to use a reasonable number of entries for debug
		cout<<"old flux has "<<oldflux->GetEntries()<<" entries"<<endl;
		cout<<"setting oldflux branch addresses"<<endl;
		oldflux->SetBranchAddress("flux",&gnumipassthruentry);
		oldflux->GetBranch("flux")->SetAutoDelete(kTRUE);
#if FLUX_STAGE==1
		oldflux->SetBranchAddress("gmcrec",&genieintx);
		oldflux->GetBranch("gmcrec")->SetAutoDelete(kTRUE);
#endif // FLUX_STAGE==1
	}
	
	TChain* newflux;
	genie::flux::GSimpleNtpNuMI* gsimplenumientry = nullptr;
	if(fluxstage==0){
		newflux = new TChain("flux");
		newflux->Add("/annie/data/flux/gsimple_bnb/gsimple_beammc_annie_*.root");
		cout<<"new flux has "<<newflux->GetEntries()<<" entries"<<endl;
		cout<<"setting branch address"<<endl;
		newflux->SetBranchAddress("numi",&gsimplenumientry);
		newflux->GetBranch("numi")->SetAutoDelete(kTRUE);
	} else {
		newflux = new TChain("gtree");
		newflux->Add("/pnfs/annie/persistent/users/moflaher/genie/BNB_World_10k_11-03-18_gsimpleflux/gntp.*.ghep.root");
		cout<<"new flux has "<<newflux->GetEntries()<<" entries"<<endl;
		cout<<"setting branch address"<<endl;
		newflux->SetBranchAddress("numi",&gsimplenumientry);
		newflux->GetBranch("numi")->SetAutoDelete(kTRUE);
#if FLUX_STAGE==1
		newflux->SetBranchAddress("gmcrec",&genieintx);
		newflux->GetBranch("gmcrec")->SetAutoDelete(kTRUE);
#endif // FLUX_STAGE==1
	}
	
	// first copy over information from original flux files
	fluxver=0;
	TFile* curf=nullptr;
	TFile* curflast=nullptr;
	long long int oldfluxentries = oldflux->GetEntries();
	cout<<"copying "<<oldfluxentries<<" oldflux entries"<<endl;
	for(int i=0; i<std::min(oldfluxentries,maxentriestoanalyse); i++){
		currentevtnum = oldflux->LoadTree(i);
		curf = oldflux->GetCurrentFile();
		if(curf!=curflast || curflast==nullptr){
			TString curftstring = curf->GetName();
			currentfilestring = std::string(curftstring.Data());
			curflast=curf;
		}
		if((i%1000)==0) cout<<"i="<<i<<" = "<<(((double)i/(double)oldfluxentries)*100.)<<"%"<<endl;
		oldflux->GetEntry(i);
		
		if(fluxstage==0){
			//nothing to do, variables are read directly
		} else if(fluxstage==1){
			// if using gntp files we need to extract the flux from the gnumifluxpassthrough object
			parentpdg = gnumipassthruentry->ptype;
			parentdecaymode = gnumipassthruentry->ndecay;
			parentdecayvtx_x = gnumipassthruentry->vx;
			parentdecayvtx_y = gnumipassthruentry->vy;
			parentdecayvtx_z = gnumipassthruentry->vz;
			parentdecaymom_x = gnumipassthruentry->pdpx;
			parentdecaymom_y = gnumipassthruentry->pdpy;
			parentdecaymom_z = gnumipassthruentry->pdpz;
			parentprodmom_x = gnumipassthruentry->ppdxdz;
			parentprodmom_y = gnumipassthruentry->ppdydz;
			parentprodmom_z = gnumipassthruentry->pppz;
			//parentprodmedium = gnumipassthruentry->ppmedium;
			parentpdgattgtexit = gnumipassthruentry->tptype;
			parenttgtexitmom_x = gnumipassthruentry->tpx;
			parenttgtexitmom_y = gnumipassthruentry->tpy;
			parenttgtexitmom_z = gnumipassthruentry->tpz;
		}
		
		// convenience type conversions
		parentdecayvtx = TVector3(parentdecayvtx_x,parentdecayvtx_y,parentdecayvtx_z);
		parentdecaymom = TVector3(parentdecaymom_x,parentdecaymom_y,parentdecaymom_z);
		parentprodmom = TVector3(parentprodmom_x,parentprodmom_y,parentprodmom_z);
		parenttgtexitmom = TVector3(parenttgtexitmom_x,parenttgtexitmom_y,parenttgtexitmom_z);
		parenttypestring = (fluxstage==0) ? GnumiToString(parentpdg) : PdgToString(parentpdg);
		parenttypestringattgtexit = (fluxstage==0) ? GnumiToString(parentpdgattgtexit) : PdgToString(parentpdgattgtexit);
		parentdecaystring = DecayTypeToString(parentdecaymode);
		//parentprodmediumstring = MediumToString(parentprodmedium);
		
#if FLUX_STAGE==1
		// neutrino interaction info
		genie::EventRecord* gevtRec = genieintx->event;
		genie::Interaction* genieint = gevtRec->Summary();
		//cout<<"scraping interaction info"<<endl;
		GenieInfo thegenieinfo;
		GetGenieEntryInfo(gevtRec, genieint, thegenieinfo);  // fill convenience struct with nu intx info
		//cout<<"done interaction info"<<endl;
		
		// retrieve info from the struct
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
		interactiontypestring=thegenieinfo.interactiontypestring;
		neutcode=thegenieinfo.neutinteractioncode; // currently disabled to prevent excessive verbosity
		
		eventq2=thegenieinfo.Q2;
		eventEnu=thegenieinfo.probeenergy;
		neutrinopdg=thegenieinfo.probepdg;
		muonenergy=thegenieinfo.fsleptonenergy;
		muonangle=thegenieinfo.fslangle;
		
		nuIntxVtx_X=thegenieinfo.Intx_x; // cm
		nuIntxVtx_Y=thegenieinfo.Intx_y; // cm
		nuIntxVtx_Z=thegenieinfo.Intx_z; // cm
		nuIntxVtx_T=thegenieinfo.Intx_t; // ns
		// check in tank
		if( ( sqrt( pow(nuIntxVtx_X, 2) + pow(nuIntxVtx_Z-MRDSpecs::tank_start-MRDSpecs::tank_radius,2) )
			  < MRDSpecs::tank_radius ) && 
			  ( abs(nuIntxVtx_Y-MRDSpecs::tank_yoffset) < MRDSpecs::tank_halfheight) ){
			isintank=true;
		} else { isintank=false; }
		// check in fiducial volume
		if( isintank &&
		  ( sqrt (pow(nuIntxVtx_X, 2) + pow(nuIntxVtx_Z-MRDSpecs::tank_start-MRDSpecs::tank_radius,2)) 
		  < MRDSpecs::fidcutradius ) && 
		  ( abs(nuIntxVtx_Y-MRDSpecs::tank_yoffset) < MRDSpecs::fidcuty ) && 
		  ( (nuIntxVtx_Z-MRDSpecs::tank_start-MRDSpecs::tank_radius) < MRDSpecs::fidcutz) ){
			isinfiducialvol=true;
		} else { isinfiducialvol = false; }
		
		fsleptonname = std::string(thegenieinfo.fsleptonname.Data());
		// this data does not appear to be populated...
		numfsprotons = thegenieinfo.numfsprotons = genieint->ExclTag().NProtons();
		numfsneutrons = thegenieinfo.numfsneutrons = genieint->ExclTag().NNeutrons();
		numfspi0 = thegenieinfo.numfspi0 = genieint->ExclTag().NPi0();
		numfspiplus = thegenieinfo.numfspiplus = genieint->ExclTag().NPiPlus();
		numfspiminus = thegenieinfo.numfspiminus = genieint->ExclTag().NPiMinus();
		
		genieintx->Clear(); // REQUIRED TO PREVENT MEMORY LEAK
#endif // FLUX_STAGE==1
		
		treeout->Fill();
	}
	int copiedoldentries=treeout->GetEntries();
	cout<<"copied "<<copiedoldentries<<" oldflux entries"<<endl;
	
	// now copy over information from new flux files
	// exactly the same, but this time we have to extract the variables from the GSimpleNtpNuMI
	fluxver=1;
	curf=nullptr;
	curflast=nullptr;
	long long int newfluxentries = newflux->GetEntries();
	cout<<"copying newflux entries"<<endl;
	for(int i=0; i<std::min(newfluxentries,maxentriestoanalyse); i++){
		currentevtnum = newflux->LoadTree(i);
		curf = newflux->GetCurrentFile();
		if(curf!=curflast || curflast==nullptr){
			TString curftstring = curf->GetName();
			currentfilestring = std::string(curftstring.Data());
		}
		if((i%1000)==0) cout<<"j="<<i<<" = "<<(((double)i/(double)newfluxentries)*100.)<<"%"<<endl;
		newflux->GetEntry(i);
		
		if(fluxstage==0||fluxstage==1){
			// extract the variables from the GSimpleNtpNuMI object
			parentpdg = gsimplenumientry->ptype;
			parentdecaymode = gsimplenumientry->ndecay;
			parentdecayvtx_x = gsimplenumientry->vx;
			parentdecayvtx_y = gsimplenumientry->vy;
			parentdecayvtx_z = gsimplenumientry->vz;
			parentdecaymom_x = gsimplenumientry->pdpx;
			parentdecaymom_y = gsimplenumientry->pdpy;
			parentdecaymom_z = gsimplenumientry->pdpz;
			parentprodmom_x = gsimplenumientry->pppx/gsimplenumientry->pppz; // ??? is this ppdxdz?
			parentprodmom_y = gsimplenumientry->pppy/gsimplenumientry->pppz;
			parentprodmom_z = gsimplenumientry->pppz;
			//parentprodmedium = gsimplenumientry->ppmedium;
			parentpdgattgtexit = gsimplenumientry->tptype;
			parenttgtexitmom_x = gsimplenumientry->tpx;
			parenttgtexitmom_y = gsimplenumientry->tpy;
			parenttgtexitmom_z = gsimplenumientry->tpz;
		}
		
		// convenience type conversions
		parentdecayvtx = TVector3(parentdecayvtx_x,parentdecayvtx_y,parentdecayvtx_z);
		parentdecaymom = TVector3(parentdecaymom_x,parentdecaymom_y,parentdecaymom_z);
		parentprodmom = TVector3(parentprodmom_x,parentprodmom_y,parentprodmom_z);
		parenttgtexitmom = TVector3(parenttgtexitmom_x,parenttgtexitmom_y,parenttgtexitmom_z);
		parenttypestring = PdgToString(parentpdg);
		parenttypestringattgtexit = PdgToString(parentpdgattgtexit);
		parentdecaystring = DecayTypeToString(parentdecaymode);
		//parentprodmediumstring = MediumToString(parentprodmedium);
		
#if FLUX_STAGE==1
		// neutrino interaction info
		genie::EventRecord* gevtRec = genieintx->event;
		genie::Interaction* genieint = gevtRec->Summary();
		//cout<<"scraping interaction info"<<endl;
		GenieInfo thegenieinfo;
		GetGenieEntryInfo(gevtRec, genieint, thegenieinfo);  // fill convenience struct with nu intx info
		//cout<<"done interaction info"<<endl;
		
		// retrieve info from the struct
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
		interactiontypestring=thegenieinfo.interactiontypestring;
		neutcode=thegenieinfo.neutinteractioncode; // currently disabled to prevent excessive verbosity
		
		eventq2=thegenieinfo.Q2;
		eventEnu=thegenieinfo.probeenergy;
		neutrinopdg=thegenieinfo.probepdg;
		muonenergy=thegenieinfo.fsleptonenergy;
		muonangle=thegenieinfo.fslangle;
		
		nuIntxVtx_X=thegenieinfo.Intx_x; // cm
		nuIntxVtx_Y=thegenieinfo.Intx_y; // cm
		nuIntxVtx_Z=thegenieinfo.Intx_z; // cm
		nuIntxVtx_T=thegenieinfo.Intx_t; // ns
		// check in tank
		if( ( sqrt( pow(nuIntxVtx_X, 2) + pow(nuIntxVtx_Z-MRDSpecs::tank_start-MRDSpecs::tank_radius,2) )
			  < MRDSpecs::tank_radius ) && 
			  ( abs(nuIntxVtx_Y-MRDSpecs::tank_yoffset) < MRDSpecs::tank_halfheight) ){
			isintank=true;
		} else { isintank=false; }
		// check in fiducial volume
		if( isintank &&
		  ( sqrt (pow(nuIntxVtx_X, 2) + pow(nuIntxVtx_Z-MRDSpecs::tank_start-MRDSpecs::tank_radius,2)) 
		  < MRDSpecs::fidcutradius ) && 
		  ( abs(nuIntxVtx_Y-MRDSpecs::tank_yoffset) < MRDSpecs::fidcuty ) && 
		  ( (nuIntxVtx_Z-MRDSpecs::tank_start-MRDSpecs::tank_radius) < MRDSpecs::fidcutz) ){
			isinfiducialvol=true;
		} else { isinfiducialvol = false; }
		
		fsleptonname = std::string(thegenieinfo.fsleptonname.Data());
		// this data does not appear to be populated...
		numfsprotons = thegenieinfo.numfsprotons = genieint->ExclTag().NProtons();
		numfsneutrons = thegenieinfo.numfsneutrons = genieint->ExclTag().NNeutrons();
		numfspi0 = thegenieinfo.numfspi0 = genieint->ExclTag().NPi0();
		numfspiplus = thegenieinfo.numfspiplus = genieint->ExclTag().NPiPlus();
		numfspiminus = thegenieinfo.numfspiminus = genieint->ExclTag().NPiMinus();
		
		genieintx->Clear(); // REQUIRED TO PREVENT MEMORY LEAK
#endif // FLUX_STAGE==1
		
		treeout->Fill();
	}
	int copiednewentries = treeout->GetEntries()-copiedoldentries;
	cout<<"copied "<<copiednewentries<<" newflux entries"<<endl;
	cout<<"treeout has "<<treeout->GetEntries()<<" combined entries"<<endl;
	
	treeout->Write();
	treeout->ResetBranchAddresses();
	
	// make any desired histograms
	TCanvas* c1 = new TCanvas("c1", "c1",13,50,1000,600);
	c1->cd();
	
	// list all parameters to make histograms of here
	std::vector<std::string> paramnames{"ParentPdg", "ParentTypeString", "ParentDecayMode", "ParentDecayString", "ParentDecayVtx_X", "ParentDecayVtx_Y", "ParentDecayVtx_Z", "ParentDecayMom_X", "ParentDecayMom_Y", "ParentDecayMom_Z", "ParentProdMom_X", "ParentProdMom_Y", "ParentProdMom_Z", "ParentPdgAtTgtExit", "ParentTypeAtTgtExitString", "ParentTgtExitMom_X", "ParentTgtExitMom_Y", "ParentTgtExitMom_Z"};
	
	cout<<"creating histograms"<<endl;
	for(fluxver=0; fluxver<2; fluxver++){
		for(std::string paramname : paramnames){
			c1->Clear();
			treeout->Draw(TString::Format("%s>>h%s%d",paramname.c_str(),paramname.c_str(),fluxver),TString::Format("fluxver==%d",fluxver));
			TH1F* temphist = (TH1F*)gROOT->FindObject(TString::Format("h%s%d",paramname.c_str(),fluxver));
			
			bool isstringhisto = ((paramname.find("String"))!=std::string::npos);
			if(isstringhisto){
				// to overlay the histos, we need to ensure they have the same bin labels and bins are in
				// the same order - which is not necessarily the case for string bins
				// to ensure they have the same bin labels, we need to add bins for all empty string labels
				// (as tree->Draw() will not create these bins automatically) and then sort the axis
				
				// first get the map of all possible labels, and note any populated (i.e. existing) bins
				std::map<std::string,double> binlabelsandconts; // even for TH1F bin contents returned by
				std::map<int,std::string>* amap;      // GetBinContent are returned as DOUBLE
				std::size_t found = paramname.find("Decay");  // then SetBinContent takes a double
				if(found!=std::string::npos) amap = GenerateDecayMap();  // and internally casts it to float!
				else amap = GenerateGnumiMap();
				for(auto apair : (*amap)) binlabelsandconts.emplace(apair.second,0.0);
				
				// scan existing bins and note which ones are present
				if(verbose) cout<<"recording contents of "<<temphist->GetName()
												<<", which has "<<temphist->GetNbinsX()<<"bins"<<endl;
				
				// XXX NOTE: USED BINS ARE NUMBERED FROM 1 (0 is underflow)
				// AND MAX USED BIN NUMBER IS = GETNBINSX() 
				// SO TO SCAN USED BINS, LOOP OVER bini=0; bini<GetNbinsX(), but ALWAYS USE (BINI+1)
				for(int bini=0; bini<temphist->GetNbinsX(); bini++){
					TString binlabel = temphist->GetXaxis()->GetBinLabel(bini+1);
					if(binlabelsandconts.count(binlabel.Data())==0){
						cerr<<"missing alphanumeric histogram label "<<binlabel.Data()<<endl;
						binlabelsandconts.emplace(binlabel.Data(),temphist->GetBinContent(bini+1));
					} else {
						if(binlabelsandconts.at(binlabel.Data())!=0) 
							cerr<<"Repeated Set call for bin "<<binlabel.Data()<<endl;
						binlabelsandconts.at(binlabel.Data())=temphist->GetBinContent(bini+1);
					}
					if(verbose) cout<<"putting in map bin "<<(bini+1)<<", "<<binlabel.Data()<<" has "
													<<temphist->GetBinContent(bini+1)<<" entries"<<endl;
				}
				
				// now add any labels that aren't already present
				temphist->SetCanExtend(TH1::kAllAxes); // << allow SetBinContent to create new bins
				int bini=temphist->GetNbinsX()+1;  // doesn't set the overflow bin conts when above is set
				for(auto abin : binlabelsandconts){
					if(verbose) cout<<"bin label "<<abin.first.c_str()<<" was recorded to have "
													<<abin.second<<" entries."<<endl;
					if(abin.second!=0) continue;
					temphist->SetBinContent(bini, abin.second); // call first to create the new bin
					temphist->GetXaxis()->SetBinLabel(bini, abin.first.c_str()); // set the new bin's label
					if(verbose) cout<<"Setting bin "<<bini<<" label to "<<abin.first.c_str()
													<<" and contents to "<<abin.second<<endl;
					bini++;
				}
				// when we ask for one new bin to be created, the num bins is doubled.
				// this means we end up with many empty bins that we need to trim by deflating the axis:
				temphist->LabelsDeflate("X");
				
				if(verbose){
					for(int bini=1; bini<temphist->GetNbinsX()+1; bini++){
						TString binlabel = temphist->GetXaxis()->GetBinLabel(bini);
						double binconts = temphist->GetBinContent(bini);
						cout<<"bin "<<bini<<" has label "<<binlabel.Data()<<" and conts "<<binconts<<endl;
					}
				}
				
				// finally sort the bins alphabetically by label to ensure same ordering
				temphist->LabelsOption("a", "X");
				
			}
			
			// line colours don't seem to be saved in the histogram? but can be saved in a canvas?
			(fluxver==0) ? temphist->SetLineColor(kRed) : temphist->SetLineColor(kBlue);
			temphist->SetTitle(TString::Format("%s_%d",paramname.c_str(),fluxver));
			cout<<"writing "<<paramname<<" hist to file"<<endl;
			temphist->Write();
			
			// make pie charts from string histograms while we have them
			if(paramname.find("String")!=std::string::npos){
				TPie* thepie = GeneratePieFromHisto(temphist, verbose);
				cout<<"writing "<<paramname<<" pie chart to file"<<endl;
				thepie->Write();
				
				delete thepie;
			}
		}
	}
	
	// not all are same format / comparable
	std::vector<std::string> compnames{"ParentTypeString", "ParentDecayMode", "ParentDecayString", "ParentDecayVtx_X", "ParentDecayVtx_Y", "ParentDecayVtx_Z", "ParentDecayMom_X", "ParentDecayMom_Y", "ParentDecayMom_Z", "ParentProdMom_X", "ParentProdMom_Y", "ParentProdMom_Z", "ParentTypeAtTgtExitString", "ParentTgtExitMom_X", "ParentTgtExitMom_Y", "ParentTgtExitMom_Z"};
	cout<<"saving overlaid normalized plots"<<endl;
	for(std::string paramname : compnames){
		c1->Clear();
		
		cout<<"making comparison plot for "<<paramname<<endl;
		
		fluxver=0;
		// set line colours, normalise and overlay
		TH1F* temphist0 = (TH1F*)gROOT->FindObject(TString::Format("h%s%d",paramname.c_str(),fluxver));
		fluxver=1;
		TH1F* temphist1 = (TH1F*)gROOT->FindObject(TString::Format("h%s%d",paramname.c_str(),fluxver));
		
		temphist0->Scale(1./(temphist0->Integral()));
		temphist0->SetLineColor(kRed);
		temphist1->Scale(1./(temphist1->Integral()));
		temphist1->SetLineColor(kBlue);
		
		if( temphist0->GetMaximum() > temphist1->GetMaximum() ){
			temphist0->Draw(); temphist1->Draw("same");
		} else {
			temphist1->Draw(); temphist0->Draw("same");
		}
		
		TLegend* aleg = c1->BuildLegend();
		aleg->SetFillStyle(0); // no background colour to legend
		c1->SetTitle(TString::Format("c%s",paramname.c_str()));
		c1->SetName(TString::Format("c%s",paramname.c_str()));
		c1->Write();
		
		// if it's a string histo, we have a corresponding pie chart:
		bool ispie = paramname.find("String")!=std::string::npos;
		if(ispie){
			// Place pies side-by-side instead of overlay
			fluxver=0;
			TH1F* temphist0 = (TH1F*)gROOT->FindObject(TString::Format("h%s%dPie",paramname.c_str(),fluxver));
			fluxver=1;
			TH1F* temphist1 = (TH1F*)gROOT->FindObject(TString::Format("h%s%dPie",paramname.c_str(),fluxver));
			
			c1->Clear();
			c1->Divide(2,1); // n.b. divisions will be removed when we clear the canvas
			gStyle->SetOptTitle(1); // in place of legends, show histogram titles. maybe need to set at view
			c1->cd(1);
			temphist0->Draw();
			c1->cd(2);
			temphist1->Draw();
			
			c1->SetTitle(TString::Format("c%sPie",paramname.c_str()));
			c1->SetName(TString::Format("c%sPie",paramname.c_str()));
			c1->Write();
			
			// while we have both pies, make a normal string histo with bins containing
			// the minimal set of required labels
			// first get the set of required bin labels and contents from the two pie charts
			TPie* temphist0pie = (TPie*)temphist0; // need to covert to TPie to use appropriate methods
			TPie* temphist1pie = (TPie*)temphist1;
			std::map<std::string,std::pair<double,double>> binlabelsandconts;
			for(int bini=0; bini<temphist0pie->GetEntries(); bini++){
				std::string binlabel = std::string(temphist0pie->GetEntryLabel(bini));
				double binconts = temphist0pie->GetEntryVal(bini);
				if(binlabelsandconts.count(binlabel.c_str())==0){
					binlabelsandconts.emplace(binlabel,std::pair<double,double>{binconts,0});
				} else {
					cerr<<"pie conversion loop 1 read label "<<binlabel<<" more than once"<<endl;
				}
			}
			for(int bini=0; bini<temphist1pie->GetEntries(); bini++){
				std::string binlabel = std::string(temphist1pie->GetEntryLabel(bini));
				double binconts = temphist1pie->GetEntryVal(bini);
				if(binlabelsandconts.count(binlabel.c_str())==0){
					binlabelsandconts.emplace(binlabel,std::pair<double,double>{0,binconts});
				} else {
					if(binlabelsandconts.at(binlabel).second==0){
						binlabelsandconts.at(binlabel).second=binconts;
					} else {
						cerr<<"pie conversion loop 2 read label "<<binlabel<<" more than once"<<endl;
					}
				}
			}
			
			// now construct two string histos, using the combined bin set
			std::vector<TH1F*> histos;
			for(fluxver=0; fluxver<2; fluxver++){
				TString hname = TString::Format("h%s%dMinimal",paramname.c_str(),fluxver);
				TH1F* histo2 = new TH1F(hname,hname,binlabelsandconts.size(),0,binlabelsandconts.size());
				// set the bin contents
				int bini=0;
				for(auto apair : binlabelsandconts){
					std::string thebinlabel = apair.first;
					double thebincontents = (fluxver==0) ? apair.second.first : apair.second.second;
					histo2->SetBinContent(bini+1,thebincontents);
					histo2->GetXaxis()->SetBinLabel(bini+1,thebinlabel.c_str());
					bini++;
				}
				histo2->LabelsOption("a", "X"); // sort the labels alphabetically for consistency
				
				(fluxver==0) ? histo2->SetLineColor(kRed) : histo2->SetLineColor(kBlue);
				cout<<"writing "<<hname.Data()<<" hist to file"<<endl;
				histo2->Write();
				
				// keep hold of them and we'll overlay them on a canvas next
				histos.push_back(histo2);
			}
			
			// now overlay the two on a canvas and save that too
			c1->Clear();
			c1->cd();
			temphist0 = histos.at(0); temphist1 = histos.at(1);
			temphist0->Scale(1./(temphist0->Integral()));
			temphist0->SetLineColor(kRed);
			temphist1->Scale(1./(temphist1->Integral()));
			temphist1->SetLineColor(kBlue);
			
			if( temphist0->GetMaximum() > temphist1->GetMaximum() ){
				temphist0->Draw(); temphist1->Draw("same");
			} else {
				temphist1->Draw(); temphist0->Draw("same");
			}
			
			TLegend* aleg = c1->BuildLegend();
			aleg->SetFillStyle(0); // no background colour to legend
			TString hname = TString::Format("c%sMinimal",paramname.c_str());
			c1->SetTitle(hname);
			c1->SetName(hname);
			c1->Write();
			
			// delete the temporary histos made
			for(TH1F* histo2 : histos) delete histo2;
		}
	}
	
	// TODO: ADD STATISTICAL UNCERTAINTIES so we can check they're consistent
	// TODO flux uncertainties vs energy, stacked for parent type, decay mode
	
	// make some plots placing cuts on neutrino flavour and decay mode
	// ==============
	cout<<"saving stacked plots"<<endl;
	std::map<std::string,int> nutypes;
	nutypes.emplace("Nue",12);
	nutypes.emplace("AntiNue",-12);
	nutypes.emplace("Numu",14);
	nutypes.emplace("AntiNumu",-14);
	for(fluxver=0; fluxver<2; fluxver++){
		// loop over flux versions
		int totalnuintx = (fluxver) ? copiednewentries : copiedoldentries;
		//double totalpots = (fluxver) ? newfluxpots : oldfluxpots;
		//double totalflux = (double)totalnuintx / (totalpots*detectortonnes); // flux per POT per tonne
		double fluxscaling = 1.; // disable for the moment
		
		//	nu flux (energy), stacked by neutrino type. total, and for CC / NC interactions?
		THStack* flavourstack = new THStack(TString::Format("NuFluxByFlavour%d",fluxver),"Nu Flux by Flavour");
		ColourWheel flavourcolourwheel = ColourWheel();
		THStack* ccflavourstack = new THStack(TString::Format("NuFluxByFlavour%dCC",fluxver), "Nu Flux by Flavour for CC Intx");
		ColourWheel ccflavourcolourwheel = ColourWheel();
		THStack* ncflavourstack = new THStack(TString::Format("NuFluxByFlavour%dNC",fluxver), "Nu Flux by Flavour for NC Intx");
		ColourWheel ncflavourcolourwheel = ColourWheel();
		TLegend* leg;
		
		// loop over the flavours, make a histo and add it to the stack
		for(auto nutype : nutypes){
			std::string nustring = nutype.first;
			int nupdg = nutype.second;
			
			// version with no CC / NC splitting
			c1->Clear();
			treeout->Draw(TString::Format("NeutrinoEnergy>>h%sFlux%d(100,0.,5.)",nustring.c_str(),fluxver), TString::Format("fluxver==%d&&NeutrinoPDG==%d",fluxver,nupdg));
			TH1F* temphist = (TH1F*)gROOT->FindObject(TString::Format("h%sFlux%d",nustring.c_str(),fluxver));
			temphist->Scale(fluxscaling); // XXX XXX XXX XXX XXX XXX XXX XXX
			temphist->SetFillColor(flavourcolourwheel.GetNextColour());
			temphist->SetFillStyle(3002); // semi dense spots
			flavourstack->Add(temphist);
			
			// make version with CC cut
			c1->Clear();
			treeout->Draw(TString::Format("NeutrinoEnergy>>h%sCCFlux%d(100,0.,5.)",nustring.c_str(),fluxver), TString::Format("fluxver==%d&&NeutrinoPDG==%d&&IsWeakCC==1",fluxver,nupdg));
			temphist = (TH1F*)gROOT->FindObject(TString::Format("h%sCCFlux%d",nustring.c_str(),fluxver));
			temphist->Scale(fluxscaling); // XXX XXX XXX XXX XXX XXX XXX XXX
			temphist->SetFillColor(ccflavourcolourwheel.GetNextColour());
			temphist->SetFillStyle(3002); // semi dense spots
			ccflavourstack->Add(temphist);
			
			// make version with NC cut
			c1->Clear();
			treeout->Draw(TString::Format("NeutrinoEnergy>>h%sNCFlux%d(100,0.,5.)",nustring.c_str(),fluxver), TString::Format("fluxver==%d&&NeutrinoPDG==%d&&IsWeakNC==1",fluxver,nupdg));
			temphist = (TH1F*)gROOT->FindObject(TString::Format("h%sNCFlux%d",nustring.c_str(),fluxver));
			temphist->Scale(fluxscaling); // XXX XXX XXX XXX XXX XXX XXX XXX
			temphist->SetFillColor(ncflavourcolourwheel.GetNextColour());
			temphist->SetFillStyle(3002); // semi dense spots
			ncflavourstack->Add(temphist);
			
			// while looping over flavours, we'll also make a stack for each flavour
			// broken down by it's parent type / decay mode
			THStack* parentstack = new THStack(TString::Format("%sParents%d",nustring.c_str(),fluxver),TString::Format("%s Flux by Parent Decay Mode",nustring.c_str()));
			ColourWheel parentcolourwheel = ColourWheel();
			
			// loop over the decay modes, make a histo and add it to the stack
			std::map<int,std::string>* amap = GenerateDecayMap();
			for(auto&& apair : (*amap)){
				int decaymode = apair.first;
				std::string decaymodestring = apair.second;
				
				c1->Clear();
				treeout->Draw(TString::Format("NeutrinoEnergy>>h%sFlux_%s%d(100,0.,5.)",nustring.c_str(),decaymodestring.c_str(),fluxver), TString::Format("fluxver==%d&&NeutrinoPDG==%d&&ParentDecayMode==%d",fluxver,nupdg,decaymode));
				TH1F* temphist = (TH1F*)gROOT->FindObject(TString::Format("h%sFlux_%s%d",nustring.c_str(),decaymodestring.c_str(),fluxver));
				temphist->Scale(fluxscaling); // XXX XXX XXX XXX XXX XXX XXX per POT per tonne per whatever
				temphist->SetFillColor(parentcolourwheel.GetNextColour());
				temphist->SetFillStyle(3002); // semi dense spots
				parentstack->Add(temphist);
				
			} // end loop over decay modes
			
			// save the stack (of decay types for this flavour)
			c1->Clear();
			parentstack->Draw(); // must draw before we can get the axes!
			parentstack->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
			parentstack->GetYaxis()->SetTitle(TString::Format("%s Flux",nustring.c_str()));
			//CenterTitles(parentstack);
			parentstack->SetTitle("");
			leg = c1->BuildLegend();
			leg->SetFillStyle(0);
			//Preliminary("Simulation"); // must get called AFTER BuildLegend, else it gets added to it!
			//gPad->Modified(); gPad->Update(); // make sure it's really (re)drawn
			//c1.SaveAs(TString::Format("%s/reconstructed_neutrino_energy.png",outdir));
			parentstack->Write();
			delete parentstack;
			
		}  // end loop over nu flavours
		
		// save the stack (of neutrino flavours)
		flavourstack->Draw(); //"nostack"
		flavourstack->GetXaxis()->SetTitle("Neutrino Energy [MeV]");
		flavourstack->GetYaxis()->SetTitle("Neutrino Flux");
		//CenterTitles(flavourstack);
		flavourstack->SetTitle("");
		leg = c1->BuildLegend();
		leg->SetFillStyle(0);
		//Preliminary("Simulation");
		//gPad->Modified(); gPad->Update(); // make sure it's really (re)drawn
		//c1.SaveAs(TString::Format("%s/reconstructed_neutrino_energy.png",outdir));
		flavourstack->Write();
		delete flavourstack;
		
		// save the cc version
		ccflavourstack->Draw(); //"nostack"
		ccflavourstack->GetXaxis()->SetTitle("Neutrino Energy [MeV]");
		ccflavourstack->GetYaxis()->SetTitle("CC Neutrino Flux");
		//CenterTitles(ccflavourstack);
		ccflavourstack->SetTitle("");
		leg = c1->BuildLegend();
		leg->SetFillStyle(0);
		//Preliminary("Simulation");
		//gPad->Modified(); gPad->Update(); // make sure it's really (re)drawn
		//c1.SaveAs(TString::Format("%s/reconstructed_neutrino_energy.png",outdir));
		ccflavourstack->Write();
		delete ccflavourstack;
		
		// save the nc version
		ncflavourstack->Draw(); //"nostack"
		ncflavourstack->GetXaxis()->SetTitle("Neutrino Energy [MeV]");
		ncflavourstack->GetYaxis()->SetTitle("NC Neutrino Flux");
		//CenterTitles(ncflavourstack);
		ncflavourstack->SetTitle("");
		leg = c1->BuildLegend();
		leg->SetFillStyle(0);
		//Preliminary("Simulation");
		//gPad->Modified(); gPad->Update(); // make sure it's really (re)drawn
		//c1.SaveAs(TString::Format("%s/reconstructed_neutrino_energy.png",outdir));
		ncflavourstack->Write();
		delete ncflavourstack;
	}
	
//		TODO tables:
//		============
//		average multiplicity, momentum and angle of each particle type in primary proton-target intx
//		average multiplicity, momentum and angle of each particle type in neutrino-detector intx
//		flux per POT per tonne of (anti)neutrino flavours
//		% of total above from decay channels
	
	// make a 2D profile of flux density at the tank (integrated over z)
	for(fluxver=0; fluxver<2; fluxver++){
		c1->Clear();
		treeout->Draw(TString::Format("NuIntxVtx_X:NuIntxVtx_Y>>hNuFluxDensity%d",fluxver), TString::Format("(fluxver==%d)&&(NuVtxInTank==1)",fluxver),"colz");
		TH2F* fluxprofile = (TH2F*)gROOT->FindObject(TString::Format("hNuFluxDensity%d",fluxver));
		fluxprofile->Write();
		
		// also make versions broken down by flavour
		for(auto nutype : nutypes){
			std::string nustring = nutype.first;
			int nupdg = nutype.second;
			c1->Clear();
			treeout->Draw(TString::Format("NuIntxVtx_X:NuIntxVtx_Y>>h%sFluxDensity%d",nustring.c_str(),fluxver), TString::Format("(fluxver==%d)&&(NuVtxInTank==1)&&(NeutrinoPDG==%d)",fluxver,nupdg),"colz");
			fluxprofile = (TH2F*)gROOT->FindObject(TString::Format("h%sFluxDensity%d",nustring.c_str(),fluxver));
			fluxprofile->Write();
			cout<<"writing "<<fluxprofile->GetName()<<" to file (fluxver="<<fluxver<<")"<<endl;
		}
	}
	
	// FIXME:
	// parentprodmom_x and y need fixing for fluxver 0, fluxstage 0. Comparable in fluxstage 1 so its ok!
	// histogram colours are not saved neither in canvases nor histograms. why not?
	// comparison canvases for types: remove bins that have no entries in either
	// TODO:
	// add colours to the pie charts (although, currently colours aren't saved to histos....)
	
	cout<<"cleaning up"<<endl;
	oldflux->ResetBranchAddresses();
	delete oldflux;
	newflux->ResetBranchAddresses();
	delete newflux;
	
	fileout->Close();
	c1->Clear();
	cout<<"done, goodbye"<<endl;
	return 0;
}


// type conversion functions:
std::string GnumiToString(int code){
	if(gnumicodetoname.size()==0) GenerateGnumiMap();
	if(gnumicodetoname.count(code)!=0){
		return gnumicodetoname.at(code);
	} else {
		cerr<<"unknown gnumi code "<<code<<endl;
		return std::to_string(code);
	}
}

std::string PdgToString(int code){
	if(pdgcodetoname.size()==0) GeneratePdgMap();
	if(pdgcodetoname.count(code)!=0){
		return pdgcodetoname.at(code);
	} else {
		cerr<<"unknown pdg code "<<code<<endl;
		return std::to_string(code);
	}
}

std::string DecayTypeToString(int code){
	if(decaymap.size()==0) GenerateDecayMap();
	if(decaymap.count(code)!=0){
		return decaymap.at(code);
	} else {
		cerr<<"unknown decay code "<<code<<endl;
		return std::to_string(code);
	}
}

std::string MediumToString(int code){
	return std::to_string(code); // TODO fill this out
}

std::map<int,std::string>* GenerateGnumiMap(){
	if(gnumicodetoname.size()!=0) return &gnumicodetoname;
	gnumicodetoname.emplace(14 ,"Proton");
	gnumicodetoname.emplace(15 ,"Anti Proton");
	gnumicodetoname.emplace(3 ,"Electron");
	gnumicodetoname.emplace(2 ,"Positron");
	gnumicodetoname.emplace(53 ,"Electron Neutrino");
	gnumicodetoname.emplace(52 ,"Anti Electron Neutrino");
	gnumicodetoname.emplace(1 ,"Photon");
	gnumicodetoname.emplace(13 ,"Neutron");
	gnumicodetoname.emplace(25 ,"Anti Neutron");
	gnumicodetoname.emplace(5 ,"Muon+");
	gnumicodetoname.emplace(6 ,"Muon-");
	gnumicodetoname.emplace(10 ,"Kaonlong");
	gnumicodetoname.emplace(8 ,"Pion+");
	gnumicodetoname.emplace(9 ,"Pion-");
	gnumicodetoname.emplace(11 ,"Kaon+");
	gnumicodetoname.emplace(12 ,"Kaon-");
	gnumicodetoname.emplace(18 ,"Lambda");
	gnumicodetoname.emplace(26 ,"Antilambda");
	gnumicodetoname.emplace(16 ,"Kaonshort");
	gnumicodetoname.emplace(21 ,"Sigma-");
	gnumicodetoname.emplace(19 ,"Sigma+");
	gnumicodetoname.emplace(20 ,"Sigma0");
	gnumicodetoname.emplace(7 ,"Pion0");
	gnumicodetoname.emplace(99,"Kaon0");  // gnumi particle code for Kaon0 and Antikaon0
	gnumicodetoname.emplace(98,"Antikaon0");  // are both listed as "10 & 16" ... 
	gnumicodetoname.emplace(56 ,"Muon Neutrino");
	gnumicodetoname.emplace(55 ,"Anti Muon Neutrino");
	gnumicodetoname.emplace(27 ,"Anti Sigma-");
	gnumicodetoname.emplace(28 ,"Anti Sigma0");
	gnumicodetoname.emplace(29 ,"Anti Sigma+");
	gnumicodetoname.emplace(22 ,"Xsi0");
	gnumicodetoname.emplace(30 ,"Anti Xsi0");
	gnumicodetoname.emplace(23 ,"Xsi-");
	gnumicodetoname.emplace(31 ,"Xsi+");
	gnumicodetoname.emplace(24 ,"Omega-");
	gnumicodetoname.emplace(32 ,"Omega+");
	gnumicodetoname.emplace(33 ,"Tau+");
	gnumicodetoname.emplace(34 ,"Tau-");
	return &gnumicodetoname;
}

std::map<int,std::string>* GeneratePdgMap(){
	if(pdgcodetoname.size()!=0) return &pdgcodetoname;
	pdgcodetoname.emplace(2212,"Proton");
	pdgcodetoname.emplace(-2212,"Anti Proton");
	pdgcodetoname.emplace(11,"Electron");
	pdgcodetoname.emplace(-11,"Positron");
	pdgcodetoname.emplace(12,"Electron Neutrino");
	pdgcodetoname.emplace(-12,"Anti Electron Neutrino");
	pdgcodetoname.emplace(22,"Photon");
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
	return &pdgcodetoname;
}

std::map<int,std::string>* GenerateDecayMap(){
	if(decaymap.size()!=0) return &decaymap;
	decaymap.emplace(1,"K0L -> nue, pi-, e+");  //→
	decaymap.emplace(2,"K0L -> nuebar, pi+, e-");
	decaymap.emplace(3,"K0L -> numu, pi-, mu+");
	decaymap.emplace(4,"K0L -> numubar, pi+, mu-");
	decaymap.emplace(5,"K+  -> numu, mu+");
	decaymap.emplace(6,"K+  -> nue, pi0, e+");
	decaymap.emplace(7,"K+  -> numu, pi0, mu+");
	decaymap.emplace(8,"K-  -> numubar, mu-");
	decaymap.emplace(9,"K-  -> nuebar, pi0, e-");
	decaymap.emplace(10,"K-  -> numubar, pi0, mu-");
	decaymap.emplace(11,"mu+ -> numubar, nue, e+");
	decaymap.emplace(12,"mu- -> numu, nuebar, e-");
	decaymap.emplace(13,"pi+ -> numu, mu+");
	decaymap.emplace(14,"pi- -> numubar, mu-");
	return &decaymap;
}

// Produce pie chart of particles that produced neutrinos
// =======================================================
TPie* GeneratePieFromHisto(std::string histoname, int verbose){
	TH1F* histo = (TH1F*)gROOT->FindObject(histoname.c_str());
	if(histo==nullptr) cerr<<"GeneratePieFromHisto could not find histo "<<histoname<<endl;
	TPie* thepie = GeneratePieFromHisto(histo, verbose);
	return thepie;
}

TPie* GeneratePieFromHisto(TH1F* histo, int verbose){
	std::string histoname = std::string(histo->GetName());
	if(verbose) cout<<"creating pie chart from histo "<<histoname<<", which has "
									<<histo->GetNbinsX()<<" bins with contents: "<<endl;
	std::vector< std::pair<std::string,float> > histbins;
	for(int bini=0; bini<histo->GetNbinsX(); bini++){
		TString binlabel = histo->GetXaxis()->GetBinLabel(bini+1);
		double binconts = histo->GetBinContent(bini+1);
		if(binconts<0.01) binconts = 0.0f;  // round floats. useful if the histo has been scaled.
		if(verbose && binconts!=0.0f) cout<<binlabel.Data()<<" : "<<binconts<<endl;
		if(binconts<0) cerr<<"error converting "<<histoname<<" to pie chart: bin "<<binlabel.Data()
			<<" has "<<binconts<<" entries!"<<endl;
		if(binconts!=0) histbins.emplace_back(binlabel.Data(),binconts);
	}
	
	TPie* thepie = new TPie(TString::Format("%sPie",histoname.c_str()), TString::Format("%s Pie Chart",histoname.c_str()), histbins.size());
	
	for(int bini=0; bini<histbins.size(); bini++){
		std::pair<std::string,float> abin = histbins.at(bini);
		std::string thebinlabel = abin.first;
		float thebincontents = abin.second;
		if(thebincontents<0) cerr<<"error converting "<<histoname<<" to pie chart: bin "<<thebinlabel
			<<" has "<<thebincontents<<" entries!"<<endl;
		thepie->SetEntryVal(bini,thebincontents);  // NO +1 - TPie's have no underflow bin!
		thepie->SetEntryLabel(bini,thebinlabel.c_str());
		
		//TPieSlice pieslice = thepie->GetSlice(bini+1);
		//pieslice->SetValue(thebincontents);
		//pieslice->SetName(thebinlabel.c_str()); // doesn't seem to work?
	}
	
	// tuning to try to make it look less shitty
	//thepie->SetAngularOffset(333); 
	//thepie->SetLabelFormat("#splitline{%txt}{#splitline{%val}{(%perc)}}");
	//thepie->SetValueFormat("%4.0f");
	//thepie->SetPercentFormat("%3.0f");
	//thepie->SetCircle(0.5, 0.4702026, 0.3302274);
	//thepie->SetTextSize(0.03455766);
	//thepie->SetLabelsOffset(0.05);	// shift labels from pie centre so they don't mash together.
	
	//TCanvas* pieCanv = new TCanvas("pieCanv","Pie Charts",385,110,700,867);
	//pieCanv->Divide(1,2);
	//thepie->Draw("rsc");
/*
  "R" Print the labels along the central "R"adius of slices.
  "T" Print the label in a direction "T"angent to circle that describes the TPie.
  "SC" Paint the the labels with the "S"ame "C"olor as the slices.
  "3D" Draw the pie-chart with a pseudo 3D effect.
  "NOL" No OutLine: Don't draw the slices' outlines, any property over the slices' line is ignored.
  ">" Sort the slices in increasing order.
  "<" Sort the slices in decreasing order.
*/
	return thepie;
}

#if FLUX_STAGE==1
void GetGenieEntryInfo(genie::EventRecord* gevtRec, genie::Interaction* genieint, GenieInfo &thegenieinfo, bool printneutrinoevent){
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
	/*Int_t*/ //thegenieinfo.neutinteractioncode = genie::utils::ghep::NeutReactionCode(gevtRec);
	/*Int_t*/ thegenieinfo.nuanceinteractioncode  = genie::utils::ghep::NuanceReactionCode(gevtRec);
	/*TLorentzVector**/ thegenieinfo.IntxVtx = gevtRec->Vertex();
	/*Double_t*/ thegenieinfo.Intx_x = thegenieinfo.IntxVtx->X() * 100.;   // same info as nuvtx in g4dirt file
	/*Double_t*/ thegenieinfo.Intx_y = thegenieinfo.IntxVtx->Y() * 100.;   // GENIE uses meters
	/*Double_t*/ thegenieinfo.Intx_z = thegenieinfo.IntxVtx->Z() * 100.;   // GENIE uses meters
	/*Double_t*/ thegenieinfo.Intx_t = thegenieinfo.IntxVtx->T() * 1000000000; // GENIE uses seconds
	
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
//	/*Double_t*/ thegenieinfo.Q2 = genieint->Kine().Q2();  // not set in our GENIE files!
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
		cout<<"This was a "<< thegenieinfo.procinfostring <<" (neut code "<<thegenieinfo.neutinteractioncode
			<<") interaction of a "
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
		cout<< " N(p) = "   << thegenieinfo.numfsprotons
			<< " N(n) = "   << thegenieinfo.numfsneutrons
			<< endl
			<< " N(pi^0) = "  << thegenieinfo.numfspi0
			<< " N(pi^+) = "  << thegenieinfo.numfspiplus
			<< " N(pi^-) = "  << thegenieinfo.numfspiminus
			<<endl;
	}
}
#endif // FLUX_STAGE==1
