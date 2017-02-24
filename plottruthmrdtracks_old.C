/* vim:set noexpandtab tabstop=4 wrap */

#ifndef PLOTVERBOSE
//#define PLOTVERBOSE
#endif
#ifndef WCSIMDEBUG
//#define WCSIMDEBUG
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
#include "TLorentzVector.h"
#include "TLegend.h"
#include <regex>
#include <exception>	// for stdexcept
#include <vector>
#include <map>
#include <string>
#include <algorithm>	// remove and remove_if
#include <iostream>
#include <iomanip>
#include <fstream> 		//std::ofstream
#include <stdlib.h> /* atoi */

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

const Float_t tank_start = 15.70;          // front face of the tank in cm
const Float_t tank_radius = 152.4;         // tank radius in cm
/* from WCSimDetectorConfigs.cc
tankouterRadius= 1.524*m;
tankzoffset = 15.70*cm;
*/

const char* dirtpath="/pnfs/annie/persistent/users/moflaher/g4dirt";
const char* geniepath="/pnfs/annie/persistent/users/rhatcher/genie";
const char* wcsimpath="/pnfs/annie/persistent/users/moflaher/wcsim";
const char* wcsimlibrarypath="/annie/app/users/moflaher/wcsim/wcsim/libWCSimRoot.so";
const char* outpath="/annie/app/users/moflaher/wcsim/root_work";

const Bool_t printneutrinoevent=false;

void truthtracks(){
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
	TFile* wcsimfile;
	TString wcsimfilepath;
	TTree* wcsimT;
	Int_t numwcsimentries;
	
	TFile* dirtfile;
	TTree* tankflux;
	TTree* tankmeta;
	
	TFile* geniefile;
	TString geniefilepath;
	TTree* gtree;
	Int_t numgenietentries;
	
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
	Double_t totalpots=0;	// count of POTs in all processed files
	
	// gtree
	cout<<"creating new genie::NtpMCEventRecord"<<endl;
	genie::NtpMCEventRecord* genierecordval = new genie::NtpMCEventRecord;
	TBranch* genierecordBranch=0;
	
	// wcsimT
	WCSimRootEvent* b=0, *m=0, *v=0;
	TBranch* bp=0, *mp=0, *vp=0;
	WCSimRootTrigger* atrigt=0, *atrigm=0, *atrigv=0;
	
	// information from genie:
	gROOT->cd();
	TH1D* incidentneutrinoenergiesall = new TH1D("incidentneutrinoenergiesall","Distribution of Probe Neutrino Energies:Energy (GeV):Num Events",100,0,0.);
	TH1D* incidentneutrinoenergiesaccepted = new TH1D("incidentneutrinoenergiesaccepted","Distribution of Accepted Probe Neutrino Energies:Energy (GeV):Num Events",100,0,0.);
	TH1D* fslanglesall = new TH1D("fslanglesall","Distribution of Final State Lepton Angles:Angle (rads):Num Events",100,0.,TMath::Pi());
	TH1D* fslenergiesall = new TH1D("fslenergiesall","Distribution of Final State Lepton Energies (GeV):Energy (GeV):Num Events",100,0.,3.);
	TH1D* fslanglesaccepted = new TH1D("fslanglesaccepted","Distribution of Accepted Final State Lepton Angles:Angle (rads):Num Events",100,0.,TMath::Pi());
	TH1D* fslenergiesaccepted = new TH1D("fslenergiesaccepted","Distribution of Accepted Final State Lepton Energies (GeV):Energy (GeV):Num Events",100,0.,3.);
	
	// information from WCSim: 
	TH1D* fslanglesacceptedwcsim = new TH1D("fslanglesacceptedwcsim","Distribution of Accepted Final State Lepton Angles:Angle (rads):Num Events",100,0.,TMath::Pi());
	TH1D* fslenergiesacceptedwcsim = new TH1D("fslenergiesacceptedwcsim","Distribution of Accepted Final State Lepton Energies (GeV):Energy (GeV):Num Events",100,0.,3.);
	
//	// debugging:
	TH1D* fsltruetracklength = new TH1D("fsltruetracklength", "Distribution of True Track Lengths", 100, 0., 1500.);
#ifdef WCSIMDEBUG
	TH3D* tankstartvertices = new TH3D("tankstartvertices", "Distribution of Tank Starting Vertices", 100, -150.,150., 100, -150., 150., 100, -100., 800.);
	TH3D* vetostartvertices = new TH3D("vetostartvertices", "Distribution of Veto Starting Vertices", 100, -150.,150., 100, -150., 150., 100, -100., 800.);
	TH3D* mrdstartvertices = new TH3D("mrdstartvertices", "Distribution of MRD Starting Vertices", 100, -150.,150., 100, -150., 150., 100, -100., 800.);
	TH3D* tankstopvertices = new TH3D("tankstopvertices", "Distribution of Tank Stopping Vertices", 100, -150.,150., 100, -150., 150., 100, -100., 800.);
	TH3D* vetostopvertices = new TH3D("vetostopvertices", "Distribution of Veto Stopping Vertices", 100, -150.,150., 100, -150., 150., 100, -100., 800.);
	TH3D* mrdstopvertices = new TH3D("mrdstopvertices", "Distribution of MRD Stopping Vertices", 100, -150.,150., 100, -150., 150., 100, -100., 800.);
#endif
	
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
	for(Int_t inputEntry=0; inputEntry<numents; inputEntry++){
		/* 	1. Load next g4dirt entry */ 
		Long64_t localEntry = c->LoadTree(inputEntry);
		if( localEntry<0){ cout<<"end of tchain"<<endl; break; }
		Int_t nextTreeNumber = c->GetTreeNumber();
		if(treeNumber!=nextTreeNumber){
			cout<<"new tree: "<<nextTreeNumber<<endl;
			treeNumber=nextTreeNumber;
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
			geniefile = TFile::Open(geniefilepath);
			if(!geniefile){cout<<"this genie file doesn't exist!"<<endl; break; }
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
			wcsimfile = TFile::Open(wcsimfilepath);
			if(!wcsimfile){cout<<"this wcsimfile doesn't exist!"<<endl; break; }
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
			if(totalpots<lasttotpots){ cout<<"ERROR! TOTAL POTS CAME DOWN FROM "<<lasttotpots<<" TO "<<totalpots<<endl; return; }
			cout<<"adding "<<pots<<" POTs to the running total, making "<<totalpots<<" POTs so far"<<endl;
			
			// wcsimT:
			// wcsim trigger classes
			wcsimT->SetBranchAddress("wcsimrootevent",&b, &bp);
			wcsimT->SetBranchAddress("wcsimrootevent_mrd",&m, &mp);
			wcsimT->SetBranchAddress("wcsimrootevent_facc",&v, &vp);
			if(bp==0||mp==0||vp==0){ cout<<"branches are zombies!"<<endl; }
			
			// gtree:
			gtree->SetBranchAddress("gmcrec",&genierecordval,&genierecordBranch);
		}
#ifdef PLOTVERBOSE
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
#ifdef PLOTVERBOSE
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
//		TLorentzVector* genieVtx = gevtRec->Vertex();
//		G4double x = genieVtx->X() * m;         // same info as nuvtx in g4dirt file
//		G4double y = genieVtx->Y() * m;         // GENIE uses meters
//		G4double z = genieVtx->Z() * m;         // GENIE uses meters
//		G4double t = genieVtx->T() * second;    // GENIE uses seconds for time
		
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
		TParticlePDG * targetnucleus = genie::PDGLibrary::Instance()->Find( targetnucleuspdg );
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
		
#ifdef PLOTVERBOSE
		cout<<"filling incident neutrino histograms"<<endl;
#endif
		incidentneutrinoenergiesall->Fill(probeenergy);
		//incidentneutrinoanglesall->Fill(probeangle);
		fslanglesall->Fill(fslanglegenie);
		fslenergiesall->Fill(fsleptonenergy);
		
		/* 4. Check if interaction is QE - if not, continue */
		if (!(isQE && isWeakCC)){ continue; }
		
		/*5. If so, load wcsim detector response. */
		// read only first subtrigger; delayed decay detector response is not interesting for primary FSL tracks
#ifdef PLOTVERBOSE
		cout<<"getting wcsim entry"<<endl;
#endif
		if(localEntry>numwcsimentries){ cout<<"can't load localEntry "<<localEntry
				<<"from "<<wcsimfilepath<<" wcsimT - not enough entries!"<<endl; continue; }
		wcsimT->GetEntry(localEntry);
		atrigt = b->GetTrigger(0);
		atrigm = m->GetTrigger(0);
		atrigv = v->GetTrigger(0);
		
		Int_t numtracks = atrigt->GetNtrack();
#ifdef PLOTVERBOSE
		cout<<"wcsim event had "<<numtracks<<" truth tracks"<<endl;
#endif
		
		std::vector<Float_t> primaryenergiesvector;
		std::vector<Double_t> scatteringanglesvector;
		std::vector<Int_t> acceptedtrackids;
		
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
#ifdef PLOTVERBOSE
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
			
#ifdef PLOTVERBOSE
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
#ifdef PLOTVERBOSE
			cout<<"primarystopvol is "<<primarystopvol<<endl;
#endif
			
			TVector3 primarystartvertex(nextrack->GetStart(0),nextrack->GetStart(1),nextrack->GetStart(2));
			TVector3 primarystopvertex(nextrack->GetStop(0), nextrack->GetStop(1), nextrack->GetStop(2));
			TVector3 differencevector  = (primarystopvertex-primarystartvertex);
			
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
			
			Float_t opp = primarystopvertex.X() - primarystartvertex.X();
			Float_t adj = primarystopvertex.Z() - primarystartvertex.Z();
			Float_t avgtrackanglex = TMath::ATan(opp/adj);
			
			Float_t xatmrd = primarystartvertex.X() + (MRD_layer2-primarystartvertex.Z())*TMath::Tan(avgtrackanglex);
#ifdef WCSIMDEBUG
			cout<<"primary start vertex: ("<<primarystartvertex.X()<<", "<<primarystartvertex.Y()
				<<","<<primarystartvertex.Z()<<")"<<endl;
			cout<<"primary stop vertex: ("<<primarystopvertex.X()<<", "<<primarystopvertex.Y()
				<<","<<primarystopvertex.Z()<<")"<<endl;
			cout<<"opp = "<<opp<<endl;
			cout<<"adj = "<<adj<<endl;
			cout<<"angle = "<<avgtrackanglex<<endl;
			cout<<"tan(angle) = "<<TMath::Tan(avgtrackanglex)<<endl;
			cout<<"projected x at z="<<MRD_layer2<<" is "<<xatmrd<<endl;
			fsltruetracklength->Fill(differencevector.Mag());
#endif
			
			opp = primarystopvertex.Y() - primarystartvertex.Y();
			Float_t avgtrackangley = TMath::ATan(opp/adj);
			Float_t yatmrd = primarystartvertex.Y() + (MRD_layer2-primarystartvertex.Z())*TMath::Tan(avgtrackangley);
#ifdef WCSIMDEBUG
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
			
			// Add this primary information to a vector. We _should_ only find one suitable primary
			// and can break once it is found, but... let's just see. Keep wcsim trackID in case.
			Int_t primarytrackid = nextrack->GetId();
			
#ifdef PLOTVERBOSE
			cout<<"found a suitable primary track"<<endl;
#endif
			scatteringanglesvector.push_back(scatteringangle);
			primaryenergiesvector.push_back(primaryenergy);
			acceptedtrackids.push_back(primarytrackid);
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
			fslanglesacceptedwcsim->Fill(scatteringanglesvector.at(indextouse));
			fslenergiesacceptedwcsim->Fill(primaryenergiesvector.at(indextouse));
		} else if (scatteringanglesvector.size()==1) { // just one matched wcsim track. 
			fslanglesacceptedwcsim->Fill(scatteringanglesvector.at(0));
			fslenergiesacceptedwcsim->Fill(primaryenergiesvector.at(0));
			
			cout<<"wcsim lepton angle = "<<scatteringanglesvector.at(0)<<endl;
			cout<<"wcsim lepton energy = "<<primaryenergiesvector.at(0)<<endl;
			cout<<"wcsim lepton trackid = "<<acceptedtrackids.at(0)<<endl;
			cout<<"compare to:"<<endl;
			cout<<"genie lepton angle = "<<fslanglegenie<<endl;
			cout<<"genie lepton energy = "<<fsleptonenergy<<endl;
		} else {
			cout<<"no accepted wcsim track"<<endl;
		}
		if(scatteringanglesvector.size()!=0){
			cout<<"MATCHED A TRACK YEEEEAH"<<endl; /*return;*/
			// by getting this far we've found a primary muon that penetrated the MRD.
			// fill genie 'accepted' histograms.
			incidentneutrinoenergiesaccepted->Fill(probeenergy);
			//incidentneutrinoanglesaccepted->Fill(probeangle);
			fslanglesaccepted->Fill(fslanglegenie);
			fslenergiesaccepted->Fill(fsleptonenergy);
		}
		
	}
	
	cout<<"generating scaled histograms"<<endl;
	Double_t numbeamspills = totalpots/(4.0 * TMath::Power(10.,12.));
	Double_t numbeamspillsperday = (24.*60.*60.*1000.)/133.3333;	// 24 hours in ms / 133.33 ms between spills
	Double_t numdays = numbeamspills/numbeamspillsperday;
	cout<<"Results based on "<<totalpots<<" POTs, or "<<numbeamspills<<" beam spills, or "<<numdays<<" days of data"<<endl;
	
	gROOT->cd();
	// create scaled histograms with bin contents scaled to the number of POTs in input files:
	TH1D* incidentneutrinoenergiesallscaled = new TH1D(*incidentneutrinoenergiesall);
	TH1D* incidentneutrinoenergiesacceptedscaled = new TH1D(*incidentneutrinoenergiesaccepted);
	TH1D* fslanglesallscaled = new TH1D(*fslanglesall);
	TH1D* fslenergiesallscaled = new TH1D(*fslenergiesall);
	TH1D* fslanglesacceptedscaled = new TH1D(*fslanglesaccepted);
	TH1D* fslenergiesacceptedscaled = new TH1D(*fslenergiesaccepted);
	TH1D* fslanglesacceptedwcsimscaled = new TH1D(*fslanglesacceptedwcsim);
	TH1D* fslenergiesacceptedwcsimscaled = new TH1D(*fslenergiesacceptedwcsim);
	
	TH1D* thehistopointersarr[]{
	incidentneutrinoenergiesallscaled,
	incidentneutrinoenergiesacceptedscaled,
	fslanglesallscaled,
	fslanglesacceptedscaled,
	fslanglesacceptedwcsimscaled,
	fslenergiesallscaled,
	fslenergiesacceptedscaled,
	fslenergiesacceptedwcsimscaled};
	//TODO: modify size if you add more histograms!!!
	std::vector<TH1D*> scaledhistopointers(thehistopointersarr,thehistopointersarr+8);
	
	for(int i=0; i<scaledhistopointers.size(); i++){
		TH1D* temp = scaledhistopointers.at(i);
		for(int j=0; j<temp->GetNbinsX(); j++){
			Double_t thecontent = temp->GetBinContent(j);
			if(thecontent!=0){ thecontent /= numdays; }
			temp->SetBinContent(j,thecontent);
		}
	}
	
	TH1D* thehistopointersarr2[]{
	incidentneutrinoenergiesall,
	incidentneutrinoenergiesaccepted,
	fslanglesall,
	fslanglesaccepted,
	fslanglesacceptedwcsim,
	fslenergiesall,
	fslenergiesaccepted,
	fslenergiesacceptedwcsim};
	//TODO: modify size if you add more histograms!!!
	std::vector<TH1D*> histopointers(thehistopointersarr2,thehistopointersarr2+8);
	
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
	for(int i=0; i<histopointers.size(); i++){
		TLegend* leg;
		TCanvas* canv;
		TH1D* hist = histopointers.at(i);
		if(i==0){
			canv = new TCanvas(TString::Format("c%d",i),TString::Format("c%d",i));
			canvaspointers.push_back(canv);
			canv->cd();
			leg = new TLegend(0.5,0.78,0.7,0.88);
			legendpointers.push_back(leg);
			leg->AddEntry(hist,"incident","l");
			hist->SetLineColor(kRed);
			hist->Draw();
		} else if(i==1){
			hist->SetLineColor(kBlue);
			leg->AddEntry(hist,"accepted (truth)","l");
			hist->Draw("same");
			leg->Draw();
		} else if((i-2)%3==0){
			canv = new TCanvas(TString::Format("c%d",i),TString::Format("c%d",i));
			canvaspointers.push_back(canv);
			canv->cd();
			leg = new TLegend(0.5,0.78,0.7,0.88);
			legendpointers.push_back(leg);
			leg->AddEntry(hist,"incident","l");
			hist->SetLineColor(kRed);
			hist->Draw();
		} else if((i-2)%3==1){
			hist->SetLineColor(kBlue);
			leg->AddEntry(hist,"accepted (truth)","l");
			hist->Draw("same");
		} else {
			hist->SetLineColor(kViolet);
			hist->SetLineStyle(2);
			leg->AddEntry(hist,"accepted (reco)","l");
			hist->Draw("same");
			leg->Draw();
		}
	}
	
	// draw the scaled histograms
	// ==========================
	std::vector<TCanvas*> scaledcanvaspointers;
	std::vector<TLegend*> scaledlegendpointers;
	for(int i=0; i<scaledhistopointers.size(); i++){
		TLegend* leg;
		TCanvas* canv;
		TH1D* hist = scaledhistopointers.at(i);
		if(i==0){
			canv = new TCanvas(TString::Format("cc%d",i),TString::Format("cc%d",i));
			scaledcanvaspointers.push_back(canv);
			canv->cd();
			leg = new TLegend(0.5,0.78,0.7,0.88);
			scaledlegendpointers.push_back(leg);
			leg->AddEntry(hist,"incident","l");
			hist->SetLineColor(kRed);
			hist->Draw();
		} else if(i==1){
			hist->SetLineColor(kBlue);
			leg->AddEntry(hist,"accepted (truth)","l");
			hist->Draw("same");
			leg->Draw();
		} else if((i-2)%3==0){
			canv = new TCanvas(TString::Format("cc%d",i),TString::Format("cc%d",i));
			scaledcanvaspointers.push_back(canv);
			canv->cd();
			leg = new TLegend(0.5,0.78,0.7,0.88);
			scaledlegendpointers.push_back(leg);
			leg->AddEntry(hist,"incident","l");
			hist->SetLineColor(kRed);
			hist->Draw();
		} else if((i-2)%3==1){
			hist->SetLineColor(kBlue);
			leg->AddEntry(hist,"accepted (truth)","l");
			hist->Draw("same");
		} else {
			hist->SetLineColor(kViolet);
			hist->SetLineStyle(2);
			leg->AddEntry(hist,"accepted (reco)","l");
			hist->Draw("same");
			leg->Draw();
		}
	}
	
//	 draw other canvases
	TCanvas* debug1 = new TCanvas("debug1","debug1");
	debug1->cd();
	fsltruetracklength->Draw();
	debug1->SaveAs("muon track lengths.png");
	if(fsltruetracklength) delete fsltruetracklength; fsltruetracklength=0;
	if(debug1) delete debug1; debug1=0;
	
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
	for(int i=0; i<histopointers.size(); i++){
		TH1D* temp = histopointers.at(i);
		delete temp;
	}
	
	// delete scaled histograms
	//cout<<"deleting scaled histograms"<<endl;
	for(int i=0; i<scaledhistopointers.size(); i++){
		TH1D* temp = scaledhistopointers.at(i);
		delete temp;
	}
	cout<<"end"<<endl;
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
