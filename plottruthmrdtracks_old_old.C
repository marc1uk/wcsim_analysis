/* vim:set noexpandtab tabstop=4 wrap */
#include <stdlib.h> /* atoi */
#include <regex>

// GENIE headers
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepUtils.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Interaction/Interaction.h"

const Float_t MRD_width = (305./2.);      // half width of steel in cm
const Float_t MRD_height = (274./2.);     // half height of steel in cm
const Float_t MRD_layer2 = 2789.45;       // position in wcsim coords of second scint layer in cm
const Float_t MRD_start = 3255.;          // position in wcsim coord of MRD front face in cm
const Float_t MRD_depth = 139.09;         // total depth of the MRD in cm
/* output from WCSim:
########## MRD front face: 3255                       ##########
########## MRD total Z length: 139.09                 ########## 
########## MRD scintillator layer 0  (H) at z=2668.35 ##########
########## MRD scintillator layer 1  (V) at z=2789.45 ##########
########## MRD scintillator layer 2  (H) at z=2910.55 ##########
########## MRD scintillator layer 3  (V) at z=3031.65 ##########
########## MRD scintillator layer 4  (H) at z=3152.75 ##########
########## MRD scintillator layer 5  (V) at z=3273.85 ##########
########## MRD scintillator layer 6  (H) at z=3394.95 ##########
########## MRD scintillator layer 7  (V) at z=3516.05 ##########
########## MRD scintillator layer 8  (H) at z=3637.15 ##########
########## MRD scintillator layer 9  (V) at z=3758.25 ##########
########## MRD scintillator layer 10 (H) at z=3879.35 ########## */

const char* dirtpath="/home/marc/LinuxSystemFiles/WChSandBox/build/fluxesandtables";
const char* geniepath="/home/marc/LinuxSystemFiles/WChSandBox/build/fluxesandtables";
const char* wcsimpath="/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in";
const char* wcsimlibrarypath="/home/marc/LinuxSystemFiles/WCSim/gitver/wcsim/libWCSimRoot.so";

const Bool_t printneutrinoevent=true;

void truthtracks(){
	// load WCSim library for reading WCSim files
	gSystem->Load(wcsimlibrarypath);
	// load genie scripts for reading genie files
	TString script_dir = gSystem->Getenv("GENIE");
	script_dir += "/src/scripts/gcint/";
	TString curr_dir = gSystem->pwd();
	gSystem->cd(script_dir.Data());
	gROOT->ProcessLine(".x loadincs.C");
	gROOT->ProcessLine(".x loadlibs.C");
	gSystem->cd(curr_dir.Data());

	// Declare paths, files, trees...
	TFile* wcsimfile;
	TTree* wcsimT;

	TFile* dirtfile;
	TTree* tankflux;
	TTree* tankmeta;

	TFile* geniefile;
	TTree* gtree;

	// TChain for dirt files - this will be the main driver of the loop - all it's events will be processed.
	TChain* c =  new TChain("tankflux");
	c->Add(TString::Format("%s/annie_tank_flux.*.root");
	Int_t numents = c->GetEntries();

	// tankflux
	Int_t genieentry=0;
	Branch* genieentrybranch=0;
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
	Int_t pots;
	Int_t totalpots=0;	// count of POTs in all processed files

	// gtree
	genie::NtpMCEventRecord* genierecordval = new genie::NtpMCEventRecord;
	TBranch* genierecordBranch=0;
	geniedata->SetBranchAddress("gmcrec",&genierecordval,&genierecordBranch);

	// wcsimT
	WCSimRootEvent* b=0, *m=0, *v=0;
	TBranch* bp=0, *mp=0, *vp=0;
	WCSimRootTrigger* atrigt=0, *atrigm=0, *atrigv=0;

	// information from genie:
	TH1F* incidentneutrinoenergies = new TH1F("incidentneutrinoenergies","Distribution of Probe Neutrino Energies:Energy (GeV):Num Events",20,0,0.);
	TH1F* finalstateleptonanglesall = new TH1F("finalstateleptonanglesall","Distribution of Final State Lepton Angles:Angle (rads):Num Events",20,-TMath::Pi(),TMath::Pi());
	TH1F* finalstateleptonenergiesall = new TH1F("finalstateleptonenergiesall","Distribution of Final State Lepton Energies (GeV):Energy (GeV):Num Events",20,0.,3.);
	TH1F* finalstateleptonanglesaccepted = new TH1F("finalstateleptonanglesaccepted","Distribution of Accepted Final State Lepton Angles:Angle (rads):Num Events",20,-TMath::Pi(),TMath::Pi());
	TH1F* finalstateleptonenergiesaccepted = new TH1F("finalstateleptonenergiesaccepted","Distribution of Accepted Final State Lepton Energies (GeV):Energy (GeV):Num Events",20,0.,3.);

	// information from WCSim: 
	TH1F* finalstateleptonanglesacceptedwcsim = new TH1F("finalstateleptonanglesacceptedwcsim","Distribution of Accepted Final State Lepton Angles:Angle (rads):Num Events",20,-TMath::Pi(),TMath::Pi());
	TH1F* finalstateleptonenergiesacceptedwcsim = new TH1F("finalstateleptonenergiesacceptedwcsim","Distribution of Accepted Final State Lepton Energies (GeV):Energy (GeV):Num Events",20,0.,3.);

	c->LoadTree(0);
	Int_t treeNumber = c->GetTreeNumber();
	tankflux = c->GetTree();
	Int_t thistreesentries = tankflux->GetEntries();

	/*
	1. Load next g4dirt entry
	2. Check if genie primary, and volume is in tank - if not, continue
	3. If so, load genie entry.
	4. Check if interaction is QE - if not, continue
	5. If so, load wcsim detector response. 
	6. Load tracks, look for primary mu track through the MRD, record interaction details. 
	*/

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
			
			// retrieve genie filename from the meta tree, and open the corresponding genie file
			tankmeta = (TTree*)dirtfile->Get("tankmeta");
			tankmeta->SetBranchAddress("inputFluxName", geniefilename, &geniefilenamebranch);
			geniefilenamebranch->GetEntry(0);
			geniefile = TFile::Open(TString::Format("%s/%s",geniepath,geniefilename);
			gtree = (TTree*)geniefile->Get("gtree");
			
			// use regexp to pull out the file number needed for identifying the corresponding wcsim file
			std::match_results<string::const_iterator> submatches;
			// filename is of the form "annie_tank_flux.####.root"
			// #### is input file num. Need this to match against genie/wcsim file names
			std::regex theexpression ("[^0-9]+.([0-9]+).*");
			std::regex_match (nextdirtfilename, submatches, theexpression);
			std::string submatch = (std::string)submatches[0];	// match 0 is 'whole match' or smthg
			if(submatch==""){ cout<<"unrecognised input file pattern: "<<nextdirtfilename<<endl; }
			submatch = (std::string)submatches[1];
			int filenum = atoi(submatch.c_str());
			
			// use filenum to open the corresponding wcsim file
			wcsimfile = TFile::Open(TString::Format("%s/wcsim_0.%d.root",wcsimpath,filenum);
			wcsimT = (TTree*)wcsimfile->Get("wcsimT");
			
			/* Set the branch addresses for the new trees */
			// tankflux:
			// genie file entry number for each entry, to get the genie intx info
			c->SetBranchAddress("entry",&genieentry,&genieentrybranch);
			// number of primaries (so we can create appropriately sized array)
			c->SetBranchAddress("ntank",&ntankbranchval, nTankBranch);
			// material of vertex - identify as 'TankWater' to pull only primaries in the tank
			c->SetBranchAddress("vtxmat",&vertexmaterial, &vertexmaterialbranch);
			// array of whether particle is a genie primary
			nuprimaryBranch=c->GetBranch("primary");
			// tankmeta:
			// POTs in this genie file, so we can normalise wcsim event frequency
			tankmeta->SetBranchAddress("inputTotalPOTs", pots, &potsbranch);
			potsbranch->GetEntry(0);
			totalpots += pots;
			
			// wcsimT:
			// wcsim trigger classes
			wcsimT->SetBranchAddress("wcsimrootevent",&b, &bp);
			wcsimT->SetBranchAddress("wcsimrootevent_mrd",&m, &mp);
			wcsimT->SetBranchAddress("wcsimrootevent_facc",&v, &vp);
			if(bp==0||mp==0||vp==0){ cout<<"branches are zombies!"<<endl; }
			
			// gtree:
			gtree->SetBranchAddress("gmcrec",&genierecordval,&genierecordBranch);
		}
		cout<<"entry "<<inputEntry<<": localEntry "<<localEntry
		    <<" of "<<thistreesentries<<" in tree "<<treeNumber<<endl;
		
		/* 2. Check if genie primary, and volume is in tank - if not, continue */
		nTankBranch->GetEntry(localEntry);
		vertexmaterialbranch->GetEntry(localEntry);
		if(vertexmaterial!="TankWater"){ continue; }	// neutrino intx wasn't in tank water
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
		genieentrybranch->GetEntry(localEntry);
		genierecordBranch->GetEntry(genieentry);
		genie::EventRecord* gevtRec = genierecordval->event;
		genie::Interaction* genieint = gevtRec->Summary();
		
		// process information:
		TString procinfostring = genieint->ProcInfo().AsString();
		TString scatteringtypestring = genieint->ScatteringTypeAsString();
		TString interactiontypestring = genieint->InteractionTypeAsString();
		Bool_t isQE = interaction.ProcInfo().IsQuasiElastic();
		Double_t neutrinoq2 = genieint->Kine().Q2();
		TLorentzVector& k1 = *(gevtRec->Probe()->P4());
		TLorentzVector& k2 = *(gevtRec->FinalStatePrimaryLepton()->P4());
		Double_t costhfsl = TMath::Cos( k2.Vect().Angle(k1.Vect()) );
//		TLorentzVector* genieVtx = gevtRec->Vertex();
//		G4double x = genieVtx->X() * m;         // same info as nuvtx in g4dirt file
//		G4double y = genieVtx->Y() * m;         // GENIE uses meters
//		G4double z = genieVtx->Z() * m;         // GENIE uses meters
//		G4double t = genieVtx->T() * second;    // GENIE uses seconds for time
		Int_t neutinteractioncode = genie::utils::ghep::NeutReactionCode(gevtRec);
		Int_t nuanceinteractioncode  = genie::utils::ghep::NuanceReactionCode(gevtRec);
		
		// neutrino information:
		Double_t probeenergy = genieint->InitState().ProbeE(genie::kRfLab);	// GeV
		Int_t probepdg = genieint->InitState().Probe()->PdgCode();
		TSring probepartname = genieint->InitState().Probe()->GetName();
		TLorentzVector* probemomentum = gevtRec->Probe()->P4();
		TVector3 probethreemomentum = probemomentum->Vect();
		TVector3 probemomentumdir = probethreemomentum.Unit();
		Double_t probeanglex = TMath::ATan(probethreemomentum.X()/probethreemomentum.Z());
		Double_t probeangley = TMath::ATan(probethreemomentum.Y()/probethreemomentum.Z());
		Double_t probeangle = TMath::Max(probeanglex,probeanglexy);
		// n.b.  genieint->InitState().Probe != gevtRec->Probe()
		
		// target nucleon:
		int targetnucleonpdg = genieint->InitState().Tgt().HitNucPdg();
		TString targetnucleonname;
		if ( genie::pdg::IsNeutronOrProton(targetnucleonpdg) ) {
			TParticlePDG * p = genie::PDGLibrary::Instance()->Find(targetnucleonpdg);
			targetnucleonname = p->GetName();
		} else {
			targetnucleonname = targetnucleonpdg;
		}
		TLorentzVector* targetnucleonmomentum = gevtRec->HitNucleon()->P4();
		TVector3 targetnucleonthreemomentum = targetnucleonmomentum->Vect();
		Double_t targetnucleonenergy = targetnucleonmomentum->Energy(); //GeV
		
		// target nucleus:
		Int_t targetnucleuspdg = genieint->InitState().Tgt().Pdg();
		TParticlePDG * targetnucleus = genie::PDGLibrary::Instance()->Find( targetnucleuspdg );
		TString targetnucleusname = "unknown";
		if(targetnucleus){ targetnucleusname = nucleartarget->GetName(); }
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
			fsleptonname = gevtRec->Particle(ipos)->Name();
			fsleptonenergy = gevtRec->Particle(ipos)->Energy();
		}
		
		// other remnants:
		Int_t numfsprotons = genieint->ExclTag().NProtons();
		Int_t numfsneutrons = genieint->ExclTag().NNeutrons();
		Int_t numfspi0 = genieint->ExclTag().NPi0();
		Int_t numfspiplus = genieint->ExclTag().NPiPlus();
		Int_t numfspiminus = genieint->ExclTag().NPiMinus();
		
		if(printneutrinoevent){
			cout<<"This was a "<< procinfoasstring <<" interaction of a "<<probeenergy<<"GeV " 
				<< probepartname << " on a "; 
			
			if( targetnucleonpdg==2212 || targetnucleonpdg==2122 ){ cout<<targetnucleonname<<" in a "; }
			else { cout<<"PDG-Code " << targetnucleonpdg<<" in a "; }
			
			if( targetnucleusname!="unknown"){ cout<<targetnucleusname<<" nucleus, "; }
			else { cout<<"Z=["<<targetnucleusZ<<","<<targetnucleusA<<"] nucleus, "; }
			
			if(remnucpos>-1){ cout<<"producing a "<<remnantnucleusenergy<<"GeV "<<remnantnucleusname; }
			else { cout<<"with no remnant nucleus"; }
			
			if(fsleppos>-1){ cout<<" and a "<<fsleptonenergy<<"GeV "<<fsleptonname<<endl; }
			else{ cout<<" and no final state leptons"<<endl; }
			
			cout<<endl<<"Q^2 was "<<neutrinoq2<<"XX, with final state lepton ejected at Cos(θ)="<<costhfsl<<endl;
			cout<<"Additional final state particles included "<<endl;
			cout<< " N(p) = "       << numfsprotons
				<< " N(n) = "       << numfsneutrons
				<< endl
				<< " N(pi^0) = "    << numfspi0
				<< " N(pi^+) = "    << numfspiplus
				<< " N(pi^-) = "    << numfspiminus
				<<endl;
		}
		
		incidentneutrinoenergies->Fill(probeenergy);
		incidentneutrinoangles->Fill(probeangle);
		finalstateleptonanglesall->Fill(TMath::ACos(costhfsl));
		finalstateleptonenergiesall->Fill(fsleptonenergy);
		
		/* 4. Check if interaction is QE - if not, continue */
		if (!(genieint->ProcInfo().IsQuasiElastic() && genieint->ProcInfo().IsWeakCC())){ continue; }
		
		/*5. If so, load wcsim detector response. */
		// read only first subtrigger; delayed decay detector response is not interesting for primary FSL tracks
		atrigt = b->GetTrigger(0);
		atrigm = m->GetTrigger(0);
		atrigv = v->GetTrigger(0);
		
		Int_t numtracks = atrigt->GetNtrack();
		
		std::vector<Float_t> primaryenergiesvector;
		std::vectr<Double_t> scatteringanglesvector;
		
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
			Int_t primaryparentpdg = GetParenttype();
			if(primaryparentpdg!=0) continue; 				// not a primary
			
			// is it a (anti)muon?
			Int_t primarypdg = nextrack->GetIpnu();
			if(TMath::Abs(primarypdg)!=13) continue; 		// not a muon
			
			// does it start in the tank?
			Int_t primarystartvol = nextrack->GetStartvol();
			if(primarystartvol!=10) continue;				// start volume is not the tank
			
			// does it stop in the mrd?
			Int_t primarystopvol = nextrack->GetStopvol();
			
			// or does it pass completely through the MRD? 
			TVector3 primarystopvertex;
			primarystopvertex.X() = nextrack->GetStop(0);
			primarystopvertex.Y() = nextrack->GetStop(1);
			primarystopvertex.Z() = nextrack->GetStop(2);
			
			// continue if neither
			if(!(primarystopvol==30 /*|| primarystopvertex.Z()>MRD_start+MRD_depth)*/) continue;
			// TODO: disabling check for complete penetration. 
			// We _could_ have a track pass thru multiple MRD layers, then exit MRD and stop with z<MRD_end
			
			/* Requiring z_stop > MRD_end is alone not sufficient - need to check it intercepts the MRD.
			The mrd is as wide as the tank. We can ensure a track enters the mrd by projecting the track
			forward from the start vertex, at the angle between start and stop vertices, and requiring:
			1) at z=MRD_start, x is between +/-(MRD_width/2);
			2) z endpoint is > MRD_start
			To require at least 2-layers of penetration, as above but with z=MRD_layer2 */
			
			TVector3 primarystartvertex;
			primarystartvertex.X() = nextrack->GetStart(0);
			primarystartvertex.Y() = nextrack->GetStart(1);
			primarystartvertex.Z() = nextrack->GetStart(2);
			
			Float_t opp = primarystopvertex.X() - primarystartvertex.X();
			Float_t adj = primarystopvertex.Z() - primarystartvertex.Z();
			Float_t avgtrackanglex = TMath::ATan(opp/adj);
			Float_t xatmrd = primarystartvertex.X() + (MRD_layer2-primarystartvertex.Z())*TMath::Tan(avgtrackanglex);
			
			opp = primarystopvertex.Y() - primarystartvertex.Y();
			Float_t avgtrackangley = TMath::ATan(opp/adj);
			Float_t yatmrd = primarystartvertex.Y() + (MRD_layer2-primarystartvertex.Z())*TMath::Tan(avgtrackangley);
			
			if((TMath::Abs(xatmrd)>MRD_width)||(TMath::Abs(yatmrd)>MRD_height)) 
				continue;	// track does not meet MRD penetration requirement
			
			// if we got to here, we have a muon that stops in, or penetrates at least 2 MRD layers!
			// retrieve any remaining information and record the event
			Float_t primaryenergy = nextrack->GetE();			// starting energy
			Float_t primarymomentummag = nextrack->GetP();		// starting momentum
			TVector3 primarymomentumdir;
			primarymomentumdir.X() = nextrack->GetPdir(0);
			primarymomentumdir.Y() = nextrack->GetPdir(1);
			primarymomentumdir.Z() = nextrack->GetPdir(2);
			Float_t starttrackanglex = TMath::ATan(primarymomentumdir.X()/primarymomentumdir.Z());
			Float_t starttrackangley = TMath::ATan(primarymomentumdir.Y()/primarymomentumdir.Y());
			Float_t starttrackangle = TMath::Max(starttrackanglex,starttrackangley);
			Double_t scatteringangle = k1.Angle(primarymomentumdir);
			// TMath::Cos(scatteringangle) should be the same as costhfsl if neutrino dir == z axis
			
			// Add this primary information to a vector. We _should_ only find one suitable primary
			// and can break once it is found, but... let's just see. Keep wcsim trackID in case.
			Int_t primarytrackid = nextrack->GetId();
			
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
			cout<<"genie lepton angle = "<<TMath::ACos(costhfsl)<<endl;
			cout<<"genie lepton energy = "<<fsleptonenergy<<endl;
			
			std::vector<Int_t>::iterator minit = std::min_element(acceptedtrackids.begin(), acceptedtrackids.end());		
			Int_t indextouse = std::distance(acceptedtrackids.begin(),minit);
			finalstateleptonanglesacceptedwcsim->Fill(TMath::ACos(scatteringanglesvector.at(indextouse));
			finalstateleptonenergiesacceptedwcsim->Fill(primaryenergiesvector.at(indextouse));
		} else if (scatteringanglesvector.size()==1) { // just one matched wcsim track. 
			finalstateleptonanglesacceptedwcsim->Fill(TMath::ACos(scatteringangle));
			finalstateleptonenergiesacceptedwcsim->Fill(primaryenergy);
		} else {
			//cout<<"no accepted wcsim track"<<endl;
		}
		
		finalstateleptonanglesaccepted->Fill(TMath::ACos(costhfsl));
		finalstateleptonenergiesaccepted->Fill(fsleptonenergy);
		
	}
	
	// create scaled histograms with bin contents scaled to the number of POTs in input files:
	TH1F* incidentneutrinoenergiesscaled = new TH1F(incidentneutrinoenergies);
	TH1F* finalstateleptonanglesallscaled = new TH1F(finalstateleptonanglesall);
	TH1F* finalstateleptonenergiesallscaled = new TH1F(finalstateleptonenergiesall);
	TH1F* finalstateleptonanglesacceptedscaled = new TH1F(finalstateleptonanglesaccepted);
	TH1F* finalstateleptonenergiesacceptedscaled = new TH1F(finalstateleptonenergiesaccepted);
	TH1F* finalstateleptonanglesacceptedwcsimscaled = new TH1F(finalstateleptonanglesacceptedwcsim);
	TH1F* finalstateleptonenergiesacceptedwcsimscaled = new TH1F(finalstateleptonenergiesacceptedwcsim);
	
	TH1F* thehistopointersarr[]{
	incidentneutrinoenergiesscaled,
	finalstateleptonanglesallscaled,
	finalstateleptonanglesacceptedscaled,
	finalstateleptonanglesacceptedwcsimscaled,
	finalstateleptonenergiesallscaled,
	finalstateleptonenergiesacceptedscaled,
	finalstateleptonenergiesacceptedwcsimscaled};
	//TODO: modify size if you add more histograms!!!
	std::vector<TH1F*> scaledhistopointers(thehistopointersarr,thehistopointersarr+7);
	
	for(int i=0; i<thehistopointers.size(); i++){
		TH1F* temp = thehistopointers.at(i);
		for(int j=0; j<temp->GetNbinsX(); j++){
			Float_t thecontent = temp->GetBinContent(j);
			thecontent /= totalpots;
			temp->SetBinContent(thecontent);
		}
	}
	
	TH1F* thehistopointersarr2[]{
	incidentneutrinoenergiesscaled,
	finalstateleptonanglesall,
	finalstateleptonanglesaccepted,
	finalstateleptonanglesacceptedwcsim,
	finalstateleptonenergiesall,
	finalstateleptonenergiesaccepted,
	finalstateleptonenergiesacceptedwcsim};
	//TODO: modify size if you add more histograms!!!
	std::vector<TH1F*> histopointers.assign(thehistopointersarr2,thehistopointersarr2+7);
	
	// draw the original histograms
	// ============================
	TLegend* leg;
	TCanvas* canv;
	
	std::vector<TCanvas*> canvaspointers;
	std::vector<TLegend*> legendpointers;
	for(int i=0; i<histopointers.size(); i++){
		TH1F* hist = histopointers.at(i);
		if(i==0){
			canv = new TCanvas();
			scaledcanvaspointers.push_back(canv);
			canv->cd();
			hist->SetLineColor(kRed);
			hist->Draw();
		else if((i-1)%3==0){
			canv = new TCanvas();
			scaledcanvaspointers.push_back(canv);
			canv->cd();
			leg = new TLegend(0.1,0.7,0.48,0.9);
			legendpointers.push_back(leg);
			leg->AddEntry(hist,"incident","l");
			hist->SetLineColor(kRed);
			hist->Draw();
		} else if((i-1)%3==1){
			hist->SetLineColor(kBlue);
			leg->AddEntry(hist,"accepted (truth)","l");
			hist->Draw("same");
		} else {
			hist->SetLineColor(kIndigo);
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
		TH1F* hist = scaledhistopointers.at(i);
		if(i==0){
			canv = new TCanvas();
			scaledcanvaspointers.push_back(canv);
			canv->cd();
			hist->SetLineColor(kRed);
			hist->Draw();
		else if((i-1)%3==0){
			canv = new TCanvas();
			scaledcanvaspointers.push_back(canv);
			canv->cd();
			leg = new TLegend(0.1,0.7,0.48,0.9);
			legendpointers.push_back(leg);
			leg->AddEntry(hist,"incident","l");
			hist->SetLineColor(kRed);
			hist->Draw();
		} else if((i-1)%3==1){
			hist->SetLineColor(kBlue);
			leg->AddEntry(hist,"accepted (truth)","l");
			hist->Draw("same");
		} else {
			hist->SetLineColor(kIndigo);
			leg->AddEntry(hist,"accepted (reco)","l");
			hist->Draw("same");
			leg->Draw();
		}
	}
	
	// cleanup
	// =======
	if(c) c->ResetBranchAddresses();
	if(tankflux) tankflux->ResetBranchAddresses();
	if(tankmeta) tankmeta->ResetBranchAddresses();
	if(dirtfile) dirtfile->Close(); // do we need to do this with a TChain? 
	if(c) delete c;					// ??
	
	if(gtree) gtree->ResetBranchAddresses();
	if(geniefile) geniefile->Close();
	// should clean up gtree
	if(genierecordval) delete genierecordval;

	if(wcsimT) wcsimT->ResetBranchAddresses();
	if(wcsimfile) wcsimfile->Close();
	// should clean up wcsimT
	if(nuprimarybranchval){delete[] nuprimarybranchval;}	// ? branch array
	
	// delete histograms
	delete incidentneutrinoenergies;
	delete finalstateleptonanglesall;
	delete finalstateleptonenergiesall;
	delete finalstateleptonanglesaccepted;
	delete finalstateleptonenergiesaccepted;
	delete finalstateleptonanglesacceptedwcsim;
	delete finalstateleptonenergiesacceptedwcsim;
	
	// delete scaled histograms
	delete incidentneutrinoenergiesscaled;
	delete finalstateleptonanglesallscaled;
	delete finalstateleptonenergiesallscaled;
	delete finalstateleptonanglesacceptedscaled;
	delete finalstateleptonenergiesacceptedscaled;
	delete finalstateleptonanglesacceptedwcsimscaled;
	delete finalstateleptonenergiesacceptedwcsimscaled;
	
	// delete canvases and legends
	std::vector<TCanvas*> canvaspointers;
	for(int i=0; i<canvaspointers.size(); i++){
		TCanvas* canv = canvaspointers.at(i);
		canv->SaveAs(TString::Format("histograms_%d.png",i));
		if(canv) delete canv;
		TLegend* leg = legendpointers.at(i);
		if(leg) delete leg;
	}
	std::vector<TCanvas*> canvaspointers;
	
	// delete scaled canvases and legends
	for(int i=0; i<scaledcanvaspointers.size(); i++){
		TCanvas* canv = scaledcanvaspointers.at(i);
		canv->SaveAs(TString::Format("scaled_histograms_%d.png",i));
		if(canv) delete canv;
		TLegend* leg = scaledlegendpointers.at(i);
		if(leg) delete leg;
	}

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
