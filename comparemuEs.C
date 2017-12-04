/* vim:set noexpandtab tabstop=4 wrap */
#include "TROOT.h"
#include "TSystem.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TApplication.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THStack.h"
#include "TNtuple.h"
#include "TParticlePDG.h"
#include "TDatabasePDG.h"
#include "TVector3.h"
#include "TMath.h"
#include "TCanvas.h"
//#include "/annie/app/users/moflaher/annie_rootlogon.C"
//void annie_rootlogon();
//void Preliminary(TString data="");
//void Preliminary(TString data="");
//void Simulation();
//void SimulationSide();
//void CenterTitles(TH1* histo);
//void CenterTitles(THStack* histo);
#include <string>
#include <iostream>

#ifndef VERBOSE
#define VERBOSE 1
#endif

const Float_t tank_start = 15.70;           // front face of the tank in cm
const Float_t tank_radius = 152.4;          // tank radius in cm
const Float_t tank_halfheight = 198.;       // tank half height in cm
const Float_t tank_yoffset = -14.46;        // tank y offset in cm

double CalculateEventQ2(double recoMuonEnergy, double recoNeutrinoEnergy, double recoMuonAngle);

void comparemuEs(){
	const char* pwd = gSystem->pwd();
	//XXX SET BEFORE RUNNING  XXX
	//XXX XXX XXX XXX XXX XXX XXX
	std::string outstring = std::string(pwd)+"/out/q2comparison";
	outstring+="/30ps_proper_split/pmts_only";
	const char* outfilename="q2comparison.root";
	const char* outdir=outstring.c_str();
#ifdef VERBOSE
	cout<<"setting output directory to "<<outdir<<endl;
#endif
	
#ifdef VERBOSE
	cout<<"loading muon reco file"<<endl;
#endif
	// muon file:
	//TFile *f1 = TFile::Open("/annie/data/users/edrakopo/OUTrecomuonE.root");
	//TFile *f1 = TFile::Open("/annie/data/users/edrakopo/OUTrecomuonENEW.root");
	//TFile *f1 = TFile::Open("/annie/data/users/edrakopo/OUTrecomuonENEWB.root");
	//TFile *f1 = TFile::Open("/annie/data/users/edrakopo/OUTrecomuonENEWC2.root");
	//TFile *f1 = TFile::Open("/annie/data/users/edrakopo/OUTrecomuonENEW2rndsplit.root");
	//TFile *f1 = TFile::Open("/annie/data/users/edrakopo/OUTrecomuonENEW2rndsplit_PMTsonly.root");
	TFile *f1 = TFile::Open("/pnfs/annie/scratch/users/edrakopo/NEWOUTRecoMuonMRD_PMTsonly.root");
	//TFile *f1 = TFile::Open("/pnfs/annie/scratch/users/edrakopo/NEWOUTRecoMuonMRD.root");
	// neutrino file:
	//TFile *f2 = TFile::Open("/annie/data/users/edrakopo/OUTreconeutrinoE.root");
	//TFile *f2 = TFile::Open("/annie/data/users/edrakopo/OUTrecomuonEneuNEWB.root");
	//TFile *f2 = TFile::Open("/annie/data/users/edrakopo/OUTreconeutrinoENEWC2.root");
	//TFile *f2 = TFile::Open("/annie/data/users/edrakopo/OUTreconeutrinoENEW2rndsplit.root");
	//TFile *f2 = TFile::Open("/annie/data/users/edrakopo/OUTreconeutrinoENEW2rndsplit_PMTsonly.root");
	TFile *f2 = TFile::Open("/pnfs/annie/scratch/users/edrakopo/NEWOUTRecoNeutrinoMRD_PMTsonly.root");
	//TFile *f2 = TFile::Open("/pnfs/annie/scratch/users/edrakopo/NEWOUTRecoNeutrinoMRD.root");
	
	// retrieve trees from files
	TTree* t1 = (TTree*)f1->Get("tuple");
//	t1->SetBranchStatus("*",0);
//	t1->SetBranchStatus("trueKE",1);
//	t1->SetBranchStatus("recoKE",1);
//	t1->SetBranchStatus("dirX",1);
//	t1->SetBranchStatus("dirY",1);
//	t1->SetBranchStatus("dirZ",1);
//	t1->SetBranchStatus("TrueMomentumTransfer",1);
//	t1->SetBranchStatus("TrueMuonAngle",1);
	TVector3 recomuondirection;
	float truemuE, recomuE, recomudirx, recomudiry, recomudirz, trueq2, trueangle;
	t1->SetBranchAddress("trueKE",&truemuE);
	t1->SetBranchAddress("recoKE",&recomuE);
	t1->SetBranchAddress("dirX",&recomudirx);
	t1->SetBranchAddress("dirY",&recomudiry);
	t1->SetBranchAddress("dirZ",&recomudirz);
	t1->SetBranchAddress("TrueMomentumTransfer",&trueq2);
	t1->SetBranchAddress("TrueMuonAngle",&trueangle);
	
#ifdef VERBOSE
	cout<<"loading neutrino reco file"<<endl;
#endif
	TTree* t2 = (TTree*)f2->Get("tuple");
//	t2->SetBranchStatus("*",0);
//	t2->SetBranchStatus("trueKE",1);
//	t2->SetBranchStatus("recoKE",1);
//	t2->SetBranchStatus("vtxX",1);
//	t2->SetBranchStatus("vtxY",1);
//	t2->SetBranchStatus("vtxZ",1);
	float truenuE, reconuE, recovtxx, recovtxy, recovtxz;
	t2->SetBranchAddress("trueKE",&truenuE);
	t2->SetBranchAddress("recoKE",&reconuE);
	t2->SetBranchAddress("vtxZ",&recovtxx);
	t2->SetBranchAddress("vtxY",&recovtxy);
	t2->SetBranchAddress("vtxX",&recovtxz);
	
#ifdef VERBOSE
	cout<<"loading MC truth file"<<endl;
#endif
	// pull the true Q2 information for comparison. 
	// We can only compare distributions for now, not event-by-event, as event # or true Q2 were not saved
	TFile* f3 = TFile::Open("/pnfs/annie/scratch/users/moflaher/trueQEvertexinfo_extv2.root");
	TTree* t3 = (TTree*)f3->Get("vertextreenocuts");
	t3->SetBranchStatus("*",0);
	t3->SetBranchStatus("NeutrinoEnergy",1);
	t3->SetBranchStatus("MuonEnergy",1);
	t3->SetBranchStatus("MuonAngle",1);
	t3->SetBranchStatus("MuonDirection",1);
	t3->SetBranchStatus("MomentumTransfer",1);
	t3->SetBranchStatus("TrackLengthInMrd",1);
	t3->SetBranchStatus("MuonStartVertex",1);
//	Event by event is not needed until we have event number for reconstructed events
//	double truenuEall, truemuEall, trueangleall, trueQ2all, trackleninmrd;
//	TVector3 truevertexall;
//	t3->SetBranchAddress("NeutrinoEnergy",&truenuEall);
//	t3->SetBranchAddress("MuonEnergy",&truemuEall);
//	t3->SetBranchAddress("MuonAngle",&trueangleall);
//	t3->SetBranchAddress("MomentumTransfer",&trueQ2all);
//	t3->SetBranchAddress("TrackLengthInMrd",&trackleninmrd);
//	t3->SetBranchAddress("MuonStartVertex",&truevertexall);
	
	// output file
	TFile* fileout = new TFile(TString::Format("%s/%s",outdir,outfilename),"RECREATE");
#ifdef VERBOSE
	cout<<"creating output file "<<fileout->GetName()<<endl;
#endif
	TTree* treeout = new TTree("q2tree","Reco vs True Comparison of Q2");
	double outtruenuE, outtruemuE, outreconuE, outrecomuE, outtrueangle, outrecoangle, outtrueq2, outrecoq2, 
	outvtxerr, outangerr, outmuEerr, outnuEerr, outq2err, outrecodirx, outrecodiry, outrecodirz;
	bool outinfidvol, outhasmrdtrack;
	TVector3 outtruevtx, outrecovtx, outtruedir, outrecodir;
	treeout->Branch("TrueNuE",&outtruenuE);
	treeout->Branch("RecoNuE",&outreconuE);
	treeout->Branch("TrueMuE",&outtruemuE);
	treeout->Branch("RecoMuE",&outrecomuE);
	treeout->Branch("TrueAngle",&outtrueangle);
	treeout->Branch("RecoAngle",&outrecoangle);
	treeout->Branch("TrueQ2",&outtrueq2);
	treeout->Branch("RecoQ2",&outrecoq2);
	treeout->Branch("VtxErr",&outvtxerr);
	treeout->Branch("AngErr",&outangerr);
	treeout->Branch("MuEerr",&outmuEerr);
	treeout->Branch("NuEerr",&outnuEerr);
	treeout->Branch("Q2Err",&outq2err);
	treeout->Branch("InFidVol",&outinfidvol);
	treeout->Branch("HasMrdTrack",&outhasmrdtrack);
	treeout->Branch("TrueVtx",&outtruevtx);
	treeout->Branch("RecoVtx",&outrecovtx);
	treeout->Branch("TrueDir",&outtruedir);
	treeout->Branch("RecoDir",&outrecodir);
	treeout->Branch("RecoDirX",&outrecodirx);
	treeout->Branch("RecoDirY",&outrecodiry);
	treeout->Branch("RecoDirZ",&outrecodirz);
	
#ifdef VERBOSE
	cout<<"creating truth histograms"<<endl;
#endif
	// histograms to compare distributions
	/* true sample histograms */
	TH1D truenuKEhist("truenuKEhist","True Neutrino KE",50,0,2500);
	t2->Draw("trueKE>>truenuKEhist");
	Preliminary("Simulation"); CenterTitles(&truenuKEhist);
	TH1D truemuKEhist("truemuKEhist","True Muon KE",50,0,2500);
	t1->Draw("trueKE>>truemuKEhist");
	Preliminary("Simulation"); CenterTitles(&truemuKEhist);
	TH1D trueanglehist("trueanglehist","True Lepton Angle",50,0,TMath::Pi());
	t1->Draw("TrueMuonAngle>>trueanglehist");
	Preliminary("Simulation"); CenterTitles(&trueanglehist);
	TH1D trueQ2hist("trueQ2hist","True Momentum Transfer",50,0,2);
	t1->Draw("(TrueMomentumTransfer)>>trueQ2hist");
	Preliminary("Simulation"); CenterTitles(&trueQ2hist);
	TH3D truevertexhist("truevertexhist","True Vertex",50,-150,150,50,-220,220,50,-150,150);
	// do not yet have XXX remember to use 'goff' option
	Preliminary("Simulation"); CenterTitles(&truevertexhist);
	TH3D truedirhist("truedirhist","True Direction",50,-1.0,1.0,50,-1.0,1.0,50,-1.0,1.0);
	// do not yet have XXX remember to use 'goff' option
	Preliminary("Simulation"); CenterTitles(&truedirhist);
	
#ifdef VERBOSE
	cout<<"creating reco histograms"<<endl;
#endif
	/* reco sample histograms */
	TH1D reconuKEhist("reconuKEhist","Reco Neutrino KE",50,0,2500);
	t2->Draw("recoKE>>reconuKEhist");
	Preliminary("Simulation"); CenterTitles(&reconuKEhist);
	TH1D recomuKEhist("recomuKEhist","Reco Muon KE",50,0,2500);
	t1->Draw("recoKE>>recomuKEhist");
	Preliminary("Simulation"); CenterTitles(&recomuKEhist);
	TH1D recoanglehist("recoanglehist","Reco Lepton Angle",50,0,TMath::Pi());
	// fill in loop
	Preliminary("Simulation"); CenterTitles(&recoanglehist);
	TH1D recoQ2hist("recoQ2hist","Reco Momentum Transfer",50,0,2);
	// fill in loop
	Preliminary("Simulation"); CenterTitles(&recoQ2hist);
	TH3D recovertexhist("recovertexhist","Reco Vertex",50,-150,150,50,-220,220,50,-150,150);
	t2->Draw("vtxX:vtxY:vtxZ>>recovertexhist","","goff");
	Preliminary("Simulation"); CenterTitles(&recovertexhist);
	TH3D recodirhist("recodirhist","Reco Direction",50,-1.0,1.0,50,-1.0,1.0,50,-1.0,1.0);
	t1->Draw("dirX:dirY:dirZ>>recodirhist","","goff");
	Preliminary("Simulation"); CenterTitles(&recodirhist);
	
#ifdef VERBOSE
	cout<<"creating input histograms"<<endl;
#endif
	/* true all histograms */
	TH1D truenuKEallhist("truenuKEallhist","True Neutrino KE (All)",50,0,2500);
	t3->Draw("NeutrinoEnergy>>(truenuKEallhist*1000.)");
	Preliminary("Simulation"); 
	TH1D truemuKEallhist("truemuKEallhist","True Muon KE (All)",50,0,2500);
	t3->Draw("(MuonEnergy*1000.)>>truemuKEallhist");
	Preliminary("Simulation"); CenterTitles(&truemuKEallhist);
	TH1D trueangleallhist("trueangleallhist","True Lepton Angle (All)",50,0,TMath::Pi());
	t3->Draw("MuonAngle>>trueangleallhist");
	Preliminary("Simulation"); CenterTitles(&trueangleallhist);
	TH1D trueQ2allhist("trueQ2allhist","True Momentum Transfer (All)",50,0,2);
	t3->Draw("(MomentumTransfer)>>trueQ2allhist");
	Preliminary("Simulation"); CenterTitles(&trueQ2allhist);
	TH3D truevertexallhist("truevertexallhist","True Vertex (All)",50,-150,150,50,-220,220,50,-150,150);
	t3->Draw("MuonStartVertex.X():MuonStartVertex.Y():MuonStartVertex.Z()>>truevertexallhist","","goff");
	Preliminary("Simulation"); CenterTitles(&truevertexallhist);
	TH3D truedirallhist("truedirallhist","True Direction (All)",50,-1.0,1.0,50,-1.0,1.0,50,-1.0,1.0);
	t3->Draw("MuonDirection.X():MuonDirection.Y():MuonDirection.Z()>>truedirallhist","","goff");
	Preliminary("Simulation"); CenterTitles(&truedirallhist);
	
#ifdef VERBOSE
	cout<<"creating fractional error histograms"<<endl;
#endif
	/* fractional errors */
	TH1D nuEfracerrhist("nuEfracerrhist","Relative Neutrino Energy Error",50,-1.0,1.0);
	Preliminary("Simulation"); CenterTitles(&nuEfracerrhist);
	TH1D muEfracerrhist("muEfracerrhist","Relative Muon Energy Error",50,-1.0,1.0);
	Preliminary("Simulation"); CenterTitles(&muEfracerrhist);
	TH1D angfracerrhist("angfracerrhist","Relative Lepton Angle Error",50,-1.0,1.0);
	Preliminary("Simulation"); CenterTitles(&angfracerrhist);
	TH1D q2fracerrhist("q2fracerrhist","Relative Q^{2} Error",50,-2.0,2.0);
	Preliminary("Simulation"); CenterTitles(&q2fracerrhist);
	
#ifdef VERBOSE
	cout<<"creating 2D histograms"<<endl;
#endif
	/* true val vs reco val */
	TH2D nuEtruevsreco("nuEtruevsreco","True vs Reco Neutrino Energy",50,0,2500,50,0,2500);
	Preliminary("Simulation"); CenterTitles(&nuEtruevsreco);
	TH2D muEtruevsreco("muEtruevsreco","True vs Reco Muon Energy",50,0,2500,50,0,2500);
	Preliminary("Simulation"); CenterTitles(&muEtruevsreco);
	TH2D angtruevsreco("angtruevsreco","True vs Reco Lepton Angle",50,0,TMath::Pi(),50,0,TMath::Pi());
	Preliminary("Simulation"); CenterTitles(&angtruevsreco);
	TH2D q2truevsreco("q2truevsreco","True vs Reco Q^{2}",50,0,1.0,50,0,1.0);
	Preliminary("Simulation"); CenterTitles(&q2truevsreco);
	
	/* fractional errors vs value */
	TH2D nuEfracerrhist2("nuEfracerrhist2","Relative Neutrino Energy Error vs True Neutrino Energy",50,0,2500,50,-1.0,1.0);
	Preliminary("Simulation"); CenterTitles(&nuEfracerrhist2);
	TH2D muEfracerrhist2("muEfracerrhist2","Relative Muon Energy Error vs True Muon Energy",50,0,2500,50,-1.0,1.0);
	Preliminary("Simulation"); CenterTitles(&muEfracerrhist2);
	TH2D angfracerrhist2("angfracerrhist2","Relative Lepton Angle Error vs True Angle",50,0,TMath::Pi(),50,-1.0,1.0);
	Preliminary("Simulation"); CenterTitles(&angfracerrhist2);
	TH2D q2fracerrhist2("q2fracerrhist2","Relative Q^{2} Error vs True Q^{2}",50,0,2,50,-2.0,2.0);
	Preliminary("Simulation"); CenterTitles(&q2fracerrhist2);
	
#ifdef VERBOSE
	cout<<"looping over reconstructed events"<<endl;
#endif
	// loop over reco trees and combine them:
	for(int i=0; i<t2->GetEntries();i++){
		t1->GetEntry(i);
		t2->GetEntry(i);
		outtruenuE=truenuE;
		outreconuE=reconuE;
		outtruemuE=truemuE;
		outrecomuE=recomuE;
		outtruevtx=TVector3(-1,-1,-1); // do not yet have
		outrecovtx=TVector3(recovtxx, recovtxy, recovtxz);
		outtruedir=TVector3(-1,-1,-1); // do not yet have
		outrecodir=TVector3(recomudirx, recomudiry, recomudirz);
		outrecodirx=recomudirx; outrecodiry=recomudiry; outrecodirz=recomudirz;
		outtrueangle=trueangle;
		outrecoangle=outrecodir.Angle(TVector3(0,0,1));
		recoanglehist.Fill(outrecoangle);
		outtrueq2=trueq2;
		outrecoq2=CalculateEventQ2(recomuE, reconuE, outrecoangle)/1000000.;
		recoQ2hist.Fill(outrecoq2);
		outnuEerr=outtruenuE-outreconuE;
		outmuEerr=outtruemuE-outrecomuE;
		outangerr=outtruedir.Angle(outrecodir);
		outvtxerr=(outtruevtx-outrecovtx).Mag();
		outq2err=outtrueq2-outrecoq2;
		outinfidvol=((sqrt(pow(outtruevtx.X(),2.)+pow(outtruevtx.Z()/*-tank_start-tank_radius*/,2.))<tank_radius)&&
		             (abs(outtruevtx.Y()/*-tank_yoffset*/)<tank_halfheight)&&((outtruevtx.Z()/*-tank_start-tank_radius*/)<0));
		outhasmrdtrack=false;          // do not yet have (trackleninmrd>0);
		
		// true vs reco histos
		nuEtruevsreco.Fill(truenuE,reconuE);
		muEtruevsreco.Fill(truemuE,recomuE);
		angtruevsreco.Fill(trueangle,outrecoangle);
		q2truevsreco.Fill(outtrueq2,outrecoq2);
		
		// fractional error histos
		nuEfracerrhist.Fill((truenuE-reconuE)/truenuE);
		muEfracerrhist.Fill((truemuE-recomuE)/truemuE);
		angfracerrhist.Fill((outtrueangle-outrecoangle)/outtrueangle);
		q2fracerrhist.Fill((outtrueq2-outrecoq2)/outtrueq2);
		// fractional error histos vs true val
		nuEfracerrhist2.Fill(truenuE,(truenuE-reconuE)/truenuE);
		muEfracerrhist2.Fill(truemuE,(truemuE-recomuE)/truemuE);
		angfracerrhist2.Fill(outtrueangle,(outtrueangle-outrecoangle)/outtrueangle);
		q2fracerrhist2.Fill(outtrueq2,(outtrueq2-outrecoq2)/outtrueq2);
		
		treeout->Fill();
	}
	
#ifdef VERBOSE
	cout<<"saving truth histograms"<<endl;
#endif
	// Save true histos
	truenuKEhist.SetLineColor(kRed);
	truenuKEhist.GetXaxis()->SetTitle("Neutrino Energy [MeV]");
	truenuKEhist.GetYaxis()->SetTitle("Fraction of Events");
	CenterTitles(&truenuKEhist);
	truenuKEhist.Write();
	
	truemuKEhist.SetLineColor(kRed);
	truemuKEhist.GetXaxis()->SetTitle("Muon Energy [MeV]");
	truemuKEhist.GetYaxis()->SetTitle("Fraction of Events");
	CenterTitles(&truemuKEhist);
	truemuKEhist.Write();
	
	trueanglehist.SetLineColor(kRed);
	trueanglehist.GetXaxis()->SetTitle("Lepton Angle [rads]");
	trueanglehist.GetYaxis()->SetTitle("Fraction of Events");
	CenterTitles(&trueanglehist);
	trueanglehist.Write();
	
	trueQ2hist.SetLineColor(kRed);
	trueQ2hist.GetXaxis()->SetTitle("Q^{2} [(MeV/c)^{2}]");
	trueQ2hist.GetYaxis()->SetTitle("Fraction of Events");
	CenterTitles(&trueQ2hist);
	trueQ2hist.Write();
	
	truevertexhist.Write();
	truedirhist.Write();
	
#ifdef VERBOSE
	cout<<"saving reco histograms"<<endl;
#endif
	// Save reco histos
	reconuKEhist.SetLineColor(kBlue);
	reconuKEhist.GetXaxis()->SetTitle("Neutrino Energy [MeV]");
	reconuKEhist.GetYaxis()->SetTitle("Fraction of Events");
	CenterTitles(&reconuKEhist);
	reconuKEhist.Write();
	
	recomuKEhist.SetLineColor(kBlue);
	recomuKEhist.GetXaxis()->SetTitle("Muon Energy [MeV]");
	recomuKEhist.GetYaxis()->SetTitle("Fraction of Events");
	CenterTitles(&recomuKEhist);
	recomuKEhist.Write();
	
	recoanglehist.SetLineColor(kBlue);
	recoanglehist.GetXaxis()->SetTitle("Lepton Angle [rads]");
	recoanglehist.GetYaxis()->SetTitle("Fraction of Events");
	CenterTitles(&recoanglehist);
	recoanglehist.Write();
	
	recoQ2hist.SetLineColor(kBlue);
	recoQ2hist.GetXaxis()->SetTitle("Q^{2} [(MeV/c)^{2}]");
	recoQ2hist.GetYaxis()->SetTitle("Fraction of Events");
	CenterTitles(&recoQ2hist);
	recoQ2hist.Write();
	
	recovertexhist.Write();
	recodirhist.Write();
	
#ifdef VERBOSE
	cout<<"saving input histograms"<<endl;
#endif
	// Save true histos without selection cuts
	truenuKEallhist.SetLineColor(kRed);
	truenuKEallhist.GetXaxis()->SetTitle("Neutrino Energy [MeV]");
	truenuKEallhist.GetYaxis()->SetTitle("Fraction of Events");
	CenterTitles(&truenuKEallhist);
	truenuKEallhist.Write();
	
	truemuKEallhist.SetLineColor(kRed);
	truemuKEallhist.GetXaxis()->SetTitle("Muon Energy [MeV]");
	truemuKEallhist.GetYaxis()->SetTitle("Fraction of Events");
	CenterTitles(&truemuKEallhist);
	truemuKEallhist.Write();
	
	trueangleallhist.SetLineColor(kRed);
	trueangleallhist.GetXaxis()->SetTitle("Lepton Angle [rads]");
	trueangleallhist.GetYaxis()->SetTitle("Fraction of Events");
	CenterTitles(&trueangleallhist);
	trueangleallhist.Write();
	
	trueQ2allhist.SetLineColor(kRed);
	trueQ2allhist.GetXaxis()->SetTitle("Q^{2} [(MeV/c)^{2}]");
	trueQ2allhist.GetYaxis()->SetTitle("Fraction of Events");
	CenterTitles(&trueQ2allhist);
	trueQ2allhist.Write();
	
	truevertexallhist.Write();
	truedirallhist.Write();
	
#ifdef VERBOSE
	cout<<"saving error histograms"<<endl;
#endif
	// Save fractional error histos
	nuEfracerrhist.GetXaxis()->SetTitle("Fractional Neutrino Energy Error");
	nuEfracerrhist.GetYaxis()->SetTitle("Fraction of Events");
	CenterTitles(&nuEfracerrhist);
	nuEfracerrhist.Write();
	
	muEfracerrhist.GetXaxis()->SetTitle("Fractional Muon Energy Error");
	muEfracerrhist.GetYaxis()->SetTitle("Fraction of Events");
	CenterTitles(&muEfracerrhist);
	muEfracerrhist.Write();
	
	angfracerrhist.GetXaxis()->SetTitle("Fractional Lepton Angle Error");
	angfracerrhist.GetYaxis()->SetTitle("Fraction of Events");
	CenterTitles(&angfracerrhist);
	angfracerrhist.Write();
	
	q2fracerrhist.GetXaxis()->SetTitle("Fractional Q^{2} Error");
	q2fracerrhist.GetYaxis()->SetTitle("Fraction of Events");
	CenterTitles(&q2fracerrhist);
	q2fracerrhist.Write();
	
#ifdef VERBOSE
	cout<<"saving 2D histograms"<<endl;
#endif
	// True vs Reco histos
	nuEtruevsreco.GetXaxis()->SetTitle("True Neutrino Energy [MeV]");
	nuEtruevsreco.GetYaxis()->SetTitle("Reconstructed Neutrino Energy [MeV]");
	CenterTitles(&nuEtruevsreco);
	nuEtruevsreco.Write();

	muEtruevsreco.GetXaxis()->SetTitle("True Muon Energy [MeV]");
	muEtruevsreco.GetYaxis()->SetTitle("Reconstructed Muon Energy [MeV]");
	CenterTitles(&muEtruevsreco);
	muEtruevsreco.Write();

	angtruevsreco.GetXaxis()->SetTitle("True Lepton Angle [rads]");
	angtruevsreco.GetYaxis()->SetTitle("Reconstructed Lepton Angle [rads]");
	CenterTitles(&angtruevsreco);
	angtruevsreco.Write();

	q2truevsreco.GetXaxis()->SetTitle("True Q^{2} [(GeV/c)^{2}]");
	q2truevsreco.GetYaxis()->SetTitle("Reconstructed Q^{2} [(GeV/c)^{2}]");
	CenterTitles(&q2truevsreco);
	q2truevsreco.Write();
	
	// Error histos
	nuEfracerrhist2.GetYaxis()->SetTitle("Fractional Neutrino Energy Error");
	nuEfracerrhist2.GetXaxis()->SetTitle("True Neutrino Energy [MeV]");
	CenterTitles(&nuEfracerrhist2);
	nuEfracerrhist2.Write();
	
	muEfracerrhist2.GetYaxis()->SetTitle("Fractional Muon Energy Error");
	muEfracerrhist2.GetXaxis()->SetTitle("True Muon Energy [MeV]");
	CenterTitles(&muEfracerrhist2);
	muEfracerrhist2.Write();
	
	angfracerrhist2.GetYaxis()->SetTitle("Fractional Scattering Error");
	angfracerrhist2.GetXaxis()->SetTitle("True Lepton Angle");
	CenterTitles(&angfracerrhist2);
	angfracerrhist2.Write();
	
	q2fracerrhist2.GetYaxis()->SetTitle("Fractional Q^{2} Error");
	q2fracerrhist2.GetXaxis()->SetTitle("True Q^{2} [(GeV/c)^{2}]");
	CenterTitles(&q2fracerrhist2);
	q2fracerrhist2.Write();
	
#ifdef VERBOSE
	cout<<"Saving PNGs of stacked plots"<<endl;
#endif
	// Save PNGs of the true and reco histograms overlaid
	// XXX N.B: alternative to Scale(norm/hist.Integral()) is histo.DrawNormalized(option,norm)
	// this draws a normalized copy of the histo. You either need to call pad->Clear(), or delete the returned pointer!
	TCanvas c1;
	THStack* thestack = new THStack("thestack","Stacked Histos");
	TLegend* leg;
	thestack->SetTitle("");
	// thestack->GetHists()->Clear(); <<< clears the histograms but doesn't allow the X-axis range to rescale with new histos.
	truenuKEhist.Scale(1/truenuKEhist.Integral());
	truenuKEhist.SetTitle("True");
	thestack->Add(&truenuKEhist);
	reconuKEhist.Scale(1/reconuKEhist.Integral());
	reconuKEhist.SetTitle("Reco");
	thestack->Add(&reconuKEhist);
	thestack->Draw("nostack");
	thestack->GetXaxis()->SetTitle("Neutrino Energy [MeV]");
	thestack->GetYaxis()->SetTitle("Fraction of Events");
	CenterTitles(thestack);
	//thestack->SetTitle("Neutrino Energy");
	leg = c1.BuildLegend();
	leg->SetFillStyle(0);
	Preliminary("Simulation"); // must get called AFTER BuildLegend, else it gets added to it!
	gPad->Modified(); gPad->Update(); // make sure it's really (re)drawn
	c1.SaveAs(TString::Format("%s/reconstructed_neutrino_energy.png",outdir));
	
#ifdef VERBOSE
	cout<<"2"<<endl;
#endif
	delete thestack; thestack = new THStack("thestack","Stacked Histos");
	truemuKEhist.Scale(1/truemuKEhist.Integral());
	truemuKEhist.SetTitle("True");
	thestack->Add(&truemuKEhist);
	recomuKEhist.Scale(1/recomuKEhist.Integral());
	recomuKEhist.SetTitle("Reco");
	thestack->Add(&recomuKEhist);
	thestack->Draw("nostack");
	gPad->Modified(); gPad->Update(); // make sure it's really (re)drawn
	thestack->GetXaxis()->SetTitle("Muon Energy [MeV]");
	thestack->GetYaxis()->SetTitle("Fraction of Events");
	CenterTitles(thestack);
	thestack->Draw("nostack");
	leg = c1.BuildLegend();
	leg->SetFillStyle(0);
	Preliminary("Simulation");
	gPad->Modified(); gPad->Update(); // make sure it's really (re)drawn
	c1.SaveAs(TString::Format("%s/reconstructed_muon_energy.png",outdir));
	
#ifdef VERBOSE
	cout<<"3"<<endl;
#endif
	/*c1.Clear();*/ thestack->GetHists()->Clear(); delete thestack; thestack = new THStack("thestack","Stacked Histos");
	trueanglehist.Scale(1/trueanglehist.Integral());
	trueanglehist.SetTitle("True");
	thestack->Add(&trueanglehist);
	//obselete: compare w/ distribution of all events
	//trueangleallhist.Scale(1/trueangleallhist.Integral());
	//thestack->Add(&trueangleallhist);
	recoanglehist.Scale(1/recoanglehist.Integral());
	recoanglehist.SetTitle("Reco");
	thestack->Add(&recoanglehist);
	thestack->Draw("nostack");
	gPad->Modified(); gPad->Update(); // make sure it's really (re)drawn
	thestack->GetXaxis()->SetTitle("Lepton Angle [rads]");
	thestack->GetYaxis()->SetTitle("Fraction of Events");
	CenterTitles(thestack);
	thestack->Draw("nostack");
	leg = c1.BuildLegend();
	leg->SetFillStyle(0);
	Preliminary("Simulation");
	gPad->Modified(); gPad->Update(); // make sure it's really (re)drawn
	c1.SaveAs(TString::Format("%s/reconstructed_scattering_angle.png",outdir));
	
#ifdef VERBOSE
	cout<<"4"<<endl;
#endif
	delete thestack; thestack = new THStack("thestack","Stacked Histos");
	trueQ2hist.Scale(1/trueQ2hist.Integral());
	trueQ2hist.SetTitle("True");
	thestack->Add(&trueQ2hist);
	//obselete: compare all events
	//trueQ2allhist.Scale(1/trueQ2allhist.Integral());
	//thestack->Add(&trueQ2allhist);
	recoQ2hist.Scale(1/recoQ2hist.Integral());
	recoQ2hist.SetTitle("Reco");
	thestack->Add(&recoQ2hist);
	thestack->Draw("nostack");
	gPad->Modified(); gPad->Update(); // make sure it's really (re)drawn
	thestack->GetXaxis()->SetTitle("Q^{2} [(GeV/c)^{2}]");
	thestack->GetYaxis()->SetTitle("Fraction of Events");
	CenterTitles(thestack);
	thestack->Draw("nostack");
	leg = c1.BuildLegend();
	leg->SetFillStyle(0);
	Preliminary("Simulation");
	gPad->Modified(); gPad->Update(); // make sure it's really (re)drawn
	c1.SaveAs(TString::Format("%s/reconstructed_momentum_transfer.png",outdir));
	
#ifdef VERBOSE
	cout<<"saving PNGs of error histos"<<endl;
#endif
	// fractional error histograms
	c1.Clear();
	nuEfracerrhist.Scale(1/nuEfracerrhist.Integral());
	nuEfracerrhist.Draw();
	Preliminary("Simulation");
	c1.SaveAs(TString::Format("%s/reconstructed_nue_error.png",outdir));
	
#ifdef VERBOSE
	cout<<"2"<<endl;
#endif
	c1.Clear();
	muEfracerrhist.Scale(1/muEfracerrhist.Integral());
	muEfracerrhist.Draw();
	Preliminary("Simulation");
	c1.SaveAs(TString::Format("%s/reconstructed_mue_error.png",outdir));
	
#ifdef VERBOSE
	cout<<"3"<<endl;
#endif
	c1.Clear();
	angfracerrhist.Scale(1/angfracerrhist.Integral());
	angfracerrhist.Draw();
	Preliminary("Simulation");
	c1.SaveAs(TString::Format("%s/reconstructed_angle_error.png",outdir));
	
#ifdef VERBOSE
	cout<<"4"<<endl;
#endif
	c1.Clear();
	q2fracerrhist.Scale(1/q2fracerrhist.Integral());
	q2fracerrhist.Draw();
	Preliminary("Simulation");
	c1.SaveAs(TString::Format("%s/reconstructed_q2_error.png",outdir));
	
#ifdef VERBOSE
	cout<<"Saving PNGs of True vs Reco histos"<<endl;
#endif
	c1.Clear();
	//nuEtruevsreco.Scale(1/nuEtruevsreco.Integral());
	nuEtruevsreco.Draw("colz");
	Preliminary("Simulation");
	c1.SaveAs(TString::Format("%s/nue_reco_vs_true.png",outdir));
	
#ifdef VERBOSE
	cout<<"2"<<endl;
#endif
	c1.Clear();
	//muEtruevsreco.Scale(1/muEtruevsreco.Integral());
	muEtruevsreco.Draw("colz");
	Preliminary("Simulation");
	c1.SaveAs(TString::Format("%s/mue_reco_vs_true.png",outdir));
	
#ifdef VERBOSE
	cout<<"3"<<endl;
#endif
	c1.Clear();
	//angtruevsreco.Scale(1/angtruevsreco.Integral());
	angtruevsreco.Draw("colz");
	Preliminary("Simulation");
	c1.SaveAs(TString::Format("%s/angle_reco_vs_true.png",outdir));
	
#ifdef VERBOSE
	cout<<"4"<<endl;
#endif
	c1.Clear();
	//q2truevsreco.Scale(1/q2truevsreco.Integral());
	q2truevsreco.Draw("colz");
	Preliminary("Simulation");
	c1.SaveAs(TString::Format("%s/reconstructed_q2_reco_vs_true.png",outdir));
	
#ifdef VERBOSE
	cout<<"Saving PNGs of Fractional Error vs True histos"<<endl;
#endif
	// fractional error vs true parameter value
	c1.Clear();
	//nuEfracerrhist2.Scale(1/nuEfracerrhist2.Integral());
	nuEfracerrhist2.Draw("colz");
	Preliminary("Simulation");
	c1.SaveAs(TString::Format("%s/reconstructed_nue_error_vs_true.png",outdir));
	
#ifdef VERBOSE
	cout<<"2"<<endl;
#endif
	c1.Clear();
	//muEfracerrhist2.Scale(1/muEfracerrhist2.Integral());
	muEfracerrhist2.Draw("colz");
	Preliminary("Simulation");
	c1.SaveAs(TString::Format("%s/reconstructed_mue_error_vs_true.png",outdir));
	
#ifdef VERBOSE
	cout<<"3"<<endl;
#endif
	c1.Clear();
	//angfracerrhist2.Scale(1/angfracerrhist2.Integral());
	angfracerrhist2.Draw("colz");
	Preliminary("Simulation");
	c1.SaveAs(TString::Format("%s/reconstructed_angle_error_vs_true.png",outdir));
	
#ifdef VERBOSE
	cout<<"4"<<endl;
#endif
	c1.Clear();
	//q2fracerrhist2.Scale(1/q2fracerrhist2.Integral());
	q2fracerrhist2.Draw("colz");
	Preliminary("Simulation");
	c1.SaveAs(TString::Format("%s/reconstructed_q2_error_vs_true.png",outdir));
	
//	// loop over trueQevertexinfo file, because draw is acting weird.
//	for(int i=0; i<t3->GetEntries();i++){
//		mrdtlb->GetEntry(i);
//		cout<<"("<<mrdtracklen<<", "
//		cout<<"("<<muE;
//		//if(mrdtracklen>0)
//		ahist.Fill(muE*1000);
//		cout<<", "<<i<<") ";
//	}
	
	treeout->Write();
	
#ifdef VERBOSE
	cout<<"Cleaning up"<<endl;
#endif
	cout<<"cleaning up t1"<<endl;
	t1->ResetBranchAddresses();
	f1->Close();
	cout<<"cleaning up t2"<<endl;
	t2->ResetBranchAddresses();
	f2->Close();
	cout<<"cleaning up t3"<<endl;
	t3->ResetBranchAddresses();
	f3->Close();
	cout<<"cleaning up treeout"<<endl;
	treeout->ResetBranchAddresses();
	fileout->Close();
	delete fileout;
}

double CalculateEventQ2(double recoMuonEnergy, double recoNeutrinoEnergy, double recoMuonAngle){
	TDatabasePDG db;
	Double_t neutronmass = (db.GetParticle(2112)->Mass())*1000.;  // converted to MeV
	Double_t protonmass = (db.GetParticle(2212)->Mass())*1000.;   // converted to MeV
	Double_t muonmass = (db.GetParticle(13)->Mass())*1000.;       // converted to MeV
	Double_t O16bindingEnergy = 7.9762086875;                     // MeV (per nucleon), from http://tinyurl.com/y8m9s4z6
	Double_t boundneutronmass = neutronmass-O16bindingEnergy;

	double part1 = recoMuonEnergy - (sqrt(pow(recoMuonEnergy,2.)-pow(muonmass,2.))*TMath::Cos(recoMuonAngle));
	double eventq2 = -pow(muonmass,2.) + 2.*recoNeutrinoEnergy*part1;
	return eventq2;
}
