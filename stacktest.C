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
#include <string>
#include <iostream>

#ifndef VERBOSE
#define VERBOSE 1
#endif

void stacktest(){
	const char* outdir="/tmp";
	
	// histograms to compare distributions
	TH1D truenuKEhist("truenuKEhist","True Neutrino KE",50,0,2000);
	truenuKEhist.FillRandom("gaus",1000);
	TH1D truemuKEhist("truemuKEhist","True Muon KE",50,0,2000);
	truemuKEhist.FillRandom("landau",1000);
	
	TH1D reconuKEhist("reconuKEhist","Reco Neutrino KE",50,0,2000);
	reconuKEhist.FillRandom("gaus",1000);
	TH1D recomuKEhist("recomuKEhist","Reco Muon KE",50,0,2000);
	recomuKEhist.FillRandom("landau",1000);
	
	// Save PNGs of the true and reco histograms overlaid
	TCanvas c1;
	THStack* thestack = new THStack("thestack","Stacked Histos");
	thestack->Add(&truenuKEhist);
	thestack->Add(&reconuKEhist);
	thestack->Draw("nostack");
	thestack->GetXaxis()->SetTitle("Neutrino Energy [MeV]");
	thestack->GetYaxis()->SetTitle("Num Events");
	c1.SaveAs(TString::Format("%s/reconstructed_neutrino_energy.png",outdir));
	
#ifdef VERBOSE
	cout<<"2"<<endl;
#endif
	delete thestack; thestack = new THStack("thestack","Stacked Histos");
	thestack->Add(&truemuKEhist);
	thestack->Add(&recomuKEhist);
	//thestack->GetXaxis()->SetTitle("Muon Energy [MeV]");
	//thestack->GetYaxis()->SetTitle("Num Events");
	thestack->Draw("nostack");
	c1.SaveAs(TString::Format("%s/reconstructed_muon_energy.png",outdir));
	
}
