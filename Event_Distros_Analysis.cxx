// call with root /home/marc/LinuxSystemFiles/WCSim/gitver/root_work/Event_Distros_Analysis.cxx+
// note: NEED THE +
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

using namespace std;

void Event_Distros_Analysis(const char* file="/home/marc/anniegpvm/wcsim/root_work/EventDistributions.root"){
	//disable stats box:
	gStyle->SetOptStat(0);
	//disable histogram title
	gStyle->SetOptTitle(0);
	
	TFile* f=TFile::Open(file);
	TTree* t = (TTree*)f->Get("treeout");
	std::map<std::string,bool> eventtypes;
	std::map<std::string,bool>* eventtypesp=&eventtypes;
	double q2, nuE, mustartE, muangle, mupencm, mutotcm, chargeincone, totmuq, mutanklen;
	int mupenlayers;
	bool infidvol, hasmu, muenters, mustops, muranges;
	t->SetBranchAddress("TypesMap",&eventtypesp);
	t->SetBranchAddress("EventQ2", &q2);                               // [0,9]
	t->SetBranchAddress("NeutrinoEnergy",&nuE);                        // [0,7]
	t->SetBranchAddress("MuonStartEnergy",&mustartE);                  // [0,5000]
	t->SetBranchAddress("MuonAngle",&muangle);                         // [0,Pi]
	t->SetBranchAddress("MuonMrdPenetrationInCm",&mupencm);            // [1,150]
	t->SetBranchAddress("MuonMrdPenetrationLayers",&mupenlayers);      // [0,12]
	t->SetBranchAddress("MuonTrackLengthInMRD",&mutotcm);              // [1,500]
	t->SetBranchAddress("FractionOfMuonChargeInCone",&chargeincone);   // [0,1]
	t->SetBranchAddress("TankChargeFromMuon",&totmuq);                 // XXX // because parents wrong (fixd)
	t->SetBranchAddress("MuonTrackLengthInTank",&mutanklen);           // [0,500]
	t->SetBranchAddress("NuVtxInFidVol",&infidvol);
	t->SetBranchAddress("EventHasMuon",&hasmu);
	t->SetBranchAddress("MuonEntersMRD",&muenters);
	t->SetBranchAddress("MuonStopsInMRD",&mustops);
	t->SetBranchAddress("MuonRangesOutMRD",&muranges);
	
	
	TH1F* enuall = new TH1F("hist1","All Events",100,0,7);
	TH1F* enumuenters = new TH1F("hist2","Fiducial Events with Muon Entering MRD",100,0,9);
	TH1F* enumustops = new TH1F("hist3","Fiducial Events with Muon Stopping in MRD",100,0,7);
	TH1F* enumupenetrates = new TH1F("hist4","Fiducial Events with Muon Fully Penetrating MRD",100,0,9);
	TH1F* enusideexit = new TH1F("hist9","Fiducial Events with Muon Exiting Side of MRD",100,0,7);
	TH1F* q2all = new TH1F("hist5","All Events",100,0,3);
	TH1F* q2muenters = new TH1F("hist6","Fiducial Events with Muon Entering MRD",100,0,3);
	TH1F* q2mustops = new TH1F("hist7","Fiducial Events with Muon Stopping in MRD",100,0,3);
	TH1F* q2mupenetrates = new TH1F("hist8","Fiducial Events with Muon Fully Penetrating MRD",100,0,3);
	TH1F* q2sideexit = new TH1F("hist10","Fiducial Events with Muon Exiting Side of MRD",100,0,3);
	TH2F* penecmvsemu = new TH2F("hist11","Penetration Vs Muon Energy",100,1,150,100,0,3000);
	TH2F* penelayersvsemu = new TH2F("hist12","Penetration vs Muon Energy",12,0,12,100,0,3000);
	TH2F* mrdlenvsemu = new TH2F("hist13","MRD Track Length vs Muon Energy",100,0,500,100,1,3000);
	std::vector<TH1F*> thehistos{enuall,enumuenters,enumustops,enumupenetrates,enusideexit,q2all,q2muenters,q2mustops,q2mupenetrates,q2sideexit};
	std::string et="IsQuasiElastic";
	std::string et2="IsWeakCC";    // everything is CC. 
//	int numCCQEstoppingmus=0, numCCQEnonstoppingmus=0, numCCotherstoppingmus=0, numCCothernonstoppingmus=0, numstoppingmus=0, numnonstoppingmus=0;
	
	for(int i=0; i<t->GetEntries();i++){
		t->GetEntry(i);
//		cout<<"Event "<<i<<" is a (";
//		eventtypes.at(et2) ? cout<<"CC" : cout<<"NC";
//		(eventtypes.at(et)) ? cout<<"QE" : cout<<"Other";
//		cout<<"nuE="<<nuE<<", Q2="<<q2<<", muangle="<<muangle;
//		cout<<") event"<<endl;
		//cout<<"muE="<<mustartE<<", mupencm="<<mupencm<<", mupenlayers="<<mupenlayers<<endl;
//		if(isCC&&isQE){ int i=0; /* TODO: do we wanna do CCQE vs CC-Other? */ }
//		if(isCC&&!isQE){ int i=0; /* TODO: do we wanna do CCQE vs CC-Other? */ }
		
		enuall->Fill(nuE);	// all events
		q2all->Fill(q2);
		
		if(infidvol){
			if(muenters){
				enumuenters->Fill(nuE);
				q2muenters->Fill(q2);
			}
		
			if(muenters&&mustops){
				enumustops->Fill(nuE);
				q2mustops->Fill(q2);
			}
		
			if(muenters&&muranges){
				enumupenetrates->Fill(nuE);
				q2mupenetrates->Fill(q2);
			}
		
			if(muenters&&!mustops&&!mupenlayers){	// events that exit the side of the MRD
				enusideexit->Fill(nuE);
				q2sideexit->Fill(q2);
			}
			penelayersvsemu->Fill(mupenlayers,mustartE);
			penecmvsemu->Fill(mupencm,mustartE);
			mrdlenvsemu->Fill(mutotcm,mustartE);
		}
	}
	// entry, stopping, penetration and side-exit counts
	cout<<"there were "<<enuall->GetEntries()<<" total entries analysed"<<endl;
	cout<<q2muenters->GetEntries()<<" entries entered the MRD, "
	    <<enumustops->GetEntries()<<" entries stopped in the MRDmrd, and "
	    <<enumupenetrates->GetEntries()<<" entries fully penetrated the MRD, leaving "
	    <<enusideexit->GetEntries()<<" entries that exited the sides of the MRD"<<endl;
	
	// this sets the canvas title but prevents BuildLegend giving the correct titles
	for(auto ahist : thehistos){
		int thecount=ahist->GetEntries();
		ahist->Scale(1./thecount);
	}
	
	TCanvas* c1=new TCanvas("c1","c1"); c1->cd();
	enuall->SetLineColor(kRed);
	enuall->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
	enuall->GetXaxis()->SetLabelFont(42);
	enuall->GetXaxis()->SetTitleSize(0.05);
	enuall->GetXaxis()->SetTitleFont(42);
	enuall->GetYaxis()->SetTitle("Num Events");
	enuall->GetYaxis()->SetLabelFont(42);
	enuall->GetYaxis()->SetTitleSize(0.05);
	enuall->GetYaxis()->SetTitleFont(42);
	enuall->Draw();
//	enumuenters->SetLineColor(kBlue);
//	enumuenters->Draw("same");
	enumustops->SetLineColor(kGreen);
	enumustops->Draw("same");
//	enumupenetrates->SetLineColor(kMagenta);
//	enumupenetrates->Draw("same");
	c1->BuildLegend();
	TCanvas* c2=new TCanvas("c2","c2"); c2->cd();
	q2all->SetLineColor(kRed);
	q2all->GetXaxis()->SetTitle("Events Q^{2} [GeV/c^{2}]");
	q2all->GetYaxis()->SetTitle("Num Events");
	q2all->GetXaxis()->SetLabelFont(42);
	q2all->GetXaxis()->SetTitleSize(0.05);
	q2all->GetXaxis()->SetTitleFont(42);
	q2all->GetYaxis()->SetTitle("Num Events");
	q2all->GetYaxis()->SetLabelFont(42);
	q2all->GetYaxis()->SetTitleSize(0.05);
	q2all->GetYaxis()->SetTitleFont(42);
	q2all->Draw();
//	q2muenters->SetLineColor(kBlue);
//	q2muenters->Draw("same");
	q2mustops->SetLineColor(kGreen);
	q2mustops->Draw("same");
//	q2mupenetrates->SetLineColor(kMagenta);
//	q2mupenetrates->Draw("same");
	c2->BuildLegend();
	
//	
//	TCanvas* c3 = new TCanvas("c3","c3"); c3->cd();
//	penelayersvsemu->Draw("box");
//	TCanvas* c4 = new TCanvas("c4","c4"); c4->cd();
//	penecmvsemu->Draw("box");
//	TCanvas* c5 = new TCanvas("c5","c5"); c5->cd();
//	mrdlenvsemu->Draw("box");
	
}
