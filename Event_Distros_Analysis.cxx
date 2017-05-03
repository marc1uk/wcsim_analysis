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
	t->SetBranchAddress("TotalChargeFromMuon",&totmuq);                // XXX // because parents wrong (fixd)
	t->SetBranchAddress("MuonTrackLengthInTank",&mutanklen);           // [0,500]
	t->SetBranchAddress("NuVtxInFidVol",&infidvol);
	t->SetBranchAddress("EventHasMuon",&hasmu);
	t->SetBranchAddress("MuonEntersMRD",&muenters);
	t->SetBranchAddress("MuonStopsInMRD",&mustops);
	t->SetBranchAddress("MuonRangesOutMRD",&muranges);
	
	
	TH1F* hist1 = new TH1F("hist1","hist1",100,0,7);
	TH1F* hist2 = new TH1F("hist2","hist2",100,0,9);
	TH1F* hist3 = new TH1F("hist3","hist3",100,0,7);
	TH1F* hist4 = new TH1F("hist4","hist4",100,0,9);
	TH1F* hist5 = new TH1F("hist5","hist5",100,0,7);
	TH1F* hist6 = new TH1F("hist6","hist6",100,0,9);
	TH1F* hist7 = new TH1F("hist7","hist7",100,0,7);
	TH1F* hist8 = new TH1F("hist8","hist8",100,0,9);
//	TH2F* hist2 = new TH2F("hist2","hist2",12,0,12,100,0,5000);
//	TH2F* hist3 = new TH2F("hist3","hist3",100,1,500,100,0,5000);
//	TH2F* hist4 = new TH2F("hist4","hist4",12,0,12,100,1,150);
	std::string et="IsQuasiElastic";
	std::string et2="IsWeakCC";    // everything is CC. 
//	int numCCQEstoppingmus=0, numCCQEnonstoppingmus=0, numCCotherstoppingmus=0, numCCothernonstoppingmus=0, numstoppingmus=0, numnonstoppingmus=0;
	
	for(int i=0; i<t->GetEntries();i++){
		t->GetEntry(i);
//		cout<<"Event "<<i<<" is a (";
//		eventtypes.at(et2) ? cout<<"CC" : cout<<"NC";
//		(eventtypes.at(et)) ? cout<<"QE" : cout<<"Other";
//		for(std::map<std::string,bool>::iterator it=eventtypes.begin(); it!=eventtypes.end(); it++)
//			cout<<it->second<<", ";
//		cout<<"nuE="<<nuE<<", Q2="<<q2<<", muangle="<<muangle;
//		cout<<") event"<<endl;
		//cout<<"muE="<<mustartE<<", mupencm="<<mupencm<<", mupenlayers="<<mupenlayers<<endl;
/*		if(mustops) numstoppingmus++;
		else numnonstoppingmus++;
		if((!(eventtypes.at(et)))&&eventtypes.at(et2)){
			if(mustops) numCCotherstoppingmus++;
			else numCCothernonstoppingmus++;
		}
*/
		hist1->Fill(nuE);	// all events
		hist2->Fill(q2);
		//if(!(((eventtypes.at(et)))&&eventtypes.at(et2))) continue;
		if(!muenters) continue;
//		if(mustops) numCCQEstoppingmus++;
//		else numCCQEnonstoppingmus++;
		hist3->Fill(nuE);	// all events with a muon that enters the MRD
		hist4->Fill(q2);	// regardless of stopping, rangeout or side exit
		if(muranges){
			hist5->Fill(nuE);	// all events with a muon that ranges out the MRD
			hist6->Fill(q2);
		} else if(mustops){
			hist7->Fill(nuE);	// all events with a muon that stops witin the MRD
			hist8->Fill(q2);
		}	// remainder are particles that exit the sides of the MRD
//		hist2->Fill(mupenlayers,mustartE);
//		hist3->Fill(mutotcm,mustartE);
//		hist4->Fill(mupenlayers,mupencm);
	}
//	cout<<"there were "<<numCCQEstoppingmus<<"CCQE events with stopped muons vs "
//	    <<numCCQEnonstoppingmus<<" non stopping muons"<<endl;
//	cout<<"there were "<<numCCotherstoppingmus<<" CC-Other events with stopped muons vs "
//	    <<numCCothernonstoppingmus<<" non stopping muons"<<endl;
//	cout<<"there were "<<numstoppingmus<<" total stopped muons vs "<<numnonstoppingmus
//	    <<" non-stopping muons"<<endl;
	TCanvas* c1=new TCanvas("c1","c1"); c1->cd();
	hist1->SetLineColor(kRed);
	hist1->Draw();
	hist3->SetLineColor(kBlue);
	hist3->Draw("same");
	hist5->SetLineColor(kGreen);
	hist5->Draw("same");
	hist7->SetLineColor(kMagenta);
	hist7->Draw("same");
	TCanvas* c2=new TCanvas("c2","c2"); c2->cd();
	hist2->SetLineColor(kRed);
	hist2->Draw();
	hist4->SetLineColor(kBlue);
	hist4->Draw("same");
	hist6->SetLineColor(kGreen);
	hist6->Draw("same");
	hist8->SetLineColor(kMagenta);
	hist8->Draw("same");
}
