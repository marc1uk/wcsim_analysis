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
#include "TEveLine.h"	// << this ALSO needs to be included by gROOT->ProcessLine BEFORE compiling
#include "TStyle.h"
#include "TGeoManager.h"
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

#ifndef TESTRACKSVERBOSE
#define TESTRACKSVERBOSE 1
#endif

void checkmrdtracks(){

//	gSystem->IgnoreSignal(kSigSegmentationViolation, true); ???!

	TChain* c = new TChain("mrdtree");
	//c->Add("/pnfs/annie/persistent/users/moflaher/wcsim_tankonly_17-06-17_ana/mrdtrackfile.1000.root");
	TString pwd = gSystem->Getenv("PWD");
	TString chainpattern = pwd+"/mrdtrackfile.1000.root";
	c->Add(chainpattern.Data());
	int numents = c->GetEntries();
	cout<<"loaded "<<numents<<" entries in the chain"<<endl;
	if(numents<1){return;}
	c->LoadTree(0);
	
	int numsubevs;
	TBranch* numsubevsb=0;
	int numtracksinev;
	TBranch* numtracksinevb=0;
	TClonesArray* thesubeventarray = new TClonesArray("cMRDSubEvent");  // string is class name
	TBranch* thesubevsb=0;
	//TClonesArray* theTrackArray = new TClonesArray("cMRDTrack");      // this is a class member
	
	c->SetBranchAddress("nummrdsubeventsthisevent",&numsubevs,&numsubevsb);
	c->SetBranchAddress("nummrdtracksthisevent",&numtracksinev,&numtracksinevb);
	c->SetBranchAddress("subeventsinthisevent",&thesubeventarray,&thesubevsb);
	if(numsubevsb==0||numtracksinevb==0||thesubevsb==0){
		cerr<<"numsubevsb="<<numsubevsb<<", numtracksinevb="<<numtracksinevb<<", thesubevsb="<<thesubevsb<<endl;
		cerr<<"zombie branch, exiting."<<endl;
		return;
	}
	thesubevsb->SetAutoDelete(kFALSE);
	
	int totnumtracks=0;
	int numstopped=0;
	int numpenetrated=0;
	int numsideexit=0;
	int numtankintercepts=0;
	int numtankmisses=0;
	
	TH1F* hnumsubevs = new TH1F("numsubevs","Number of MRD SubEvents",20,0,10);
	TH1F* hnumtracks = new TH1F("numtracks","Number of MRD Tracks",10,0,10);
	
	TH1F* hrun = new TH1F("run","Run Number Histogram",20,0,9);
	TH1F* hevent = new TH1F("event","Event Number Histogram",100,0,2000);
	TH1F* hmrdsubev = new TH1F("subevent","MRD SubEvent Number Histogram",20,0,19);
	TH1F* htrigger = new TH1F("trigg","Trigger Number Histogram",10,0,9);
	TH1F* hnumhclusters = new TH1F("numhclusts","Number of H Clusters in MRD Track",20,0,10);
	TH1F* hnumvclusters = new TH1F("numvclusts","Number of V Clusters in MRD Track",20,0,10);
	TH1F* hnumhcells = new TH1F("numhcells","Number of H Cells in MRD Track",20,0,10);
	TH1F* hnumvcells = new TH1F("numvcells","Number of V Cells in MRD Track",20,0,10);
	TH1F* hpaddleids = new TH1F("paddleids","IDs of MRD Paddles Hit",500,0,499);
	TH1F* hpaddleinlayeridsh = new TH1F("inlayerpaddleidsh","IDs of MRD Paddles Hit Within H Layer",60,0,30);
	TH1F* hpaddleinlayeridsv = new TH1F("inlayerpaddleidsv","IDs of MRD Paddles Hit Within V Layer",60,0,30);
	TH1D* hdigittimes = new TH1D("digittimes","Digit Times",500,900,1300);
	
	TH1F* hhangle = new TH1F("hangle","Track Angle in Top View",100,-TMath::Pi(),TMath::Pi());
	TH1F* hhangleerr = new TH1F("hangleerr","Error in Track Angle in Top View",100,0,TMath::Pi());
	TH1F* hvangle = new TH1F("vangle","Track Angle in Side View",100,-TMath::Pi(),TMath::Pi());
	TH1F* hvangleerr = new TH1F("vangleerr","Error in Track Angle in Side View",100,0,TMath::Pi());
	TH1F* htotangle = new TH1F("totangle","Track Angle from Beam Axis",100,0,TMath::Pi());
	TH1F* htotangleerr = new TH1F("totangleerr","Error in Track Angle from Beam Axis",100,0,TMath::Pi());
	TH1F* henergyloss = new TH1F("energyloss","Track Energy Loss in MRD",100,0,1000);
	TH1F* henergylosserr = new TH1F("energylosserr","Error in Track Energy Loss in MRD",100,0,2000);
	TH1F* htracklength = new TH1F("tracklen","Total Track Length in MRD",100,0,220);
	TH1F* htrackpen = new TH1F("trackpen","Track Penetration in MRD",100,0,200);
	TH2F* htrackpenvseloss = new TH2F("trackpenvseloss","Track Penetration vs E Loss",100,0,220,100,0,1000);
	TH2F* htracklenvseloss = new TH2F("tracklenvseloss","Track Length vs E Loss",100,0,220,100,0,1000);
	
	TH3D* htrackstart = new TH3D("trackstart","MRD Track Start Vertices",100,-170,170,100,300,480,100,-230,220);
	TH3D* htrackstop = new TH3D("trackstop","MRD Track Stop Vertices",100,-170,170,100,300,480,100,-230,220);
	TH3D* hpep = new TH3D("pep","Back Projected Tank Exit",100,-500,500,100,0,480,100,-330,320);
	
	// Import the gdml geometry for the detector:
	TString gdmlpath = pwd+"/annie_v04_aligned.gdml";
	TGeoManager::Import(gdmlpath);
	// make all the materials semi-transparent so we can see the tracks going through them:
	TList* matlist = gGeoManager->GetListOfMaterials();
	TIter nextmaterial(matlist);
	while(TGeoMixture* amaterial=(TGeoMixture*)nextmaterial()) amaterial->SetTransparency(90);
	gGeoManager->GetVolume("WORLD2_LV")->Draw("ogl");
	TCanvas* gdmlcanv = (TCanvas*)gROOT->FindObject("c1");
	// maybe we can use these:
//   TEveViewer *ev = gEve->GetDefaultViewer();
//   TGLViewer  *gv = ev->GetGLViewer();
//   gv->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);

//   gEve->Redraw3D(kTRUE);
//   gSystem->ProcessEvents();

//   gv->CurrentCamera().RotateRad(-0.5, 1.4);
//   gv->RequestDraw();
	
	std::vector<TEveLine*> thiseventstracks;
	int numevents=c->GetEntries();
	// loop over events
	//numevents=1;
	int numtracksdrawn=0;
	bool earlyexit=false;
	for(int evi=0; evi<numevents; evi++){
		numsubevsb->GetEntry(evi);
		hnumsubevs->Fill(numsubevs);
		thesubeventarray->Clear();
		thesubevsb->GetEntry(evi);
		assert(thesubeventarray->GetEntriesFast()==numsubevs&&"Num subevents is different from TClonesArray size!");
		numtracksinevb->GetEntry(evi); // compare this later with sum of cMRDTracks found in all SubEvents
		
		int thetracki=-1;
		// loop over subevents (collections of hits on the MRD within a narrow time window)
		for(int subevi=0; subevi<numsubevs; subevi++){
			cMRDSubEvent* thesubevent = (cMRDSubEvent*)thesubeventarray->At(subevi);
			thesubevent->Print();
			
			std::vector<cMRDTrack>* tracksthissubevent = thesubevent->GetTracks();
			int numtracks = tracksthissubevent->size();
			hnumtracks->Fill(numtracks);
			totnumtracks+=numtracks;
			//cout<<"event "<<evi<<", subevent "<<subevi<<" had "<<numtracks<<" tracks"<<endl;
			// loop over tracks (collections of hits in a line within a subevent)
			for(auto&& thetrack : *tracksthissubevent){
				thetracki++;
				thetrack.Print2();  // finding tracks needs to be re-run with new class def
				/* the same as Print2.
				cout<<"wcsimfile: "<<thetrack.wcsimfile<<endl
				<<"run: "<<thetrack.run_id<<", event: "<<thetrack.event_id<<", trigger: "<<thetrack.trigger<<", subevent: "<<thetrack.mrdsubevent_id
				<<", track: "<<thetrack.MRDtrackID<<endl
				<<"num digits: "<<thetrack.digi_ids.size()<<endl
				<<"num pmts hit: "<<thetrack.pmts_hit.size()<<endl
				<<"digit times: ";
				for(auto atime : thetrack.digi_ts) cout<<atime<<", ";
			cout<<endl<<"digit charges: ";
				for(auto acharge : thetrack.digi_qs) cout<<acharge<<", ";
			cout<<endl<<"layers hit: ";
				for(auto alayer : thetrack.layers_hit) cout<<alayer<<", ";
			cout<<endl<<"energy deposited: ";
				for(auto aqdeposit : thetrack.eDepsInLayers) cout<<aqdeposit<<", ";
			cout<<endl<<"track ";
				if(thetrack.ispenetrating) cout<<"fully penetrates mrd"<<endl;
				if(thetrack.isstopped) cout<<"stops within mrd"<<endl;
				if(thetrack.sideexit) cout<<"exits side of mrd"<<endl;
			cout<<"side view clusters: "<<thetrack.htrackclusters.size()<<endl
				<<"top view clusters: "<<thetrack.vtrackclusters.size()<<endl
				<<"total track length: "<<thetrack.mutracklengthinMRD<<endl
				<<"penetration depth: "<<thetrack.penetrationdepth<<endl
				<<"energy loss: "<<thetrack.EnergyLoss<<", error: "<<thetrack.EnergyLossError<<endl
				<<"track start: ("<<thetrack.trackfitstart.X()<<", "<<thetrack.trackfitstart.Y()<<", "<<thetrack.trackfitstart.Z()
				<<"), track end: ("<<thetrack.trackfitstop.X()<<", "<<thetrack.trackfitstop.Y()<<", "<<thetrack.trackfitstop.Z()<<")"<<endl
				<<"h origin: "<<thetrack.htrackorigin<<", error: "<<thetrack.htrackoriginerror
				<<", gradient: "<<thetrack.htrackgradient<<", error: "<<thetrack.htrackgradienterror<<endl
				<<"v origin: "<<thetrack.vtrackorigin<<", error: "<<thetrack.vtrackoriginerror
				<<", gradient: "<<thetrack.vtrackgradient<<", error: "<<thetrack.vtrackgradienterror<<endl
				<<"fit chi2: "<<thetrack.htrackfitchi2<<", "<<thetrack.vtrackfitchi2<<endl
				<<"angle from z: "<<thetrack.trackangle<<", error: "<<thetrack.trackangleerror<<endl
				<<"back projection ";
				if(thetrack.interceptstank) cout<<"intercepts the tank"<<endl;
				else cout<<"does not intercept the tank"<<endl
				<<endl;
				*/
				
				if(thetrack.ispenetrating) numpenetrated++;
				else if(thetrack.isstopped) numstopped++;
				else numsideexit++;
				if(thetrack.interceptstank) numtankintercepts++;
				else numtankmisses++;
				
				hrun->Fill(thetrack.run_id);
				hevent->Fill(thetrack.event_id);
				hmrdsubev->Fill(thetrack.mrdsubevent_id);
				htrigger->Fill(thetrack.trigger);
				
				hnumhclusters->Fill(thetrack.htrackclusters.size());
				hnumvclusters->Fill(thetrack.vtrackclusters.size());
				for(auto&& acluster : thetrack.htrackclusters){
					hpaddleinlayeridsh->Fill(acluster.GetCentreIndex());
					Int_t uptubetopid = (2*acluster.xmaxid) + layeroffsets.at(acluster.layer);
					Int_t uptubebottomid = (2*acluster.xminid) + layeroffsets.at(acluster.layer);
					for(int i=uptubebottomid; i<=uptubebottomid; i++) hpaddleids->Fill(i);
					for(auto&& adigitime : acluster.digittimes) hdigittimes->Fill(adigitime);
				}
				for(auto&& acluster : thetrack.vtrackclusters){
					hpaddleinlayeridsv->Fill(acluster.GetCentreIndex());
					Int_t uptubetopid = (2*acluster.xmaxid) + layeroffsets.at(acluster.layer);
					Int_t uptubebottomid = (2*acluster.xminid) + layeroffsets.at(acluster.layer);
					for(int i=uptubebottomid; i<=uptubebottomid; i++) hpaddleids->Fill(i);
					for(auto&& adigitime : acluster.digittimes) hdigittimes->Fill(adigitime);
				}
				hnumhcells->Fill(thetrack.htrackcells.size());
				hnumvcells->Fill(thetrack.vtrackcells.size());
				
				hhangle->Fill(thetrack.htrackgradient);
				hhangleerr->Fill(thetrack.htrackgradienterror);
				hvangle->Fill(thetrack.htrackgradient);
				hvangleerr->Fill(thetrack.vtrackgradienterror);
				htotangle->Fill(thetrack.trackangle);
				htotangleerr->Fill(thetrack.trackangleerror);
				henergyloss->Fill(thetrack.EnergyLoss);
				henergylosserr->Fill(thetrack.EnergyLossError);
				htracklength->Fill(thetrack.mutracklengthinMRD);
				htrackpen->Fill(thetrack.penetrationdepth);
				htrackpenvseloss->Fill(thetrack.penetrationdepth,thetrack.EnergyLoss);
				htracklenvseloss->Fill(thetrack.mutracklengthinMRD,thetrack.EnergyLoss);
				
				TVector3* sttv = &thetrack.trackfitstart;
				TVector3* stpv = &thetrack.trackfitstop;
				TVector3* pep = &thetrack.projectedtankexitpoint;
				htrackstart->Fill(sttv->X(),sttv->Z(),sttv->Y());
				htrackstop->Fill(stpv->X(),stpv->Z(),stpv->Y());
				//cout<<"track "<<thetracki<<" started at ("<<sttv->X()<<", "<<sttv->Y()<<", "<<sttv->Z()<<")"
				//	<<" and ended at ("<<stpv->X()<<", "<<stpv->Y()<<", "<<stpv->Z()<<")"<<endl;
				
				TEveLine* evl = new TEveLine("track1",2);
				evl->SetLineWidth(4);
				evl->SetLineStyle(1);
				evl->SetMarkerColor(kRed);
				evl->SetRnrPoints(kTRUE);  // enable rendering of points
				//int icolour = thetracki;
				//while(icolour>=cMRDSubEvent::trackcolours.size()) icolour-=cMRDSubEvent::trackcolours.size();
				//evl->SetLineColor(cMRDSubEvent::trackcolours.at(icolour));
				evl->SetPoint(0,sttv->X(),sttv->Y(),sttv->Z());
				evl->SetPoint(1,stpv->X(),stpv->Y(),stpv->Z());
				if(!(pep->X()==0&&pep->Y()==0&&pep->Z()==0)){
					evl->SetPoint(2,pep->X(),pep->Y(),pep->Z());
					hpep->Fill(pep->X(),pep->Z(),pep->Y());
					cout<<"back projected point intercepts the tank at ("
						<<pep->X()<<", "<<pep->Y()<<", "<<pep->Z()<<")"<<endl;
					if(abs(pep->X())>450.) earlyexit=true;
				}
				thiseventstracks.push_back(evl);
				numtracksdrawn++;
				//cout<<"eveline constructed at "<<evl<<endl;
				//if(numtracksdrawn>100) earlyexit=true;
				
				
				/*
				// draw the track
				int trackcolourindex=thetrack.MRDtrackID+1; // element 0 is black
				while(trackcolourindex+1>=cMRDSubEvent::trackcolours.size()) 
					trackcolourindex-=cMRDSubEvent::trackcolours.size();
				EColor thistrackscolour = cMRDSubEvent::trackcolours.at(trackcolourindex);
				EColor fittrackscolour = cMRDSubEvent::trackcolours.at(trackcolourindex+1);
				thetrack->DrawReco(imgcanvas, trackarrows, thistrackscolour, paddlepointers);
				thetrack->DrawFit(imgcanvas, trackfitarrows, fittrackscolour);
				*/
				if(earlyexit) break;
			} // end loop over tracks
			
			gdmlcanv->cd();
			for(auto&& aline : thiseventstracks){
				//cout<<"drawing track at "<<aline<<" from ("<<aline->GetLineStart().fX
				//<<", "<<aline->GetLineStart().fY<<", "<<aline->GetLineStart().fZ<<") to ("
				//<<aline->GetLineEnd().fX<<", "<<aline->GetLineEnd().fY<<", "<<aline->GetLineEnd().fZ<<")"
				//<<endl;
				aline->Draw();
			}
			
			
			thesubevent->DrawMrdCanvases();  // creates the canvas with the digits
			thesubevent->DrawTracks();       // adds the CA tracks and their fit
			thesubevent->DrawTrueTracks();   // draws true tracks over the event
			thesubevent->imgcanvas->SaveAs(TString::Format("checkmrdtracks_%d.png",subevi+evi));
			
			//thesubevent->RemoveArrows();
			if(earlyexit) break;
		} // end loop over subevents
		gdmlcanv->Update();
		//std::this_thread::sleep_for (std::chrono::seconds(2));
		//for(auto aline : thiseventstracks) delete aline;
		//thiseventstracks.clear();
		if(earlyexit) break;
	} // end loop over events
	
	cout<<"Analysed "<<numevents<<" events, found "<<totnumtracks<<" MRD tracks, of which "
		<<numstopped<<" stopped in the MRD, "<<numpenetrated<<" fully penetrated and the remaining "
		<<numsideexit<<" exited the side."<<endl
		<<"Back projection suggests "<<numtankintercepts<<" tracks would have originated in the tank, while "
		<<numtankmisses<<" would not intercept the tank through back projection."<<endl;
	
	TCanvas c2;
	hnumsubevs->Draw();
	TCanvas c3;
	hnumtracks->Draw();
	TCanvas c4;
	hrun->Draw();
	TCanvas c5;
	hevent->Draw();
	TCanvas c6;
	hmrdsubev->Draw();
	TCanvas c7;
	htrigger->Draw();
	TCanvas c8;
	hnumhclusters->Draw();
	TCanvas c9;
	hnumvclusters->Draw();
	TCanvas c10;
	hnumhcells->Draw();
	TCanvas c11;
	hnumvcells->Draw();
	TCanvas c12;
	hpaddleids->Draw();
	TCanvas c13;
	hdigittimes->Draw();
	TCanvas c14;
	hhangle->Draw();
	TCanvas c15;
	hhangleerr->Draw();
	TCanvas c16;
	hvangle->Draw();
	TCanvas c17;
	hvangleerr->Draw();
	TCanvas c18;
	htotangle->Draw();
	TCanvas c19;
	htotangleerr->Draw();
	TCanvas c20;
	henergyloss->Draw();
	TCanvas c21;
	henergylosserr->Draw();
	TCanvas c22;
	htracklength->Draw();
	TCanvas c23;
	htrackpen->Draw();
	TCanvas c24;
	htrackpenvseloss->Draw();
	TCanvas c25;
	htracklenvseloss->Draw();
	TCanvas c26;
	htrackstart->SetMarkerStyle(20);
	htrackstart->SetMarkerColor(kRed);
	htrackstart->Draw();
	TCanvas c27;
	htrackstop->SetMarkerStyle(20);
	htrackstop->SetMarkerColor(kBlue);
	htrackstop->Draw();
	TCanvas c28;
	hpep->Draw();
	
	gPad->WaitPrimitive();
	
} // end

