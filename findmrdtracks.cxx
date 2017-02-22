/* vim:set noexpandtab tabstop=4 wrap */
// #######################################################################
#include "MRDTrackClass.hh"				// a class for defining MRD tracks

#ifndef MRDTRACKVERBOSE
//#define MRDTRACKVERBOSE 1
#endif

// CREATE+OPEN OUTPUT FILE
// =======================
void WCSimAnalysis::OpenMRDtrackOutfile(){
	cout<<"opening mrd output file"<<endl;
	mrdtrackfile = new TFile("mrdtrackfile.root","RECREATE","MRD Tracks file");
	mrdtrackfile->cd();
	recotree = new TTree("recotree","Tree for reconstruction data");
	gROOT->cd();
	
	nummrdtracksthisevent=0;
	aTrack = new TClonesArray("cMRDTrack");	// string is class name
	
	nummrdtracksthiseventb = recotree->Branch("nummrdtracksthisevent",&nummrdtracksthisevent);
	tracksinthiseventb = recotree->Branch("tracksinthisevent",&aTrack, MAXTRACKSPEREVENT);
	
}

// #######################################################################

// SPLIT HITS BY TIME
// ==================
void WCSimAnalysis::FindMRDtracksInEvent(){
#ifdef MRDTRACKVERBOSE
	cout<<"eventnum is "<<eventnum<<endl;
#endif
	//aTracka.Clear();
/* 
if your class contains pointers, use aTrack.Clear("C"). You MUST then provide a Clear() method in your class that properly performs clearing and memory freeing. (or "implements the reset procedure for pointer objects")
 see https://root.cern.ch/doc/master/classTClonesArray.html#a025645e1e80ea79b43a08536c763cae2
*/
	
	if(mrddigittimesthisevent.size()<4){
		// =====================================================================================
		// NO TRACKS IN THIS EVENT
		// ======================================
		new((*aTrack)[0]) cMRDTrack();
		nummrdtracksthisevent=0;
		nummrdtracksthiseventb->Fill();
		tracksinthiseventb->Fill();
		//recotree->Fill();						// fill even though there are not tracks so the entries align.
		return;
		// skip remainder
		// ======================================================================================
	}
	
	// MEASURE EVENT DURATION TO DETERMINE IF THERE IS MORE THAN ONE TRACK
	// ===================================================================
	std::vector<Int_t> digitidsinatrack;
	std::vector<Int_t> tubeidsinatrack;
	std::vector<Double_t> digitqsinatrack;
	std::vector<Double_t> digittimesinatrack;
	std::vector<Int_t> digitnumtruephots;
	std::vector<Int_t> particleidsinatrack;
	std::vector<Double_t> photontimesinatrack;
	
	// first check: are all hits within a 30ns window (maxtrackduration) If so, just one track. 
#ifdef MRDTRACKVERBOSE
	cout<<"mrddigittimesthisevent.size()="<<mrddigittimesthisevent.size()<<endl;
#endif
	Double_t eventendtime = *std::max_element(mrddigittimesthisevent.begin(),mrddigittimesthisevent.end());
	Double_t eventstarttime = *std::min_element(mrddigittimesthisevent.begin(),mrddigittimesthisevent.end());
	Double_t eventduration = (eventendtime - eventstarttime);
#ifdef MRDTRACKVERBOSE
	cout<<"event start: "<<eventstarttime<<", end : "<<eventendtime<<", duration : "<<eventduration<<endl;
#endif
	
	if(eventduration<maxtrackduration){
	// JUST ONE TRACK
	// ==============
#ifdef MRDTRACKVERBOSE
		cout<<"all hits this event within one track."<<endl;
#endif
		
		// loop over digits and convert them into cMRDdigit objects that contain all their information
		for(Int_t thisdigit=0;thisdigit<mrddigittimesthisevent.size();thisdigit++){
			digitidsinatrack.push_back(thisdigit);
			WCSimRootCherenkovDigiHit* thedigihit = (WCSimRootCherenkovDigiHit*)atrigm->GetCherenkovDigiHits()->At(thisdigit);
			Int_t thisdigitstubeid = thedigihit->GetTubeId();
			tubeidsinatrack.push_back(thisdigitstubeid);
			Int_t thisdigitsq = thedigihit->GetQ();
			digitqsinatrack.push_back(thisdigitsq);
			double thisdigitstime = thedigihit->GetT();
			digittimesinatrack.push_back(thisdigitstime);
			// add all the unique parent ID's for digits contributing to this track (truth level info)
			std::vector<int> truephotonindices = thedigihit->GetPhotonIds();
			digitnumtruephots.push_back(truephotonindices.size());
			for(int truephoton=0; truephoton<truephotonindices.size(); truephoton++){
				int thephotonsid = truephotonindices.at(truephoton);
				WCSimRootCherenkovHitTime *thehittimeobject = (WCSimRootCherenkovHitTime*)atrigm->GetCherenkovHitTimes()->At(thephotonsid);
				int thephotonsparenttrack = thehittimeobject->GetParentID();
				int checkcount = std::count(particleidsinatrack.begin(), particleidsinatrack.end(), thephotonsparenttrack);
				if(checkcount==0){ particleidsinatrack.push_back(thephotonsparenttrack); }
				double thephotonstruetime = thehittimeobject->GetTruetime();
				photontimesinatrack.push_back(thephotonstruetime);
			}
			// append the digit - nah, just construct track with all digits
			//aTracka.AppendDigit(thisdigitstime, thisdigitsq, thisdigitstubeid,
			//											photontimesinatrack,particleidsinatrack);
		}
		
		// construct the track from all the digits
		if(digitidsinatrack.size()>3){
#ifdef MRDTRACKVERBOSE
			cout<<"constructing a single track for this event"<<endl;
#endif
			new((*aTrack)[0]) cMRDTrack(0, currentfilestring, runnum, eventnum, subtriggernum, digitidsinatrack, tubeidsinatrack, digitqsinatrack, digittimesinatrack, digitnumtruephots, photontimesinatrack, particleidsinatrack);
		}
		// can also use 'cMRDTrack* = (cMRDTrack*)aTrack.ConstructedAt(0);' followed by a bunch of 'Set' calls
		// to set all relevant fields. This bypasses the constructor, calling it only when necessary, 
		// saving time. In that case, we do not need to call aTracka.Clear();
		
		nummrdtracksthisevent=1;
		nummrdtracksthiseventb->Fill();
		tracksinthiseventb->Fill();
		//recotree->Fill();
		
	} else {
	// MORE THAN ONE TRACK
	// ===================
		// this event has multiple tracks. Need to split hits into which track they belong to.
		// scan over the times and look for gaps where no digits lie, using these to delimit 'tracks'
		std::vector<Float_t> trackhittimesv;	// a vector of the starting times of a given 'track'
		std::vector<double> sorteddigittimes(mrddigittimesthisevent);
		std::sort(sorteddigittimes.begin(), sorteddigittimes.end());
		trackhittimesv.push_back(sorteddigittimes.at(0));
		for(Int_t i=0;i<sorteddigittimes.size()-1;i++){
			Float_t timetonextdigit = sorteddigittimes.at(i+1)-sorteddigittimes.at(i);
			if(timetonextdigit>maxtrackduration){
				trackhittimesv.push_back(sorteddigittimes.at(i+1));
#ifdef MRDTRACKVERBOSE
				cout<<"Setting track time threshold at "<<trackhittimesv.back()<<endl;
//				cout<<"this digit is at "<<sorteddigittimes.at(i)<<endl;
//				cout<<"next digit is at "<<sorteddigittimes.at(i+1)<<endl;
//				try{
//					cout<<"next next digit is at "<<sorteddigittimes.at(i+2)<<endl;
//				} catch (...) { int i=0; }
#endif
			}
		}
#ifdef MRDTRACKVERBOSE
		cout<<trackhittimesv.size()<<" tracks this event"<<endl;
#endif
		
		// a vector to record the track number for each hit, to know if we've allocated it yet.
		std::vector<Int_t> tracknumthisevent(mrddigittimesthisevent.size(),-1);
		
		// now we need to sort the digits into the tracks they belong to:
		// loop over tracks
		Int_t actualtrackcounter=0;	// not all time groups will have enough digits
		for(Int_t thistrack=0; thistrack<trackhittimesv.size(); thistrack++){
#ifdef MRDTRACKVERBOSE
			cout<<"Digits in MRD at = "<<trackhittimesv.at(thistrack)<<"ns in event "<<eventnum<<endl;
#endif
			// don't need to worry about lower bound as we start from lowest t peak and 
			// exclude already allocated hits
			
			Float_t endtime = (thistrack<(trackhittimesv.size()-1)) ? trackhittimesv.at(thistrack+1) : (eventendtime+1.);
#ifdef MRDTRACKVERBOSE
			cout<<"endtime for track "<<thistrack<<" is "<<endtime<<endl;
#endif
			// times are not ordered, so scan through all digits for each track
			for(Int_t thisdigit=0;thisdigit<mrddigittimesthisevent.size();thisdigit++){
				if(tracknumthisevent.at(thisdigit)<0 && mrddigittimesthisevent.at(thisdigit)< endtime ){
					// thisdigit is in thistrack
#ifdef MRDTRACKVERBOSE
					cout<<"adding digit at "<<mrddigittimesthisevent.at(thisdigit)<<" to track "<<thistrack<<endl;
#endif
					digitidsinatrack.push_back(thisdigit);
					tracknumthisevent.at(thisdigit)=thistrack;
					WCSimRootCherenkovDigiHit* thedigihit = (WCSimRootCherenkovDigiHit*)atrigm->GetCherenkovDigiHits()->At(thisdigit);
					Int_t thisdigitstubeid = thedigihit->GetTubeId();
					tubeidsinatrack.push_back(thisdigitstubeid);
					Int_t thisdigitsq = thedigihit->GetQ();
					digitqsinatrack.push_back(thisdigitsq);
					double thisdigitstime = thedigihit->GetT();
					digittimesinatrack.push_back(thisdigitstime);
					// add all the unique parent ID's for digits contributing to this track (truth level info)
					std::vector<int> truephotonindices = thedigihit->GetPhotonIds();
					digitnumtruephots.push_back(truephotonindices.size());
					for(int truephoton=0; truephoton<truephotonindices.size(); truephoton++){
						int thephotonsid = truephotonindices.at(truephoton);
						WCSimRootCherenkovHitTime *thehittimeobject = (WCSimRootCherenkovHitTime*)atrigm->GetCherenkovHitTimes()->At(thephotonsid);
						int thephotonsparenttrack = thehittimeobject->GetParentID();
						int checkcount = std::count(particleidsinatrack.begin(), particleidsinatrack.end(), thephotonsparenttrack);
						if(checkcount==0){ particleidsinatrack.push_back(thephotonsparenttrack); }
						double thephotonstruetime = thehittimeobject->GetTruetime();
						photontimesinatrack.push_back(thephotonstruetime);
					}
				}
			}
			
			// construct the track from all the digits
			if(digitidsinatrack.size()>3){	// must have enough for a track
#ifdef MRDTRACKVERBOSE
				cout<<"constructing track "<<actualtrackcounter<<" with "<<digitidsinatrack.size()<<" digits for this event"<<endl;
#endif
				new((*aTrack)[actualtrackcounter]) cMRDTrack(actualtrackcounter, currentfilestring, runnum, eventnum, subtriggernum, digitidsinatrack, tubeidsinatrack, digitqsinatrack, digittimesinatrack, digitnumtruephots, photontimesinatrack, particleidsinatrack);
				actualtrackcounter++;
			}
			// can also use 'cMRDTrack* = (cMRDTrack*)aTrack.ConstructedAt(0);' followed by a bunch of 'Set' calls
			// to set all relevant fields. This bypasses the constructor, calling it only when necessary, 
			// saving time. In that case, we do not need to call aTracka.Clear();
			
			// clear the vectors and loop to the next track
			digitidsinatrack.clear();
			tubeidsinatrack.clear();
			digitqsinatrack.clear();
			digittimesinatrack.clear();
			particleidsinatrack.clear();
			photontimesinatrack.clear();
			digitnumtruephots.clear();
			
		}
		
		// quick scan to check for any unallocated hits
		for(Int_t k=0;k<tracknumthisevent.size();k++){
			if(tracknumthisevent.at(k)==-1){cout<<"*****unbinned hit!"<<k<<" "<<mrddigittimesthisevent.at(k)<<endl;}
		}
		
		nummrdtracksthisevent=trackhittimesv.size();
		nummrdtracksthiseventb->Fill();
		tracksinthiseventb->Fill();
		//recotree->Fill();
		
	}	// end multiple tracks case
	
	// WRITE+CLOSE OUTPUT FILES
	// ========================
	mrdtrackfile->cd();
	recotree->SetEntries(nummrdtracksthiseventb->GetEntries());
	recotree->Write("",kOverwrite);
	gROOT->cd();
}

// #######################################################################
