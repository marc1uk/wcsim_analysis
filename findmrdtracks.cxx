/* vim:set noexpandtab tabstop=4 wrap */
// #######################################################################
#include "MRDSubEventClass.hh"			// a class for defining subevents
#include "MRDTrackClass.hh"				// a class for defining MRD tracks

#ifndef MRDSPLITVERBOSE
//#define MRDSPLITVERBOSE 1
#endif

// CREATE+OPEN OUTPUT FILE
// =======================
void WCSimAnalysis::OpenMRDtrackOutfile(int filenum){
	TString filenameout = TString::Format("%s/mrdtrackfile.%d.root",outputdir,filenum);
	cout<<"opening mrd output file "<<filenameout.Data()<<endl;
	if(mrdtrackfile) mrdtrackfile->Close();
	mrdtrackfile = new TFile(filenameout.Data(),"RECREATE","MRD Tracks file");
	mrdtrackfile->cd();
	mrdtree = new TTree("mrdtree","Tree for reconstruction data");
	gROOT->cd();
	
	nummrdtracksthisevent=0;
	aSubEvent = new TClonesArray("cMRDSubEvent");	// string is class name
	
	mrdeventnumb = mrdtree->Branch("EventID",&eventnum);
	mrdtriggernumb = mrdtree->Branch("TriggerID",&triggernum);
	nummrdsubeventsthiseventb = mrdtree->Branch("nummrdsubeventsthisevent",&nummrdsubeventsthisevent);
	subeventsinthiseventb = mrdtree->Branch("subeventsinthisevent",&aSubEvent, nummrdsubeventsthisevent);
	nummrdtracksthiseventb = mrdtree->Branch("nummrdtracksthisevent",&nummrdtracksthisevent);
	
}

// #######################################################################

// SPLIT HITS BY TIME
// ==================
void WCSimAnalysis::FindMRDtracksInEvent(){

	const Int_t minimumdigits =4;
#ifdef MRDSPLITVERBOSE
	cout<<"Searching for MRD tracks in event "<<eventnum<<endl;
	cout<<"mrddigittimesthisevent.size()="<<mrddigittimesthisevent.size()<<endl;
#endif
	aSubEvent->Clear("C");
/* 
if your class contains pointers, use aTrack.Clear("C"). You MUST then provide a Clear() method in your class that properly performs clearing and memory freeing. (or "implements the reset procedure for pointer objects")
 see https://root.cern.ch/doc/master/classTClonesArray.html#a025645e1e80ea79b43a08536c763cae2
*/
	Int_t mrdeventcounter=0;
	if(mrddigittimesthisevent.size()<minimumdigits){
		// =====================================================================================
		// NO DIGITS IN THIS EVENT
		// ======================================
		//new((*aSubEvent)[0]) cMRDSubEvent();
		nummrdsubeventsthisevent=mrdeventcounter;
		nummrdtracksthisevent=0;
		nummrdsubeventsthiseventb->Fill();
		nummrdtracksthiseventb->Fill();
		subeventsinthiseventb->Fill();
		//mrdtree->Fill();						// fill the branches so the entries align.
		mrdtrackfile->cd();
		mrdtree->SetEntries(nummrdtracksthiseventb->GetEntries());
		mrdtree->Write("",kOverwrite);
		gROOT->cd();
		return;
		// skip remainder
		// ======================================================================================
	}
	
	// MEASURE EVENT DURATION TO DETERMINE IF THERE IS MORE THAN ONE MRD SUB EVENT
	// ===========================================================================
	std::vector<Int_t> digitidsinasubevent;
	std::vector<Int_t> tubeidsinasubevent;
	std::vector<Double_t> digitqsinasubevent;
	std::vector<Double_t> digittimesinasubevent;
	std::vector<Int_t> digitnumtruephots;
	std::vector<Int_t> particleidsinasubevent;
	std::vector<Double_t> photontimesinasubevent;
	
	// first check: are all hits within a 30ns window (maxsubeventduration) If so, just one subevent. 
#ifdef MRDSPLITVERBOSE
	cout<<"mrddigittimesthisevent.size()="<<mrddigittimesthisevent.size()<<endl;
#endif
	Double_t eventendtime = *std::max_element(mrddigittimesthisevent.begin(),mrddigittimesthisevent.end());
	Double_t eventstarttime = *std::min_element(mrddigittimesthisevent.begin(),mrddigittimesthisevent.end());
	Double_t eventduration = (eventendtime - eventstarttime);
#ifdef MRDSPLITVERBOSE
	cout<<"mrd event start: "<<eventstarttime<<", end : "<<eventendtime<<", duration : "<<eventduration<<endl;
#endif
	
	if((eventduration<maxsubeventduration)&&(mrddigittimesthisevent.size()>=minimumdigits)){
	// JUST ONE SUBEVENT
	// =================
#ifdef MRDSPLITVERBOSE
		cout<<"all hits this event within one subevent."<<endl;
#endif
		// loop over digits and extract info
		for(Int_t thisdigit=0;thisdigit<mrddigittimesthisevent.size();thisdigit++){
			digitidsinasubevent.push_back(thisdigit);
			WCSimRootCherenkovDigiHit* thedigihit =
				(WCSimRootCherenkovDigiHit*)atrigm->GetCherenkovDigiHits()->At(thisdigit);
			Int_t thisdigitstubeid = thedigihit->GetTubeId()-1;  //XXX IMPORTANT: Tube index = TubeId()-1
			tubeidsinasubevent.push_back(thisdigitstubeid);
			Int_t thisdigitsq = thedigihit->GetQ();
			digitqsinasubevent.push_back(thisdigitsq);
			double thisdigitstime = thedigihit->GetT();
			digittimesinasubevent.push_back(thisdigitstime);
			// add all the unique parent ID's for digits contributing to this subevent (truth level info)
			std::vector<int> truephotonindices = thedigihit->GetPhotonIds();
			digitnumtruephots.push_back(truephotonindices.size());
			for(int truephoton=0; truephoton<truephotonindices.size(); truephoton++){
				int thephotonsid = truephotonindices.at(truephoton);
				WCSimRootCherenkovHitTime *thehittimeobject = (WCSimRootCherenkovHitTime*)atrigm->GetCherenkovHitTimes()->At(thephotonsid);
				int thephotonsparentsubevent = thehittimeobject->GetParentID();
				particleidsinasubevent.push_back(thephotonsparentsubevent);
				double thephotonstruetime = thehittimeobject->GetTruetime();
				photontimesinasubevent.push_back(thephotonstruetime);
			}
			// append the digit - nah, just construct subevent with all digits
			//aSubEventa.AppendDigit(thisdigitstime, thisdigitsq, thisdigitstubeid,
			//										photontimesinasubevent,particleidsinasubevent);
		}
		
		// we store the truth tracks in each cMRDSubEvent, so we need to pass the ones within the subevent
		// time window to the constructor
		Int_t numtracks = atrigt->GetNtrack();
		std::vector<WCSimRootTrack*> truetrackpointers;
		for(int truetracki=0; truetracki<numtracks; truetracki++){
			WCSimRootTrack* nextrack = (WCSimRootTrack*)atrigt->GetTracks()->At(truetracki);
			if((nextrack->GetFlag()==0)/*&&(nextrack->GetIpnu()!=11)*/) truetrackpointers.push_back(nextrack);
		}
		// construct the subevent from all the digits
#ifdef MRDSPLITVERBOSE
		cout<<"constructing a single subevent for this event"<<endl;
#endif
		cMRDSubEvent* currentsubevent = new((*aSubEvent)[0]) cMRDSubEvent(0, currentfilestring, runnum, eventnum, triggernum, digitidsinasubevent, tubeidsinasubevent, digitqsinasubevent, digittimesinasubevent, digitnumtruephots, photontimesinasubevent, particleidsinasubevent, truetrackpointers);
		mrdeventcounter++;
		// can also use 'cMRDSubEvent* = (cMRDSubEvent*)aSubEvent.ConstructedAt(0);' followed by a bunch of
		// 'Set' calls to set all relevant fields. This bypasses the constructor, calling it only when 
		// necessary, saving time. In that case, we do not need to call aSubEventa.Clear();
		
		int tracksthissubevent=currentsubevent->GetTracks()->size();
#ifdef MRDSPLITVERBOSE
		cout<<"the only subevent this event found "<<tracksthissubevent<<" tracks"<<endl;
#endif
		nummrdsubeventsthisevent=1;
		nummrdtracksthisevent=tracksthissubevent;
		nummrdsubeventsthiseventb->Fill();
		nummrdtracksthiseventb->Fill();
		subeventsinthiseventb->Fill();
		
	} else {
	// MORE THAN ONE MRD SUBEVENT
	// ===========================
		// this event has multiple subevents. Need to split hits into which subevent they belong to.
		// scan over the times and look for gaps where no digits lie, using these to delimit 'subevents'
		std::vector<Float_t> subeventhittimesv;	// a vector of the starting times of a given subevent
		std::vector<double> sorteddigittimes(mrddigittimesthisevent);
		std::sort(sorteddigittimes.begin(), sorteddigittimes.end());
		subeventhittimesv.push_back(sorteddigittimes.at(0));
		for(Int_t i=0;i<sorteddigittimes.size()-1;i++){
			Float_t timetonextdigit = sorteddigittimes.at(i+1)-sorteddigittimes.at(i);
			if(timetonextdigit>maxsubeventduration){
				subeventhittimesv.push_back(sorteddigittimes.at(i+1));
#ifdef MRDSPLITVERBOSE
				cout<<"Setting subevent time threshold at "<<subeventhittimesv.back()<<endl;
//				cout<<"this digit is at "<<sorteddigittimes.at(i)<<endl;
//				cout<<"next digit is at "<<sorteddigittimes.at(i+1)<<endl;
//				try{
//					cout<<"next next digit is at "<<sorteddigittimes.at(i+2)<<endl;
//				} catch (...) { int i=0; }
#endif
			}
		}
#ifdef MRDSPLITVERBOSE
		cout<<subeventhittimesv.size()<<" subevents this event"<<endl;
#endif
		
		// a vector to record the subevent number for each hit, to know if we've allocated it yet.
		std::vector<Int_t> subeventnumthisevent(mrddigittimesthisevent.size(),-1);
		// another for true tracks
		Int_t numtracks = atrigt->GetNtrack();
		std::vector<int> subeventnumthisevent2(numtracks,-1);
		
		// now we need to sort the digits into the subevents they belong to:
		// loop over subevents
		Int_t mrdtrackcounter=0;	// not all the subevents will have a track
		for(Int_t thissubevent=0; thissubevent<subeventhittimesv.size(); thissubevent++){
#ifdef MRDSPLITVERBOSE
			cout<<"Digits in MRD at = "<<subeventhittimesv.at(thissubevent)<<"ns in event "<<eventnum<<endl;
#endif
			// don't need to worry about lower bound as we start from lowest t peak and 
			// exclude already allocated hits
			
			Float_t endtime = (thissubevent<(subeventhittimesv.size()-1)) ? 
				subeventhittimesv.at(thissubevent+1) : (eventendtime+1.);
#ifdef MRDSPLITVERBOSE
			cout<<"endtime for subevent "<<thissubevent<<" is "<<endtime<<endl;
#endif
			// times are not ordered, so scan through all digits for each subevent
			for(Int_t thisdigit=0;thisdigit<mrddigittimesthisevent.size();thisdigit++){
				if(subeventnumthisevent.at(thisdigit)<0 && mrddigittimesthisevent.at(thisdigit)< endtime ){
					// thisdigit is in thissubevent
#ifdef MRDSPLITVERBOSE
					cout<<"adding digit at "<<mrddigittimesthisevent.at(thisdigit)<<" to subevent "<<thissubevent<<endl;
#endif
					digitidsinasubevent.push_back(thisdigit);
					subeventnumthisevent.at(thisdigit)=thissubevent;
					WCSimRootCherenkovDigiHit* thedigihit = (WCSimRootCherenkovDigiHit*)atrigm->GetCherenkovDigiHits()->At(thisdigit);
					Int_t thisdigitstubeid = thedigihit->GetTubeId()-1; //XXX Tube index = TubeId()-1
					tubeidsinasubevent.push_back(thisdigitstubeid);
					Int_t thisdigitsq = thedigihit->GetQ();
					digitqsinasubevent.push_back(thisdigitsq);
					double thisdigitstime = thedigihit->GetT();
					digittimesinasubevent.push_back(thisdigitstime);
					// add all the unique parent ID's for digits contributing to this subevent (truth level info)
					std::vector<int> truephotonindices = thedigihit->GetPhotonIds();
					digitnumtruephots.push_back(truephotonindices.size());
					for(int truephoton=0; truephoton<truephotonindices.size(); truephoton++){
						int thephotonsid = truephotonindices.at(truephoton);
						WCSimRootCherenkovHitTime *thehittimeobject = (WCSimRootCherenkovHitTime*)atrigm->GetCherenkovHitTimes()->At(thephotonsid);
						int thephotonsparentsubevent = thehittimeobject->GetParentID();
						particleidsinasubevent.push_back(thephotonsparentsubevent);
						double thephotonstruetime = thehittimeobject->GetTruetime();
						photontimesinasubevent.push_back(thephotonstruetime);
					}
				}
			}
			
			// construct the subevent from all the digits
			if(digitidsinasubevent.size()>=minimumdigits){	// must have enough for a subevent
				// first scan through all truth tracks to find those within this subevent time window
				std::vector<WCSimRootTrack*> truetrackpointers;
				for(int truetracki=0; truetracki<numtracks; truetracki++){
					WCSimRootTrack* nextrack = (WCSimRootTrack*)atrigt->GetTracks()->At(truetracki);
					if((subeventnumthisevent2.at(truetracki)<0)&&(nextrack->GetTime()<endtime)
						&&(nextrack->GetFlag()==0)&&(nextrack->GetIpnu()!=11)){
						truetrackpointers.push_back(nextrack);
						subeventnumthisevent2.at(truetracki)=thissubevent;
					}
				}
				
#ifdef MRDSPLITVERBOSE
				cout<<"constructing subevent "<<mrdeventcounter<<" with "<<digitidsinasubevent.size()<<" digits"<<endl;
#endif
				cMRDSubEvent* currentsubevent = new((*aSubEvent)[mrdeventcounter]) cMRDSubEvent(mrdeventcounter, currentfilestring, runnum, eventnum, triggernum, digitidsinasubevent, tubeidsinasubevent, digitqsinasubevent, digittimesinasubevent, digitnumtruephots, photontimesinasubevent, particleidsinasubevent, truetrackpointers);
				mrdeventcounter++;
				mrdtrackcounter+=currentsubevent->GetTracks()->size();
#ifdef MRDSPLITVERBOSE
				cout<<"subevent "<<thissubevent<<" found "<<currentsubevent->GetTracks()->size()<<" tracks"<<endl;
#endif
			}
			
			// clear the vectors and loop to the next subevent
			digitidsinasubevent.clear();
			tubeidsinasubevent.clear();
			digitqsinasubevent.clear();
			digittimesinasubevent.clear();
			particleidsinasubevent.clear();
			photontimesinasubevent.clear();
			digitnumtruephots.clear();
			
		}
		
		// quick scan to check for any unallocated hits
		for(Int_t k=0;k<subeventnumthisevent.size();k++){
			if(subeventnumthisevent.at(k)==-1){cout<<"*****unbinned hit!"<<k<<" "<<mrddigittimesthisevent.at(k)<<endl;}
		}
		
		nummrdsubeventsthisevent=mrdeventcounter;
		nummrdtracksthisevent=mrdtrackcounter;
		nummrdsubeventsthiseventb->Fill();
		nummrdtracksthiseventb->Fill();
		subeventsinthiseventb->Fill();
		//mrdtree->Fill();
		
	}	// end multiple subevents case
	
	//if(eventnum==735){ assert(false); }
	//if(nummrdtracksthisevent) std::this_thread::sleep_for (std::chrono::seconds(5));
	
	// WRITE+CLOSE OUTPUT FILES
	// ========================
#ifdef MRDSPLITVERBOSE
	cout<<"writing output files, end of finding MRD SubEvents and Tracks in this event"<<endl;
#endif
	mrdtrackfile->cd();
	mrdtree->SetEntries(nummrdtracksthiseventb->GetEntries());
	mrdtree->Write("",kOverwrite);
	if(cMRDSubEvent::imgcanvas) cMRDSubEvent::imgcanvas->Update();
	gROOT->cd();
}

// #######################################################################
