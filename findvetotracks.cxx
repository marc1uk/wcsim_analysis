/* vim:set noexpandtab tabstop=4 wrap */
// #######################################################################
#include "VetoEventClass.hh"			// a class for defining subevents

#ifndef VETOSPLITVERBOSE
//#define VETOSPLITVERBOSE 1
#endif

// CREATE+OPEN OUTPUT FILE
// =======================
void WCSimAnalysis::OpenFACCtrackOutfile(){
	cout<<"opening facc output file"<<endl;
	vetotrackfile = new TFile("vetotrackfile.root","RECREATE","Veto Events file");
	vetotrackfile->cd();
	vetotree = new TTree("vetotree","Tree for event data");
	gROOT->cd();
	
	numvetoeventsthisevent=0;
	FaccSubEvents = new TClonesArray("cVetoEvent");	// string is class name 
	
	numvetoeventsthiseventb = vetotree->Branch("NumVetoEvents",&numvetoeventsthisevent);
	vetoeventsinthiseventb = vetotree->Branch("VetoEvents",&FaccSubEvents, numvetoeventsthisevent);
	
}

// #######################################################################

// SPLIT HITS BY TIME
// ==================
void WCSimAnalysis::FindVetoTracksInEvent(){

const Int_t minimumdigits =1;
#ifdef VETOSPLITVERBOSE
	cout<<"eventnum is "<<eventnum<<endl;
	cout<<"vetodigittimesthisevent.size()="<<vetodigittimesthisevent.size()<<endl;
#endif
	FaccSubEvents->Clear("C");
/* 
if your class contains pointers, use aTrack.Clear("C"). You MUST then provide a Clear() method in your class that properly performs clearing and memory freeing. (or "implements the reset procedure for pointer objects")
 see https://root.cern.ch/doc/master/classTClonesArray.html#a025645e1e80ea79b43a08536c763cae2
*/
	Int_t vetoeventcounter=0;
	if(vetodigittimesthisevent.size()<minimumdigits){
		// =====================================================================================
		// NO DIGITS IN THIS EVENT
		// ======================================
		//new((*FaccSubEvents)[0]) cVetoEvent();
		numvetoeventsthisevent=vetoeventcounter;
		numvetoeventsthiseventb->Fill();
		vetoeventsinthiseventb->Fill();
		//vetotree->Fill();						// fill the branches so the entries align.
		vetotrackfile->cd();
		vetotree->SetEntries(numvetoeventsthiseventb->GetEntries());
		vetotree->Write("",kOverwrite);
		gROOT->cd();
		return;
		// skip remainder
		// ======================================================================================
	}
	
	// MEASURE EVENT DURATION TO DETERMINE IF THERE IS MORE THAN ONE MRD SUB EVENT
	// ===========================================================================
	std::vector<Int_t> digitidsinFaccSubEvents;
	std::vector<Int_t> tubeidsinFaccSubEvents;
	std::vector<Double_t> digitqsinFaccSubEvents;
	std::vector<Double_t> digittimesinFaccSubEvents;
	std::vector<Int_t> digitnumtruephots;
	std::vector<Int_t> particleidsinFaccSubEvents;
	std::vector<Double_t> photontimesinFaccSubEvents;
	
	// first check: are all hits within a 30ns window (maxsubeventduration) If so, just one subevent. 
	Double_t eventendtime = *std::max_element(vetodigittimesthisevent.begin(),vetodigittimesthisevent.end());
	Double_t eventstarttime = *std::min_element(vetodigittimesthisevent.begin(),vetodigittimesthisevent.end());
	Double_t eventduration = (eventendtime - eventstarttime);
#ifdef VETOSPLITVERBOSE
	cout<<"veto event start: "<<eventstarttime<<", end : "<<eventendtime<<", duration : "<<eventduration<<endl;
#endif
	
	if((eventduration<maxsubeventduration)&&(vetodigittimesthisevent.size()>=minimumdigits)){
	// JUST ONE SUBEVENT
	// =================
#ifdef VETOSPLITVERBOSE
		cout<<"all veto hits this event within one veto event."<<endl;
#endif
		// loop over digits and extract info
		for(Int_t thisdigit=0;thisdigit<vetodigittimesthisevent.size();thisdigit++){
			digitidsinFaccSubEvents.push_back(thisdigit);
			WCSimRootCherenkovDigiHit* thedigihit =
				(WCSimRootCherenkovDigiHit*)atrigv->GetCherenkovDigiHits()->At(thisdigit);
			Int_t thisdigitstubeid = thedigihit->GetTubeId()-1;  //XXX IMPORTANT: Tube index = TubeId()-1
			tubeidsinFaccSubEvents.push_back(thisdigitstubeid);
			Int_t thisdigitsq = thedigihit->GetQ();
			digitqsinFaccSubEvents.push_back(thisdigitsq);
			double thisdigitstime = thedigihit->GetT();
			digittimesinFaccSubEvents.push_back(thisdigitstime);
			// add all the unique parent ID's for digits contributing to this subevent (truth level info)
			std::vector<int> truephotonindices = thedigihit->GetPhotonIds();
			digitnumtruephots.push_back(truephotonindices.size());
			for(int truephoton=0; truephoton<truephotonindices.size(); truephoton++){
				int thephotonsid = truephotonindices.at(truephoton);
				WCSimRootCherenkovHitTime *thehittimeobject = (WCSimRootCherenkovHitTime*)atrigv->GetCherenkovHitTimes()->At(thephotonsid);
				int thephotonsparentsubevent = thehittimeobject->GetParentID();
				particleidsinFaccSubEvents.push_back(thephotonsparentsubevent);
				double thephotonstruetime = thehittimeobject->GetTruetime();
				photontimesinFaccSubEvents.push_back(thephotonstruetime);
			}
			// append the digit - nah, just construct subevent with all digits
			//FaccSubEventsa.AppendDigit(thisdigitstime, thisdigitsq, thisdigitstubeid,
			//										photontimesinFaccSubEvents,particleidsinFaccSubEvents);
		}
		
		// we store the truth tracks in each cVetoEvent, so we need to pass the ones within the subevent
		// time window to the constructor
		Int_t numtracks = atrigt->GetNtrack();
		std::vector<WCSimRootTrack*> truetrackpointers;
		for(int truetracki=0; truetracki<numtracks; truetracki++){
			WCSimRootTrack* nextrack = (WCSimRootTrack*)atrigt->GetTracks()->At(truetracki);
			if((nextrack->GetFlag()==0)&&(nextrack->GetIpnu()!=11)) truetrackpointers.push_back(nextrack);
		}
		// construct the subevent from all the digits
#ifdef VETOSPLITVERBOSE
		cout<<"constructing a single veto event for this event"<<endl;
#endif
		cVetoEvent* currentsubevent = new((*FaccSubEvents)[0]) cVetoEvent(0, currentfilestring, runnum, eventnum, subtriggernum, digitidsinFaccSubEvents, tubeidsinFaccSubEvents, digitqsinFaccSubEvents, digittimesinFaccSubEvents, digitnumtruephots, photontimesinFaccSubEvents, particleidsinFaccSubEvents, truetrackpointers);
		vetoeventcounter++;
		// can also use 'cVetoEvent* = (cVetoEvent*)FaccSubEvents.ConstructedAt(0);' followed by a bunch of
		// 'Set' calls to set all relevant fields. This bypasses the constructor, calling it only when 
		// necessary, saving time. In that case, we do not need to call FaccSubEventsa.Clear();
		
		numvetoeventsthisevent=1;
		numvetoeventsthiseventb->Fill();
		vetoeventsinthiseventb->Fill();
		
	} else {
	// MORE THAN ONE MRD SUBEVENT
	// ===========================
		// this event has multiple subevents. Need to split hits into which subevent they belong to.
		// scan over the times and look for gaps where no digits lie, using these to delimit 'subevents'
		std::vector<Float_t> subeventhittimesv;	// a vector of the starting times of a given subevent
		std::vector<double> sorteddigittimes(vetodigittimesthisevent);
		std::sort(sorteddigittimes.begin(), sorteddigittimes.end());
		subeventhittimesv.push_back(sorteddigittimes.at(0));
		for(Int_t i=0;i<sorteddigittimes.size()-1;i++){
			Float_t timetonextdigit = sorteddigittimes.at(i+1)-sorteddigittimes.at(i);
			if(timetonextdigit>maxsubeventduration){
				subeventhittimesv.push_back(sorteddigittimes.at(i+1));
#ifdef VETOSPLITVERBOSE
				cout<<"Setting veto event time threshold at "<<subeventhittimesv.back()<<endl;
//				cout<<"this digit is at "<<sorteddigittimes.at(i)<<endl;
//				cout<<"next digit is at "<<sorteddigittimes.at(i+1)<<endl;
//				try{
//					cout<<"next next digit is at "<<sorteddigittimes.at(i+2)<<endl;
//				} catch (...) { int i=0; }
#endif
			}
		}
#ifdef VETOSPLITVERBOSE
		cout<<subeventhittimesv.size()<<" veto events this event"<<endl;
#endif
		
		// a vector to record the subevent number for each hit, to know if we've allocated it yet.
		std::vector<Int_t> subeventnumthisevent(vetodigittimesthisevent.size(),-1);
		// another for true tracks
		Int_t numtracks = atrigt->GetNtrack();
		std::vector<int> subeventnumthisevent2(numtracks,-1);
		
		// now we need to sort the digits into the subevents they belong to:
		// loop over subevents
		for(Int_t thissubevent=0; thissubevent<subeventhittimesv.size(); thissubevent++){
#ifdef VETOSPLITVERBOSE
			cout<<"Digits in Veto at = "<<subeventhittimesv.at(thissubevent)<<"ns in event "<<eventnum<<endl;
#endif
			// don't need to worry about lower bound as we start from lowest t peak and 
			// exclude already allocated hits
			
			Float_t endtime = (thissubevent<(subeventhittimesv.size()-1)) ? 
				subeventhittimesv.at(thissubevent+1) : (eventendtime+1.);
#ifdef VETOSPLITVERBOSE
			cout<<"endtime for veto event "<<thissubevent<<" is "<<endtime<<endl;
#endif
			// times are not ordered, so scan through all digits for each subevent
			for(Int_t thisdigit=0;thisdigit<vetodigittimesthisevent.size();thisdigit++){
				if(subeventnumthisevent.at(thisdigit)<0 && vetodigittimesthisevent.at(thisdigit)< endtime ){
					// thisdigit is in thissubevent
#ifdef VETOSPLITVERBOSE
					cout<<"adding digit at "<<vetodigittimesthisevent.at(thisdigit)<<" to veto event "<<thissubevent<<endl;
#endif
					digitidsinFaccSubEvents.push_back(thisdigit);
					subeventnumthisevent.at(thisdigit)=thissubevent;
					WCSimRootCherenkovDigiHit* thedigihit = (WCSimRootCherenkovDigiHit*)atrigv->GetCherenkovDigiHits()->At(thisdigit);
					Int_t thisdigitstubeid = thedigihit->GetTubeId()-1; //XXX Tube index = TubeId()-1
					tubeidsinFaccSubEvents.push_back(thisdigitstubeid);
					Int_t thisdigitsq = thedigihit->GetQ();
					digitqsinFaccSubEvents.push_back(thisdigitsq);
					double thisdigitstime = thedigihit->GetT();
					digittimesinFaccSubEvents.push_back(thisdigitstime);
					// add all the unique parent ID's for digits contributing to this subevent (truth level info)
					std::vector<int> truephotonindices = thedigihit->GetPhotonIds();
					digitnumtruephots.push_back(truephotonindices.size());
					for(int truephoton=0; truephoton<truephotonindices.size(); truephoton++){
						int thephotonsid = truephotonindices.at(truephoton);
						WCSimRootCherenkovHitTime *thehittimeobject = (WCSimRootCherenkovHitTime*)atrigv->GetCherenkovHitTimes()->At(thephotonsid);
						int thephotonsparentsubevent = thehittimeobject->GetParentID();
						particleidsinFaccSubEvents.push_back(thephotonsparentsubevent);
						double thephotonstruetime = thehittimeobject->GetTruetime();
						photontimesinFaccSubEvents.push_back(thephotonstruetime);
					}
				}
			}
			
			// construct the subevent from all the digits
			if(digitidsinFaccSubEvents.size()>=minimumdigits){	// must have enough for a subevent
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
				
#ifdef VETOSPLITVERBOSE
				cout<<"constructing veto event "<<vetoeventcounter<<" with "<<digitidsinFaccSubEvents.size()<<" digits"<<endl;
#endif
				cVetoEvent* currentsubevent = new((*FaccSubEvents)[vetoeventcounter]) cVetoEvent(vetoeventcounter, currentfilestring, runnum, eventnum, subtriggernum, digitidsinFaccSubEvents, tubeidsinFaccSubEvents, digitqsinFaccSubEvents, digittimesinFaccSubEvents, digitnumtruephots, photontimesinFaccSubEvents, particleidsinFaccSubEvents, truetrackpointers);
				vetoeventcounter++;
			}
			
			// clear the vectors and loop to the next subevent
			digitidsinFaccSubEvents.clear();
			tubeidsinFaccSubEvents.clear();
			digitqsinFaccSubEvents.clear();
			digittimesinFaccSubEvents.clear();
			particleidsinFaccSubEvents.clear();
			photontimesinFaccSubEvents.clear();
			digitnumtruephots.clear();
			
		}
		
		// quick scan to check for any unallocated hits
		for(Int_t k=0;k<subeventnumthisevent.size();k++){
			if(subeventnumthisevent.at(k)==-1){cout<<"*****unbinned hit!"<<k<<" "<<vetodigittimesthisevent.at(k)<<endl;}
		}
		
		numvetoeventsthisevent=vetoeventcounter;
		numvetoeventsthiseventb->Fill();
		vetoeventsinthiseventb->Fill();
		//vetotree->Fill();
		
	}	// end multiple subevents case
	
	//if(eventnum==735){ assert(false); }
	//if(numvetoeventsthisevent) std::this_thread::sleep_for (std::chrono::seconds(5));
	
	// WRITE+CLOSE OUTPUT FILES
	// ========================
#ifdef VETOSPLITVERBOSE
	cout<<"writing output files, there were "<<vetoeventcounter<<" veto Events in this event"<<endl;
#endif
	vetotrackfile->cd();
	vetotree->SetEntries(numvetoeventsthiseventb->GetEntries());
	vetotree->Write("",kOverwrite);
	gROOT->cd();
}

// #######################################################################
