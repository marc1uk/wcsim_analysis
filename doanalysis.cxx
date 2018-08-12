/* vim:set noexpandtab tabstop=4 wrap */
//############################################################################################

#define FILE_VERSION 2
#ifndef VERBOSE
#define VERBOSE 1
#endif
/*
Version 1:
wcsim_wdirt_07-02-17, 200 PMTs + 200 LAPPDs, including dirt intx.

Version 2:
wcsim_tankonly_03-05-17, 200 PMTs + 200 LAPPDs, tank only, with bug fixes, bad lappd pulse timing resoln

Version 3:
wcsim_tankonly_17-06-17, 120 PMTs of 3 different types (LUX, Watchboy, LBNE, 8inHQE), no LAPPDs.

Version 4:
???????????????????????, 120 PMTs of 3 different types, no LAPPDs?, ANNIE Hefty triggering.
*/

// LOOP OVER EVENTS TO DO ANALYSIS
// ===============================
void WCSimAnalysis::DoAnalysis(){

	if(!create_emulated_output){  // ensure all other elements are turned off
		add_emulated_ccdata=false;
		add_emulated_triggerdata=false;
		add_emulated_pmtdata=false;
	}

	InitEnvironment();			// things like loading the libraries and class header locations
	LoadInputFiles();			// open input tchain
	if(create_emulated_output)	// if raw emulation
	LoadOutputFiles();			// create and open output files - must go after loading first entry
	MakePMTmap(); 				// map pmt positions from geometry file (only uses geotree)
	GetTreeData(); 				// get branches from tree, get triggers from the first entry

	// Declare loop locals
	// ===================
	eventnum=0;
	sequence_id=0;
	minibuffer_id=0;
	firstfilenum=0;
	treeNumber=-1;
	
	// Perform Pre-Loop actions
#ifdef VERBOSE
	cout<<"performing Tank PreEventLoop"<<endl;
#endif
	DoTankPreEventLoop();
#ifdef VERBOSE
	cout<<"performing MRD PreEventLoop"<<endl;
#endif
	DoMRDpreEventLoop();
#ifdef VERBOSE
	cout<<"performing FACC PreEventLoop"<<endl;
#endif
	DoVetoPreEventLoop();
	
	if(create_emulated_output) FillEmulatedRunInformation();
	
	// Loop over events
	// ================
	cout<<"Looping over entries"<<endl;
	int breakearlyat=-1;
	int breakearlysequenceid=-1;
	int maxdigits=0;
	do {
		// load next entry, including new trees and setting branch addresses when necessary
#ifdef VERBOSE
		cout<<"loading entry "<<eventnum<<endl;
#endif
		int entryvalid = LoadTchainEntry(eventnum); // note eventnum may be modified within this function!
#ifdef VERBOSE
		cout<<"analyzing event "<<eventnum<<endl;
#endif
		if( (entryvalid==0) || (eventnum>=breakearlyat&&breakearlyat>0)
			|| (sequence_id>breakearlysequenceid&&breakearlysequenceid>0) ){ break; }
		
#ifdef VERBOSE
		cout<<"Doing PreTriggerLoops"<<endl;
#endif
		DoTankPreTriggerLoop();
		DoMRDpreTriggerLoop();
		DoVetoPreTriggerLoop();
		
		// store event wide counters outside the trigger loop.
		int numtankhits=0, numtankdigits=0;
		int nummrdhits=0, nummrddigits=0;
		int numvetohits=0, numvetodigits=0;
		
		triggernum=0;
		for(triggernum=0; triggernum< (b->GetNumberOfEvents()); triggernum++){
#ifdef VERBOSE
			cout<<"Getting Triggers "<<triggernum<<endl;
#endif
			atrigt = b->GetTrigger(triggernum);
			atrigm = m->GetTrigger(triggernum);
			atrigv = v->GetTrigger(triggernum);
			
#ifdef VERBOSE
			cout<<"Getting Header info"<<endl;
#endif
			header = atrigt->GetHeader();
			//int event_id = header->GetEvtNum();               <<< not sure this works
			runnum = header->GetRun();
			//int trigger_id = header->GetSubEvtNumber();
			
			// TANK ANALYSIS
			// =============
#ifdef VERBOSE
			cout<<"Doing Tank Trigger Analaysis"<<endl;
#endif
			DoTankTrigger(numtankhits, numtankdigits);
			// pre hit loop actions
#ifdef VERBOSE
			cout<<"Doing Tank PreHitLoop Analaysis"<<endl;
#endif
			DoTankPreHitLoop();
			// loop over true hits and digits internally
#ifdef VERBOSE
			cout<<"Doing Tank True Hit Analaysis"<<endl;
#endif
			DoTankTrueHits();
#ifdef VERBOSE
			cout<<"Doing Tank Digit Analaysis"<<endl;
#endif
			DoTankDigitHits();
			// post hit loop actions
#ifdef VERBOSE
			cout<<"Doing Tank PostHitLoop Analaysis"<<endl;
#endif
			DoTankPostHitLoop();
			//if(numtankdigits>maxdigits){ maxdigits=numtankdigits; cout<<"maxdigits now "<<maxdigits<<endl; }
			//if(numtankdigits==79){ cout<<numtankdigits<<" digits in this event"<<endl; break; }
			
			// MRD ANALYSIS
			// ============
#ifdef VERBOSE
			cout<<"Doing MRD Analaysis"<<endl;
#endif
			DoMRDtrigger(nummrdhits, nummrddigits);
			// pre hit loop actions
			DoMRDpreHitLoop();
			// loop over true hits and digits internally
			DoMRDtrueHits();
			DoMRDdigitHits();
			// post hit loop actions
			DoMRDpostHitLoop();
			
			// VETO ANALYSIS
			// =============
#ifdef VERBOSE
			cout<<"Doing FACC Analaysis"<<endl;
#endif
			DoVetoTrigger(numvetohits, numvetodigits);
			// pre hit loop actions
			DoVetoPreHitLoop();
			// loop over true hits and digits internally
			DoVetoTrueHits();
			DoVetoDigitHits();
			// post hit loop actions
			DoVetoPostHitLoop();
			
			if(create_emulated_output){
				// advance the counter of triggers (minibuffers)
				minibuffer_id++;
				if(minibuffer_id==minibuffers_per_fullbuffer){
#ifdef VERBOSE
					cout<<"#########################"<<endl;
					cout<<"Filling Emulated Data"<<endl;
#endif
					if(add_emulated_pmtdata) FillEmulatedPMTData();
					if(add_emulated_triggerdata) FillEmulatedTrigData();
					sequence_id++;
					minibuffer_id=0;
				}
			}
			// LOOP TO NEXT TRIGGER
			// ====================
#ifdef VERBOSE
			cout<<"Looping to next trigger"<<endl;
#endif
			b->ReInitialize();
			m->ReInitialize();
			v->ReInitialize();
		}
		
#ifdef VERBOSE
			cout<<"Doing PostTriggerLoops"<<endl;
#endif
		DoTankPostTriggerLoop(numtankhits, numtankdigits);
		DoMRDpostTriggerLoop(nummrdhits, nummrddigits);
		DoVetoPostTriggerLoop(numvetohits, numvetodigits);
		
		//std::this_thread::sleep_for (std::chrono::seconds(5));	// a little wait so we can look at histos
		eventnum++;
		
		cout<<"minibuffer_id="<<minibuffer_id<<", sequence_id="<<sequence_id<<endl;
		
		// LOOP TO NEXT EVENT
		// ==================
	} while (1);
	cout<<"Reached end of TChain: analysed "<<eventnum<<" events"<<endl;
	
#ifdef VERBOSE
	cout<<"doing PostEventLoops"<<endl;
#endif
	DoTankPostEventLoop();
	DoMRDpostEventLoop();
	DoVetoPostEventLoop();
	
	if(rawfileout) rawfileout->Write();
	
//	DrawGlobalHistos(); 	// doesn't fall into any other category.... 
}

