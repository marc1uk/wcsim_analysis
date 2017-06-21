/* vim:set noexpandtab tabstop=4 wrap */
//############################################################################################

// LOOP OVER EVENTS TO DO ANALYSIS
// ===============================
void WCSimAnalysis::DoAnalysis(){

	InitEnvironment();			// things like loading the libraries and class header locations
	LoadInputFiles();			// open input tchain
	MakePMTmap(); 				// map pmt positions from geometry file (only uses geotree)
	GetTreeData(); 				// get branches from tree, get triggers from the first entry

	// Declare loop locals
	// ===================
	eventnum=0;
	treeNumber=0;
	
	// Perform Pre-Loop actions
	DoTankPreLoop();
	DoMRDpreLoop();
	DoVetoPreLoop();
	
	// Loop over events
	// ================
	cout<<"Looping over entries"<<endl;
	int breakearlyat=10;
	int maxdigits=0;
	do {
		// load next entry, including new trees and setting branch addresses when necessary
		int entryvalid = LoadTchainEntry(eventnum);
		if(entryvalid==0 || eventnum==breakearlyat){ break; }
		
		// TODO: should include a loop over subtriggers here
		// TODO: TO BE ABLE TO DO THIS MRDTRACKCLASS AND VETOTRACKCLASS NEED TO SUPPORT
		// MULTIPLE TRIGGERS PER TREE ENTRY
		Int_t subtrigger=0;
		// for(Int_t subtrigger=0; subtrigger< (b->GetNumberOfEvents()); subtrigger++){
		atrigt = b->GetTrigger(subtrigger);
		atrigm = m->GetTrigger(subtrigger);
		atrigv = v->GetTrigger(subtrigger);
		
		header = atrigt->GetHeader();
		//eventnum = header->GetEvtNum();					<<< not sure this works
		runnum = header->GetRun();
		subtriggernum = header->GetSubEvtNumber();
		
		// TANK ANALYSIS
		// =============
		int numtanktruehits=0, numtankdigits=0;
		// pre hit loop actions
		DoTankEventwide(numtanktruehits, numtankdigits);
		// loop over true hits and digits internally
		DoTankPreHitLoop();
		DoTankTrueHits();
		DoTankDigitHits();
		// post hit loop actions
		DoTankPostHitLoop();
		//if(numtankdigits>maxdigits){ maxdigits=numtankdigits; cout<<"maxdigits now "<<maxdigits<<endl; }
//		if(numtankdigits==79){ cout<<numtankdigits<<" digits in this event"<<endl; break; }
		
		// MRD ANALYSIS
		// ============
		int nummrdtruehits=0, nummrddigits=0;
		// pre hit loop actions
		DoMRDeventwide(nummrdtruehits, nummrddigits);
		// loop over true hits and digits internally
		DoMRDpreHitLoop();
		DoMRDtrueHits();
		DoMRDdigitHits();
		// post hit loop actions
		DoMRDpostHitLoop();
		
		// VETO ANALYSIS
		// =============
		int numvetotruehits=0, numvetodigits=0;
		// pre hit loop actions
		DoVetoEventwide(numvetotruehits, numvetodigits);
		// loop over true hits and digits internally
		DoVetoPreHitLoop();
		DoVetoTrueHits();
		DoVetoDigitHits();
		// post hit loop actions
		DoVetoPostHitLoop();
		
		// LOOP TO NEXT EVENT
		// ==================
		//std::this_thread::sleep_for (std::chrono::seconds(5));	// a little wait so we can look at histos
		eventnum++;
	} while (1);
	cout<<"Reached end of TChain"<<endl;
	
	DoTankPostLoop();
	DoMRDpostLoop();
	DoVetoPostLoop();
	
//	DrawGlobalHistos(); 	// doesn't fall into any other category.... 
}

