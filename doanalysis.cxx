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
	Int_t eventnum=0;
	treeNumber=0;
	maxtrackduration=30.;
	
	// Perform Pre-Loop actions
	DoTankPreLoop();
	DoMRDpreLoop();
	DoVetoPreLoop();
	
	// Loop over events
	// ================
	cout<<"Looping over entries"<<endl;
	//for(eventnum=0; eventnum<bp->GetEntries();eventnum++){
	do {
		// load next entry, including new trees and setting branch addresses when necessary
		int entryvalid = LoadTchainEntry(eventnum);
		if(entryvalid==0){ break; }
		
		// TODO: should include a loop over subtriggers here
		Int_t subtrigger=0;
		// for(Int_t subtrigger=0; subtrigger< (b->GetNumberOfEvents()); subtrigger++){
		atrigt = b->GetTrigger(subtrigger);
		atrigm = m->GetTrigger(subtrigger);
		atrigv = v->GetTrigger(subtrigger);
		
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
		eventnum++;
	} while (1);
	cout<<"Reached end of TChain"<<endl;
	
	DoTankPostLoop();
	DoMRDpostLoop();
	DoVetoPostLoop();
	
	DrawGlobalHistos(); 	// doesn't fall into any other category.... 
}

