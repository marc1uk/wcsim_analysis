//############################################################################################

// LOOP OVER EVENTS TO DO ANALYSIS
// ===============================
void WCSimAnalysis::DoAnalysis(){

	InitEnvironment();	// things like loading the libraries and class header locations
	LoadInputFiles();	// open input tchain
	MakePMTmap(); 		// map pmt positions from geometry file (only uses geotree)
	GetTreeData(); 		// get hit & digit data from file
	// declare histograms to be filled in loop
	DefineTankHistos();
	DefineMRDhistos();
	DefineVetoHistos();
	//OpenMRDtrackOutfile();	// open file for writing mrd tracks

	// Declare loop locals
	// ===================
	Int_t eventnum=0;
	treeNumber=0;
	maxtrackduration=30.;
	
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
		
		// > TANK ANALYSIS
		// > =============
//		int numtanktruehits=0, numtankdigits=0;
//		DoTankEventwide(numtanktruehits, numtankdigits);
//		// loop over true hits and digits internally
//		DoTankTrueHits(numtanktruehits);
//		DoTankDigitHits(numtankdigits);
		
		// > MRD ANALYSIS
		// > ============
		int nummrdtruehits=0, nummrddigits=0;
		DoMRDeventwide(nummrdtruehits, nummrddigits);
		// loop over true hits and digits internally
		DoMRDtrueHits(nummrdtruehits);
		DoMRDdigitHits(nummrddigits);
		
		// > VETO ANALYSIS
		// > =============
//		int numvetotruehits=0, numvetodigits=0;
//		DoVetoEventwide(numvetotruehits, numvetodigits);
//		// loop over true hits and digits internally
//		DoVetoTrueHits(numvetotruehits);
//		DoVetoDigitHits(numvetodigits);
		
		// > LOOP TO NEXT EVENT
		// > ==================
		eventnum++;
	} while (1);
	
//	DrawGlobalHistos();
//	DrawTankHistos();
//	DrawMRDhistos();
//	DrawVetoHistos();
	cout<<"Reached end of TChain"<<endl;
}

