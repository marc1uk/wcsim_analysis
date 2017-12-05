/* vim:set noexpandtab tabstop=4 wrap */
//############################################################################################

// VETO PRE-EVENT-LOOP ACTIONS
// ===============================
void WCSimAnalysis::DoVetoPreEventLoop(){
	DefineVetoHistos();
	//OpenFACCtrackOutfile(wcsimfilenum);
}

//############################################################################################

// VETO PRE-TRIGGER LOOP ACTIONS
// ===============================
void WCSimAnalysis::DoVetoPreTriggerLoop(){
}

//############################################################################################

// VETO TRIGGER ACTIONS
// ===============================
void WCSimAnalysis::DoVetoTrigger(Int_t &numtruehits, Int_t &numdigits){
	numtruehits += atrigv->GetCherenkovHits()->GetEntries();
	numdigits += atrigv->GetCherenkovDigiHits()->GetEntries();
}

//############################################################################################

// VETO PRE-HIT-LOOP ACTIONS
// ===============================
void WCSimAnalysis::DoVetoPreHitLoop(){
	vetodigittimesthisevent.clear();
}

//############################################################################################

// VETO TRUE HIT ACTIONS
// ===============================
void WCSimAnalysis::DoVetoTrueHits(){
	Int_t numtruehits = atrigv->GetCherenkovHits()->GetEntries();
	for(Int_t i=0; i<numtruehits; i++){
		// retrieve the hit information
		WCSimRootCherenkovHit* hit = (WCSimRootCherenkovHit*)atrigv->GetCherenkovHits()->At(i);
		//WCSimRootCherenkovHit has methods GetTubeId(), GetTotalPe(int). only really need int=0

		// HitTimes and Hits are not 1:1 lists; use hit->TotalPe(0) to map from 1 to the other.
		WCSimRootCherenkovHitTime* hittime = (WCSimRootCherenkovHitTime*)atrigv->GetCherenkovHitTimes()->At(hit->GetTotalPe(0));
		// WCSimRootCherenkovHitTime has methods GetTruetime() and GetParentID(); 
		
		// call functions that use this information
		FillVetoTrueHitsHist(hit, hittime);	// relies on the maps generated in MakePMTmap
	}
}

//############################################################################################

// VETO DIGIT ACTIONS
// ===============================
void WCSimAnalysis::DoVetoDigitHits(){
	Int_t numdigits = atrigv->GetCherenkovDigiHits()->GetEntries();
	for(Int_t i=0; i<numdigits; i++){
		// retrieve the digit information
		// ============================
		WCSimRootCherenkovDigiHit* digihit = (WCSimRootCherenkovDigiHit*)atrigv->GetCherenkovDigiHits()->At(i);
		//WCSimRootChernkovDigiHit has methods GetTubeId(), GetT(), GetQ()
		vetodigittimesthisevent.push_back(digihit->GetT());

		// call functions that use this information
		// ========================================
		FillVetoDigiHitsHist(digihit);
		
		if(add_emulated_ccdata) AddCCDataEntry(digihit);
	}
}

//############################################################################################

// VETO POST-HIT-LOOP ACTIONS
// ===============================
void WCSimAnalysis::DoVetoPostHitLoop(){
	//FindVetoTracksInEvent();
	//FillEmulatedCCData(); Fill() is called per digit.
}

//############################################################################################

// VETO POST-TRIGGER LOOP ACTIONS
// ===============================
void WCSimAnalysis::DoVetoPostTriggerLoop(Int_t &numtruehits, Int_t &numdigits){
	FillVetoEventWideHists(numtruehits, numdigits);
}

//############################################################################################

// VETO POST-EVENT-LOOP ACTIONS
// ===============================
void WCSimAnalysis::DoVetoPostEventLoop(){
	DrawVetoHistos();
}
