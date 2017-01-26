/* vim:set noexpandtab tabstop=4 wrap */
//############################################################################################

// TANK PRE-EVENT-LOOP ACTIONS
// ===============================
void WCSimAnalysis::DoTankPreLoop(){
	DefineTankHistos();
}

//############################################################################################

// TANK EVENT-WIDE ACTIONS
// ===============================
void WCSimAnalysis::DoTankEventwide(Int_t &numtruehits, Int_t &numdigits){
	numtruehits = atrigt->GetCherenkovHits()->GetEntries();
	numdigits = atrigt->GetCherenkovDigiHits()->GetEntries();
	FillTankEventWideHists(numtruehits, numdigits);
}

//############################################################################################

// TANK PRE-HIT-LOOP ACTIONS
// ===============================
void WCSimAnalysis::DoTankPreHitLoop(){
// nothing here yet
}

//############################################################################################

// TANK TRUE HIT ACTIONS
// ===============================
void WCSimAnalysis::DoTankTrueHits(){
	Int_t numtruehits = atrigt->GetCherenkovHits()->GetEntries();
	for(Int_t i=0; i<numtruehits; i++){
		// retrieve the hit information
		// ============================
		WCSimRootCherenkovHit* hit = (WCSimRootCherenkovHit*)atrigt->GetCherenkovHits()->At(i);
		//WCSimRootCherenkovHit has methods GetTubeId(), GetTotalPe(int)
		
		// how is HitTimes related? Is it an event-wide? digit wide?....
		WCSimRootCherenkovHitTime* hittime = (WCSimRootCherenkovHitTime*)atrigt->GetCherenkovHitTimes()->At(i);
		// WCSimRootCherenkovHitTime has methods GetTruetime() and GetParentID(); 

		// call functions that use this information
		// ========================================
		FillTankTrueHitsHist(hit);	// relies on the maps generated in MakePMTmap
	}
}

//############################################################################################

// TANK DIGIT ACTIONS
// ===============================
void WCSimAnalysis::DoTankDigitHits(){
	Int_t numdigits = atrigt->GetCherenkovDigiHits()->GetEntries();
	for(Int_t i=0; i<numdigits; i++){
		// retrieve the digit information
		// ============================
		WCSimRootCherenkovDigiHit* digihit = (WCSimRootCherenkovDigiHit*)atrigt->GetCherenkovDigiHits()->At(i);
		//WCSimRootChernkovDigiHit has methods GetTubeId(), GetT(), GetQ()

		// call functions that use this information
		// ========================================
		FillTankDigiHitsHist(digihit);
	}
}

//############################################################################################

// TANK POST-HIT-LOOP ACTIONS
// ===============================
void WCSimAnalysis::DoTankPostHitLoop(){
// nothing here yet
}

//############################################################################################

// TANK POST-EVENT-LOOP ACTIONS
// ===============================
void WCSimAnalysis::DoTankPostLoop(){
	DrawTankHistos();
}
