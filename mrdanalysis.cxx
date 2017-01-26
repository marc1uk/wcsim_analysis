//############################################################################################

// MRD EVENT-WIDE ACTIONS
// ===============================
void WCSimAnalysis::DoMRDeventwide(Int_t &numtruehits, Int_t &numdigits){
	numtruehits = atrigm->GetCherenkovHits()->GetEntries();
	numdigits = atrigm->GetCherenkovDigiHits()->GetEntries();
	FillMRDeventWideHists(numtruehits, numdigits);
	
	FindMRDtracksInEvent();
}

//############################################################################################

// MRD TRUE HIT ACTIONS
// ===============================
void WCSimAnalysis::DoMRDtrueHits(Int_t numtruehits){
	for(Int_t i=0; i<numtruehits; i++){
		// retrieve the hit information
		WCSimRootCherenkovHit* hit = (WCSimRootCherenkovHit*)atrigm->GetCherenkovHits()->At(i);
		//WCSimRootCherenkovHit has methods GetTubeId(), GetTotalPe(int)

		// how is HitTimes related? Is it an event-wide? digit wide?....
		WCSimRootCherenkovHitTime* hittime = (WCSimRootCherenkovHitTime*)atrigm->GetCherenkovHitTimes()->At(i);
		// WCSimRootCherenkovHitTime has methods GetTruetime() and GetParentID(); 
		
		// call functions that use this information
		FillMRDtrueHitsHist(hit, hittime);	// relies on the maps generated in MakePMTmap
	}
}

//############################################################################################

// MRD DIGIT ACTIONS
// ===============================
void WCSimAnalysis::DoMRDdigitHits(Int_t numdigits){
	for(Int_t i=0; i<numdigits; i++){
		// retrieve the digit information
		// ============================
		WCSimRootCherenkovDigiHit* digihit = (WCSimRootCherenkovDigiHit*)atrigm->GetCherenkovDigiHits()->At(i);
		//WCSimRootChernkovDigiHit has methods GetTubeId(), GetT(), GetQ()

		// call functions that use this information
		// ========================================
		FillMRDdigiHitsHist(digihit);
	}
}
