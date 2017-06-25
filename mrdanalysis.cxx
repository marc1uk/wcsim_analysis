/* vim:set noexpandtab tabstop=4 wrap */
//############################################################################################

// MRD PRE-EVENT-LOOP ACTIONS
// ===============================
void WCSimAnalysis::DoMRDpreLoop(){
	DefineMRDhistos();
	//OpenMRDtrackOutfile(wcsimfilenum);	// open file for writing mrd tracks
	//TClonesArray* aTrack; clear this?
	//TClonesArray &aTracka = *aTrack; reset this?
}

//############################################################################################

// MRD EVENT-WIDE ACTIONS
// ===============================
void WCSimAnalysis::DoMRDeventwide(Int_t &numtruehits, Int_t &numdigits){
	numtruehits = atrigm->GetCherenkovHits()->GetEntries();
	numdigits = atrigm->GetCherenkovDigiHits()->GetEntries();
	FillMRDeventWideHists(numtruehits, numdigits);
}

//############################################################################################

// MRD PRE-HIT-LOOP ACTIONS
// ===============================
void WCSimAnalysis::DoMRDpreHitLoop(){
	mrddigittimesthisevent.clear();
}

//############################################################################################

// MRD TRUE HIT ACTIONS
// ===============================
void WCSimAnalysis::DoMRDtrueHits(){
	Int_t numtruehits = atrigm->GetCherenkovHits()->GetEntries();
	for(Int_t i=0; i<numtruehits; i++){
		// retrieve the hit information
		WCSimRootCherenkovHit* hit = (WCSimRootCherenkovHit*)atrigm->GetCherenkovHits()->At(i);
		//WCSimRootCherenkovHit has methods GetTubeId(), GetTotalPe(int). only really need int=0

		// HitTimes and Hits are not 1:1 lists; use hit->TotalPe(0) to map from 1 to the other.
		WCSimRootCherenkovHitTime* hittime = (WCSimRootCherenkovHitTime*)atrigm->GetCherenkovHitTimes()->At(hit->GetTotalPe(0));
		// WCSimRootCherenkovHitTime has methods GetTruetime() and GetParentID(); 
		
		// call functions that use this information
		FillMRDtrueHitsHist(hit, hittime);	// relies on the maps generated in MakePMTmap
	}
}

//############################################################################################

// MRD DIGIT ACTIONS
// ===============================
void WCSimAnalysis::DoMRDdigitHits(){
	Int_t numdigits = atrigm->GetCherenkovDigiHits()->GetEntries();
	for(Int_t i=0; i<numdigits; i++){
		// retrieve the digit information
		// ============================
		WCSimRootCherenkovDigiHit* digihit = (WCSimRootCherenkovDigiHit*)atrigm->GetCherenkovDigiHits()->At(i);
		//WCSimRootChernkovDigiHit has methods GetTubeId(), GetT(), GetQ()
		mrddigittimesthisevent.push_back(digihit->GetT());
		
		// call functions that use this information
		// ========================================
		FillMRDdigiHitsHist(digihit);
	}
}

//############################################################################################

// MRD POST-HIT-LOOP ACTIONS
// ===============================
void WCSimAnalysis::DoMRDpostHitLoop(){
	FindMRDtracksInEvent();
}

//############################################################################################

// MRD POST-EVENT-LOOP ACTIONS
// ===============================
void WCSimAnalysis::DoMRDpostLoop(){
	DrawMRDhistos();
}
