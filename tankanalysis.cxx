/* vim:set noexpandtab tabstop=4 wrap */
//############################################################################################

// TANK PRE-EVENT-LOOP ACTIONS
// ===============================
void WCSimAnalysis::DoTankPreEventLoop(){
	DefineTankHistos();
}

//############################################################################################

// TANK PRE-TRIGGER LOOP ACTIONS
// ===============================
void WCSimAnalysis::DoTankPreTriggerLoop(){
}

//############################################################################################

// TANK TRIGGER ACTIONS
// ====================
void WCSimAnalysis::DoTankTrigger(Int_t &numtruehits, Int_t &numdigits){
	numtruehits += atrigt->GetCherenkovHits()->GetEntries();
	numdigits += atrigt->GetCherenkovDigiHits()->GetEntries();
	
	if(add_emulated_pmtdata){
		// if starting a fresh readout, set the top level (non-minibuffer level) readout details
		if(minibuffer_id==0) ConstructEmulatedPmtDataReadout();
		AddMinibufferStartTime();
	}
	//cout<<endl<<"ok"<<endl;
}

//############################################################################################

// TANK PRE-HIT-LOOP ACTIONS
// ===============================
void WCSimAnalysis::DoTankPreHitLoop(){
	// these histograms become event-scope.
	if((wallhist->GetEntries()+topcaphist->GetEntries()+bottomcaphist->GetEntries())!=0){
		topcaphist->Reset();
		bottomcaphist->Reset();
		wallhist->Reset();
	}
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
		//WCSimRootCherenkovHit has methods GetTubeId(), GetTotalPe(int). only really need int=0
		
		// HitTimes and Hits are not 1:1 lists; use hit->TotalPe(0) to map from 1 to the other.
		WCSimRootCherenkovHitTime* hittime = (WCSimRootCherenkovHitTime*)atrigt->GetCherenkovHitTimes()->At(hit->GetTotalPe(0));
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
		if(add_emulated_pmtdata) AddPMTDataEntry(digihit);
	}
}

//############################################################################################

// TANK POST-HIT-LOOP ACTIONS
// ===============================
void WCSimAnalysis::DoTankPostHitLoop(){
	if(drawtankhistos){ 
		//DrawTankHistos();
		if((wallhist->GetEntries()+topcaphist->GetEntries()+bottomcaphist->GetEntries())!=0){
			wallmapcanv->cd();
			// for no weighting, no range setting necessary? maybe SetMinimum(0)? all bins always 1???
			// for weighting by Q, no range setting necessary
			// for weighting by timing, need to account for 'trigger' offset of 950ns. 
//			wallhist->SetMaximum(970.);			// for timing - maybe only min necessary
			wallhist->SetMinimum(950.);
			wallhist->SetBit(TH1::kNoStats);
			wallhist->Draw("colz");
			topcapmapcanv->cd();
			topcaphist->SetMinimum(950.);
			topcaphist->SetBit(TH1::kNoStats);
			topcaphist->Draw("colz");
			bottomcapmapcanv->cd();
			bottomcaphist->SetMinimum(950.);
			bottomcaphist->SetBit(TH1::kNoStats);
			bottomcaphist->Draw("colz");
		}
	}
}

//############################################################################################

// TANK POST-TRIGGER LOOP ACTIONS
// ===============================
void WCSimAnalysis::DoTankPostTriggerLoop(Int_t &numtruehits, Int_t &numdigits){
	FillTankEventWideHists(numtruehits, numdigits);
}

//############################################################################################

// TANK POST-EVENT-LOOP ACTIONS
// ===============================
void WCSimAnalysis::DoTankPostEventLoop(){
	DrawTankHistos();
}
