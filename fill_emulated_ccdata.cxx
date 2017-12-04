/* vim:set noexpandtab tabstop=4 wrap */
// #######################################################################

// called in DoMRDdigitHits() loop over MRD digits.
void WCSimAnalysis::AddCCDataEntry(WCSimRootCherenkovDigiHit* digihit){
	//WCSimRootChernkovDigiHit has methods GetTubeId(), GetT(), GetQ()
	std::vector<unsigned int> placeholdervector{};
	fileout_Trigger = 0;                       // FIXME
	fileout_OutNumber = 0;                     // FIXME
	fileout_Type = std::vector<string>{};      // FIXME
	fileout_Value = placeholdervector;         // FIXME
	fileout_Slot = placeholdervector;          // FIXME
	fileout_Channel = placeholdervector;       // FIXME
	fileout_TimeStamp = digihit->GetT();       // FIXME
	tCCData->Fill();
}

void WCSimAnalysis::FillEmulatedCCData(){
	// called in post mrd hit loop - if needed, add code here.
}

