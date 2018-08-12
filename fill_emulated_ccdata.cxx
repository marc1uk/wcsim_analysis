/* vim:set noexpandtab tabstop=4 wrap */
// #######################################################################

// called in DoMRDdigitHits() loop over MRD digits.
void WCSimAnalysis::AddCCDataEntry(WCSimRootCherenkovDigiHit* digihit){
	//WCSimRootChernkovDigiHit has methods GetTubeId(), GetT(), GetQ()
	
#if FILE_VERSION<4
	int digits_time_ns = digihit->GetT() + pre_trigger_window_ns - 950;
#else
	int digits_time_ns = digihit->GetT() + pre_trigger_window_ns;
#endif
	// MRD TDC common stop timout: 85 x 50ns 
	if(digits_time_ns>MRD_TIMEOUT_NS) return;
	
	// remember tubeids start from 1
	int channelnum = (digihit->GetTubeId()-1)%channels_per_tdc_card;
	int cardid = ((digihit->GetTubeId()-1)-channelnum)/channels_per_tdc_card;
	
	// convert to 'value' = number of TDC ticks
	int digits_time_ticks = digits_time_ns / MRD_NS_PER_SAMPLE;
	fileout_Value.push_back(digits_time_ticks);
	fileout_Slot.push_back(cardid);
	fileout_Channel.push_back(channelnum);
}

void WCSimAnalysis::FillEmulatedCCData(){
	// called in post mrd hit loop
	// the following are all per-readout (per trigger).
	if(fileout_Value.size()==0) return; // no hits, no entries.
	
	// Timestamp is applied by the PC post-readout so is actually delayed from the trigger!
	unsigned long long timestamp_ms = (static_cast<unsigned long long>(header->GetDate()) 
		+ placeholder_date_ns + MRD_TIMESTAMP_DELAY) / 1000000.;
	
	fileout_Trigger = eventnum+triggernum; // TDC readout number
	fileout_TimeStamp = timestamp_ms; // UTC MS since unix epoch
	fileout_OutNumber = fileout_Value.size();  //number of hits in this event
	fileout_Type.assign(fileout_OutNumber,"TDC"); // all cards are TDC cards for now.
	tCCData->Fill();
}

