/* vim:set noexpandtab tabstop=4 wrap */
// #######################################################################

namespace {
	// Used to convert between seconds and nanoseconds
	constexpr unsigned long long BILLION = 1000000000;
	// The cards have clock frequencies of 125 MHz, so they take samples every 8 ns.
	constexpr unsigned long long CLOCK_TICK = 8; // ns
}

void WCSimAnalysis::AddPMTDataEntry(WCSimRootCherenkovDigiHit* digihit){
	//WCSimRootChernkovDigiHit has methods GetTubeId(), GetT(), GetQ()
	
	/*
	   annie readout: when a (ndigits) trigger happens, 2us is read out. 
	   if another trigger happens (during this)<<<??? (how does trigger detection during ndigits trigger work?)
	   another minibuffer is appended to it and readout.
	   
	   WE SHOULD HAVE PRE-TRIGGER WINDOW ??? AND POST-TRIGGERWINDOW ??? (SUM =2US),
	   CURRENT FILES HAVE PRE-TRIGGER WINDOW 400ns AND POST-TRIGGER WINDOW 950ns
	   (or pre-trigger 1us and post-trigger 2us for MRD+VETO)
	   current files have digit times relative to the trigger time + 950:
	   i.e. after removing 950, all negative time digits are in pre-trigger.
	   WE ALSO NEED TO SHIFT UPPERLIMIT IF WE TRUNCATE THE LOWER LIMIT when combining triggers.
	   in current files if the next trigger window starts before the last one ends, the start time
	   of the next trigger is delayed, but the end is not: so the trigger window is shorter.
	   note that when scanning for triggers, when a trigger is found the scan is resumed from the
	   triggertime + post-trigger window. This means only the pre-trigger-window of the subsequent
	   trigger could overlap, so the tructation of the second trigger is somewhat limited - 
	   the second trigger size will be at least post-trigger-window in length.
	   
	   trigger integration time (ndigitsWindow) should represent max light travel time + scattering, 
	   so that light from an event has time to arrive at opposite side of detector.
	   this shuld also be shortened for ANNIE
	*/
	
	int channelnum = digihit->GetTubeId()%channels_per_adc_card;
	int cardid = (digihit->GetTubeId()-channelnum)/channels_per_adc_card;
	
	/* digit time is relative to the trigger time (ndigits threshold crossing).
	   we need to convert this to position of the digit within the minibuffer data array 
	   and construct a waveform centred at that location with suitable integral to represent its charge. */
	
	// first get the digit time relative to start of minibuffer
#if FILE_VERSION<4
	int digits_time_ns = digihit->GetT() + pre_trigger_window_ns - 950;
#else
	int digits_time_ns = digihit->GetT() + pre_trigger_window_ns;
#endif
	// convert to index of this time in the minibuffer waveform array
	int digits_time_index = 
		(digits_time_ns * minibuffer_datapoints_per_channel) / (pre_trigger_window_ns+post_trigger_window_ns);
	
	// construct a suitable pulse waveform using the position and charge of the digit
	pulsevector.assign(minibuffer_datapoints_per_channel,0);
	GenerateMinibufferPulse(digits_time_index, digihit->GetQ(), pulsevector);
	
	// calculate the offset of the minibuffer, for this channel, within the full buffer for the readout
	int channeloffset = (channelnum * full_buffer_size) / num_adc_cards;
	int minibufferoffset = minibuffer_id*minibuffer_datapoints_per_channel;
	
	// get the full data array for this card, and an iterator to the appropriate start point
	auto& thiscards_fullbuffer = temporary_databuffers.at(cardid);
	auto minibuffer_start = thiscards_fullbuffer.begin() + channeloffset + minibufferoffset;
	
	// add the generated pulse to the minibuffer at the appropriate location
	std::transform(minibuffer_start, minibuffer_start+minibuffer_datapoints_per_channel, pulsevector.begin(),
		minibuffer_start, std::plus<uint16_t>() );
	
}

void WCSimAnalysis::GenerateMinibufferPulse(int digit_index, double digit_charge, std::vector<uint16_t> &pulsevector){
	// need to construct a waveform which crosses a Hefty threshold at digit_index
	// and has integral digit_charge. A landau function has approximately the right shape
	if(fLandau==nullptr){
		fLandau = new TF1("fLandau","[2]*TMath::Landau(x,[0],[1],1)",-10,100);
		fLandau->SetRange(-10,100);
	}
	// [0] is ~position of maximum, [1] is width (sigma), [2] is a scaling.
	// the last parameter (fixed to 1 above) is 'normalise by dividing by sigma'.
	// with 1, the integral is fixed to 1 and the maximum is varied
	// (with 0 the maximum is fixed to ~0.18 and the integral varies)
	// by choosing 1 we can always get the desired integral (Q). How should we vary height vs width?
	// that's given by the typical aspect ratio of a PMT pulse: 
	// Looking at data: with X scale in samples (2ns) fitting a landau gives a sigma of ~2
	
	fLandau->SetParameters(digit_charge, digit_index, 2.);
	// landau function is interesting in region -5*sigma -> 50*sigma, or for sigma=2, -10 to 100
	//fLandau->SetRange(-10,100); only set in creation as fixed
	
	for(int i=(digit_index-10); i<(digit_index+100); i++){
		// pulses very close to the front/end of the minibuffer: get tructated.
		if( (i<0) || (i>pulsevector.size()) ) continue;
		pulsevector.at(digit_index)=fLandau->Eval(i);
	}
}

void WCSimAnalysis::AddMinibufferStartTime(){
	
	// we need to record in the readout the start time of this minibuffer relative to the card 
	// initialization time (StartTimeSec+StartTimeNSec, or StartCount)
	// since each card gets initialized at a slightly different time, the minibuffer start counts
	// are slightly different in each card, even though they represent the same trigger.
	for(auto&& acard : emulated_pmtdata_readout){
		// get start time of this card in abs ns
		int readout_starttime_s = acard.StartTimeSec;
		int readout_starttime_ns = acard.StartTimeNSec;
		unsigned long long readout_time = static_cast<unsigned long long>(readout_starttime_s)
			+ static_cast<unsigned long long>(readout_starttime_ns);
		
		// calculate start time of this minibuffer (trigger) in abs ns
		unsigned long long minibuffer_time = static_cast<unsigned long long>(header->GetDate()) 
			+ placeholder_date_ns;
		// calculate the clock ticks since that time and when the card was initialized
		int relative_start_time_clockticks = (minibuffer_time - readout_time) / CLOCK_TICK;
		
		acard.TriggerCounts.at(minibuffer_id) = relative_start_time_clockticks;
	}
}

void WCSimAnalysis::ConstructEmulatedPmtDataReadout(){
	// fill all the non-minibuffer info and reset the minibuffers
	
	// going to use the time of the first trigger in the readout as StartTime, for all cards.
	// this isn't true for real data... Should we do something else? FIXME
	unsigned long long triggertime = static_cast<unsigned long long>(header->GetDate()) + placeholder_date_ns;
	int triggersecs = triggertime % BILLION;
	int triggernsecs = triggertime - triggersecs;
	
	if(emulated_pmtdata_readout.size()<num_adc_cards){
		emulated_pmtdata_readout= std::vector<CardData>(num_adc_cards);
	}
	
	for(int cardi=0; cardi<num_adc_cards; cardi++){
		auto&& acard = emulated_pmtdata_readout.at(cardi);
		acard.Reset();
		//acard.LastSync = default BOGUS_INT until i know what to put here.
		acard.SequenceID = sequence_id;
		acard.StartTimeSec = triggersecs;
		acard.StartTimeNSec = triggernsecs;
		//acard.StartCount = default BOGUS_INT until i know what to put here
		acard.CardID = cardi;
		//acard.Data.assign(full_buffer_size,0); not necessary since we're using the temporary buffers
		temporary_databuffers.at(cardi).assign(full_buffer_size,0);
	}
	
	//  Data to be filled while processing triggers:
	//  --------------------------------------------
	//  TriggerCounts.assign(TriggerNumber,BOGUS_UINT64);   // start counts of each minibuffer
	//  Rates.assign(Channels,BOGUS_UINT32);                // avg pulse rate on each PMT
	//  Data.assign(FullBufferSize,BOGUS_UINT16);           // 160k datapoints
	
}

void WCSimAnalysis::FillEmulatedPMTData(){
	//fileout_TriggerCounts = new ULong64_t[minibuffers_per_fullbuffer];  // must be >= fileout_TriggerNumber
	//fileout_Rates         = new Int_t[channels_per_adc_card];           // must be >= fileout_Channels
	//fileout_Data          = new UShort_t[full_buffer_size];             // must be >= fileout_FullBufferSize
	////*............................................................................*
	//TBranch *bLastSync       = tPMTData->Branch("LastSync", &fileout_LastSync);
	//TBranch *bSequenceID     = tPMTData->Branch("SequenceID", &fileout_SequenceID);
	//TBranch *bStartTimeSec   = tPMTData->Branch("StartTimeSec", &fileout_StartTimeSec);
	//TBranch *bStartTimeNSec  = tPMTData->Branch("StartTimeNSec", &fileout_StartTimeNSec);
	//TBranch *bStartCount     = tPMTData->Branch("StartCount", &fileout_StartCount);
	//TBranch *bTriggerNumber  = tPMTData->Branch("TriggerNumber", &fileout_TriggerNumber);
	//TBranch *bTriggerCounts  = tPMTData->Branch("TriggerCounts", &fileout_TriggerCounts);
	//TBranch *bCardID         = tPMTData->Branch("CardID", &fileout_CardID);
	//TBranch *bChannels       = tPMTData->Branch("Channels", &fileout_Channels);
	//TBranch *bRates          = tPMTData->Branch("Rates", &fileout_Rates);
	//TBranch *bBufferSize     = tPMTData->Branch("BufferSize", &fileout_BufferSize);
	//TBranch *bEventsize      = tPMTData->Branch("Eventsize", &fileout_Eventsize);
	//TBranch *bFullBufferSize = tPMTData->Branch("FullBufferSize", &fileout_FullBufferSize);
	//TBranch *bData           = tPMTData->Branch("Data", &fileout_Data);
	
	// PMTData tree has: one entry per VME card, 16 entries (cards) per readout, 40 minibuffers per readout.
	// A single trigger may span multiple minibuffers, but there's nothing to indicate it - just the 
	// timestamps of the starts of each minibuffer will be separated by 2us. 
	// The Data[] array in each entry concatenates all channels on that card:
	// ([Card 0: {minibuffer 0}{minibuffer 1}...{minibuffer 39}][Card 1: {minibuffer 0}...]...)
	// Entries are ordered according to the card position in the vme crate, so are consistent but not 
	// 0,1,2...  (there are 16 cards, numbered up to 21). Each readout has a unique SequenceID.
	// Data[] arrays are waveforms of 40,000 datapoints per minibuffer. 
	
	// first, shuffle your library. this also copies temporary_databuffers into emulated_pmtdata_readout
	RiffleShuffle();
	
	// loop over all cards and fill the PMTData tree with the constructed data
	for(auto& acard : emulated_pmtdata_readout){
		fileout_LastSync          = acard.LastSync;
		fileout_SequenceID        = acard.SequenceID;
		fileout_StartTimeSec      = acard.StartTimeSec;
		fileout_StartTimeNSec     = acard.StartTimeNSec;
		fileout_StartCount        = acard.StartCount;
		fileout_TriggerNumber     = acard.TriggerNumber;
		fileout_CardID            = acard.CardID;
		fileout_Channels          = acard.Channels;
		fileout_BufferSize        = acard.BufferSize;
		fileout_Eventsize         = acard.Eventsize;
		fileout_FullBufferSize    = acard.FullBufferSize;
		fileout_TriggerCounts     = acard.TriggerCounts.data();
		fileout_Rates             = acard.Rates.data();
		fileout_Data              = acard.Data.data();
		
		tPMTData->Fill();
	}
	
}

void WCSimAnalysis::RiffleShuffle(){
/*
	4 channels are separate and isolated within the full buffer:
	[chan 1][chan 2][chan 3][chan 4]
	within each channel, there are 40 minibuffers concatenated:
	[{ch1:mb1}{ch1:mb2}...{ch1:mb40}][{ch2:mb1}{ch2:mb2}...{ch2:mb40}] --- [{ch4:mb1}{ch4:mb2}...{ch4:mb40}]
	but wait! within each channel section, each subarray is split in half and pairs are interleaved:
	so if a channel subarray had indices:
	[0, 1, 2, 3, 4 ... 40k]
	then the actual sample ordering should be:
	[0, 1, 4, 5, 8, 9 ... 20k, 2, 3, 6, 7, 10, 11 ... 40k]
	
	So let's shuffle!
*/
	
	int channel_buffer_size = (full_buffer_size / channels_per_adc_card);
	for(int cardi=0; cardi<num_adc_cards; cardi++){
		auto& acard = emulated_pmtdata_readout.at(cardi);
		auto& card_buffer = acard.Data;
		auto& temp_card_buffer = temporary_databuffers.at(cardi);
		
		// split the deck into four piles
		for(int channeli=0; channeli<channels_per_adc_card; channeli++){
			auto channel_buffer_start = card_buffer.begin() + (channel_buffer_size * channeli);
			auto temp_channel_buffer_start = temp_card_buffer.begin() + (channel_buffer_size * channeli);
			// split each pile into two, then interleave pairs in the two piles
			for(int samplei=0, samplej=0; samplei<channel_buffer_size; samplei+=4, samplej+=2){
				//0 = 0
				*(channel_buffer_start + samplej) = 
					*(temp_channel_buffer_start + samplei);
				//1 = 1
				*(channel_buffer_start + samplej + 1) = 
					*(temp_channel_buffer_start + samplei + 1);
				//20,000 = 2
				*(channel_buffer_start + (channel_buffer_size/2) + samplej ) 
					= *(temp_channel_buffer_start + samplei + 2);
				//20,001 = 3
				*(channel_buffer_start + (channel_buffer_size/2) + samplej + 1) 
					= *(temp_channel_buffer_start + samplei + 3);
			}
		}
	}
	// finally, ask a player to cut the deck
}
