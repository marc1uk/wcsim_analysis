/* vim:set noexpandtab tabstop=4 wrap */
// #######################################################################
#ifndef EMULATED_OUT_VERBOSE
#define EMULATED_OUT_VERBOSE 1
#endif

void WCSimAnalysis::AddPMTDataEntry(WCSimRootCherenkovDigiHit* digihit){
#ifdef EMULATED_OUT_VERBOSE
	//cout<<"adding emulated PMT data digit"<<endl;
#endif
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
	
	// remember tubeids start from 1
	int channelnum = (digihit->GetTubeId()-1)%channels_per_adc_card;
	int cardid = ((digihit->GetTubeId()-1)-channelnum)/channels_per_adc_card;
	
	/* digit time is relative to the trigger time (ndigits threshold crossing).
	   we need to convert this to position of the digit within the minibuffer data array 
	   and construct a waveform centred at that location with suitable integral to represent its charge. */
	
	// first get the digit time relative to start of minibuffer
	int digits_time_ns = digihit->GetT() + pre_trigger_window_ns - triggeroffset;
	// convert to index of this time in the minibuffer waveform array
	int digits_time_index = digits_time_ns / ADC_NS_PER_SAMPLE;
	// int digits_time_index =
	//	(digits_time_ns * minibuffer_datapoints_per_channel) / (pre_trigger_window_ns+post_trigger_window_ns);
	//cout<<"digit time is "<<digits_time_ns<<"ns, total trigger window time is "
	//	<<(pre_trigger_window_ns+post_trigger_window_ns)<<", so digit lies at index "<<digits_time_index
	//	<<" in the trigger window."<<endl;
	
	// need to scale area so that amplitude is in ADC counts and time is in samples
	//convert digihit->GetQ() pe's to Coulombs by * e.
	//C corresponds to y axis I[A], x axis s. Y AXIS (I) IS CALCULATED AUTOMATICALLY.
	//then change x axis from s to (ns/8) => area = area * 1e9 * (1/8)
	//then change y axis from I to ADC => ADC = (V/ADC_TO_VOLT[V/ct]) = (I*R/ADC_TO_VOLT) = I * (R / ADC_TO_VOLT)
	//so total scaling area = digihit->GetQ() * e * 1e9 * (1/8) * 50 * (1/ADC_TO_VOLT)
	// ADC_TO_VOLT ~0.0006
	
	// digit Q is ~0-30... units? few pe * gain 10^7 = few 10^7 pe's... 
	// charge on electron = 10^-19C = few 10^-12 C... where 30 from?
	double adjusted_digit_q = digihit->GetQ() * (1./ADC_NS_PER_SAMPLE) /* * pow(10.,9.) */
		* (ADC_INPUT_RESISTANCE/ADC_TO_VOLT) * PULSE_HEIGHT_FUDGE_FACTOR;
	
	// construct a suitable pulse waveform using the position and charge of the digit
	pulsevector.assign(minibuffer_datapoints_per_channel,0);
	GenerateMinibufferPulse(digits_time_index, adjusted_digit_q, pulsevector);
	
	// calculate the offset of the minibuffer, for this channel, within the full buffer for the readout
	int channeloffset = channelnum * (full_buffer_size / channels_per_adc_card);
	int minibufferoffset = minibuffer_id*minibuffer_datapoints_per_channel;
	//cout<<"this event+channel's minbuffer starts at index "<<channeloffset + minibufferoffset
	//	<<" within the full buffer"<<endl;
	
	// get the full data array for this card, and an iterator to the appropriate start point
	auto& thiscards_fullbuffer = temporary_databuffers.at(cardid);
	auto minibuffer_start = thiscards_fullbuffer.begin() + channeloffset + minibufferoffset;
	
	// add the generated pulse to the minibuffer at the appropriate location
	std::transform(minibuffer_start, minibuffer_start+minibuffer_datapoints_per_channel, pulsevector.begin(),
		minibuffer_start, std::plus<uint16_t>() );
	
	//int totalindex=channeloffset + minibufferoffset + digits_time_index;
	//cout<<"NS: "<<digits_time_index<<", I: "<<totalindex<<", E: "<<eventnum<<", T: "<<triggernum<<", S: "<<sequence_id<<", M: "<<minibuffer_id<<endl;
	
	// draw for debug
//	//cout<<"making graph of this channel's buffer of size "<<(full_buffer_size/channels_per_adc_card)<<endl;
//	static TGraph* buffergraph = new TGraph(full_buffer_size/channels_per_adc_card);
//	static TCanvas* buffercanvas = new TCanvas();
//	int i=0;
//	cout<<"adding hit at "<<digits_time_ns<<endl;
//	for(auto it=thiscards_fullbuffer.begin()+channeloffset; it!=thiscards_fullbuffer.end(); it++){
//		auto theval = *it;
//		buffergraph->SetPoint(i,i,theval); i++;
//		if(theval!=0) cout<<"("<<i<<","<<theval<<")"<<endl;
//	}
//	buffercanvas->cd();
//	buffergraph->Draw("alp");
//	buffercanvas->Update();
//	//std::this_thread::sleep_for (std::chrono::seconds(5));
//	gPad->WaitPrimitive();
	
}

void WCSimAnalysis::GenerateMinibufferPulse(int digit_index, double adjusted_digit_q, std::vector<uint16_t> &pulsevector){
//#ifdef EMULATED_OUT_VERBOSE
//	cout<<"generating emulated digit pulse"<<endl;
//#endif
	// need to construct a waveform which crosses a Hefty threshold at digit_index
	// and has integral digit_charge. A landau function has approximately the right shape
	if(fLandau==nullptr){
		fLandau = new TF1("fLandau","[2]*TMath::Landau(x,[0],[1],1)",-10,100);
	}
	// [0] is ~position of maximum, [1] is width (sigma), [2] is a scaling.
	// the last parameter (fixed to 1 above) is 'normalise by dividing by sigma'.
	// with 1, the integral is fixed to 1 and the maximum is varied
	// (with 0 the maximum is fixed to ~0.18 and the integral varies)
	// by choosing 1 we can always get the desired integral (Q). How should we vary height vs width?
	// that's given by the typical aspect ratio of a PMT pulse: 
	// Looking at data: with X scale in samples (8ns) fitting a landau gives a sigma of ~2
	// typical digit Qs are 0-30. (PEs?)
	
	// TODO improve this by trading off the width vs height based on the time between the first and last
	// photons within the digit XXX XXX XXX XXX
	
	fLandau->SetParameters(digit_index, 2., adjusted_digit_q);
	// landau function is interesting in region -5*sigma -> 50*sigma, or for sigma=2, -10 to 100
	//fLandau->SetRange(-10,100); only set in creation as fixed
	
	for(int i=(digit_index-10); i<(digit_index+100); i++){
		// pulses very close to the front/end of the minibuffer: get tructated.
		if( (i<0) || (i>(pulsevector.size()-1)) ) continue;
		pulsevector.at(i)=fLandau->Eval(i);
	}
	
//	// draw for debug
//	static TGraph* pulsegraph = new TGraph(pulsevector.size());
//	static TCanvas* pulsecanvas = new TCanvas();
//	int i=0;
//	for(i=0; i<pulsevector.size(); i++) pulsegraph->SetPoint(i,(float)i,(float)pulsevector.at(i));
//	i=0;
//	double px, py;
//	for(i=0; i<pulsegraph->GetN();i++){pulsegraph->GetPoint(i,px,py); cout<<" "<<i<<"("<<px<<","<<py<<")";} cout<<endl;
//	pulsecanvas->cd();
//	pulsegraph->Draw("alp");
//	pulsecanvas->Update();
//	//std::this_thread::sleep_for (std::chrono::seconds(5));
//	//gPad->WaitPrimitive();
}

void WCSimAnalysis::AddMinibufferStartTime(bool droppingremainingsubtriggers){
#ifdef EMULATED_OUT_VERBOSE
	cout<<"Adding a new minibuffer start time to emulated PMT data"<<endl;
#endif
	
	// we need to make the timing variables correctly line up to give the correct TSinceBeam.
	// It's probably possible to use TriggerCounts, which defines the start time of the minibuffer
	// in clock ticks, but this is not relative to the beam arrival. In fact the relationship between
	// TSinceBeam (which we have) and the timing variables in PMTData and TrigData is non-trivial.
	// For now, though, we can just skip proper simulation of remaining timing variables.
	// This is because ToolAnalysis doesn't actually read the raw timing variables, but instead
	// reads TSinceBeam from a pre-processed 'Hefty Timing file', which we can also generate.
	
	if(triggernum==0){
		// First, to replicate the beam structure and timing, we need to add some running time offsets.
		// We expect ~0.03 events per spill, so throw from a Poisson distribution with mean of 1/0.03
		// to see how many beam spills there were since the last event.
		Int_t numspillssincelast = R.Poisson(1./0.03);
		// spills occur at an (average) rate of 7.5Hz. Add time to the running timer since last spill.
		runningeventtime += (numspillssincelast * (1./7.5) * SEC_TO_NS);
		
		// now calculate time from start of the spill to the first trigger (start of this event)
		// each spill lasts 1.6us, but is made up of 84 bunches, each separated by 19ns
		int starttimeofbunch = floor(R.Uniform(0,84)) * 19;
		// each bunch is a gaussian with width +-2.5ns
		int timewithinbunch = R.Gaus(0.,2.5); // add fractional ns for to subsequent digit times? TODO
		// this gives the time to the start of the event within the beam spill
		currenteventtime = starttimeofbunch + timewithinbunch;
	}
	
	// get start time of this minibuffer (trigger) since event start in ns
	// FIXME do we need to subtract the ~6670ns offset from upstream?
	int minibuffer_start_time_ns = header->GetDate();
	
	// and for now put it into the TriggerCounts variable
	timefileout_TSinceBeam[minibuffer_id] = minibuffer_start_time_ns;
	timefileout_More[minibuffer_id]=droppingremainingsubtriggers;
	timefileout_Label[minibuffer_id] = (triggernum==0) ? (0x1<<4) : (0x1<<24); // Beam or Window trigger
	timefileout_Time[minibuffer_id]  = fileout_StartTimeSec*SEC_TO_NS            // run start time
		 + currenteventtime + minibuffer_start_time_ns; // this beam trigger offset + this minibuffer offset
	
}

void WCSimAnalysis::ConstructEmulatedPmtDataReadout(){
#ifdef EMULATED_OUT_VERBOSE
	cout<<"constructing emulated PMT data"<<endl;
#endif
	// fill all the non-minibuffer info and reset the minibuffers
	if(emulated_pmtdata_readout.size()<num_adc_cards){
		cout<<"constructing vector of "<<num_adc_cards<<" CardData objects"
			<<" and temp databuffer vector with "<<full_buffer_size<<" element data arrays"<<endl;
		emulated_pmtdata_readout= std::vector<CardData>(num_adc_cards);
		temporary_databuffers = std::vector< std::vector<uint16_t> >(num_adc_cards);
	}
	
	// Calculate the timing values
	///////////////////////////////////
	
	// the first LastSync value is ~ one step above the first StartCount value
	// (StartCount values are pre-calculated in createrawfile.cxx)
	// at the start of each new sequence_id, LastSync steps up by a large step
	// to calculate a suitable step size, we need more details:
	
	// LastSync is dominated by a series of steps, each step being around 1e9 (700-1200 x 10*6)
	// in size and with steps occurring at each new sequence_id.
	// Step sizes are distributed in discrete regions, separated by 125e6,
	// spread roughly gaussian around a centre of 950e6.
	// Each discrete region contains step sizes varying by about ~1000.
	// There's also a group of small steps, which generally occur in close proximity to the large ones,
	// of around ~120e6 (121e6 - 127e6) in size.
	// After subtracting StartCount the size of small steps reduces in spread to ~10k (124.98e6 to 124.99e6).
	// The spread within the discrete regions is unchanged:
	// -- e.g. without StartCount removed: (876.8082 - 876.8086 e6)
	// -- e.g. with    StartCount removed: (875.2140 - 875.2145 e6)
	// Finally between steps there are the regular decrements between steps of <9e6, equal
	// to the steps in StartTimeNSec, on top of that the same noise as in StartCount, 
	// and finally some small residual noise ~+-500.
	// (These regular step sizes are negative once StartCount is subtracted (~-15000+-1500))
	
	// for simplicity we'll ignore the presence of small steps for now.
	TRandom3 R;
	uint64_t largespread = static_cast<int64_t>(R.Gaus(950e6,150e6));  // calculate main step size
	uint64_t roundedval = 125e6*(largespread/int64_t(125e6));          // discretize the step size to a region
	int64_t smallspread = static_cast<int64_t>(R.Gaus(0,250));         // calculate spread within the region
	uint64_t stepsize = roundedval + smallspread;
	runningsteps += stepsize;                                          // add to the running total
	
	// TriggerCounts has 40 entries per readout (one per minibuffer).
	// Within each entry, the 40 values have a linear rising component and a series of ~4-5 steps,
	// which are irregular both in amplitude and spacing.
	// At each new SequenceID, all values step up by roughly the same amount as LastSync.
	// Although the step sizes are not exactly the same, over many sequence_ids they track
	// well enough that (TriggerCounts - LastSync) has no long term rise. The positions of steps
	// also changes and the amplitude of the sawtooth changes.
	int numsteps = R.Gaus(4,1);
	std::map<int,int> steps;
	for(int stepi=0; stepi<numsteps; stepi++){
		int asteppos = R.Gaus((stepi+1)*(minibuffers_per_fullbuffer/(numsteps+1)),3);
		int astepsize = R.Uniform(50e6,250e6);  // step sizes are dominated by several randomly placed peaks
		steps.emplace(asteppos,astepsize); // in this range of varying frequency, but we'll ignore for now
	}
	std::vector<ULong64_t> theTriggerCounts(minibuffers_per_fullbuffer);
	int runningval=0;
	for(int mbi=0; mbi<minibuffers_per_fullbuffer; mbi++){
		runningval = runningval + R.Gaus(8332e3,6e3);
		if(steps.count(mbi)) runningval+=steps.at(mbi);
		theTriggerCounts.at(mbi) = runningval;
	}
	
	/////////////////////////////
	// fill the PMTData entries for each card.
	
	for(int cardi=0; cardi<num_adc_cards; cardi++){
		auto&& acard = emulated_pmtdata_readout.at(cardi);
		
		// fixed values
		// these need to be set first as they're used by CardData::Reset
		acard.Channels = channels_per_adc_card;
		acard.TriggerNumber = minibuffers_per_fullbuffer;
		acard.FullBufferSize = full_buffer_size;           // trigger window sizes and sampling rate
		acard.Eventsize = emulated_event_size;             // wcsimanalysis::emulated_event_size etc.
		acard.BufferSize = buffer_size;                    // calculated in utilityfuncs based on
		
		acard.Reset();
		acard.CardID = cardi;
		
		// mostly randomly generated timing references
		acard.StartTimeSec = fileout_StartTimeSec;         // start time of run (read from config file)
		acard.StartTimeNSec = StartTimeNSecVals.at(cardi); // steady decreases
		acard.StartCount = StartCountVals.at(cardi);       // periodic spikes
		int64_t smallresidual = static_cast<int64_t>(R.Uniform(-500,500)); // varies per value
		acard.LastSync = acard.StartCount + runningsteps - StartTimeNSecVals.at(cardi)/8 + smallresidual;
		acard.LastSync = (acard.LastSync/4)*4;             // always a mutiple of 4
		
		// randomly generated TriggerCounts
		// In each of the 16 entries with the same sequence_id (one for each card) the array of
		// TriggerCounts is roughly the same, creating a sawtooth with only small variations
		std::vector<ULong64_t> TriggerCountsWithNoise(minibuffers_per_fullbuffer);
		for(int mbi=0; mbi<minibuffers_per_fullbuffer; mbi++){
			int thenoise = R.Gaus(5e6,2e6);
			TriggerCountsWithNoise.at(mbi) = theTriggerCounts.at(mbi) + thenoise;
		}
		acard.TriggerCounts = TriggerCountsWithNoise;
		
		// actually set per event
		acard.SequenceID = sequence_id;
		//acard.Data.assign(full_buffer_size,0); not necessary since we're using the temporary buffers
		temporary_databuffers.at(cardi).assign(full_buffer_size,0);
		
		//  Remaining event data to be filled while processing triggers:
		//  -----------------------------------------------------------
		//  TriggerCounts.assign(TriggerNumber,BOGUS_UINT64);   // start counts of each minibuffer
		//  Rates.assign(Channels,BOGUS_UINT32);                // avg pulse rate on each PMT
		//  Data.assign(FullBufferSize,BOGUS_UINT16);           // 160k datapoints
	}
	
}

void WCSimAnalysis::FillEmulatedPMTData(){
#ifdef EMULATED_OUT_VERBOSE
	cout<<"Filling Emulated PMT Data"<<endl;
#endif
	//fileout_TriggerCounts = new ULong64_t[minibuffers_per_fullbuffer];  // must be >= fileout_TriggerNumber
	//fileout_Rates         = new UInt_t[channels_per_adc_card];           // must be >= fileout_Channels
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
	
	// first, add noise to the temporary waveforms
	AddNoiseToWaveforms();
	
	// next, shuffle your library. this also copies temporary_databuffers into emulated_pmtdata_readout
	RiffleShuffle(false); // XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX 
	
#ifdef EMULATED_OUT_VERBOSE
	cout<<"Filling PMTData tree"<<endl;
#endif
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
		tPMTData->SetBranchAddress("Data",fileout_Data);
		// ^ why?! without this first set of full buffers all read garbage... 
		
		// draw for debug
		if(full_buffer_size!=acard.Data.size()){
			cout<<"full_buffer_size="<<full_buffer_size
				<<", fileout_Data.size()="<<acard.Data.size()<<endl;
		}
		static TGraph* buffergraph = new TGraph(full_buffer_size);
		static TCanvas* buffercanvas = new TCanvas();
		int i=0;
		for(auto it=acard.Data.begin(); it!=acard.Data.end(); it++){
			auto theval = *it;
			buffergraph->SetPoint(i,i,theval); i++;
			//if(theval!=0) cout<<"("<<i<<","<<theval<<")"<<", ";
		}
		//cout<<endl;
		buffercanvas->cd();
		buffergraph->SetTitle(TString::Format("Card %d",fileout_CardID));
		buffergraph->Draw("alp");
		buffercanvas->Update();
		//std::this_thread::sleep_for (std::chrono::seconds(5));
		//gPad->WaitPrimitive();
		
		tPMTData->Fill();
	}
	
	// Fill the hefty timing file. Branches are populated already
	cout<<"Filling heftdb tree: More branch has values: ";
	for(int i=0; i<minibuffers_per_fullbuffer; i++){
		cout<<timefileout_More[i]<<", ";
	}
	cout<<endl;
	theftydb->Fill();
	
//	// draw for debug
//	int samples_per_channel = full_buffer_size/channels_per_adc_card;
//	static TGraph* buffergraph = new TGraph(samples_per_channel);
//	static TCanvas* buffercanvas = new TCanvas();
//	for(int i=0; i<samples_per_channel; i++){
//		buffergraph->SetPoint(i,(float)i,(float)fileout_Data[i]);
//	}
//	buffergraph->SetLineColor(kRed);
//	buffercanvas->cd();
//	buffergraph->Draw("alp");
//	buffercanvas->Update();
//	gPad->WaitPrimitive();
//	
//	// retrieve the version in the tree
//	std::vector<uint16_t> z(full_buffer_size,0);
//	tPMTData->SetBranchAddress("Data",z.data());
//	cout<<"retrieving Data[] from tree entry "<<(tPMTData->GetEntries()-1)<<endl;
//	tPMTData->GetEntry(tPMTData->GetEntries()-1);
//	for(int i=0; i<samples_per_channel; i++){
//		buffergraph->SetPoint(i,(float)i,(float)z.at(i));
//	}
//	buffergraph->SetLineColor(kBlue);
//	buffercanvas->cd();
//	buffergraph->Draw("alp");
//	buffercanvas->Update();
//	gPad->WaitPrimitive();
//	tPMTData->SetBranchAddress("Data",fileout_Data);
//	z.clear();
	
}

void WCSimAnalysis::AddNoiseToWaveforms(){
#ifdef EMULATED_OUT_VERBOSE
	cout<<"Adding noise to waveforms"<<endl;
#endif
//	// add noise to the temporary buffers
//	for(int cardi=0; cardi<num_adc_cards; cardi++){
//		// get the temporary Data array for this card, which currently only has pulses in
//		auto& temp_card_buffer = temporary_databuffers.at(cardi);
//		// Data arrays for channels on this card are concatenated: loop over channels
//		for(int channelnum=0; channelnum<channels_per_adc_card; channelnum++){
//			int channeloffset = channelnum * (full_buffer_size / channels_per_adc_card);
//			// minibuffers for this channel are concatenated: loop over minibuffers
//			for(int minibufi=0; minibufi<minibuffers_per_fullbuffer; minibufi++){
//				// calculate the offset of the minibuffer, for this channel, within the temporary buffer
//				int minibufferoffset = minibufi*minibuffer_datapoints_per_channel;
//				
//				// now generate a vector of noise to add
//				std::vector<uint16_t> noisevector(minibuffer_datapoints_per_channel);
//				// Each minibuffer has an offset of ~330 +- 20 ADC counts, distributed... well..
//				// there may be an underlying sine wave of varying amplitude, which is sorta kinda uniform..?
//				uint16_t mboffset = static_cast<uint16_t>(R.Uniform(310,350));
//				// then within the minibuffer a spread of gaussian noise +-5 ADC counts
//				for(int avali=0; avali<minibuffer_datapoints_per_channel; avali++){
//					noisevector.at(avali)=mboffset+static_cast<uint16_t>(R.Gaus(0,2));
//				}
//				
//				if(cardi==0&&channelnum==0&&(minibufi==11||minibufi==12)){
//					cout<<"card "<<cardi<<", channel "<<channelnum<<", minibuffer "<<minibufi
//						<<" noise values: "<<endl;
//					for(int noiseval : noisevector){
//						cout<<noiseval<<", ";
//					}
//					cout<<endl;
//				}
//				
//				// get the data array for this card, and an iterator to the appropriate start point
//				auto minibuffer_start = temp_card_buffer.begin() + channeloffset + minibufferoffset;
//				
//				if(cardi==0&&channelnum==0&&(minibufi==11||minibufi==12)){
//					cout<<"card "<<cardi<<", channel "<<channelnum<<", minibuffer "<<minibufi
//						<<", channeloffset = "<<channeloffset
//						<<", minibufferoffset = "<<minibufferoffset
//						<<", mb_dp_p_ch = "<<minibuffer_datapoints_per_channel
//						<<endl;
//				}
//				
////				cout<<"card "<<cardi<<", channel "<<channelnum<<", minibuffer "<<minibufi<<" without noise: "<<endl;
////				for(auto it = minibuffer_start; it!=(minibuffer_start+minibuffer_datapoints_per_channel); it++){
////					cout<<(*(it))<<", ";
////				}
//				
//				// add the noise to the correct region of the full data array
//				std::transform(minibuffer_start, minibuffer_start+minibuffer_datapoints_per_channel,
//					 noisevector.begin(), minibuffer_start, std::plus<uint16_t>() );
//				
////				cout<<"card "<<cardi<<", channel "<<channelnum<<", minibuffer "<<minibufi<<" with noise: "<<endl;
////				for(auto it = minibuffer_start; it!=(minibuffer_start+minibuffer_datapoints_per_channel); it++){
////					cout<<(*(it))<<", ";
////				}
////				cout<<endl;
//			}
//		}
//	}
	
//////
//	Alt method: not significantly faster
//	// add noise to the temporary buffers
	for(int cardi=0; cardi<num_adc_cards; cardi++){
		// get the temporary Data array for this card, which currently only has pulses in
		auto& temp_card_buffer = temporary_databuffers.at(cardi);
		
		// loop over all the minibuffers in the full buffer
		uint16_t mboffset=0;
		for(int samplei=0; samplei<temp_card_buffer.size(); samplei++){
			// Each minibuffer has an offset of ~330 +- 20 ADC counts, distributed... well..
			// there may be an underlying sine wave of varying amplitude, which is sorta kinda uniform..?
			if((samplei%(full_buffer_size/minibuffers_per_fullbuffer))==0){
				// if new minibuffer, generate a new minibuffer offset
				mboffset = static_cast<uint16_t>(R.Uniform(310,350));
			}
			
			// within the minibuffer there's a spread of gaussian noise +-5 ADC counts
			//temp_card_buffer.at(samplei) += mboffset + static_cast<uint16_t>(R.Gaus(0,2));
		}
	}
}

void WCSimAnalysis::RiffleShuffle(bool do_shuffle){
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
#ifdef EMULATED_OUT_VERBOSE
	cout<<"Shuffling Data arrays"<<endl;
#endif
	
	int channel_buffer_size = (full_buffer_size / channels_per_adc_card);
	for(int cardi=0; cardi<num_adc_cards; cardi++){
		auto& acard = emulated_pmtdata_readout.at(cardi);
		auto& card_buffer = acard.Data;
		auto& temp_card_buffer = temporary_databuffers.at(cardi);
		
		if(do_shuffle){
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
		} else {
			for(int i=0; i<temp_card_buffer.size(); i++) card_buffer.at(i)=temp_card_buffer.at(i);
		}
	}
	// finally, ask a player to cut the deck
}
