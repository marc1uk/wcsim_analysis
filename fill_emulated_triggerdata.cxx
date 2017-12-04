/* vim:set noexpandtab tabstop=4 wrap */
// #######################################################################

void WCSimAnalysis::FillEmulatedTrigData(){
	
	// fileout_EventIDs        = new UShort_t[MAXEVENTSIZE];    // must be >= EventSize
	// fileout_EventTimes      = new ULong64_t[MAXEVENTSIZE];   // must be >= EventSize
	// fileout_TriggerMasks    = new UInt_t[MAXTRIGGERSIZE];    // must be >= TriggerSize
	// fileout_TriggerCounters = new UInt_t[MAXTRIGGERSIZE];    // muse be >= TriggerSize
	//*............................................................................*
	// TBranch *bFirmwareVersion = tTrigData->Branch("FirmwareVersion", &fileout_FirmwareVersion);
	// TBranch *bSequenceID      = tTrigData->Branch("SequenceID", &fileout_SequenceID);
	// TBranch *bEventSize       = tTrigData->Branch("EventSize", &fileout_Eventsize);
	// TBranch *bTriggerSize     = tTrigData->Branch("TriggerSize", &fileout_TriggerSize);
	// TBranch *bFIFOOverflow    = tTrigData->Branch("FIFOOverflow", &fileout_FIFOOverflow);
	// TBranch *bDriverOverfow   = tTrigData->Branch("DriverOverfow", &fileout_DriverOverfow);
	// TBranch *bEventIDs        = tTrigData->Branch("EventIDs", &fileout_EventIDs);
	// TBranch *bEventTimes      = tTrigData->Branch("EventTimes", &fileout_EventTimes);
	// TBranch *bTriggerMasks    = tTrigData->Branch("TriggerMasks", &fileout_TriggerMasks);
	// TBranch *bTriggerCounters = tTrigData->Branch("TriggerCounters", &fileout_TriggerCounters);
	
	fileout_FirmwareVersion = 0;     // FIXME represent WCSim version instead
	fileout_SequenceID = 0;          // FIXME
	fileout_FIFOOverflow = 0;        // FIXME
	fileout_DriverOverfow = 0;       // FIXME
	
	//fileout_EventSize = ???;       // FIXME there's already a fileout_EventSize for PMTData. Same?
	fileout_TriggerSize = 0;         // FIXME
	
	if(fileout_Eventsize>MAXEVENTSIZE){
		MAXEVENTSIZE=fileout_Eventsize;
		if(fileout_EventIDs) delete[] fileout_EventIDs;
		fileout_EventIDs = new UShort_t[MAXEVENTSIZE];
		if(fileout_EventTimes) delete[] fileout_EventTimes;
		fileout_EventTimes = new ULong64_t[MAXEVENTSIZE];
	}
	if(fileout_TriggerSize>MAXTRIGGERSIZE){
		MAXTRIGGERSIZE=fileout_TriggerSize;
		if(fileout_TriggerMasks) delete[] fileout_TriggerMasks;
		fileout_TriggerMasks = new UInt_t[MAXTRIGGERSIZE];
		if(fileout_TriggerCounters) delete[] fileout_TriggerCounters;
		fileout_TriggerCounters = new UInt_t[MAXTRIGGERSIZE];
	}
	
	if(fileout_EventIDs==nullptr) fileout_EventIDs = new UShort_t[0];
	if(fileout_EventTimes==nullptr) fileout_EventTimes = new ULong64_t[0];
	if(fileout_TriggerMasks==nullptr) fileout_TriggerMasks = new UInt_t[0];
	if(fileout_TriggerCounters==nullptr) fileout_TriggerCounters = new UInt_t[0];
	
	//fileout_EventIDs = ???;          // FIXME
	//fileout_EventTimes = ???;        // FIXME
	//fileout_TriggerMasks = ???;      // FIXME
	//fileout_TriggerCounters = ???;   // FIXME
	
	tTrigData->Fill();
}
