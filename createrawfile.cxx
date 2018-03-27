/* vim:set noexpandtab tabstop=4 wrap */
// #######################################################################
#include <regex>
#include <stdlib.h> 	/* atoi */

// CREATE+OPEN RAW FORMAT OUTPUT FILE
// ==================================
void WCSimAnalysis::LoadOutputFiles(){
	TTree* currenttree = t->GetTree();
	TFile* currentfile = currenttree->GetCurrentFile();
	std::string wcsimfilename = currentfile->GetName(); // wcsim_0.AAA.B.root
	
	// use regexp to pull out the file numbers - we'll call this 'run' and 'subrun'
	std::match_results<string::const_iterator> submatches;
	// filename is of the form "wcsim_0.AAA.B.root"
	//std::regex theexpression (".*/[^0-9]+\\.([0-9]+)\\.([0-9]+)\\.root");
	std::regex theexpression (".*/?[^\\.]+\\.([0-9]+)\\.?([0-9]+)?\\.root");
	//cout<<"matching regex for filename "<<wcsimfilename<<endl;
	std::regex_match (wcsimfilename, submatches, theexpression);
	std::string submatch = (std::string)submatches[0]; // match 0 is whole match
	int rawfilerun, rawfilesubrun, rawfilepart;
	if(submatch==""){ 
		cerr<<"unrecognised input file pattern: "<<wcsimfilename
			<<", will set rawfilerun=0, rawfilesubrun=0, rawfilepart=0"<<endl;
		//return;
		rawfilerun=0;
		rawfilesubrun=0;
		rawfilepart = 0;
	} else {
		submatch = (std::string)submatches[1];
		//cout<<"extracted submatch is "<<submatch<<endl;
		rawfilerun = atoi(submatch.c_str());
		submatch = (std::string)submatches[2];
		rawfilesubrun = atoi(submatch.c_str());
		rawfilepart = 0;
	}
	placeholder_date_ns = 1511965741;
	
	std::string rawfilename="RAWDataR"+to_string(rawfilerun)+"S"+to_string(rawfilesubrun)+"p"+to_string(rawfilepart)+".root";
	cout<<"creating RAW output file "<<rawfilename<<endl;
	rawfileout = new TFile(rawfilename.c_str(),"RECREATE");
	//*----------------------------------------------------------------------------*
	tPMTData      = new TTree("PMTData","");
	fileout_TriggerCounts = new ULong64_t[minibuffers_per_fullbuffer];     // must be >= TriggerNumber
	fileout_Rates         = new UInt_t[channels_per_adc_card];              // must be >= Channels
	fileout_Data          = new UShort_t[full_buffer_size];                // must be >= FullBufferSize
	//*............................................................................*
	TBranch *bLastSync       = tPMTData->Branch("LastSync", &fileout_LastSync);
	TBranch *bSequenceID     = tPMTData->Branch("SequenceID", &fileout_SequenceID);
	TBranch *bStartTimeSec   = tPMTData->Branch("StartTimeSec", &fileout_StartTimeSec);
	TBranch *bStartTimeNSec  = tPMTData->Branch("StartTimeNSec", &fileout_StartTimeNSec);
	TBranch *bStartCount     = tPMTData->Branch("StartCount", &fileout_StartCount);
	TBranch *bTriggerNumber  = tPMTData->Branch("TriggerNumber", &fileout_TriggerNumber);
	TBranch *bTriggerCounts  = tPMTData->Branch("TriggerCounts", &fileout_TriggerCounts, 
		"fileout_TriggerCounts[TriggerNumber]/l");
	TBranch *bCardID         = tPMTData->Branch("CardID", &fileout_CardID);
	TBranch *bChannels       = tPMTData->Branch("Channels", &fileout_Channels);
	TBranch *bRates          = tPMTData->Branch("Rates", &fileout_Rates, 
		"fileout_Rates[Channels]/i");
	TBranch *bBufferSize     = tPMTData->Branch("BufferSize", &fileout_BufferSize);
	TBranch *bEventsize      = tPMTData->Branch("Eventsize", &fileout_Eventsize);
	TBranch *bFullBufferSize = tPMTData->Branch("FullBufferSize", &fileout_FullBufferSize);
	TBranch *bData           = tPMTData->Branch("Data", &fileout_Data, "fileout_Data[FullBufferSize]/s");
	//*----------------------------------------------------------------------------*
	tRunInformation = new TTree("RunInformation","");
	//*............................................................................*
	TBranch *bInfoTitle = tRunInformation->Branch("InfoTitle", &fileout_InfoTitle);
	TBranch *bInfoMessage = tRunInformation->Branch("InfoMessage", &fileout_InfoMessage);
	//*----------------------------------------------------------------------------*
	tTrigData       = new TTree("TrigData","");
	fileout_EventIDs        = new UShort_t[MAXEVENTSIZE];    // must be >= Eventsize
	fileout_EventTimes      = new ULong64_t[MAXEVENTSIZE];   // must be >= Eventsize
	fileout_TriggerMasks    = new UInt_t[MAXTRIGGERSIZE];    // must be >= TriggerSize
	fileout_TriggerCounters = new UInt_t[MAXTRIGGERSIZE];    // muse be >= TriggerSize
	//*............................................................................*
	TBranch *bFirmwareVersion = tTrigData->Branch("FirmwareVersion", &fileout_FirmwareVersion);
	TBranch *bSequenceID2     = tTrigData->Branch("SequenceID", &fileout_SequenceID);
	TBranch *bEventsize2      = tTrigData->Branch("Eventsize", &fileout_Eventsize);
	TBranch *bTriggerSize     = tTrigData->Branch("TriggerSize", &fileout_TriggerSize);
	TBranch *bFIFOOverflow    = tTrigData->Branch("FIFOOverflow", &fileout_FIFOOverflow);
	TBranch *bDriverOverfow   = tTrigData->Branch("DriverOverfow", &fileout_DriverOverfow);
	TBranch *bEventIDs        = tTrigData->Branch("EventIDs", &fileout_EventIDs, 
		"fileout_EventIDs[Eventsize]/s");
	TBranch *bEventTimes      = tTrigData->Branch("EventTimes", &fileout_EventTimes, 
		"fileout_EventTimes[Eventsize]/l");
	TBranch *bTriggerMasks    = tTrigData->Branch("TriggerMasks", &fileout_TriggerMasks, 
		"fileout_TriggerMasks[TriggerSize]/I");
	TBranch *bTriggerCounters = tTrigData->Branch("TriggerCounters", &fileout_TriggerCounters, "fileout_TriggerCounters[TriggerSize]/I");
	//*----------------------------------------------------------------------------*
	tCCData = new TTree("CCData","");
	//*............................................................................*
	TBranch *bTrigger   = tCCData->Branch("Trigger", &fileout_Trigger);         // readout number
	TBranch *bOutNumber = tCCData->Branch("OutNumber", &fileout_OutNumber);     // num hits this event/readout
	TBranch *bType      = tCCData->Branch("Type", &fileout_Type);               // card type string, "TDC", "ADC"
	TBranch *bValue     = tCCData->Branch("Value", &fileout_Value);             // see below
	TBranch *bSlot      = tCCData->Branch("Slot", &fileout_Slot);               // card position in crate
	TBranch *bChannel   = tCCData->Branch("Channel", &fileout_Channel);         // channel in card
	TBranch *bTimeStamp = tCCData->Branch("TimeStamp", &fileout_TimeStamp);     // see below
	//*............................................................................*
	/*
	The MRD process is:
		1. Trigger card sends common start to MRD cards
		2. A timer is started on all channels.
		3. When a channel receives a pulse, the timer stops. XXX: only the first pulse is recorded!XXX
		4. After all channels either record hits or time out (currently 4.2us) everything is read out.
		   A timestamp is created at time of readout. Note this is CLOSE TO, BUT NOT EQUAL TO the trigger
		   time. (Probably ~ trigger time + timeout...)
		5. Channels with a hit will have an entry created with 'Value' = clock ticks between the common
		   start and when the hit arrived. Channels that timed out have no entry.
	
	Timestamp is a UTC [ns] timestamp of when the readout ended. To correctly map Value to an actual time
	one would need to match the MRD timestamp to the trigger card timestamp (which will be more accurate)
	then add Value * MRD_CLOCK_TICK_NS to the Trigger time. 
	
	Only the first pulse will be recorded .... IS THIS OK? XXX
	What about pre-trigger? - common start issued by trigger is from Beam not on NDigits, so always pre-beam...?
	TDC records with resolution 4ns = 1 sample per ns? or round times to nearest/round down to 4ns?
	TDCRes https://github.com/ANNIEDAQ/ANNIEDAQ/blob/master/configfiles/TDCreg
	*/
	
	gROOT->cd();
	
//	RunInformation: fixed 11 entries per run;
//	PMTData: one entry per trigger readout (6128 entries);
//	TrigData: one entry per trigger (383 entries);
//	CCData: one entry per MRD trigger (3758 entries);
}

void WCSimAnalysis::FillEmulatedRunInformation(){
	// RunInformation tree in RAW file is essentially a map of strings.
	// It stores (at present) 11 entries, each with a 'key' in the InfoTitle branch, 
	// and a corresponding xml object string in the InfoMessage branch - 
	// the xml string may represent multiple variables.
	// Since this is largely irrelevant, we'll just use template values and modify if necessary.
	
	std::vector<std::string> InfoTitles{"ToolChainVariables","InputVariables","PostgresVariables",
	"SlackBotVariables","HVComs","TriggerVariables","NetworkReceiveDataVariables","LoggerVariables",
	"RootDataRecorderVariables","MonitoringVariables","MRDVariables"};
	std::vector<std::string>* TemplateInfoMessages = GetTemplateRunInfo();
	
	for(int entryi=0; entryi<InfoTitles.size(); entryi++){
		fileout_InfoTitle=InfoTitles.at(entryi);
		fileout_InfoMessage=TemplateInfoMessages->at(entryi);
		tRunInformation->Fill();
	}
}
