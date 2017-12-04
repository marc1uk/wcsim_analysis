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
	if(submatch==""){ cerr<<"unrecognised input file pattern: "<<wcsimfilename<<endl; return; }
	submatch = (std::string)submatches[1];
	//cout<<"extracted submatch is "<<submatch<<endl;
	int rawfilerun = atoi(submatch.c_str());
	submatch = (std::string)submatches[2];
	int rawfilesubrun = atoi(submatch.c_str());
	int rawfilepart = 0;
	placeholder_date_ns = 1511965741;
	
	std::string rawfilename="RAWDataR"+to_string(rawfilerun)+"S"+to_string(rawfilesubrun)+"p"+to_string(rawfilepart)+".root";
	cout<<"creating RAW output file "<<rawfilename<<endl;
	rawfileout = new TFile(rawfilename.c_str(),"RECREATE");
	//*----------------------------------------------------------------------------*
	tPMTData      = new TTree("PMTData","");
	fileout_TriggerCounts = new ULong64_t[minibuffers_per_fullbuffer];     // must be >= TriggerNumber
	fileout_Rates         = new Int_t[channels_per_adc_card];              // must be >= Channels
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
		"fileout_Rates[Channels]/I");
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
	TBranch *bEventsize2       = tTrigData->Branch("Eventsize", &fileout_Eventsize);
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
	TBranch *bTrigger   = tCCData->Branch("Trigger", &fileout_Trigger);
	TBranch *bOutNumber = tCCData->Branch("OutNumber", &fileout_OutNumber);
	TBranch *bType      = tCCData->Branch("Type", &fileout_Type);
	TBranch *bValue     = tCCData->Branch("Value", &fileout_Value);
	TBranch *bSlot      = tCCData->Branch("Slot", &fileout_Slot);
	TBranch *bChannel   = tCCData->Branch("Channel", &fileout_Channel);
	TBranch *bTimeStamp = tCCData->Branch("TimeStamp", &fileout_TimeStamp);
	//*............................................................................*
	
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
