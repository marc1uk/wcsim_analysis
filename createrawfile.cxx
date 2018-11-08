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
	
	std::string rawfilename="RAWDataR"+to_string(rawfilerun)+"S"+to_string(rawfilesubrun)+"p"+to_string(rawfilepart)+".root";
	cout<<"creating RAW output file "<<rawfilename<<endl;
	rawfileout = new TFile(rawfilename.c_str(),"RECREATE");
	//*----------------------------------------------------------------------------*
	tPMTData                 = new TTree("PMTData","");
	fileout_TriggerCounts    = new ULong64_t[minibuffers_per_fullbuffer];     // must be >= TriggerNumber
	fileout_Rates            = new UInt_t[channels_per_adc_card];             // must be >= Channels
	fileout_Data             = new UShort_t[full_buffer_size];                // must be >= FullBufferSize
	//*............................................................................*
	TBranch *bLastSync       = tPMTData->Branch("LastSync", &fileout_LastSync);
	TBranch *bSequenceID     = tPMTData->Branch("SequenceID", &fileout_SequenceID);
	TBranch *bStartTimeSec   = tPMTData->Branch("StartTimeSec", &fileout_StartTimeSec);
	TBranch *bStartTimeNSec  = tPMTData->Branch("StartTimeNSec", &fileout_StartTimeNSec);
	TBranch *bStartCount     = tPMTData->Branch("StartCount", &fileout_StartCount);
	TBranch *bTriggerNumber  = tPMTData->Branch("TriggerNumber", &fileout_TriggerNumber);
	TBranch *bTriggerCounts  = tPMTData->Branch("TriggerCounts", &fileout_TriggerCounts, 
		"TriggerCounts[TriggerNumber]/l");
	TBranch *bCardID         = tPMTData->Branch("CardID", &fileout_CardID);
	TBranch *bChannels       = tPMTData->Branch("Channels", &fileout_Channels);
	TBranch *bRates          = tPMTData->Branch("Rates", &fileout_Rates, 
		"Rates[Channels]/i");
	TBranch *bBufferSize     = tPMTData->Branch("BufferSize", &fileout_BufferSize);
	TBranch *bEventsize      = tPMTData->Branch("Eventsize", &fileout_Eventsize);
	TBranch *bFullBufferSize = tPMTData->Branch("FullBufferSize", &fileout_FullBufferSize);
	TBranch *bData           = tPMTData->Branch("Data", &fileout_Data, "Data[FullBufferSize]/s");
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
	TBranch *bEventsize2      = tTrigData->Branch("EventSize", &fileout_Eventsize);
	TBranch *bTriggerSize     = tTrigData->Branch("TriggerSize", &fileout_TriggerSize);
	TBranch *bFIFOOverflow    = tTrigData->Branch("FIFOOverflow", &fileout_FIFOOverflow);
	TBranch *bDriverOverfow   = tTrigData->Branch("DriverOverfow", &fileout_DriverOverfow);
	TBranch *bEventIDs        = tTrigData->Branch("EventIDs", &fileout_EventIDs, 
		"EventIDs[EventSize]/s");
	TBranch *bEventTimes      = tTrigData->Branch("EventTimes", &fileout_EventTimes, 
		"EventTimes[EventSize]/l");
	TBranch *bTriggerMasks    = tTrigData->Branch("TriggerMasks", &fileout_TriggerMasks, 
		"TriggerMasks[TriggerSize]/I");
	TBranch *bTriggerCounters = tTrigData->Branch("TriggerCounters", &fileout_TriggerCounters, "TriggerCounters[TriggerSize]/I");
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
	
	Timestamp is a UTC [MILLISECONDS] timestamp of when the readout ended. 
	To correctly map Value to an actual time one would need to match the MRD timestamp 
	to the trigger card timestamp (which will be more accurate)
	then add Value * MRD_NS_PER_SAMPLE to the Trigger time. 
	
	Only the first pulse will be recorded .... IS THIS OK? XXX
	What about pre-trigger? - common start issued by trigger is from Beam not on NDigits, so always pre-beam...?
	TDC records with resolution 4ns = 1 sample per ns? or round times to nearest/round down to 4ns?
	TDCRes https://github.com/ANNIEDAQ/ANNIEDAQ/blob/master/configfiles/TDCreg
	*/
	
	//*----------------------------------------------------------------------------*
	// Since we can't properly synthesize all the timing variables we'll directly create
	// the simplified output Hefty Timing file.
	std::string timingfilename="DataR"+to_string(rawfilerun)+"S"+to_string(rawfilesubrun)+"p"+to_string(rawfilepart)+"_timing.root";
	cout<<"creating RAW timing file "<<timingfilename<<endl;
	timingfileout = new TFile(timingfilename.c_str(),"RECREATE");
	timingfileout->cd();
	//*----------------------------------------------------------------------------*
	theftydb                  = new TTree("heftydb","");
	timefileout_Time          = new ULong_t[minibuffers_per_fullbuffer];
	timefileout_Label         = new Int_t[minibuffers_per_fullbuffer];
	timefileout_TSinceBeam    = new Long_t[minibuffers_per_fullbuffer];
	timefileout_More          = new Int_t[minibuffers_per_fullbuffer];
	//*----------------------------------------------------------------------------*
	theftydb->Branch("SequenceID", &fileout_SequenceID);
	theftydb->Branch("Time", &timefileout_Time, TString::Format("Time[%d]/l",minibuffers_per_fullbuffer));
	theftydb->Branch("Label", &timefileout_Label, TString::Format("Label[%d]/I",minibuffers_per_fullbuffer));
	theftydb->Branch("TSinceBeam", &timefileout_TSinceBeam, TString::Format("TSinceBeam[%d]/L",minibuffers_per_fullbuffer));
	theftydb->Branch("More", &timefileout_More, TString::Format("More[%d]/I",minibuffers_per_fullbuffer));
	
	gROOT->cd();
	
	//===============================
	// We can derive some of the output values here, as they're just randomly generated
	
	// First, StartTimeSec represents the unix seconds of the start of the run.
	// This is read from the options file in utilityfuncs - we need to convert date string to unixns
	// first convert config file string to parts
	int hh, mm, ss, MM, DD, YYYY;
	sscanf(startDate.c_str(), "%d/%d/%d %d:%d:%d", &DD, &MM, &YYYY, &hh, &mm, &ss);
	// combine parts into a time structure
	struct std::tm runstart;
	runstart.tm_year = YYYY;
	runstart.tm_mon = MM;
	runstart.tm_mday = DD;
	runstart.tm_hour = hh;
	runstart.tm_min = mm;
	runstart.tm_sec = ss;
	// finally convert time structure to unix seconds
	time_t runstarttime = mktime(&runstart);
	// use runstarttime = time(NULL);  to get current time.
	fileout_StartTimeSec = static_cast<Int_t>(runstarttime);
	
	// StartCount is a periodic, spiky kinda cosine ish thing, which repeats with each sequence_id
	// amplitude of variation ~10e6, offset ~1e15, first entry seems to be the highest.
	// let's generate a random periodic signal once off for this run:
	int64_t StartCountOffset;
	int64_t anum;
	int64_t aref=100000000000000;  // offset should be of this order
	do{                            // generate a suitable random number
		StartCountOffset=static_cast<int64_t>(R.Gaus(5e14,2e14));
		anum=StartCountOffset/aref;
		//printf("%lu, %lu\n",StartCountOffset, anum);
	} while ((anum==0)||(anum>9));
	
	StartCountVals.resize(num_adc_cards);
	for(int i=0; i<num_adc_cards; i++){
		double randp = static_cast<double>(R.Gaus(-0.1,0.08));
		double cosval = -0.9*pow(cos(TMath::Pi()*(double(i)/10.)),2.);
		if((randp+cosval)>0||(randp+cosval)<-1.) randp*=-1.;  // constrain to -1->0
		StartCountVals.at(i)= static_cast<int64_t>(randp+cosval*10000000.) + StartCountOffset;
	}
	
	// StartTimeNSec is a rising sawtooth, monotonically stepping with a size of ~117+-3 x10*3
	// step sizes are discretized with separation 500.
	// even the deviations from equal step sizes appear to be periodic, so we can prepare the values now
	// starting values appear to be ~1-10 e8
	StartTimeNSecVals.resize(num_adc_cards);
	int64_t firstval = static_cast<int64_t>(R.Gaus(5e8,1e8));
	StartTimeNSecVals.at(0) = (firstval/500)*500;
	for(int i=1; i<num_adc_cards; i++){
		uint64_t randp = static_cast<uint64_t>(R.Gaus(117e3,1.5e3));
		uint64_t roundedval=(randp/500)*500;
		StartTimeNSecVals.at(i) = StartTimeNSecVals.at(i-1) + roundedval;
	}
	
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
