/* vim:set noexpandtab tabstop=4 wrap */
#ifndef WCSimAnalysisClass
#define WCSimAnalysisClass

#define FILE_VERSION 2

#include "TROOT.h"
#include "TRint.h"
#include "TSystem.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TApplication.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TString.h"
#include "TMath.h"
#include "TColor.h"
#include "TStyle.h"
#include "TF1.h"
#include <exception>	// for stdexcept
#include <vector>
#include <map>
#include <string>
#include <algorithm>	// remove and remove_if
#include <iostream>
#include <iomanip>
#include <fstream> 		//std::ofstream
#include <stdlib.h>
#include <regex>
// #######################################################################
// we need to #include all the WCSim headers.
// we need to have the path to these headers exported to $ROOT_INCLUDE_PATH
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimPmtInfo.hh"
#include "WCSimLAPPDInfo.hh"
#include "WCSimEnumerations.hh"
#include "WCSimRootLinkDef.hh"
#if FILE_VERSION>1
#include "WCSimRootOptions.hh"
#endif
//#include "../wcsim/include/WCSimRootEvent.hh"
//#include "../wcsim/include/WCSimRootGeom.hh"
//#include "../wcsim/include/WCSimPmtInfo.hh"
//#include "../wcsim/include/WCSimLAPPDInfo.hh"
//#include "../wcsim/include/WCSimEnumerations.hh"
//#include "../wcsim/include/WCSimRootLinkDef.hh"
//#if FILE_VERSION>2
//#include "../wcsim/include/WCSimRootOptions.hh"
//#endif
// #######################################################################
// from ANNIEDAQ
#include "CardData.h"

// CLASS DEF
// =========
class WCSimAnalysis : public TObject {
	private:
	// CONSTS - put in a namespace not the class. combine with MRDspecs.hh. Pull from geo/opts if possible.
	const Int_t numtankpmts=128+2*(26);     // 26 pmts and lappds on each cap
	const Int_t nummrdpmts=307;             //
	const Int_t numvetopmts=26;             //
	const Int_t caparraysize=8;             // pmts on the cap form an nxn grid where caparraysize=n
	const Int_t pmtsperring=16;             // pmts around each ring of the main walls
	const Int_t numpmtrings=8;              // num rings around the main walls
	const Int_t MAXTRACKSPEREVENT=50;       //
	int pre_trigger_window_ns;              // read from options file
	int post_trigger_window_ns;             // read from options file
	
	int treeNumber=-1;
	Double_t maxsubeventduration=30.;  // in ns?
	// canvas sizes
	float win_scale;
	int n_wide;
	int n_high;
	
	// Input file read variables
	const char* inputdir;
	TChain* t; 
	WCSimRootGeom* geo=0;
	WCSimRootOptions *opt=0;
	WCSimRootEvent* b=0;
	WCSimRootEvent* m=0;
	WCSimRootEvent* v=0;
	TBranch* bp=0, *mp=0, *vp=0;
	WCSimRootTrigger* atrigt=0;
	WCSimRootTrigger* atrigm=0;
	WCSimRootTrigger* atrigv=0;
	WCSimRootEventHeader* header=0;
	int firstfilenum=0; // skips analysis of files in a TChain with number less than this
	TTree* currenttree;
	TFile* currentfile;
	std::string currentfilestring;
	Int_t eventnum;
	Int_t runnum;
	Int_t triggernum;
	
	// Variables for output files replicating ANNIE Raw data format
	std::string rawfilename;
	TFile* rawfileout=nullptr;
	// pmt data tree
	// ~~~~~~~~~~~~~
	TTree* tPMTData;
	ULong64_t fileout_LastSync, fileout_StartCount;
	Int_t fileout_SequenceID, fileout_StartTimeSec, fileout_StartTimeNSec, fileout_TriggerNumber,
		fileout_CardID, fileout_Channels, fileout_BufferSize, fileout_Eventsize, fileout_FullBufferSize;
	ULong64_t* fileout_TriggerCounts=nullptr;
	UInt_t* fileout_Rates=nullptr;
	UShort_t* fileout_Data=nullptr;
	// run info tree
	// ~~~~~~~~~~~~~
	TTree* tRunInformation;
	std::string fileout_InfoTitle, fileout_InfoMessage;
	// trigger data tree
	// ~~~~~~~~~~~~~~~~~
	TTree* tTrigData;
	Int_t fileout_FirmwareVersion, /*fileout_EventSize,*/ fileout_TriggerSize,
		fileout_FIFOOverflow, fileout_DriverOverfow;
	UShort_t* fileout_EventIDs=nullptr;
	ULong64_t* fileout_EventTimes=nullptr;
	UInt_t* fileout_TriggerMasks=nullptr, *fileout_TriggerCounters=nullptr;
	int MAXEVENTSIZE=10, MAXTRIGGERSIZE=10;
	// mrd data tree
	// ~~~~~~~~~~~~~
	TTree* tCCData;
	UInt_t fileout_Trigger, fileout_OutNumber;
	ULong64_t fileout_TimeStamp;
	std::vector<string> fileout_Type;
	std::vector<unsigned int> fileout_Value, fileout_Slot, fileout_Channel;
	// filling the output file
	// ~~~~~~~~~~~~~~~~~~~~~~~
	// methods
	void LoadOutputFiles();
	void FillEmulatedCCData();
	void FillEmulatedTrigData();
	void FillEmulatedPMTData();
	void FillEmulatedRunInformation();
	void AddCCDataEntry(WCSimRootCherenkovDigiHit* digihit);
	void AddPMTDataEntry(WCSimRootCherenkovDigiHit* digihit);
	void ConstructEmulatedPmtDataReadout();
	void AddMinibufferStartTime();
	void GenerateMinibufferPulse(int digit_index, double adjusted_digit_q, std::vector<uint16_t> &pulsevector);
	void RiffleShuffle(bool do_shuffle=true);
	std::vector<std::string>* GetTemplateRunInfo();
	// variables: TODO put in a namespace rather than the WCSimAnalysis class?
	Int_t sequence_id;
	Int_t minibuffer_id;
	// switches to allow turning this conversion on/off:
	bool add_emulated_ccdata=false, add_emulated_triggerdata=true, add_emulated_pmtdata=true;
	int full_buffer_size;                  // = CardData::FullBufferSize; read in utilityfuncs from options
	int minibuffer_datapoints_per_channel; // = CardData::BufferSize / CardData::TriggerNumber;
	int emulated_event_size;               // = CardData::Eventsize; 
	const int minibuffers_per_fullbuffer = CardData::TriggerNumber;
	const int channels_per_adc_card = CardData::Channels;
	const int channels_per_tdc_card = 32;
	const int num_adc_cards = (numtankpmts-1)/channels_per_adc_card + 1; // round up, requiring numtankpmts!=0
	// n.b. variant that doesn't require !=0, but could overflow:
	// num_adc_cards = (numtankpmts + channels_per_adc_card - 1) / channels_per_adc_card;
	unsigned long long placeholder_date_ns;
	std::vector<CardData> emulated_pmtdata_readout;
	std::vector<std::vector<uint16_t>> temporary_databuffers;
	std::vector<uint16_t> pulsevector;
	TF1* fLandau{nullptr};
	const int ADC_NS_PER_SAMPLE=2;
	const int MRD_NS_PER_SAMPLE=4;
	const int MRD_TIMEOUT_NS=4200;
	const unsigned long long MRD_TIMESTAMP_DELAY = static_cast<unsigned long long>(MRD_TIMEOUT_NS);
	const double ADC_INPUT_RESISTANCE = 50.;  // Ohm
	const double ADC_TO_VOLT = 2.415 / std::pow(2., 12);// * by this constant converts ADC counts to Volts
	const double PULSE_HEIGHT_FUDGE_FACTOR = (1./300.); // WHAT UNITS ARE DIGIT Q's IN?!?
	
	// PMT MAP VARIABLES
	// needed for drawing tank histograms
	std::map<int, std::pair<int,int> > topcappositionmap;
	std::map<int, std::pair<int,int> > bottomcappositionmap;
	std::map<int, std::pair<int,int> > wallpositionmap;
	
	
	// HISTOGRAMS
	// tank histograms
	// ~~~~~~~~~~~~~~~
	TH1F* tubeshithisttank=0;
	TH1F* truehitcounthisttank=0;
	TH1F* digihitcounthisttank=0;
	TH1D* PEdisttank=0;
	TH1D* PMThitfrequencytank=0;
	TH2D* QvsTtank=0;
	// TH1D* hitTimeDisttank=0;
	// TH1D* digitTimeDisttank=0;
	TH2D* topcaphist=0;
	TH2D* bottomcaphist=0;
	TH2D* wallhist=0;
	TH1F* AvgQvsTtankHist=0;
	// mrd histograms
	// ~~~~~~~~~~~~~~
	TH1F* tubeshithistmrd=0;
	TH1F* truehitcounthistmrd=0;
	TH1F* digihitcounthistmrd=0;
	TH1D* PEdistmrd=0;
	TH1D* PMThitfrequencymrd=0;
	TH2D* QvsTmrd=0;
	TH1D* hitTimeDistmrd=0;
	TH1D* digitTimeDistmrd=0;
	TH2D* mrdhist=0;
	// veto histograms
	// ~~~~~~~~~~~~~~~
	TH1F* tubeshithistveto=0;
	TH1F* truehitcounthistveto=0;
	TH1F* digihitcounthistveto=0;
	TH1D* PEdistveto=0;
	TH1D* PMThitfrequencyveto=0;
	TH2D* QvsTveto=0;
	TH1D* hitTimeDistveto=0;
	TH1D* digitTimeDistveto=0;
	TH2D* facchist=0;
	TH2D* PMTsvDigitTimeveto=0;
	
	// CANVASES
	// global canvases
	// ~~~~~~~~~~~~~~~
	TCanvas* hitCountCanv=0;
	// tank canvases
	// ~~~~~~~~~~~~~~~
	TCanvas* wallmapcanv=0;
	TCanvas* topcapmapcanv=0;
	TCanvas* bottomcapmapcanv=0;
	TCanvas *QvsTtankCanv=0;
	// mrd canvases
	// ~~~~~~~~~~~~~~~
	TCanvas *QvsTmrdCanv=0;
	// veto canvases
	// ~~~~~~~~~~~~~~~
	TCanvas *QvsTvetoCanv=0;
	TCanvas* vetoLastTimevsNumHitsCanv=0;
	
	// MRD TRACK RECONSTRUCTION
	// ~~~~~~~~~~~~~~~~~~~~~~~~
	// variables for file writing
	int wcsimfilenum;
	const char* outputdir="";
	TFile* mrdtrackfile=0, *vetotrackfile=0;
	TTree* mrdtree=0; // mrd track reconstruction tree
	TTree* vetotree=0; // veto track reconstruction tree
	std::vector<Double_t> mrddigittimesthisevent;
	Int_t nummrdsubeventsthisevent;
	Int_t nummrdtracksthisevent;
	TBranch* mrdeventnumb;
	TBranch* mrdtriggernumb;
	TBranch* nummrdsubeventsthiseventb=0;
	TBranch* nummrdtracksthiseventb=0;
	TBranch* subeventsinthiseventb=0;
	TClonesArray* aTrack=0;
	TClonesArray* aSubEvent=0;
	std::vector<double> vetodigittimesthisevent;
	Int_t numvetoeventsthisevent=0;
	TBranch* vetoeventnumb;
	TBranch* vetotriggernumb;
	TBranch* numvetoeventsthiseventb;
	TBranch* vetoeventsinthiseventb;
	TClonesArray* FaccSubEvents=0;
	
	// DISABLING STUFF
	// ~~~~~~~~~~~~~~
	Bool_t drawtankhistos=false;
	Bool_t drawmrdhistos=false;
	Bool_t drawvetohistos=false;
	
	public:
	// constructor + destructor
	WCSimAnalysis(const char* indir="/home/marc/anniegpvm/stats10k", const char* outdir=".");
	~WCSimAnalysis();
	
	// functions - initialization and utility
	void InitEnvironment();
	void LoadInputFiles();
	void MakePMTmap();
	void GetTreeData();
	int LoadTchainEntry(Int_t &eventnum);
	
	// functions - pre-event-loop analysis initializations
	void DoTankPreEventLoop();
	void DoMRDpreEventLoop();
	void DoVetoPreEventLoop();
	
	// functions - event and hit wide loops
	void DoTankPreHitLoop();
	void DoTankPreTriggerLoop();
	void DoTankTrigger(int &numtruehits, int &numdigits);
	void DoTankPostTriggerLoop(int &numtruehits, int &numdigits);
	void DoTankTrueHits();
	void DoTankDigitHits();
	void DoTankPostHitLoop();
	void DoMRDpreHitLoop();
	void DoMRDpreTriggerLoop();
	void DoMRDtrigger(int &numtruehits, int &numdigits);
	void DoMRDpostTriggerLoop(int &numtruehits, int &numdigits);
	void DoMRDtrueHits();
	void DoMRDdigitHits();
	void DoMRDpostHitLoop();
	void DoVetoPreHitLoop();
	void DoVetoPreTriggerLoop();
	void DoVetoTrigger(int &numtruehits, int &numdigits);
	void DoVetoPostTriggerLoop(int &numtruehits, int &numdigits);
	void DoVetoTrueHits();
	void DoVetoDigitHits();
	void DoVetoPostHitLoop();
	
	// functions - post-event-loop analysis initializations
	void DoTankPostEventLoop();
	void DoMRDpostEventLoop();
	void DoVetoPostEventLoop();
	
	// functions - histogram declaration
	void DefineTankHistos();
	void DefineMRDhistos();
	void DefineVetoHistos();
	
	// functions - histogram filling
	void FillTankEventWideHists(Int_t numtruehits, Int_t numdigits);
	void FillMRDeventWideHists(Int_t numtruehits, Int_t numdigits);
	void FillVetoEventWideHists(Int_t numtruehits, Int_t numdigits);
	
	void FillTankTrueHitsHist(WCSimRootCherenkovHit* hit);
	void FillTankDigiHitsHist(WCSimRootCherenkovDigiHit* digihit);
	void FillMRDtrueHitsHist(WCSimRootCherenkovHit* hit, WCSimRootCherenkovHitTime *hittime);
	void FillMRDdigiHitsHist(WCSimRootCherenkovDigiHit* digihit);
	void FillVetoTrueHitsHist(WCSimRootCherenkovHit* hit, WCSimRootCherenkovHitTime *hittime);
	void FillVetoDigiHitsHist(WCSimRootCherenkovDigiHit* digihit);
	
	// functions - histogram drawing
	void DrawTankHistos();
	void DrawMRDhistos();
	void DrawVetoHistos();
	void DrawGlobalHistos();
	
	// functions - plot styling
	void ColourPlotStyle();
	
	// functions - mrd track finding
	void OpenMRDtrackOutfile(int filenum);
	void FindMRDtracksInEvent();
	void OpenFACCtrackOutfile(int filenum);
	void FindVetoTracksInEvent();
	
	// the one that calls all the others
	void DoAnalysis();
	
	ClassDef(WCSimAnalysis,1);	// INCREMENT VERSION NUM EVERY TIME CLASS MEMBERS CHANGE
};

// #######################################################################

// CONSTRUCTOR
// ===========
WCSimAnalysis::WCSimAnalysis(const char* indir, const char* outdir) : inputdir(indir), outputdir(outdir) {

}

// ########################################################################

// DESTRUCTOR
// ===========
WCSimAnalysis::~WCSimAnalysis(){
	cout<<"Resetting chain branch addresses"<<endl;
	t->ResetBranchAddresses();
	//f->Close();
	
	cout<<"deleting geometry, branches and triggers"<<endl;
	// wcsim classes
	if( geo ) { delete geo; geo=0; }
	if( opt ) { delete opt; opt=0; }
	cout<<"deleting branches"<<endl;
	if( b ) { delete b; b=0; }
	if( m ) { delete m; m=0; }
	if( v ) { delete v; v=0; }
// doesn't like this?
//	cout<<"deleting triggers"<<endl;
//	if( atrigt ) { delete atrigt; atrigt=0; }
//	if( atrigm ) { delete atrigm; atrigm=0; }
//	if( atrigv ) { delete atrigv; atrigv=0; }
	
	// HISTOGRAMS
	cout<<"deleting histograms"<<endl;
	// tank histograms
	// ~~~~~~~~~~~~~~
	if( tubeshithisttank ) { delete tubeshithisttank; tubeshithisttank=0; }
	if( truehitcounthisttank ) { delete truehitcounthisttank; truehitcounthisttank=0; }
	if( digihitcounthisttank ) { delete digihitcounthisttank; digihitcounthisttank=0; }
	if( PEdisttank ) { delete PEdisttank; PEdisttank=0; }
	if( PMThitfrequencytank ) { delete PMThitfrequencytank; PMThitfrequencytank=0; }
	if( QvsTtank ) { delete QvsTtank; QvsTtank=0; }
	// if( hitTimeDisttank ) { delete hitTimeDisttank; hitTimeDisttank=0; }
	// if( digitTimeDisttank ) { delete digitTimeDisttank; digitTimeDisttank=0; }
	if( topcaphist ) { delete topcaphist; topcaphist=0; }
	if( bottomcaphist ) { delete bottomcaphist; bottomcaphist=0; }
	if( wallhist ) { delete wallhist; wallhist=0; }
	if( AvgQvsTtankHist ) { delete AvgQvsTtankHist; AvgQvsTtankHist=0; }
	// mrd histograms
	// ~~~~~~~~~~~~~~
	if( tubeshithistmrd ) { delete tubeshithistmrd; tubeshithistmrd=0; }
	if( truehitcounthistmrd ) { delete truehitcounthistmrd; truehitcounthistmrd=0; }
	if( digihitcounthistmrd ) { delete digihitcounthistmrd; digihitcounthistmrd=0; }
	if( PEdistmrd ) { delete PEdistmrd; PEdistmrd=0; }
	if( PMThitfrequencymrd ) { delete PMThitfrequencymrd; PMThitfrequencymrd=0; }
	if( QvsTmrd ) { delete QvsTmrd; QvsTmrd=0; }
	if( hitTimeDistmrd ) { delete hitTimeDistmrd; hitTimeDistmrd=0; }
	if( digitTimeDistmrd ) { delete digitTimeDistmrd; digitTimeDistmrd=0; }
	if( mrdhist ) { delete mrdhist; mrdhist=0; }
	// veto histograms
	// ~~~~~~~~~~~~~~~
	if( tubeshithistveto ) { delete tubeshithistveto; tubeshithistveto=0; }
	if( truehitcounthistveto ) { delete truehitcounthistveto; truehitcounthistveto=0; }
	if( digihitcounthistveto ) { delete digihitcounthistveto; digihitcounthistveto=0; }
	if( PEdistveto ) { delete PEdistveto; PEdistveto=0; }
	if( PMThitfrequencyveto ) { delete PMThitfrequencyveto; PMThitfrequencyveto=0; }
	if(QvsTveto ) { delete QvsTveto; QvsTveto=0; }
	if( hitTimeDistveto ) { delete hitTimeDistveto; hitTimeDistveto=0; }
	if( digitTimeDistveto ) { delete digitTimeDistveto; digitTimeDistveto=0; }
	if( facchist ) { delete facchist; facchist=0; }
	if( PMTsvDigitTimeveto ) { delete PMTsvDigitTimeveto; PMTsvDigitTimeveto=0; }
	
	// CANVASES
	cout<<"deleting canvases"<<endl;
	// global canvases
	if( hitCountCanv ) { delete hitCountCanv; hitCountCanv=0; }
	// tank canvases
	if( wallmapcanv ) { delete wallmapcanv; wallmapcanv=0; }
	if( topcapmapcanv ) { delete topcapmapcanv; topcapmapcanv=0; }
	if( bottomcapmapcanv ) { delete bottomcapmapcanv; bottomcapmapcanv=0; }
	if( QvsTtankCanv ) { delete QvsTtankCanv; QvsTtankCanv=0; }
	// mrd canvases
	if( QvsTmrdCanv ) { delete QvsTmrdCanv; QvsTmrdCanv=0; }
	// veto canvases
	if( QvsTvetoCanv ) { delete QvsTvetoCanv; QvsTvetoCanv=0; }
	if( vetoLastTimevsNumHitsCanv ) { delete vetoLastTimevsNumHitsCanv; vetoLastTimevsNumHitsCanv=0; }
	
	// EMULATED OUTPUT STUFF
	cout<<"clearing up after emulated file conversion"<<endl;
	if( fLandau ){ delete fLandau; fLandau=0; }
	if( rawfileout ){ rawfileout->Close(); delete rawfileout; rawfileout=0; }
	
	// Close output files
	// TChain* t doesn't need closing...
	cout<<"closing file"<<endl;
	if( mrdtrackfile ) { mrdtrackfile->Close(); delete mrdtrackfile; mrdtrackfile=0; }	// deletes member branches too
	if( vetotrackfile ) { vetotrackfile->Close(); delete vetotrackfile; vetotrackfile=0; }
	cout<<"done with destructor"<<endl;
}

// ########################################################################

// FUNCTION DEFINITIONS
// ===========
#include "utilityfuncs.cxx"
#include "makepmtmaps.cxx"
#include "tankhists.cxx"
#include "vetohists.cxx"
#include "mrdhists.cxx"
#include "globalhists.cxx"
#include "findmrdtracks.cxx"
#include "doanalysis.cxx"		// main analysis: calls all the above and loops over events
#include "tankanalysis.cxx"
#include "mrdanalysis.cxx"
#include "vetoanalysis.cxx"
#include "findvetotracks.cxx"
#include "createrawfile.cxx"
#include "fill_emulated_ccdata.cxx"
#include "fill_emulated_triggerdata.cxx"
#include "fill_emulated_pmtdata.cxx"
#include "gettemplateruninfo.cxx"

//TODO std::string intxnumtotype(gst genieeventasclass){


#ifdef __CINT__
#pragma link C++ class WCSimAnalysis+;
#endif

// ###########################################################################

// WCSim class summary:

/* 
WCSimRootTrack has methods: 
Int_t     GetIpnu()             pdg
Int_t     GetFlag()             -1 = probe nu, -2 = target, 1 = fsl, 2 = most energetic fs nucleon
Float_t   GetM()                mass
Float_t   GetP()                momentum magnitude
Float_t   GetE()                energy (inc rest mass^2)
Int_t     GetStartvol()         starting volume
Int_t     GetStopvol()          stopping volume
Float_t   GetDir(Int_t i=0)     momentum unit vector
Float_t   GetPdir(Int_t i=0)    momentum vector
Float_t   GetStop(Int_t i=0)    stopping vertex x,y,z for i=0-2
Float_t   GetStart(Int_t i=0)   starting vertex x,y,z for i=0-2
Int_t     GetParenttype()       parent pdg
Float_t   GetTime()             trj->GetGlobalTime(); stopping(?) time of particle
Int_t     GetId()               wcsim trackid
*/

#endif

