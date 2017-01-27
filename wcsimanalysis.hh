/* vim:set noexpandtab tabstop=4 wrap */
#include "TROOT.h"
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
#include "TString.h"
#include "TMath.h"
#include <exception>	// for stdexcept
#include <vector>
#include <map>
#include <string>
#include <algorithm>	// remove and remove_if
#include <iostream>
#include <iomanip>
#include <fstream> 		//std::ofstream
#include <stdlib.h>
// #######################################################################
// we need to #include all the WCSim headers.
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimPmtInfo.hh"
#include "WCSimLAPPDInfo.hh"
#include "WCSimEnumerations.hh"
#include "WCSimRootLinkDef.hh"
// #######################################################################

#ifndef WCSimAnalysisClass
#define WCSimAnalysisClass

// CLASS DEF
// =========
class WCSimAnalysis : public TObject {
	private:
	// CONSTS
	const Int_t numtankpmts=128+2*(26);	// 26 pmts and lappds on each cap
	const Int_t nummrdpmts=336;
	const Int_t numvetopmts=26;
	const Int_t caparraysize=8;    // pmts on the cap form an nxn grid where caparraysize=n
	const Int_t pmtsperring=16;    // pmts around each ring of the main walls
	const Int_t numpmtrings=8;     // num rings around the main walls
	const Int_t MAXTRACKSPEREVENT=50;

	int treeNumber;
	Double_t maxtrackduration;  // in ns?
	// canvas sizes
	float win_scale;
	int n_wide;
	int n_high;

	// Input file read variables
	const char* inputdir;
	TChain* t; 
	WCSimRootGeom* geo=0;
	WCSimRootEvent* b=0;
	WCSimRootEvent* m=0;
	WCSimRootEvent* v=0;
	TBranch* bp=0, *mp=0, *vp=0;
	WCSimRootTrigger* atrigt=0;
	WCSimRootTrigger* atrigm=0;
	WCSimRootTrigger* atrigv=0;

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
	TFile* mrdtrackfile=0;
	TTree* recotree=0;
	Int_t nummrddigitsthisevent;
	Int_t nummrdtracksthisevent;
	std::vector<std::vector<Int_t>> tubeidsinthistrack;
	std::vector<std::vector<Double_t>> digitqsinthistrack;
	std::vector<std::vector<Double_t>> digittimesinthistrack;
	std::vector<std::vector<Int_t>> particleidsinthistrack;
	TBranch* nummrddigitsthiseventb=0;
	TBranch* nummrdtracksthiseventb=0;
	TBranch* tubeidsinthistrackb=0;
	TBranch* digitqsinthistrackb=0;
	TBranch* digittimesinthistrackb=0;
	TBranch* particleidsinthistrackb=0;
	// vectors filled in histogram functions, used for track finding
	std::vector<int> mrddigittubesthisevent;
	std::vector<double> mrddigittimesthisevent;
	
	// DISABLING STUFF
	// ~~~~~~~~~~~~~~
	Bool_t drawtankhistos=false;
	Bool_t drawmrdhistos=false;
	Bool_t drawvetohistos=false;

	public:
	// constructor + destructor
	WCSimAnalysis(const char* indir="/home/marc/anniegpvm/stats10k");
	~WCSimAnalysis();
	
	// functions - initialization and utility
	void InitEnvironment();
	void LoadInputFiles();
	void MakePMTmap();
	void GetTreeData();
	int LoadTchainEntry(Int_t eventnum);
	
	// functions - pre-event-loop analysis initializations
	void DoTankPreLoop();
	void DoMRDpreLoop();
	void DoVetoPreLoop();
	
	// functions - event and hit wide loops
	void DoTankPreHitLoop();
	void DoTankEventwide(Int_t &numtruehits, Int_t &numdigits);
	void DoTankTrueHits();
	void DoTankDigitHits();
	void DoTankPostHitLoop();
	void DoMRDpreHitLoop();
	void DoMRDeventwide(Int_t &numtruehits, Int_t &numdigits);
	void DoMRDtrueHits();
	void DoMRDdigitHits();
	void DoMRDpostHitLoop();
	void DoVetoPreHitLoop();
	void DoVetoEventwide(Int_t &numtruehits, Int_t &numdigits);
	void DoVetoTrueHits();
	void DoVetoDigitHits();
	void DoVetoPostHitLoop();
	
	// functions - post-event-loop analysis initializations
	void DoTankPostLoop();
	void DoMRDpostLoop();
	void DoVetoPostLoop();
	
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
	void OpenMRDtrackOutfile();
	void SplitMrdTracks();			//TODO: store tubeids as well
	void FindMRDtracksInEvent();	//TODO: write this
	
	// the one that calls all the others
	void DoAnalysis();

	ClassDef(WCSimAnalysis,1);	// INCREMENT VERSION NUM EVERY TIME CLASS MEMBERS CHANGE
};

// #######################################################################

// CONSTRUCTOR
// ===========
WCSimAnalysis::WCSimAnalysis(const char* indir) : inputdir(indir) {

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

	// Close output files
	// TChain* t doesn't need closing...
	cout<<"closing file"<<endl;
	if( mrdtrackfile ) { mrdtrackfile->Close(); delete mrdtrackfile; mrdtrackfile=0; }	// deletes member branches too
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
#include "MRDStrikeClass.hh"	// basically an extension of a digit, with MRD specific info
#include "MRDTrackClass.hh"		// a class for defining MRD tracks

//TODO std::string intxnumtotype(gst genieeventasclass){

#endif

#ifdef __CINT__
#pragma link C++ class WCSimAnalysis+;
#endif

