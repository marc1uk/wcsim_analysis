/* vim:set noexpandtab tabstop=4 wrap */
#ifndef _MRDSubEvent_Class_
#define _MRDSubEvent_Class_ 1

#ifndef _MRDSubEvent_VERBOSE_
//#define _MRDSubEvent_VERBOSE_ 1
#endif

#include <TObject.h>
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TText.h"
#include "TLine.h"
#include "TArrow.h"
#include "TColor.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include "Math/GenVector/LorentzVector.h"
#include <exception>	// for stdexcept
#include <vector>
#include <algorithm>
#include <map>
#include "WCSimRootEvent.hh"
#include "MRDspecs.hh"
#include "MRDSubEvent_ReconstructionClasses.hh"	// defines classes used in DoReconstruction() function
#include "MRDTrackClass.hh"

class cMRDSubEvent : public TObject {
	
	// Private members
	// ===============
	private:
	Int_t mrdsubevent_id;					// ID of this track within the trigger
	
	// Raw Info:
	std::string wcsimfile;					// which wcsim file this was in
	Int_t run_id;							// which run this file was in  FIXME is same as mrdsubevent_id
	Int_t event_id;							// which event this track was in
	Int_t trigger;							// which (sub)trigger this track was in FIXME loads -1s/1's, no 0s
	std::vector<Int_t> digi_ids;			// vector of digi ids: GetCherenkovDigiHits()->At(digi_ids.at(i))
	std::vector<Int_t> pmts_hit;			// vector of PMT IDs hit
	std::vector<Double_t> digi_qs;			// vector of digit charges
	std::vector<Double_t> digi_ts;			// vector of digit times
	std::vector<Int_t> digi_numphots;		// number of true photons for each digit
	std::vector<Double_t> digi_phot_ts;		// true hit times of photons in a digit FIXME some up to -2000???
	std::vector<Int_t> digi_phot_parents;	// wcsim track IDs of parents that provided photons for a digit
//	std::vector<WCSimRootCherenkovDigiHit> digits;
	std::vector<WCSimRootTrack> truetracks;	// true WCSim tracks within this event time window
	// XXX FIXME XXX we could also use digi_phot_parents to sort hits in a track by their parent for truth comp!
	
	// Calculated/Reconstructed Info
	std::vector<cMRDTrack> tracksthissubevent;	// tracks created this SubEvent
	std::vector<Int_t> layers_hit;				// vector of layers hit TODO currently empty
	std::vector<Double_t> eDepsInLayers;		// fixed len vector of energy deposition in each layer TODO
	
	// Involved in drawing
	std::pair<double, double> xupcorner1, xupcorner2, xdowncorner1, xdowncorner2, yupcorner1, yupcorner2, ydowncorner1, ydowncorner2;  // TODO remove me, probably can just be defined in makemrdimage.cxx
	std::vector<TArrow*> trackfitarrows;  //!  stores TLines which are associated with track boundaries
	std::vector<TArrow*> trackarrows;     //!  stores TArrows associated with CA reconstructed tracks
	std::vector<TArrow*> truetrackarrows; //!  stores TArrows associated with true tracks
	
	public:
	// SubEvent Level Getters
	// ======================
	// Locate the subevent in file>run>event>trigger hierarchy
	Int_t GetSubEventID(){return mrdsubevent_id;}
	std::string GetFile(){return wcsimfile;}
	Int_t GetRunID(){return run_id;}
	Int_t GetEventID(){return event_id;}
	Int_t GetTrigger(){return trigger;}
	
	// Top level information about the subevent
	Int_t GetNumDigits(){return digi_ids.size();}
	Int_t GetNumLayersHit(){return layers_hit.size();}
	Int_t GetNumPMTsHit(){return pmts_hit.size();}
	std::vector<Int_t> GetDigitIds(){return digi_ids;}
	std::vector<Double_t> GetDigitQs(){return digi_qs;}
	std::vector<Double_t> GetDigitTs(){return digi_ts;}
	std::vector<Int_t> GetDigiNumPhots(){return digi_numphots;}
	std::vector<Double_t> GetDigiPhotTs(){return digi_phot_ts;}
	std::vector<Int_t> GetDigiPhotParents(){return digi_phot_parents;}
	std::vector<Int_t> GetLayersHit(){return layers_hit;}
	std::vector<Int_t> GetPMTsHit(){return pmts_hit;}
	std::vector<WCSimRootTrack> GetTrueTracks(){return truetracks;}
	
	// Reconstructed Variables
	std::vector<cMRDTrack>* GetTracks(){ return &tracksthissubevent;}
	std::vector<Double_t> GetEdeps(){return eDepsInLayers;}
	std::vector<TArrow*> GetTrackArrows(){return trackarrows;}
	std::vector<TArrow*> GetTrueTrackArrows(){return truetrackarrows;}
	std::vector<TArrow*> GetTrackFitArrows(){return trackfitarrows;}
	
	void Print(); // print the subevent info.
	
//	// Digit Level Getters - all the obtainable information about a digit. Or just return the digit?
//	// ===============================
//	Int_t GetPMTNumber(Int_t digitnum){return PMTnum.at(digitnum);}
//	Double_t GetTime(Int_t digitnum){return DigitTs.at(digitnum);}
//	Double_t GetCharge(Int_t digitnum){return DigitQs.at(digitnum);}
//	Double_t GetEdeposited(Int_t digitnum){ return DigitQs.at(digitnum)*3;}	//TODO derive from above. units?
//	Int_t GetNumTrueHits(Int_t digitnum){return NumTrueHits.at(digitnum);}	//size of PhotonIds()
//	std::vector<Int_t> GetPhotonIds(Int_t digitnum){return PhotonIds.at(digitnum);}
//	std::vector<Int_t> GetPhotonParents(Int_t digitnum){return PhotonParents.at(digitnum);}
//	std::vector<Double_t> GetPhotonTrueTimes(Int_t digitnum){return PhotonTrueTimes.at(digitnum);}
//	WCSimRootCherenkovDigiHit* GetDigit(Int_t i){
//		try{return &(digits.at(i));}
//		catch(const std::out_of_range& oor){return 0;}
//	}
	
	//---------------------------
	// digit information not in WCSimRootCherenkovDigiHit
	// Geometric information about a given digit
//	std::vector<Int_t> MRDlayers;
//	std::vector<Int_t> MRDpaddles;
//	std::vector<std::pair<Double_t, Double_t> > xranges;	// from width of panel(s) hit
//	std::vector<std::pair<Double_t, Double_t> > yranges;	// 
//	std::vector<std::pair<Double_t, Double_t> > zranges;	// from depth of panel
//	std::vector<std::pair<Double_t, Double_t> > tranges; 	// from uncertainty in PMT timing resoluton
	//---------------------------
//	Int_t GetLayerNum(Int_t digitnum){return MRDlayers.at(digitnum);}
//	Int_t GetPaddleNum(Int_t digitnum){return MRDpaddles.at(digitnum);}	// number of paddle within this panel
//	std::pair<Double_t, Double_t> GetXrange(Int_t digitnum){return xranges.at(digitnum);}
//	std::pair<Double_t, Double_t> GetYrange(Int_t digitnum){return yranges.at(digitnum);}
//	std::pair<Double_t, Double_t> GetZrange(Int_t digitnum){return zranges.at(digitnum);}
//	std::pair<Double_t, Double_t> GetTrange(Int_t digitnum){return tranges.at(digitnum);}
	//---------------------------
	
	// "Setters"
	// =========
//	void AppendDigit(WCSimRootCherenkovDigiHit digitin){digits.push_back(digitin);}
//	void AppendDigit(thisdigitstime, thisdigitsq, thisdigitstubeid, photontimesinatrack, particleidsinatrack);
	
	// Functions to do reconstruction
	// ==============================
	private:
	// Main track reconstruction code. Groups paddles in a line into a MRDTrack
	void DoReconstruction(bool printtracks, bool drawcells, bool drawfit);
	// Used within DoReconstruction (CA version) 
	void LeastSquaresMinimizer(Int_t numdatapoints, Double_t datapointxs[], Double_t datapointys[], Double_t datapointweights[], Double_t errorys[], Double_t &fit_gradient, Double_t &fit_offset, Double_t &chi2);
	bool BridgeSearch(const std::vector<mrdcell*> &tracktotest, const std::vector<std::pair<int,int> > &matchedtracks, const std::vector<std::vector<mrdcell*> > &allpaddletracks, const std::string horv);
	bool SearchForClusterInTracks(const std::vector<std::pair<int,int> > &matchedtracks, const std::vector<std::vector<mrdcell*> > &allpaddletracks, const std::vector<mrdcell*> tracktotest, const std::string horv);
	
	void FillStaticMembers(){
//		// fill static members
//		std::vector<Int_t> temp(aspectrum, aspectrum+19);
//		//aspectrumv.assign(temp.rbegin(), temp.rend());
//		aspectrumv.assign(temp.begin(), temp.end());
		for(int i=0; i<19; i++){
			Int_t colorindex = TColor::GetColor(colorhexes.at(i).c_str());
			aspectrumv.at(i)=colorindex;
		}
		fillstaticmembers=false;
	}
	
	// Default Constructor
	// ====================
	public:
	// Default constructor that initialises all private members required for ROOT classes
	cMRDSubEvent() : mrdsubevent_id(-1), wcsimfile(""), run_id(-1), event_id(-1), trigger(-1), digi_ids(), pmts_hit(), digi_qs(), digi_ts(), digi_numphots(), digi_phot_ts(), digi_phot_parents(), layers_hit(), eDepsInLayers(), tracksthissubevent(), truetracks() {};
	
	// destructor
	~cMRDSubEvent(){
#ifdef _MRDSubEvent_VERBOSE_
		cout<<"cMRDSubEvent destructor"<<endl;
#endif
		RemoveArrows();
	}
	
	// Actual Constructor
	// ==================
	cMRDSubEvent(Int_t mrdsubevent_idin, std::string wcsimefilein, Int_t runidin, Int_t eventidin,
	Int_t triggerin, std::vector<Int_t> digitidsin, std::vector<Int_t> digittubesin, std::vector<Double_t>
	digitqsin, std::vector<Double_t> digittimesin, std::vector<Int_t> digitnumphotsin, std::vector<Double_t> 
	digitstruetimesin, std::vector<Int_t> digitsparentsin, std::vector<WCSimRootTrack*> truetracksin) :
	/* information retrieved when creating the track: initialize with input */
	mrdsubevent_id(mrdsubevent_idin), wcsimfile(wcsimefilein), run_id(runidin), event_id(eventidin),
	trigger(triggerin), digi_ids(digitidsin), pmts_hit(digittubesin), digi_qs(digitqsin),
	digi_ts(digittimesin), digi_numphots(digitnumphotsin), digi_phot_ts(digitstruetimesin),
	digi_phot_parents(digitsparentsin),
	/* information calculated: initialize to default */
	layers_hit(), tracksthissubevent(), trackarrows(), truetrackarrows(), trackfitarrows() {
		eDepsInLayers.assign(MRDSpecs::numpanels, 0.);	// can't assign the size in the class def. 
		// we receive a set of pointers to truth tracks: to store them in the MRDSubEventClass
		// we need to clone them into a vector of objects here:
		for( WCSimRootTrack* atrack : truetracksin){
			truetracks.push_back(WCSimRootTrack(*atrack));  // copy constructor
		}
		
#ifdef _MRDSubEvent_VERBOSE_
		cout<<endl<<"constructing a subevent with "<<digi_ids.size()<<" digits"<<endl;
#endif
		if(fillstaticmembers) FillStaticMembers();
		
		Bool_t printtracks=false;
		Bool_t drawcells=true;
		Bool_t drawfit=true;
		Bool_t drawtruetracks=true;
		Bool_t saveimage=false;
		
		if(drawcells||drawfit||drawtruetracks) DrawMrdCanvases();  // creates the canvas with the digits
		if(drawtruetracks) DrawTrueTracks();   // draws true tracks over the event
		DoReconstruction(printtracks, drawcells, drawfit); // adds the tracks to the canvas
		if(saveimage) imgcanvas->SaveAs(TString::Format("mrdtracks_%d.png",event_id));
		//cout<<"sleeping for 5 seconds to analyse output"<<endl;
		//if(tracksthissubevent.size()) std::this_thread::sleep_for (std::chrono::seconds(15));
		//cout<<"moving to next event"<<endl;
		gPad->WaitPrimitive();
		RemoveArrows();		// removes true and reco track arrows so the canvas can be re-used
		//assert(false);
	}
	
	// Drawing
	// =======
	void DrawMrdCanvases();
	void DrawTrueTracks();
	void DrawTracks();
	static Bool_t fillstaticmembers;
	static TCanvas* imgcanvas;
	static TText* titleleft;
	static TText* titleright;
	static std::vector<TBox*> paddlepointers;
	void ComputePaddleTransformation (const Int_t copyNo, TVector3 &origin, Bool_t &ishpaddle);
	//Int_t aspectrum[19] = {kYellow, kOrange, (kOrange-3), (kOrange+8), (kOrange+10), kRed, (kRed+1), (kPink+4), (kMagenta+2), (kMagenta+1), kMagenta, (kViolet-2), (kViolet-3), (kViolet+7), (kViolet+9), (kBlue+2), (kBlue+1), kAzure, (kAzure+7)};
	static std::vector<Int_t> aspectrumv;
	static std::vector<std::string> colorhexes;
	static std::vector<EColor> trackcolours;
	
	void RemoveArrows(){	// sometimes need to clear the arrows even before deleting the subevent.
		for(auto anarrow : trackarrows){
			delete anarrow;
		}
		trackarrows.clear();
		
		for(auto anarrow : truetrackarrows){
			delete anarrow;
		}
		truetrackarrows.clear();
		
		for(auto anarrow : trackfitarrows){
			delete anarrow;
		}
		trackfitarrows.clear();
	}
	
	// Required by ROOT
	// ================
	void Clear(){
		cout<<"calling clear on cMRDSubEvent "<<mrdsubevent_id<<endl;
		mrdsubevent_id=-1;
		wcsimfile="";
		run_id=-1;
		event_id=-1;
		trigger=-1;
		digi_ids.clear();
		pmts_hit.clear();
		digi_qs.clear();
		digi_ts.clear();
		digi_numphots.clear();
		digi_phot_ts.clear();
		digi_phot_parents.clear();
//		digits.clear();
		layers_hit.clear();
		eDepsInLayers.assign(MRDSpecs::numpanels,0.);
		truetracks.clear();
		// heap allocated objects for drawing the event on a canvas
		RemoveArrows();
	}
	
	// End class definition
	// ====================
	ClassDef(cMRDSubEvent,1);					// INCREMENT VERSION NUM EVERY TIME CLASS MEMBERS CHANGE
};

Bool_t cMRDSubEvent::fillstaticmembers=true;
TCanvas* cMRDSubEvent::imgcanvas=0;
TText* cMRDSubEvent::titleleft=0;
TText* cMRDSubEvent::titleright=0;
// allocate paddle vector now: they'll be filled in first call to DrawMrdCanvases
std::vector<TBox*> cMRDSubEvent::paddlepointers(MRDSpecs::nummrdpmts+(2*MRDSpecs::numpanels));
//std::vector<Int_t> cMRDSubEvent::aspectrumv(19);
//std::vector<Int_t> cMRDSubEvent::aspectrumv = ( []()->std::vector<Int_t> { std::vector<Int_t> temp {kYellow, kOrange, (kOrange-3), (kOrange+8), (kOrange+10), kRed, (kRed+1), (kPink+4), (kMagenta+2), (kMagenta+1), kMagenta, (kViolet-2), (kViolet-3), (kViolet+7), (kViolet+9), (kBlue+2), (kBlue+1), kAzure, (kAzure+7)}; return temp; }() );
//std::vector<Int_t> cMRDSubEvent::aspectrumv{kYellow, kOrange, (kOrange-3), (kOrange+8), (kOrange+10), kRed, (kRed+1), (kPink+4), (kMagenta+2), (kMagenta+1), kMagenta, (kViolet-2), (kViolet-3), (kViolet+7), (kViolet+9), (kBlue+2), (kBlue+1), kAzure, (kAzure+7)};
std::vector<Int_t> cMRDSubEvent::aspectrumv{1435, 1436, 1437, 1438, 1439, 1440, 1441, 1442, 1443, 1443, 1444, 1445, 885, 1446, 1447, 1448, 879, 1449, 1450};
std::vector<std::string> cMRDSubEvent::colorhexes{"#21ffff", "#20deea", "#1fbcd5", "#21a8cd", "#269bcb", "#2b8fca", "#367fcb", "#416fcb", "#4965cd", "#505dcf", "#5855d0", "#6247cf", "#6d3ace", "#782acd", "#851acc", "#910dcc", "#9e07cd", "#aa00ce", "#bf00d7"};
std::vector<EColor> cMRDSubEvent::trackcolours{kBlack, kBlue, (EColor)TColor::GetColorDark(kGreen), kRed, kViolet, kOrange, kMagenta,(EColor)(kAzure+2),(EColor)(kOrange+4),(EColor)(kViolet-6),(EColor)(kTeal-6)};

#include "MRDSubEvent_DoReconstruction.cxx"			// contains reconstruction function definitions
#include "makemrdimage.cxx"							// functions to draw the MRD top and side views
#include "MRDSubEvent_DrawTruthTracks.cxx"			// contains the function to draw truth tracks
#include "MRDSubEvent_Draw_Print.cxx"				// contains the print function

#endif

#ifdef __CINT__
#pragma link C++ class cMRDSubEvent+;
//#pragma link C++ class ROOT::Math::XYZTVector+;
//#pragma link C++ class std::vector<ROOT::Math::XYZTVector>+;
#endif
