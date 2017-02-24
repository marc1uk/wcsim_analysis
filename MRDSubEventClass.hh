/* vim:set noexpandtab tabstop=4 wrap */
#ifndef _MRDSubEvent_VERBOSE_
//#define _MRDSubEvent_VERBOSE_ 1
#endif

#ifndef _MRDSubEvent_Class_
#define _MRDSubEvent_Class_ 1

#include <TObject.h>
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TVector3.h"
#include "TText.h"
#include "TLine.h"
#include "TArrow.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include "Math/GenVector/LorentzVector.h"
#include <exception>	// for stdexcept
#include <vector>
#include <algorithm>
#include <map>
#include "MRDspecs.hh"
#include "MRDSubEvent_ReconstructionClasses.hh"	// defines classes used in DoReconstruction() function
#include "MRDTrackClass.hh"

class cMRDSubEvent : public TObject {
	
	// Private members
	// ===============
	private:
	Int_t MRDSubEventID;						// ID of this track within the subtrigger
	
	// Raw Info:
	std::string wcsimfile;					// which wcsim file this was in
	Int_t run_id;							// which run this file was in
	Int_t event_id;							// which event this track was in
	Int_t subtrigger;						// which (sub)trigger this track was in
	std::vector<Int_t> digi_ids;			// vector of digi ids: GetCherenkovDigiHits()->At(digi_ids.at(i))
	std::vector<Int_t> pmts_hit;			// vector of PMT IDs hit
	std::vector<Double_t> digi_qs;			// vector of digit charges
	std::vector<Double_t> digi_ts;			// vector of digit times
	std::vector<Int_t> digi_numphots;		// number of true photons for each digit
	std::vector<Double_t> digi_phot_ts;		// true hit times of photons in a digit
	std::vector<Int_t> digi_phot_parents;	// wcsim track IDs of parents that provided photons for a digit
//	std::vector<WCSimRootCherenkovDigiHit> digits;
	
	// Calculated/Reconstructed Info
	std::vector<cMRDTrack*> tracksthissubevent;	//! tracks created this SubEvent (don't delete the '!')
	std::vector<Int_t> layers_hit;			// vector of layers hit
	std::vector<Double_t> eDepsInLayers;	// (fixed length) vector of net energy deposition in each layer
	
	// Involved in drawing
	std::pair<double, double> xupcorner1, xupcorner2, xdowncorner1, xdowncorner2, yupcorner1, yupcorner2, ydowncorner1, ydowncorner2;
	std::vector<TLine*> trackboundaries;
	std::vector<TArrow*> trackarrows;
	
	public:
	// SubEvent Level Getters
	// ======================
	// Locate the subevent in file>run>event>trigger hierarchy
	std::string GetFile(){return wcsimfile;}
	Int_t GetRunID(){return run_id;}
	Int_t GetEventID(){return event_id;}
	Int_t GetSubTrigger(){return subtrigger;}
	
	// Top level information about the subevent
	Int_t GetNumDigits(){return digi_ids.size();}
	Int_t GetNumLayersHit(){return layers_hit.size();}
	Int_t GetNumPMTsHit(){return pmts_hit.size();}
	std::vector<Int_t> GetDigitIds(){return digi_ids;}
	std::vector<Int_t> GetLayersHit(){return layers_hit;}
	std::vector<Int_t> GetPMTsHit(){return pmts_hit;}
	
	// Reconstructed Variables
	std::vector<cMRDTrack*> GetTracks(){return tracksthissubevent;}
	std::vector<Double_t> GetEdeps(){return eDepsInLayers;}
	
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
	void DoReconstruction();
	// Used within DoReconstruction (CA version) 
	void LeastSquaresMinimizer(Int_t numdatapoints, Double_t datapointxs[], Double_t datapointys[], Double_t datapointweights[], Double_t errorys[], Double_t &fit_gradient, Double_t &fit_offset, Double_t &chi2);
	
	// Default Constructor
	// ====================
	public:
	// Default constructor that initialises all private members required for ROOT classes
	cMRDSubEvent() : MRDSubEventID(-1), wcsimfile(""), run_id(-1), event_id(-1), subtrigger(-1), digi_ids(), pmts_hit(), digi_qs(), digi_ts(), digi_numphots(), digi_phot_ts(), digi_phot_parents(), layers_hit(), eDepsInLayers(), tracksthissubevent() {};
	
	// destructor
	~cMRDSubEvent(){
		// tracksthissubevent is a vector of pointers to cMRDTrack objects created on the heap
		for(auto atrack : tracksthissubevent){
			delete atrack;
		}
		tracksthissubevent.clear();
		
		// trackboundaries vector stores TLines on the heap which are associated with drawing tracks
		// delete them when the SubEvent gets destroyed. 
		for(auto aline : trackboundaries){
			delete aline;
		}
		trackboundaries.clear();
		for(auto anarrow : trackarrows){
			delete anarrow;
		}
		trackarrows.clear();
	}
	
	// Actual Constructor
	// ==================
	cMRDSubEvent(Int_t mrdsubeventidin, std::string wcsimefilein, Int_t runidin, Int_t eventidin,
	Int_t subtriggerin, std::vector<Int_t> digitidsin, std::vector<Int_t> digittubesin, std::vector<Double_t>
	digitqsin, std::vector<Double_t> digittimesin, std::vector<Int_t> digitnumphotsin, std::vector<Double_t> 
	digitstruetimesin, std::vector<Int_t> digitsparentsin) :
	/* information retrieved when creating the track: initialize with input */
	MRDSubEventID(mrdsubeventidin), wcsimfile(wcsimefilein), run_id(runidin), event_id(eventidin),
	subtrigger(subtriggerin), digi_ids(digitidsin), pmts_hit(digittubesin), digi_qs(digitqsin),
	digi_ts(digittimesin), digi_numphots(digitnumphotsin), digi_phot_ts(digitstruetimesin),
	digi_phot_parents(digitsparentsin),
	/* information calculated: initialize to default */
	layers_hit(), tracksthissubevent() {
		eDepsInLayers.assign(numpanels, 0.);	// can't assign the size in the class def. 
		
#ifdef _MRDSubEvent_VERBOSE_
		cout<<endl<<"constructing a subevent with "<<digi_ids.size()<<" digits"<<endl;
#endif
//		if(fillstaticmembers){
//			// fill static members
//			std::vector<Int_t> temp(aspectrum, aspectrum+19);
//			//aspectrumv.assign(temp.rbegin(), temp.rend());
//			aspectrumv.assign(temp.begin(), temp.end());
//			fillstaticmembers=false;
//		}
		DrawMrdCanvases();	// creates the canvas with the digits
		DoReconstruction();	// adds the tracks to the canvas
		imgcanvas->SaveAs(TString::Format("mrdtracks_%d.png",event_id));
		//if(tracksthissubevent.size()) std::this_thread::sleep_for (std::chrono::seconds(5));
		RemoveArrows();		// removes arrows so the canvas can be re-used
		//assert(false);
	}
	
	// Drawing
	// =======
	void DrawMrdCanvases();
	void AddTrackLines();
	static Bool_t fillstaticmembers;
	static TCanvas* imgcanvas;
	static TText* titleleft;
	static TText* titleright;
	static std::vector<TBox*> paddlepointers;
	void ComputePaddleTransformation (const Int_t copyNo, TVector3 &origin, Bool_t &ishpaddle, Bool_t paddleishit);
	//Int_t aspectrum[19] = {kYellow, kOrange, (kOrange-3), (kOrange+8), (kOrange+10), kRed, (kRed+1), (kPink+4), (kMagenta+2), (kMagenta+1), kMagenta, (kViolet-2), (kViolet-3), (kViolet+7), (kViolet+9), (kBlue+2), (kBlue+1), kAzure, (kAzure+7)};
	static std::vector<Int_t> aspectrumv;
	
	void RemoveArrows(){	// sometimes need to clear the arrows even before deleting the subevent.
		for(auto anarrow : trackarrows){
			delete anarrow;
		}
		trackarrows.clear();
	}
	
	// Required by ROOT
	// ================
	void Clear(){
		MRDSubEventID=-1;
		wcsimfile="";
		run_id=-1;
		event_id=-1;
		subtrigger=-1;
		digi_ids.clear();
		pmts_hit.clear();
		digi_qs.clear();
		digi_ts.clear();
		digi_numphots.clear();
		digi_phot_ts.clear();
		digi_phot_parents.clear();
//		digits.clear();
		// tracksthissubevent is a vector of cMRDTrack objects on the heap.
		for(auto atrack : tracksthissubevent){
			delete atrack;
		}
		tracksthissubevent.clear();
		layers_hit.clear();
		eDepsInLayers.assign(numpanels,0.);
		// heap allocated objects for drawing the event on a canvas
		for(auto aline : trackboundaries){
			delete aline;
		}
		trackboundaries.clear();
		for(auto anarrow : trackarrows){
			delete anarrow;
		}
		trackarrows.clear();
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
std::vector<TBox*> cMRDSubEvent::paddlepointers(nummrdpmts+(2*numpanels));
//std::vector<Int_t> cMRDSubEvent::aspectrumv(19);
std::vector<Int_t> cMRDSubEvent::aspectrumv = ( []()->std::vector<Int_t> { std::vector<Int_t> temp {kYellow, kOrange, (kOrange-3), (kOrange+8), (kOrange+10), kRed, (kRed+1), (kPink+4), (kMagenta+2), (kMagenta+1), kMagenta, (kViolet-2), (kViolet-3), (kViolet+7), (kViolet+9), (kBlue+2), (kBlue+1), kAzure, (kAzure+7)}; return temp; }() );

#include "MRDSubEvent_DoReconstruction.cxx"			// contains reconstruction function definitions
#include "makemrdimage.cxx"							// functions to draw the MRD top and side views

#endif

#ifdef __CINT__
#pragma link C++ class cMRDSubEvent+;
//#pragma link C++ class ROOT::Math::XYZTVector+;
//#pragma link C++ class std::vector<ROOT::Math::XYZTVector>+;
//#pragma link C++ class cMRDStrike+;
#pragma link C++ class std::vector<cMRDStrike>+;
#endif
