/* vim:set noexpandtab tabstop=4 wrap */
#ifndef _VetoEvent_VERBOSE_
//#define _VetoEvent_VERBOSE_ 1
#endif

#ifndef _VetoEvent_Class_
#define _VetoEvent_Class_ 1

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

class cVetoEvent : public TObject {
	
	// Private members
	// ===============
	private:
	Int_t VetoEventID;						// ID of this track within the subtrigger
	
	// Raw Info:
	std::string wcsimfile;					// which wcsim file this was in
	Int_t run_id;							// which run this file was in  FIXME is same as VetoEventID
	Int_t event_id;							// which event this track was in FIXME loads of 0's
	Int_t subtrigger;						// which (sub)trigger this track was in FIXME loads -1s/1's, no 0s
	std::vector<Int_t> digi_ids;			// vector of digi ids: GetCherenkovDigiHits()->At(digi_ids.at(i))
	std::vector<Int_t> pmts_hit;			// vector of PMT IDs hit
	std::vector<Double_t> digi_qs;			// vector of digit charges
	std::vector<Double_t> digi_ts;			// vector of digit times
	std::vector<Int_t> digi_numphots;		// number of true photons for each digit
	std::vector<Double_t> digi_phot_ts;		// true hit times of photons in a digit FIXME some up to -2000???
	std::vector<Int_t> digi_phot_parents;	// wcsim track IDs of parents that provided photons for a digit
//	std::vector<WCSimRootCherenkovDigiHit> digits;
	std::vector<WCSimRootTrack> truetracks;	// true WCSim tracks within this event time window
	
	// Calculated/Reconstructed Info
	std::vector<Int_t> layers_hit;			// vector of layers hit TODO currently empty
	std::vector<Double_t> eDepsInLayers;	// (fixed length) vector of net energy deposition in each layer ^TODO
	
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
	std::vector<WCSimRootTrack> GetTrueTracks(){return truetracks;}
	
	// Reconstructed Variables
	std::vector<Double_t> GetEdeps(){return eDepsInLayers;}
	
	// Functions to do reconstruction
	// ==============================
	private:
	
	// Default Constructor
	// ====================
	public:
	// Default constructor that initialises all private members required for ROOT classes
	cVetoEvent() : VetoEventID(-1), wcsimfile(""), run_id(-1), event_id(-1), subtrigger(-1), digi_ids(), pmts_hit(), digi_qs(), digi_ts(), digi_numphots(), digi_phot_ts(), digi_phot_parents(), layers_hit(), eDepsInLayers(), truetracks() {cout<<"cVetoEvent default constructor"<<endl;};
	
	// destructor
	~cVetoEvent(){
		//cout<<"cVetoEvent destructor"<<endl;
	}
	
	// Actual Constructor
	// ==================
	cVetoEvent(Int_t VetoEventidin, std::string wcsimefilein, Int_t runidin, Int_t eventidin,
	Int_t subtriggerin, std::vector<Int_t> digitidsin, std::vector<Int_t> digittubesin, std::vector<Double_t>
	digitqsin, std::vector<Double_t> digittimesin, std::vector<Int_t> digitnumphotsin, std::vector<Double_t> 
	digitstruetimesin, std::vector<Int_t> digitsparentsin, std::vector<WCSimRootTrack*> truetracksin) :
	/* information retrieved when creating the track: initialize with input */
	VetoEventID(VetoEventidin), wcsimfile(wcsimefilein), run_id(runidin), event_id(eventidin),
	subtrigger(subtriggerin), digi_ids(digitidsin), pmts_hit(digittubesin), digi_qs(digitqsin),
	digi_ts(digittimesin), digi_numphots(digitnumphotsin), digi_phot_ts(digitstruetimesin),
	digi_phot_parents(digitsparentsin),
	/* information calculated: initialize to default */
	layers_hit() {
		eDepsInLayers.assign(MRDSpecs::numpanels, 0.);	// can't assign the size in the class def. 
		// we receive a set of pointers to truth tracks: to store them in the VetoEventClass
		// we need to clone them into a vector of objects here:
		for( WCSimRootTrack* atrack : truetracksin){
			truetracks.push_back(WCSimRootTrack(*atrack));  // copy constructor
		}
		
#ifdef _VetoEvent_VERBOSE_
		cout<<endl<<"constructing a veto event with "<<digi_ids.size()<<" digits"<<endl;
#endif
		if(fillstaticmembers){
//			// fill static members
//			// placeholder
			fillstaticmembers=false;
		}
		
		bool drawhits=false;
		//if(drawhits) DrawHits();   // TODO colour veto paddles on the global event image (doesn't exist yet)
		//cout<<"sleeping for 5 seconds to analyse output"<<endl;
		//if(tracksthissubevent.size()) std::this_thread::sleep_for (std::chrono::seconds(15));
		//cout<<"moving to next event"<<endl;
		//gPad->WaitPrimitive();
	}
	
	// Drawing
	// =======
	// void DrawHits(); // TODO: add coloring to veto paddles on MRD to indicate hits.
	// 					// TODO: i suppose this requires a complete event-wide image canvas... 
	static Bool_t fillstaticmembers;
	
	// Required by ROOT
	// ================
	void Clear(){
		cout<<"calling clear on cVetoEvent "<<VetoEventID<<endl;
		VetoEventID=-1;
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
		truetracks.clear();
		layers_hit.clear();
		eDepsInLayers.assign(MRDSpecs::numpanels,0.);
	}
	
	// End class definition
	// ====================
	ClassDef(cVetoEvent,1);					// INCREMENT VERSION NUM EVERY TIME CLASS MEMBERS CHANGE
};
Bool_t cVetoEvent::fillstaticmembers=true;


#endif

#ifdef __CINT__
#pragma link C++ class cVetoEvent+;
//#pragma link C++ class ROOT::Math::XYZTVector+;
//#pragma link C++ class std::vector<ROOT::Math::XYZTVector>+;
#endif
