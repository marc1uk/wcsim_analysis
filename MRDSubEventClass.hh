/* vim:set noexpandtab tabstop=4 wrap */
#ifndef _MRDSubEvent_VERBOSE_
#define _MRDSubEvent_VERBOSE_ 1
#endif

#ifndef _MRDSubEvent_Class_
#define _MRDSubEvent_Class_ 1

#include <TObject.h>
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TVector3.h"
#include "TText.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include "Math/GenVector/LorentzVector.h"
#include <exception>	// for stdexcept
#include <vector>
#include <algorithm>
#include <map>
#include "MRDspecs.hh"
#include "MRDSubEventClass_ReconstructionClasses.hh"	// defines classes used in DoReconstruction() function

class cMRDSubEvent : public TObject {
	
	// Private members
	// ===============
	private:
	Int_t MRDSubEventID;						// ID of this track within the subtrigger
	Int_t tanktrackID;						// correlated tank track within the subtrigger
	
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
	std::vector<mrdcluster*> theclusters;	// vector of the clusters in the track
	std::vector<mrdcell*> thecells;			// vector of the cells in the track
	std::vector<Int_t> layers_hit;			// vector of MRD layers hit
	std::vector<Double_t> eDepsInLayers;	// fixed length vector of energy deposited in each layer 
	Double_t KEStart;						// from depth or estimate if fully penetrating
	Double_t KEEnd;							// 0 for stopping, or estimate if fully penetrating
	Int_t particlePID;						// estimated PID
	Int_t tracktype;						// are all strikes consistent with one track? 
											// -1: unchecked, 0: inconsistent, 1: consistent
	
	// Truth information:
	Int_t trueTrackID;						// index in the WCSimRootTrack clones array
//	WCSimRootTrack* trueTrack;				// from WCSim GetTracks. Should we keep this?
	
//	// Do we need these?
//	std::vector<ROOT::Math::XYZTVector> mrdpoints;	// estimated times and positions of interaction points
//	std::vector<ROOT::Math::XYZTVector> tankpoints;	// projected points from tank track, if available
	
	// These are needed for drawing the track in root
	std::map<const char*, double> recovalshoriz;	// keys: "xycrossing", "zcrossing", "angmin" and "angmax"
	std::map<const char*, double> recovalsvert;		// keys: "xycrossing", "zcrossing", "angmin" and "angmax"
	
	public:
	// Track Level Getters
	// ===================
	// Locate the track in file>run>event>trigger hierarchy
	std::string GetFile(){return wcsimfile;}
	Int_t GetRunID(){return run_id;}
	Int_t GetEventID(){return event_id;}
	Int_t GetSubTrigger(){return subtrigger;}
	
	// Top level information about the track
	Int_t GetNumDigits(){return digi_ids.size();}
	Int_t GetNumLayersHit(){return layers_hit.size();}
	Int_t GetNumPMTsHit(){return pmts_hit.size();}
	std::vector<Int_t> GetDigitIds(){return digi_ids;}
	std::vector<Int_t> GetLayersHit(){return layers_hit;}
	std::vector<Int_t> GetPMTsHit(){return pmts_hit;}
	
	// Reconstructed Variables
	std::vector<Double_t> GetEdeps(){return eDepsInLayers;}
	Double_t GetKEStart(){return KEStart;}
	Double_t GetKEEnd(){return KEEnd;}
	Double_t GetTotalEdep(){return KEStart-KEEnd;}
	Int_t GetparticlePID(){return particlePID;}
	Int_t GetTrackType(){return tracktype;}
	std::map<const char*, double> GetTrackXStats(){ return recovalshoriz; }
	std::map<const char*, double> GetTrackYStats(){ return recovalsvert; }
	
	// Truth Level Info
	Int_t GetTrueTrackID(){return trueTrackID;}
//	WCSimRootTrack* GetTrueTrack(){return trueTrack;}
	
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
	// Geometric information about a given digit
//	Int_t GetLayerNum(Int_t digitnum){return MRDlayers.at(digitnum);}
//	Int_t GetPaddleNum(Int_t digitnum){return MRDpaddles.at(digitnum);}	// number of paddle within this panel
//	std::pair<Double_t, Double_t> GetXrange(Int_t digitnum){return xranges.at(digitnum);}
//	std::pair<Double_t, Double_t> GetYrange(Int_t digitnum){return yranges.at(digitnum);}
//	std::pair<Double_t, Double_t> GetZrange(Int_t digitnum){return zranges.at(digitnum);}
//	std::pair<Double_t, Double_t> GetTrange(Int_t digitnum){return tranges.at(digitnum);}
	//---------------------------
	// digit information not in WCSimRootCherenkovDigiHit
//	std::vector<Int_t> MRDlayers;
//	std::vector<Int_t> MRDpaddles;
//	std::vector<std::pair<Double_t, Double_t> > xranges;	// from width of panel(s) hit
//	std::vector<std::pair<Double_t, Double_t> > yranges;	// 
//	std::vector<std::pair<Double_t, Double_t> > zranges;	// from depth of panel
//	std::vector<std::pair<Double_t, Double_t> > tranges; 	// from uncertainty in PMT timing resoluton
	//---------------------------
	
	// Geometric Variables
//	std::vector<ROOT::Math::XYZTVector> GetMRDPoints(){return mrdpoints;} 			// do we need this?
//	ROOT::Math::XYZTVector* GetPoint(Int_t i){										// 
//		try{return &(mrdpoints.at(i));}
//		catch(const std::out_of_range& oor){return 0;}
//	}
//	ROOT::Math::XYZTVector* GetFirstPoint(){return &(mrdpoints.front());}			// 
//	ROOT::Math::XYZTVector* GetLastPoint(){return &(mrdpoints.back());}				// 
//	Double_t GetPenetration(){return layers_hit.back().Z()-layers_hit.front().Z();}	// 
	
	
	// "Setters"
	// =========
//	void AppendDigit(WCSimRootCherenkovDigiHit digitin){digits.push_back(digitin);}
//	void AppendDigit(thisdigitstime, thisdigitsq, thisdigitstubeid, photontimesinatrack, particleidsinatrack);
	void SetTankTrack(Int_t trackidin){tanktrackID = trackidin;}
	
	// Functions to do reconstruction
	// ==============================
	private:
	// Main track reconstruction code
	void DoReconstruction();
	// Used within DoReconstruction (CA version) 
	void LeastSquaresMinimizer(Int_t numdatapoints, Double_t datapointxs[], Double_t datapointys[], Double_t datapointweights[], Double_t errorys[], Double_t &fit_gradient, Double_t &fit_offset, Double_t &chi2);
	// check a prospective trajectory hits all layers
	Bool_t AngValid(Int_t layerstart, Double_t angle, Int_t MaxMin, Int_t axis);
	// check a trajectory hits a given layer
	Bool_t CheckIntersection(Int_t layerstart, Int_t layerend, Double_t angle, Int_t MaxMin, Int_t axis);
	
	void CalculateKEstart();				// based on ? depth?
	void CalculateKEend();					// fully penetrating: estimate from dE/dx? Tank_E - E_loss?
	void CalculateParticlePID();			// based on rate of loss? num tracks..? penetration? tank? 
	
	// Default Constructor
	// ====================
	public:
	// Default constructor that initialises all private members required for ROOT classes
	cMRDSubEvent() : MRDSubEventID(-1), wcsimfile(""), run_id(-1), event_id(-1), subtrigger(-1), digi_ids(), pmts_hit(), digi_qs(), digi_ts(), digi_numphots(), digi_phot_ts(), digi_phot_parents(), tanktrackID(-1), layers_hit(), eDepsInLayers(), KEStart(-1.), KEEnd(-1.), particlePID(-1), tracktype(-1), trueTrackID(-1), recovalshoriz(), recovalsvert() {};
	
	// destructor
	~cMRDSubEvent(){}
	
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
	tanktrackID(-1), layers_hit(), KEStart(-1.), KEEnd(-1.), particlePID(-1),
	tracktype(-1), trueTrackID(-1), recovalshoriz(), recovalsvert() {
		
#ifdef _MRDSubEvent_VERBOSE_
		cout<<"constructing a track with "<<digi_ids.size()<<" digits"<<endl;
#endif
//		if(fillstaticmembers){
//			// fill static members
//			std::vector<Int_t> temp(aspectrum, aspectrum+19);
//			//aspectrumv.assign(temp.rbegin(), temp.rend());
//			aspectrumv.assign(temp.begin(), temp.end());
//			fillstaticmembers=false;
//		}
		eDepsInLayers.assign(numpanels, 0.);	// can't assign the size in the class def. 
		DoReconstruction();
		DrawMrdCanvases();
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
	
	std::pair<double, double> xupcorner1, xupcorner2, xdowncorner1, xdowncorner2, yupcorner1, yupcorner2, ydowncorner1, ydowncorner2;
	
	
	// End class definition
	// ====================
	ClassDef(cMRDSubEvent,1);					// INCREMENT VERSION NUM EVERY TIME CLASS MEMBERS CHANGE
};

Bool_t cMRDSubEvent::fillstaticmembers=true;
TCanvas* cMRDSubEvent::imgcanvas=0;
TText* cMRDSubEvent::titleleft=0;
TText* cMRDSubEvent::titleright=0;
std::vector<TBox*> cMRDSubEvent::paddlepointers(nummrdpmts);
//std::vector<Int_t> cMRDSubEvent::aspectrumv(19);
std::vector<Int_t> cMRDSubEvent::aspectrumv = ( []()->std::vector<Int_t> { std::vector<Int_t> temp {kYellow, kOrange, (kOrange-3), (kOrange+8), (kOrange+10), kRed, (kRed+1), (kPink+4), (kMagenta+2), (kMagenta+1), kMagenta, (kViolet-2), (kViolet-3), (kViolet+7), (kViolet+9), (kBlue+2), (kBlue+1), kAzure, (kAzure+7)}; return temp; }() );

#include "MRDSubEvent_DoReconstruction2.cxx"	// contains reconstruction function definitions
#include "makemrdimage.cxx"					// functions to draw the MRD top and side views

#endif

#ifdef __CINT__
#pragma link C++ class cMRDSubEvent+;
//#pragma link C++ class ROOT::Math::XYZTVector+;
//#pragma link C++ class std::vector<ROOT::Math::XYZTVector>+;
//#pragma link C++ class cMRDStrike+;
#pragma link C++ class std::vector<cMRDStrike>+;
#endif
