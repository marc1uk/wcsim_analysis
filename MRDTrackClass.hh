#ifndef _MRDTrack_VERBOSE_
#define _MRDTrack_VERBOSE_ 1
#endif

#ifndef _MRDTrack_Class_
#define _MRDTrack_Class_ 1

// TODO: what should we do about multiple strikes in the same track on the same pmt? 
// this would screw with reconstruction...  should we just add those hits to the same strike?

#include <TObject.h>
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TVector3.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include "Math/GenVector/LorentzVector.h"
#include <exception>	// for stdexcept
#include <vector>
#include <map>

//#include "MRDStrikeClass.hh"
#include "MRDspecs.hh"
extern Int_t numpanels;
extern Int_t numpaddlesperpanelh;
extern Int_t numpaddlesperpanelv;

class cMRDTrack : public TObject {
	
	// Private members
	// ===============
	private:
	Int_t MRDtrackID;		// an ID to uniquely identify this track within the MRD.
	Int_t tanktrackID;	// if this track has been correlated with a track in the tank, store the tank track ID here.
	
	// Raw Info:
	std::vector<cMRDStrike> strikes;			//vector of cMRDStrike class objects. 
	// strike class objects contain all information derivable from just one strike object - layer number, number and vector of 
	// PMT hits, which PMTs were hit, hit times, energy depositions, positional error from struck PMTs... 
	std::vector<Double_t> eDepsInLayers;		// fixed length vector of energy deposited in each layer 
	Int_t valid;	// are all strikes consistent with one track? -1: u nchecked, 0: inconsistent, 1: consistent
	Int_t tankconsistency;	// are strikes consistent with associated tank track, if there is one? (-2, no tank track)
	
	// Calculated Info
	Int_t particlePID;	// particle type code, based on rate of energy loss, penetration... other?
	Double_t KEStart;		// projected from energy loss? what about from tank info?
	Double_t KEEnd;			// projected from energy loss, if possible. (-1; not enough info)
	std::vector<ROOT::Math::XYZTVector> mrdpoints;	// estimated times and positions of interaction points, or 'strikes'.
	std::vector<ROOT::Math::XYZTVector> tankpoints;	// projected crossing points from tank track, if available
	std::map<const char*, double> recovalshoriz, recovalsvert; // has keys "xycrossing", "zcrossing", "angmin" and "angmax"
	
	public:
	// Setters
	// =======
	void AppendStrike(cMRDStrike strikein){strikes.push_back(strikein);}
//void SetTankEntryPoint(ROOT::Math::XYZTVector tankStartin){tankStart = tankStartin;} // entry point derived from a matched tank track
	void SetparticlePID(Int_t particlePIDin){particlePID = particlePIDin;}		// could also be determined from tank info
	// void SetKEEnd(Double_t KEEndin){KEEnd = KEEndin;}											// should there be setters for other properties?
	// void SetKEStart(Double_t KEStartin){KEStart = KEStartin;}
	
	// Getters
	// =======
	Int_t GetparticlePID(){return particlePID;}				// reconstruct from rate of ionisation...? likelihood? farm out to a function.
	Double_t GetKEStart(){return KEStart;}						// reconstruct from rate of energy deposition
	Double_t GetKEEnd(){return KEEnd;}
	Double_t GetTotalEdep(){return KEStart-KEEnd;}
	std::vector<Double_t> GetEdeps(){return eDepsInLayers;}
	
	Int_t GetNumStrikes(){return strikes.size();}
	std::vector<cMRDStrike> GetStrikes(){ return strikes; }
	cMRDStrike* GetFirstStrike(){return &(strikes.front());}	// should check we have at least one.
	cMRDStrike* GetLastStrike(){return &(strikes.back());}		// should check we have at least one.
	cMRDStrike* GetStrike(Int_t i){
		try{return &(strikes.at(i));}
		catch(const std::out_of_range& oor){return 0;}
	}
	
	std::vector<ROOT::Math::XYZTVector> GetMRDPoints(){return mrdpoints;}
	Double_t GetPenetration(){return mrdpoints.back().Z()-mrdpoints.front().Z();}
	ROOT::Math::XYZTVector* GetFirstPoint(){return &(mrdpoints.front());}		// should check we have at least one.
	ROOT::Math::XYZTVector* GetLastPoint(){return &(mrdpoints.back());}			// should check we have at least one.
	ROOT::Math::XYZTVector* GetPoint(Int_t i){
		try{return &(mrdpoints.at(i));}
		catch(const std::out_of_range& oor){return 0;}
	}
	
	std::map<const char*, double> GetTrackXStats(){ return recovalshoriz; }
	std::map<const char*, double> GetTrackYStats(){ return recovalsvert; }
	
	public:
	// Constructors & Destructors
	// ==========================
	// ##### Must provide a default constructor that initialises all private members!! #####
	cMRDTrack() : MRDtrackID(-1), tanktrackID(-1), valid(-1), tankconsistency(-2), particlePID(-1), KEStart(0.), KEEnd(0.), strikes(), eDepsInLayers(numpanels, 0.), mrdpoints(), tankpoints(), recovalshoriz(), recovalsvert() {}
	~cMRDTrack(){}	// vector container automatically handles resource deallocation: destruction of points and strikes
	cMRDTrack(std::vector<cMRDStrike> strikesin) : MRDtrackID(-1), tanktrackID(-1), valid(-1), tankconsistency(-2), particlePID(-1), KEStart(0.), KEEnd(0.), strikes(strikesin), eDepsInLayers(numpanels, 0.), mrdpoints(), tankpoints(), recovalshoriz(), recovalsvert() {
#ifdef _MRDTrack_VERBOSE_
		cout<<"constructing a track with "<<strikes.size()<<" strikes"<<endl;
#endif
		DoReconstruction();	// automatically do reconstruction on creation, since passed a vector of strikes.
	}
	
	// Private functions for calculating particle/track properties
	private:
	void DoReconstruction();	// Main track reconstruction code
	Bool_t AngValid(Int_t layerstart, Double_t angle, Int_t MaxMin, Int_t axis);	// check a prospective trajectory hits all layers
	Bool_t CheckIntersection(Int_t layerstart, Int_t layerend, Double_t angle, Int_t MaxMin, Int_t axis);	// check a trajectory hits a given layer

	ClassDef(cMRDTrack,1);	// INCREMENT VERSION NUM EVERY TIME CLASS MEMBERS CHANGE
};

// ********* trajectory reconstruction: main *****************
	void cMRDTrack::DoReconstruction(){
#ifdef _MRDTrack_VERBOSE_
	cout<<"doing track reconstruction"<<endl;
#endif
	/* OK, this is going to generate two maps, containing keys: (xycrossing, zcrossing, angmax, angmin). 
	Each map contains the elements that define a region of acceptance in one plane - starting from the zcrossing and xycrossing
	values, (defining a line in the x/y plane) and projecting in the y/x plane at angles angmax and angmin defines a wedge.
	This wedge, projected back to the tank, defines a region of possible starting points and angles (coupled - tracks must lead
	to the vertex) within the tank.
	The two sets of constraints in x and y are independent and define a non-trivial geometry when combined.*/

	std::vector<TVector3> bottompoints;
	std::vector<TVector3> toppoints;
	std::vector<Int_t> strucklayers;
	for(Int_t i=0;i<strikes.size();i++){
		cMRDStrike thisstrike=strikes.at(i);
		// bottompoints stores most negative values: smallest z, most negative x, y. toppoints the opposite.
		bottompoints.push_back(TVector3(thisstrike.GetXrange().first,thisstrike.GetYrange().first,thisstrike.GetZrange().first));
		toppoints.push_back(TVector3(thisstrike.GetXrange().second,thisstrike.GetYrange().second,thisstrike.GetZrange().second));
		Int_t strucklayer = thisstrike.GetLayerNum();
		strucklayers.push_back(strucklayer);
		eDepsInLayers.at(strucklayer)+=thisstrike.GetEdeposited();
	}
	Double_t xangmax = (0.0/0.0);	// set to NaN angmax represents angle of uppermost limit of valid trajectories
	Double_t xangmin = (0.0/0.0);	// set to NaN. angmin represents angle of lowermost limit of valid trajectories 
	Double_t yangmax = (0.0/0.0);
	Double_t yangmin = (0.0/0.0);
	Double_t xlayermax = -1;
	Double_t xlayermin = -1;
	Double_t ylayermax = -1;
	Double_t ylayermin = -1;
	
	// calculating span of possible angles that the particle could have traversed to hit all paddles.
	// this assumes paddles are infinitely thin - 0.6mm so good approximation. 
	// angmin is the angle of the trajectory that defines the lower bound of the backward projection
	// it may be either upgoing or downgoing, depending on the paddles hit. 
	// angmax defines the trajectory of the upper bound of the backward projected 'cone' 
	
	// angmin:
	// scan through from all bottom points to all top points to find the greatest upward angle a particle may
	// have had and still hit all recorded layers (note: this may be negative, indicating a particle was
	// going downward - in this case angmin is the angle at which it went downward _least_ steeply)
	// angle is measured from the positive z axis
	//for(std::vector<Int_t>::iterator startlayerit=strucklayers.begin(); startlayerit!=strucklayers.end(); ++startlayerit){
	//	Int_t startlayer=(*startlayerit);
	//	for(std::vector<Int_t>::iterator nextlayerit=++startlayerit; nextlayerit!=strucklayers.end(); ++nextlayerit){
	//		Int_t nextlayer=(*nextlayerit);
	for(int startlayer=0; startlayer<strikes.size(); startlayer++){
		for(int nextlayer=(startlayer+1); nextlayer<strikes.size(); nextlayer++){
			Double_t opp = toppoints.at(nextlayer).X()-bottompoints.at(startlayer).X();
			Double_t adj = bottompoints.at(nextlayer).Z()-bottompoints.at(startlayer).Z();	// really it should be an average for Z points
			Double_t ang = TMath::ATan(opp/adj);
#ifdef _MRDTrack_VERBOSE_
			cout<<"checking most upgoing angle from bottom of layer "<<strucklayers.at(startlayer)<<" at x="
					<<bottompoints.at(startlayer).X()<<", z="<<bottompoints.at(startlayer).Z()<<", to top of layer "
					<<strucklayers.at(nextlayer)<<" at x="<<toppoints.at(nextlayer).X()<<", z="<<bottompoints.at(nextlayer).Z()
					<<", at angle= "<<TMath::RadToDeg()*ang<<"degs"<<endl;
#endif
			if(TMath::IsNaN(xangmin)||(ang>xangmin)){
				if(AngValid(startlayer,ang,0, 0)){	// check this angle actually hits all paddles
#ifdef _MRDTrack_VERBOSE_
					cout<<"setting this x angle as the new most upgoing"<<endl;
#endif
					xangmin=ang;										// set this as the lower bound angle
					xlayermin=startlayer;					// to define a trajectory we also need a starting point
				}
			}
		}
	}
#ifdef _MRDTrack_VERBOSE_
	cout<<"finished calculating most upgoing x angle"<<endl;
#endif
	if(TMath::IsNaN(xangmin)){cout<<" ####### NO UPGOING X ANGLE FOUND ########## "<<endl;}
	
	// angmax: same, this time scan from top of each layer to bottom of next layer, for the uppermost entering trajectory
	for(Int_t startlayer=0;startlayer<strikes.size();startlayer++){
		for(Int_t nextlayer=(startlayer+1);nextlayer<strikes.size(); nextlayer++){
			Double_t opp = toppoints.at(startlayer).X() - bottompoints.at(nextlayer).X();	// negative if upgoing
			Double_t adj = bottompoints.at(nextlayer).Z()-bottompoints.at(startlayer).Z();
			Double_t ang = TMath::ATan(opp/adj);										// negative if downgoing
#ifdef _MRDTrack_VERBOSE_
			cout<<"checking most downgoing angle from top of layer "<<strucklayers.at(startlayer)<<" at x="
			<<toppoints.at(startlayer).X()<<", z="<<bottompoints.at(startlayer).Z()<<" to bottom of layer "
			<<strucklayers.at(nextlayer)<<", at x="<<bottompoints.at(nextlayer).X()<<", z="<<bottompoints.at(nextlayer).Z()
			<<", at angle= "<<TMath::RadToDeg()*ang<<endl;
#endif
			if(TMath::IsNaN(xangmax)||(ang>xangmax)){
				if(AngValid(startlayer,ang,1, 0)){							// project both forward and back from bottom vertex of
					xangmax=ang;																// nextlayer to all other layers and check within limits
					xlayermax=startlayer;
#ifdef _MRDTrack_VERBOSE_
					cout<<"setting this x angle as the new most downgoing"<<endl;
#endif
				}
			}
		}
	}
#ifdef _MRDTrack_VERBOSE_
	cout<<"finished calculating most downgoing x angle"<<endl;
#endif
	if(TMath::IsNaN(xangmax)){cout<<" ####### NO DOWNGOING X ANGLE FOUND ########## "<<endl;}
	
	//now find angles and layers for y projections
	for(Int_t startlayer=0;startlayer<strikes.size();startlayer++){
		for(Int_t nextlayer=startlayer+1;nextlayer<strikes.size(); nextlayer++){
			Double_t opp = toppoints.at(nextlayer).Y()-bottompoints.at(startlayer).Y();
			Double_t adj = bottompoints.at(nextlayer).Z()-bottompoints.at(startlayer).Z();	// really it should be an average for Z points
			Double_t ang = TMath::ATan(opp/adj);
#ifdef _MRDTrack_VERBOSE_
			cout<<"checking most upgoing angle from bottom of layer "<<strucklayers.at(startlayer)<<" at y="
					<<bottompoints.at(startlayer).Y()<<", z="<<bottompoints.at(startlayer).Z()<<", to top of layer "
					<<strucklayers.at(nextlayer)<<" at y="<<toppoints.at(nextlayer).Y()<<", z="<<bottompoints.at(nextlayer).Z()
					<<", at angle= "<<TMath::RadToDeg()*ang<<"degs"<<endl;
#endif
			if(TMath::IsNaN(yangmin)||(ang>yangmin)){
				if(AngValid(startlayer,ang,0, 1)){	// check this angle actually hits all paddles
					yangmin=ang;										// set this as the lower bound angle
					ylayermin=startlayer;					// to define a trajectory we also need a starting point
#ifdef _MRDTrack_VERBOSE_
					cout<<"setting this y angle as the new most upgoing"<<endl;
#endif
				}
			}
		}
	}
#ifdef _MRDTrack_VERBOSE_
	cout<<"finished calculating most upgoing y angle"<<endl;
#endif
	if(TMath::IsNaN(yangmin)){cout<<" ####### NO UPGOING Y ANGLE FOUND ########## "<<endl;}
	
	//now find angles and layers for y projections
	for(Int_t startlayer=0;startlayer<strikes.size();startlayer++){
		for(Int_t nextlayer=startlayer+1;nextlayer<strikes.size(); nextlayer++){
			// angmax: same, this time scan from top of each layer to bottom of next layer, for the uppermost entering trajectory
			Double_t opp = toppoints.at(startlayer).Y() - bottompoints.at(nextlayer).Y();	// negative if upgoing
			Double_t adj = bottompoints.at(nextlayer).Z()-bottompoints.at(startlayer).Z();
			Double_t ang = TMath::ATan(opp/adj);																					// negative if upgoing
#ifdef _MRDTrack_VERBOSE_
			cout<<"checking most downgoing angle from top of layer "<<strucklayers.at(startlayer)<<" at y="
			<<toppoints.at(startlayer).Y()<<", z="<<bottompoints.at(startlayer).Z()<<" to bottom of layer "
			<<strucklayers.at(nextlayer)<<", at y="<<bottompoints.at(nextlayer).Y()<<", z="<<bottompoints.at(nextlayer).Z()
			<<", at angle= "<<TMath::RadToDeg()*ang<<endl;
#endif
			if(TMath::IsNaN(yangmax)||(ang>yangmax)){
				if(AngValid(startlayer,ang,1, 1)){							// project both forward and back from bottom vertex of
					yangmax=ang;																	// nextlayer to all other layers and check within limits
					ylayermax=startlayer;
#ifdef _MRDTrack_VERBOSE_
					cout<<"setting this y angle as the new maximum"<<endl;
#endif
				}
			}
		}
	}
#ifdef _MRDTrack_VERBOSE_
	cout<<"finished calculating most downgoing y angle"<<endl;
#endif
	if(TMath::IsNaN(yangmax)){cout<<" ####### NO DOWNGOING Y ANGLE FOUND ########## "<<endl;}

	recovalshoriz.emplace("angminzval",bottompoints.at(xlayermin).Z());
	recovalshoriz.emplace("angminxval",bottompoints.at(xlayermin).X());
	recovalshoriz.emplace("angmin", xangmin);
	recovalshoriz.emplace("angmaxzval",bottompoints.at(xlayermax).Z());
	recovalshoriz.emplace("angmaxxval",toppoints.at(xlayermax).X());
	recovalshoriz.emplace("angmax", xangmax);

	recovalsvert.emplace("angminzval",bottompoints.at(ylayermin).Z());
	recovalsvert.emplace("angminyval",bottompoints.at(ylayermin).Y());
	recovalsvert.emplace("angmax", yangmax);
	recovalsvert.emplace("angmaxzval",bottompoints.at(ylayermax).Z());
	recovalsvert.emplace("angmaxyval",toppoints.at(ylayermax).Y());
	recovalsvert.emplace("angmin", yangmin);
	
	
}

// ********* trajectory reconstruction: scan projections to all layers *****************
Bool_t cMRDTrack::AngValid(Int_t layerstart, Double_t angle, Int_t MaxMin, Int_t axis){
#ifdef _MRDTrack_VERBOSE_
	cout<<"checking if ";
	if(axis){cout<<"y ";} else {cout<<"x ";}
	cout<<"trajectory, going ";
	if(MaxMin){cout<<"downwards ";} else {cout<<"upwards ";}
	cout<<"from layer "<<layerstart<<" at angle "<<TMath::RadToDeg()*angle<<" hits other layers"<<endl;
#endif
	Bool_t angvalid=true;
	for(Int_t layerend=0;layerend<strikes.size();layerend++){
		if(layerstart==layerend){continue;}	// don't check from a layer to itself
		angvalid = CheckIntersection(layerstart, layerend, angle, MaxMin, axis);
		if(!angvalid){
#ifdef _MRDTrack_VERBOSE_
			cout<<"trajectory does not hit layer "<<layerend<<endl;
#endif
			return false;
		}
	}
#ifdef _MRDTrack_VERBOSE_
	cout<<"trajectory hits all layers"<<endl;
#endif
	if(axis){cout<<"y ";} else {cout<<"x ";}
	cout<<"trajectory, going ";
	if(MaxMin){cout<<"downwards ";} else {cout<<"upwards ";}
	cout<<"from layer "<<layerstart<<" at angle "<<TMath::RadToDeg()*angle<<" hits all other layers"<<endl;
	return true;
}

// ********* trajectory reconstruction: check single projection *****************
Bool_t cMRDTrack::CheckIntersection(Int_t layerstart, Int_t layerend, Double_t angle, Int_t MaxMin, Int_t axis){
	Double_t z1 = strikes.at(layerstart).GetZrange().first;
	Double_t z2 = strikes.at(layerend).GetZrange().first;
	Double_t zstart;
	Double_t zend;
	Int_t LHlayer;
	Int_t RHlayer;
	if(z1<z2){
		LHlayer = layerstart;
		RHlayer = layerend;
		zstart=z1; 
		zend=z2;
	} else {
		LHlayer = layerend;
		RHlayer = layerstart;
		zstart=z2; 
		zend=z1;
	}
#ifdef _MRDTrack_VERBOSE_
	cout<<"checking if projection from zstart="<<z1;
	if(axis){
		cout<<", ystart=";
		if(MaxMin==1){
			cout<<strikes.at(layerstart).GetYrange().second; 
			if(z2<z1){cout<<", going backwards and upwards to zend=";} else {cout<<", going downwards to zend=";}
		} else {
			cout<<strikes.at(layerstart).GetYrange().first; 
			if(z2<z1){cout<<" going backwards and downwards to zend=";} else {cout<<", going upwards to zend=";}
		}
		cout<<z2<<" at angle "<<TMath::RadToDeg()*angle<<" lies within range "<<strikes.at(layerend).GetYrange().first
				<<"<y<"<<strikes.at(layerend).GetYrange().second<<endl;
	} else {
		cout<<", xstart=";
		if(MaxMin==1){
			cout<<strikes.at(layerstart).GetXrange().second; 
			if(z2<z1){cout<<", going backwards and upwards to zend=";} else {cout<<", going downwards to zend=";}
		} else {
			cout<<strikes.at(layerstart).GetXrange().first;
			if(z2<z1){cout<<", going backwards and downwards to zend=";} else {cout<<", going upwards to zend=";}
		}
		cout<<z2<<" at angle "<<TMath::RadToDeg()*angle<<" lies within range "<<strikes.at(layerend).GetXrange().first
				<<"<x<"<<strikes.at(layerend).GetXrange().second<<endl;
	}
#endif
	Double_t ydiff = (zend-zstart)*TMath::Tan(angle);
	Double_t yproj;
	if(MaxMin==1){	//downgoing angle
	// layerstart defines a maximum angle - check projection from it's top to layerend at 'angle'
	// downwards (or at 'angle' upwards if projecting backward) strikes layerend
		if(z1==zstart){	// layerstart is on LHS - project forwards (-'ve y displacement)
			if(axis){	// y axis
				yproj = strikes.at(layerstart).GetYrange().second - ydiff;
			} else {	// x axis
				yproj = strikes.at(layerstart).GetXrange().second - ydiff;
			}
		} else {		// layerstart is on RHS - project backwards (+'ve y displacement)
			if(axis){
				yproj = strikes.at(layerstart).GetYrange().second + ydiff;
			} else {
				yproj = strikes.at(layerstart).GetXrange().second + ydiff;
			}
		}
	} else {	//upgoing angle
	// layerstart defines a min angle - check projection from it's bottom at angle upwards 
	// (if projecting forwards, or downwards if backwards) to layerend, strikes within 
	// layerend's y boundaries
		if(z1==zstart){	// layerstart is on LHS - project forwards (+'ve y displacement)
			if(axis){
				yproj = strikes.at(layerstart).GetYrange().first + ydiff;
			} else {
				yproj = strikes.at(layerstart).GetXrange().first + ydiff;
			}
		} else {		// layerstart is on RHS - project backwards (-'ve y displacement)
			if(axis){
				yproj = strikes.at(layerstart).GetYrange().first - ydiff;
			} else {
				yproj = strikes.at(layerstart).GetXrange().first - ydiff;
			}
		}
	}
	Double_t upperbound, lowerbound;
	if(axis){
		upperbound=strikes.at(layerend).GetYrange().first;
		lowerbound=strikes.at(layerend).GetYrange().second;
	} else {
		upperbound=strikes.at(layerend).GetXrange().first;
		lowerbound=strikes.at(layerend).GetXrange().second;
	}
	
	if( ((yproj+0.05) >= upperbound)&&((yproj-0.05) <= lowerbound)){
#ifdef _MRDTrack_VERBOSE_ 
		cout<<"projected point is "<<yproj<<"; trajectory valid"<<endl;
#endif
		return true; 
	}	else { 
#ifdef _MRDTrack_VERBOSE_ 
	cout<<"projected point is "<<yproj<<"; trajectory does not cross"<<endl;
#endif
	return false; 
	}
}

#endif

#ifdef __CINT__
#pragma link C++ class cMRDTrack+;
//#pragma link C++ class ROOT::Math::XYZTVector+;
//#pragma link C++ class std::vector<ROOT::Math::XYZTVector>+;
//#pragma link C++ class cMRDStrike+;
#pragma link C++ class std::vector<cMRDStrike>+;
#endif
