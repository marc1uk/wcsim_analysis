/* vim:set noexpandtab tabstop=2 wrap */
#ifndef _MRDDigit_VERBOSE_
//#define _MRDDigit_VERBOSE_ 1
#endif

#ifndef _MRDDigit_Class_
#define _MRDDigit_Class_ 1

#include <TObject.h>
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include "Math/GenVector/LorentzVector.h"
#include <exception>	// for stdexcept
#include <vector>
#include <utility>
#include "TMath.h"

#include "MRDspecs.hh"
extern Int_t numpanels;
extern Int_t numpaddlesperpanelh;
extern Int_t numpaddlesperpanelv;
extern Double_t scintzedges[24];

class cMRDdigit : public TObject {
// a cMRDStrike' is old sandbox terminology; this class is being repurposed to hold useful information about a 'digit' in the MRD. 
	private: 	
	//Double_t xEstimate;			// this information is stored in the points vector in MRDTrackClass
	//Double_t yEstimate;
	//Double_t zEstimate;
	//Double_t tEstimate;
	Int_t MRDlayer;
	Int_t MRDpaddle;
	Int_t PMTnum;				// assume photons are confined to a paddle.
	Double_t eDeposited;			// use wavelength of photon for each to obtain energy
	std::pair<Double_t, Double_t> xrange;	// from width of panel and whether multiple PMTs were hit.
	std::pair<Double_t, Double_t> yrange;
	std::pair<Double_t, Double_t> zrange;	// from depth of panel
	std::pair<Double_t, Double_t> trange; 	// from uncertainty in PMT timing resoluton
	std::vector<Double_t> PMThitTimes;	// vector of photon hit times.
	
	public:

	
	// Constructors:
	cMRDdigit() : MRDlayer(-1), MRDpaddle(-1), PMTnum(-1), eDeposited(0), xrange(), yrange(), zrange(), trange(), PMThitTimes() {}
	// default constructor.
	~cMRDdigit(){}	// nothing allocated with new, nothing to delete.
	
	//TODO: cMRDdigit(thisdigitstime, thisdigitsq, thisdigitstubeid, photontimesinatrack,particleidsinatrack);
	cMRDdigit(std::vector<Double_t> PMThitTimesin, Int_t PMTnumin) : MRDlayer(-1), MRDpaddle(-1), eDeposited(0), xrange(), yrange(), 
	zrange(), trange(), PMThitTimes(PMThitTimesin), PMTnum(PMTnumin) {
	// actual constructor. 
#ifdef _MRDDigit_VERBOSE_
	cout<<"constructing a mrd digit"<<endl;
#endif
		std::sort(PMThitTimes.begin(), PMThitTimes.end());
//	std::vector<Int_t> hitpmts = GetNumPMTsHit();			// should only be one PMT hit per mrd digit.
//	for(std::vector<Int_t>::iterator it=hitpmts.begin();it!=hitpmts.end();it++){
#ifdef _MRDDigit_VERBOSE_
		cout<<"there are "<<PMThitTimes.size()<<" hits in this mrd digit"<<endl;
#endif
		if(PMThitTimes.size()!=0){
			Int_t panelpairnum = TMath::Floor(PMTnum/(numpaddlesperpanelh+numpaddlesperpanelv));
			Int_t pmtnumrem = PMTnum-((numpaddlesperpanelh+numpaddlesperpanelv)*panelpairnum);
			// first layer is a vertical layer (check this!) //TODO
			if(pmtnumrem<numpaddlesperpanelv){
#ifdef _MRDDigit_VERBOSE_
				cout<<"this is a vertical paddle"<<endl;
#endif
				MRDlayer=(panelpairnum*2);
			} else {
				MRDlayer=(panelpairnum*2)+1; 
				pmtnumrem = pmtnumrem - numpaddlesperpanelv;
#ifdef _MRDDigit_VERBOSE_
				cout<<"this is a horizontal paddle"<<endl;
#endif
			}
			MRDpaddle = pmtnumrem;
#ifdef _MRDDigit_VERBOSE_
			cout<<"the hit PMT was PMT number "<<PMTnum<<" which is paddle "<<MRDpaddle<<" of layer "<<MRDlayer<<endl;
#endif
			// note within a panel PMTs {0,1}, {2,3} etc are pairs stacked end-to-end
			if(MRDlayer%2!=0){	//horizontal panels, use PMT number to determine y position
				if((pmtnumrem%2)==0){	// row 1 - RHS looking downstream
					pmtnumrem=pmtnumrem/2;	// merge rows
					xrange.first = 0;
					xrange.second = scintvfullylen;
					yrange.first = (scintfullxlen+scintbordergap)*(pmtnumrem-(numpaddlesperpanelh/4)-0.5);
					yrange.second = (scintfullxlen+scintbordergap)*(pmtnumrem+1-(numpaddlesperpanelh/4)-0.5)-scintalugap;
				} else {
					pmtnumrem=TMath::Floor(pmtnumrem/2);	// merge rows
					xrange.first = -scintvfullylen;
					xrange.second = 0;
					yrange.first = (scintfullxlen+scintbordergap)*(pmtnumrem-(numpaddlesperpanelh/4)-0.5);
					yrange.second = (scintfullxlen+scintbordergap)*(pmtnumrem+1-(numpaddlesperpanelh/4)-0.5)-scintalugap;
				}			
			} else {
				if((pmtnumrem%2)==0){	// row 1 or 2?
					pmtnumrem=pmtnumrem/2;	// merge rows
					yrange.first = -scinthfullylen;
					yrange.second = 0;
					xrange.first = (scintfullxlen+scintbordergap)*(pmtnumrem-(numpaddlesperpanelv/4)-0.5);
					xrange.second = (scintfullxlen+scintbordergap)*(pmtnumrem+1-(numpaddlesperpanelv/4)-0.5)-scintalugap;
				} else {
					pmtnumrem=TMath::Floor(pmtnumrem/2);	// merge rows
					yrange.first = 0;
					yrange.second = scinthfullylen;
					xrange.first = (scintfullxlen+scintbordergap)*(pmtnumrem-(numpaddlesperpanelv/4)-0.5);
					xrange.second = (scintfullxlen+scintbordergap)*(pmtnumrem+1-(numpaddlesperpanelv/4)-0.5)-scintalugap;
				}
			}
//			zrange.first = scintzedges[(Int_t)(2*TMath::Floor(MRDlayer))];
//			zrange.second = scintzedges[(Int_t)(2*TMath::Floor(MRDlayer))+1];
			Double_t Zposition = MRDlayer*(steelfullzlen + alufullzlen + scintfullzlen + layergap);
			if(MRDlayer==0){Zposition = Zposition + alufullzlen + scintalugap;}		// layer 0 scints intrude into layer 1's offset
			Zposition = Zposition + steelfullzlen + steelscintgap;								// scint follows closely behind first steel
			Zposition = Zposition + (scintfullzlen/2);														// offset by half depth so we place front face not centre
			Zposition = Zposition - (mrdZlen/2);																	// offset by half total length to shift to front
			zrange.first = Zposition;
			zrange.second = Zposition + scintfullzlen;
#ifdef _MRDDigit_VERBOSE_
			cout<<"the paddle bounds are from "<<xrange.first<<"<x<"<<xrange.second<<", "<<yrange.first<<"<y<"<<yrange.second
					<<", "<<zrange.first<<"<z<"<<zrange.second<<endl;
#endif
			trange.first = PMThitTimes.front();
			trange.second = PMThitTimes.back();
#ifdef _MRDDigit_VERBOSE_
			cout<<"hit times range from "<<trange.first<<"<t<"<<trange.second<<"ns"<<endl;
#endif
		}
	}
	
	ClassDef(cMRDdigit,1);	// INCREMENT VERSION NUM EVERY TIME CLASS MEMBERS CHANGE
};
#endif

#ifdef __CINT__
#pragma link C++ class cMRDdigit+;
//#pragma link C++ class std::pair<Double_t, Double_t>+;
//#pragma link C++ class std::vector<Double_t>+;
#endif
