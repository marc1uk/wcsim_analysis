/* vim:set noexpandtab tabstop=4 wrap */
#include "TROOT.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TEveLine.h"	// << this ALSO needs to be included by gROOT->ProcessLine BEFORE compiling
#include "TStyle.h"
#include "TGeoManager.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>

void gdml_test() {
	// Import the gdml geometry for the detector:
	//TString gdmlpath = "/pnfs/annie/scratch/users/moflaher/annie_v04_aligned.gdml";
	TString gdmlpath = "/home/marc/anniegpvm/wcsim/root_work/annie_v04_aligned.gdml";
	TGeoManager::Import(gdmlpath);
	// optional: make all the materials semi-transparent so we can see the tracks going through them:
	TList* matlist = gGeoManager->GetListOfMaterials();
	TIter nextmaterial(matlist);
	while(TGeoMixture* amaterial=(TGeoMixture*)nextmaterial()) amaterial->SetTransparency(90);
	gGeoManager->GetVolume("WORLD2_LV")->Draw("ogl");
	TCanvas* gdmlcanv = (TCanvas*)gROOT->FindObject("c1");
	
	std::vector< std::pair<TVector3,TVector3> > sometracks{
	std::make_pair(TVector3(-91.35, -105.56, 336.08),TVector3(-172.55, -89.3194, 432.96)),
	std::make_pair(TVector3(0, 20.2995, 336.08),TVector3(0, -10.1492, 372.41)),
	std::make_pair(TVector3(0, -73.0802, 336.08),TVector3(0, -54.8097, 445.07)),
	std::make_pair(TVector3(-30.4558, 101.498, 336.08),TVector3(50.7547, 101.498, 432.96)),
	std::make_pair(TVector3(-101.5, 44.6602, 336.08),TVector3(-101.5, 20.2998, 432.96)),
	};
	
	std::vector<TEveLine*> thiseventstracks;
	for(auto&& apair : sometracks){
		TEveLine* evl = new TEveLine("track1",2);
		evl->SetLineWidth(4);
		evl->SetLineStyle(1);
		evl->SetMarkerColor(kRed);
		evl->SetRnrPoints(kTRUE);  // enable rendering of points
		TVector3 sttv = apair.first;
		TVector3 stpv = apair.second;
		evl->SetPoint(0,sttv.X(),sttv.Y(),sttv.Z());
		evl->SetPoint(1,stpv.X(),stpv.Y(),stpv.Z());
		thiseventstracks.push_back(evl);
	}
	
	gdmlcanv->cd();
	for(auto&& aline : thiseventstracks){
		cout<<"drawing track at "<<aline<<" from ("<<aline->GetLineStart().fX
			<<", "<<aline->GetLineStart().fY<<", "<<aline->GetLineStart().fZ<<") to ("
			<<aline->GetLineEnd().fX<<", "<<aline->GetLineEnd().fY<<", "<<aline->GetLineEnd().fZ<<")"
			<<endl;
		aline->Draw();
	}
	gdmlcanv->Update();
	gPad->WaitPrimitive();
	for(auto&& aline : thiseventstracks) delete aline; 
}
