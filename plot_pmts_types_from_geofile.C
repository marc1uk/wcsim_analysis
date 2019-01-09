/* vim:set noexpandtab tabstop=2 wrap */
//C++
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <thread>
#include <chrono>
#include <time.h>
//ROOT
#include <TSystem.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TPolyMarker3D.h>
#include <TPointSet3D.h>   // use TPointSet3D over TPolyMarker3D for ogl viewer (needed for gdml overlay)
#include <TLegend.h>
#include <TGeoManager.h>
#include <TGLViewer.h>
//CUSTOM
#include "ColourWheel.hh"
//WCSIM
#include "../wcsim/include/WCSimRootEvent.hh"
#include "../wcsim/include/WCSimRootGeom.hh"
#include "../wcsim/include/WCSimPmtInfo.hh"
#include "../wcsim/include/WCSimLAPPDInfo.hh"
#include "../wcsim/include/WCSimEnumerations.hh"
#include "../wcsim/include/WCSimRootLinkDef.hh"
#include "../wcsim/include/WCSimRootOptions.hh"

// setup static members for ColourWheel
const std::vector<EColor> ColourWheel::colours{kBlack, kBlue, (EColor)TColor::GetColorDark(kGreen), kRed, kViolet, kOrange, kMagenta,(EColor)(kAzure+2),(EColor)(kOrange+4),(EColor)(kViolet-6),(EColor)(kTeal-6)};
const std::vector<std::string> ColourWheel::colournames{"kBlack", "kBlue", "kGreen", "kRed", "kViolet", "kOrange", "kMagenta","kAzure","kOrange","kViolet","kTeal"};

void plot_pmts_types_from_geofile(){
	TFile* f = TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_0.root");
	TTree* t= (TTree*)f->Get("wcsimGeoT");
	WCSimRootGeom* geo =0;
	t->SetBranchAddress("wcsimrootgeom",&geo);
	t->GetEntry(0);
	
	// make a map of polymarkers, which will mark the PMT locations
	// each marker collection will represent a different PMT type
	std::map<std::string, TPointSet3D*> pmt_polymarkers;
	
	bool useGeomNames =true; // plot using the names from the geom, or from the PMT objects
	
	// loop over the PMTs in the geometry and make a marker for each
	std::cout<<"geometry has "<<geo->GetWCNumPMT()<<" tank PMTs"<<std::endl;
	for(int pmti=0; pmti<geo->GetWCNumPMT(); pmti++){
		WCSimRootPMT thepmt = geo->GetPMT(pmti);
		std::string pmt_type_name;
		if(useGeomNames){
			int tubeid = thepmt.GetTubeNo();
			int tubetypeindex = geo->GetTubeIndex(tubeid);
			pmt_type_name = geo->GetWCPMTNameAt(tubetypeindex);
		} else {
			pmt_type_name = thepmt.GetName();
		}
		
		// consistency check: use with useGeomNames==true
		std::string pmt_type_name2 = thepmt.GetName();
		if(pmt_type_name2.find('_') != std::string::npos)
			pmt_type_name2.erase(pmt_type_name2.begin(), pmt_type_name2.begin()+pmt_type_name2.find('_')+1);
		if(pmt_type_name!=pmt_type_name2){
			std::cerr<<"Differing PMT Type names for PMT "<<pmti
							 <<" with tubeNo "<<thepmt.GetTubeNo()
							 <<", tube type index "<<geo->GetTubeIndex(thepmt.GetTubeNo())
							 <<", type name (from PMT object) "<<pmt_type_name2
							 <<", type name (from geometry type index) "
							 <<geo->GetWCPMTNameAt(geo->GetTubeIndex(thepmt.GetTubeNo()))
							 <<std::endl;
		}
		
		// make the marker collection if it doesn't exist
		if(pmt_polymarkers.count(pmt_type_name)==0){
			std::cout<<"new PMT type: "<<pmt_type_name<<std::endl;
			if(pmt_type_name=="INDEX_OUT_OF_BOUNDS"){
				std::cerr<<"PMT "<<pmti<<" with tubeNo "<<thepmt.GetTubeNo()
								 <<", tube type index "<<geo->GetTubeIndex(thepmt.GetTubeNo())
								 <<", type name (from PMT object) "<<thepmt.GetName()
								 <<", type name (from geometry type index) "
								 <<geo->GetWCPMTNameAt(geo->GetTubeIndex(thepmt.GetTubeNo()))
								 <<" and position ("
								 <<thepmt.GetPosition(0)<<", "
								 <<thepmt.GetPosition(1)<<", "
								 <<thepmt.GetPosition(2)<<")"
								 <<" has no type!"<<std::endl;
			}
			pmt_polymarkers.emplace(pmt_type_name, new TPointSet3D());
			pmt_polymarkers.at(pmt_type_name)->SetName(pmt_type_name.c_str());
		}
		// add this pmt to the marker collection
		float pmtx = thepmt.GetPosition(0);
		float pmty = thepmt.GetPosition(1);
		float pmtz = thepmt.GetPosition(2);
		//std::cout<<"adding a pmt of type "<<pmt_type_name<<std::endl;
		pmt_polymarkers.at(pmt_type_name)->SetNextPoint(pmtx,pmty, pmtz);
		
	}
	
	TCanvas thecanvas("thecanvas");
	// Import the gdml geometry for the detector:
	TGeoManager::Import("../root_work/annie_v04_aligned.gdml");
	// make all the materials semi-transparent so we can see the PMTs
	TList* matlist = gGeoManager->GetListOfMaterials();
	TIter nextmaterial(matlist);
	while(TGeoMixture* amaterial=(TGeoMixture*)nextmaterial()) amaterial->SetTransparency(90);
	gGeoManager->GetVolume("WORLD2_LV")->Draw("ogl");
	
//	// need to create axes to add to the TPolyMarker3D... maybe not if overlaying over gdml
//	TH3F *frame3d = new TH3F("frame3d","frame3d",10,-1.5,1.5,10,0,3.3,10,-2.5,2.5);
//	frame3d->Draw();
//	frame3d->SetTitle("PMT Positions and Types");
	
	// give all polymarker sets a unique colour and draw them
	ColourWheel thecolourwheel = ColourWheel();
	for(auto& amarkerset : pmt_polymarkers){
		std::cout<<"drawing markers for collection "<<amarkerset.first<<std::endl;
		TPointSet3D* themarkers = amarkerset.second;
		themarkers->SetMarkerStyle(20);
		themarkers->SetMarkerColor(thecolourwheel.GetNextColour());
		themarkers->Draw("same");
	}
	thecanvas.BuildLegend();
	thecanvas.Modified();
	thecanvas.Update();
	
	do{
		gSystem->ProcessEvents();
		// N.B. std::chrono::(milli)seconds must take INTEGER arguments!
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
	} while(gROOT->FindObject("thecanvas")!=nullptr);
}
