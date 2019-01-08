#ifndef _MRD_SPECS_
#define _MRD_SPECS_ 1

// first panel is a HORIZONTAL LAYER.
// TODO: this should probably be a namespace
// or better, stored in the WCSimRootGeom

namespace MRDSpecs{
const Int_t numpaddlesperpanelv=30;
const Int_t numpaddlesperpanelh=26;
const Int_t numhpanels=6;
const Int_t numvpanels=5;
const Int_t numpanels=numhpanels+numvpanels;
const Int_t numplates=12;

const Double_t scintfullxlen = 20;
const Double_t scintfullzlen = 0.6;
const Double_t scinthfullylen = 147.2;
const Double_t scintvfullylen = 130.2;

const Int_t numalustructs=11;
const Int_t numvetopaddles=26;
const Int_t vetopaddlesperpanel=13;

const Double_t steelfullxlen = 305;
const Double_t steelfullylen = 274;
const Double_t steelfullzlen = 5;

const Double_t scinttapfullwidth = 17.1; 	// width of the scintilling taper at the narrow end
const Double_t scinttapfullheight = 7.8; 	// z length of tapering part of scint paddles.

const Double_t scintlgfullwidth = 5.08; 	// tapered light guides at narrow end
const Double_t scintlgfullheight = 33.3; 

const Double_t alufullxlen = steelfullxlen+15;
const Double_t alufullylen = steelfullylen+15;
const Double_t alufullzlen = 3.81;		// TODO: check this with measurement
const Double_t alufullxthickness = 2.54;	// TODO: check this with measurement
const Double_t alufullythickness = 2.54;	// TODO:  " "
const Double_t windowwidth = (steelfullxlen-(4*alufullxthickness))/3;
const Double_t windowheight= (steelfullylen-(4*alufullythickness))/3;

const Double_t scintbordergap=0.3;		// gap between each scintillator (cm) to account for cladding etc
const Double_t steelscintgap=0.5;		// gap between steel and scintillator
const Double_t scintalugap=0.2;			// gap between scintillator and alu struct
const Double_t alusteelgap=2.0; 		// gap between alu struct and subsequent steel of next layer
const Double_t layergap = steelscintgap + scintalugap + alusteelgap;	// total gaps between layers
const Double_t nothickness = 0.01;

const Double_t mrdZlen = numplates*steelfullzlen + (numpanels+1)*scintfullzlen + numalustructs*alufullzlen + numpanels*layergap + scintalugap; // TODO reconcile this with MRD_depth

const Double_t tankouterRadius= 152.4;

// calculate full extent of MRD for it's mother volume / the canvas scaling factor.
const Double_t MRDPMTExposeHeight=0;
const Double_t MRDPMTRadius=5.08;
const Double_t mrdpmtfullheight = MRDPMTExposeHeight;
const Double_t  widths[] = {2*(scinthfullylen+scinttapfullheight+scintlgfullheight+(scintbordergap/2)+mrdpmtfullheight+nothickness),((numpaddlesperpanelv/2)*(scintfullxlen+scintbordergap))};
const Double_t  heights[] = {2*(scintvfullylen+scinttapfullheight+scintlgfullheight+(scintbordergap/2)+mrdpmtfullheight+nothickness),((numpaddlesperpanelh/2)*(scintfullxlen+scintbordergap))};
const Double_t maxwidth = *std::max_element(widths,widths+(sizeof(widths)/sizeof(widths[0])))+0.1;
const Double_t maxheight = *std::max_element(heights,heights+(sizeof(heights)/sizeof(heights[0])))+0.1;

const Int_t nummrdpmts=306;
const std::vector<int> layeroffsets {0, 26, 56, 82, 112, 138, 168, 194, 224, 250, 280, 306};
// ids of the first pmt in each layer.                                   KEEP an extra ^^^

const Float_t MRD_width = ((numpaddlesperpanelv/2)*(scintfullxlen+scintbordergap))/2.;
const Float_t MRD_height = ((numpaddlesperpanelh/2)*(scintfullxlen+scintbordergap))/2.;
const Float_t MRD_layer2 = 290.755;             // position in wcsim coords of second scint layer in cm
const Float_t MRD_start = 325.5;                // position in wcsim coord of MRD front face in cm
const Float_t MRD_depth = 139.09;               // total depth of the MRD in cm
const Float_t MRD_end = 464.59;                 // end of the MRD in cm
const Float_t MRD_steel_width = (305./2.);      // half width of steel in cm
const Float_t MRD_steel_height = (274./2.);     // half height of steel in cm
/* output from WCSim:
########## MRD front face: 325.5                     ##########
########## MRD total Z length: 139.09                ##########
########## MRD scintillator layer 0  (H) at z=336.08 ##########     1
########## MRD scintillator layer 1  (V) at z=348.19 ########## 1
########## MRD scintillator layer 2  (H) at z=360.30 ##########     2
########## MRD scintillator layer 3  (V) at z=372.41 ########## 2
########## MRD scintillator layer 4  (H) at z=384.52 ##########     3
########## MRD scintillator layer 5  (V) at z=396.63 ########## 3
########## MRD scintillator layer 6  (H) at z=408.74 ##########     4
########## MRD scintillator layer 7  (V) at z=420.85 ########## 4
########## MRD scintillator layer 8  (H) at z=432.96 ##########     5
########## MRD scintillator layer 9  (V) at z=445.07 ########## 5
########## MRD scintillator layer 10 (H) at z=457.18 ##########
*/
std::vector<double> mrdscintlayers{336.080, 348.190, 360.300, 372.410, 384.520, 396.630, 408.740, 420.850, 432.960, 445.070, 457.180 };

//TODO: should retrieve this info from the geo in wcsimanalysis class
const Float_t tank_start = 15.70;           // front face of the tank in cm
const Float_t tank_radius = 152.4;          // tank radius in cm
const Float_t tank_halfheight = 198.;       // tank half height in cm
const Float_t tank_yoffset = -14.46;        // tank y offset in cm

//totMRD_box = new G4Box("totMRD",(maxwidth/2),(maxheight/2),mrdZlen/2);
}
#endif 
