
#ifndef _MRD_SPECS_
#define _MRD_SPECS_ 1
// THIS NEEDS UPDATING
Double_t scintzedges[24] = {19, 25, 75, 81, 150, 156, 225, 231, 300, 306, 375, 381, 450, 456, 525, 531, 600, 606, 675, 681, 750, 756, 825, 831};// taken from printed output generated by SteppingAction.cc at line 122
// ==================

Int_t numpanels=12;													// scintillator panels
Int_t	numplates=12;													// steel plates
Int_t numpaddlesperpanelh=26;								// paddles per h scintillator panel
Int_t numpaddlesperpanelv=30;								// paddles per v scintillator panel
Int_t numhpanels=6;													// verified
Int_t numvpanels=5;													// check?
Int_t numalustructs=13;											// number of supporting structs. We may be dropping one as we have fewer scintillators?
Int_t numvetopaddles=26;										// number of scintillator paddles in the FACC; 13 panels in 2 layers
Int_t vetopaddlesperpanel=13;								// number of scintillator paddles in each FACC panel

Double_t steelfullxlen = 305;
Double_t steelfullylen = 274;
Double_t steelfullzlen = 5;

Double_t scintfullxlen = 20;
Double_t scintfullzlen= 0.6;
Double_t scinthfullylen = 147.2; //155cm - 7.8cm tapered section		147.1 according to sciboone gdml export		//swapped???
Double_t scintvfullylen= 130.2;  //138cm - 7.8cm tapered section		129.4 according to sciboone gdml export

Double_t scinttapfullwidth = 17.1; 			// width of the tapering part of scintillator paddles at the narrow end
Double_t scinttapfullheight = 7.8; 			// z length of tapering part of scint paddles.

Double_t scintlgfullwidth = 5.08; 			// tapered light guides at narrow end
Double_t scintlgfullheight = 33.3; 			// 

Double_t alufullxlen = steelfullxlen+15;	
Double_t alufullylen = steelfullylen+15;	
Double_t alufullzlen = 3.81;							// from eye and a skim of the mrdmodule.txt file, i'm guessing depth is ~0.75 inches (1.9cm)
Double_t alufullxthickness = 2.54;					// as above, guessing frame to be 1 inch box cross-section
Double_t alufullythickness = 2.54;
Double_t windowwidth = (steelfullxlen-(4*alufullxthickness))/3;
Double_t windowheight= (steelfullylen-(4*alufullythickness))/3;

Double_t scintbordergap=0.3;								// gap between each scintillator (cm) (border to account for cladding etc)
Double_t steelscintgap=0.5;									// gap between steel and scintillator
Double_t scintalugap=0.2;										// gap between scintillator and alu struct
Double_t alusteelgap=2.0; 									// gap between alu struct and subsequent steel of next layer
Double_t layergap = steelscintgap + scintalugap + alusteelgap;	// total gaps of layers
Double_t nothickness = 0.01;

Double_t mrdZlen = numplates*steelfullzlen + (numpanels+1)*scintfullzlen + numalustructs*alufullzlen + numpanels*layergap + scintalugap; 

Double_t tankouterRadius= 152.4;

// calculate full extent of MRD for it's mother volume / the canvas scaling factor.
Double_t MRDPMTExposeHeight=0;
Double_t mrdpmtfullheight = MRDPMTExposeHeight;
Double_t  widths[] = {2*(scinthfullylen+scinttapfullheight+scintlgfullheight+(scintbordergap/2)+mrdpmtfullheight+nothickness),((numpaddlesperpanelv/2)*(scintfullxlen+scintbordergap))};
Double_t  heights[] = {2*(scintvfullylen+scinttapfullheight+scintlgfullheight+(scintbordergap/2)+mrdpmtfullheight+nothickness),((numpaddlesperpanelh/2)*(scintfullxlen+scintbordergap))};
Double_t maxwidth = *std::max_element(widths,widths+(sizeof(widths)/sizeof(widths[0])))+0.1;
Double_t maxheight = *std::max_element(heights,heights+(sizeof(heights)/sizeof(heights[0])))+0.1;

//	totMRD_box = new G4Box("totMRD",(maxwidth/2),(maxheight/2),mrdZlen/2);

#endif 
