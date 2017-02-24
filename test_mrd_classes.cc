//#ifndef _PADDLE_PLACEMENT_VERBOSE_
//#define _PADDLE_PLACEMENT_VERBOSE_ 1
//#endif

#include "MRDStrikeClass.hh"
#include "MRDTrackClass.hh"
#include "MRDspecs.hh"
#include "TRandom.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TLine.h"

void ComputePaddleTransformation(const Int_t copyNo, TVector3 &origin, Bool_t &ishpaddle, Bool_t verbose);
void ComputeSteelTransformation (const Int_t copyNo, TVector3 &origin);

void generatesegfault(){

	TRandom RandomGenerator;
	
	// cMRDTrack takes in a vector of cMRDStrikes, so first we need to generate a bunch of strikes
	std::vector<cMRDStrike> strikesthistrack;
	std::vector<int> layeroffsets {0,30,56,86,112,142,168,198,224,254,280,310};
	std::vector<int> pmtids;
	//std::vector<int> pmtidsin {22, 2, 18, 4, 14, 6, 10, 8, 6, 10, 2};		// these are layer offsets to easily generate a track
	//std::vector<int> pmtidsin {22, 2, 20, 4, 18, 6, 16, 8, 14, 10, 12};
	std::vector<int> pmtidsin {22, 2, 22, 4, 20, 6, 20, 8, 18, 10};
	// note there are TWO PADDLES ON THE SAME LEVEL (panel consists of 2 halves) so add 2 to offset by each layer!
	Int_t acounter=0, ancounter=0;
	for(std::vector<int>::iterator it=pmtidsin.begin(); it!=pmtidsin.end(); ++it){
		//if((acounter%2)!=0){ 
			pmtids.push_back((*it)+layeroffsets.at(acounter)); 
			cout<<pmtids.at(ancounter)<<", ";
			ancounter++;
		//}
		acounter++;
	}
	cout<<"set pmts"<<endl;
	// push in these strikes for this track.
	for(int i=0; i<pmtids.size(); i++){
		// Each cMRDStrike takes in a vector of hit times and an int of PMT ID - all strikes in this track can be on the same PMT for now
		Int_t PMTID=pmtids.at(i);
		std::vector<Double_t> hittimesthisstrike;
		// generate some hits in this strike. poisson the number of hits up to 100.
		int numhits = RandomGenerator.Poisson(100);
		//cout<<"this strike hits pmt "<<PMTID<<" with "<<numhits<<" photons"<<endl;
		// generate the times for the hits and add them to the hit vector
		for(int j=0; j<numhits; j++){
			// uniform distribution of hit times over 30ns
			double hittime = RandomGenerator.Uniform(0.,30.);
			hittimesthisstrike.push_back(hittime);
		}
		// make the strike
		//cout<<"making the strike"<<endl;
		cMRDStrike astrike = cMRDStrike(hittimesthisstrike, PMTID);
		// push this strike to the vector of strikes for this track
		strikesthistrack.push_back((astrike));
	}
	
	cout<<"making track with "<<strikesthistrack.size()<<" strikes"<<endl;
	cMRDTrack* aTrack = new cMRDTrack(strikesthistrack);
	cout<<"track made"<<endl;
	std::vector<cMRDStrike> thestrikes = aTrack->GetStrikes();
	std::vector<int> thestruckpmts;
	for(std::vector<cMRDStrike>::iterator it=thestrikes.begin(); it!=thestrikes.end(); ++it){
		thestruckpmts.push_back(it->GetPMTNumber());
	}
	
	cout<<"making canvas plots"<<endl;
	Double_t scintboxwidth=1;
	Double_t canvw=700., canvh=700.;
	// build a root pad that represents the MRD and its paddles. highlight the paddles that were hit, and draw the reconstructed region.
	TCanvas *topview = new TCanvas("topview", "MRD Vis Top View",canvw,canvh);
	TCanvas *sideview  = new TCanvas("sideview", "MRD Vis Side View",canvw,canvh);
	//cout<<"adding scints"<<endl;
	Bool_t firsthpaddledone=false, firstvpaddledone=false;
	std::pair<double, double> xupcorner1, xupcorner2, xdowncorner1, xdowncorner2, yupcorner1, yupcorner2, ydowncorner1, ydowncorner2;
	Bool_t rightsidehit=false, leftsidehit=false, tophit=false, bottomhit=false;
	
	// loop over all paddles and build a map of which were hit
	for(int paddle=0; paddle<((numpaddlesperpanelv+numpaddlesperpanelh)*(numpanels/2)); paddle++){
	
		// check if the paddle was hit by finding it's PMT ID in the struck PMTs
		Bool_t paddleishit;
		if(std::find(thestruckpmts.begin(), thestruckpmts.end(), paddle)!=thestruckpmts.end()){paddleishit=true;} else {paddleishit=false;}
		
		// calculate the paddle position from it's ID
		TVector3 origin;
		Double_t holderx, holdery, holderz, holderzoffset;
		Bool_t ishpaddle;
		ComputePaddleTransformation(paddle, origin, ishpaddle, paddleishit);
		
		// establish if this hits a new half of the horizontal or vertical panel.
				 if(paddleishit&&ishpaddle&&origin.Y()>0){tophit=true;}
		else if(paddleishit&&ishpaddle&&origin.Y()<0){bottomhit=true;}
		else if(paddleishit&&(!ishpaddle)&&origin.X()>0){rightsidehit=true;}
		else if(paddleishit&&(!ishpaddle)&&origin.X()<0){leftsidehit=true;}
		
		// convert from simulation 'cm' units to canvas units: 0-1, 0-1 horizontal and vertical
		holderx = ((origin.X()/(maxwidth*1.2)))+0.5;
		Double_t anoffset=0, viewoffset=(scintfullzlen+scintalugap)*12;
		if((ishpaddle&&(origin.X()<0))||((!ishpaddle)&&(origin.Y()<0))){anoffset=(scintfullzlen+scintalugap)*5;}
		// in order to view the two sets in the same graph, shift one half on the canvas
		holdery = (origin.Y()/(maxheight*1.2))+0.5;
		holderz = ((origin.Z()+anoffset)/(mrdZlen*1.2))+0.5;
		/*cout<<"origin is: ("<<origin.X()<<","<<origin.Y()<<","<<origin.Z()<<")"<<endl;
		cout<<"scaled origin is: ("<<holderx<<","<<holdery<<","<<holderz<<")"<<endl;
		cout<<"ishpaddle="<<ishpaddle<<endl;
		if(paddle==0){
			cout<<"drawing box with height "<<(scintfullxlen/(maxwidth*1.2))<<" and width "<<(scintboxwidth/(mrdZlen*1.2))<<endl;
		}*/
		if((!ishpaddle)&&(!firsthpaddledone)&&paddleishit){
			firsthpaddledone=true;
			xupcorner1=std::pair<double,double>(holderx,holderz);
			xdowncorner1=std::pair<double,double>(holderx+(scintfullxlen/(maxwidth*1.3)),holderz+(scintboxwidth/mrdZlen));
			//cout<<"first hit h paddle at ("<<holderx<<","<<holdery<<","<<holderz<<")"<<endl;
		}
		if((!ishpaddle)&&paddleishit){
			xdowncorner2=std::pair<double,double>(holderx,holderz);
			xupcorner2=std::pair<double,double>(holderx+(scintfullxlen/(maxwidth*1.3)),holderz+(scintboxwidth/mrdZlen));
			//cout<<"last h paddle at ("<<holderx<<","<<holdery<<","<<holderz<<")"<<endl;
		}
		if(ishpaddle&&(!firstvpaddledone)&&paddleishit){
			firstvpaddledone=true;
			yupcorner1=std::pair<double,double>(holdery,holderz);
			ydowncorner1=std::pair<double,double>(holdery+(scintfullxlen/(maxwidth*1.2)),holderz+(scintboxwidth/mrdZlen));
			//cout<<"first hit v paddle at ("<<holderx<<","<<holdery<<","<<holderz<<")"<<endl;
		}
		if(ishpaddle&&paddleishit){
			ydowncorner2=std::pair<double,double>(holdery,holderz);
			yupcorner2=std::pair<double,double>(holdery+(scintfullxlen/(maxwidth*1.2)),holderz+(scintboxwidth/mrdZlen));
			//cout<<"last h paddle at ("<<holderx<<","<<holdery<<","<<holderz<<")"<<endl;
		}
		
		// Draw the box
		TBox *thepaddle;
		if(!ishpaddle){
			sideview->cd();
			//TBox (Double_t leftx, Double_t bottomy, Double_t rightx, Double_t topy)
			thepaddle = new TBox(holderz,holderx,holderz+(scintboxwidth/mrdZlen),holderx+(scintfullxlen/(maxwidth*1.3)));
		} else {
			topview->cd();
			thepaddle = new TBox(holderz,holdery,holderz+(scintboxwidth/mrdZlen),holdery+(scintfullxlen/(maxwidth*1.2)));
		}
		thepaddle->SetFillStyle(1001);	// solid fill
		if(paddleishit){ thepaddle->SetFillColor(kYellow); } else { thepaddle->SetFillColor(11); }
		// if the paddle was hit light it up - color 5 (yellow) otherwise leave it dark - colour 11 (grey)
		thepaddle->SetLineColor(1);	//black
		thepaddle->SetLineWidth(2);
		thepaddle->Draw();
		
		// check if we're changing panel
		if((std::find(layeroffsets.begin(), layeroffsets.end(), (paddle+1))!=layeroffsets.end())&&paddle!=0){
			Double_t holderx1, holderx2, holdery1, holdery2;
			Double_t otherviewoffset;
			// this is the last paddle of this layer - draw the layer in the opposite view, combining all paddles
			if(!ishpaddle){
				otherviewoffset=-0.014;
				topview->cd();
				holderx1=((scintalugap*10)/(maxwidth*1.2));
				holderx2=(scinthfullylen/(maxwidth*1.2));
				thepaddle = new TBox(holderz+otherviewoffset,0.5+holderx1,holderz+(scintboxwidth/mrdZlen)+otherviewoffset, 0.5+holderx1+holderx2);
				thepaddle->SetFillStyle(1001);	// solid fill
				if(leftsidehit){ thepaddle->SetFillColor(kYellow); } else { thepaddle->SetFillColor(11); }	// yellow if hit or grey otherwise
				thepaddle->SetLineColor(1);	//black
				thepaddle->SetLineWidth(2);
				thepaddle->Draw();
				
				thepaddle = new TBox(holderz+otherviewoffset,0.5-holderx1,holderz+(scintboxwidth/mrdZlen)+otherviewoffset, 0.5-(holderx1+holderx2));
				thepaddle->SetFillStyle(1001);	// solid fill
				if(rightsidehit){ thepaddle->SetFillColor(kYellow); } else { thepaddle->SetFillColor(11); }	// yellow if hit or grey otherwise
				thepaddle->SetLineColor(1);	//black
				thepaddle->SetLineWidth(2);
				thepaddle->Draw();
				
			} else {
				otherviewoffset=-0.01;
				sideview->cd();
				holdery1=((scintalugap*10)/(maxheight*1.2));
				holdery2=(scintvfullylen/(maxheight*1.2));
				thepaddle = new TBox(holderz+otherviewoffset,holdery1+0.5,holderz+(scintboxwidth/mrdZlen)+otherviewoffset,0.5+holdery1+holdery2);
				thepaddle->SetFillStyle(1001);	// solid fill
				if(bottomhit){ thepaddle->SetFillColor(kYellow); } else { thepaddle->SetFillColor(11); }	// yellow if hit or grey otherwise
				thepaddle->SetLineColor(1);	//black
				thepaddle->SetLineWidth(2);
				thepaddle->Draw();
				
				thepaddle = new TBox(holderz+otherviewoffset,0.5+-holdery1,holderz+(scintboxwidth/mrdZlen)+otherviewoffset,0.5-(holdery1+holdery2));
				thepaddle->SetFillStyle(1001);	// solid fill
				if(tophit){ thepaddle->SetFillColor(kYellow); } else { thepaddle->SetFillColor(11); }	// yellow if hit or grey otherwise
				thepaddle->SetLineColor(1);	//black
				thepaddle->SetLineWidth(2);
				thepaddle->Draw();
			}
			
			// reset for the next layer
			tophit=false; bottomhit=false; leftsidehit=false; rightsidehit=false;
		}
		
	}
	//TODO: rather than projecting to 0.8, project to the appropriate further depth of penetration.
	// add trajectory lines 
	//=====================
	// TLine::TLine(Double_t x1, Double_t y1, Double_t x2, Double_t y2))
	std::map<const char*, double> xlines = aTrack->GetTrackXStats();
	std::map<const char*, double> ylines = aTrack->GetTrackYStats();
	// has keys "angminxval", "angminzval", "angmin", same for angmax
	Double_t holderx1, holdery1, holderz1, holderx2, holdery2, holderz2, projectedleft, projectedright;
	TLine *xupgoing, *xdowngoing, *yupgoing, *ydowngoing;
	Bool_t usebasic=false, useclass=true;
	
	// first upgoing x line
	cout<<"making basic upgoing x line"<<endl;
	if(usebasic){	// simple reconstruction: just use corners of first and last paddles hit
		holderx1=xupcorner1.first; holderz1=xupcorner1.second;
		holderx2=xupcorner2.first; holderz2=xupcorner2.second;
		projectedleft = holderx1-holderz1*((holderx2-holderx1)/(holderz2-holderz1));
		//xupgoing = new TLine(holderz1,holderx1,holderz2, holderx2);	//without projection, just join corners
		xupgoing = new TLine(0.,projectedleft,holderz2, holderx2);		// with projection to left
		//-------------------------
		xupgoing->SetLineColor(kRed);
		xupgoing->SetLineWidth(4);
		xupgoing->SetLineStyle(3);
		sideview->cd();
		xupgoing->Draw();
	} 
	cout<<"making class upgoing x line"<<endl;
	if (useclass){			// track class has a more thorough method that ensures all paddles are hit
		holderx1=((xlines.at("angminxval"))/(maxwidth*1.2))+0.5;
		holderz1=((xlines.at("angminzval"))/(mrdZlen*1.2))+0.5;
		projectedleft=holderx1-(holderz1*(TMath::Tan(xlines.at("angmin")))*(mrdZlen/maxwidth));
		projectedright=holderx1+((0.8-holderz1)*TMath::Tan(xlines.at("angmin"))*(mrdZlen/maxwidth));
		xupgoing = new TLine(0.02,projectedleft,0.8,projectedright);
				/*cout<<"xlines.at('angminxval')="<<xlines.at("angminxval")<<endl;
				cout<<"xlines.at('angminzval')="<<xlines.at("angminzval")<<endl;
				cout<<"xlines.at('angmin')="<<TMath::RadToDeg()*xlines.at("angmin")<<endl;
				cout<<"holderx1="<<holderx1<<" projected left by z="<<xlines.at("angminzval")<<"="<<holderz1<<" at angle "
				<<TMath::RadToDeg()*((mrdZlen/maxwidth)*xlines.at("angmin"))<<" giving offset of "
				<<(holderz1*TMath::Tan(xlines.at("angmin"))*(mrdZlen/maxwidth))
				<<", projectedleft="<<projectedleft<<endl
				<<"holderz1="<<holderz1
				<<", (TMath::Tan(xlines.at('angmin'))="<<((maxwidth/mrdZlen)*TMath::Tan(xlines.at("angmin")))
				<<", projectedright="<<projectedright<<endl;*/
		//--------------------------
		xupgoing->SetLineWidth(4);
		xupgoing->SetLineColor(kBlue);
		xupgoing->SetLineStyle(3);
		sideview->cd();
		xupgoing->Draw();
	}
	
	cout<<"making basic downgoing x line"<<endl;
	// now downgoing x line
	if(usebasic){
		holderx1=xdowncorner1.first; holderz1=xdowncorner1.second;
		holderx2=xdowncorner2.first; holderz2=xdowncorner2.second;
		projectedleft = holderx1+holderz1*((holderx1-holderx2)/(holderz2-holderz1));
		xdowngoing = new TLine(0.,projectedleft,holderz2, holderx2);
		//--------------------------
		xdowngoing->SetLineColor(kRed);
		xdowngoing->SetLineWidth(4);
		xdowngoing->SetLineStyle(3);
		sideview->cd();
		xdowngoing->Draw();
	}
	cout<<"making class downgoing x line"<<endl;
	if(useclass){
		holderx1=((xlines.at("angmaxxval"))/(maxwidth*1.2))+0.5;
		holderz1=((xlines.at("angmaxzval"))/(mrdZlen*1.2))+0.5;
		projectedleft=holderx1+(holderz1*(TMath::Tan(xlines.at("angmax")))*(mrdZlen/maxwidth));
		projectedright=holderx1-((0.8-holderz1)*TMath::Tan(xlines.at("angmax"))*(mrdZlen/maxwidth));
		xdowngoing = new TLine(0.02,projectedleft,0.8,projectedright);
				/*cout<<"xlines.at('angmaxxval')="<<xlines.at("angmaxxval")<<endl;
				cout<<"xlines.at('angmaxzval')="<<xlines.at("angmaxzval")<<endl;
				cout<<"xlines.at('angmax')="<<TMath::RadToDeg()*xlines.at("angmax")<<endl;
				cout<<"holderx1="<<holderx1<<" projected left by z="<<xlines.at("angmaxzval")<<"="<<holderz1<<" at angle "
				<<TMath::RadToDeg()*((mrdZlen/maxwidth)*xlines.at("angmax"))<<" giving offset of "
				<<(holderz1*TMath::Tan(xlines.at("angmax"))*(mrdZlen/maxwidth))
				<<", projectedleft="<<projectedleft<<endl
				<<"holderz1="<<holderz1
				<<", (TMath::Tan(xlines.at('angmax'))="<<((maxwidth/mrdZlen)*TMath::Tan(xlines.at("angmax")))
				<<", projectedright="<<projectedright<<endl;*/
		//--------------------------
		xdowngoing->SetLineColor(kBlue);
		xdowngoing->SetLineWidth(4);
		xdowngoing->SetLineStyle(3);
		sideview->cd();
		xdowngoing->Draw();
	}

	cout<<"making basic upgoing y line"<<endl;
	// now upgoing y line
	if(usebasic){
		holdery1=yupcorner1.first; holderz1=yupcorner1.second;
		holdery2=yupcorner2.first; holderz2=yupcorner2.second;
		projectedleft = holdery1-holderz1*((holdery2-holdery1)/(holderz2-holderz1));
		yupgoing = new TLine(0.,projectedleft,holderz2, holdery2);
		//--------------------------
		yupgoing->SetLineWidth(4);
		yupgoing->SetLineColor(kRed);
		yupgoing->SetLineStyle(3);
		topview->cd();
		yupgoing->Draw();
	}
	cout<<"making class upgoing y line"<<endl;
	if(useclass){
		holdery1=((ylines.at("angminyval"))/(maxheight*1.2))+0.5;
		holderz1=((ylines.at("angminzval"))/(mrdZlen*1.2))+0.5;
		projectedleft=holdery1-(holderz1*(TMath::Tan(ylines.at("angmin")))*(mrdZlen/maxheight));
		projectedright=holdery1+((0.8-holderz1)*TMath::Tan(ylines.at("angmin"))*(mrdZlen/maxheight));
		yupgoing = new TLine(0.02,projectedleft,0.8,projectedright);
		//yupgoing = new TLine(holderz1,holdery1,0.8,projectedright);
				/*cout<<"ylines.at('angminyval')="<<ylines.at("angminyval")<<endl;
				cout<<"ylines.at('angminzval')="<<ylines.at("angminzval")<<endl;
				cout<<"ylines.at('angmin')="<<TMath::RadToDeg()*ylines.at("angmin")<<endl;
				cout<<"holdery1="<<holdery1<<" projected left by z="<<ylines.at("angminzval")<<"="<<holderz1<<" at angle "
				<<TMath::RadToDeg()*((mrdZlen/maxwidth)*ylines.at("angmin"))<<" giving offset of "
				<<(holderz1*TMath::Tan(ylines.at("angmin"))*(mrdZlen/maxwidth))
				<<", projectedleft="<<projectedleft<<endl
				<<"holderz1="<<holderz1
				<<", (TMath::Tan(ylines.at('angmin'))="<<((maxwidth/mrdZlen)*TMath::Tan(ylines.at("angmin")))
				<<", projectedright="<<projectedright<<endl;*/
		//--------------------------
		yupgoing->SetLineColor(kBlue);
		yupgoing->SetLineWidth(4);
		yupgoing->SetLineStyle(3);
		topview->cd();
		yupgoing->Draw();
	}
	
	// finally downgoing y line
	cout<<"making basic downgoing y line"<<endl;
	if(usebasic){
		holdery1=ydowncorner1.first; holderz1=ydowncorner1.second;
		holdery2=ydowncorner2.first; holderz2=ydowncorner2.second;
		projectedleft = holdery1+holderz1*((holdery1-holdery2)/(holderz2-holderz1));
		ydowngoing = new TLine(0.,projectedleft,holderz2, holdery2);
		//--------------------------
		ydowngoing->SetLineWidth(4);
		ydowngoing->SetLineColor(kRed);
		ydowngoing->SetLineStyle(3);
		topview->cd();
		ydowngoing->Draw();
	}
	cout<<"making class downgoing y line"<<endl;
	if(useclass){
		holdery1=((ylines.at("angmaxyval"))/(maxheight*1.2))+0.5;
		holderz1=((ylines.at("angmaxzval"))/(mrdZlen*1.2))+0.5;
		projectedleft=holdery1+(holderz1*(TMath::Tan(ylines.at("angmax")))*(mrdZlen/maxheight));
		projectedright=holdery1-((0.8-holderz1)*TMath::Tan(ylines.at("angmax"))*(mrdZlen/maxheight));
		ydowngoing = new TLine(0.02,projectedleft,0.8,projectedright);
		//ydowngoing = new TLine(holderz1,holdery1,0.8,projectedright);
		//--------------------------
		ydowngoing->SetLineColor(kBlue);
		ydowngoing->SetLineWidth(4);
		ydowngoing->SetLineStyle(3);
		topview->cd();
		ydowngoing->Draw();
	}
	
	cout<<"done"<<endl;
		

}



/*(	// informational only
 	static G4bool doonce = true;
  	if(doonce){
		scintzedges.push_back(alufullzlen+scintalugap);									// first panel is offset by depth of one alu struct. 
		scintzedges.push_back(alufullzlen+scintalugap+scintfullzlen);		// rear face fo first panel
		G4cout<<"layer one at z="<<(scintzedges.at(0)+tankouterRadius+2)<<" to "<<(scintzedges.at(1)+tankouterRadius+2)<<G4endl;
		// all others are at even intervals thereafter.
		layerthickness = steelfullzlen + alufullzlen + scintfullzlen + layergap;
		for(G4int i=1;i<numpanels;i++){
			scintzedges.push_back(i*layerthickness);		// front face
			scintzedges.push_back(i*layerthickness+scintfullzlen);	// rear face
			G4cout<<"layer "<<i<<" at z="<<(scintzedges.at(i*2)+tankouterRadius+2)<<" to "<<(scintzedges.at((i*2)+1)+tankouterRadius+2)<<G4endl;
		}
		doonce=false;
	}
*/

// code taken from MRDDetectorConstruction.cc, WCSimDetectorConstruction::ComputePaddleTransformation
//============================= 
void ComputePaddleTransformation(const Int_t copyNo, TVector3 &origin, Bool_t &ishpaddle, Bool_t verbose){
		Double_t Xposition=0, Yposition=0, Zposition=0;
		Int_t panelpairnum = TMath::Floor(copyNo/(numpaddlesperpanelv+numpaddlesperpanelh));	// which pair of panels
		Int_t panelnumrem = copyNo - panelpairnum*(numpaddlesperpanelv+numpaddlesperpanelh);	// copy num within a pair of panels
#ifdef _PADDLE_PLACEMENT_VERBOSE_
if(verbose){
		cout<<"computing transformation for PMT "<<copyNo<<", panelpairnum="<<panelpairnum<<", panelnumrem="<<panelnumrem<<endl;
}
#endif
		Int_t panelnum;
		Int_t paddlenum;
		//Bool_t ishpaddle=false;	// passed as arg in this 
		if(panelnumrem>(numpaddlesperpanelv-1)){					// first layer is a VERTICAL layer (leading horizontal layer removed)
			panelnum = (panelpairnum*2) +1;
			ishpaddle = true;
			paddlenum = panelnumrem%numpaddlesperpanelv; //copyNo%numpaddlesperpanelh;
		} else {
			panelnum = (panelpairnum*2);
			ishpaddle = false;
			paddlenum = panelnumrem; //copyNo%numpaddlesperpanelv;
		}
		Int_t pairnum = floor(paddlenum/2);					// Paddles 0&1 are a pair, then Y offset is the same for every pair
		Zposition = panelnum*(steelfullzlen + alufullzlen + scintfullzlen + layergap);	// fixed layer separation (except first)
		if(panelnum==0){Zposition = Zposition + alufullzlen + scintalugap;}							// layer 0 scints intrude into layer 1's offset
		Zposition = Zposition + steelfullzlen + steelscintgap;								// scint follows closely behind first steel
		Zposition = Zposition + (scintfullzlen/2);														// offset by half depth so we place front face not centre
		Zposition = Zposition - (mrdZlen/2);																	// offset by half total length to shift to front
		
#ifdef _PADDLE_PLACEMENT_VERBOSE_
if(verbose){
		cout<<"panel num "<<panelnum<<", paddle num "<<paddlenum;
		if(ishpaddle){cout<<", this is a horizontal paddle"<<endl;} else {cout<<", this is a vertical paddle"<<endl;}
//		if(paddlenum==0){cout <<"paddle placed at z = " 
//														<< (Zposition + (mrdZlen/2) - (scintfullzlen/2) + tankouterRadius + 2)
//														<< " cm to "<< (Zposition + (mrdZlen/2) + (scintfullzlen/2) + tankouterRadius + 2)
//														<< " cm." << endl;} 
}
#endif
		// Y position is offset by half length, so that one end is at the origin. This needs to be the correct half-length
		// for the appropriate paddle (H or V type).
		if (ishpaddle){
			if (paddlenum%2==0){
				Xposition=((scinthfullylen+scintbordergap)/2); 										// offset by +half length so one end is at x=0
			} else {
				Xposition=-((scinthfullylen+scintbordergap)/2);										// offset by -half length so one end is at x=0
			}
			Yposition = pairnum*(scintfullxlen+scintbordergap); 								// individual offset by pair number
			Yposition = Yposition - 0.5*(((scintfullxlen+scintbordergap)/2)*numpaddlesperpanelh);//-(scintfullxlen/2);
			// shift whole set by 1/2 total X extent to shift center back to X=0: HalfLength cancels doubed num of paddles
		} else {	// vertical panel
			if (paddlenum%2==0){
				Yposition=((scintvfullylen+scintbordergap)/2); 
			} else {
				Yposition=-((scintvfullylen+scintbordergap)/2);
			}
			Xposition = pairnum*(scintfullxlen+scintbordergap); 	// individual offset by pair number
			Xposition = Xposition - 0.5*(((scintfullxlen+scintbordergap)/2)*numpaddlesperpanelv);//-(scintfullxlen/2); 
			// shift whole set by 1/2 total Y extent to shift center back to Y=0: HalfLength cancels doubed num of paddles
		}
		
		TVector3 position(Xposition,Yposition,Zposition);
#ifdef _PADDLE_PLACEMENT_VERBOSE_
if(verbose){
		//cout<<"paddle positioned at ("<<Xposition<<","<<Yposition<<","<<Zposition<<")"<<endl;
		cout<<"paddle positioned at x="<<Xposition<<" to "<<Xposition+scintfullxlen;
		if(ishpaddle){cout<<"; y="<<Yposition<<" to "<<Yposition+scinthfullylen;}
		else {cout<<"; y="<<Yposition<<" to "<<Yposition+scintvfullylen;}
		cout<<"; z="<<Zposition<<" to "<<Zposition+scintfullzlen<<endl;
}
#endif
		origin=position;
		
}

void ComputeSteelTransformation (const Int_t copyNo, TVector3 &origin) {
		Double_t Xposition=0, Yposition=0, Zposition=0;
		Zposition=(copyNo)*(steelfullzlen + scintfullzlen + alufullzlen + layergap);	// layer width offset is always constant
		Zposition=Zposition + (steelfullzlen/2);													// offset by half depth so we are placing front face not centre
		Zposition=Zposition - (mrdZlen/2);																// offset by half total length to shift to front.
		
		TVector3 theposition(Xposition,Yposition,Zposition);
		origin=theposition;
	}
