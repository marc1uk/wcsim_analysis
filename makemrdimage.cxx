/* vim:set noexpandtab tabstop=4 wrap */
#include "TRandom.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TLine.h"
#include "TArrow.h"

#ifndef DRAWVERBOSE
//#define DRAWVERBOSE 1
#endif

#ifndef DRAWSUPERVERBOSE
//#define DRAWSUPERVERBOSE 1
#endif

void cMRDSubEvent::DrawMrdCanvases(){
#ifdef DRAWVERBOSE
	cout<<"making canvas plots"<<endl;
#endif
	if(fillstaticmembers) FillStaticMembers();  // fills the canvas colours
	//debug: highlight the corner paddles for checking
//	std::vector<int> holdme(pmts_hit); // remember to disable this also at the end of this method
//	std::vector<int> pmtshit_debug{
//	0,1,
//	(MRDSpecs::numpaddlesperpanelh/2)-1,(MRDSpecs::numpaddlesperpanelh/2),
//	MRDSpecs::numpaddlesperpanelh-1,MRDSpecs::numpaddlesperpanelh-2,
//	MRDSpecs::nummrdpmts-MRDSpecs::numpaddlesperpanelh,MRDSpecs::nummrdpmts-MRDSpecs::numpaddlesperpanelh+1,
//	MRDSpecs::nummrdpmts-(MRDSpecs::numpaddlesperpanelh/2)-1,MRDSpecs::nummrdpmts-(MRDSpecs::numpaddlesperpanelh/2),
//	MRDSpecs::nummrdpmts-2, MRDSpecs::nummrdpmts-1};
// test set 1:
//	0,2,4,6,8,10,12,14,16,18,20,22,24};  // all paddles in the RIGHT HAND (x>0) side of the first H layer.
// test set 2:
//	for(auto&& apmt : pmtshit_debug) apmt++; // all left hand (x<0) paddles of first H layer. 
// test set 3:
//	0,2,4,6,8,10,12,14,16,18,20,22,24,26,28};  // in conjunction with below:
//	for(auto&& apmt : pmtshit_debug) apmt+=MRDSpecs::numpaddlesperpanelh; // all TOP (y>0) paddles of first V layer.
//// test set 4:
//	for(auto&& apmt : pmtshit_debug) apmt++;  // all BOTTOM (y<0) paddles of first V layer.
	
//	pmts_hit = pmtshit_debug; 
//	std::vector<double> holdmetoo(digi_ts);
//	std::vector<double> digi_tsdebug
//	{
//	0,220,
//	220,220,
//	220,220,
//	220,220,
//	220,220,
//	220,220};
//	(pmtshit_debug.size(),100.);
//	digi_ts=digi_tsdebug;
	
	Double_t scintboxwidth=1;
	Double_t topboxheight = (MRDSpecs::scintfullxlen/(MRDSpecs::maxwidth*1.2));
	Double_t sideboxheight = (MRDSpecs::scintfullxlen/(MRDSpecs::maxwidth*1.3));
	Int_t paddlecolour;
	Double_t maxtime = (*std::max_element(digi_ts.begin(), digi_ts.end()));
	Double_t mintime = (*std::min_element(digi_ts.begin(), digi_ts.end()));
	static TText *texttop;
	static TText *textbottom;
#ifdef DRAWVERBOSE
	cout<<"min time is "<<mintime<<endl;
	cout<<"max time is "<<maxtime+mintime<<endl;
#endif
	
	// build a root pad that represents the MRD and its paddles. 
	// highlight the paddles that were hit, and draw the reconstructed region.
	Double_t canvw=1200., canvh=700.;
	if(imgcanvas==0){
#ifdef DRAWVERBOSE
	cout<<"making canvases"<<endl;
#endif
		imgcanvas = new TCanvas("imgcanvas","MRD Digit Visualiser",canvw,canvh); 
		imgcanvas->Divide(2,1);
		imgcanvas->cd(1);
		titleleft = new TText(.32,.9,TString::Format("Side View, Event %d",event_id));
		titleleft->Draw();
		imgcanvas->cd(2);
		titleright = new TText(.32,.9,TString::Format("Top View, Event %d",event_id));
		titleright->Draw();
		imgcanvas->cd(2);
		for(int i=0; i<aspectrumv.size(); i++){
			TBox* colourbox = new TBox(0.02,0.1+(i*0.04),0.05,0.1+(i*0.04)+0.03);
			colourbox->SetFillColor(aspectrumv.at(i));
			colourbox->Draw();
		}
		char titlebuf[50];
		snprintf(titlebuf,50,"%d",(int)(maxtime*2));
		texttop = new TText(0.0,0.88,titlebuf);
		texttop->SetTextSize(0.04);
		texttop->Draw();
		snprintf(titlebuf,50,"%d",(int)mintime);
		textbottom = new TText(0.0,0.05,titlebuf);
		textbottom->SetTextSize(0.04);
		textbottom->Draw();
	} else {
#ifdef DRAWVERBOSE
		cout<<"re-drawing titles and time labels on canvases"<<endl;
#endif
		imgcanvas->cd(1);
//		gPad->Clear();
//		titleleft->Draw();
		titleleft->SetTitle(TString::Format("Side View, Event %d",event_id));
		titleleft->Draw();
		imgcanvas->cd(2);
//		gPad->Clear();
//		titleright->Draw();
		titleright->SetTitle(TString::Format("Top View, Event %d",event_id));
		titleright->Draw();
		char titlebuf[50];
		snprintf(titlebuf,50,"%d",(int)(maxtime*2));
		texttop->SetTitle(titlebuf);
		texttop->Draw();
		snprintf(titlebuf,50,"%d",(int)mintime);
		textbottom->SetTitle(titlebuf);
		textbottom->Draw();
	}

#ifdef DRAWVERBOSE
	cout<<"adding scints"<<endl;
#endif
	Bool_t firsthpaddledone=false, firstvpaddledone=false;
	//std::pair<double, double> xupcorner1, xupcorner2, xdowncorner1, xdowncorner2, yupcorner1, yupcorner2, ydowncorner1, ydowncorner2;	// part of cMRDSubEvent class now
	Bool_t rightsidehit=false, leftsidehit=false, tophit=false, bottomhit=false;
	
	// loop over all paddles and build a map of which were hit
	for(int paddle=0; paddle<((MRDSpecs::numpaddlesperpanelv*MRDSpecs::numvpanels)+(MRDSpecs::numpaddlesperpanelh*MRDSpecs::numhpanels)); paddle++){
	
		// check if the paddle was hit by finding it's PMT ID in the struck PMTs
		Bool_t paddleishit;
		Double_t thetime;
		std::vector<Int_t>::iterator theit = std::find(pmts_hit.begin(), pmts_hit.end(), paddle);
		if(theit!=pmts_hit.end()){
			paddleishit=true;
			// grab the time the digit on this PMT. (what to do with more than one?)
			Int_t theindex = std::distance(pmts_hit.begin(),theit);
			thetime = digi_ts.at(theindex);
			Double_t relatime = (maxtime!=mintime) ? (thetime-mintime)/(maxtime-mintime) : 0.;
			Int_t colorindex = TMath::Floor((aspectrumv.size()-1)*(relatime)); //relatime/2.
			paddlecolour = aspectrumv.at(colorindex);
			// kRed is an EColor, kRed+2 is an Int_t, representing generically a TColor, TBox requires a Color_t!
		} else {
			paddleishit=false;
		}
		
		// calculate the paddle position from it's ID
		TVector3 origin;
		Double_t holderx, holdery, holderz, holderzoffset;
		Bool_t ishpaddle;
		ComputePaddleTransformation(paddle, origin, ishpaddle);
		origin=TVector3(origin.X(),origin.Y(),origin.Z());
#ifdef DRAWSUPERVERBOSE
		cout<<"paddle "<<paddle<<" at ("<<origin.X()<<", "<<origin.Y()<<", "<<origin.Z()<<") was ";
		(paddleishit) ? cout<<" struck at "<<thetime : cout<<" not hit"; cout<<endl;
#endif
		
		// establish if this hits a new half of the horizontal or vertical panel.
		if(paddleishit&&origin.Y()>0){tophit=true;}
		else if(paddleishit&&origin.Y()<0){bottomhit=true;}
		if(paddleishit&&origin.X()>0){rightsidehit=true;}
		else if(paddleishit&&origin.X()<0){leftsidehit=true;}
#ifdef DRAWSUPERVERBOSE
		cout<<"updated stats are: "<<"tophit: "<<tophit<<", bottomhit: "<<bottomhit
			<<", rightsidehit: "<<rightsidehit<<", leftsidehit: "<<leftsidehit<<endl;
#endif
		
		// convert from simulation 'cm' units to canvas units: 0-1, 0-1 horizontal and vertical
		// flip x sign, a "top view" with z from Left->Right has 'RHS' for x>0, but that corresponds
		// to the LOWER (negative) half of the canvas -> i.e. x>0 should maps to canvas_y<0. 
		holderx = 0.5+((-origin.X()/(MRDSpecs::maxwidth*1.2)));
		Double_t anoffset=0;
		if((ishpaddle&&(origin.X()>0))||((!ishpaddle)&&(origin.Y()<0))){anoffset=(MRDSpecs::scintfullzlen+MRDSpecs::scintalugap)*5;}
		// in order to view the two sets in the same graph, shift one half on the canvas
		holdery = (origin.Y()/(MRDSpecs::maxheight*1.2))+0.5;
		holderz = ((origin.Z()+anoffset)/(MRDSpecs::mrdZlen*1.2))+0.5;
//		cout<<"origin is: ("<<origin.X()<<","<<origin.Y()<<","<<origin.Z()<<")"<<endl;
//		cout<<"scaled origin is: ("<<holderx<<","<<holdery<<","<<holderz<<")"<<endl;
//		cout<<"ishpaddle="<<ishpaddle<<endl;
//		if(paddle==0){
//			cout<<"drawing box with height "<<(MRDSpecs::scintfullxlen/(MRDSpecs::maxwidth*1.2))
//				<<" and width "<<(scintboxwidth/(MRDSpecs::mrdZlen*1.2))<<endl;
//		}
		
		/*
		// make a note of the first and last paddles hit in each view. These are used to do a very simple
		// track fit, just by drawing a straight line between them. --------Currently not used---------
		if((!ishpaddle)&&(!firsthpaddledone)&&paddleishit){	//TODO: why !ishpaddle instead of ishpaddle??
			firsthpaddledone=true;
			xupcorner1=std::pair<double,double>(holderx-(sideboxheight/2.),holderz);
			xdowncorner1=std::pair<double,double>(holderx+(sideboxheight/2.),holderz+(scintboxwidth/MRDSpecs::mrdZlen));
			//cout<<"first hit h paddle at ("<<holderx<<","<<holdery<<","<<holderz<<")"<<endl;
		}
		if((!ishpaddle)&&paddleishit){
			xdowncorner2=std::pair<double,double>(holderx-(sideboxheight/2.),holderz);
			xupcorner2=std::pair<double,double>(holderx+(sideboxheight/2.),holderz+(scintboxwidth/MRDSpecs::mrdZlen));
			//cout<<"last h paddle at ("<<holderx<<","<<holdery<<","<<holderz<<")"<<endl;
		}
		if(ishpaddle&&(!firstvpaddledone)&&paddleishit){
			firstvpaddledone=true;
			yupcorner1=std::pair<double,double>(holdery-(topboxheight/2.),holderz);
			ydowncorner1=std::pair<double,double>(holdery+(topboxheight/2.),holderz+(scintboxwidth/MRDSpecs::mrdZlen));
			//cout<<"first hit v paddle at ("<<holderx<<","<<holdery<<","<<holderz<<")"<<endl;
		}
		if(ishpaddle&&paddleishit){
			ydowncorner2=std::pair<double,double>(holdery-(topboxheight/2.),holderz);
			yupcorner2=std::pair<double,double>(holdery+(topboxheight/2.),holderz+(scintboxwidth/MRDSpecs::mrdZlen));
			//cout<<"last h paddle at ("<<holderx<<","<<holdery<<","<<holderz<<")"<<endl;
		}
		// -----------------------------
		*/
		
		// Draw the box
		TBox *thepaddle;
		if(paddlepointers.at(paddle)!=0){
#ifdef DRAWSUPERVERBOSE
			cout<<"getting box pointer for paddle "<<paddle<<endl;
#endif
			thepaddle = paddlepointers.at(paddle);
		} else {
#ifdef DRAWSUPERVERBOSE
			cout<<"making box pointer for paddle "<<paddle<<endl;
#endif
			if(!ishpaddle){
				//TBox (Double_t leftx, Double_t bottomy, Double_t rightx, Double_t topy)
				thepaddle = new TBox(holderz,holderx-(sideboxheight/2.),holderz+(scintboxwidth/MRDSpecs::mrdZlen),holderx+(sideboxheight/2.));
			} else {
				thepaddle = new TBox(holderz,holdery-(topboxheight/2.),holderz+(scintboxwidth/MRDSpecs::mrdZlen),holdery+(topboxheight/2.));
			}
			paddlepointers.at(paddle) = thepaddle;
		}
		thepaddle->SetFillStyle(1001);	// solid fill
		if(paddleishit){ thepaddle->SetFillColor(paddlecolour); } else { thepaddle->SetFillColor(11); } 
		// if the paddle was hit light it up - color 5 (kYellow) otherwise leave it dark - colour 11 (grey)
		thepaddle->SetLineColor(1);	//black
		thepaddle->SetLineWidth(2);
		(ishpaddle) ? imgcanvas->cd(1) : imgcanvas->cd(2);
		thepaddle->Draw();
		
		// check if we're changing panel
		if( paddle!=0 && (std::count(MRDSpecs::layeroffsets.begin(), MRDSpecs::layeroffsets.end(), (paddle+1))!=0) ){
#ifdef DRAWSUPERVERBOSE
			cout<<"drawing paddles in alternate view"<<endl;
#endif
			Double_t holderx1, holderx2, holdery1, holdery2;
			Double_t otherviewoffset;
			// this is the last paddle of this layer - draw the layer in the opposite view, combining all paddles
			if(!ishpaddle){
#ifdef DRAWSUPERVERBOSE
				cout<<"drawing paddles in side view"<<endl;
#endif
				otherviewoffset=-0.014;
				imgcanvas->cd(1);
				holderx1=((MRDSpecs::scintalugap*10)/(MRDSpecs::maxwidth*1.2));
				holderx2=(MRDSpecs::scinthfullylen/(MRDSpecs::maxwidth*1.2)); // TODO: wrong length, should be MRDSpecs::scintvfullylen
				// some paddles are shown in both views, so we have more TBoxes than real PMTs
				Int_t overflowindex = MRDSpecs::nummrdpmts + (2*mrdcluster::paddle_layers.at(paddle));	// RH paddle
				if(paddlepointers.at(overflowindex)==0){
					thepaddle = new TBox(holderz+otherviewoffset,0.5+holderx1,holderz+(scintboxwidth/MRDSpecs::mrdZlen)+otherviewoffset, 0.5+holderx1+holderx2);
					paddlepointers.at(overflowindex) = thepaddle;
				} else {
					thepaddle = paddlepointers.at(overflowindex);
				}
				thepaddle->SetFillStyle(1001);	// solid fill
				if(tophit){ thepaddle->SetFillColor(paddlecolour); } else { thepaddle->SetFillColor(11); }
#ifdef DRAWSUPERVERBOSE
				(tophit) ? cout<<"top paddle is hit"<<endl : cout<<"top paddle not hit"<<endl;
#endif
				// kYellow if hit or grey otherwise
				thepaddle->SetLineColor(1);	//black
				thepaddle->SetLineWidth(2);
				thepaddle->Draw();
				overflowindex++;	// left side paddle
				if(paddlepointers.at(overflowindex)==0){
					thepaddle = new TBox(holderz+otherviewoffset,0.5-holderx1,holderz+(scintboxwidth/MRDSpecs::mrdZlen)+otherviewoffset, 0.5-(holderx1+holderx2));
					paddlepointers.at(overflowindex) = thepaddle;
				} else {
					thepaddle = paddlepointers.at(overflowindex);
				}
				thepaddle->SetFillStyle(1001);	// solid fill
				if(bottomhit){ thepaddle->SetFillColor(paddlecolour); } else { thepaddle->SetFillColor(11); }
#ifdef DRAWSUPERVERBOSE
				(bottomhit) ? cout<<"bottom paddle is hit"<<endl : cout<<"bottom paddle not hit"<<endl;
#endif
				// kYellow if hit or grey otherwise
				thepaddle->SetLineColor(1);	//black
				thepaddle->SetLineWidth(2);
				thepaddle->Draw();
				
			} else {
#ifdef DRAWSUPERVERBOSE
				cout<<"drawing paddles in top view"<<endl;
#endif
				otherviewoffset=-0.01;
				imgcanvas->cd(2);
				holdery1=((MRDSpecs::scintalugap*10)/(MRDSpecs::maxheight*1.2));
				holdery2=(MRDSpecs::scintvfullylen/(MRDSpecs::maxheight*1.2));
				Int_t overflowindex = MRDSpecs::nummrdpmts + (2*mrdcluster::paddle_layers.at(paddle));	// top paddle
				if(paddlepointers.at(overflowindex)==0){
					thepaddle = new TBox(holderz+otherviewoffset,holdery1+0.5,holderz+(scintboxwidth/MRDSpecs::mrdZlen)+otherviewoffset,0.5+holdery1+holdery2);
					paddlepointers.at(overflowindex) = thepaddle;
				} else {
					thepaddle = paddlepointers.at(overflowindex);
				}
				thepaddle->SetFillStyle(1001);	// solid fill
				if(leftsidehit){ thepaddle->SetFillColor(paddlecolour); }
				else { thepaddle->SetFillColor(11); }
#ifdef DRAWSUPERVERBOSE
				(leftsidehit) ? cout<<"left paddle is hit"<<endl : cout<<"left paddle not hit"<<endl;
#endif
				// kYellow if hit or grey otherwise
				thepaddle->SetLineColor(1);	//black
				thepaddle->SetLineWidth(2);
				thepaddle->Draw();
				overflowindex++;	// bottom paddle
				if(paddlepointers.at(overflowindex)==0){
					thepaddle = new TBox(holderz+otherviewoffset,0.5-holdery1,holderz+(scintboxwidth/MRDSpecs::mrdZlen)+otherviewoffset,0.5-(holdery1+holdery2));
					paddlepointers.at(overflowindex) = thepaddle;
				} else {
					thepaddle = paddlepointers.at(overflowindex);
				}
				thepaddle->SetFillStyle(1001);	// solid fill
				if(rightsidehit){ thepaddle->SetFillColor(paddlecolour); }
				else { thepaddle->SetFillColor(11); }
#ifdef DRAWSUPERVERBOSE
				(rightsidehit) ? cout<<"right paddle is hit"<<endl : cout<<"right paddle not hit"<<endl;
#endif
				// kYellow if hit or grey otherwise
				thepaddle->SetLineColor(1);	//black
				thepaddle->SetLineWidth(2);
				thepaddle->Draw();
			}
			
			// reset for the next layer
			tophit=false; bottomhit=false; leftsidehit=false; rightsidehit=false;
#ifdef DRAWVERBOSE
			cout<<"layer finished, moving to next one"<<endl;
#endif
		}
	}
	//gPad->WaitPrimitive();	// wait for user to click
	// imgcanvas->SaveAs(TString::Format("mrdtracks_%d.png",event_id)); // invoke after adding tracks
//	pmts_hit = holdme;
//	digi_ts = holdmetoo;
}

// code taken from MRDDetectorConstruction.cc, WCSimDetectorConstruction::ComputePaddleTransformation
// with comments stripped, G4int->Int_t, G4bool->Bool_t, G4ThreeVector->TVector3, remove G4RotationMatrix
//============================= 
void cMRDSubEvent::ComputePaddleTransformation (const Int_t copyNo, TVector3 &origin, Bool_t &ishpaddle) {
	Double_t Xposition=0, Yposition=0, Zposition=0;
	Int_t panelpairnum = floor(copyNo/(MRDSpecs::numpaddlesperpanelv+MRDSpecs::numpaddlesperpanelh));
	Int_t panelnumrem = copyNo - panelpairnum*(MRDSpecs::numpaddlesperpanelv+MRDSpecs::numpaddlesperpanelh);
	Int_t panelnum;
	Int_t paddlenum;
	ishpaddle=false;
	if(panelnumrem>(MRDSpecs::numpaddlesperpanelh-1)){
		panelnum = (panelpairnum*2) +1;
		paddlenum = panelnumrem-MRDSpecs::numpaddlesperpanelh;
	} else {
		panelnum = (panelpairnum*2);
		ishpaddle = true;
		paddlenum = panelnumrem;
	}
	Int_t pairnum = floor(paddlenum/2);
	Zposition = panelnum*(MRDSpecs::steelfullzlen + MRDSpecs::alufullzlen + MRDSpecs::scintfullzlen + MRDSpecs::layergap);
	Zposition = Zposition + MRDSpecs::steelfullzlen + MRDSpecs::steelscintgap;
	Zposition = Zposition + (MRDSpecs::scintfullzlen/2);
	Zposition = Zposition + MRDSpecs::MRDPMTRadius - (MRDSpecs::mrdZlen/2);
	
	if (ishpaddle){
		if (paddlenum%2==0){
			Xposition=((MRDSpecs::scinthfullylen+MRDSpecs::scintbordergap)/2);
		} else {
			Xposition=-((MRDSpecs::scinthfullylen+MRDSpecs::scintbordergap)/2);
		}
		Yposition = pairnum*(MRDSpecs::scintfullxlen+MRDSpecs::scintbordergap);
		Yposition = Yposition - 0.5*(((MRDSpecs::scintfullxlen+MRDSpecs::scintbordergap)/2)*MRDSpecs::numpaddlesperpanelh)+(MRDSpecs::scintfullxlen/2);
	} else {
		if (paddlenum%2==0){
			Yposition=((MRDSpecs::scintvfullylen+MRDSpecs::scintbordergap)/2); 
		} else {
			Yposition=-((MRDSpecs::scintvfullylen+MRDSpecs::scintbordergap)/2);
		}
		Xposition = pairnum*(MRDSpecs::scintfullxlen+MRDSpecs::scintbordergap);
		Xposition = Xposition - 0.5*(((MRDSpecs::scintfullxlen+MRDSpecs::scintbordergap)/2)*MRDSpecs::numpaddlesperpanelv)+(MRDSpecs::scintfullxlen/2); 
	}
	TVector3 theposition(Xposition,Yposition,Zposition);
	origin=theposition;
}
