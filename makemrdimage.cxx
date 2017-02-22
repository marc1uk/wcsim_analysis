/* vim:set noexpandtab tabstop=4 wrap */
#include "TRandom.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TLine.h"

#ifndef DRAWVERBOSE
//#define DRAWVERBOSE 1
#endif

void cMRDTrack::DrawMrdCanvases(){
// TODO: need to handle cleanup
#ifdef DRAWVERBOSE
	cout<<"making canvas plots"<<endl;
#endif
	Double_t scintboxwidth=1;
	Double_t topboxheight = (scintfullxlen/(maxwidth*1.2));
	Double_t sideboxheight = (scintfullxlen/(maxwidth*1.3));
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
		imgcanvas = new TCanvas("imgcanvas","MRD Digit Visualiser",canvw,canvh); 
		imgcanvas->Divide(2,1);
		imgcanvas->cd(1);
		titleleft = new TText(.32,.9,TString::Format("Top View, Event %d",event_id));
		titleleft->Draw();
		imgcanvas->cd(2);
		titleright = new TText(.32,.9,TString::Format("Side View, Event %d",event_id));
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
//		imgcanvas->cd(1);
//		gPad->Clear();
//		titleleft->Draw();
//		imgcanvas->cd(2);
//		gPad->Clear();
//		titleright->Draw();
		titleleft->SetTitle(TString::Format("Top View, Event %d",event_id));
		titleleft->Draw();
		titleright->SetTitle(TString::Format("Side View, Event %d",event_id));
		titleright->Draw();
		char titlebuf[50];
		snprintf(titlebuf,50,"%d",(int)(maxtime*2));
		texttop->SetTitle(titlebuf);
		texttop->Draw();
		snprintf(titlebuf,50,"%d",(int)mintime);
		textbottom->SetTitle(titlebuf);
		textbottom->Draw();
	}

	//cout<<"adding scints"<<endl;
	Bool_t firsthpaddledone=false, firstvpaddledone=false;
	//std::pair<double, double> xupcorner1, xupcorner2, xdowncorner1, xdowncorner2, yupcorner1, yupcorner2, ydowncorner1, ydowncorner2;	// part of cMRDTrack class now
	Bool_t rightsidehit=false, leftsidehit=false, tophit=false, bottomhit=false;
	
	// loop over all paddles and build a map of which were hit
	for(int paddle=0; paddle<((numpaddlesperpanelv*numvpanels)+(numpaddlesperpanelh*numhpanels)); paddle++){
	
		// check if the paddle was hit by finding it's PMT ID in the struck PMTs
		Bool_t paddleishit;
		if(std::find(pmts_hit.begin(), pmts_hit.end(), paddle)!=pmts_hit.end()){
			paddleishit=true;
			// grab the time the digit on this PMT. (what to do with more than one?)
			std::vector<Int_t>::iterator theit = std::find(pmts_hit.begin(), pmts_hit.end(), paddle);
			Int_t theindex = std::distance(pmts_hit.begin(),theit);
			Double_t thetime = digi_ts.at(theindex);
			Double_t relatime = (thetime-mintime)/(maxtime-mintime);
			Int_t colorindex = TMath::Floor((aspectrumv.size()-1)*(relatime/2));
			paddlecolour = aspectrumv.at(colorindex);
			// kRed is an EColor, kRed+2 is an Int_t, representing generically a TColor, TBox requires a Color_t!
		} else {
			paddleishit=false;
		}
		
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
			xupcorner1=std::pair<double,double>(holderx-(sideboxheight/2.),holderz);
			xdowncorner1=std::pair<double,double>(holderx+(sideboxheight/2.),holderz+(scintboxwidth/mrdZlen));
			//cout<<"first hit h paddle at ("<<holderx<<","<<holdery<<","<<holderz<<")"<<endl;
		}
		if((!ishpaddle)&&paddleishit){
			xdowncorner2=std::pair<double,double>(holderx-(sideboxheight/2.),holderz);
			xupcorner2=std::pair<double,double>(holderx+(sideboxheight/2.),holderz+(scintboxwidth/mrdZlen));
			//cout<<"last h paddle at ("<<holderx<<","<<holdery<<","<<holderz<<")"<<endl;
		}
		if(ishpaddle&&(!firstvpaddledone)&&paddleishit){
			firstvpaddledone=true;
			yupcorner1=std::pair<double,double>(holdery-(topboxheight/2.),holderz);
			ydowncorner1=std::pair<double,double>(holdery+(topboxheight/2.),holderz+(scintboxwidth/mrdZlen));
			//cout<<"first hit v paddle at ("<<holderx<<","<<holdery<<","<<holderz<<")"<<endl;
		}
		if(ishpaddle&&paddleishit){
			ydowncorner2=std::pair<double,double>(holdery-(topboxheight/2.),holderz);
			yupcorner2=std::pair<double,double>(holdery+(topboxheight/2.),holderz+(scintboxwidth/mrdZlen));
			//cout<<"last h paddle at ("<<holderx<<","<<holdery<<","<<holderz<<")"<<endl;
		}
		
		// Draw the box
		TBox *thepaddle;
		if(!ishpaddle){
			imgcanvas->cd(2);
			//TBox (Double_t leftx, Double_t bottomy, Double_t rightx, Double_t topy)
			if(paddlepointers.at(paddle)==0){
				thepaddle = new TBox(holderz,holderx-(sideboxheight/2.),holderz+(scintboxwidth/mrdZlen),holderx+(sideboxheight/2.));
			} else {
				thepaddle = paddlepointers.at(paddle);
			}
		} else {
			imgcanvas->cd(1);
			if(paddlepointers.at(paddle)==0){
				thepaddle = new TBox(holderz,holdery-(topboxheight/2.),holderz+(scintboxwidth/mrdZlen),holdery+(topboxheight/2.));
			} else {
				thepaddle = paddlepointers.at(paddle);
			}
		}
		thepaddle->SetFillStyle(1001);	// solid fill
		if(paddleishit){ thepaddle->SetFillColor(paddlecolour); } else { thepaddle->SetFillColor(11); } 
		// if the paddle was hit light it up - color 5 (kYellow) otherwise leave it dark - colour 11 (grey)
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
				imgcanvas->cd(1);
				holderx1=((scintalugap*10)/(maxwidth*1.2));
				holderx2=(scinthfullylen/(maxwidth*1.2));
				if(paddlepointers.at(paddle)==0){
					thepaddle = new TBox(holderz+otherviewoffset,0.5+holderx1,holderz+(scintboxwidth/mrdZlen)+otherviewoffset, 0.5+holderx1+holderx2);
				} else {
					thepaddle = paddlepointers.at(paddle);
				}
				thepaddle->SetFillStyle(1001);	// solid fill
				if(rightsidehit){ thepaddle->SetFillColor(paddlecolour); } else { thepaddle->SetFillColor(11); }
				// kYellow if hit or grey otherwise
				thepaddle->SetLineColor(1);	//black
				thepaddle->SetLineWidth(2);
				thepaddle->Draw();
				if(paddlepointers.at(paddle)==0){
					thepaddle = new TBox(holderz+otherviewoffset,0.5-holderx1,holderz+(scintboxwidth/mrdZlen)+otherviewoffset, 0.5-(holderx1+holderx2));
				} else {
					thepaddle = paddlepointers.at(paddle);
				}
				thepaddle->SetFillStyle(1001);	// solid fill
				if(leftsidehit){ thepaddle->SetFillColor(paddlecolour); } else { thepaddle->SetFillColor(11); }
				// kYellow if hit or grey otherwise
				thepaddle->SetLineColor(1);	//black
				thepaddle->SetLineWidth(2);
				thepaddle->Draw();
				
			} else {
				otherviewoffset=-0.01;
				imgcanvas->cd(2);
				holdery1=((scintalugap*10)/(maxheight*1.2));
				holdery2=(scintvfullylen/(maxheight*1.2));
				if(paddlepointers.at(paddle)==0){
					thepaddle = new TBox(holderz+otherviewoffset,holdery1+0.5,holderz+(scintboxwidth/mrdZlen)+otherviewoffset,0.5+holdery1+holdery2);
				} else {
					thepaddle = paddlepointers.at(paddle);
				}
				thepaddle->SetFillStyle(1001);	// solid fill
				if(tophit){ thepaddle->SetFillColor(paddlecolour); } else { thepaddle->SetFillColor(11); }
				// kYellow if hit or grey otherwise
				thepaddle->SetLineColor(1);	//black
				thepaddle->SetLineWidth(2);
				thepaddle->Draw();
				if(paddlepointers.at(paddle)==0){
					thepaddle = new TBox(holderz+otherviewoffset,0.5+-holdery1,holderz+(scintboxwidth/mrdZlen)+otherviewoffset,0.5-(holdery1+holdery2));
				} else {
					thepaddle = paddlepointers.at(paddle);
				}
				thepaddle->SetFillStyle(1001);	// solid fill
				if(bottomhit){ thepaddle->SetFillColor(paddlecolour); } else { thepaddle->SetFillColor(11); }
				// kYellow if hit or grey otherwise
				thepaddle->SetLineColor(1);	//black
				thepaddle->SetLineWidth(2);
				thepaddle->Draw();
			}
			
			// reset for the next layer
			tophit=false; bottomhit=false; leftsidehit=false; rightsidehit=false;
		}
		
	}
	//gPad->WaitPrimitive();	// wait for user to click
	imgcanvas->SaveAs(TString::Format("mrdtracks_%d.png",event_id));
}

void cMRDTrack::AddTrackLines(){
//TODO: rather than projecting to 0.8, project to the appropriate further depth of penetration.
	// add trajectory lines 
	//=====================
	// TLine::TLine(Double_t x1, Double_t y1, Double_t x2, Double_t y2))
	std::map<const char*, double> xlines = this->GetTrackXStats();
	std::map<const char*, double> ylines = this->GetTrackYStats();
	// has keys "angminxval", "angminzval", "angmin", same for angmax
	Double_t holderx1, holdery1, holderz1, holderx2, holdery2, holderz2, projectedleft, projectedright;
	TLine *xupgoing, *xdowngoing, *yupgoing, *ydowngoing;
	Bool_t usebasic=false, useclass=true;
	
	// first upgoing x line
#ifdef DRAWVERBOSE
	cout<<"making basic upgoing x line"<<endl;
#endif
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
		imgcanvas->cd(2);
		xupgoing->Draw();
	} 
#ifdef DRAWVERBOSE
	cout<<"making class upgoing x line"<<endl;
#endif
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
		imgcanvas->cd(2);
		xupgoing->Draw();
	}
	
#ifdef DRAWVERBOSE
	cout<<"making basic downgoing x line"<<endl;
#endif
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
		imgcanvas->cd(2);
		xdowngoing->Draw();
	}
#ifdef DRAWVERBOSE
	cout<<"making class downgoing x line"<<endl;
#endif
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
		imgcanvas->cd(2);
		xdowngoing->Draw();
	}
#ifdef DRAWVERBOSE
	cout<<"making basic upgoing y line"<<endl;
#endif
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
		imgcanvas->cd(1);
		yupgoing->Draw();
	}
#ifdef DRAWVERBOSE
	cout<<"making class upgoing y line"<<endl;
#endif
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
		imgcanvas->cd(1);
		yupgoing->Draw();
	}
	
	// finally downgoing y line
#ifdef DRAWVERBOSE
	cout<<"making basic downgoing y line"<<endl;
#endif
	if(usebasic){
		holdery1=ydowncorner1.first; holderz1=ydowncorner1.second;
		holdery2=ydowncorner2.first; holderz2=ydowncorner2.second;
		projectedleft = holdery1+holderz1*((holdery1-holdery2)/(holderz2-holderz1));
		ydowngoing = new TLine(0.,projectedleft,holderz2, holdery2);
		//--------------------------
		ydowngoing->SetLineWidth(4);
		ydowngoing->SetLineColor(kRed);
		ydowngoing->SetLineStyle(3);
		imgcanvas->cd(1);
		ydowngoing->Draw();
	}
#ifdef DRAWVERBOSE
	cout<<"making class downgoing y line"<<endl;
#endif
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
		imgcanvas->cd(1);
		ydowngoing->Draw();
	}
}

// code taken from MRDDetectorConstruction.cc, WCSimDetectorConstruction::ComputePaddleTransformation
// with comments stripped, G4int->Int_t, G4bool->Bool_t, G4ThreeVector->TVector3, remove G4RotationMatrix
//============================= 
void cMRDTrack::ComputePaddleTransformation (const Int_t copyNo, TVector3 &origin, Bool_t &ishpaddle, Bool_t paddleishit) {
	Double_t Xposition=0, Yposition=0, Zposition=0;
	Int_t panelpairnum = floor(copyNo/(numpaddlesperpanelv+numpaddlesperpanelh));
	Int_t panelnumrem = copyNo - panelpairnum*(numpaddlesperpanelv+numpaddlesperpanelh);
	Int_t panelnum;
	Int_t paddlenum;
	ishpaddle=false;
	if(panelnumrem>(numpaddlesperpanelh-1)){
		panelnum = (panelpairnum*2) +1;
		paddlenum = panelnumrem-numpaddlesperpanelh;
	} else {
		panelnum = (panelpairnum*2);
		ishpaddle = true;
		paddlenum = panelnumrem;
	}
	Int_t pairnum = floor(paddlenum/2);
	Zposition = panelnum*(steelfullzlen + alufullzlen + scintfullzlen + layergap);
	Zposition = Zposition + steelfullzlen + steelscintgap;
	Zposition = Zposition + (scintfullzlen/2);
	Zposition = Zposition + MRDPMTRadius - (mrdZlen/2);
	
	if (ishpaddle){
		if (paddlenum%2==0){
			Xposition=((scinthfullylen+scintbordergap)/2);
		} else {
			Xposition=-((scinthfullylen+scintbordergap)/2);
		}
		Yposition = pairnum*(scintfullxlen+scintbordergap);
		Yposition = Yposition - 0.5*(((scintfullxlen+scintbordergap)/2)*numpaddlesperpanelh)+(scintfullxlen/2);
	} else {
		if (paddlenum%2==0){
			Yposition=((scintvfullylen+scintbordergap)/2); 
		} else {
			Yposition=-((scintvfullylen+scintbordergap)/2);
		}
		Xposition = pairnum*(scintfullxlen+scintbordergap);
		Xposition = Xposition - 0.5*(((scintfullxlen+scintbordergap)/2)*numpaddlesperpanelv)+(scintfullxlen/2); 
	}
	TVector3 theposition(Xposition,Yposition,Zposition);
	origin=theposition;
}
