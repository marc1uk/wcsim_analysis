/* vim:set noexpandtab tabstop=4 wrap */
#ifndef MRDSUBEVENTDRAWTRUETRACKS
#define MRDSUBEVENTDRAWTRUETRACKS

#ifndef DRAWTRUEVERBOSE
//#define DRAWTRUEVERBOSE
#endif

bool CheckLineBox( TVector3 L1, TVector3 L2, TVector3 B1, TVector3 B2, TVector3 &Hit, TVector3 &Hit2, bool &error);
int inline InBox( TVector3 Hit, TVector3 B1, TVector3 B2, const int Axis);
int inline GetIntersection( float fDst1, float fDst2, TVector3 P1, TVector3 P2, TVector3 &Hit);

void cMRDSubEvent::DrawTrueTracks(){
	// to visualize the accuracy and efficiency of track finding, as well as diagnose problems
	// it is very helpful to draw the truth tracks as arrows on the reconstruction plot
	// truth tracks are assumed to be straight lines from their start to finish (only info we have)
	// we have std::vector<WCSimRootTrack> truetracks; to pull the track info from
	// we have std::vector<TArrow*> truetrackarrows; to add the arrows to
		
#ifdef DRAWTRUEVERBOSE
	cout<<"MRD width is "<<MRDSpecs::MRD_width<<", MRDSpecs::maxwidth is "<<(MRDSpecs::maxwidth/2.)<<", MRD height is "<<MRDSpecs::MRD_height<<", MRDSpecs::maxheight is "<<(MRDSpecs::maxheight/2.)<<", MRD start is "<<MRDSpecs::MRD_start<<", MRD end is "<<(MRDSpecs::MRD_start+MRDSpecs::MRD_depth)<<", MRD depth is "<<(MRDSpecs::MRD_depth/2.)<<", MRDSpecs::mrdZlen is "<<(MRDSpecs::mrdZlen/2.)<<endl;
#endif
	
#ifdef DRAWTRUEVERBOSE
	cout<<"Adding lines for "<<truetracks.size()<<" true tracks"<<endl;
#endif

/////////////////////////////////////////////////////////////////
//	std::vector<double> placeholderx;   // start points in the front H/V layer of the MRD
//	std::vector<double> placeholdery;
//	std::vector<double> placeholderz;
//	std::vector<bool> placeholderorient;
//	std::vector<double> placeholderx2;  // end points in the last H/V layer of the MRD
//	std::vector<double> placeholdery2;
//	std::vector<double> placeholderz2;
//	for(int paddle=numpaddlesperpanelh-1; paddle>0; paddle-=2){
//		TVector3 origin;
//		Bool_t ishpaddle;
//		ComputePaddleTransformation(/*paddle*/numpaddlesperpanelh-1, origin, ishpaddle);
//		cout<<"origin "<<paddle<<" is ("<<origin.X()<<", "<<origin.Y()<<", "<<origin.Z()<<")"<<endl;
//		placeholderx.push_back(origin.X());
//		placeholdery.push_back(origin.Y());
//		placeholderz.push_back(origin.Z());
//		placeholderorient.push_back(ishpaddle);
//		
//		int paddle2 = (ishpaddle) ? paddle+(5*numpaddlesperpanelv)+(5*numpaddlesperpanelh) +3
//		                     : paddle+(4*numpaddlesperpanelv)+(4*numpaddlesperpanelh) +3;
//		paddle2=numpaddlesperpanelh+numpaddlesperpanelv+paddle;
//		ComputePaddleTransformation(paddle2, origin, ishpaddle);
//		cout<<"origin2 "<<paddle2<<" is ("<<origin.X()<<", "<<origin.Y()<<", "<<origin.Z()<<")"<<endl;
//		placeholderx2.push_back(origin.X());
//		placeholdery2.push_back(origin.Y());
//		placeholderz2.push_back(origin.Z());
//		//break;
//	}
/////////////////////////////////////////////////////////////////
	//gStyle->SetOptStat(1);    // show stats / statistics box
	//gStyle->SetOptStat(111111); //show overflow & underflow bins
//	if(firsttime==true){
//		cout<<"creating histograms:"<<endl;
//		cout<<"startxhist="<<startxhist<<endl
//			<<"startyhist="<<startyhist<<endl
//			<<"startzhist="<<startzhist<<endl;
//		firsttime=false;
//	} else {
//		cout<<"retrieving histograms:"<<endl;
//		cout<<"startxhist="<<startxhist<<endl
//			<<"startyhist="<<startyhist<<endl
//			<<"startzhist="<<startzhist<<endl;
//	}
	
	int trackit=0;
	std::vector<int> placeholder(2,0);
	for(auto nextrack : truetracks){ //truetracks, placeholderx
		
//		// TODO maybe we should colour the track by time?
//		Double_t thetime = digi_ts.at(theindex);
//		Double_t relatime = (thetime-mintime)/(maxtime-mintime);
//		Int_t colorindex = TMath::Floor((aspectrumv.size()-1)*(relatime/2));
//		paddlecolour = aspectrumv.at(colorindex);
		
////////////////////////////////////////////////////////////////////////////////////
//// REAPLCEMENT STARTS HERE
		TLorentzVector primarystartvertex(  nextrack.GetStart(0),
											nextrack.GetStart(1),
											nextrack.GetStart(2),
											nextrack.GetTime());
#if FILE_VERSION>2
		TLorentzVector primarystopvertex(   nextrack.GetStop(0),
											nextrack.GetStop(1),
											nextrack.GetStop(2),
											nextrack.GetStopTime());
#else
		TLorentzVector primarystopvertex(   nextrack.GetStop(0),
											nextrack.GetStop(1),
											nextrack.GetStop(2),
											-1); // not stored prior to this
#endif
		// Draw vector indicating the true track's intial direction, for better idea of trajectory
		double dirx = nextrack.GetDir(0);
		double diry = nextrack.GetDir(1);
		double dirz = nextrack.GetDir(2);
		
		// skip the easy ones
		if(primarystopvertex.Z()<MRDSpecs::MRD_start){
#ifdef DRAWTRUEVERBOSE
			cout<<"skipping track "<<trackit<<" as it doesn't make it to mrd: trackstopz="
				<<primarystopvertex.Z()<<", MRDSpecs::MRD_start="<<MRDSpecs::MRD_start<<endl;
#endif
			continue;
		}  // the track never makes it to the MRD
		if(primarystopvertex.Z()>(MRDSpecs::MRD_start+MRDSpecs::MRD_depth)){
#ifdef DRAWTRUEVERBOSE
			cout<<"skipping track "<<trackit<<" as it starts after mrd: trackstartz="
				<<primarystopvertex.Z()<<", MRDSpecs::MRD_start+MRDSpecs::MRD_depth="<<MRDSpecs::MRD_start+MRDSpecs::MRD_depth<<endl;
#endif
			continue;
		} // the track is somehow created after the MRD

#ifdef DRAWTRUEVERBOSE
		cout<<"Track start: ("<<primarystartvertex.X()<<", "<<primarystartvertex.Y()
			<<", "<<primarystartvertex.Z()<<"), track end: ("<<primarystopvertex.X()
			<<", "<<primarystopvertex.Y()<<", "<<primarystopvertex.Z()<<")"<<endl;
#endif
		//new MRD entry/exit point calculation based on external function calls
		///////////////////////////////////////////////////////////////////////
		// bool CheckLineBox( TVector3 L1, TVector3 L2, TVector3 B1, TVector3 B2, 
		//					  TVector3 &Hit, TVector3 &Hit2, bool &error)
		// returns true if line (L1, L2) intersects with the box (B1, B2)
		// returns intersection with smaller Z in Hit
		// if 2 interceptions are found, returns interception with larger Z,
		// if 1 interception is found, returns L2 (stopping point).
		// error returns true if >2 intercepts are found, or other error.
		TVector3 MRDentrypoint(0,0,0), MRDexitpoint(0,0,0), MuTrackInMRD(0,0,0);
		double mutracklengthinMRD, mrdpenetrationcm;
		int mrdpenetrationlayers;
		// check for intercept and record entry point
		bool checkboxlinerror=false, muonentersMRD=false, muonstopsinMRD, muonrangesoutMRD;
		muonentersMRD  =  CheckLineBox( primarystartvertex.Vect(), primarystopvertex.Vect(), 
										TVector3(-MRDSpecs::MRD_width,-MRDSpecs::MRD_height,MRDSpecs::MRD_start), 
										TVector3(MRDSpecs::MRD_width,MRDSpecs::MRD_height,MRDSpecs::MRD_end),
										MRDentrypoint, MRDexitpoint, checkboxlinerror );
		// sanity check: XXX DISABLE TO ALLOW TRACKS STARTING IN THE MRD
		//assert(MRDentrypoint!=primarystartvertex.Vect()&&"track starts in MRD!?");
		// check if MRD stops in the MRD
		muonstopsinMRD = ( abs(primarystopvertex.X())<MRDSpecs::MRD_width&&
						   abs(primarystopvertex.Y())<MRDSpecs::MRD_height&&
						   primarystopvertex.Z()>MRDSpecs::MRD_start&&
						   primarystopvertex.Z()<(MRDSpecs::MRD_start+MRDSpecs::MRD_depth) );
		if(muonentersMRD){
			muonrangesoutMRD = ((MRDentrypoint.Z()==MRDSpecs::MRD_start)&&(MRDexitpoint.Z()==MRDSpecs::MRD_end));
			MuTrackInMRD = (MRDexitpoint-MRDentrypoint);
			mutracklengthinMRD = MuTrackInMRD.Mag();
			mrdpenetrationcm = MuTrackInMRD.Z();
			mrdpenetrationlayers=0;
			for(auto layerzval : MRDSpecs::mrdscintlayers){
				if(MRDexitpoint.Z()<layerzval) break;
				mrdpenetrationlayers++;
			}
		} else {
			muonrangesoutMRD=false;
			mrdpenetrationcm=0.;
			mrdpenetrationlayers=0;
			mutracklengthinMRD=0.;
		}
		
		// old stuff needed for backward compatibility with remainder of code
		bool trackisbackwardgoing=false;
		double frontx,fronty,backz,backx,backy;
		double frontz=min(primarystartvertex.Z(),primarystopvertex.Z());
		if(frontz==primarystartvertex.Z()){
			frontx=primarystartvertex.X();
			fronty=primarystartvertex.Y();
			backx=primarystopvertex.X();
			backy=primarystopvertex.Y();
			backz=primarystopvertex.Z();
		} else {
			frontx=primarystopvertex.X();
			fronty=primarystopvertex.Y();
			backx=primarystartvertex.X();
			backy=primarystartvertex.Y();
			backz=primarystartvertex.Z();
			trackisbackwardgoing=true;
		}
		//assert(frontx==trackstartx);
		double avgtrackanglex=/*TMath::ATan*/((backx-frontx)/(backz-frontz));
		double avgtrackangley=/*TMath::ATan*/((backy-fronty)/(backz-frontz));
		
		bool drawinitdir=false;
		double mrdentryx,mrdentryy,mrdentryz, mrdexitx, mrdexity, mrdexitz;
		mrdentryx=MRDentrypoint.X();
		mrdentryy=MRDentrypoint.Y();
		mrdentryz=MRDentrypoint.Z();
		mrdexitx=MRDexitpoint.X();
		mrdexity=MRDexitpoint.Y();
		mrdexitz=MRDexitpoint.Z();
#ifdef DRAWTRUEVERBOSE
		cout<<"MRD entry and exit points are:"<<endl
			<<"("<<mrdentryx<<", "<<mrdentryy<<", "<<mrdentryz<<") to "
			<<"("<<mrdexitx<<", "<<mrdexity<<", "<<mrdexitz<<")"<<endl;
#endif
		
////////////////////////////////////////////////////////////////////////////////////
//// REPLACED FROM HERE
//		// these are the true track starting & stopping vertices
//		double trackstartx=nextrack.GetStart(0);
//		double trackstarty=nextrack.GetStart(1);
//		double trackstartz=nextrack.GetStart(2);
//		double trackstopx=nextrack.GetStop(0);
//		double trackstopy=nextrack.GetStop(1);
//		double trackstopz=nextrack.GetStop(2);
//		
//		// Draw vector indicating the true track's intial direction, for better idea of trajectory
//		double dirx = nextrack.GetDir(0);
//		double diry = nextrack.GetDir(1);
//		double dirz = nextrack.GetDir(2);
//		
//		// don't bother drawing short tracks?
//////		double tracklength=TMath::Sqrt(TMath::Power((trackstopx-trackstartx),2)+
//////			TMath::Power((trackstopy-trackstarty),2)+TMath::Power((trackstopz-trackstartz),2));
//////		if(tracklength<20) continue;
//		
//		// for testing:
///////////////////////////////////////////////////////////////////
////		double trackstartx=-MRDSpecs::MRD_width;
////		double trackstarty=-MRDSpecs::MRD_height;
////		double trackstartz=MRDSpecs::MRD_start-10.;
////		double trackstopx=MRDSpecs::MRD_width;
////		double trackstopy=MRDSpecs::MRD_height;
////		double trackstopz=MRDSpecs::MRD_start+MRDSpecs::MRD_depth+10.;
///////////////////////////////////////////////////////////////////
//		
///////////////////////////////////////////////////////////////////
////		double trackstartx, trackstarty, trackstartz, trackstopx, trackstopy, trackstopz;
////		switch (trackit) {
////		case 0:
////			trackstartx=MRDSpecs::MRD_width*1.3;
////			trackstarty=MRDSpecs::MRD_height*1.3;
////			trackstartz=MRDSpecs::MRD_start+(MRDSpecs::MRD_depth/2.)-(MRDSpecs::MRD_depth*(1.3/2.));
////			trackstopx=-MRDSpecs::MRD_width*1.3;
////			trackstopy=-MRDSpecs::MRD_height*1.3;
////			trackstopz=MRDSpecs::MRD_start+(MRDSpecs::MRD_depth/2.)+(MRDSpecs::MRD_depth*(1.3/2.));
////			break;
////		case 1:
////			trackstartx=-MRDSpecs::MRD_width*0.8;
////			trackstarty=-MRDSpecs::MRD_height*0.8;
////			trackstartz=MRDSpecs::MRD_start+(MRDSpecs::MRD_depth/2.)-(MRDSpecs::MRD_depth*0.4);
////			trackstopx=MRDSpecs::MRD_width*0.8;
////			trackstopy=MRDSpecs::MRD_height*0.8;
////			trackstopz=MRDSpecs::MRD_start+(MRDSpecs::MRD_depth/2.)+(MRDSpecs::MRD_depth*0.4);
////			break;
////		}
///////////////////////////////////////////////////////////////////
//		
///////////////////////////////////////////////////////////////////
////		cout<<"origin is ("<<placeholderx.at(trackit)<<", "<<placeholdery.at(trackit)
////			<<", "<<placeholderz.at(trackit)<<")"<<endl;
////		
////		double trackstartx, trackstarty, trackstartz, trackstopx, trackstopy, trackstopz;
////		double scintfullylen = (placeholderorient.at(trackit)) ? scinthfullylen : scintvfullylen;
////		if(placeholderorient.at(trackit)){   // horizontal paddle. these are rotated; need to swap x,y
////			trackstartx=placeholderx.at(trackit);
////			trackstarty=placeholdery.at(trackit);
////		
////			trackstopx=placeholderx2.at(trackit);
////			trackstopy=placeholdery2.at(trackit);
////		} else {
////			trackstartx=placeholderx.at(trackit);
////			trackstarty=placeholdery.at(trackit);
////		
////			trackstopx=placeholderx2.at(trackit);
////			trackstopy=placeholdery2.at(trackit);
////		}
////		trackstartz=placeholderz.at(trackit) + MRDSpecs::MRD_start + (MRDSpecs::MRD_depth/2.);
////		trackstopz= placeholderz2.at(trackit) + MRDSpecs::MRD_start + (MRDSpecs::MRD_depth/2.);
///////////////////////////////////////////////////////////////////
//		
///////////////////////////////////////////////////////////////////
////		double trackstartx=0;
////		double trackstarty=0;
////		double trackstartz=MRDSpecs::MRD_start-10.;
////		double trackstopx=0;
////		double trackstopy=0;
////		double trackstopz=MRDSpecs::MRD_start+MRDSpecs::MRD_depth+10.;
///////////////////////////////////////////////////////////////////
//		
//		// skip the easy ones
//		if(trackstopz<MRDSpecs::MRD_start){
//#ifdef DRAWTRUEVERBOSE
//			cout<<"skipping track "<<trackit<<" as it doesn't make it to mrd: trackstopz="
//				<<trackstopz<<", MRDSpecs::MRD_start="<<MRDSpecs::MRD_start<<endl;
//#endif
//			continue;
//		}  // the track never makes it to the MRD
//		if(trackstartz>(MRDSpecs::MRD_start+MRDSpecs::MRD_depth)){
//#ifdef DRAWTRUEVERBOSE
//			cout<<"skipping track "<<trackit<<" as it starts after mrd: trackstartz="
//				<<trackstartz<<", MRDSpecs::MRD_start+MRDSpecs::MRD_depth="<<MRDSpecs::MRD_start+MRDSpecs::MRD_depth<<endl;
//#endif
//			continue;
//		} // the track is somehow created after the MRD

//#ifdef DRAWTRUEVERBOSE
//		cout<<"Track start: ("<<trackstartx<<", "<<trackstarty<<", "<<trackstartz
//			<<"), track end: ("<<trackstopx<<", "<<trackstopy<<", "<<trackstopz<<")"<<endl;
//#endif
//		
//		// tracks may be forward or backward going: get the front and back, irrespective of direction
//		// while doing this, remove the MRD offset
//		double frontz=min(trackstartz,trackstopz);
//		double frontx,fronty,backz,backx,backy;
//		bool trackisbackwardgoing=false;
//		if(frontz==trackstartz){
//			frontx=trackstartx;
//			fronty=trackstarty;
//			backx=trackstopx;
//			backy=trackstopy;
//			backz=trackstopz;
//		} else {
//			frontx=trackstopx;
//			fronty=trackstopy;
//			backx=trackstartx;
//			backy=trackstarty;
//			backz=trackstartz;
//			trackisbackwardgoing=true;
//		}
//		//assert(frontx==trackstartx);
//		double avgtrackanglex=/*TMath::ATan*/((backx-frontx)/(backz-frontz));
//		double avgtrackangley=/*TMath::ATan*/((backy-fronty)/(backz-frontz));
//		
//		bool drawinitdir=false;
//		double mrdentryx,mrdentryy,mrdentryz, mrdexitx, mrdexity, mrdexitz;
//		// we have 2 cases: a track starts before the MRD or in the MRD
//		bool sidentry=false;
//		if(frontz<MRDSpecs::MRD_start){  // case 1: track starts before the MRD.
//#ifdef DRAWTRUEVERBOSE
//			cout<<"track starts before the MRD";
//#endif
//			// project the track to the z plane of the MRD start
//			mrdentryx=frontx+(MRDSpecs::MRD_start-frontz)*/*TMath::Tan*/(avgtrackanglex);
//			mrdentryy=fronty+(MRDSpecs::MRD_start-frontz)*/*TMath::Tan*/(avgtrackangley);
//			if((abs(mrdentryx)<MRDSpecs::MRD_width)&&(abs(mrdentryy)<MRDSpecs::MRD_height)){
//				mrdentryz=MRDSpecs::MRD_start;
//#ifdef DRAWTRUEVERBOSE
//				cout<<" and projects within X and Y bounds - zstart is MRDSpecs::MRD_start"<<endl;
//#endif
//			} else {
//				sidentry=true;  // projecting to the front face doesn't put it within MRD bounds
//#ifdef DRAWTRUEVERBOSE
//				cout<<" but projects outside X and Y bounds - side entry"<<endl;
//#endif
//			}
//		}
//		if((!(frontz<MRDSpecs::MRD_start)) || sidentry){           // case 2: z entry is within MRD bounds. 
//			// first the simple case: the track start was within the MRD
//			if((abs(frontx)<MRDSpecs::MRD_width)&&(abs(fronty)<MRDSpecs::MRD_height)){
//#ifdef DRAWTRUEVERBOSE
//				cout<<"track starts within the MRD"<<endl;
//#endif
//				mrdentryx=frontx; mrdentryy=fronty; mrdentryz=frontz;
//				drawinitdir=true;
//			} else {
//#ifdef DRAWTRUEVERBOSE
//				cout<<"side entry"<<endl;
//#endif
//				// the track must have entered from a side.
//				// we need to find the z at which the projected track enters both x and y bounds
//				// this is the greater of projected z entry point for x or y.
//				double projectedentryzx, projectedentryzy;
//				if(abs(frontx)<MRDSpecs::MRD_width){
//#ifdef DRAWTRUEVERBOSE
//					cout<<"start x is within bounds"<<endl;
//#endif
//					projectedentryzx = -1.; // it's already within bounds
//				} else {
//					double traverselength = (frontx>0) ? (frontx-MRDSpecs::MRD_width) : (frontx+MRDSpecs::MRD_width);
//					projectedentryzx = frontz - traverselength/ /*TMath::Tan*/(avgtrackanglex);
//#ifdef DRAWTRUEVERBOSE
//					cout<<"using signed projectedentryzx="<<projectedentryzx
//						<<", using unsigned ="<<(frontz+((abs(frontx)-MRDSpecs::MRD_width)/abs(avgtrackanglex)))<<endl;
//					//TODO i forget what this is about but they aren't always the same.
//					cout<<"projected x entry point is at z="<<projectedentryzx<<endl;
//#endif
//				}
//				if(abs(fronty)<MRDSpecs::MRD_height){
//#ifdef DRAWTRUEVERBOSE
//					cout<<"start y is within bounds"<<endl;
//#endif
//					projectedentryzy = -1.;
//				} else {
//					double traverselength = (fronty>0) ? (fronty-MRDSpecs::MRD_height) : (fronty+MRDSpecs::MRD_height);
//					projectedentryzy = frontz - traverselength/ /*TMath::Tan*/(avgtrackangley);
//#ifdef DRAWTRUEVERBOSE
//					cout<<"using signed projectedentryzy="<<projectedentryzy
//						<<", using unsigned ="<<(frontz+((abs(fronty)-MRDSpecs::MRD_height)/abs(avgtrackangley)))<<endl;
//					cout<<"projected y entry point is at z="<<projectedentryzy<<endl;
//#endif
//				}
//				mrdentryz=max(projectedentryzy,projectedentryzx);
//				mrdentryx=frontx+(mrdentryz-frontz)*/*TMath::Tan*/(avgtrackanglex);
//				mrdentryy=fronty+(mrdentryz-frontz)*/*TMath::Tan*/(avgtrackangley);
//			}
//		}
//#ifdef DRAWTRUEVERBOSE
//		cout<<"mrd entry point is ("<<mrdentryx<<", "<<mrdentryy<<", "<<mrdentryz<<")"<<endl;
//#endif
//		
//		// repeat the logic for the exit point
//		// we have 2 cases: a track ends after the MRD or in the MRD
//		bool sidexit=false;
//		if(backz>(MRDSpecs::MRD_start+MRDSpecs::MRD_depth)){  // case 1: track ends after the MRD.
//#ifdef DRAWTRUEVERBOSE
//			cout<<"track ends after the MRD";
//#endif
//			// project the track to the z plane of the MRD end
//			mrdexitx=frontx+(MRDSpecs::MRD_start+MRDSpecs::MRD_depth-frontz)*/*TMath::Tan*/(avgtrackanglex);
//			mrdexity=fronty+(MRDSpecs::MRD_start+MRDSpecs::MRD_depth-frontz)*/*TMath::Tan*/(avgtrackangley);
//			if((abs(mrdexitx)<MRDSpecs::MRD_width)&&(abs(mrdexity)<MRDSpecs::MRD_height)){
//#ifdef DRAWTRUEVERBOSE
//				cout<<" and within MRD bounds - mrdzexit = MRDSpecs::MRD_end"<<endl;
//#endif
//				mrdexitz=MRDSpecs::MRD_start+MRDSpecs::MRD_depth;
//			} else {
//				sidexit=true;  // projecting to the back face doesn't put it within MRD bounds
//#ifdef DRAWTRUEVERBOSE
//				cout<<" but projects outside X and Y bounds - side exit"<<endl;
//#endif
//			}
//		}
//		if((!(backz>(MRDSpecs::MRD_start+MRDSpecs::MRD_depth)))||sidexit){  // case 2: z exit is within MRD bounds. 
//			// first the simple case: the track end is within the MRD
//			if((abs(backx)<MRDSpecs::MRD_width)&&(abs(backy)<MRDSpecs::MRD_height)){
//#ifdef DRAWTRUEVERBOSE
//				cout<<"track ends within the MRD"<<endl;
//#endif
//				mrdexitx=backx; mrdexity=backy; mrdexitz=backz;
//			} else {
//#ifdef DRAWTRUEVERBOSE
//				cout<<"side exit"<<endl;
//#endif
//				// the track must have exited from a side.
//				// we need to find the z at which the projected track exits either x or y bounds
//				// this is the lesser of projected z exit point for x or y.
//				double projectedexitzx, projectedexitzy;
//				if(abs(backx)<MRDSpecs::MRD_width){
//					projectedexitzx = 999.; // it's already within bounds
//				} else {
//					double traverselength = (backx>0) ? (backx-MRDSpecs::MRD_width) : (backx+MRDSpecs::MRD_width);
//					projectedexitzx = backz - traverselength/ /*TMath::Tan*/(avgtrackanglex);
//#ifdef DRAWTRUEVERBOSE
//					cout<<"using signed projectedexitzx="<<projectedexitzx<<", unsigned ="
//						<<(backz-((abs(backx)-MRDSpecs::MRD_width)/abs(avgtrackanglex)))<<endl;
//#endif
//				}
//				if(abs(backy)<MRDSpecs::MRD_height){
//					projectedexitzy = 999.;
//				} else {
//					double traverselength = (backy>0) ? (backy-MRDSpecs::MRD_height) : (backy+MRDSpecs::MRD_height);
//					projectedexitzy = backz - traverselength/ /*TMath::Tan*/(avgtrackangley);
//#ifdef DRAWTRUEVERBOSE
//					cout<<"using signed projectedexitzy="<<projectedexitzy<<", unsigned ="
//						<<(backz-((abs(backy)-MRDSpecs::MRD_height)/abs(avgtrackangley)))<<endl;
//#endif
//				}
//				mrdexitz=min(projectedexitzy,projectedexitzx);
//				mrdexitx=frontx+(mrdexitz-frontz)*/*TMath::Tan*/(avgtrackanglex);
//				mrdexity=fronty+(mrdexitz-frontz)*/*TMath::Tan*/(avgtrackangley);
//			}
//		}
//#ifdef DRAWTRUEVERBOSE
//		cout<<"mrd exit point is ("<<mrdexitx<<", "<<mrdexity<<", "<<mrdexitz<<")"<<endl;
//#endif
// END OF REPLACEMENT
////////////////////////////////////////////////////////////////////////////////////
		
		if(trackisbackwardgoing){
#ifdef DRAWTRUEVERBOSE
			cout<<"switching entry and exit points as track is backward going"<<endl;
#endif
			double tempx = mrdentryx;
			double tempy = mrdentryy;
			double tempz = mrdentryz;
			mrdentryx = mrdexitx;
			mrdentryy = mrdexity;
			mrdentryz = mrdexitz;
			mrdexitx = tempx;
			mrdexity = tempy;
			mrdexitz = tempz;
		}
		
		// up to now all measurements are in WCSim absolute coordinates. Shift z axis so that
		// the MRD is centered on (0,0);
		mrdentryz -= (MRDSpecs::MRD_start+(MRDSpecs::MRD_depth/2.));
		mrdexitz -= (MRDSpecs::MRD_start+(MRDSpecs::MRD_depth/2.));
#ifdef DRAWTRUEVERBOSE
		cout<<"shifting z axis; new entry and exit points are "<<mrdentryz<<" and "<<mrdexitz<<endl;
		cout<<"entry and exit points in terms of mrd width, height and depth are: ("
			<<(mrdentryx/MRDSpecs::maxwidth)<<", "<<(mrdentryy/MRDSpecs::maxheight)<<", "<<(mrdentryz/MRDSpecs::mrdZlen)<<") -> ("
			<<(mrdexitx/MRDSpecs::maxwidth)<<", "<<(mrdexity/MRDSpecs::maxheight)<<", "<<(mrdexitz/MRDSpecs::mrdZlen)<<")"<<endl;
#endif
		
		//============================================================
		// convert from simulation 'cm' units to canvas units: 0-1, 0-1 horizontal and vertical
		// we also need to know the 'side' (left or right) of the track at it's start, end, and where
		// it crosses if applicable, for both views, so that we can shift the track in both views
		// to account for the paddle shift used to visualize two halves on the same canvas
		double anoffset=(MRDSpecs::scintfullzlen+MRDSpecs::scintalugap)*5.;  // this offset accuonts for the half shift
		
		// add the starting x and y, scaled to canvas size and offset to start of MRD diagram
		/*  ✩ ✨ Magic Numbers! ✨ ✩ */
//		double topscalefactor=1.5;           // compress canvas width to paddle diagram height
//		double sidescalefactor=1.5;          //   "         "      "       "       "    width
//		double topdepthscalefactor=1.35;     //   "         "      "       "       "    depth (top view)
//		double sidedepthscalefactor=1.4;     // compress canvas depth to paddle diagram depth (side view)
//		double xscalefactor=(0.5/0.403825);  // correct differences in definition of MRD width and height
//		double yscalefactor=(0.5/0.384671);  // between this method and that for paddle placements
//		double topzoffset=0.03;              // shifts the track arrows +z          (top  view)
//		double sidezoffset=0.005;            // to account for centering of diagram (side view)

		double topscalefactor=1.5;           // compress canvas width to paddle diagram height
		double sidescalefactor=1.55;         //   "         "      "       "       "    width
		double topdepthscalefactor=1.18;     //   "         "      "       "       "    depth (top view)
		double sidedepthscalefactor=1.2;     // compress canvas depth to paddle diagram depth (side view)
		double xscalefactor=(0.5/0.403825);  // correct differences in definition of MRD width and height
		double yscalefactor=(0.5/0.384671);  // between this method and that for paddle placements
		double topzoffset=0.0;               // shifts the track arrows +z          (top  view)
		double sidezoffset=0.0;              // to account for centering of diagram (side view)
		
		
		mrdentryx*=xscalefactor;
		mrdexitx*=xscalefactor;
		mrdentryy*=yscalefactor;
		mrdexity*=yscalefactor;
#ifdef DRAWTRUEVERBOSE
		cout<<"scaled entry and exit points in terms of mrd width, height and depth are: ("
			<<(mrdentryx/MRDSpecs::maxwidth)<<", "<<(mrdentryy/MRDSpecs::maxheight)<<", "<<(mrdentryz/MRDSpecs::mrdZlen)<<") -> ("
			<<(mrdexitx/MRDSpecs::maxwidth)<<", "<<(mrdexity/MRDSpecs::maxheight)<<", "<<(mrdexitz/MRDSpecs::mrdZlen)<<")"<<endl;
#endif
		
		// one last thing: the beam comes from the left. In the top view, right-hand-side (x>0)
		// needs to map to the bottom of the canvas (canvas_y<0) - so, let's swap the signs of all
		// x points
		mrdentryx*=-1.;
		mrdexitx*=-1.;
		avgtrackanglex*=-1.;
		
		std::vector<double> xstarts, ystarts, zstartsx, zstartsy, xstops, ystops, zstopsx, zstopsy;
		xstarts.push_back((mrdentryx/(MRDSpecs::maxwidth*topscalefactor))+0.5);
		ystarts.push_back((mrdentryy/(MRDSpecs::maxheight*sidescalefactor))+0.5);
		// starting z may need a shift depending on the appropriate half
		// in top view
		if(mrdentryy>0){
			zstartsx.push_back((mrdentryz/(MRDSpecs::mrdZlen*topdepthscalefactor))+0.5+topzoffset);
		} else {
			zstartsx.push_back(((mrdentryz+anoffset)/(MRDSpecs::mrdZlen*topdepthscalefactor))+0.5+topzoffset);
		}
		// in side view
		if(mrdentryx>0){
			zstartsy.push_back((mrdentryz/(MRDSpecs::mrdZlen*sidedepthscalefactor))+0.5+sidezoffset);
		} else {
			zstartsy.push_back(((mrdentryz+anoffset)/(MRDSpecs::mrdZlen*sidedepthscalefactor))+0.5+sidezoffset);
		}
		// check if we cross sides, and if so, create a middle stop and start set
		// top view
		if((mrdentryy*mrdexity)<0){
			// we'll need two lines with a bit of a disconnect. find the crossing point.
			double crossingz = mrdentryz-((mrdentryy/yscalefactor) / /*TMath::Tan*/(avgtrackangley));
			double crossingx = (mrdentryx/xscalefactor) + ((crossingz-mrdentryz)*avgtrackanglex);
			xstops.push_back(((crossingx*xscalefactor)/(MRDSpecs::maxwidth*topscalefactor))+0.5);
			if(mrdentryy>0){
				zstopsx.push_back((crossingz/(MRDSpecs::mrdZlen*topdepthscalefactor))+0.5+topzoffset);
			} else {
				zstopsx.push_back(((crossingz+anoffset)/(MRDSpecs::mrdZlen*topdepthscalefactor))+0.5+topzoffset);
			}
			xstarts.push_back(((crossingx*xscalefactor)/(MRDSpecs::maxwidth*topscalefactor))+0.5);
			if(mrdexity>0){
				zstartsx.push_back((crossingz/(MRDSpecs::mrdZlen*topdepthscalefactor))+0.5+topzoffset);
			} else {
				zstartsx.push_back(((crossingz+anoffset)/(MRDSpecs::mrdZlen*topdepthscalefactor))+0.5+topzoffset);
			}
		}
		// side view
		if((mrdentryx*mrdexitx)<0){
			// we'll need two lines with a bit of a disconnect. find the crossing point.
			double crossingz = mrdentryz-((mrdentryx/xscalefactor) / /*TMath::Tan*/(avgtrackanglex));
			double crossingy = (mrdentryy/yscalefactor) + ((crossingz-mrdentryz)*avgtrackangley);
			ystops.push_back(((crossingy*yscalefactor)/(MRDSpecs::maxheight*sidescalefactor))+0.5);
			if(mrdentryx>0){
				zstopsy.push_back((crossingz/(MRDSpecs::mrdZlen*sidedepthscalefactor))+0.5+sidezoffset);
			} else {
				zstopsy.push_back(((crossingz+anoffset)/(MRDSpecs::mrdZlen*sidedepthscalefactor))+0.5+sidezoffset);
			}
			ystarts.push_back(((crossingy*yscalefactor)/(MRDSpecs::maxheight*sidescalefactor))+0.5);
			if(mrdexitx>0){
				zstartsy.push_back((crossingz/(MRDSpecs::mrdZlen*sidedepthscalefactor))+0.5+sidezoffset);
			} else {
				zstartsy.push_back(((crossingz+anoffset)/(MRDSpecs::mrdZlen*sidedepthscalefactor))+0.5+sidezoffset);
			}
		}
		// finally add the endpoint values, with offset for z according to ending MRD half.
		xstops.push_back((mrdexitx/(MRDSpecs::maxwidth*topscalefactor))+0.5);
		ystops.push_back((mrdexity/(MRDSpecs::maxheight*sidescalefactor))+0.5);
		// top view
		if(mrdexity>0){
			zstopsx.push_back((mrdexitz/(MRDSpecs::mrdZlen*topdepthscalefactor))+0.5+topzoffset);
		} else {
			zstopsx.push_back(((mrdexitz+anoffset)/(MRDSpecs::mrdZlen*topdepthscalefactor))+0.5+topzoffset);
		}
		// side view
		if(mrdexitx>0){
			zstopsy.push_back((mrdexitz/(MRDSpecs::mrdZlen*sidedepthscalefactor))+0.5+sidezoffset);
		} else {
			zstopsy.push_back(((mrdexitz+anoffset)/(MRDSpecs::mrdZlen*sidedepthscalefactor))+0.5+sidezoffset);
		}
		
		// OK done.
		// now loop over the pairs and make the arrows
		// top view
		for(int i=0; i < zstartsx.size(); i++){
			// Draw arrow representing "true" (assumed straight) trajectory in top view
			std::string arrowdir = (!trackisbackwardgoing) ? ">" : "<";
			TArrow* myarrow = 
				new TArrow(zstartsx.at(i), xstarts.at(i), zstopsx.at(i), xstops.at(i), 0.005, arrowdir.c_str());
			myarrow->SetLineWidth(1);
			myarrow->SetLineColor(kRed);
			imgcanvas->cd(2);  // top view for x positions
			myarrow->Draw();
#ifdef DRAWTRUEVERBOSE
			cout<<"drawing top view truth track arrow from "<<myarrow->GetX1()<<", "<<myarrow->GetY1()
				<<" to "<<myarrow->GetX2()<<", "<<myarrow->GetY2()<<endl;
#endif
			truetrackarrows.push_back(myarrow);
			
			if(zstartsx.size()==2&&i==0){
				// add an intermediate link dashed line to link the arrows
				myarrow = 
					new TArrow(zstopsx.at(i), xstops.at(i), zstartsx.at(i+1), xstarts.at(i+1),0.0,">");
					// an arrow size of 0.0 gives no arrow head (just a line)
				myarrow->SetLineWidth(1);
				myarrow->SetLineColor(kRed);
				myarrow->SetLineStyle(7);  //dashed
				myarrow->Draw();
#ifdef DRAWTRUEVERBOSE
				cout<<"drawing line from "<<myarrow->GetX1()<<", "<<myarrow->GetY1()<<" to "
					<<myarrow->GetX2()<<", "<<myarrow->GetY2()<<endl;
#endif
				truetrackarrows.push_back(myarrow);
			}
			
			// draw the arrow indicating the true initial direction
			if(i==0&&drawinitdir){
				TArrow* myarrow = new TArrow( zstartsx.at(i), xstarts.at(i), 
						zstartsx.at(i)+(dirz*10./(MRDSpecs::mrdZlen*topdepthscalefactor)),
						xstarts.at(i)+(diry*10.*xscalefactor/(MRDSpecs::maxwidth*topscalefactor)),
						0.005, arrowdir.c_str() );
				myarrow->SetLineWidth(2);
				myarrow->SetLineColor(kMagenta);
				myarrow->Draw();
				truetrackarrows.push_back(myarrow);
			}
		}
		
		// need to do top and side views separately as they may have different sizes
		// side view
		for(int i=0; i < zstartsy.size(); i++){
			// Draw arrow representing "true" (assumed straight) trajectory in top view
			std::string arrowdir = (!trackisbackwardgoing) ? ">" : "<";
			TArrow* myarrow = 
				new TArrow(zstartsy.at(i), ystarts.at(i), zstopsy.at(i), ystops.at(i), 0.005, arrowdir.c_str());
			myarrow->SetLineWidth(1);
			myarrow->SetLineColor(kRed);
			imgcanvas->cd(1);  // side view for y positions
			myarrow->Draw();
#ifdef DRAWTRUEVERBOSE
			cout<<"drawing side view truth track arrow from "<<myarrow->GetX1()<<", "
			<<myarrow->GetY1()<<" to "<<myarrow->GetX2()<<", "<<myarrow->GetY2()<<endl;
#endif
			truetrackarrows.push_back(myarrow);
			
			if(zstartsy.size()==2&&i==0){
				// add an intermediate link dashed line to link the arrows
				myarrow = 
					new TArrow(zstopsy.at(i), ystops.at(i), zstartsy.at(i+1), ystarts.at(i+1),0.0,">");
					// an arrow size of 0.0 gives no arrow head (just a line)
				myarrow->SetLineWidth(1);
				myarrow->SetLineColor(kRed);
				myarrow->SetLineStyle(7);  //dashed
				myarrow->Draw();
#ifdef DRAWTRUEVERBOSE
				cout<<"drawing line from "<<myarrow->GetX1()<<", "<<myarrow->GetY1()<<" to "
					<<myarrow->GetX2()<<", "<<myarrow->GetY2()<<endl;
#endif
				truetrackarrows.push_back(myarrow);
			}
			
			// draw the arrow indicating the true initial direction
			if(i==0&&drawinitdir){
				TArrow* myarrow = new TArrow( zstartsy.at(i), ystarts.at(i), 
						zstartsy.at(i)+(dirz*10./(MRDSpecs::mrdZlen*sidedepthscalefactor)),
						ystarts.at(i)+(diry*10.*yscalefactor/(MRDSpecs::maxheight*sidescalefactor)),
						0.005, arrowdir.c_str() );
				myarrow->SetLineWidth(2);
				myarrow->SetLineColor(kMagenta);
				myarrow->Draw();
				truetrackarrows.push_back(myarrow);
			}
		}
		trackit++;
	}
}

//============================================================================

// a test to see if a projected point in a plane is within a box in that plane
int inline InBox( TVector3 Hit, TVector3 B1, TVector3 B2, const int Axis) {
	if ( Axis==1 && Hit.Z() > B1.Z() && Hit.Z() < B2.Z() && Hit.Y() > B1.Y() && Hit.Y() < B2.Y()) return 1;
	if ( Axis==2 && Hit.Z() > B1.Z() && Hit.Z() < B2.Z() && Hit.X() > B1.X() && Hit.X() < B2.X()) return 1;
	if ( Axis==3 && Hit.X() > B1.X() && Hit.X() < B2.X() && Hit.Y() > B1.Y() && Hit.Y() < B2.Y()) return 1;
	return 0;
}

// projects the hitpoint by adding a scaled vector to the start point
int inline GetIntersection( float fDst1, float fDst2, TVector3 P1, TVector3 P2, TVector3 &Hit) {
	if ( (fDst1 * fDst2) >= 0.0f) return 0;
	if ( fDst1 == fDst2) return 0; 
	Hit = P1 + (P2-P1) * ( -fDst1/(fDst2-fDst1) );
	return 1;
}

// returns true if line (L1, L2) intersects with the box (B1, B2)
// returns intersection point in Hit
bool CheckLineBox( TVector3 L1, TVector3 L2, TVector3 B1, TVector3 B2, TVector3 &Hit, TVector3 &Hit2, bool &error){
	error=false;
	bool startsinmrd=false, stopsinmrd=false;
	// check if it misses the box entirely by being on one side of a plane over entire track
	if (L2.X() < B1.X() && L1.X() < B1.X()) return false;
	if (L2.X() > B2.X() && L1.X() > B2.X()) return false;
	if (L2.Y() < B1.Y() && L1.Y() < B1.Y()) return false;
	if (L2.Y() > B2.Y() && L1.Y() > B2.Y()) return false;
	if (L2.Z() < B1.Z() && L1.Z() < B1.Z()) return false;
	if (L2.Z() > B2.Z() && L1.Z() > B2.Z()) return false;
	// check if it's inside the box to begin with (classed as an interception at start vtx)
	if (L1.X() > B1.X() && L1.X() < B2.X() &&
		L1.Y() > B1.Y() && L1.Y() < B2.Y() &&
		L1.Z() > B1.Z() && L1.Z() < B2.Z())
		{startsinmrd=true; /*Hit = L1; return true;*/}
		// XXX DISABLE RETURN FOR TRACKS STARTING IN MRD - NEED TO FIND EXIT POINT
	// check if it's inside the box to end with (classed as an interception at stop vtx)
	if (L2.X() > B1.X() && L2.X() < B2.X() &&
		L2.Y() > B1.Y() && L2.Y() < B2.Y() &&
		L2.Z() > B1.Z() && L2.Z() < B2.Z())
		{stopsinmrd=true; /*Hit2 = L2; return true;*/}
	if(startsinmrd&&stopsinmrd){
		// no point looking for interceptions
		Hit=L1; Hit2=L2; return true;
	}
	// going forward, either one or both track endpoints are outside the MRD, and the MRD is hit.
	// we need to find the interception(s) with the MRD volumetric surface.

	// check for an interception in X, Y then Z.
	//if ( (GetIntersection( L1.X()-B1.X(), L2.X()-B1.X(), L1, L2, Hit) && InBox( Hit, B1, B2, 1 ))
	//  || (GetIntersection( L1.Y()-B1.Y(), L2.Y()-B1.Y(), L1, L2, Hit) && InBox( Hit, B1, B2, 2 ))
	//  || (GetIntersection( L1.Z()-B1.Z(), L2.Z()-B1.Z(), L1, L2, Hit) && InBox( Hit, B1, B2, 3 ))
	//  || (GetIntersection( L1.X()-B2.X(), L2.X()-B2.X(), L1, L2, Hit) && InBox( Hit, B1, B2, 1 ))
	//  || (GetIntersection( L1.Y()-B2.Y(), L2.Y()-B2.Y(), L1, L2, Hit) && InBox( Hit, B1, B2, 2 ))
	//  || (GetIntersection( L1.Z()-B2.Z(), L2.Z()-B2.Z(), L1, L2, Hit) && InBox( Hit, B1, B2, 3 )))
	//	return true;

	// Above seems to assume there will only be one interception!!
	// e.g. if X has an interception, there are no checks for Z interception - if it enters
	// the front face and exits the side, only the side exit will be returned. 
	// Instead, note all interception points and return the first (and second if it exists)
	std::vector<TVector3> interceptions;
	bool anyinterception=false;
	bool thisinterception;

	thisinterception=
	GetIntersection( L1.X()-B1.X(), L2.X()-B1.X(), L1, L2, Hit) && InBox(Hit, B1, B2, 1);
	if(thisinterception){ interceptions.push_back(Hit); anyinterception=true; }
	thisinterception=
	GetIntersection( L1.Y()-B1.Y(), L2.Y()-B1.Y(), L1, L2, Hit) && InBox( Hit, B1, B2, 2 );
	if(thisinterception){ interceptions.push_back(Hit); anyinterception=true; }
	thisinterception=
	GetIntersection( L1.Z()-B1.Z(), L2.Z()-B1.Z(), L1, L2, Hit) && InBox( Hit, B1, B2, 3 );
	if(thisinterception){ interceptions.push_back(Hit); anyinterception=true; }
	thisinterception=
	GetIntersection( L1.X()-B2.X(), L2.X()-B2.X(), L1, L2, Hit) && InBox( Hit, B1, B2, 1 );
	if(thisinterception){ interceptions.push_back(Hit); anyinterception=true; }
	thisinterception=
	GetIntersection( L1.Y()-B2.Y(), L2.Y()-B2.Y(), L1, L2, Hit) && InBox( Hit, B1, B2, 2 );
	if(thisinterception){ interceptions.push_back(Hit); anyinterception=true; }
	thisinterception=
	GetIntersection( L1.Z()-B2.Z(), L2.Z()-B2.Z(), L1, L2, Hit) && InBox( Hit, B1, B2, 3 );
	if(thisinterception){ interceptions.push_back(Hit); anyinterception=true; }

	if((interceptions.size()>2)||(startsinmrd&&interceptions.size()==2)){
		cerr<<"CheckLineBox found more than two intercepts?! Or 2 intercepts but starts in MRD?? Intercepts at:"<<endl;
		for(auto&& avec : interceptions)
			cerr<<"("<<avec.X()<<", "<<avec.Y()<<", "<<avec.Z()<<")"<<endl;
		error=true;
		//assert(false); // leave for later so we can print debug info.
		return false;
	} else if(interceptions.size()==2){
		auto vec1 = interceptions.at(0);
		auto vec2 = interceptions.at(1);
		if(vec1.Z()<vec2.Z()){
			Hit=vec1;
			Hit2=vec2;
		} else {
			Hit=vec2;
			Hit2=vec1;
		}
		return true;
	} else if(interceptions.size()==1) {
		if(stopsinmrd){
			Hit=interceptions.at(0);
			Hit2=L2; // Hit2 is 'mrd exit' point - return track end. 
		} else if (startsinmrd){
			Hit = L1;
			Hit=interceptions.at(0);
		} else {
			cerr<<"only one intercept found, but track neither starts nor stops in mrd??"<<endl;
			error=true;
			return false;
		}
		return true;
	} else {
		return false;
	}
}

#endif
