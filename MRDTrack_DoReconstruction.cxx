/* vim:set noexpandtab tabstop=4 wrap */
// =====================================================================
// ********* trajectory reconstruction: main *****************
/* information received during construction */
// event information: MRDtrackID, wcsimfile, run_id, event_id, subtrigger, digi_ids, 
// digit information: pmts_hit, digi_qs, digi_ts, 
// truth information: digi_numphots, digi_phot_ts, digi_phot_parents
/* make functions to retrieve these? */
// digits(), trueTrackID(-1)
// ------------------------
// 
/* information to calculate here: */
// TODO: KEStart(), KEEnd(), particlePID(), tracktype() (one track, 2 tracks, inconsistent hits...),
// Done: layers_hit(), eDepsInLayers(), recovalshoriz(), recovalsvert()
// ^^^^^^^^^^^^^^^^^^^^^^^^

#ifndef MRDTrack_RECO_VERBOSE
#define MRDTrack_RECO_VERBOSE
#endif

void cMRDTrack::DoReconstruction(){
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"doing track reconstruction"<<endl;
#endif
	/* 
	Reconstruction generates two maps (recovalshoriz, recovalsvert) each with keys: 
	(xycrossing, zcrossing, angmax, angmin)
	These values define the region of acceptance in that plane - projecting from the crossing points,
	the angles define a wedge in the x-z / y-z plane.
	This wedge, projected back to the tank, defines the region of possible starting points and angles
	within the tank 
	Tracks starting from a point in this region must have the appropriate angle to hit the crossing vertex.
	The two sets of constraints in x and y are independent and define a non-trivial geometry when combined.*/
	
	int numdigits = digi_ids.size();
	std::vector<TVector3> bottompoints;
	std::vector<TVector3> toppoints;
	// sort the hits in digi_ids, pmts_hit (and correspondingly layers_hit) by tube id:
	// useful for debugging but shouldn't be necessary for functionality - better not to.
	std::vector<std::pair<int,int>> digisandtubeids;
	for(int i=0; i<numdigits; i++){
		digisandtubeids.emplace_back(pmts_hit.at(i),i);
	}
	std::sort(digisandtubeids.begin(), digisandtubeids.end()); // sorts by first of pair by default
	std::vector<int> digi_ids_dup(digi_ids), pmts_hit_dup(pmts_hit);
	for(int i=0; i<numdigits; i++){
		int sortedindex = digisandtubeids.at(i).second;
		pmts_hit.at(i) = pmts_hit_dup.at(sortedindex);
		digi_ids.at(i) = digi_ids_dup.at(sortedindex);
	}
	
	// FIXME: CURRENTLY NOT WORKING: implemented counting number of hit layers to do a 'best fit'
	// selection for events where no straight line will hit all paddles in the track.
	// but projection signing(?) is not working - when using MRD_muon_sample, turn up verbosity and see
	// event 0.
	// But this method is being abandoned for a better chi2 fit method using the cluster extents
	// as errors.
	
	for(Int_t i=0;i<numdigits;i++){
		// get the extent of the corresponding paddle
		Int_t tube_id = pmts_hit.at(digi_ids.at(i));
#ifdef MRDTrack_RECO_VERBOSE
		cout<<"getting ranges for tube "<<tube_id<<endl;
#endif
		std::pair<Double_t, Double_t> thexrange = mrdcluster::paddle_extentsx.at(tube_id);
		std::pair<Double_t, Double_t> theyrange = mrdcluster::paddle_extentsy.at(tube_id);
		std::pair<Double_t, Double_t> thezrange = mrdcluster::paddle_extentsz.at(tube_id);
		
		// bottompoints stores most negative values: smallest z, most negative x, y. toppoints the opposite.
		bottompoints.push_back(TVector3(thexrange.first, theyrange.first, thezrange.first));
		toppoints.push_back(TVector3(thexrange.second, theyrange.second, thezrange.second));
#ifdef MRDTrack_RECO_VERBOSE
		cout<<"adding bottom point ("<<thexrange.first<<", "<<theyrange.first<<", "<<thezrange.first<<")"<<endl;
		cout<<"and top point ("<<thexrange.second<<", "<<theyrange.second<<", "<<thezrange.second<<")"<<endl;
#endif
		Int_t strucklayer = mrdcluster::paddle_layers.at(tube_id);
		if(std::count(layers_hit.begin(), layers_hit.end(), strucklayer)==0){
			layers_hit.push_back(strucklayer);
		}
		eDepsInLayers.at(strucklayer)+= (digi_qs.at(i));	// TODO: convert q to energy?
	}
	Double_t xangmax = (0.0/0.0);	// angmax represents angle of uppermost limit of valid trajectories
	Double_t xangmin = (0.0/0.0);	// angmin represents angle of lowermost limit of valid trajectories 
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
	int maxmatches=0; // the most layers hit by one trajectory so far
	for(int startlayer=0; startlayer<numdigits; startlayer++){
		for(int nextlayer=(startlayer+1); nextlayer<numdigits; nextlayer++){
			Double_t opp = toppoints.at(nextlayer).X()-bottompoints.at(startlayer).X();
			Double_t adj = bottompoints.at(nextlayer).Z()-bottompoints.at(startlayer).Z();
			// really it should be an average for Z points
			Double_t ang = TMath::ATan(opp/adj);
#ifdef MRDTrack_RECO_VERBOSE
			cout<<"checking most upgoing angle from bottom of layer "
				<<layers_hit.at(startlayer)<<" at x="<<bottompoints.at(startlayer).X()
				<<", z="<<bottompoints.at(startlayer).Z()<<", to top of layer "
				<<layers_hit.at(nextlayer)<<" at x="<<toppoints.at(nextlayer).X()
				<<", z="<<bottompoints.at(nextlayer).Z()
				<<", at angle= "<<TMath::RadToDeg()*ang<<"degs"<<endl;
#endif
/*
			if(TMath::IsNaN(xangmin)||(ang>xangmin)){
				if(AngValid(startlayer,ang,0, 0)){	// check this angle actually hits all paddles
#ifdef MRDTrack_RECO_VERBOSE
					cout<<"setting this x angle as the new most upgoing"<<endl;
#endif
					xangmin=ang;					// set this as the lower bound angle
					xlayermin=startlayer;			// to define a trajectory we also need a starting point
				}
			}
*/
			if(TMath::IsNaN(xangmin)||(ang>xangmin)||(maxmatches<(numdigits-1))){
				// instead count for how many paddles this trajectory is valid:
				int nummatches = AngMatchCount(startlayer,ang,0,0);
				if( (nummatches>maxmatches) ||
					( (nummatches==maxmatches)&&(ang>xangmin) ) ){
#ifdef MRDTrack_RECO_VERBOSE
					cout<<"setting this x angle as the new most upgoing"<<endl;
#endif
					xangmin=ang;					// set this as the lower bound angle
					xlayermin=startlayer;			// to define a trajectory we also need a starting point
				}
				if(nummatches>maxmatches) maxmatches=nummatches;
			}
		}
	}
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"finished calculating most upgoing x angle"<<endl;
#endif
	if(TMath::IsNaN(xangmin)){cout<<" ####### NO UPGOING X ANGLE FOUND ########## "<<endl; return;}
	
	// angmax: this time scan from top of each layer to bottom of next layer, looking 
	// for the uppermost entering trajectory
	maxmatches=0;
	for(Int_t startlayer=0;startlayer<numdigits;startlayer++){
		for(Int_t nextlayer=(startlayer+1);nextlayer<numdigits; nextlayer++){
			Double_t opp = toppoints.at(startlayer).X() - bottompoints.at(nextlayer).X();
			// negative if upgoing
			Double_t adj = bottompoints.at(nextlayer).Z()-bottompoints.at(startlayer).Z();
			Double_t ang = TMath::ATan(opp/adj);
			// negative if downgoing
#ifdef MRDTrack_RECO_VERBOSE
			cout<<"checking most downgoing angle from top of layer "
				<<layers_hit.at(startlayer)<<" at x="<<toppoints.at(startlayer).X()
				<<", z="<<bottompoints.at(startlayer).Z()<<" to bottom of layer "
				<<layers_hit.at(nextlayer)<<", at x="<<bottompoints.at(nextlayer).X()
				<<", z="<<bottompoints.at(nextlayer).Z()
				<<", at angle= "<<TMath::RadToDeg()*ang<<endl;
#endif
/*
			if(TMath::IsNaN(xangmax)||(ang>xangmax)){
				if(AngValid(startlayer,ang,1, 0)){		// project both forward and back from 
					xangmax=ang;						// bottom vertex of nextlayer to all other layers 
					xlayermax=startlayer;				// and check within limits
#ifdef MRDTrack_RECO_VERBOSE
					cout<<"setting this x angle as the new most downgoing"<<endl;
#endif
				}
			}
*/
			if(TMath::IsNaN(xangmax)||(ang>xangmax)||(maxmatches<(numdigits-1))){
				// instead count for how many paddles this trajectory is valid:
				int nummatches = AngMatchCount(startlayer,ang,1,0);
				if( (nummatches>maxmatches) ||
					( (nummatches==maxmatches)&&(ang>xangmax) ) ){
#ifdef MRDTrack_RECO_VERBOSE
					cout<<"setting this x angle as the new most downgoing"<<endl;
#endif
					xangmax=ang;					// set this as the lower bound angle
					xlayermax=startlayer;			// to define a trajectory we also need a starting point
				}
				if(nummatches>maxmatches) maxmatches=nummatches;
			}
		}
	}
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"finished calculating most downgoing x angle"<<endl;
#endif
	if(TMath::IsNaN(xangmax)){cout<<" ####### NO DOWNGOING X ANGLE FOUND ########## "<<endl; return;}
	
	//now find angles and layers for y projections
	maxmatches=0;
	for(Int_t startlayer=0;startlayer<numdigits;startlayer++){
		for(Int_t nextlayer=startlayer+1;nextlayer<numdigits; nextlayer++){
			Double_t opp = toppoints.at(nextlayer).Y()-bottompoints.at(startlayer).Y();
			Double_t adj = bottompoints.at(nextlayer).Z()-bottompoints.at(startlayer).Z();
			// really it should be an average for Z points
			Double_t ang = TMath::ATan(opp/adj);
#ifdef MRDTrack_RECO_VERBOSE
			cout<<"checking most upgoing angle from bottom of layer "
				<<layers_hit.at(startlayer)<<" at y="<<bottompoints.at(startlayer).Y()
				<<", z="<<bottompoints.at(startlayer).Z()<<", to top of layer "
				<<layers_hit.at(nextlayer)<<" at y="<<toppoints.at(nextlayer).Y()
				<<", z="<<bottompoints.at(nextlayer).Z()
				<<", at angle= "<<TMath::RadToDeg()*ang<<"degs"<<endl;
#endif
/*
			if(TMath::IsNaN(yangmin)||(ang>yangmin)){
				if(AngValid(startlayer,ang,0, 1)){			// check this angle actually hits all paddles
					yangmin=ang;							// set this as the lower bound angle
					ylayermin=startlayer;					// a trajectory also needs a starting point
#ifdef MRDTrack_RECO_VERBOSE
					cout<<"setting this y angle as the new most upgoing"<<endl;
#endif
				}
			}
*/
			if(TMath::IsNaN(yangmin)||(ang>yangmin)||(maxmatches<(numdigits-1))){
				// instead count for how many paddles this trajectory is valid:
				int nummatches = AngMatchCount(startlayer,ang,0,1);
				if( (nummatches>maxmatches) ||
					( (nummatches==maxmatches)&&(ang>yangmin) ) ){
#ifdef MRDTrack_RECO_VERBOSE
					cout<<"setting this y angle as the new most upgoing"<<endl;
#endif
					yangmin=ang;					// set this as the lower bound angle
					ylayermin=startlayer;			// to define a trajectory we also need a starting point
				}
				if(nummatches>maxmatches) maxmatches=nummatches;
			}
		}
	}
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"finished calculating most upgoing y angle"<<endl;
#endif
	if(TMath::IsNaN(yangmin)){cout<<" ####### NO UPGOING Y ANGLE FOUND ########## "<<endl; return;}
	
	//now find angles and layers for y projections
	maxmatches=0;
	for(Int_t startlayer=0;startlayer<numdigits;startlayer++){
		for(Int_t nextlayer=startlayer+1;nextlayer<numdigits; nextlayer++){
			// angmax: same, this time scan from top of each layer 
			// to bottom of next layer, for the uppermost entering trajectory
			Double_t opp = toppoints.at(startlayer).Y() - bottompoints.at(nextlayer).Y();
			// negative if upgoing
			Double_t adj = bottompoints.at(nextlayer).Z()-bottompoints.at(startlayer).Z();
			Double_t ang = TMath::ATan(opp/adj);
			// negative if upgoing
#ifdef MRDTrack_RECO_VERBOSE
			cout<<"checking most downgoing angle from top of layer "
				<<layers_hit.at(startlayer)<<" at y="<<toppoints.at(startlayer).Y()
				<<", z="<<bottompoints.at(startlayer).Z()<<" to bottom of layer "
				<<layers_hit.at(nextlayer)<<", at y="<<bottompoints.at(nextlayer).Y()
				<<", z="<<bottompoints.at(nextlayer).Z()
				<<", at angle= "<<TMath::RadToDeg()*ang<<endl;
#endif
/*
			if(TMath::IsNaN(yangmax)||(ang>yangmax)){
				if(AngValid(startlayer,ang,1, 1)){		// project both forward and back from
					yangmax=ang;						// bottom vertex of nextlayer to all other layers
					ylayermax=startlayer;				// and check within limits
#ifdef MRDTrack_RECO_VERBOSE
					cout<<"setting this y angle as the new maximum"<<endl;
#endif
				}
			}
*/
			if(TMath::IsNaN(yangmax)||(ang>yangmax)||(maxmatches<(numdigits-1))){
				// instead count for how many paddles this trajectory is valid:
				int nummatches = AngMatchCount(startlayer,ang,1,1);
				if( (nummatches>maxmatches) ||
					( (nummatches==maxmatches)&&(ang>yangmax) ) ){
#ifdef MRDTrack_RECO_VERBOSE
					cout<<"setting this y angle as the new most downgoing"<<endl;
#endif
					yangmax=ang;					// set this as the lower bound angle
					ylayermax=startlayer;			// to define a trajectory we also need a starting point
				}
				if(nummatches>maxmatches) maxmatches=nummatches;
			}
		}
	}
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"finished calculating most downgoing y angle"<<endl;
#endif
	if(TMath::IsNaN(yangmax)){cout<<" ####### NO DOWNGOING Y ANGLE FOUND ########## "<<endl; return;}

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

// =====================================================================
// ********* trajectory reconstruction: scan projections to all layers *****************
Bool_t cMRDTrack::AngValid(Int_t layerstart, Double_t angle, Int_t MaxMin, Int_t axis){
#ifdef MRDTrack_RECO_VERBOSE
//	cout<<"checking if ";
//	if(axis){cout<<"y ";} else {cout<<"x ";}
//	cout<<"trajectory, going ";
//	if(MaxMin){cout<<"downwards ";} else {cout<<"upwards ";}
//	cout<<"from layer "<<layerstart<<" at angle "<<TMath::RadToDeg()*angle<<" hits other layers"<<endl;
#endif
	Bool_t angvalid=true;
	for(Int_t layerend=0;layerend<digi_ids.size();layerend++){
		if(layerstart==layerend){continue;}	// don't check from a layer to itself
		angvalid = CheckIntersection(layerstart, layerend, angle, MaxMin, axis);
		if(!angvalid){
#ifdef MRDTrack_RECO_VERBOSE
			cout<<"trajectory does not hit layer "<<layerend<<endl;
#endif
			return false;
		}
	}
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"trajectory hits all layers"<<endl;
	if(axis){cout<<"y ";} else {cout<<"x ";}
	cout<<"trajectory, going ";
	if(MaxMin){cout<<"downwards ";} else {cout<<"upwards ";}
	cout<<"from layer "<<layerstart<<" at angle "<<TMath::RadToDeg()*angle<<" hits all other layers"<<endl;
#endif
	return true;
}

// =====================================================================
// ********* trajectory reconstruction: scan projections to all layers *****************
Int_t cMRDTrack::AngMatchCount(Int_t layerstart, Double_t angle, Int_t MaxMin, Int_t axis){
#ifdef MRDTrack_RECO_VERBOSE
//	cout<<"checking if ";
//	if(axis){cout<<"y ";} else {cout<<"x ";}
//	cout<<"trajectory, going ";
//	if(MaxMin){cout<<"downwards ";} else {cout<<"upwards ";}
//	cout<<"from layer "<<layerstart<<" at angle "<<TMath::RadToDeg()*angle<<" hits other layers"<<endl;
#endif
	Int_t matchcount=0;
	for(Int_t layerend=0;layerend<digi_ids.size();layerend++){
		if(layerstart==layerend){continue;}  // don't check from a layer to itself
		bool angvalid = CheckIntersection(layerstart, layerend, angle, MaxMin, axis);
		if(angvalid){ matchcount++;
#ifdef MRDTrack_RECO_VERBOSE
			cout<<"trajectory hits layer "<<layerend<<endl;
#endif
		}
	}
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"trajectory hits "<<matchcount<<" layers of "<<digi_ids.size()<<endl;
	if(matchcount==(digi_ids.size()-1)){
		if(axis){cout<<"y ";} else {cout<<"x ";}
		cout<<"trajectory, going ";
		if(MaxMin){cout<<"downwards ";} else {cout<<"upwards ";}
		cout<<"from layer "<<layerstart<<" at angle "<<TMath::RadToDeg()*angle
			<<" hits all other layers"<<endl;
	}
#endif
	return matchcount;
}

// =====================================================================
// ********* trajectory reconstruction: check single projection *****************
Bool_t cMRDTrack::CheckIntersection(Int_t layerstart, Int_t layerend, Double_t angle, Int_t MaxMin, Int_t axis){

	Double_t z1 = mrdcluster::paddle_extentsz.at(pmts_hit.at(layerstart)).first;
	Double_t z2 = mrdcluster::paddle_extentsz.at(pmts_hit.at(layerend)).first;
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
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"checking if projection from zstart="<<z1;
	if(axis){
		cout<<", ystart=";
		if(MaxMin==1){
			cout<<mrdcluster::paddle_extentsy.at(pmts_hit.at(layerstart)).second; 
			if(z2<z1){cout<<", going backwards and upwards to zend=";} else {cout<<", going downwards to zend=";}
		} else {
			cout<<mrdcluster::paddle_extentsy.at(pmts_hit.at(layerstart)).first; 
			if(z2<z1){cout<<" going backwards and downwards to zend=";} else {cout<<", going upwards to zend=";}
		}
		cout<<z2<<" at angle "<<TMath::RadToDeg()*angle<<" lies within range "
			<<mrdcluster::paddle_extentsy.at(pmts_hit.at(layerend)).first
			<<"<y<"<<mrdcluster::paddle_extentsy.at(pmts_hit.at(layerend)).second<<endl;
	} else {
		cout<<", xstart=";
		if(MaxMin==1){
			cout<<mrdcluster::paddle_extentsx.at(pmts_hit.at(layerstart)).second; 
			if(z2<z1){cout<<", going backwards and upwards to zend=";} else {cout<<", going downwards to zend=";}
		} else {
			cout<<mrdcluster::paddle_extentsx.at(pmts_hit.at(layerstart)).first;
			if(z2<z1){cout<<", going backwards and downwards to zend=";} else {cout<<", going upwards to zend=";}
		}
		cout<<z2<<" at angle "<<TMath::RadToDeg()*angle<<" lies within range "
			<<mrdcluster::paddle_extentsx.at(pmts_hit.at(layerend)).first
			<<"<x<"<<mrdcluster::paddle_extentsx.at(pmts_hit.at(layerend)).second<<endl;
	}
#endif
	Double_t ydiff = (zend-zstart)*TMath::Tan(angle);
	Double_t yproj;
	if(MaxMin==1){	//downgoing angle
	// layerstart defines a maximum angle - check projection from it's top to layerend at 'angle'
	// downwards (or at 'angle' upwards if projecting backward) strikes layerend
		if(z1==zstart){	// layerstart is on LHS - project forwards (-'ve y displacement)
			if(axis){	// y axis
				yproj = mrdcluster::paddle_extentsy.at(pmts_hit.at(layerstart)).second - ydiff;
			} else {	// x axis
				yproj = mrdcluster::paddle_extentsx.at(pmts_hit.at(layerstart)).second - ydiff;
			}
		} else {		// layerstart is on RHS - project backwards (+'ve y displacement)
			if(axis){
				yproj = mrdcluster::paddle_extentsy.at(pmts_hit.at(layerstart)).second + ydiff;
			} else {
				yproj = mrdcluster::paddle_extentsx.at(pmts_hit.at(layerstart)).second + ydiff;
			}
		}
	} else {	//upgoing angle
	// layerstart defines a min angle - check projection from it's bottom at angle upwards 
	// (if projecting forwards, or downwards if backwards) to layerend, strikes within 
	// layerend's y boundaries
		if(z1==zstart){	// layerstart is on LHS - project forwards (+'ve y displacement)
			if(axis){
				yproj = mrdcluster::paddle_extentsy.at(pmts_hit.at(layerstart)).first + ydiff;
			} else {
				yproj = mrdcluster::paddle_extentsx.at(pmts_hit.at(layerstart)).first + ydiff;
			}
		} else {		// layerstart is on RHS - project backwards (-'ve y displacement)
			if(axis){
				yproj = mrdcluster::paddle_extentsy.at(pmts_hit.at(layerstart)).first - ydiff;
			} else {
				yproj = mrdcluster::paddle_extentsx.at(pmts_hit.at(layerstart)).first - ydiff;
			}
		}
	}
	Double_t upperbound, lowerbound;
	if(axis){
		upperbound=mrdcluster::paddle_extentsy.at(pmts_hit.at(layerend)).first;
		lowerbound=mrdcluster::paddle_extentsy.at(pmts_hit.at(layerend)).second;
	} else {
		upperbound=mrdcluster::paddle_extentsx.at(pmts_hit.at(layerend)).first;
		lowerbound=mrdcluster::paddle_extentsx.at(pmts_hit.at(layerend)).second;
	}
	
	if( ((yproj+0.05) >= upperbound)&&((yproj-0.05) <= lowerbound)){
#ifdef MRDTrack_RECO_VERBOSE 
		cout<<"projected point is "<<yproj<<"; trajectory valid"<<endl;
#endif
		return true; 
	}	else { 
#ifdef MRDTrack_RECO_VERBOSE 
	cout<<"projected point is "<<yproj<<"; trajectory does not cross"<<endl;
#endif
	return false; 
	}
}
