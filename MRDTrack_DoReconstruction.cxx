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
//#define MRDTrack_RECO_VERBOSE
#endif

//#include "Math/Polynomial.h"  // can be used to find ROOTs of polynomials

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
	
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"getting hit layers and energy depositions"<<endl;
#endif
	int numdigits = digi_ids.size();
	for(Int_t i=0;i<numdigits;i++){
		// get the extent of the corresponding paddle
		int digiindex = digi_ids.at(i);
		Int_t tube_id = pmts_hit.at(i);
		Int_t strucklayer = mrdcluster::paddle_layers.at(tube_id);
		if(std::count(layers_hit.begin(), layers_hit.end(), strucklayer)==0){
			layers_hit.push_back(strucklayer);
		}
		eDepsInLayers.at(strucklayer)+= (digi_qs.at(i));	// TODO: convert q to energy?
	}
	
	DoTGraphErrorsFit();
	
	// calculate track fit start and endpoints from frontmost and backmost cell z values,
	// and fit formulae to determine x and y
	Double_t trackfitstartz = TMath::Min(mrdscintlayers.at(vtrackclusters.back().layer),
										mrdscintlayers.at(htrackclusters.back().layer));
	Double_t trackfitstarty = htrackorigin + htrackgradient*trackfitstartz;
	Double_t trackfitstartx = vtrackorigin + vtrackgradient*trackfitstartz;
	Double_t trackfitstopz = TMath::Max(mrdscintlayers.at(vtrackclusters.front().layer),
										mrdscintlayers.at(htrackclusters.front().layer));
	Double_t trackfitstopy = htrackorigin + htrackgradient*trackfitstopz;
	Double_t trackfitstopx = vtrackorigin + vtrackgradient*trackfitstopz;
	
	trackfitstart = TVector3(trackfitstartx,trackfitstarty,trackfitstartz);
	trackfitstop = TVector3(trackfitstopx,trackfitstopy,trackfitstopz);
	
	CheckIfStopping();
	
	// Projection to MRD front face: XXX do we need this? 
//	Double_t mrdentryz = MRD_start;
//	Double_t mrdentryx = vtrackorigin + vtrackgradient*mrdentryz;
//	Double_t mrdentryy = htrackorigin + htrackgradient*mrdentryz;
//	mrdentrypoint=TVector3(mrdentryx,mrdentryy,mrdentryz);
//	// and find allowed area at the MRD front face, for checking compatibility with tank forward projections
//	double mrdentryxmax = vtrackorigin + vtrackoriginerror + (vtrackgradient*mrdentryz);
//	double mrdentryymax = htrackorigin + htrackoriginerror + (htrackgradient*mrdentryz);
//	double mrdentryxmin = vtrackorigin - vtrackoriginerror - (vtrackgradient*mrdentryz);
//	double mrdentryymin = htrackorigin - htrackoriginerror - (htrackgradient*mrdentryz);
//	mrdentryxbounds=std::make_pair(mrdentryxmin, mrdentryxmax);
//	mrdentryybounds=std::make_pair(mrdentryymin, mrdentryymax);
	
	// Projection to tank, see if it's compatible with an interception.
	// to account for errors, use the shallowest angles allowed to give best chance of interception
	interceptstank = CheckTankIntercept(&projectedtankexitpoint,0,-1);
	// finally if compatible, note best fit interception point
	if(true||interceptstank){ // calculate it anyway, it may be useful during analysis/debugging
		// get the best fit interception point
		CheckTankIntercept(&projectedtankexitpoint,0,0);
	}
	
	CalculateEnergyLoss();
}

void cMRDTrack::CalculateEnergyLoss(){
	// determine energy of particle from track length and rate of energy loss.
	TVector3 differencevector  = trackfitstop - trackfitstart;
	TVector3 azaxisvector(0,0,1);
	trackangle = differencevector.Angle(azaxisvector);
	// TODO: MRDenergyvspenetration assumes a muon. Should we have functions for other particles?
	// TODO: For error on energy, need to fit the upper and lower bounds of dEdxVsAng. Then:
	// min possible E = value of fit to lower bounds of dEdxVsAng evaluated @ highest angle
	// max possible E = value of fit to upper bounds of dEdxVsAng evaluated @ lowest possible angle
	// error on angle for shortest tracks (40cm) is ~1rad, then falls to ~0.2 rads @ longest tracks (140cm)
	// for small angles (<0.5) & long tracks (>100cm ⇒ ang error 0.3 ) error  is ~+/-20%
	// for large angles (>1.5) & short tracks (<80cm ⇒ ang error 1.0 ) error is ~-50%/+300%
	// ~5MeV/cm  (5 to 6) @ large angle 
	// ~16MeV/cm (10 to >40) @ small angles
	// calculate the total track length in cm
	penetrationdepth=trackfitstopz-MRD_start;
	double muXdistanceinMRD=trackfitstopx-trackfitstartx;
	double muYdistanceinMRD=trackfitstopy-trackfitstarty;
	mutracklengthinMRD=
		TMath::Sqrt(TMath::Power(muXdistanceinMRD,2)+TMath::Power(muYdistanceinMRD,2)
		+TMath::Power(penetrationdepth,2));
	// calculate the energy loss
	double dEdx = MRDenergyvspenetration.Eval(trackangle);
	EnergyLoss = mutracklengthinMRD*dEdx;
	
	// calculate energy by recalculating with dEdx for maximum and minimum gradients
	double htrackgradientllim = htrackgradient+(htrackgradienterror*((htrackgradient>0) ? -1. : 1.));
	double vtrackgradientllim = vtrackgradient+(vtrackgradienterror*((vtrackgradient>0) ? -1. : 1.));
	double htrackgradienthlim = htrackgradient+(htrackgradienterror*((htrackgradient>0) ? 1. : -1.));
	double vtrackgradienthlim = vtrackgradient+(vtrackgradienterror*((vtrackgradient>0) ? 1. : -1.));
	
	double trackanglemax=TMath::Sqrt(TMath::Power(htrackgradienthlim,2)
									+TMath::Power(vtrackgradienthlim,2));
	double trackanglemin=TMath::Sqrt(TMath::Power(htrackgradientllim,2)
									+TMath::Power(vtrackgradientllim,2));
	double dEdxmax = MRDenergyvspenetration.Eval(trackanglemax); // TODO evaluate @ fit min when availabe
	double dEdxmin = MRDenergyvspenetration.Eval(trackanglemin); // TODO evaluate @ fit max when available
	double EnergyLossMax = mutracklengthinMRD*dEdxmax;
	double EnergyLossMin = mutracklengthinMRD*dEdxmin;
	EnergyLossError = max(EnergyLossMax-EnergyLoss,EnergyLoss-EnergyLossMin); // XXX asymmetric errors!
	
}

void cMRDTrack::DoTGraphErrorsFit(){
	
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"creating TGraphErrors arrays"<<endl;
#endif
	// Reconstruction based on a fit to a TGraphErrors of paddles, using extent of clusters as errors
	//std::vector<mrdcluster> htrackclusters;	// generated by CA algorithm
	//std::vector<mrdcluster> vtrackclusters;	// 
	std::vector<double> hclusterzpositions, hclusterxpositions, hclusterzerrors, hclusterxerrors;
	for(auto acluster : htrackclusters){
		hclusterzpositions.push_back(mrdscintlayers.at(acluster.layer));
		hclusterxpositions.push_back(acluster.GetCentre()/10.);
		hclusterzerrors.push_back(scintfullzlen);
		hclusterxerrors.push_back(TMath::Abs(acluster.GetXmax()-acluster.GetXmin())/10.);
	}
	std::vector<double> vclusterzpositions, vclusterxpositions, vclusterzerrors, vclusterxerrors;
	for(auto acluster : vtrackclusters){
		vclusterzpositions.push_back(mrdscintlayers.at(acluster.layer));
		vclusterxpositions.push_back(acluster.GetCentre()/10.);
		vclusterzerrors.push_back(scintfullzlen);
		vclusterxerrors.push_back(TMath::Abs(acluster.GetXmax()-acluster.GetXmin())/10.);
	}
	
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"creating TGraphErrors hclustergraph"<<endl;
#endif
	
	TGraphErrors hclustergraph = TGraphErrors(htrackclusters.size(), &hclusterzpositions[0], &hclusterxpositions[0], &hclusterzerrors[0], &hclusterxerrors[0]);
	
	// add additional points to the TGraphErrors that represent bonsai vertex, and/or possibly veto hit
	for(int pointi=0; pointi<extrahpoints.size(); pointi++){
		hclustergraph.SetPoint(htrackclusters.size()+i, extrahpoints.at(i), extrazpoints.at(i));
		hclustergraph.SetPointError(htrackclusters.size()+i, extrahpointerrors.at(i), extrazpointerrors.at(i));
	}
	
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"creating TGraphErrors vclustergraph"<<endl;
#endif
	TGraphErrors vclustergraph = TGraphErrors(vtrackclusters.size(), &vclusterzpositions[0], &vclusterxpositions[0], &vclusterzerrors[0], &vclusterxerrors[0]);
	
	// add additional points to the TGraphErrors that represent bonsai vertex, and/or possibly veto hit
	for(int pointi=0; pointi<extrahpoints.size(); pointi++){
		vclustergraph.SetPoint(vtrackclusters.size()+i, extravpoints.at(i), extrazpoints.at(i));
		vclustergraph.SetPointError(vtrackclusters.size()+i, extravpointerrors.at(i), extrazpointerrors.at(i));
	}
	
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"Drawing and getting fit parameters"<<endl;
#endif
	TCanvas c1;
	hclustergraph.Draw("AP");
	TF1 htrackfit("htrackfit","pol1",MRD_start,(MRD_start+MRD_depth));
	htrackfit.SetParameters(0,0);
#ifdef MRDTrack_RECO_VERBOSE
	TFitResultPtr htrackfitresult = hclustergraph.Fit(&htrackfit,"SR");
#else
	TFitResultPtr htrackfitresult = hclustergraph.Fit(&htrackfit,"SRQ");
#endif
	htrackorigin = htrackfitresult->Value(0);
	htrackoriginerror = htrackfitresult->ParError(0);
	htrackgradient = htrackfitresult->Value(1);
	htrackgradienterror = htrackfitresult->ParError(1);
	htrackfitchi2 = htrackfitresult->Chi2();
	
	//c1.SaveAs(TString::Format("htrackfit_%d.png",MRDtrackID));
	
	c1.Clear();
	vclustergraph.Draw("AP");
	TF1 vtrackfit("vtrackfit","pol1",MRD_start,(MRD_start+MRD_depth));
	vtrackfit.SetParameters(0,0);
#ifdef MRDTrack_RECO_VERBOSE
	TFitResultPtr vtrackfitresult = vclustergraph.Fit(&vtrackfit,"SR");
#else
	TFitResultPtr vtrackfitresult = vclustergraph.Fit(&vtrackfit,"SRQ");
#endif
	vtrackorigin = vtrackfitresult->Value(0);
	vtrackoriginerror = vtrackfitresult->ParError(0);
	vtrackgradient = vtrackfitresult->Value(1);
	vtrackgradienterror = vtrackfitresult->ParError(1);
	vtrackfitchi2 = vtrackfitresult->Chi2();
	//c1.SaveAs(TString::Format("vtrackfit_%d.png",MRDtrackID));
	
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"htrack fit had angle "<<((180/TMath::Pi())*TMath::ATan(htrackgradient))<<endl;
	cout<<"fit parameters were "<<htrackorigin<<", "<<htrackgradient
		<<" with errors "<<htrackoriginerror<<", "<<htrackgradienterror<<endl;
	cout<<"vtrack fit had angle "<<((180/TMath::Pi())*TMath::ATan(vtrackgradient))<<endl;
	cout<<"fit parameters were "<<vtrackorigin<<", "<<vtrackgradient
		<<" with errors "<<vtrackoriginerror<<", "<<vtrackgradienterror<<endl;
#endif
	
}

bool cMRDTrack::CheckTankIntercept(TVector3* entrypoint, TVector3* exitpoint=0, int tracktype=0){
	//entrypoint is entry BACK PROJECTING from mrd. Likewise exit point is exit point in upstream direction.
	double xgradient, ygradient, xoffset, yoffset;
	switch(tracktype){
	case 0:
		xgradient=vtrackgradient;
		ygradient=htrackgradient;
		xoffset=vtrackorigin;
		yoffset=htrackorigin;
		break;
	case 1:
		// steepest angles
		xgradient = vtrackgradient+(vtrackgradienterror*((vtrackgradient>0) ? -1. : 1.));
		ygradient = htrackgradient+(htrackgradienterror*((htrackgradient>0) ? -1. : 1.));
		xoffset=vtrackorigin+vtrackoriginerror;
		yoffset=htrackorigin+htrackoriginerror;
		break;
	case -1:
		// shallowest angles
		xgradient = vtrackgradient+(vtrackgradienterror*((vtrackgradient>0) ? 1. : -1.));
		ygradient = htrackgradient+(htrackgradienterror*((htrackgradient>0) ? 1. : -1.));
		xoffset=vtrackorigin-vtrackoriginerror;
		yoffset=htrackorigin-htrackoriginerror;
		break;
	}
	
	return CheckTankIntercept(ygradient, xgradient, yoffset, xoffset, entrypoint, exitpoint);
}

bool cMRDTrack::CheckTankIntercept(double htrackgradientin, double vtrackgradientin, double htrackoriginin, 
						double vtrackoriginin, TVector3* solution1, TVector3* solution2=0){
	// first as it's easiest, check if projection vertically actually enters tank height.
	double projectedtankexity = htrackoriginin + htrackgradientin*(tank_start+(2*tank_radius));
	double projectedtankexitz, projectedtankexitx;
	if(TMath::Abs(projectedtankexity-tank_yoffset)>tank_halfheight){
		return false;
	} else {
		// we know the track at least has height within the tank.
		// we now need to try to find an intercept in X-Z plane.
		// from simultaneous equations:
		// 1)  x^2 + z^2 = r^2
		// 2)  x = x0 + m*(z - z0) = m*z + (x0 - m*z0) = m*z + c (<< defines 'c')
		// we obtain:
		// z_intercept = [ -mc +/- Sqrt(m^2*c^2 + (1+m^2)*(r^2-c^2)) ] / (2*(1+m^2))
		// (with z axis shifted to centre of tank)
		double z0 = mrdscintlayers.at(vtrackclusters.back().layer) - tank_start - tank_radius;
		double linec = vtrackoriginin - (vtrackgradientin*z0);
		double coeffa = 1.+pow(vtrackgradientin,2);
		double coeffb = (2*vtrackgradientin*linec);
		double coeffc = -(pow(tank_radius,2) - pow(linec,2));
		// first check if there is an intercept
		double determinnt = pow(coeffb,2) - 4*coeffa*coeffc;
		if(determinnt>=0){
			projectedtankexitz = ( -coeffb + TMath::Sqrt(determinnt) ) / (2*coeffa);
			projectedtankexitz += tank_start + tank_radius; // remove tank-centering offset
			projectedtankexitx = vtrackoriginin + (vtrackgradientin*projectedtankexitz);
			// we need to recalculate this as the relevant z value will have changed
			projectedtankexity = htrackoriginin + htrackgradientin*projectedtankexitz;
			if(abs(projectedtankexity-tank_yoffset)>tank_halfheight) return false;
			// ^ bailing here neglects the possibility of an upgoing line exiting from the tank cap
			*solution1=TVector3(projectedtankexitx,projectedtankexity,projectedtankexitz);
			if(solution2!=0){
				projectedtankexitz = ( -coeffb - TMath::Sqrt(determinnt) ) / (2*coeffa);
				projectedtankexitz += tank_start + tank_radius; // remove tank-centering offset
				projectedtankexitx = vtrackoriginin + (vtrackgradientin*projectedtankexitz);
				projectedtankexity = htrackoriginin + htrackgradientin*projectedtankexitz;
				if(abs(projectedtankexity-tank_yoffset)>tank_halfheight){
					// the track exits a tank cap. 
					if((projectedtankexity-tank_yoffset)>tank_halfheight){
						// exits top cap
						projectedtankexity = tank_halfheight;
					} else {
						// exits bottom cap
						projectedtankexity = -tank_halfheight;
					}
					projectedtankexitz = (projectedtankexity-htrackoriginin)/htrackgradientin;
					projectedtankexitx = vtrackoriginin + (vtrackgradientin*projectedtankexitz);
				}
				*solution2=TVector3(projectedtankexitx,projectedtankexity,projectedtankexitz);
			}
			return true;
		} else {
			return false;
		}
	}
}

void cMRDTrack::AddTrackPoint(TVector3 pointposition, TVector3 pointerror){
	extravpoints.push_back(pointposition.X());
	extravpointerrors.push_back(pointerror.X());
	extrahpoints.push_back(pointposition.Y());
	extrahpointerrors.push_back(pointerror.Y());
	extrazpoints.push_back(pointposition.Z());
	extrazpointerrors.push_back(pointerror.Z());
}

void cMRDTrack::GetProjectionLimits(double zplane, double &xmax, double &xmin, double &ymax, double &ymin){
	// get the allowed region at a given z
	
	double htrackgradientllim = htrackgradient+(htrackgradienterror*((htrackgradient>0) ? -1. : 1.));
	double vtrackgradientllim = vtrackgradient+(vtrackgradienterror*((vtrackgradient>0) ? -1. : 1.));
	double htrackgradienthlim = htrackgradient+(htrackgradienterror*((htrackgradient>0) ? 1. : -1.));
	double vtrackgradienthlim = vtrackgradient+(vtrackgradienterror*((vtrackgradient>0) ? 1. : -1.));
	
	// first in y plane
	ymax= (htrackgradienthlim/htrackgradient)*htrackorigin + htrackoriginerror + htrackgradienthlim*zplane;
	ymin= (htrackgradientllim/htrackgradient)*htrackorigin - htrackoriginerror + htrackgradientllim*zplane;
	// in x plane
	xmax= (vtrackgradienthlim/vtrackgradient)*vtrackorigin + vtrackoriginerror + vtrackgradienthlim*zplane;
	xmin= (vtrackgradientllim/vtrackgradient)*vtrackorigin - vtrackoriginerror + vtrackgradientllim*zplane;
}

void cMRDTrack::GetProjectedPoint(double zplane){
	// get the track best fit projection 
	double yval = htrackorigin + htrackgradient*zplane;
	double xval = vtrackorigin + vtrackgradient*zplane;
	return TVector3(xval, yval, zplane);
}

void cMRDTrack::CheckIfStopping(){
	// define 'fiducial' MRD volume to call tracks that stop sufficiently far from MRD edges as 
	// 'stopping'
	double depthfidfrac = 0.9;
	double widthfidfrac = 0.9;
	double heightfidfrac = 0.9;
	if( TMath::Abs(trackfitstopz)>((MRD_start+MRD_depth)*depthfidfrac)){
		double projectedxexit = vtrackorigin + vtrackgradient*MRD_end;
		double projectedyexit = htrackorigin + htrackgradient*MRD_end;
		if( (TMath::Abs(projectedxexit)<MRD_width) && (TMath::Abs(projectedyexit)<MRD_height) ){
			ispenetrating=true;
		} else {
			sideexit=true;
		}
	} else if ( (TMath::Abs(trackfitstopx)>(MRD_width*widthfidfrac))   ||
				(TMath::Abs(trackfitstopy)>(MRD_height*heightfidfrac)) ){
		sideexit=true;
	} else {
		isstopped=true;
	}
}

double cMRDTrack::GetClosestApproach(TVector3 pointin, int tracktype){
	// tracktype =  0: best fit mrd track
	// tracktype =  1: upper limit of allowed track
	// tracktype = -1: lower limit of allowed track
	
	double xgradient, ygradient, xoffset, yoffset;
	switch(tracktype){
	case 0:
		xgradient=vtrackgradient;
		ygradient=htrackgradient;
		xoffset=vtrackorigin;
		yoffset=htrackorigin;
		break;
	case 1:
		xgradient = vtrackgradient+(vtrackgradienterror*((vtrackgradient>0) ? -1. : 1.));
		ygradient = htrackgradient+(htrackgradienterror*((htrackgradient>0) ? -1. : 1.));
		xoffset=vtrackorigin+vtrackoriginerror;
		yoffset=htrackorigin+htrackoriginerror;
		break;
	case -1:
		xgradient = vtrackgradient+(vtrackgradienterror*((vtrackgradient>0) ? 1. : -1.));
		ygradient = htrackgradient+(htrackgradienterror*((htrackgradient>0) ? 1. : -1.));
		xoffset=vtrackorigin-vtrackoriginerror;
		yoffset=htrackorigin-htrackoriginerror;
		break;
	}
	
	double xvala, yvala, zvala;
	zvala = pointa.Z()-500.;
	xvala = xoffset + xgradient*zvala;
	yvala = yoffset + ygradient*zvala;
	TVector3 pointa(xvala, yvala, zvala);
	double xvalb, yvalb, zvalb;
	zvalb = pointa.Z()+500.;
	xvalb = xoffset + xgradient*zvalb;
	yvalb = yoffset + ygradient*zvalb;
	TVector3 pointa(xvalb, yvalb, zvalb);
	TVector3 linea = pointin-pointa;
	TVector3 lineb = pointin-pointb;
	TVector3 acrossb = linea.Cross(lineb);
	TVector3 abdiff = pointb-pointa;
	
	double closestapp = acrossb.Mag() / abdiff.Mag();
	// x value first
}

TVector3 cMRDTrack::GetClosestPoint(TVector3 origin, int tracktype){
	// no idea how to calculate this for tracktypes other than zero - that's an open
	// issue in combinedanalysis on finding closest point for the projection of square-based 
	// pyramid representing possible mrd tracks to a bonsai vertex.
	if(tracktype!=0) return TVector3(0,0,0);
	
	// for the best fit, this should be possible
	double numerator =  vtrackgradient*(origin.X()-vtrackorigin) + 
						htrackgradient*(origin.Y()-htrackorigin) + origin.Z();
	double denominator = pow(vtrackgradient,2.) + pow(htrackgradient,2.) + 1;
	double closestz = numerator / denominator;
	double closestx = vtrackorigin + vtrackgradient*closestz;
	double closesty = htrackorigin + htrackgradient*closestz;
	return TVector3(closestx,closesty,closestz);
	
}
