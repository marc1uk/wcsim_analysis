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

//#include "Math/Polynomial.h"  // can be used to find ROOTs of polynomials

double xgradmin,ygradmin,xgradmax,ygradmax,xoffmin,yoffmin,xoffmax,yoffmax;

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
	
	//TVirtualFitter::SetDefaultFitter("Minuit");
	DoTGraphErrorsFit();
	
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"calculating track start and ends"<<endl;
#endif
	// calculate track fit start and endpoints from frontmost and backmost cell z values,
	// and fit formulae to determine x and y
	Double_t trackfitstartz = TMath::Min(MRDSpecs::mrdscintlayers.at(vtrackclusters.back().layer),
										MRDSpecs::mrdscintlayers.at(htrackclusters.back().layer));
	Double_t trackfitstarty = htrackorigin + htrackgradient*trackfitstartz;
	Double_t trackfitstartx = vtrackorigin + vtrackgradient*trackfitstartz;
	Double_t trackfitstopz = TMath::Max(MRDSpecs::mrdscintlayers.at(vtrackclusters.front().layer),
										MRDSpecs::mrdscintlayers.at(htrackclusters.front().layer));
	Double_t trackfitstopy = htrackorigin + htrackgradient*trackfitstopz;
	Double_t trackfitstopx = vtrackorigin + vtrackgradient*trackfitstopz;
	
	trackfitstart = TVector3(trackfitstartx,trackfitstarty,trackfitstartz);
	trackfitstop = TVector3(trackfitstopx,trackfitstopy,trackfitstopz);
	
	CheckIfStopping();
	if((isstopped&&ispenetrating)||(isstopped&&sideexit)||(ispenetrating&&sideexit)){
		cerr<<"isstopped="<<isstopped<<", ispenetrating="<<ispenetrating<<", sideexit="<<sideexit<<endl;
		assert(false);
	}
	
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
	// To account for errors, use the shallowest angles allowed to give best chance of interception
	//interceptstank = CheckTankIntercept(&projectedtankexitpoint,0,-1);
	// if compatible, note best fit interception point
	if(true||interceptstank){ // calculate it anyway, it may be useful during analysis/debugging
		// get the best fit interception point
		interceptstank = CheckTankIntercept(&projectedtankexitpoint,0,0); // FIXME - intercept bool should be based on shallowest...
	}
	//if(!interceptstank) assert(false&&"track did not come from tank?");
	
	CalculateEnergyLoss();
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"Reconstruction done"<<endl;
#endif
}

void cMRDTrack::CalculateEnergyLoss(){
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"Calculating energy loss"<<endl;
#endif
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
	penetrationdepth=trackfitstop.Z()-MRDSpecs::MRD_start;
	double muXdistanceinMRD=trackfitstop.X()-trackfitstart.X();
	double muYdistanceinMRD=trackfitstop.Y()-trackfitstart.Y();
	mutracklengthinMRD = 
		TMath::Sqrt(TMath::Power(muXdistanceinMRD,2)+TMath::Power(muYdistanceinMRD,2)
		+TMath::Power(penetrationdepth,2));
	// calculate the energy loss
	double dEdx = MRDenergyvspenetration.Eval(trackangle);
	EnergyLoss = mutracklengthinMRD*dEdx;
	
	// calculate energy by recalculating with dEdx for maximum and minimum gradients
	double htrackgradsteepest = htrackgradient+(htrackgradienterror*((htrackgradient>0) ? 1. : -1.));
	double vtrackgradsteepest = vtrackgradient+(vtrackgradienterror*((vtrackgradient>0) ? 1. : -1.));
	double htrackgradshallowest = htrackgradient+(htrackgradienterror*((htrackgradient>0) ? -1. : 1.));
	double vtrackgradshallowest = vtrackgradient+(vtrackgradienterror*((vtrackgradient>0) ? -1. : 1.));
	
	double trackanglemax=TMath::Sqrt(TMath::Power(htrackgradsteepest,2)
									+TMath::Power(vtrackgradsteepest,2));
	double trackanglemin=TMath::Sqrt(TMath::Power(htrackgradshallowest,2)
									+TMath::Power(vtrackgradshallowest,2));
	double dEdxmax = MRDenergyvspenetration.Eval(trackanglemax); // TODO evaluate @ fit min when availabe
	double dEdxmin = MRDenergyvspenetration.Eval(trackanglemin); // TODO evaluate @ fit max when available
	double EnergyLossMax = mutracklengthinMRD*dEdxmax;
	double EnergyLossMin = mutracklengthinMRD*dEdxmin;
	EnergyLossError = max(EnergyLossMax-EnergyLoss,EnergyLoss-EnergyLossMin); // XXX asymmetric errors!
	
}

void cMRDTrack::DoTGraphErrorsFit(){
/*	A function that constructs a TGraphErrors from the hit clusters and performs a straight line fit through them.
	These are used to obtain the h/vtrackorigin, h/vtrackgradient, along with their respective errors and covariances. */
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"creating TGraphErrors arrays"<<endl;
#endif
	// Reconstruction based on a fit to a TGraphErrors of paddles, using extent of clusters as errors
	//std::vector<mrdcluster> htrackclusters;	// generated by CA algorithm
	//std::vector<mrdcluster> vtrackclusters;	// 
	std::vector<double> hclusterzpositions, hclusterxpositions, hclusterzerrors, hclusterxerrors;
	for(auto acluster : htrackclusters){
		hclusterzpositions.push_back(MRDSpecs::mrdscintlayers.at(acluster.layer));
		hclusterxpositions.push_back(acluster.GetCentre()/10.);
		hclusterzerrors.push_back(MRDSpecs::scintfullzlen);
		hclusterxerrors.push_back(TMath::Abs(acluster.GetXmax()-acluster.GetXmin())/10.);
	}
	std::vector<double> vclusterzpositions, vclusterxpositions, vclusterzerrors, vclusterxerrors;
	for(auto acluster : vtrackclusters){
		vclusterzpositions.push_back(MRDSpecs::mrdscintlayers.at(acluster.layer));
		vclusterxpositions.push_back(acluster.GetCentre()/10.);
		vclusterzerrors.push_back(MRDSpecs::scintfullzlen);
		vclusterxerrors.push_back(TMath::Abs(acluster.GetXmax()-acluster.GetXmin())/10.);
	}
	
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"creating TGraphErrors hclustergraph"<<endl;
#endif
	
	TGraphErrors hclustergraph = TGraphErrors(htrackclusters.size(), &hclusterzpositions[0], &hclusterxpositions[0], &hclusterzerrors[0], &hclusterxerrors[0]);
	
	// XXX: note additional points must specify complete (x,y,z) points, as they (currently) share Z values. 
	// add additional points to the TGraphErrors that represent bonsai vertex, and/or possibly veto hit
	for(int i=0; i<extrahpoints.size(); i++){
		hclustergraph.SetPoint(htrackclusters.size()+i, extrahpoints.at(i), extrazpoints.at(i));
		hclustergraph.SetPointError(htrackclusters.size()+i, extrahpointerrors.at(i), extrazpointerrors.at(i));
	}
	
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"creating TGraphErrors vclustergraph"<<endl;
#endif
	TGraphErrors vclustergraph = TGraphErrors(vtrackclusters.size(), &vclusterzpositions[0], &vclusterxpositions[0], &vclusterzerrors[0], &vclusterxerrors[0]);
	
	// add additional points to the TGraphErrors that represent bonsai vertex, and/or possibly veto hit
	for(int i=0; i<extravpoints.size(); i++){
		vclustergraph.SetPoint(vtrackclusters.size()+i, extravpoints.at(i), extrazpoints.at(i));
		vclustergraph.SetPointError(vtrackclusters.size()+i, extravpointerrors.at(i), extrazpointerrors.at(i));
	}
	
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"Drawing and getting fit parameters"<<endl;
#endif
	// horizontal cluster fit
	// ----------------------
	TCanvas c1;
	hclustergraph.Draw("AP");
	TF1 htrackfit("htrackfit","pol1",MRDSpecs::MRD_start,(MRDSpecs::MRD_start+MRDSpecs::MRD_depth));
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
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"htrackfitchi2="<<htrackfitchi2<<endl;
	//auto returnedcovmatrix = htrackfitresult->GetCovarianceMatrix();
	//std::cout << type_name<decltype(returnedcovmatrix)>() << endl;
	cout<<"htrackfitcov has "<<htrackfitcov.GetNrows()<<" rows and "<<htrackfitcov.GetNcols()<<" columns"<<endl;
#endif
	htrackfitcov = htrackfitresult->GetCovarianceMatrix();
	
	//c1.SaveAs(TString::Format("htrackfit_%d.png",MRDtrackID));
	//TestCovarianceErrorCalc();  // it works!
	
	// vertical cluster fit
	// --------------------
	c1.Clear();
	vclustergraph.Draw("AP");
	TF1 vtrackfit("vtrackfit","pol1",MRDSpecs::MRD_start,(MRDSpecs::MRD_start+MRDSpecs::MRD_depth));
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
	vtrackfitcov = vtrackfitresult->GetCovarianceMatrix();
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"vtrackfitchi2="<<vtrackfitchi2<<endl;
	cout<<"vtrackfitcov has "<<vtrackfitcov.GetNrows()<<" rows and "<<vtrackfitcov.GetNcols()<<" columns"<<endl;
	//c1.SaveAs(TString::Format("vtrackfit_%d.png",MRDtrackID));
		cout<<"vtrack fit had angle "<<((180/TMath::Pi())*TMath::ATan(vtrackgradient))<<endl;
		cout<<"fit parameters were "<<vtrackorigin<<", "<<vtrackgradient
			<<" with errors "<<vtrackoriginerror<<", "<<vtrackgradienterror<<endl;
#endif
	
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"htrack fit had angle "<<((180/TMath::Pi())*TMath::ATan(htrackgradient))<<endl;
	cout<<"fit parameters were "<<htrackorigin<<", "<<htrackgradient
		<<" with errors "<<htrackoriginerror<<", "<<htrackgradienterror<<endl;
	cout<<"htrack fit covariance matrix: "<<endl;
	htrackfitcov.Print();
	cout<<"vtrack fit had angle "<<((180/TMath::Pi())*TMath::ATan(vtrackgradient))<<endl;
	cout<<"fit parameters were "<<vtrackorigin<<", "<<vtrackgradient
		<<" with errors "<<vtrackoriginerror<<", "<<vtrackgradienterror<<endl;
	cout<<"vtrack fit covariance matrix: "<<endl;
	vtrackfitcov.Print();
#endif
	
	/*
	// double check the steepest and shallowest angles by also drawing these on the TGraphErrors
	// steepest angles
	double xgraderrorsteep = (vtrackgradienterror*((vtrackgradient>0) ? 1. : -1.));
	double xgradientsteep = vtrackgradient+(vtrackgradienterror*((vtrackgradient>0) ? 1. : -1.));
	double ygraderrorsteep = (htrackgradienterror*((htrackgradient>0) ? 1. : -1.));
	double ygradientsteep = htrackgradient+(htrackgradienterror*((htrackgradient>0) ? 1. : -1.));
	double xoffseterrorsteep = CalculateCovariantError(xgraderrorsteep, 1, vtrackfitcov);
	double xoffsetsteep = vtrackorigin + xoffseterrorsteep;
	double yoffseterrorsteep = CalculateCovariantError(ygraderrorsteep, 1, htrackfitcov);
	double yoffsetsteep = htrackorigin + yoffseterrorsteep;
	// shallowest angles
	double xgraderrorshallow = (vtrackgradienterror*((vtrackgradient>0) ? -1. : 1.));
	double xgradientshallow = vtrackgradient+(vtrackgradienterror*((vtrackgradient>0) ? -1. : 1.));
	double ygraderrorshallow = (htrackgradienterror*((htrackgradient>0) ? -1. : 1.));
	double ygradientshallow = htrackgradient+(htrackgradienterror*((htrackgradient>0) ? -1. : 1.));
	double xoffseterrorshallow = CalculateCovariantError(xgraderrorshallow, 1, vtrackfitcov);
	double xoffsetshallow = vtrackorigin + xoffseterrorshallow;
	double yoffseterrorshallow = CalculateCovariantError(ygraderrorshallow, 1, htrackfitcov);
	double yoffsetshallow = htrackorigin + yoffseterrorshallow;
	
	cout<<"best fit x line has gradient "<<vtrackgradient<<", error "<<vtrackgradienterror<<", steepest gradient "<<xgradientsteep<<", shallowest gradient "<<xgradientshallow<<endl;
	cout<<"best fit offset is "<<vtrackorigin<<", error "<<vtrackoriginerror<<", steepest offset "<<xoffsetsteep<<", shallowest offset "<<xoffsetshallow<<endl;
	c2.cd();
	TF1* xsteepestline = new TF1("xsteep","[0]*x+[1]",MRD_start,(MRD_start+MRD_depth));
	xsteepestline->SetParameters(xgradientsteep,xoffsetsteep);
	xsteepestline->SetLineColor(kBlue);
	xsteepestline->Draw("same");
	TF1* xshallowestline = new TF1("xshallow","[0]*x+[1]",MRD_start,(MRD_start+MRD_depth));
	xshallowestline->SetParameters(xgradientshallow,xoffsetshallow);
	xshallowestline->SetLineColor(kRed);
	xshallowestline->Draw("same");
	c2.Modified();
	
	c1.cd();
	TF1* ysteepestline = new TF1("ysteep","[0]*x+[1]",MRD_start,(MRD_start+MRD_depth));
	ysteepestline->SetParameters(ygradientsteep,yoffsetsteep);
	ysteepestline->SetLineColor(kBlue);
	ysteepestline->Draw("same");
	TF1* yshallowestline = new TF1("yshallow","[0]*x+[1]",MRD_start,(MRD_start+MRD_depth));
	yshallowestline->SetParameters(ygradientshallow,yoffsetshallow);
	yshallowestline->SetLineColor(kRed);
	yshallowestline->Draw("same");
	c1.Modified();
	
	gPad->WaitPrimitive();
	*/
}

std::pair<double, double> cMRDTrack::CalculateCovariantError(double errorx, int indexgiven, TMatrixDSym covariancematrix, bool flag){
	/* calculating error in x based on error in y given covariance.
	for a distribution (position) described by a function of two parameters (angle, offset),
	the parameters best describing the distribution may be found my minimizing the log likelihood
	of the parameter values based on estimators which derive estimates of the parameter values from 
	samples of the distribution. 
	The variation of the log-likelihood with parameter values, about the maximum, is described by:
	ln(L) = ln(Lmax) - \frac{1}{2(1-rho^2)}*[ (\frac{(x-xmax)^2}{var_x}) + (\frac{(y-ymax)^2}{var_y}) 
					 - 2*(\frac{x-xmax}{sqrt(var_x)})*(\frac{y-ymax}{sqrt(var_y)}) ]
	where var_x is the variance in x and rho = \frac{cov_xy}{var_x*var_y} is the correlation coefficient of x-y. 
	For a single parameter, a difference of 1 std_dev from the maximized log-likelihood parameter value results
	in a log-likelihood reduced by 1/2 (absolute!). For two parameters, a reduction of log-likelihood by 1/2
	describes an ellipse in the x-y plane, whose equation is found from (ln(L)-ln(Lmax))=1/2, with the formula above.
	The ellipse is centred on (xmax,ymax), has maximum and minimum in x and y axes of (x-var_x) and (y-var_y),
	and is inclined at angle alpha given by:
	tan(2*alpha) = \frac{2*rho*sqrt(var_x)*sqrt(var_y)}{var_x-var_y}
	
	To find the value of the parameter y, when the log-likelihood is reduce by 1/2 and the other parameter x 
	is at it's minimum (x_max - var_x), we can solve this equation for the difference (y-y_max)/sqrt(var_y).
	Defining A=\frac{x-xmax}{sqrt(var_x)} and B=\frac{y-ymax}{sqrt(var_y)} we can rearrange for B to obtain:
	B = rho*A +/- Sqrt[ (rho^2-1)*(A^2-1) ]
	*/
	//cout<<"calculating offset error for gradient "<<errorx<<endl;
	int indextocalculate=0; // which row/column of the covariance matrix is the passed error.
	if(indexgiven==0) indextocalculate=1;
	double myA = errorx/sqrt(covariancematrix(indexgiven,indexgiven)); // note sigma, one SD, is the SQRT of variance
	double rho = covariancematrix(0,1)/(sqrt(covariancematrix(0,0)*covariancematrix(1,1)));
	//cout<<"myA="<<myA<<", rho="<<rho<<endl;
	double firstterm = rho*myA;
	double secondterm = (pow(rho,2.)-1)*(pow(myA,2.)-1);
	//cout<<"first="<<firstterm<<", second="<<secondterm<<", third="<<thirdterm<<endl;
	double solution1 = firstterm + sqrt(secondterm);
	double solution2 = firstterm - sqrt(secondterm);
	//cout<<"sol1="<<solution1<<", 2="<<solution2<<endl;
	// solutions are values of B, which is error, in units of variance. 
	// We multiply by variance -  BUT NEED TO ***add best fit val to get actual parameter estimate***.
	double value1 = solution1*sqrt(covariancematrix(indextocalculate,indextocalculate));
	double value2 = solution2*sqrt(covariancematrix(indextocalculate,indextocalculate));
	//cout<<"pair=("<<value1<<", "<<value2<<")"<<endl;
	return std::make_pair(value1,value2);
}

double cMRDTrack::CalculateCovariantError(double errorx, int indexgiven, TMatrixDSym covariancematrix){
	std::pair<double,double> bothresults = CalculateCovariantError(errorx, indexgiven, covariancematrix, true);
	return bothresults.first;
}

void cMRDTrack::TestCovarianceErrorCalc(){
	// to check the calculation of error in offset given error in gradient, we'll draw the ellipse by varying
	// the gradient error, then see if it matches the ellipse from minuit. parater 0=offset, 1=gradient
	//cout<<"covariancematrix(1,1) = "<<vtrackfitcov(1,1)<<", vtrackgradienterror^2="<<pow(vtrackgradienterror,2)<<endl;
	std::vector<double> xpoints;
	std::vector<double> ypoints(200,0);
	double xpointstep = vtrackgradienterror/50;
	for(int pointi=0; pointi<100; pointi++){  // 100 points around the ellipse
		double nextgraderror = - vtrackgradienterror + pointi*xpointstep; // spans - vtrackgraderr to + vtrackgraderr
		xpoints.push_back(vtrackgradient+nextgraderror);
		double offseterror = CalculateCovariantError(nextgraderror, 1, vtrackfitcov);
		ypoints.at(pointi)= offseterror + vtrackorigin;
		ypoints.at(ypoints.size()-pointi-1)= offseterror + vtrackorigin;
	}
	xpoints.insert(xpoints.end(),xpoints.rbegin(),xpoints.rend());
	cout<<"making testgraph from points:"<<endl;
	for(int i=0; i<xpoints.size(); i++){ cout<<"("<<xpoints.at(i)<<", "<<ypoints.at(i)<<")"<<endl; }
	TGraph* testgraph = new TGraph(xpoints.size(), &ypoints[0], &xpoints[0]);
	TCanvas* c1 = new TCanvas();
	cout<<"drawing manual ellipse"<<endl;
	testgraph->SetMarkerStyle(kOpenCircle);
	testgraph->Draw("AP");
	testgraph->GetYaxis()->SetTitle("parameter 1 (gradient)");
	testgraph->GetXaxis()->SetTitle("parameter 0 (intercept)");
	testgraph->SetTitle("1-sigma uncertainties on fit parameters");
	
	TMatrixD covmat(2,2);
	gMinuit->mnemat(covmat.GetMatrixArray(),2);
	cout<<"printing minuit covmat"<<endl;
	covmat.Print(); 
	cout<<"printed. Printing correlation matrix:"<<endl;
	gMinuit->mnmatu(1);
	//TCanvas* c2 = new TCanvas();
	cout<<"drawing minuit ellipse"<<endl;
	cout<<"setting errordef"<<endl;
	gMinuit->SetErrorDef(1);
	cout<<"getting graph"<<endl;
	TGraph* mingraph = (TGraph*)gMinuit->Contour(50,0,1);
	cout<<"graph at "<<mingraph<<", drawing"<<endl;
	mingraph->SetLineColor(kRed);
	mingraph->SetMarkerStyle(2);
	mingraph->Draw("l same");
	mingraph->SetTitle("Minuit Plot");
	gPad->WaitPrimitive();
	
}

bool cMRDTrack::CheckTankIntercept(TVector3* exitpoint, TVector3* entrypoint=0, int tracktype=0){
#ifdef MRDTrack_RECO_VERBOSE
	if(tracktype==0) cout<<"Checking best fit track tank interception"<<endl;
	if(tracktype==-1) cout<<"Checking minimum track tank interception"<<endl;
	if(tracktype==1) cout<<"Checking maximum track tank interception"<<endl;
#endif
	//entrypoint is entry BACK PROJECTING from mrd. Likewise exit point is exit point in upstream direction.
	double xgradient, ygradient, xoffset, yoffset;
	GetAngleAndOffset(tracktype, xgradient, ygradient, xoffset, yoffset);
	
	return CheckTankIntercept(ygradient, xgradient, yoffset, xoffset, exitpoint, entrypoint);
}

bool cMRDTrack::CheckTankIntercept(double htrackgradientin, double vtrackgradientin, double htrackoriginin, 
						double vtrackoriginin, TVector3* solution1, TVector3* solution2=0){
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"Checking for and calculating track tank interception points"<<endl;
#endif
	// first as it's easiest, check if projection vertically actually enters tank height.
	double projectedtankexity = htrackoriginin + 
		htrackgradientin*(MRDSpecs::tank_start+(2*MRDSpecs::tank_radius));
	double projectedtankexitz, projectedtankexitx;
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"initial back projection of track to tank back edge places projected height at "
		<<projectedtankexity-MRDSpecs::tank_yoffset
		<<" relative to tank half-height "<<MRDSpecs::tank_halfheight<<endl;
#endif
	if(TMath::Abs(projectedtankexity-MRDSpecs::tank_yoffset)>MRDSpecs::tank_halfheight){
#ifdef MRDTrack_RECO_VERBOSE
		cout<<"tank misses tank height"<<endl;
#endif
		return false;
	} else {
		// we know the track at least has height within the tank.
		// we now need to try to find an intercept in X-Z plane.
		// from simultaneous equations:
		// 1)  x^2 + z^2 = r^2
		// 2)  x = x0 + m*(z - z0) = m*z + (x0 - m*z0) = m*z + c (<< defines 'c') <<<< wrong?
		// 3)  x = mz + c = m(z'+tank_start+tank_radius) + c = mz' + c'; c' = c + m(tank_start+tank_radius)
		// we obtain:
		// z_intercept = [ -mc +/- Sqrt(m^2*c^2 - (1+m^2)*(c^2-r^2)) ] / (1+m^2)
		// (with z axis shifted to centre of tank)
		double linec = vtrackoriginin + vtrackgradientin*(MRDSpecs::tank_start+MRDSpecs::tank_radius);
		double coeffa = 1.+pow(vtrackgradientin,2);
		double coeffb = 2*vtrackgradientin*linec;
		double coeffc = pow(linec,2) - pow(MRDSpecs::tank_radius,2);
		// first check if there is an intercept
		double determinnt = pow(coeffb,2) - 4*coeffa*coeffc;
#ifdef MRDTrack_RECO_VERBOSE
		cout<<"x-z plane interception equation has determinant "<<determinnt<<endl;
#endif
		if(determinnt>=0){
			projectedtankexitz = ( -coeffb + TMath::Sqrt(determinnt) ) / (2*coeffa);
#ifdef MRDTrack_RECO_VERBOSE
			cout<<"x-z plane interception occurs at z= "<<projectedtankexitz
				<<" in coordinates with origin centred on the tank"<<endl;
#endif
			projectedtankexitz += MRDSpecs::tank_start + MRDSpecs::tank_radius; // remove tank-centering offset
			projectedtankexitx = vtrackoriginin + (vtrackgradientin*projectedtankexitz);
			// we need to recalculate this as the relevant z value will have changed
			projectedtankexity = htrackoriginin + htrackgradientin*projectedtankexitz;
#ifdef MRDTrack_RECO_VERBOSE
			cout<<"interception point in global coords is ("<<projectedtankexitx<<", "
				<<projectedtankexity<<", "<<projectedtankexitz<<"). Rechecking y in range."<<endl;
#endif
			if(abs(projectedtankexity-MRDSpecs::tank_yoffset)>MRDSpecs::tank_halfheight)
				{ cout<<"y range fail"<<endl; return false; }
			*solution1=TVector3(projectedtankexitx,projectedtankexity,projectedtankexitz);
#ifdef MRDTrack_RECO_VERBOSE
			cout<<"y value at interception to tank cylinder is within tank height, doing sanity checks"<<endl;
#endif
			// sanity check:
			if( (abs(projectedtankexitx)>MRDSpecs::tank_radius)                            ||
				(abs(projectedtankexity-MRDSpecs::tank_yoffset)>MRDSpecs::tank_halfheight) ||
				(projectedtankexitz<MRDSpecs::tank_start)                                  ||
				(projectedtankexitz>(MRDSpecs::tank_start+(2*MRDSpecs::tank_radius)))      ||
				(sqrt(pow(projectedtankexitz-MRDSpecs::tank_start-MRDSpecs::tank_radius,2)
					 +pow(projectedtankexitx,2))>MRDSpecs::tank_radius) ){
				cerr<<"projected exit point is OUTSIDE the tank!"<<endl;
				cerr<<"x > tank radius           : "<<(abs(projectedtankexitx)>MRDSpecs::tank_radius)<<endl;
				cerr<<"y > tank height           : "<<(abs(projectedtankexity
				-MRDSpecs::tank_yoffset)>MRDSpecs::tank_halfheight)<<endl;
				cerr<<"z < tank start            : "<<(projectedtankexitz<MRDSpecs::tank_start)<<endl;
				cerr<<"z > tank end              : "<<(projectedtankexitz>
					(MRDSpecs::tank_start+(2*MRDSpecs::tank_radius)))<<endl;
				cerr<<"(x^2 + z^2) > tank radius : "<<
				(sqrt(pow(projectedtankexitz-MRDSpecs::tank_start-MRDSpecs::tank_radius,2)+
					  pow(projectedtankexitx,2))>MRDSpecs::tank_radius)<<endl;
				assert(false);
			}
			// second sanity check:
			if( (abs((projectedtankexity-MRDSpecs::tank_yoffset)-MRDSpecs::tank_halfheight)>10.) &&
				(abs(sqrt(pow(projectedtankexitz-MRDSpecs::tank_start-MRDSpecs::tank_radius,2)
						 +pow(projectedtankexitx,2))-MRDSpecs::tank_radius)>10.) ){
				cerr<<"projected exit point doesn't lie on the tank surface!"<<endl;
				assert(false);
			}
#ifdef MRDTrack_RECO_VERBOSE
			cout<<"sanity checks passed";
#endif
			
			// if requested, also project backwards to the tank entry point
			if(solution2!=0){
#ifdef MRDTrack_RECO_VERBOSE
				cout<<", calculating corresponding tank entry point"<<endl;
#endif
				projectedtankexitz = ( -coeffb - TMath::Sqrt(determinnt) ) / (2*coeffa);
				projectedtankexitz += MRDSpecs::tank_start + MRDSpecs::tank_radius; // remove tank-center offset
				projectedtankexitx = vtrackoriginin + (vtrackgradientin*projectedtankexitz);
				projectedtankexity = htrackoriginin + htrackgradientin*projectedtankexitz;
#ifdef MRDTrack_RECO_VERBOSE
			cout<<"second solution (tank exit) in global coords is ("<<projectedtankexitx<<", "
				<<projectedtankexity<<", "<<projectedtankexitz<<"). Rechecking y in range."<<endl;
#endif
				// since we project tracks back from the MRD, the exit point must be the barrel
				// (cap exits would not produce an MRD track to back-project), but for the tank entry point
				// we need to check whether the projected interception point is outside a cap. If so, recalculate.
				if(abs(projectedtankexity-MRDSpecs::tank_yoffset)>MRDSpecs::tank_halfheight){
					// the track enters via a tank cap. 
					if((projectedtankexity-MRDSpecs::tank_yoffset)>MRDSpecs::tank_halfheight){
						// enters top cap
#ifdef MRDTrack_RECO_VERBOSE
						cout<<"tank enters top cap, recalculating"<<endl;
#endif
						projectedtankexity = MRDSpecs::tank_halfheight;
					} else {
						// enters bottom cap
#ifdef MRDTrack_RECO_VERBOSE
						cout<<"tank enters bottom cap, recalculating"<<endl;
#endif
						projectedtankexity = -MRDSpecs::tank_halfheight;
					}
					projectedtankexitz = (projectedtankexity-htrackoriginin)/htrackgradientin;
					projectedtankexitx = vtrackoriginin + (vtrackgradientin*projectedtankexitz);
				}
#ifdef MRDTrack_RECO_VERBOSE
				cout<<"second solution is ("<<projectedtankexitx<<", "
					<<projectedtankexity<<", "<<projectedtankexitz<<") "<<endl;
#endif
				*solution2=TVector3(projectedtankexitx,projectedtankexity,projectedtankexitz);
			}
#ifdef MRDTrack_RECO_VERBOSE
			cout<<", tank entry point not requested, returning true"<<endl;
#endif
			return true;
		} else {
#ifdef MRDTrack_RECO_VERBOSE
			cout<<"track projection misses tank barrel."<<endl;
#endif
			return false;
		}
	}
}

void cMRDTrack::AddTrackPoint(TVector3 pointposition, TVector3 pointerror){
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"Adding external track constraint point"<<endl;
#endif
	extravpoints.push_back(pointposition.X());
	extravpointerrors.push_back(pointerror.X());
	extrahpoints.push_back(pointposition.Y());
	extrahpointerrors.push_back(pointerror.Y());
	extrazpoints.push_back(pointposition.Z());
	extrazpointerrors.push_back(pointerror.Z());
}

void cMRDTrack::GetProjectionLimits(double zplane, double &xmax, double &xmin, double &ymax, double &ymin){
	// get the allowed region at a given z
	
	double htrackgradientupgoing = htrackgradient-htrackgradienterror;
	double vtrackgradientupgoing = vtrackgradient-vtrackgradienterror;
	double htrackgradientdowngoing = htrackgradient+htrackgradienterror;
	double vtrackgradientdowngoing = vtrackgradient+vtrackgradienterror;
	
	double htrackupgoingoff = CalculateCovariantError(htrackgradient-htrackgradientupgoing, 1, htrackfitcov)+htrackorigin;
	double htrackdowngoingoff = CalculateCovariantError(htrackgradient-htrackgradientdowngoing, 1, htrackfitcov)+htrackorigin;
	double vtrackupgoingoff = CalculateCovariantError(vtrackgradient-vtrackgradientupgoing, 1, vtrackfitcov)+vtrackorigin;
	double vtrackdowngoingoff = CalculateCovariantError(vtrackgradient-vtrackgradientdowngoing, 1, vtrackfitcov)+vtrackorigin;
	
	// assumes projection is backward, so that highest y is given by most downgoing gradient
	// first in y plane
	ymin= htrackgradientupgoing*zplane + htrackupgoingoff;
	ymax= htrackgradientdowngoing*zplane + htrackdowngoingoff;
	// in x plane
	xmin= vtrackgradientupgoing*zplane + vtrackupgoingoff;
	xmax= vtrackgradientdowngoing*zplane + vtrackdowngoingoff;
}

TVector3 cMRDTrack::GetProjectedPoint(double zplane){
	// get the track best fit projection 
	double yval = htrackorigin + htrackgradient*zplane;
	double xval = vtrackorigin + vtrackgradient*zplane;
	return TVector3(xval, yval, zplane);
}

void cMRDTrack::CheckIfStopping(){
/*	Test whether the track fit stop vertex lies within the fiducial volume of the MRD. */
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"checking if stopping"<<endl;
#endif
	// define 'fiducial' MRD volume to call tracks that stop sufficiently far from MRD edges as 'stopping'
	double depthfidfrac = 0.9;
	double widthfidfrac = 0.9;
	double heightfidfrac = 0.9;
	if( TMath::Abs(trackfitstop.Z())>(MRDSpecs::MRD_start+(MRDSpecs::MRD_depth*depthfidfrac)) ){
		// if the last point is sufficiently close to the MRD end, we consider it as penetrated
		// and take the 'zstop' as MRD_end. But, we can check from the angle if it would have exited through the sides first.
		double projectedxexit = vtrackorigin + vtrackgradient*MRDSpecs::MRD_end;
		double projectedyexit = htrackorigin + htrackgradient*MRDSpecs::MRD_end;
		if( (TMath::Abs(projectedxexit)<MRDSpecs::MRD_width) && 
			(TMath::Abs(projectedyexit)<MRDSpecs::MRD_height) ){
			ispenetrating=true;
		} else {
			sideexit=true;
		}
	} else if ( (TMath::Abs(trackfitstop.X())>(MRDSpecs::MRD_width*widthfidfrac))   ||
				(TMath::Abs(trackfitstop.Y())>(MRDSpecs::MRD_height*heightfidfrac)) ){
		sideexit=true;
	} else {
		isstopped=true;
	}
}

double cMRDTrack::GetClosestApproach(TVector3 pointin, int tracktype){
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"Getting closest approach distance"<<endl;
#endif
	// tracktype =  0: best fit mrd track
	// tracktype =  1: upper limit of allowed track, steepest angle
	// tracktype = -1: lower limit of allowed track, shallowest angle
	
	double xgradient, ygradient, xoffset, yoffset;
	GetAngleAndOffset(tracktype, xgradient, ygradient, xoffset, yoffset);
	
	double xvala, yvala, zvala;
	zvala = pointin.Z()-500.;
	xvala = xoffset + xgradient*zvala;
	yvala = yoffset + ygradient*zvala;
	TVector3 pointa(xvala, yvala, zvala);
	double xvalb, yvalb, zvalb;
	zvalb = pointin.Z()+500.;
	xvalb = xoffset + xgradient*zvalb;
	yvalb = yoffset + ygradient*zvalb;
	TVector3 pointb(xvalb, yvalb, zvalb);
	TVector3 linea = pointin-pointa;
	TVector3 lineb = pointin-pointb;
	TVector3 acrossb = linea.Cross(lineb);
	TVector3 abdiff = pointb-pointa;
	
	double closestapp = acrossb.Mag() / abdiff.Mag();
	return closestapp;
}

TVector3 cMRDTrack::GetClosestPoint(TVector3 origin, int tracktype){
#ifdef MRDTrack_RECO_VERBOSE
	cout<<"Getting closest point of approach"<<endl;
#endif
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

void cMRDTrack::GetAngleAndOffset(int tracktype, double &xgradient, double &ygradient, double &xoffset, double &yoffset){
	/* Get the gradient and offset of the steepest or shallowest track */
	// tracktype =  0: best fit track
	// tracktype =  1: steepest angle
	// tracktype = -1: shallowest angle
	
	switch(tracktype){
		case 0: {
			xgradient=vtrackgradient;
			ygradient=htrackgradient;
			xoffset=vtrackorigin;
			yoffset=htrackorigin;
			break;
		}
		case 1: { // you need to enclose the case block in { ... } to be able to declare variables within a block
			// steepest
			xgradient = vtrackgradient+(vtrackgradienterror*((vtrackgradient>0) ? 1. : -1.));
			ygradient = htrackgradient+(htrackgradienterror*((htrackgradient>0) ? 1. : -1.));
			// since angle and offset are correlated, we need to find the corresponding error in the offset
			// when specifying the angle:
			double xoffseterror = CalculateCovariantError(vtrackgradienterror, 1, vtrackfitcov);
			xoffset=vtrackorigin + xoffseterror;
			double yoffseterror = CalculateCovariantError(htrackgradienterror, 1, htrackfitcov);
			yoffset=htrackorigin + yoffseterror;
			// XXX first or second - for the extremities it shouldn't matter as at these points the two converge.
			// however this would need to be considered for other points - there are two possible values of error in A
			// that give the same total error, for a given error in B.
			break;
		}
		case -1: { // shallowest
			xgradient = vtrackgradient+(vtrackgradienterror*((vtrackgradient>0) ? -1. : 1.));
			ygradient = htrackgradient+(htrackgradienterror*((htrackgradient>0) ? -1. : 1.));
			double xoffseterror = CalculateCovariantError(vtrackgradienterror, 1, vtrackfitcov);
			xoffset=vtrackorigin + xoffseterror;
			double yoffseterror = CalculateCovariantError(htrackgradienterror, 1, htrackfitcov);
			yoffset=htrackorigin + yoffseterror;
			break;
		}
	}
	
}
