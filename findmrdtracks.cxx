// #######################################################################

// CREATE+OPEN OUTPUT FILE
// =======================
void WCSimAnalysis::OpenMRDtrackOutfile(){
	outfile = new TFile("outfile","outfile");
	outfile->cd();
	recotree = new TTree("recotree","Tree for reconstruction data");
	nummrddigitsthisevent=0;
	nummrdtracksthisevent=0;

	nummrddigitsthiseventb = recotree->Branch("nummrddigitsthisevent",&nummrddigitsthisevent);
	nummrdtracksthiseventb = recotree->Branch("nummrdtracksthisevent",&nummrdtracksthisevent);
	tubeidsinthistrackb = recotree->Branch("tubeidsinthistrack",&tubeidsinthistrack);
	digitqsinthistrackb = recotree->Branch("digitqsinthistrack",&digitqsinthistrack);
	digittimesinthistrackb = recotree->Branch("digittimesinthistrack",&digittimesinthistrack);
	particleidsinthistrackb = recotree->Branch("particleidsinthistrack",&particleidsinthistrack);
}

// #######################################################################

// SPLIT HITS BY TIME
// ==================
void WCSimAnalysis::SplitMrdTracks(std::vector<double> mrddigittimesthisevent){
	
	// MEASURE EVENT DURATION TO DETERMINE IF THERE IS MORE THAN ONE TRACK
	// ===================================================================
	std::vector<Int_t> tubeidsinatrack;
	std::vector<Double_t> digitqsinatrack;
	std::vector<Double_t> digittimesinatrack;
	std::vector<Int_t> particleidsinatrack;
	
	// first check: are all hits within a 30ns window (maxtrackduration) If so, just one track. 
	Double_t eventendtime = *std::max_element(mrddigittimesthisevent.begin(),mrddigittimesthisevent.end());
	Double_t eventstarttime = *std::min_element(mrddigittimesthisevent.begin(),mrddigittimesthisevent.end());
	Double_t eventduration = (eventendtime - eventstarttime);
	//cout<<"event start: "<<eventstarttime<<", end : "<<eventendtime<<", duration : "<<eventduration<<endl;
	
	if(eventduration<maxtrackduration){
	// JUST ONE TRACK
	// ==============
		//cout<<"all hits this event within one track."<<endl;
		
		for(Int_t thisdigit=0;thisdigit<mrddigittimesthisevent.size();thisdigit++){
			WCSimRootCherenkovDigiHit* thedigihit = (WCSimRootCherenkovDigiHit*)atrigm->GetCherenkovDigiHits()->At(thisdigit);
			Int_t thisdigitstubeid = thedigihit->GetTubeId();
			tubeidsinatrack.push_back(thisdigitstubeid);
			Int_t thisdigitsq = thedigihit->GetQ();
			digitqsinatrack.push_back(thisdigitsq);
			double thisdigitstime = thedigihit->GetT();
			digittimesinatrack.push_back(thisdigitstime);
			// add all the unique parent ID's for digits contributing to this track (truth level info)
			std::vector<int> truephotonindices = thedigihit->GetPhotonIds();
			for(int truephoton=0; truephoton<truephotonindices.size(); truephoton++){
				int thephotonsid = truephotonindices.at(truephoton);
				WCSimRootCherenkovHitTime *thehittimeobject = (WCSimRootCherenkovHitTime*)atrigm->GetCherenkovHitTimes()->At(thephotonsid);
				int thephotonsparenttrack = thehittimeobject->GetParentID();
				int checkcount = std::count(particleidsinatrack.begin(), particleidsinatrack.end(), thephotonsparenttrack);
				if(checkcount==0){ particleidsinatrack.push_back(thephotonsparenttrack); }
			}
		}
		
		nummrddigitsthisevent = mrddigittimesthisevent.size();
		nummrdtracksthisevent=1;
		tubeidsinthistrack.push_back(tubeidsinatrack);
		digitqsinthistrack.push_back(digitqsinatrack);
		digittimesinthistrack.push_back(digittimesinatrack);
		particleidsinthistrack.push_back(particleidsinatrack);
		
		nummrddigitsthiseventb->Fill();
		nummrdtracksthiseventb->Fill();
		tubeidsinthistrackb->Fill();
		digitqsinthistrackb->Fill();
		digittimesinthistrackb->Fill();
		particleidsinthistrackb->Fill();
		
		// end of event with all hits in one mrd 'track'
		
	} else {
	// MORE THAN ONE TRACK
	// ===================
		// this event has multiple tracks. Need to split hits into which track they belong to.
		// scan over the times and look for gaps where no digits lie, using these to delimit 'tracks'
		std::vector<Float_t> trackhittimesv;	// a vector of the starting times of a given 'track'
		std::vector<double> sorteddigittimes(mrddigittimesthisevent);
		std::sort(sorteddigittimes.begin(), sorteddigittimes.end());
		trackhittimesv.push_back(sorteddigittimes.at(0));
		for(Int_t i=0;i<sorteddigittimes.size()-1;i++){
			Float_t timetonextdigit = sorteddigittimes.at(i+1)-sorteddigittimes.at(i);
			if(timetonextdigit>maxtrackduration){
				trackhittimesv.push_back(sorteddigittimes.at(i+1));
				//cout<<"Setting track time threshold at "<<trackhittimesv.back()<<endl;
			}
		}
		//cout<<trackhittimesv.size()<<" tracks this event"<<endl;
		
		// a vector to record the track number for each hit, to know if we've allocated it yet.
		std::vector<Int_t> tracknumthisevent(mrddigittimesthisevent.size(),-1);
		
		// now we need to sort the digits into the tracks they belong to
		for(Int_t j=0; j<trackhittimesv.size(); j++){
			//cout<<"Track struck MRD at = "<<trackhittimesv.at(j)<<"ns in event "<<entry<<endl;
			// don't need to worry about lower bound as we start from lowest t peak and exclude already allocated hits
			
			// hit times are not ordered, so scan through them all 
			for(Int_t thisdigit=0;thisdigit<mrddigittimesthisevent.size();thisdigit++){
				Float_t endtime = (j<(trackhittimesv.size()-1)) ? trackhittimesv.at(j+1) : (eventendtime+1.);
				if(tracknumthisevent.at(thisdigit)<0 && mrddigittimesthisevent.at(thisdigit)< endtime ){
					tracknumthisevent.at(thisdigit)=j;
					WCSimRootCherenkovDigiHit* thedigihit = (WCSimRootCherenkovDigiHit*)atrigm->GetCherenkovDigiHits()->At(thisdigit);
					Int_t thisdigitstubeid = thedigihit->GetTubeId();
					tubeidsinatrack.push_back(thisdigitstubeid);
					Int_t thisdigitsq = thedigihit->GetQ();
					digitqsinatrack.push_back(thisdigitsq);
					double thisdigitstime = thedigihit->GetT();
					digittimesinatrack.push_back(thisdigitstime);
					// add all the unique parent ID's for digits contributing to this track (truth level info)
					std::vector<int> truephotonindices = thedigihit->GetPhotonIds();
					for(int truephoton=0; truephoton<truephotonindices.size(); truephoton++){
						int thephotonsid = truephotonindices.at(truephoton);
						WCSimRootCherenkovHitTime *thehittimeobject = (WCSimRootCherenkovHitTime*)atrigm->GetCherenkovHitTimes()->At(thephotonsid);
						int thephotonsparenttrack = thehittimeobject->GetParentID();
						int checkcount = std::count(particleidsinatrack.begin(), particleidsinatrack.end(), thephotonsparenttrack);
						if(checkcount==0){ particleidsinatrack.push_back(thephotonsparenttrack); }
					}
				}
			}
			tubeidsinthistrack.push_back(tubeidsinatrack);
			digitqsinthistrack.push_back(digitqsinatrack);
			digittimesinthistrack.push_back(digittimesinatrack);
			particleidsinthistrack.push_back(particleidsinatrack);
			
			tubeidsinatrack.clear();
			digitqsinatrack.clear();
			digittimesinatrack.clear();
			particleidsinatrack.clear();
		}
		
		// quick scan to check for any unallocated hits
		for(Int_t k=0;k<tracknumthisevent.size();k++){
			if(tracknumthisevent.at(k)==-1){cout<<"*****unbinned hit!"<<k<<" "<<mrddigittimesthisevent.at(k)<<endl;}
		}
		
		nummrddigitsthisevent = mrddigittimesthisevent.size();
		nummrdtracksthisevent=trackhittimesv.size();
		
		nummrddigitsthiseventb->Fill();
		nummrdtracksthiseventb->Fill();
		tubeidsinthistrackb->Fill();
		digitqsinthistrackb->Fill();
		digittimesinthistrackb->Fill();
		particleidsinthistrackb->Fill();
		
	}
	
	// WRITE+CLOSE OUTPUT FILES
	// ========================
	outfile->cd();
	recotree->SetEntries(nummrddigitsthiseventb->GetEntries());
	recotree->Write("",kOverwrite);
	outfile->Close();
	delete outfile;
}

// #######################################################################

// USE COLLECTIONS OF HITS TO FORM TRACKS
// ======================================
void WCSimAnalysis::FindMRDtracksInEvent(){
// TODO: fill this in based on mrd track reco.
}
