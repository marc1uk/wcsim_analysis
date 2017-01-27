/* vim:set noexpandtab tabstop=2 wrap */
// #######################################################################

// CREATE+OPEN OUTPUT FILE
// =======================
void WCSimAnalysis::OpenMRDtrackOutfile(){
	cout<<"opening mrd output file"<<endl;
	mrdtrackfile = new TFile("mrdtrackfile","RECREATE","MRD Tracks file");
	mrdtrackfile->cd();
	recotree = new TTree("recotree","Tree for reconstruction data");
	gROOT->cd();
	nummrddigitsthisevent=0;
	nummrdtracksthisevent=0;

	nummrddigitsthiseventb = recotree->Branch("nummrddigitsthisevent",&nummrddigitsthisevent);
	nummrdtracksthiseventb = recotree->Branch("nummrdtracksthisevent",&nummrdtracksthisevent);
	tracksinthiseventb = recotree->Branch("tracksinthisevent",&tracksinthisevent);
	
	// ?????????????
	TClonesArray* aDigit = new TClonesArray("cMRDdigit");	// string is class name
	TClonesArray* aTrack = new TClonesArray("cMRDTrack");	// string is class name
	TClonesArray &aTracka = *aTrack;
	TClonesArray &aDigita = *aDigit;
}

// #######################################################################

// SPLIT HITS BY TIME
// ==================
void WCSimAnalysis::FindMRDtracksInEvent(){

	if(numhitsthisevent==0){
		// =====================================================================================
		// NO TRACKS IN THIS EVENT
		// ======================================
		numtracksthisevent=0;
		new(aTracka[0]) cMRDTrack();
		//recotree->Fill();	// fill the tree anyway so entries align
		pmtnumtracksthiseventb->Fill();
		pmtnumdigitsthiseventb->Fill();
		pmtTracksBranch->Fill();
		continue;
	}	// skip remainder
	// ==========================================================================================
	

	aTracka.Clear();
	aDigita.Clear();
/* if you class contains pointers, use aTrack.Clear("C"). You MUST then provide a Clear() method in your class that properly performs clearing and memory freeing. (or "implements the reset procedure for pointer objects")
 see https://root.cern.ch/doc/master/classTClonesArray.html#a025645e1e80ea79b43a08536c763cae2
*/
	
	// MEASURE EVENT DURATION TO DETERMINE IF THERE IS MORE THAN ONE TRACK
	// ===================================================================
	std::vector<Int_t> tubeidsinatrack;
	std::vector<Double_t> digitqsinatrack;
	std::vector<Double_t> digittimesinatrack;
	std::vector<Int_t> particleidsinatrack;
	std::vector<Double_t> photontimesinatrack;
	
	// first check: are all hits within a 30ns window (maxtrackduration) If so, just one track. 
	Double_t eventendtime = *std::max_element(mrddigittimesthisevent.begin(),mrddigittimesthisevent.end());
	Double_t eventstarttime = *std::min_element(mrddigittimesthisevent.begin(),mrddigittimesthisevent.end());
	Double_t eventduration = (eventendtime - eventstarttime);
	//cout<<"event start: "<<eventstarttime<<", end : "<<eventendtime<<", duration : "<<eventduration<<endl;
	
	if(eventduration<maxtrackduration){
	// JUST ONE TRACK
	// ==============
		//cout<<"all hits this event within one track."<<endl;
		
		// loop over digits and convert them into cMRDdigit objects that contain all their information
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
				double thephotonstruetime = thehittimeobject->GetTruetime();
				photontimesinatrack.push_back(thephotonstruetime);
			}
			// construct the digit
			new(aDigita[thisdigit]) cMRDdigit(thisdigitstime, thisdigitsq, thisdigitstubeid, photontimesinatrack,particleidsinatrack);
		}
		
		// construct a track from all these digits
		new(aTracka[0]) cMRDTrack(aDigita);
		// can also use 'cMRDTrack* = (cMRDTrack*)aTrack.ConstructedAt(0);' followed by a bunch of 'Set' calls
		// to set all relevant fields. This bypasses the constructor, calling it only when necessary, 
		// saving time. In that case, we do not need to call aTracka.Clear(); 
		//recotree->Fill();
		
		nummrddigitsthisevent = mrddigittimesthisevent.size();
		nummrdtracksthisevent=1;
		tubeidsinthisevent.push_back(tubeidsinatrack);
		digitqsinthisevent.push_back(digitqsinatrack);
		digittimesinthisevent.push_back(digittimesinatrack);
		particleidsinthisevent.push_back(particleidsinatrack);
		photontimesinthisevent.push_back(photontimesinatrack);
		
		
		nummrddigitsthiseventb->Fill();
		nummrdtracksthiseventb->Fill();
		tubeidsinthiseventb->Fill();
		digitqsinthiseventb->Fill();
		digittimesinthiseventb->Fill();
		particleidsinthiseventb->Fill();
		photontimesinthiseventb->Fill();
		
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
		
		// now we need to sort the digits into the tracks they belong to:
		// loop over tracks
		for(Int_t thistrack=0; thistrack<trackhittimesv.size(); thistrack++){
			//cout<<"Track struck MRD at = "<<trackhittimesv.at(thistrack)<<"ns in event "<<entry<<endl;
			// don't need to worry about lower bound as we start from lowest t peak and exclude already allocated hits
			
			// hit times are not ordered, so scan through them all 
			for(Int_t thisdigit=0;thisdigit<mrddigittimesthisevent.size();thisdigit++){
				Float_t endtime = (thistrack<(trackhittimesv.size()-1)) ? trackhittimesv.at(thistrack+1) : (eventendtime+1.);
				if(tracknumthisevent.at(thisdigit)<0 && mrddigittimesthisevent.at(thisdigit)< endtime ){
					// thisdigit is in thistrack
					tracknumthisevent.at(thisdigit)=thistrack;
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
						double thephotonstruetime = thehittimeobject->GetTruetime();
						photontimesinatrack.push_back(thephotonstruetime);
					}
					// construct the digit
					new(aDigita[thisdigit]) cMRDdigit(thisdigitstime, thisdigitsq, thisdigitstubeid, photontimesinatrack,particleidsinatrack);
				}
			}
			tubeidsinthisevent.push_back(tubeidsinatrack);
			digitqsinthisevent.push_back(digitqsinatrack);
			digittimesinthisevent.push_back(digittimesinatrack);
			particleidsinthisevent.push_back(particleidsinatrack);
			photontimesinthisevent.push_back(photontimesinatrack);
			
			
			tubeidsinatrack.clear();
			digitqsinatrack.clear();
			digittimesinatrack.clear();
			particleidsinatrack.clear();
			photontimesinatrack.clear();
			
		}
		
		// quick scan to check for any unallocated hits
		for(Int_t k=0;k<tracknumthisevent.size();k++){
			if(tracknumthisevent.at(k)==-1){cout<<"*****unbinned hit!"<<k<<" "<<mrddigittimesthisevent.at(k)<<endl;}
		}
		
		nummrddigitsthisevent = mrddigittimesthisevent.size();
		nummrdtracksthisevent=trackhittimesv.size();
		
		nummrddigitsthiseventb->Fill();
		nummrdtracksthiseventb->Fill();
		tubeidsinthiseventb->Fill();
		digitqsinthiseventb->Fill();
		digittimesinthiseventb->Fill();
		particleidsinthiseventb->Fill();
		photontimesinthiseventb->Fill();
		
	}
	
	// WRITE+CLOSE OUTPUT FILES
	// ========================
	mrdtrackfile->cd();
	recotree->SetEntries(nummrddigitsthiseventb->GetEntries());
	recotree->Write("",kOverwrite);
	gROOT->cd();
}

// #######################################################################

// USE COLLECTIONS OF HITS TO FORM TRACKS
// ======================================
void WCSimAnalysis::FindMRDtracksInEvent(){
	// TODO: fill this in based on mrd track reco.
//	tubeidsinthistrack
//	digitqsinthistrack
//	digittimesinthistrack
//	particleidsinthistrack

Int_t numhitsthisevent=0;
Int_t numhitsthiseventtemp=0;
Int_t numtracksthisevent;

TCanvas* splittracksCanv = new TCanvas("splittracksCanv","Title");
TH1F* splittrackshist=0;

Bool_t* hitunallocated = new Bool_t[maxtracksperevent];

TClonesArray &aTracka = *aTrack;
std::vector<Double_t> hittimesthisdigit;
std::vector<cMRDdigit> digitsthistrack;

TBranch* pmtnumtracksthiseventb = recotree->Branch("pmtnumtracksthisevent",&numtracksthisevent);
TBranch* pmtnumhitsthiseventb = recotree->Branch("pmtnumhitsthisevent",&numhitsthisevent);
TBranch* pmtTracksBranch = recotree->Branch("PMTrecotracks",&aTrack, maxtracksperevent); 
//TBranch* pmtTracksBranch  = recotree->Branch("PMTrecotracks[pmtnumtracksthisevent]/cMRDTrack",&aTrack);

std::vector<Double_t>* hittimesbydigit = new std::vector<Double_t>[mrdnumlayers];
std::vector<Int_t>* copynumsbydigit = new std::vector<Int_t>[mrdnumlayers];

aTracka.Clear();
// if you class contains pointers, use aTrack.Clear("C"). You MUST then provide a Clear() method in your class that properly
// performs clearing and memory freeing. (or "implements the reset procedure for pointer objects") 
// https://root.cern.ch/doc/master/classTClonesArray.html#a025645e1e80ea79b43a08536c763cae2
numhitsbranch->GetEntry(entry); // first check if we have any hits this entry. 
numhitsthisevent = numhitsthiseventtemp;	// numhitsthiseventtemp WILL GET OVERWRITTEN BY THE DRAW() CALL BELOW!

if(numhitsthisevent==0){
	// =====================================================================================
	// NO TRACKS IN THIS EVENT
	// ======================================
	numtracksthisevent=0;
	new(aTracka[0]) cMRDTrack();
	//recotree->Fill();	// fill the tree anyway so entries align
	pmtnumtracksthiseventb->Fill();
	pmtnumhitsthiseventb->Fill();
	pmtTracksBranch->Fill();
	continue;
}	// skip remainder
// =========================================================================================

// AT LEAST ONE TRACK THIS EVENT
// ======================================
// if continued to here, we have at least some hits. 
hittimebranch->ResetAddress(); 	// so we don't have to ensure sufficient memory allocation for 'Draw()' call.
trackidsbranch->ResetAddress();
copynumbranch->ResetAddress();

// isolate each track by hit time - multiple hits from the same track generally all occur within 1ns span 
// (1 unit in G4 native)
// need to plot histogram with sufficient resolution that there will be at least one empty bin between
// hits from adjacent tracks
if(splittrackshist!=0){delete splittrackshist; splittrackshist=0;}
pmttree->Draw("pmthit_t>>splittrackshist",TString::Format("Entry$==%d",entry));
splittrackshist = (TH1F*)splittracksCanv->GetPrimitive("splittrackshist");

// first check: are all hits within a 30ns window (maxtrackduration) If so, just one track. 
Double_t eventendtime = splittrackshist->GetXaxis()->GetXmax();
Double_t eventstarttime = splittrackshist->GetXaxis()->GetXmin();
Double_t eventduration = (eventendtime - eventstarttime);
//cout<<"event start: "<<eventstarttime<<" end : "<<eventendtime<<endl;
//cout<<"event duration is "<<eventduration<<endl;
if(eventduration<maxpmttrackduration){
	// ONE TRACK THIS EVENT
	// ======================================
	//cout<<"all hits this event within one track."<<endl;
	// all hits are within a short period of time - assume they are from the same track.
	numtracksthisevent=1;
	hittimesthisdigit.clear();
	hittimes = new Double_t[numhitsthisevent];
	hittimebranch->SetAddress(hittimes);
	hittimebranch->GetEntry(entry);
	hittimesthisdigit.assign(hittimes,hittimes+numhitsthisevent);
	copynums = new Int_t[numhitsthisevent];	
	copynumbranch->SetAddress(copynums);
	copynumbranch->GetEntry(entry);
	Int_t PMTmosthit = findmode(copynums, numhitsthisevent);
	cMRDdigit aMRDdigit = cMRDdigit(hittimesthisdigit,PMTmosthit);
	digitsthistrack.clear();
	digitsthistrack.push_back(aMRDdigit);
	new(aTracka[0]) cMRDTrack(digitsthistrack);
	// can also use 'cMRDTrack* = (cMRDTrack*)aTrack.ConstructedAt(0);' followed by a bunch of 'Set' calls to 
	// set all relevant fields. This bypasses the constructor, calling it only when necessary, saving time.
	// in this case we do not need to call aTracka.Clear(); 
	//recotree->Fill();
	pmtnumtracksthiseventb->Fill();
	pmtnumhitsthiseventb->Fill();
	pmtTracksBranch->Fill();

} else {
	// =========================================================================================
	// MULTIPLE TRACKS IN THIS EVENT
	// ======================================
	// this event has multiple tracks. Need to split hits into which track they belong to.
	// first, to keep track of what's been assigned, initialize all to -1.		
	for(Int_t j=0;j<numhitsthisevent;j++){hitunallocated[j]=true;}
	
	// Now ensure binning is fine enough to have at least one empty bin between tracks.
	Int_t numbins = TMath::Floor(eventduration)/50;
	if(splittrackshist!=0){delete splittrackshist;splittrackshist=0;}
	pmttree->Draw(TString::Format("pmthit_t>>splittrackshist(%d,%f,%f)",numbins,eventstarttime,eventendtime),TString::Format("Entry$==%d",entry));
	// unless we explicitly give the start/end times ROOT decides to arbitrarily clip the range when specifying # bins
	splittrackshist = (TH1F*)splittracksCanv->GetPrimitive("splittrackshist");
	//cout<<"using "<<numbins<<", "<<splittrackshist->GetNbinsX()<<" bins of width "<<splittrackshist->GetXaxis()->GetBinWidth(1)<<endl;

	bool inpeak=false;
	std::vector<Float_t> trackhittimesv;
	for(Int_t i=0;i<numbins;i++){
		if(splittrackshist->GetBinContent(i)>0){inpeak=true;}
		else {if(inpeak==true){
			trackhittimesv.push_back(splittrackshist->GetXaxis()->GetBinLowEdge(i));
			//cout<<"Setting track time threshold at "<<trackhittimesv.back()<<endl;
			inpeak=false;}
		}
	}
	numtracksthisevent = trackhittimesv.size();
  		cout<<numtracksthisevent<<" tracks and "<<numhitsthisevent<<" hits in event "<<entry<<endl;
	
	hittimes = new Double_t[numhitsthisevent];
	hittimebranch->SetAddress(hittimes);
	hittimebranch->GetEntry(entry);
	trackids = new Int_t[numhitsthisevent];	// only used for checking
	trackidsbranch->SetAddress(trackids);
	trackidsbranch->GetEntry(entry);
	copynums = new Int_t[numhitsthisevent];
	copynumbranch->SetAddress(copynums);
	copynumbranch->GetEntry(entry);

	// for each track time threshold, scan through all hits and collect those in this track
	// collect hits in a given track into an array of vectors; each array element respresents an MRD layer/potential digit
	for(Int_t layer=0;layer<mrdnumlayers;layer++){
		hittimesbydigit[layer].clear();
		copynumsbydigit[layer].clear();
	}

	for(Int_t j=0; j<numtracksthisevent; j++){
		digitsthistrack.clear();
		//cout<<"processing track "<<j<<endl;
		// hit times are not ordered, so scan through them all 
		for(Int_t thishit=0;thishit<numhitsthisevent;thishit++){
			if(hitunallocated[thishit]&&hittimes[thishit]<trackhittimesv.at(j)){
				//cout<<"allocating hit "<<thishit<<endl;
				hitunallocated[thishit]=false;
				// timing gives the track - now check z value to determine appropriate digit
				Int_t panelnum = TMath::Floor(copynums[thishit]/mrdpaddlesperpanel);
				hittimesbydigit[panelnum].push_back(hittimes[thishit]);
				copynumsbydigit[panelnum].push_back(copynums[thishit]);
				//if(entry==94){cout<<"put hit "<<thishit<<" into track "<<j<<" of this event;";
				//cout<<" hit time: "<<hittimes[thishit]<<endl;}
			}
		}

		// we have an array of 12 potential digits; check which have hits and for those with >0 hits,
		// create a digit object.
		for(Int_t layer=0;layer<mrdnumlayers;layer++){
			//cout<<"forming digit for layer "<<layer<<endl;
			// each layer with hits in is a digit
			if(hittimesbydigit[layer].size()>0){
				cout<<hittimesbydigit[layer].size()<<" hits in digit in layer "<<layer;
				// create a digit for this layer. 
				hittimesthisdigit.clear();
				hittimesthisdigit=hittimesbydigit[layer];
				std::vector<Int_t> copynumsthisdigit = copynumsbydigit[layer];
				Int_t PMTmosthit = findmode(&copynumsthisdigit.front(),hittimesthisdigit.size());
				cout<<", PMT most hit: "<<PMTmosthit<<endl;
				cMRDdigit aMRDdigit = cMRDdigit(hittimesthisdigit,PMTmosthit);
				digitsthistrack.push_back(aMRDdigit);
			} //else {cout<<"no digits this layer"<<endl;}
		}
		new(aTracka[j]) cMRDTrack(digitsthistrack);
		cout<<"creating track of "<<digitsthistrack.size()<<" digits."<<endl;
	}
	
	// quick scan to check for any unallocated hits (i.e. in a peak not detected by TSpectrum search)
	for(Int_t k=0;k<numhitsthisevent;k++){
		if(hitunallocated[k]){cout<<"*****unbinned hit!"<<k<<" "<<hittimes[k]<<endl;}
	}
	//cout<<"filling tree event "<<entry<<" with "<<numhitsthisevent<<" hits"<<endl;
	//recotree->Fill();
	pmtnumtracksthiseventb->Fill();
	pmtnumhitsthiseventb->Fill();
	pmtTracksBranch->Fill();
	delete[] hittimes;
	delete[] trackids;
	delete[] copynums;
}	// more than one track


delete[] hittimesbydigit;	// using delete rather than delete[] results in segfault!
delete[] copynumsbydigit;
delete splittracksCanv;
pmttree->ResetBranchAddresses();
// NOTE: WHEN STORING AN OBJECT IN A BRANCH: The object must not be destroyed (i.e. be deleted) until the TTree is deleted or TTree::ResetBranchAddress is called. 
cout<<"resetting recotree"<<endl;
recotree->ResetBranchAddresses(); // so this MUST be called BEFORE aTrack gets deleted!
//cout<<"deleting aTrack"<<endl;	// doesn't seem necessary for TClonesArray? 
//delete[] aTrack;		// using delete rather than delete[] results in segfault!
cout<<"exiting"<<endl;

}

}
