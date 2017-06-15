/* vim:set noexpandtab tabstop=4 wrap */
#ifndef _cMRD_reconstruction_classes
#define _cMRD_reconstruction_classes 1

#ifdef __MAKECINT__
#pragma link C++ class std::pair<mrdcluster*,mrdcluster*>+;
#endif
// supporting classes used in DoReconstruction cellular algorithm version (cMRDTrack_DoReconstruction2.cxx)
// 1. define a cluster, by an id and it's position
class mrdcluster : public TObject{
	public: 
	static Int_t clustercounter;
	Int_t clusterid;
	Int_t xmaxid;
	Int_t xminid; 
	Int_t layer;
	Int_t altviewhalf;
	std::vector<Int_t> in_layer_tubeids;
	std::vector<Double_t> digittimes;
	std::vector<Int_t> digitindexes;  //index of the digit within the WCSimRootTrigger
	void AddDigit(Int_t digitindexin, Int_t pmtidin, Double_t timein){
		digitindexes.push_back(digitindexin);
		digittimes.push_back(timein);
		
		// calculate in-layer tube index - first subtract id of first tube in the layer
		pmtidin -= layeroffsets.at(layer);
		// correct for the fact there are two halves to each layer. 
		// tubes 0 and 1 are in opposite halves!!
		// if pmtidin%2==1, it is in the second half. Merge halves by subtracting 1, then divide by 2. 
		if(pmtidin%2==1){
			if(altviewhalf>0) altviewhalf=0; // added digit is in the other half to constructor digit
			pmtidin -= 1;
		} else {
			if(altviewhalf<0) altviewhalf=0;
		}
		pmtidin /= 2;
		
		if(pmtidin>(*std::max_element(in_layer_tubeids.begin(), in_layer_tubeids.end()))){
			// this PMT is 'higher' than the others - increase the xmax accordingly
			xmaxid = pmtidin;
		} else {
			// this PMT is lower - decrease the xmin accordingly
			xminid = pmtidin;
		}
		
		in_layer_tubeids.push_back(pmtidin);
	}
	// return an *effective in-layer tube index* of the cluster centre.
	Double_t GetCentreIndex(){
		if(xmaxid==xminid){
			//cout<<"xmaxid==xminid, centreindex is "<<xmaxid<<endl;
			return xmaxid;
		} else {
			//cout<<"xmaxid!=xminid, centreindex is "<<((xmaxid+xminid)/2.)<<endl;
			return ((xmaxid+xminid)/2.);
		}
	}
	// return an actual cm in-layer position
	Double_t GetXmax(){
		//cout<<"GetXmax for cluster "<<clusterid<<endl;
		if (layer%2==0){	// h layer: short paddle extents are those in extentsy
			//cout<<"layer is a h layer, xmaxid is "<<xmaxid
			//	<<", getting upper paddle extent of pmt id "<<(2*xmaxid) + layeroffsets.at(layer)
			//	<<" = "<<paddle_extentsy.at((2*xmaxid) + layeroffsets.at(layer)).second<<endl;
			return paddle_extentsy.at((2*xmaxid) + layeroffsets.at(layer)).second;
		} else {
			//cout<<"layer is a v layer, xmaxid is "<<xmaxid
			//	<<", getting upper paddle extent of pmt id "<<(2*xmaxid) + layeroffsets.at(layer)
			//	<<" = "<<paddle_extentsx.at((2*xmaxid) + layeroffsets.at(layer)).second<<endl;
			return paddle_extentsx.at((2*xmaxid) + layeroffsets.at(layer)).second;
		}
	}
	Double_t GetXmin(){
		//cout<<"GetXmin for cluster "<<clusterid<<endl;
		if (layer%2==0){
			//cout<<"layer is a h layer, xminid is "<<xminid
			//	<<", getting lower paddle extent of pmt id "<<(2*xminid) + layeroffsets.at(layer)
			//	<<" = "<<paddle_extentsy.at((2*xminid) + layeroffsets.at(layer)).first<<endl;
			return paddle_extentsy.at((2*xminid) + layeroffsets.at(layer)).first;
		} else {
			//cout<<"layer is a v layer, xminid is "<<xminid
			//	<<", getting lower paddle extent of pmt id "<<(2*xminid) + layeroffsets.at(layer)
			//	<<" = "<<paddle_extentsx.at((2*xminid) + layeroffsets.at(layer)).first<<endl;
			return paddle_extentsx.at((2*xminid) + layeroffsets.at(layer)).first;
		}
	}
	Double_t GetCentre(){
		return ((this->GetXmax()+this->GetXmin())/2.);
	}
	Double_t GetTime(){
		return digittimes.at(0);
	}
	// constructor
	mrdcluster(Int_t digitidin, Int_t pmtidin, Int_t layerin, Double_t timein){
		digitindexes.push_back(digitidin);
		layer=layerin;
		clusterid=clustercounter;
		clustercounter++;
		digittimes.push_back(timein);
		//cout<<"generating cluster with tube "<<pmtidin<<endl;
		
		// calculate pmt number within this layer
		pmtidin -= layeroffsets.at(layer);
		// correct for the fact there are two halves to each layer. 
		// tubes 0 and 1 are in opposite halves!!
		// if pmtidin%2==1, it is in the second half. Merge halves by subtracting 1. 
		if(pmtidin%2==1){ altviewhalf=-1; pmtidin -= 1; }
		else { altviewhalf=1; }
		pmtidin /= 2;
		//cout<<"in-layer tube id is "<<pmtidin<<endl;
		
		in_layer_tubeids.push_back(pmtidin);
		
		// as the constructor takes the first digit, this is both the max and min hit tube
		xmaxid=pmtidin;
		xminid=pmtidin;
		
		// if this is the first constructor, fill the static members
		if(fillstaticmembers){
			//StripMrdPositions("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/MRD_positions_raw");
			fillstaticmembers=false;
//			for(size_t paddlei=0; paddlei<paddle_extentsx.size(); paddlei++){
//				cout<<"paddle_extentsx.at("<<paddlei<<") = "<<paddle_extentsx.at(paddlei).first
//					<<", "<<paddle_extentsx.at(paddlei).second<<endl;
//			}
//			assert(false);
		}
	}
	
	// default constructor 
	mrdcluster(){}
	
	// destructor
	~mrdcluster(){}
	
	// static members that store information about the paddles so they can be looked up from the tube ID
	static Bool_t fillstaticmembers;
	static std::vector<Int_t> paddle_orientations;	// H is 0, V is 1
	static std::vector<Int_t> paddle_layers;
	static std::vector<Double_t> paddle_originx;
	static std::vector<Double_t> paddle_originy;
	static std::vector<Double_t> paddle_originz;
	static std::vector<std::pair<Double_t,Double_t> > paddle_extentsx;
	static std::vector<std::pair<Double_t,Double_t> > paddle_extentsy;
	static std::vector<std::pair<Double_t,Double_t> > paddle_extentsz;
	
	// Function to build the information about the paddle positions from a wcsim output file
	static int StripMrdPositions(std::string filename);	// fills the member vectors paddle_*** 
};
#include "StripMrdPositions.C"	// function to pull information about paddle positions from file
Int_t mrdcluster::clustercounter=0;
Bool_t mrdcluster::fillstaticmembers=true;
std::vector<Int_t> mrdcluster::paddle_orientations(nummrdpmts);
std::vector<Int_t> mrdcluster::paddle_layers(nummrdpmts);
std::vector<Double_t> mrdcluster::paddle_originx(nummrdpmts);
std::vector<Double_t> mrdcluster::paddle_originy(nummrdpmts);
std::vector<Double_t> mrdcluster::paddle_originz(nummrdpmts);
std::vector<std::pair<Double_t,Double_t> > mrdcluster::paddle_extentsx(nummrdpmts);
std::vector<std::pair<Double_t,Double_t> > mrdcluster::paddle_extentsy(nummrdpmts);
std::vector<std::pair<Double_t,Double_t> > mrdcluster::paddle_extentsz(nummrdpmts);
int nothing = mrdcluster::StripMrdPositions();
// FIXME: RESULTS RETURNED FROM STRIPMRDPOSITIONS ARE IN MM!
// EVERYTHING ELSE HERE IS IN CM. THIS IS CONFUSING.

// ===================================================================================================
// 2. define a cell as an object containing an ID, a line-segment (pair of clusters), an integer state (initialised to 0), and a vector of neighbour cell ids.
class mrdcell : public TObject{
	public:
	static Int_t cellcounter;
	Int_t cellid;
	std::pair<mrdcluster*, mrdcluster*> clusters; //! used in building the cell
	Bool_t isdownstreamgoing;	// is the track going upstream instead of downstream
	Int_t status;
	Int_t utneighbourcellindex;	// the up-track neighbour of this cell (*)
	Int_t dtneighbourcellindex;	// the down-track neighbour of this cell (*)
								// a down-track cell index of -2 indicates more than one
	Int_t parentcellindex;		// if the up-track cell has more than one downtrack neighbour, the track splits (*)
								// so it's upstream neighbour becomes it's parent
	Bool_t hasdaughters;		// so we can keep short tracks that have daughters (*)
	Double_t neighbourchi2;		// if we have more than one candidate, use the most aligned neighbour (*)
								// (*) denotes fields that are only useful during reconstruction
	
	void IncrementStatus(){ status +=1; }
	void SetClusterAddresses(std::vector<mrdcluster> &trackclusters){
		if(clusters.first!=nullptr) return;
		clusters=std::make_pair(&trackclusters.at(cellid), &trackclusters.at(cellid+1));
	}
	
	mrdcell(mrdcluster* startcluster, mrdcluster* endcluster) : status(0), utneighbourcellindex(-1), 
		dtneighbourcellindex(-1), parentcellindex(-1), hasdaughters(false), neighbourchi2(-1.) {
		clusters = std::make_pair(startcluster, endcluster);
		isdownstreamgoing = (clusters.first->GetTime() < clusters.second->GetTime());
		cellid=cellcounter;
		cellcounter++;
	}
	
	// 'copy' constructor (also takes a new cell id)
	mrdcell(const mrdcell& cellin, int cellidin){
		cellid=cellidin;  // id within this track
		clusters = std::pair<mrdcluster*, mrdcluster*>(nullptr,nullptr);
		isdownstreamgoing=cellin.isdownstreamgoing;
		status=cellin.status;
		utneighbourcellindex=cellin.utneighbourcellindex;
		dtneighbourcellindex=cellin.dtneighbourcellindex;
		parentcellindex=cellin.parentcellindex;
		hasdaughters=cellin.hasdaughters;
		neighbourchi2=cellin.neighbourchi2;
	}
	
	// default constructor
	mrdcell(){}
	
	// destructor
	~mrdcell(){}
};
Int_t mrdcell::cellcounter=0;
#endif
