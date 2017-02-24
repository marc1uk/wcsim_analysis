/* vim:set noexpandtab tabstop=4 wrap */
#ifndef __CINT__
#include "Riostream.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TGraphErrors.h"
#include "TDecompChol.h"
#include "TDecompSVD.h"
#include "TF1.h"
#include "TBox.h"
#endif
#include "MRDSubEventClass.hh"

#ifndef TRACKFINDVERBOSE
//#define TRACKFINDVERBOSE 1
#endif

class mrdcluster;
class mrdcell;
const Int_t mintracklength=1;	// tracks must have (mintracklength+1) cells to qualify
const Double_t chi2limit=40.0;	// max chi^2 from a linear least squares fit to all 3 clusters for those
								// to be classed as clusters in the same track

void cMRDSubEvent::DoReconstruction(){
#ifdef __CINT__
	gSystem->Load("libMatrix");
#endif
//	// get digit tubes, layers and energy depositions. (not part of track reconstruction)
//	for(Int_t i=0;i<digi_ids.size();i++){
//		Int_t tube_id = pmts_hit.at(i);
//		Int_t strucklayer = mrdcluster::paddle_layers.at(tube_id);
//		if(std::count(layers_hit.begin(), layers_hit.end(), strucklayer)==0){
//			layers_hit.push_back(strucklayer);
//		}
//		eDepsInLayers.at(strucklayer)+= (digi_qs.at(i));	// TODO: convert q to energy?
//	}
	
	// 1. define a cluster, by an id and it's position. - done in cMRDSubEvent_ReconstructionClasses
	std::vector<mrdcluster*> theclusters;
	mrdcluster::clustercounter=0;	// reset for the start of a new event
	
	//2. define a cell as an object containing an ID, a pair of clusters, an integer state (initialised to 0),
	//   and a vector of neighbour cell ids. - done in cMRDSubEvent_ReconstructionClasses
	std::vector<mrdcell*> thecells;
	mrdcell::cellcounter=0;	// reset for the start of a new event
	
	// 3. find large clusters. where there are multiple hits in adjacent paddles within the same layer, these combine into a single cluster
	// sort digits into layers
#ifdef TRACKFINDVERBOSE
	cout<<"sorting digits into layers"<<endl;
#endif
	std::vector< std::vector<int> > digits_by_layer( (numpanels) );	//note the inner parenthesis are rqd!
	for(Int_t adigiindex=0; adigiindex<digi_ids.size(); adigiindex++){
		Int_t tube_id = pmts_hit.at(adigiindex);
		Int_t strucklayer = mrdcluster::paddle_layers.at(tube_id);
		digits_by_layer.at(strucklayer).push_back(adigiindex);
	}
#ifdef TRACKFINDVERBOSE
	for(int layer=0; layer<numpanels; layer++){
		cout<<"found "<<digits_by_layer.at(layer).size()<<" digits in layer "<<layer<<endl;
	}
#endif
	
	// first make note of clusters spanning two or more adjacent pmts in the same layer
	// loop over layers, looking for multiple digits in adjacent pmts
	// VVVV this has been skipped. Particles can't hit two paddles in the MRD, showers should be separate tracks
	std::vector<Int_t> donedigits;
	/*
	for(int thislayer=0; thislayer<numpanels; thislayer++){
		if(digits_by_layer.at(thislayer).size()==0) continue;
		std::vector<Int_t> digitsinthislayer = digits_by_layer.at(thislayer);
		// it will be easier to perform a linear scan across the layer looking for adjacent pmts
		// this requires that digits in the layer are sorted by tube_id in ascending order
		
		
		// make a vector of pairs of digits in this layer, pair is <tube_id, digi_index>
		// then sort this vector - it will be sorted by tube_id.
		// then scan through the vector of pairs, comparing tube_ids. if match, use digi_index
		// part to retrieve the other required information
#ifdef TRACKFINDVERBOSE
		cout<<"building vector of pairs in layer "<<thislayer<<" for sorting"<<endl;
		cout<<"digits in this layer include:"<<endl;
#endif
		std::vector< std::pair<Int_t, Int_t> > digitssortedbytube;
		for(int j=0; j<digitsinthislayer.size(); j++){
			Int_t adigiindex = digitsinthislayer.at(j);
			Int_t atubeid = pmts_hit.at(adigiindex);
#ifdef TRACKFINDVERBOSE
			cout<<"   tube "<<atubeid<<", digit "<<adigiindex<<endl;
#endif
			digitssortedbytube.push_back(std::pair<Int_t, Int_t>(atubeid, adigiindex));
		}
#ifdef TRACKFINDVERBOSE
		cout<<"sorting pairs in this layer by tube id"<<endl;
#endif
		std::sort(digitssortedbytube.begin(), digitssortedbytube.end());
#ifdef TRACKFINDVERBOSE
		cout<<"sorted set of pairs is:"<<endl;
		for(int pairveci=0; pairveci<digitssortedbytube.size(); pairveci++){
			cout<<"   tube "<<digitssortedbytube.at(pairveci).first
				<<", digit "<<digitssortedbytube.at(pairveci).second<<endl;
		}
#endif
		
		// scan through all digits except the last one, in a forward direction
		// 1. look at the next digit's tube. If it's adjacent to this digit's tube, they're a cluster.
		// 2. if this digit is already in a cluster, add the next digit.
		// 3. otherwise, create a cluster with this digit and add the next digit.
#ifdef TRACKFINDVERBOSE
		cout<<"searching for multi-paddle clusters"<<endl;
#endif
		mrdcluster* thisdigitscluster=nullptr;
		for(Int_t thisdigit=0; thisdigit<(digitsinthislayer.size()-1); thisdigit++){
			// get digit index
			Int_t tube_id = digitssortedbytube.at(thisdigit).first;
			Int_t tube_id2 = digitssortedbytube.at(thisdigit+1).first;
			// if the next digit's tube_id is adjacent to that of this digit
			if(tube_id2==(tube_id+2)){	// offset is 2 because +1 is in the other half of the MRD
				// they're a cluster! grab the digi_ids
			Int_t adigiindex = digitssortedbytube.at(thisdigit).second;
			Int_t adigiindex2 = digitssortedbytube.at(thisdigit+1).second;
				//check if this digit is already in a cluster
				if(thisdigitscluster==nullptr){
					// it's not: make a cluster and add both digits
#ifdef TRACKFINDVERBOSE
					cout<<"creating multi-paddle cluster in layer "<<thislayer
						<<" from digits in tubes "<<tube_id<<" and "<<tube_id2<<endl;
#endif
					thisdigitscluster = new mrdcluster(adigiindex, tube_id, thislayer, digi_ts.at(adigiindex));
					theclusters.push_back(thisdigitscluster);
					thisdigitscluster->AddDigit(adigiindex2, tube_id2, digi_ts.at(adigiindex2));
					donedigits.push_back(adigiindex);
					donedigits.push_back(adigiindex2);
				} else {
					// it is: add the next digit to this digit's cluster
#ifdef TRACKFINDVERBOSE
					cout<<"adding digiti in tube "<<tube_id2<<" to this cluster"<<endl;
#endif
					thisdigitscluster->AddDigit(adigiindex2, tube_id2, digi_ts.at(adigiindex2));
					donedigits.push_back(adigiindex2);
				}
			} else {
				thisdigitscluster=nullptr;
			}
		}
#ifdef TRACKFINDVERBOSE
		cout<<"found "<<theclusters.size()<<" multi-digit clusters in this layer"<<endl;
#endif
	} */
	
	// 4. make any remaining non-grouped hits their own cluster.
#ifdef TRACKFINDVERBOSE
	cout<<"constructed "<<theclusters.size()<<" multi-paddle clusters."<<endl
		<<"Making remaining digits into clusters"<<endl;
#endif
	for(Int_t adigiindex=0; adigiindex<digi_ids.size(); adigiindex++){
		if(std::count(donedigits.begin(), donedigits.end(), adigiindex)==0){
			Int_t tube_id = pmts_hit.at(adigiindex);
			Int_t strucklayer = mrdcluster::paddle_layers.at(tube_id);
			mrdcluster* thisdigitscluster = new mrdcluster(adigiindex, tube_id, strucklayer, digi_ts.at(adigiindex));
			theclusters.push_back(thisdigitscluster);
		}
	}
#ifdef TRACKFINDVERBOSE
	cout<<"constructed "<<theclusters.size()<<" total clusters from "<<digi_ids.size()<<" digits"<<endl;
#endif
	// 5. generate a std::vector all the possible cells for each layer
	for(int thislayer=0; thislayer<numpanels; thislayer++){
	// loop over all clusters in the layer ... 
#ifdef TRACKFINDVERBOSE
		cout<<"constructing cells from clusters starting in layer "<<thislayer<<endl;
#endif
		for(int thiscluster=0; thiscluster<theclusters.size(); thiscluster++){
			mrdcluster* startcluster = theclusters.at(thiscluster);
			if(startcluster->layer==thislayer){
				std::vector<mrdcell*> thisclusterscells;
				// ... and generate a cell between it and all clusters in the next and next-next layer
				// - but separating the two views. So 'next' is 'next horizontal' or 'next vertical' as appopriate.
				for(int thiscluster2=0; thiscluster2<theclusters.size(); thiscluster2++){
					mrdcluster* endcluster = theclusters.at(thiscluster2);
					Int_t nextlayer = endcluster->layer;
					if((nextlayer>thislayer)&&((nextlayer-thislayer)<3)&&((nextlayer-thislayer)%2==0)){
						// if making a cell between two adjacent layers, make the cell - 
						// no need to check for equivalent cells with intermediate layers
						if(((nextlayer-thislayer)==2)||(thisclusterscells.size()==0)){
#ifdef TRACKFINDVERBOSE
							cout<<"generating a cell between cluster "<<thiscluster
								<<" with centre "<<startcluster->GetCentreIndex()<<" in layer "<<thislayer
								<<" and cluster "<<thiscluster2<<" with centre "<<endcluster->GetCentreIndex()
								<<" in layer "<<nextlayer<<endl;
#endif
							// generate the cell between these thiscluster and thiscluster2
							mrdcell* thiscell = new mrdcell(startcluster, endcluster);
							thecells.push_back(thiscell);
							thisclusterscells.push_back(thiscell);
						} else {
							// check if a cell already exists starting from startcluster and going to 
							// an intermediate cluster in perfect alignment with this endcluster
							// if so, we will skip making this cell, as it will overlap with the existing cell
#ifdef TRACKFINDVERBOSE
							cout<<"scanning for intermediate cells between cluster "<<thiscluster
								<<" with centre "<<startcluster->GetCentreIndex()<<" in layer "<<thislayer
								<<" and cluster "<<thiscluster2<<" with centre "<<endcluster->GetCentreIndex()
								<<" in layer "<<nextlayer<<endl;
#endif
							Bool_t foundintermediate=false;
							// scan through existing cells starting from this cluster
							for(int acelli=0; acelli<thisclusterscells.size(); acelli++){
								mrdcell* acell = thisclusterscells.at(acelli);
								// skip if the cell end cluster is in nextlayer - it's not intermediate
								if(acell->clusters.second->layer==nextlayer) {
#ifdef TRACKFINDVERBOSE
									cout<<"skipping cell to present layer (can't be intermediate)"<<endl;
#endif
									continue;
								}
								// existing cell is in the intermediate layer - 
								// does the projection to nextlayer hit this cluster?
								double startclustercentre = startcluster->GetCentreIndex();
								// need to make sure we extract the downstream (not necessarily second) cluster
								mrdcluster* secondcluster = acell->clusters.second;
								double midclustercentre = secondcluster->GetCentreIndex();
								int xdiff = midclustercentre-startclustercentre;
								double projectedendpoint = startclustercentre+(2*xdiff);
#ifdef TRACKFINDVERBOSE
								cout<<"projection of cell from tube "<<startclustercentre
									<<" to "<<midclustercentre<<" to the subsequent layer gives tube "
									<<projectedendpoint<<", compared to present cluster's centre of "
									<<endcluster->GetCentreIndex()<<endl;
#endif
								if(projectedendpoint==endcluster->GetCentreIndex()) {
									foundintermediate=true;  // same as existing cell
#ifdef TRACKFINDVERBOSE
									cout<<"this is an intermediate, skipping this endpoint cluster."<<endl;
#endif
									break;
								}
							}
							if(!foundintermediate){
#ifdef TRACKFINDVERBOSE
								cout<<"no intermediate cells found, creating a new cell"<<endl;
#endif
								// no existing cell that overlaps with this one; make the cell
								mrdcell* thiscell = new mrdcell(startcluster, endcluster);
								thecells.push_back(thiscell);
								thisclusterscells.push_back(thiscell);
							}
						}
					}
				}
			}
		}
	}
#ifdef TRACKFINDVERBOSE
	cout<<endl<<"constructed a total of "<<thecells.size()<<" cells"<<endl;
#endif
	
	//6. fill the neighbourhood for each cell. the neighbours are defined considering a pair of cells:
	//  a. the two cells must have a shared cluster (there will be two end clusters and a shared middle cluster)
	//  b. the clusters must be sufficiently aligned - the x^2 of a linear fit through all 3 cluster centres,
	// with errors defined by paddle boundaries, must be less than chi2limit.
	//  c. then there's some tricky stuff based on track splitting and converging. 
	// if an uptrack cell already has a downtrack cell, then another downtrack candidate is found, then the
	// track splits. In this case, we unset the previous downtrack candidate as a neighbour and instead make
	// both downtrack cells "daughters".
	// if a downtrack cell is found to align with multiple uptrack cells, then the uptrack cells converge.
	// in this case, we take the uptrack candidate with the best chi2 as the uptrack cell, and ignore others.
	//  d. finally, we do the scan from upstream to downstream, but the definition of "uptrack" and "downtrack"
	// must account for upstream directed tracks. We use the cell direction to account for this. 
	
	// loop over cells
	for(int acelli=0; acelli<thecells.size(); acelli++){
		mrdcell* upcell = thecells.at(acelli);
		mrdcluster* sharedcluster1 = upcell->clusters.second;
#ifdef TRACKFINDVERBOSE
		cout<<"searching for downstream neighbours of cell "<<acelli
			<<" with downstream cluster at "<<sharedcluster1<<endl;
#endif
		mrdcell* thiscellsneighbour = nullptr;
		// scan over other downstream cells
		for(int bcelli=(acelli+1); bcelli<thecells.size(); bcelli++){
			mrdcell* downcell = thecells.at(bcelli);
			mrdcluster* sharedcluster2 = downcell->clusters.first;
#ifdef TRACKFINDVERBOSE
			cout<<"checking cell with upstream cluster at "<<sharedcluster2<<endl;
#endif
			// check if they share a cluster
			if(sharedcluster1==sharedcluster2){
#ifdef TRACKFINDVERBOSE
				cout<<"shared cluster! Checking if they're aligned enough to be a track"<<endl;
#endif
				// they share a central cluster. Next check they are sufficiently 'aligned' to be considered
				// part of the same track. Do this by checking a x^2 fit to a straight line through all clusters
				mrdcluster* upcluster = upcell->clusters.first;
				mrdcluster* downcluster = downcell->clusters.second;
				// make a linear-least-square fit to the three cluster centres, with errors on the three
				// centres based on the cluster paddle limits
#ifdef TRACKFINDVERBOSE
				cout<<"checking clusters at "<<upcluster<<", "<<sharedcluster1<<", "<<downcluster<<endl;
#endif
				Double_t clusterzcentres[] = 
					{(Double_t)upcluster->layer, (Double_t)sharedcluster1->layer, (Double_t)downcluster->layer};
				Double_t clusterxcentres[] = 
					{upcluster->GetCentre(), sharedcluster1->GetCentre(), downcluster->GetCentre()};
				Double_t upclustererror = upcluster->GetXmax() - upcluster->GetXmin();
				Double_t sharedclustererror = sharedcluster1->GetXmax() - sharedcluster1->GetXmin();
				Double_t downclustererror = downcluster->GetXmax() - downcluster->GetXmin();
				Double_t clusterxerrors[] = {upclustererror, sharedclustererror, downclustererror};
				Double_t fit_gradient=0, fit_offset=0, chi2 = 0.;
				// is errors the same as weighting in chi^2 implementation?
#ifdef TRACKFINDVERBOSE
				cout<<"Performing least squares fit with centres at y values ("
					<<clusterxcentres[0]<<", "<<clusterxcentres[1]<<", "<<clusterxcentres[2]<<"), y errors ("
					<<clusterxerrors[0]<<", "<<clusterxerrors[1]<<", "<<clusterxerrors[2]<<"), z values ("
					<<clusterzcentres[0]<<", "<<clusterzcentres[1]<<", "<<clusterzcentres[2]<<")"<<endl;
#endif
				LeastSquaresMinimizer(3, clusterzcentres, clusterxcentres, 0, clusterxerrors, 
					fit_gradient, fit_offset, chi2);
#ifdef TRACKFINDVERBOSE
				cout<<"fit produces gradient "<<fit_gradient<<", offset "<<fit_offset<<", chi^2 "<<chi2<<endl;
#endif
				if(chi2<chi2limit){
#ifdef TRACKFINDVERBOSE
					cout<<"chi2 is good, checking for existing neighbours and daughters"<<endl;
#endif
					// these cells are neighbours! 
					// determine the track direction to see how to proceed
					if(upcell->isdownstreamgoing&&downcell->isdownstreamgoing){
						// both are cells of a downstream going track.
						// now we need to check if this cell is already the upstream neighbour of
						// another cell - i.e. does the track split here
						if((upcell->dtneighbourcellindex==-2)||(upcell->dtneighbourcellindex>-1)){
#ifdef TRACKFINDVERBOSE
							cout<<"this is a second downtrack candidate! Track splits!"<<endl;
#endif
							// we need to un-set the existing neighbour. Instead set both downstream
							// candidates as their own new tracks, and call this track the 'parent'
							if(upcell->dtneighbourcellindex!=-2){
								mrdcell* dtneighbourcell = thecells.at(upcell->dtneighbourcellindex);
								dtneighbourcell->parentcellindex = upcell->cellid;
								dtneighbourcell->utneighbourcellindex = -1;
								dtneighbourcell->neighbourchi2 = 0.;
								upcell->dtneighbourcellindex=-2;
							}
							downcell->parentcellindex = upcell->cellid;
							upcell->hasdaughters = true;
						} else {
#ifdef TRACKFINDVERBOSE
							cout<<"no existing downtrack candidate..."<<endl;
#endif
							// this cell does not have a downstream cell yet. If the downstream cell has no
							// existing upstream cell, use this one. If the downstream cell already has an
							// upstream neighbour defined (converging upstream tracks), use the better fit
							if( (downcell->utneighbourcellindex==-1) || 
								(downcell->utneighbourcellindex!=-1 && chi2<downcell->neighbourchi2) ){
#ifdef TRACKFINDVERBOSE
								cout<<"considering cell "<<acelli<<" the upstream neighbour of cell "<<bcelli<<endl;
#endif
								downcell->utneighbourcellindex = upcell->cellid;
								downcell->neighbourchi2 = chi2;
								upcell->dtneighbourcellindex = downcell->cellid;
							} // else this cell already has a neighbour with a better fit. Nothing to do.
						}
						// we need to check if the track splits, so continue scanning downstream cells
					} else if((!(upcell->isdownstreamgoing))&&(!(downcell->isdownstreamgoing))){
						// both are components of upstream going tracks. Check for any existing daughters
						// of this parent
						if(downcell->dtneighbourcellindex!=-1){
#ifdef TRACKFINDVERBOSE
							cout<<"this is a second downtrack candidate! Track splits!"<<endl;
#endif
							// this cell has already been given a down-track cell
							if(downcell->dtneighbourcellindex!=-2){
								mrdcell* dtneighbourcell = thecells.at(downcell->dtneighbourcellindex);
								dtneighbourcell->parentcellindex = downcell->cellid;
								dtneighbourcell->utneighbourcellindex = -1;
								dtneighbourcell->neighbourchi2 = 0.;
								downcell->dtneighbourcellindex=-2;
							}
							upcell->parentcellindex = downcell->cellid;
							downcell->hasdaughters = true;
						} else {
							// this cell has no existing downtrack cells, so proceed to make neighbours
							// (first checking if we have an upstream neighbour - choose best)
							if( (upcell->utneighbourcellindex==-1) || 
							(upcell->utneighbourcellindex!=-1 && chi2<upcell->neighbourchi2) ){
								// this cell has no upstream candidate, or this one has a better alignment
#ifdef TRACKFINDVERBOSE
								cout<<"considering cell "<<bcelli<<" the upstream neighbour of cell "<<acelli<<endl;
#endif
								upcell->utneighbourcellindex = downcell->cellid;
								upcell->neighbourchi2 = chi2;
								downcell->dtneighbourcellindex = upcell->cellid;
							} // else this cell already has a better fitting uptrack neighbour. Nothing to do.
						}
					}	// end if based on upstream/downstream
				}	// end if based on chi2 limit.
			}	// end if based on shared cluster
		}	// end loop over downstream cells
#ifdef TRACKFINDVERBOSE
			if(upcell->dtneighbourcellindex==-1) cout<<"no downstream neighbour cell found"<<endl;
			else cout<<"downstream neighbour cell found"<<endl;
#endif
	}	// end loop over upstream cells
	
	int cellswithoutupstreamneighbours=0;
	int cellswithupstreamneighbours=0;
	for(int acelli=0; acelli<thecells.size(); acelli++){
		mrdcell* thiscell = thecells.at(acelli);
		if(thiscell->utneighbourcellindex==-1) { cellswithoutupstreamneighbours++; }
		else { cellswithupstreamneighbours++; }
	}
#ifdef TRACKFINDVERBOSE
	cout<<"found "<<cellswithoutupstreamneighbours<<" cells without any upstream neighbours"
		<<" and "<<cellswithupstreamneighbours<<" cells with upstream neighbours"<<endl;
#endif
	
	// 7. now the automation [ do {...} while (...) ]
	//  a. for each iteration, loop over all cells
	//   > if any upstream neighbour cells have the same state as this cell, increment this cell's state number
	// repeat while there at least 1 cell has a neighbour with the same state number.
	if(cellswithupstreamneighbours){
#ifdef TRACKFINDVERBOSE
		cout<<"doing automation!"<<endl;
#endif
		int loopiterations=0;
		Bool_t madeachange;
		std::vector <Bool_t> newcellstatuses(thecells.size());
		do {
	#ifdef TRACKFINDVERBOSE
			cout<<"automation loop "<<loopiterations<<endl;
	#endif
			madeachange=false;
			for(int acelli=0; acelli<thecells.size(); acelli++){
				mrdcell* thiscell = thecells.at(acelli);
				if(thiscell->utneighbourcellindex!=-1){
					// the cell has a neighbour: compare statuses
					mrdcell* neighbourcell = thecells.at(thiscell->utneighbourcellindex);
					if(thiscell->status==neighbourcell->status){
						// the cell has the same status as it's upstream neighbour - increment it's status
						newcellstatuses.at(acelli) = true;
						madeachange=true;
#ifdef TRACKFINDVERBOSE
						cout<<"incrementing status of cell "<<acelli<<endl;
#endif
					} else {
						// the cell doesn't have the same status as upstream neighbour - keep same status
						newcellstatuses.at(acelli) = false;
					}
				} else {
					// no neighbour - don't increment any statuses
					newcellstatuses.at(acelli) = false;
				}
			}
			// now perform the actual status incrementing - we have to do this independantly of status checking
			for(int acelli=0; acelli<thecells.size(); acelli++){
				if(newcellstatuses.at(acelli)) thecells.at(acelli)->IncrementStatus();
			}
			loopiterations++;
		} while (madeachange&&loopiterations<30);
		assert(loopiterations!=29 && "automation loop failed to exit normally!");
		// should never have loopiterations > (max(numhlayers, numvlayers)-1) 
		//  - this would be a track longer than the MRD!
	}
	
	// 8. cells with a state of >= mintracklength and no (upstream) neighbours are defined as 
	// the downstream track ends.
	std::vector<Int_t> trackstartindices; 
	for(int acelli=0; acelli<thecells.size(); acelli++){
		mrdcell* thiscell = thecells.at(acelli);
		// check if the cell is a down-track endpoint - i.e. if it has downtrack neighbours
		if(thiscell->dtneighbourcellindex>-1) continue;
		// filter short tracks here by also not making tracks if the downtrack endpoint has a low status
		// value - i.e. few upstream neighbours. We allow short tracks if they instead have or are a parent. 
		if((thiscell->status<mintracklength)&&(thiscell->parentcellindex==-1)&&(thiscell->hasdaughters==false)){
#ifdef TRACKFINDVERBOSE
			cout<<"skipping short track with endpoint cell "<<acelli<<endl;
#endif
			continue;
		}
		trackstartindices.push_back(acelli);	// this is a valid downstream end cell!
	}
#ifdef TRACKFINDVERBOSE
	cout<<"found "<<trackstartindices.size()<<" track startpoints"<<endl;
#endif
	
	// 9. These cells are combined with upstream neighbours until the cell with state 0 (upstream edge). 
	//    This combination is the track. 
	std::vector<std::vector<mrdcell*> > alltracks;
	for(int atracki=0; atracki<trackstartindices.size(); atracki++){
		std::vector<mrdcell*> thistrack;
		mrdcell* thiscell = thecells.at(trackstartindices.at(atracki));
		thistrack.push_back(thiscell);
#ifdef TRACKFINDVERBOSE
		cout<<"track "<<atracki<<" has start cell index "<<trackstartindices.at(atracki)<<endl;
#endif
		do {
			int nextupcellindex = thiscell->utneighbourcellindex;
#ifdef TRACKFINDVERBOSE
			cout<<"cell "<<thiscell<<" has neighbour index "<<nextupcellindex<<endl;
#endif
			if(nextupcellindex==-1) {
				if(thiscell->status!=0) cerr<<"cell with nonzero status but no upstream neighbour!"<<endl;
				break;
			}
#ifdef TRACKFINDVERBOSE
			cout<<"adding cell "<<nextupcellindex<<" to track "<<atracki<<endl;
#endif
			thiscell = thecells.at(nextupcellindex);
			thistrack.push_back(thiscell);
		} while (1);
#ifdef TRACKFINDVERBOSE
		cout<<"track "<<atracki<<" has a total of "<<thistrack.size()<<" cells"<<endl;
#endif
		if(thistrack.size()) alltracks.push_back(thistrack);
	}
	
	// Make the actual cMRDTrack track objects, now that we have a coherent set of cells:
	for(int tracki=0; tracki<alltracks.size(); tracki++){
		// get the track (vector of cells)
		std::vector<mrdcell*> thistrack = alltracks.at(tracki);
		
		// declare all the things to make a track class
		std::vector<Int_t> digitids_inthistrack;
		std::vector<Int_t> tubeids_inthistrack;
		std::vector<Double_t> digitqs_inthistrack;
		std::vector<Double_t> digittimes_inthistrack;
		std::vector<Int_t> digitnumphots_inthistrack;
		std::vector<Double_t> digittruetimes_inthistrack;
		std::vector<Int_t> digittrueparents_inthistrack;
		
		// scan through all the cells, and grab the digit ids of all digits in all clusters in this track
		for(int celli=0; celli<thistrack.size(); celli++){
			mrdcell* acell = thistrack.at(celli);
			mrdcluster* upcluster = acell->clusters.first;
			mrdcluster* downcluster = acell->clusters.second;
			if(celli==0){
				// since clusters are shared between cells, take the upcluster only from the first cell
				digitids_inthistrack.insert(digitids_inthistrack.end(), upcluster->digitindexes.begin(), upcluster->digitindexes.end());
			}
			digitids_inthistrack.insert(digitids_inthistrack.end(), downcluster->digitindexes.begin(), downcluster->digitindexes.end());
		}
		
		for(int digiti=0; digiti<digitids_inthistrack.size(); digiti++){
			Int_t thedigitindex = digitids_inthistrack.at(digiti);
			tubeids_inthistrack.push_back(pmts_hit.at(thedigitindex));
			digitqs_inthistrack.push_back(digi_qs.at(thedigitindex));
			digittimes_inthistrack.push_back(digi_ts.at(thedigitindex));
			digitnumphots_inthistrack.push_back(digi_numphots.at(thedigitindex));
			digittruetimes_inthistrack.push_back(digi_phot_ts.at(thedigitindex));
			digittrueparents_inthistrack.push_back(digi_phot_parents.at(thedigitindex));
		}
		
		cMRDTrack* newtrack = new cMRDTrack(MRDSubEventID+tracki, wcsimfile, run_id, event_id, subtrigger, digitids_inthistrack, tubeids_inthistrack, digitqs_inthistrack, digittimes_inthistrack, digitnumphots_inthistrack, digittruetimes_inthistrack, digittrueparents_inthistrack);
		tracksthissubevent.push_back(newtrack);
//		tracksthissubevent.emplace_back(MRDSubEventID+tracki, wcsimfile, run_id, event_id, subtrigger, digitids_inthistrack, tubeids_inthistrack, digitqs_inthistrack, digittimes_inthistrack, digitnumphots_inthistrack, digittruetimes_inthistrack, digittrueparents_inthistrack);
	}
	
	Bool_t printtracks=true;
	if(printtracks){
	cout<<"Found "<<alltracks.size()<<" tracks this subevent."<<endl;
		for(int i=0; i<alltracks.size(); i++){
			std::vector<mrdcell*> thistrack = alltracks.at(i);
			cout<<"track "<<i<<" has "<<thistrack.size()<<" cells:"<<endl;
			for(int j=(thistrack.size()-1); j>-1; j--){
				cout<<endl;
				mrdcell* acell = thistrack.at(j);
				mrdcluster* upcluster = acell->clusters.first;
				mrdcluster* downcluster = acell->clusters.second;
				if(!(acell->isdownstreamgoing)){
					// switch the order of clusters, since the track is upstream going,
					// but cluster ordering is downstream going
					mrdcluster* temp = upcluster;
					upcluster=downcluster;
					downcluster=temp;
				}
				cout<<"   cell "<<j<<" goes from cluster "<<upcluster->clusterid<<" in layer "<<upcluster->layer
					<<" to cluster "<<downcluster->clusterid<<" in layer "<<downcluster->layer<<endl;
			
				cout<<"      uptrack cluster contains tubes ";
				for(int k=0; k<upcluster->in_layer_tubeids.size(); k++){
					(k!=0) ? cout<<", " : cout<<" ";
					cout<<upcluster->in_layer_tubeids.at(k);
				}
				cout<<endl;
				cout<<"      maximum in-layer index is "<<upcluster->xmaxid
					<<", minimum in-layer index is "<<upcluster->xminid
					<<", effective centre index is "<<upcluster->GetCentreIndex()<<endl;
				cout<<"      physical extent is "<<upcluster->GetXmin()<<" to "<<upcluster->GetXmax()<<endl;
			
				cout<<endl;
			
				cout<<"      downtrack cluster contains tubes ";
				for(int k=0; k<downcluster->in_layer_tubeids.size(); k++){
					(k!=0) ? cout<<", " : cout<<" ";
					cout<<downcluster->in_layer_tubeids.at(k);
				}
				cout<<endl;
				cout<<"      maximum in-layer index is "<<downcluster->xmaxid
					<<", minimum in-layer index is "<<downcluster->xminid
					<<", effective centre index is "<<downcluster->GetCentreIndex()<<endl;
				cout<<"      physical extent is "<<downcluster->GetXmin()<<" to "<<downcluster->GetXmax()<<endl;
			}
		}
	}
	
	// if we want to draw the cells, we need to do it here - cells are not stored in the track!
	Bool_t drawcells=true;
	if(drawcells){
		for(int i=0; i<alltracks.size(); i++){
			std::vector<mrdcell*> thistrack = alltracks.at(i);
			for(int j=0; j<thistrack.size(); j++){
				mrdcell* acell = thistrack.at(j);
				mrdcluster* upcluster = acell->clusters.first;
				mrdcluster* downcluster = acell->clusters.second;
				// cluster ordering is downstream going, regardless of track direction
				// a cluster may have multiple tubes hit; find the top and bottom tubes
				Int_t uptubetopid = (2*upcluster->xmaxid) + layeroffsets.at(upcluster->layer);
				Int_t uptubebottomid = (2*upcluster->xminid) + layeroffsets.at(upcluster->layer);
				// find the positions in the canvas of the those tubes' bottom corners
				TBox* uptopbox = cMRDSubEvent::paddlepointers.at(uptubetopid);
				TBox* upbottombox = cMRDSubEvent::paddlepointers.at(uptubebottomid);
				if(uptopbox==0||upbottombox==0){
					cout<<"null box! uptubetopid="<<uptubetopid<<",uptubebottomid="<<uptubebottomid<<endl;
				}
				// find the average x and y positions: this is the centre position on the canvas of this cluster
				Double_t upclustery = (upbottombox->GetY1()+uptopbox->GetY2())/2.;
				Double_t upclusterx = upbottombox->GetX2()+0.005;
				// repeat process with downcluster to find arrow endpoint in canvas
				Int_t downtubetopid = (2*downcluster->xmaxid) + layeroffsets.at(downcluster->layer);
				Int_t downtubebottomid = (2*downcluster->xminid) + layeroffsets.at(downcluster->layer);
				TBox* downtopbox = paddlepointers.at(downtubetopid);
				TBox* downbottombox = paddlepointers.at(downtubebottomid);
				Double_t downclustery = (downbottombox->GetY1()+downtopbox->GetY2())/2.;
				Double_t downclusterx = downbottombox->GetX2()+0.005;
				
				// choose arrow head based on track direction
				TArrow* myarrow;
				if(acell->isdownstreamgoing){
					// downcluster is the arrow head
					myarrow = new TArrow(upclusterx, upclustery, downclusterx, downclustery, 0.01, ">");
				} else {
					myarrow = new TArrow(upclusterx, upclustery, downclusterx, downclustery, 0.01, "<");
				}
				myarrow->SetLineWidth(2);
				(downcluster->layer%2==0) ? imgcanvas->cd(1) : imgcanvas->cd(2);
				myarrow->Draw();
				cout<<"drawing arrow from "<<myarrow->GetX1()<<", "<<myarrow->GetY1()<<" to "
					<<myarrow->GetX2()<<", "<<myarrow->GetY2()<<endl;
				trackarrows.push_back(myarrow);
			}
		}
	}
	
	
	// CLEANUP: cMRDTrack objects have been created: we no longer need the clusters or cells objects.
	for(auto acluster : theclusters){
		delete acluster;
	}
	theclusters.clear();
	for(auto acell : thecells){
		delete acell;
	}
	thecells.clear();
}

void cMRDSubEvent::LeastSquaresMinimizer(Int_t numdatapoints, Double_t datapointxs[], Double_t datapointys[], Double_t datapointweights[], Double_t errorys[], Double_t &fit_gradient, Double_t &fit_offset, Double_t &chi2){
	// linear least squares finder from $ROOTSYS/tutorials/matrix/solveLinear.C
	const Int_t nrVar  = 2;				// number of...? fit variables? 2 because a linear fit?
	Int_t method = 2;					// implementation method to use
	
//	// now passed as args
//	//const Int_t numdatapoints = 3;
//	Double_t datapointxs[]      = {0.0,1.0,2.0};
//	Double_t datapointys[]      = {1.4,1.5,3.7};
//	Double_t datapointweights[] = {0.5,0.2,1.0};
	
	// Make the vectors 'Use" the data : they are not copied, the vector data
	// pointer is just set appropriately
	TVectorD x; x.Use(numdatapoints,datapointxs);
	TVectorD y; y.Use(numdatapoints,datapointys);
	TVectorD e;
	if(datapointweights!=0) e.Use(numdatapoints,datapointweights);
	
	TMatrixD A(numdatapoints,nrVar);
	TMatrixDColumn(A,0) = 1.0;
	TMatrixDColumn(A,1) = x;
	
	// Two possible methods: first the 'Normal Equations' whose derivation can be found
	// here: http://mathworld.wolfram.com/LeastSquaresFitting.html
	// for root documentation of NormalEqn see https://root.cern/doc/master/TDecompChol_8cxx.html
	if(method==1){
		// Method 1. solve through 'Normal Equations'
		if(datapointweights==0){
#ifdef TRACKFINDVERBOSE
			cout<<"solving with method 1, no weightings"<<endl;
#endif
			const TVectorD c_norm = NormalEqn(A,y);
			fit_offset=c_norm(0);
			fit_gradient=c_norm(1);
	//		return;
		} else {
#ifdef TRACKFINDVERBOSE
			cout<<"solving with method 2, with weightings"<<endl;
#endif
			// weightings are passed by providing the 'e' vector to NormalEqn.
			const TVectorD c_norm = NormalEqn(A,y,e);
			fit_offset=c_norm(0);
			fit_gradient=c_norm(1);
	//		return;
		}
	}
	
	if(method==2){
		// Method 2. solve through single-value-decomposition (SVD) (numerically  preferred method)
		if(datapointweights==0){
#ifdef TRACKFINDVERBOSE
			cout<<"solving with method 2, no weightings"<<endl;
#endif
			// no weightings, proceed with unweighted datapoints
			TDecompSVD svd(A);
			Bool_t ok;
			const TVectorD c_svd = svd.Solve(y,ok);
#ifdef TRACKFINDVERBOSE
			cout<<"TDecompSVD returned ok="<<ok<<endl;
#endif
			fit_offset=c_svd(0);
			fit_gradient=c_svd(1);
	//		return;
		} else {
#ifdef TRACKFINDVERBOSE
			cout<<"solving with method 2, with weightings"<<endl;
#endif
			// weightings - fill weighted matrix/vectors
			TMatrixD Aw = A;
			TVectorD yw = y;
			for (Int_t irow = 0; irow < A.GetNrows(); irow++) {
			TMatrixDRow(Aw,irow) *= 1/e(irow);
			yw(irow) /= e(irow);
			}
			// solve with the weighted datapoints
			TDecompSVD svd(Aw);
			Bool_t ok;
			const TVectorD c_svd = svd.Solve(yw,ok);
#ifdef TRACKFINDVERBOSE
			cout<<"TDecompSVD returned ok="<<ok<<endl;
#endif
			fit_offset=c_svd(0);
			fit_gradient=c_svd(1);
	//		return;
		}
	}
	
	chi2=0.;
	for(int i=0; i<numdatapoints; i++){
		Double_t fitpred = (fit_gradient*datapointxs[i]) + fit_offset;
#ifdef TRACKFINDVERBOSE
		cout<<"predicted datapoint "<<i<<" = "<<fitpred
			<<", actual datapoint = "<<datapointys[i]<<endl;
#endif
		// calculate the error on the position from the paddle extents. Since the cluster may span
		// several paddles, we have to calculate this on a per-cluster basis
		chi2 += TMath::Power((datapointys[i]-fitpred),2)/errorys[i];
	}
	
	return;
}
