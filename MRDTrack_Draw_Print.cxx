/* vim:set noexpandtab tabstop=4 wrap */

void cMRDTrack::Print(){
	for(int j=(htrackcells.size()-1); j>-1; j--){
		cout<<endl;
		mrdcell* acell = &htrackcells.at(j);
		acell->SetClusterAddresses(htrackclusters);
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
	// ------------
	cout<<"vertical track "<<MRDtrackID<<" at "<<&vtrackcells[0]<<" has "
		<<vtrackcells.size()<<" cells:"<<endl;
	for(int j=(vtrackcells.size()-1); j>-1; j--){
		cout<<endl;
		mrdcell* acell = &vtrackcells.at(j);
		cout<<"setting cluster addresses for cell "<<j<<endl;
		acell->SetClusterAddresses(vtrackclusters);
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

void cMRDTrack::Draw(TCanvas* imgcanvas, std::vector<TArrow*> &trackarrows, EColor thistrackscolour, std::vector<TBox*> paddlepointers){
	// draw horizontal paddle track
	imgcanvas->cd(1);
	for(int j=0; j<htrackcells.size(); j++){
		mrdcell* acell = &htrackcells.at(j);
		acell->SetClusterAddresses(htrackclusters);
		mrdcluster* upcluster = acell->clusters.first;
		mrdcluster* downcluster = acell->clusters.second;
		// cluster ordering is downstream going, regardless of track direction
		// a cluster may have multiple tubes hit; find the top and bottom tubes
		Int_t uptubetopid = (2*upcluster->xmaxid) + layeroffsets.at(upcluster->layer);
		Int_t uptubebottomid = (2*upcluster->xminid) + layeroffsets.at(upcluster->layer);
		if((upcluster->altviewhalf)<0){ uptubetopid++; uptubebottomid++; }
		// find the positions in the canvas of the those tubes' bottom corners
		TBox* uptopbox = paddlepointers.at(uptubetopid);
		TBox* upbottombox = paddlepointers.at(uptubebottomid);
		if(uptopbox==0||upbottombox==0){
			cerr<<"null box! uptubetopid="<<uptubetopid<<",uptubebottomid="<<uptubebottomid<<endl;
		}
		// find the average x and y positions: this is the centre position on the canvas of this cluster
		Double_t upclustery = (upbottombox->GetY1()+uptopbox->GetY2())/2.;
		Double_t upclusterx = upbottombox->GetX2();
		// repeat process with downcluster to find arrow endpoint in canvas
		Int_t downtubetopid = (2*downcluster->xmaxid) + layeroffsets.at(downcluster->layer);
		Int_t downtubebottomid = (2*downcluster->xminid) + layeroffsets.at(downcluster->layer);
		if((downcluster->altviewhalf)<0){ downtubetopid++; downtubebottomid++; }
		TBox* downtopbox = paddlepointers.at(downtubetopid);
		TBox* downbottombox = paddlepointers.at(downtubebottomid);
		Double_t downclustery = (downbottombox->GetY1()+downtopbox->GetY2())/2.;
		Double_t downclusterx = downbottombox->GetX2();
		//if((upcluster->altviewhalf)==0){ upclusterx+=0.0025; }
		//if((downcluster->altviewhalf)==0){ downclusterx+=0.0025; }
		
		// choose arrow head based on track direction
		TArrow* myarrow;
		if(acell->isdownstreamgoing){
			// downcluster is the arrow head
			myarrow = new TArrow(upclusterx, upclustery, downclusterx, downclustery, 0.01, ">");
		} else {
			myarrow = new TArrow(upclusterx, upclustery, downclusterx, downclustery, 0.01, "<");
		}
		myarrow->SetLineWidth(2);
		myarrow->SetLineColor(thistrackscolour);
		//(downcluster->layer%2==0) ? imgcanvas->cd(1) : imgcanvas->cd(2);
		myarrow->Draw();
		imgcanvas->Update();
#ifdef MRDTrack_VERBOSE
		cout<<"drawing reconstructed track arrow from "<<myarrow->GetX1()<<", "<<myarrow->GetY1()<<" to "
			<<myarrow->GetX2()<<", "<<myarrow->GetY2()<<endl;
#endif
		//std::this_thread::sleep_for (std::chrono::seconds(3));
		trackarrows.push_back(myarrow);
	}
	// draw vertical paddle track
	imgcanvas->cd(2);
	for(int j=0; j<vtrackcells.size(); j++){
		mrdcell* acell = &vtrackcells.at(j);
		acell->SetClusterAddresses(vtrackclusters);
		mrdcluster* upcluster = acell->clusters.first;
		mrdcluster* downcluster = acell->clusters.second;
		// cluster ordering is downstream going, regardless of track direction
		// a cluster may have multiple tubes hit; find the top and bottom tubes
		Int_t uptubetopid = (2*upcluster->xmaxid) + layeroffsets.at(upcluster->layer);
		Int_t uptubebottomid = (2*upcluster->xminid) + layeroffsets.at(upcluster->layer);
		if((upcluster->altviewhalf)<0){ uptubetopid++; uptubebottomid++; }
		// find the positions in the canvas of the those tubes' bottom corners
		TBox* uptopbox = paddlepointers.at(uptubetopid);
		TBox* upbottombox = paddlepointers.at(uptubebottomid);
		if(uptopbox==0||upbottombox==0){
			cerr<<"null box! uptubetopid="<<uptubetopid<<",uptubebottomid="<<uptubebottomid<<endl;
		}
		// find the average x and y positions: this is the centre position on the canvas of this cluster
		Double_t upclustery = (upbottombox->GetY1()+uptopbox->GetY2())/2.;
		Double_t upclusterx = upbottombox->GetX2();
		// repeat process with downcluster to find arrow endpoint in canvas
		Int_t downtubetopid = (2*downcluster->xmaxid) + layeroffsets.at(downcluster->layer);
		Int_t downtubebottomid = (2*downcluster->xminid) + layeroffsets.at(downcluster->layer);
		if((downcluster->altviewhalf)<0){ downtubetopid++; downtubebottomid++; }
		TBox* downtopbox = paddlepointers.at(downtubetopid);
		TBox* downbottombox = paddlepointers.at(downtubebottomid);
		Double_t downclustery = (downbottombox->GetY1()+downtopbox->GetY2())/2.;
		Double_t downclusterx = downbottombox->GetX2();
		//if((upcluster->altviewhalf)==0){ upclusterx+=0.0025; }
		//if((downcluster->altviewhalf)==0){ downclusterx+=0.0025; }
	
		// choose arrow head based on track direction
		TArrow* myarrow;
		if(acell->isdownstreamgoing){
			// downcluster is the arrow head
			myarrow = new TArrow(upclusterx, upclustery, downclusterx, downclustery, 0.01, ">");
		} else {
			myarrow = new TArrow(upclusterx, upclustery, downclusterx, downclustery, 0.01, "<");
		}
		myarrow->SetLineWidth(2);
		myarrow->SetLineColor(thistrackscolour);
		//(downcluster->layer%2==0) ? imgcanvas->cd(1) : imgcanvas->cd(2);
		myarrow->Draw();
		imgcanvas->Update();
#ifdef MRDTrack_VERBOSE
		cout<<"drawing reconstructed track arrow from "<<myarrow->GetX1()<<", "<<myarrow->GetY1()<<" to "
			<<myarrow->GetX2()<<", "<<myarrow->GetY2()<<endl;
#endif
		//std::this_thread::sleep_for (std::chrono::seconds(3));
		trackarrows.push_back(myarrow);
	}
}
