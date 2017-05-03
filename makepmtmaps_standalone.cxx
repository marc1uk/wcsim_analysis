/* vim:set noexpandtab tabstop=4 wrap */
//############################################################################################

// PRODUCE MAP OF PMT POSITIONS
// ============================
// these maps, of wall, topcap and bottomcap, are used for 2D histograms for SK-style hit mapping
void MakePMTmap(WCSimRootGeom* geo, std::map<int, std::pair<int,int> > &topcappositionmap, std::map<int, std::pair<int,int> > &bottomcappositionmap, std::map<int, std::pair<int,int> > &wallpositionmap){
	// get the extent of z,x for the caps and x,y,z for the walls. 
	// Use the array of x,y,z values to produce a mapping of tubeID vs x-y for an unrolled view of the tank.
	std::vector<double> topcapxvalsall, topcapzvalsall, bottomcapxvalsall, bottomcapzvalsall, wallyvalsall, wallthetavalsall;
	std::vector<double> topcapxvals, topcapzvals, bottomcapxvals, bottomcapzvals, wallyvals, wallthetavals;
	std::vector<int> topcaptubeids, bottomcaptubeids, walltubeids;

	cout<<"Producing map of tank PMT positions"<<endl;
	cout<<geo->GetWCNumPMT()<<" total PMTs in geofile"<<endl;
	for(int i=0; i<geo->GetWCNumPMT();i++){
		WCSimRootPMT p = geo->GetPMT(i);
		double thexval = p.GetPosition(0);
		double theyval = p.GetPosition(1);
		double thezval = p.GetPosition(2);
		// we have some spurious entries? their positions are effectively 0,0,0.. 
		if(TMath::Abs(thexval)<1.||TMath::Abs(theyval)<1.||TMath::Abs(thezval)<1.){ cout<<"skipping spurious entry"<<endl; continue; }
		switch(p.GetCylLoc()){
			case 0: { // top cap pmt
				topcaptubeids.push_back(p.GetTubeNo());
				// make a vector of all the pmt positions
				topcapxvalsall.push_back(thexval);
				topcapzvalsall.push_back(thezval);
				// add the position to the vectors if it's a new one
				std::vector<double>::iterator it = std::find(topcapxvals.begin(), topcapxvals.end(), thexval);
				if(it==topcapxvals.end()){ topcapxvals.push_back(thexval); }
				it = std::find(topcapzvals.begin(), topcapzvals.end(), thezval);
				if(it==topcapzvals.end()){ topcapzvals.push_back(thezval); }
				break;
			}
			case 1: { // wall pmt
				walltubeids.push_back(p.GetTubeNo());
				// make a vector of all the pmt positions
				thezval -= geo->GetWCOffset(2); // z offset of the tank origin
				double thethetaval = TMath::ATan(thexval/TMath::Abs(thezval));
				if(thezval<0.){ (thexval<0.) ? thethetaval=(-TMath::Pi()+thethetaval) : thethetaval=(TMath::Pi()-thethetaval); }
				wallthetavalsall.push_back(thethetaval);
				wallyvalsall.push_back(theyval);
				// add the position to the vectors if it's a new one
				std::vector<double>::iterator it = std::find(wallyvals.begin(), wallyvals.end(), theyval);
				if(it==wallyvals.end()){ wallyvals.push_back(theyval); }
				it = std::find(wallthetavals.begin(), wallthetavals.end(), thethetaval);
				if(it==wallthetavals.end()){ wallthetavals.push_back(thethetaval); }
				break;
			}
			case 2: { // bottom cap pmt
				bottomcaptubeids.push_back(p.GetTubeNo());
				// make a vector of all the pmt positions
				bottomcapxvalsall.push_back(thexval);
				bottomcapzvalsall.push_back(thezval);
				// add the position to the vectors if it's a new one
				std::vector<double>::iterator it = std::find(bottomcapxvals.begin(), bottomcapxvals.end(), thexval);
				if(it==bottomcapxvals.end()){ bottomcapxvals.push_back(thexval); }
				it = std::find(bottomcapzvals.begin(), bottomcapzvals.end(), thezval);
				if(it==bottomcapzvals.end()){ bottomcapzvals.push_back(thezval); }
				break;
			}
		}
	}
	cout<<"    wall PMTs have "<<wallyvals.size()<<" unique y values and "
			<<wallthetavals.size()<<" unique theta vals"<<endl
			<<"    top cap PMTs have "<<topcapxvals.size()<<" unique x values and "
			<<topcapzvals.size()<<" unique z values"<<endl
			<<"    bottom cap PMTs have "<<bottomcapxvals.size()<<" unique x values and "
			<<bottomcapzvals.size()<<" unique z values"<<endl;
	// sort the unique positions into ascending order
	std::sort(topcapxvals.begin(),topcapxvals.end());
	std::sort(topcapzvals.begin(),topcapzvals.end());
	std::sort(bottomcapxvals.begin(),bottomcapxvals.end());
	std::sort(bottomcapzvals.begin(),bottomcapzvals.end());
	std::sort(wallyvals.begin(),wallyvals.end());
	std::sort(wallthetavals.begin(),wallthetavals.end());

	std::ofstream mapfile;
	mapfile.open("out/mapfile.txt", std::ios::out);
	//mapfile<<"wall map ordering:\n"<<endl;
	//for(int i=0; i<wallthetavals.size(); i++){
	//	mapfile<<wallthetavals.at(i)<<endl;
	//}
	//mapfile<<endl;
	// create a map (for each pmt, via it's tubeid) of pairs (for the 2 positions), but with positions now ints representing bin
	//std::map<int, std::pair<int,int> > topcappositionmap;	 in class def
	for(int i=0; i<topcaptubeids.size();i++){
		double thexval=topcapxvalsall.at(i);
		double thezval=topcapzvalsall.at(i);
		std::vector<double>::iterator it = std::find(topcapxvals.begin(), topcapxvals.end(), thexval);
		int xbin = std::distance(topcapxvals.begin(),it);
		it = std::find(topcapzvals.begin(), topcapzvals.end(), thezval);
		int zbin = std::distance(topcapzvals.begin(),it);
		topcappositionmap.emplace(topcaptubeids.at(i), std::pair<int,int>(xbin,zbin));
		mapfile << "topcapmap entry "<<i<<" tubeID: "<<topcaptubeids.at(i)<<", x: "<<xbin<<", z: "<<zbin<<endl;
	}
	cout<<"    top cap map has "<<topcappositionmap.size()<<" entries"<<endl;
	// repeat for bottom cap
	//std::map<int, std::pair<int,int> > bottomcappositionmap;
	for(int i=0; i<bottomcaptubeids.size();i++){
		double thexval=bottomcapxvalsall.at(i);
		double thezval=bottomcapzvalsall.at(i);
		std::vector<double>::iterator it = std::find(bottomcapxvals.begin(), bottomcapxvals.end(), thexval);
		int xbin = std::distance(bottomcapxvals.begin(),it);
		it = std::find(bottomcapzvals.begin(), bottomcapzvals.end(), thezval);
		int zbin = std::distance(bottomcapzvals.begin(),it);
		bottomcappositionmap.emplace(bottomcaptubeids.at(i), std::pair<int,int>(xbin,zbin));
		mapfile << "bottomcapmap entry "<<i<<" tubeID: "<<bottomcaptubeids.at(i)<<", x: "<<xbin<<", z: "<<zbin<<endl;
	}
	cout<<"    bottom cap map has "<<bottomcappositionmap.size()<<" entries"<<endl;
	// once more for the walls
	//std::map<int, std::pair<int,int> > wallpositionmap;
	for(int i=0; i<walltubeids.size();i++){
		double theyval=wallyvalsall.at(i);
		double thethetaval=wallthetavalsall.at(i);
		std::vector<double>::iterator it = std::find(wallthetavals.begin(), wallthetavals.end(), thethetaval);
		int thetabin = std::distance(wallthetavals.begin(),it);
		it = std::find(wallyvals.begin(), wallyvals.end(), theyval);
		int ybin = std::distance(wallyvals.begin(),it);
		wallpositionmap.emplace(walltubeids.at(i), std::pair<int,int>(thetabin,ybin));
		mapfile << "wallmap entry "<<i<<" tubeID: "<<walltubeids.at(i)<<", theta: "<<thetabin<<", y: "<<ybin<<endl;
	}
	cout<<"    wall map has "<<wallpositionmap.size()<<" entries"<<endl;
	mapfile.close();
	// incidentally, topcapxvalsall.size() should be the same as caparraysize
	//f->cd();
}
