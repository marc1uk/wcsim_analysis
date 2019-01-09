/* vim:set noexpandtab tabstop=4 wrap */
// #######################################################################

// DEFINING HISTOGRAMS
// ===================
void WCSimAnalysis::DefineTankHistos(){
	// Define Histograms
	// =================
	tubeshithisttank = new TH1F("tubeshithist","Number of Tank PMTs Hit;Number of tubes;Frequency", 100, 1., numtankpmts);
	truehitcounthisttank = new TH1F("truehitcounthisttank","Number of True Tank PMT Hits per Event;Number of true hits;Frequency", 100, 1., 200.);
	digihitcounthisttank = new TH1F("digihitcounthisttank","Number of Digitized Tank PMT Hits per Event;Number of digitized hits;Frequency", 100, 1., 200.);
	digihitcounthisttank->SetLineColor(kRed);
	PEdisttank = new TH1D("PEdisttank","Tank Num PEs Distribution;Num PEs;Frequency", 16,-0.5,30.5);
	PMThitfrequencytank = new TH1D("PMThitfrequencytank","Num Hits vs Tank PMT;PMT ID;Num Hits", numtankpmts,0.5,numtankpmts+0.5);
	QvsTtank = new TH2D("QvsTtank","Tank charge vs. time;time;charge", 40, 900, 1400, 40, -0.5, 30.5);
	topcaphist = new TH2D("topcaphist","Map of Hit Frequency for PMTs on the Top Cap",caparraysize+2,-1,caparraysize+1,caparraysize+2,-1,caparraysize+1);
	bottomcaphist = new TH2D("bottomcaphist","Map of Hit Frequency for PMTs on the Bottom Cap",caparraysize+2,-1,caparraysize+1,caparraysize+2,-1,caparraysize+1);
	wallhist = new TH2D("wallhist","Map of Hit Frequency for PMTs on the Tank Wall",pmtsperring+2,-1,pmtsperring+1,numpmtrings+2,-1,numpmtrings+1);

	// Define Canvases
	// ===============
	if(drawtankhistos){
		win_scale=1.5;
		n_wide=1;
		n_high=1;
		wallmapcanv = new TCanvas("wallmapcanv","Map of Wall Hits",700*n_wide*win_scale,500*n_high*win_scale);
		topcapmapcanv = new TCanvas("topcapmapcanv","Map of Top Cap Hits",700*n_wide*win_scale,500*n_high*win_scale);
		bottomcapmapcanv = new TCanvas("bottomcapmapcanv","Map of Bottom Cap Hits",700*n_wide*win_scale,500*n_high*win_scale);
	
		win_scale=0.7;
		n_wide=2;
		n_high=3;
		QvsTtankCanv = new TCanvas("QvsTtankCanv","QvsTtankCanv",700*n_wide*win_scale,500*n_high*win_scale);
		QvsTtankCanv->Divide(n_wide,n_high);
	}
}

// #######################################################################

// FILLING EVENT-WIDE HISTOGRAMS
// =============================
void WCSimAnalysis::FillTankEventWideHists(Int_t numtruehits, Int_t numdigits){
	truehitcounthisttank->Fill(numtruehits);
	digihitcounthisttank->Fill(numdigits);
	
	int numtankpmtshit = atrigt->GetNumTubesHit();
	tubeshithisttank->Fill(numtankpmtshit);
}

// #######################################################################

// FILLING TRUE HITS HISTOGRAMS
// =============================
void WCSimAnalysis::FillTankTrueHitsHist(WCSimRootCherenkovHit* hit){
	// ADD TO HISTOGRAM OF PMT HIT FREQUENCY
	PMThitfrequencytank->Fill(hit->GetTubeID());
	
	// ADD TO THE HISTOGRAM OF NUMBER OF PE'S FOR A TANK HIT
	PEdisttank->Fill(hit->GetTotalPe(1));	//(num 'pe's = number of true photon hits after QE)
}

// #######################################################################

// FILLING DIGITS HISTOGRAMS
// =========================
void WCSimAnalysis::FillTankDigiHitsHist(WCSimRootCherenkovDigiHit* digihit){
	QvsTtank->Fill(digihit->GetT(), digihit->GetQ());
	
	//let's do this better with a 2D map with positions
	int tubeID = digihit->GetTubeId();
	WCSimRootPMT pmt = geo->GetPMT(tubeID - 1);	// TUBE ID NEEDS -1 IN GEO FILE
	// WCSimRootPMT has members GetTubeNo(), GetCylLoc(), GetPosition(j), GetOrientation(j)
	// GetCylLoc(): 0=top cap, 2=bottom cap, 1=wall, 4=mrd, 5=veto, 3=obselete outer veto (shouldnt come up)
	// GetPosition(j), j=0..2: Returns x,y,z coordinates of the center of the sphere that forms the PMT.
	// GetOrientation(j), j=0..2: Returns the x,y,z components of a vector of the direction the PMT faces.
	// GetPMT(j) Returns a pmt object - NOT a pointer to a PMT object.
	//cout<<"Filling histogram for cylloc "<<pmt.GetCylLoc()<<" for tubeID "<<tubeID<<endl;
	switch(pmt.GetCylLoc()){
		//TODO: weighting by time doesn't make sense unless you can assure each PMT only has one hit
		case 0: {
			if(topcappositionmap.count(tubeID)){
				std::pair<int,int> thebins = topcappositionmap.at(tubeID);
				topcaphist->Fill(thebins.first, thebins.second, digihit->GetT());
				// do we weight by time? Do this in post-digit loop to find avg T?
			} else {cout<<"bad pmt: ID "<<tubeID<<" in CylLoc "<<pmt.GetCylLoc()<<endl;}
			break;
		}
		case 1: {
			if(wallpositionmap.count(tubeID)){
				std::pair<int,int> thebins = wallpositionmap.at(tubeID);
				wallhist->Fill(thebins.first+0.5, thebins.second, digihit->GetT());
			} else {cout<<"bad pmt: ID "<<tubeID<<" in CylLoc "<<pmt.GetCylLoc()<<endl;}
			break;
		}
		case 2: {
			if(bottomcappositionmap.count(tubeID)){
				std::pair<int,int> thebins = bottomcappositionmap.at(tubeID);
				bottomcaphist->Fill(thebins.first, thebins.second, digihit->GetT());
			} else {cout<<"bad pmt: ID "<<tubeID<<" in CylLoc "<<pmt.GetCylLoc()<<endl;}
			break;
		}
		case 4: {
//				std::pair<int,int> thebins = mrdpositionmap.at(tubeID);
//				mrdhist->Fill(thebins.first, thebins.second);
			break;
		}
		case 5: {
//				std::pair<int,int> thebins = faccpositionmap.at(tubeID);
//				facchist->Fill(thebins.first, thebins.second);
			break;
		}
		default: {
			//cout<<"PMT "<<tubeID<<" has unknown location "<<pmt.GetCylLoc()<<"!"<<endl; 
			break;
		}
	}
	
}

// #######################################################################

// DRAWING HISTOGRAMS
// ==================
void WCSimAnalysis::DrawTankHistos(){
	if(drawmrdhistos){
		wallmapcanv->cd();
		wallhist->Draw("colz");
		topcapmapcanv->cd();
		topcaphist->Draw("colz");
		bottomcapmapcanv->cd();
		bottomcaphist->Draw("colz");

		QvsTtankCanv->cd(1);
		QvsTtank->Draw("colz");
	
		TH1 *temp;
		QvsTtankCanv->cd(2);
		temp=QvsTtank->ProjectionY();
		temp->SetTitle("charge");
		temp->Draw();
		QvsTtankCanv->GetPad(2)->SetLogy();

		QvsTtankCanv->cd(3);
		temp=QvsTtank->ProjectionX();
		temp->SetTitle("hits vs time");
		temp->Draw();
		QvsTtankCanv->GetPad(3)->SetLogy();

		QvsTtankCanv->cd(4);
		temp=QvsTtank->ProfileX();
		temp->SetTitle("average charge vs time");
		temp->Draw();
		// fill in the lines... need to deal with memory allocation
		//AvgQvsTtankHist = new TH1F("AvgQvsTtankHist","average charge vs time",temp->GetNbinsX(),temp->GetXaxis()->GetBinLowEdge(0),temp->GetXaxis()->GetBinUpEdge(temp->GetNbinsX()));
		//for(int i=-1; i<(temp->GetNbinsX()+1); i++){AvgQvsTtankHist->SetBinContent(i,temp->GetBinContent(i)); }
		//AvgQvsTtankHist->SetTitle("average charge vs time");
		//AvgQvsTtankHist->SetLineColor(kRed);
		//AvgQvsTtankHist->Draw("same");

		QvsTtankCanv->cd(5);
		temp=PEdisttank;
		temp->Draw();
		QvsTtankCanv->GetPad(5)->SetLogy();

		//QvsTtankCanv->cd(6);
		//temp=PMThitfrequencytank;
		//temp->SetTitle("hit frequency distribution;pmt number;hit frequency");
		//temp->Draw();
	}
}
