/* vim:set noexpandtab tabstop=4 wrap */
// #######################################################################

// DEFINING HISTOGRAMS
// ===================
void WCSimAnalysis::DefineMRDhistos(){
	tubeshithistmrd = new TH1F("tubeshithistmrd","Number of MRD PMTs Hit;Number of tubes;Frequency", 100, 1., nummrdpmts);
	truehitcounthistmrd = new TH1F("truehitcounthistmrd","Number of True MRD PMT Hits per Event;Number of true hits;Frequency", 100, 1., 250.);
	digihitcounthistmrd = new TH1F("digihitcounthistmrd","Number of Digitized MRD PMT Hits per Event;Number of digitized hits;Frequency", 100, 1., 250.);
	digihitcounthistmrd->SetLineColor(kRed);
	PEdistmrd = new TH1D("PEdistmrd","MRD Num PEs Distribution;Num PEs;Frequency", 16,-0.5,30.5);
	PMThitfrequencymrd = new TH1D("PMThitfrequencymrd","Num Hits vs MRD PMT;PMT ID;Num Hits", nummrdpmts,0.5,nummrdpmts+0.5);
	QvsTmrd = new TH2D("QvsTmrd","MRD charge vs. time;time;charge", 40, 900, 1400, 40, -0.5, 30.5);
	hitTimeDistmrd = new TH1D("hitTimeDistmrd","MRD Hit Time Distribution;Hit Time (ns);Frequency", 100,900,3400);
	digitTimeDistmrd = new TH1D("digitTimeDistmrd","MRD Digitized Hit Time Distribution;Digit Time (ns);Frequency", 100,900,3400);
	mrdhist = new TH2D("mrdhist","Map of MRD Hits", 15, 0, 14, 13, 0, 12);
	
	if(drawmrdhistos){
		win_scale=0.7;
		n_wide=2;
		n_high=3;
	
		QvsTmrdCanv = new TCanvas("QvsTmrdCanv","QvsTmrdCanv",700*n_wide*win_scale,500*n_high*win_scale);
		QvsTmrdCanv->Divide(n_wide,n_high);
	}
}

// #######################################################################

// FILLING EVENT-WIDE HISTOGRAMS
// =============================
void WCSimAnalysis::FillMRDeventWideHists(Int_t numtruehits, Int_t numdigits){
	truehitcounthistmrd->Fill(numtruehits);
	digihitcounthistmrd->Fill(numdigits);
	
	int nummrdpmtshit = atrigm->GetNumTubesHit();
	tubeshithistmrd->Fill(nummrdpmtshit);
}

// #######################################################################

// FILLING TRUE HITS HISTOGRAMS
// =============================
void WCSimAnalysis::FillMRDtrueHitsHist(WCSimRootCherenkovHit* hit, WCSimRootCherenkovHitTime *hittime){
	PMThitfrequencymrd->Fill(hit->GetTubeID());
	PEdistmrd->Fill(hit->GetTotalPe(1));	//(num 'pe's = number of true photon hits after QE)
	
	hitTimeDistmrd->Fill(((hittime->GetTruetime())/500)+triggeroffset);
}

// #######################################################################

// FILLING DIGIT HISTOGRAMS
// ========================
void WCSimAnalysis::FillMRDdigiHitsHist(WCSimRootCherenkovDigiHit* digihit){
	QvsTmrd->Fill(digihit->GetT(), digihit->GetQ());
	
	digitTimeDistmrd->Fill(digihit->GetT());
}

// #######################################################################

// DRAWING HISTOGRAMS
// ==================
void WCSimAnalysis::DrawMRDhistos(){
	if(drawmrdhistos){
		QvsTmrdCanv->cd(1);
		QvsTmrd->Draw("colz");

		TH1 *temp;
		QvsTmrdCanv->cd(2);
		temp=QvsTmrd->ProjectionY();
		temp->SetTitle("charge");
		temp->Draw();
		QvsTmrdCanv->GetPad(2)->SetLogy();

		QvsTmrdCanv->cd(3);
		temp=QvsTmrd->ProjectionX();
		temp->SetTitle("hits vs time");
		temp->Draw();
		QvsTmrdCanv->GetPad(3)->SetLogy();

		QvsTmrdCanv->cd(4);
		temp=QvsTmrd->ProfileX();
		temp->SetTitle("average charge vs time");
		temp->Draw();

		QvsTmrdCanv->cd(5);
		temp=PEdistmrd;
		temp->Draw();
		QvsTmrdCanv->GetPad(5)->SetLogy();

		QvsTmrdCanv->cd(6);
		temp=hitTimeDistmrd;
		temp->Draw();
		temp=digitTimeDistmrd;
		temp->SetLineColor(kRed);
		temp->Draw("same");
		QvsTmrdCanv->GetPad(6)->SetLogy();
	}
}



