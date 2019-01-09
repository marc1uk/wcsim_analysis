/* vim:set noexpandtab tabstop=4 wrap */
// #######################################################################

// DEFINING HISTOGRAMS
// ===================
void WCSimAnalysis::DefineVetoHistos(){
	tubeshithistveto = new TH1F("tubeshithistveto","Number of Veto PMTs Hit;Number of tubes;Frequency", numvetopmts, 1., numvetopmts);
	truehitcounthistveto = new TH1F("truehitcounthistveto","Number of True Veto PMT Hits per Event;Number of true hits;Frequency", numvetopmts, 1., numvetopmts);
	digihitcounthistveto = new TH1F("digihitcounthistveto","Number of Digitized Veto PMT Hits per Event;Number of digitized hits;Frequency", numvetopmts, 1., numvetopmts);
	digihitcounthistveto->SetLineColor(kRed);
	PEdistveto = new TH1D("PEdistveto","Veto Num PEs Distribution;Num PEs;Frequency", 16,-0.5,30.5);
	PMThitfrequencyveto = new TH1D("PMThitfrequencyveto","Num Hits vs Veto PMT;PMT ID;Num Hits", numvetopmts,0.5,numvetopmts+0.5);
	QvsTveto = new TH2D("QvsTveto","Veto charge vs. time;time;charge", 40, 900, 1400, 40, -0.5, 30.5);
	hitTimeDistveto = new TH1D("hitTimeDistveto","Veto Hit Time Distribution;Hit Time (ns);Frequency", 100,900,3400);
	digitTimeDistveto = new TH1D("digitTimeDistveto","Veto Digitized Hit Time Distribution;Digit Time (ns);Frequency", 100,900,3400);
//	PMTsvDigitTimeveto = new TH2D("PMTsvDigitTimeveto","Num Veto PMTs hit vs Last Digit Time;Digit Time (ns);Num PMTs Hit", 100,900,3400, numvetopmts, 0, numvetopmts);
	facchist = new TH2D("facchist","Map of FACC Hits", 30, 0, 29, 2, 0, 1);
	
	if(drawvetohistos){
		win_scale=0.7;
		n_wide=2;
		n_high=3;
	
		// charge distributions for veto
		QvsTvetoCanv = new TCanvas("QvsTvetoCanv","QvsTvetoCanv",700*n_wide*win_scale,500*n_high*win_scale);
		QvsTvetoCanv->Divide(n_wide,n_high);
		
//		n_wide=1;
//		n_high=1;
//		vetoLastTimevsNumHitsCanv = new TCanvas("vetoLastTimevsNumHitsCanv","Last Veto Digit vs Num Veto PMTs hit",700*n_wide*win_scale,500*n_high*win_scale);
	}
}


// #######################################################################

// FILLING EVENT-WIDE HISTOGRAMS
// =============================
void WCSimAnalysis::FillVetoEventWideHists(Int_t numtruehits, Int_t numdigits){
	truehitcounthistveto->Fill(numtruehits);
	digihitcounthistveto->Fill(numdigits);
	
	int numvetopmtshit = atrigv->GetNumTubesHit();
	tubeshithistveto->Fill(numvetopmtshit);
}

// #######################################################################

// FILLING TRUE HITS HISTOGRAMS
// =============================
void WCSimAnalysis::FillVetoTrueHitsHist(WCSimRootCherenkovHit* hit, WCSimRootCherenkovHitTime *hittime){
	PMThitfrequencyveto->Fill(hit->GetTubeID());
	PEdistveto->Fill(hit->GetTotalPe(1));	//(num 'pe's = number of true photon hits after QE)
	
	hitTimeDistveto->Fill(((hittime->GetTruetime())/500)+triggeroffset);
}

// #######################################################################

// FILLING DIGITS HISTOGRAMS
// =========================
void WCSimAnalysis::FillVetoDigiHitsHist(WCSimRootCherenkovDigiHit* digihit){
	QvsTveto->Fill(digihit->GetT(), digihit->GetQ());
	digitTimeDistveto->Fill(digihit->GetT());
	//if(DigiHit->GetT()>lastdigittime){ lastdigittime=digihit->GetT(); }
}

// #######################################################################

// DRAWING HISTOGRAMS
// ==================
void WCSimAnalysis::DrawVetoHistos(){
	if(drawvetohistos){
		QvsTvetoCanv->cd(1);
		QvsTveto->Draw("colz");

		TH1 *temp;
		QvsTvetoCanv->cd(2);
		temp=QvsTveto->ProjectionY();
		temp->SetTitle("charge");
		temp->Draw();
		QvsTvetoCanv->GetPad(2)->SetLogy();

		QvsTvetoCanv->cd(3);
		temp=QvsTveto->ProjectionX();
		temp->SetTitle("hits vs time");
		temp->Draw();
		QvsTvetoCanv->GetPad(3)->SetLogy();

		QvsTvetoCanv->cd(4);
		temp=QvsTveto->ProfileX();
		temp->SetTitle("average charge vs time");
		temp->Draw();

		QvsTvetoCanv->cd(5);
		temp=PEdistveto;
		temp->Draw();
		QvsTvetoCanv->GetPad(5)->SetLogy();

		QvsTvetoCanv->cd(6);
		temp=hitTimeDistveto;
		temp->Draw();
		temp=digitTimeDistveto;
		temp->SetLineColor(kRed);
		temp->Draw("same");
		QvsTvetoCanv->GetPad(6)->SetLogy();

	//	vetoLastTimevsNumHitsCanv->cd();
	//	//temp=PMTsvDigitTimeveto->ProjectionX();
	//	//temp->Draw();
	//	//PMTsvDigitTimeveto->Draw();
	}
}
