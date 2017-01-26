// #######################################################################

// Draw Canvases with Histograms for All Detector Elements
// =======================================================
void WCSimAnalysis::DrawGlobalHistos(){
	win_scale=0.7;
	n_wide=3;
	n_high=3;
	
	// histograms combining data about all detector elements
	hitCountCanv = new TCanvas("hitCountCanv","Title",700*n_wide*win_scale,500*n_high*win_scale);
	hitCountCanv->Divide(3,3);

	//tank plots
	hitCountCanv->cd(1);
	tubeshithisttank->Draw();
	hitCountCanv->cd(2);
	truehitcounthisttank->Draw();
	//hitCountCanv->cd(3);
	digihitcounthisttank->Draw("same");
	hitCountCanv->cd(3);
	PMThitfrequencytank->Draw();

	//mrd plots
	hitCountCanv->cd(4);
	tubeshithistmrd->Draw();
	hitCountCanv->cd(5);
	truehitcounthistmrd->Draw();
	//hitCountCanv->cd(6);
	digihitcounthistmrd->Draw("same");
	hitCountCanv->cd(6);
	PMThitfrequencymrd->Draw();

	//veto plots
	hitCountCanv->cd(7);
	tubeshithistveto->Draw();
	hitCountCanv->cd(8);
	truehitcounthistveto->Draw();
	//hitCountCanv->cd(9);
	digihitcounthistveto->Draw("same");
	hitCountCanv->cd(9);
	PMThitfrequencyveto->Draw();
}
