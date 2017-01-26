{

Bool_t high_verbosity=false;
gStyle->SetOptStat(111111);	// draw overflow (and underflow? bins)
// note SetOptStat(0); to remove the stats box

TChain* c = new TChain("tankflux");
cout<<"adding chain files"<<endl;
c->Add("/home/marc/LinuxSystemFiles/WChSandBox/build/fluxesandtables/annie_tank_flux.20*");
Int_t numents = c->GetEntries();
cout<<"loaded "<<numents<<" entries"<<endl;
if(high_verbosity) cout<<"making histograms"<<endl;
TH1D* primaryoffsets = new TH1D("primaryoffsets","t0 Offsets of Primaries",100,6670,6690);
TH1D* secondaryoffsets = new TH1D("secondaryoffsets","t0 Offsets of Secondaries",100,-1,2000000000000);
TH1D* primaryranges = new TH1D("primaryranges","Time span of Event Primaries",100,-1,1);
TH1D* secondaryranges = new TH1D("secondaryranges","Time span of Event Secondaries",100,-1,2000000000000);
TH1D* allranges = new TH1D("allranges","Time span of Event Primaries+Secondaries",100,-1,2000000000000);
Int_t ntankbranchval=0;
TBranch *nTankBranch=0;
Double_t* vtxtbranchval=0;
TBranch* vtxtBranch=0;
Int_t* nuprimarybranchval=0;
TBranch* nuprimaryBranch=0;
if(high_verbosity) cout<<"loading tree 0"<<endl;
Int_t treeNumber = -1;
cout<<"looping over entries"<<endl;
for(Int_t inputEntry=0; inputEntry<numents; inputEntry++){
	if(high_verbosity) cout<<"loading entry"<<endl;
	Long64_t localEntry = c->LoadTree(inputEntry);
	if(localEntry<0){ cout<<"end of file"<<endl; break; }
	Int_t nextTreeNumber = c->GetTreeNumber();
	if(treeNumber!=nextTreeNumber){
		if(high_verbosity) cout<<"loading new tree"<<endl;
		c->SetBranchAddress("ntank",&ntankbranchval,&nTankBranch);
		vtxtBranch=c->GetBranch("vt");
		nuprimaryBranch=c->GetBranch("primary");
		if(nuprimaryBranch==0||vtxtBranch==0){cout<<"ZOMBIE BRANCHES"<<endl; exit(-1);}
		treeNumber=nextTreeNumber;
	}
	if(high_verbosity) cout<<"getting num primaries in this entry"<<endl;
	nTankBranch->GetEntry(localEntry);
	if(high_verbosity) cout<<"making arrays"<<endl;
	if(vtxtbranchval){delete[] vtxtbranchval;}
	if(nuprimarybranchval){delete[] nuprimarybranchval;}
	vtxtbranchval = new Double_t[ntankbranchval];
	nuprimarybranchval = new Int_t[ntankbranchval];
	if(high_verbosity) cout<<"setting branch addresses"<<endl;
	if(high_verbosity) cout<<"setting vtxtBranch ("<<vtxtBranch<<") address to "<<vtxtbranchval<<endl;
	vtxtBranch->SetAddress(vtxtbranchval);
	if(high_verbosity) cout<<"setting nuprimaryBranch ("<<nuprimaryBranch<<") address to "<<nuprimarybranchval<<endl;
	nuprimaryBranch->SetAddress(nuprimarybranchval);
	if(high_verbosity) cout<<"getting arrays"<<endl;
	vtxtBranch->GetEntry(localEntry);
	nuprimaryBranch->GetEntry(localEntry);
	if(high_verbosity) cout<<"looping over particles in this entry"<<endl;
	Double_t mintimethisevent=999.;
	Double_t maxtimethisevent=999.;
	Double_t mintimethiseventprims=999.;
	Double_t maxtimethiseventprims=999.;
	Double_t mintimethiseventsecs=999.;
	Double_t maxtimethiseventsecs=999.;
	for(int i=0;i<ntankbranchval;i++){
		if(nuprimarybranchval[i]==1){
			primaryoffsets->Fill(vtxtbranchval[i]);
			if(vtxtbranchval[i]<mintimethiseventprims||mintimethiseventprims==999.){
				mintimethiseventprims=vtxtbranchval[i];
			}
			if(vtxtbranchval[i]>maxtimethiseventprims||maxtimethiseventprims==999.){
				maxtimethiseventprims=vtxtbranchval[i];
			}
		}
		else {
			secondaryoffsets->Fill(vtxtbranchval[i]);
			if(vtxtbranchval[i]<mintimethiseventsecs||mintimethiseventsecs==999.){
				mintimethiseventsecs=vtxtbranchval[i];
			}
			if(vtxtbranchval[i]>maxtimethiseventsecs||maxtimethiseventsecs==999.){
				maxtimethiseventsecs=vtxtbranchval[i];
			}
		}
		if(vtxtbranchval[i]<mintimethisevent||mintimethisevent==999.){
			mintimethisevent=vtxtbranchval[i]; 
		}
		if(vtxtbranchval[i]>maxtimethisevent||maxtimethisevent==999.){
			maxtimethisevent=vtxtbranchval[i]; 
		}
		if(maxtimethisevent>100000000){ TTree* temp = c->GetTree(); temp->Show(localEntry); }
	}
	if(high_verbosity) cout<<"filling histograms"<<endl;
	allranges->Fill(maxtimethisevent-mintimethisevent);
	primaryranges->Fill(maxtimethiseventprims-mintimethiseventprims);
	secondaryranges->Fill(maxtimethiseventsecs-mintimethiseventsecs);
}
cout<<"drawing histograms"<<endl;
/*protons->SetTitle("Number of Protons and Neutrons from Genie v2.8");
protons->GetXaxis()->SetTitle("multiplicity");
protons->GetYaxis()->SetTitle("num events");
protons->Draw();
neutrons->SetLineColor(kRed);
neutrons->Draw("same");
nucleons->SetLineColor(kViolet);
nucleons->Draw("same");
TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(neutrons,"neutrons","l");
leg->AddEntry(protons,"protons","l");
leg->AddEntry(nucleons,"nucleons","l");
leg->Draw();*/

/*
TCanvas* c1 = new TCanvas("c1","c1");
c1->cd();
c1->SetLogy();
primaryoffsets->Draw();
TCanvas* c2 = new TCanvas("c2","c2");
c2->cd();
c2->SetLogy();
secondaryoffsets->Draw();
TCanvas* c3 = new TCanvas("c3","c3");
c3->cd();
c3->SetLogy();
primaryranges->Draw();
TCanvas* c4 = new TCanvas("c4","c4");
c4->cd();
c4->SetLogy();
secondaryranges->Draw();
*/

//TCanvas* c5 = new TCanvas("c5","c5");
//c5->cd();
//c5->SetLogy();
//allranges->Draw();


}
