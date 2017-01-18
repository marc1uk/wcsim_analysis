{

TChain* c = new TChain("tankflux");
c->Add("/home/marc/LinuxSystemFiles/WChSandBox/build/fluxesandtables/annie_tank_flux.20*");
Int_t numents = c->GetEntries();
TH1I* neutrons = new TH1I("neutrons","Number of Neutrons",20,0,20);
TH1I* protons = new TH1I("protons","Number of Protons",20,0,20);
TH1I* nucleons = new TH1I("nucleons","Number of Nucleons",20,0,20);
Int_t ntankbranchval=0;
TBranch *nTankBranch=0;
Int_t* pdgbranchval=0;
TBranch* pdgBranch=0;
Int_t* nuprimarybranchval=0;
TBranch* nuprimaryBranch=0;
c->LoadTree(0);
Int_t treeNumber = -1;
for(Int_t inputEntry=0; inputEntry<numents; inputEntry++){
	Long64_t localEntry = c->LoadTree(inputEntry);
	Int_t nextTreeNumber = c->GetTreeNumber();
	if(treeNumber!=nextTreeNumber){
		c->SetBranchAddress("ntank",&ntankbranchval,&nTankBranch);
		pdgBranch=c->GetBranch("pdgtank");
		nuprimaryBranch=c->GetBranch("primary");
		treeNumber=nextTreeNumber;
	}
	nTankBranch->GetEntry(localEntry);
	if(pdgbranchval){delete[] pdgbranchval;}
	if(nuprimarybranchval){delete[] nuprimarybranchval;}
	pdgbranchval = new Int_t[ntankbranchval];
	nuprimarybranchval = new Int_t[ntankbranchval];
	pdgBranch->SetAddress(pdgbranchval);
	nuprimaryBranch->SetAddress(nuprimarybranchval);
	pdgBranch->GetEntry(localEntry);
	nuprimaryBranch->GetEntry(localEntry);
	Int_t protonsthisevent=0;
	Int_t neutronsthisevent=0;
	for(int i=0;i<ntankbranchval;i++){
		if(pdgbranchval[i]==2212){ protonsthisevent++; }
		else if(pdgbranchval[i]==2112){ neutronsthisevent++; }
	}
	protons->Fill(protonsthisevent);
	neutrons->Fill(neutronsthisevent);
	nucleons->Fill(protonsthisevent+neutronsthisevent);
}
protons->SetTitle("Number of Protons and Neutrons from Genie v2.8");
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
leg->Draw();


}
