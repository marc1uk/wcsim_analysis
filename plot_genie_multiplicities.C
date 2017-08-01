{

TChain* c = new TChain("tankflux");
TChain* cm = new TChain("tankmeta");
c->Add("/home/marc/LinuxSystemFiles/WChSandBox/build/fluxesandtables/annie_tank_flux.20*");
cm->Add("/home/marc/LinuxSystemFiles/WChSandBox/build/fluxesandtables/annie_tank_flux.20*");
Int_t numents = c->GetEntries();
TH1D* neutrons = new TH1D("neutrons","Number of Neutrons",20,0,20);
TH1D* protons = new TH1D("protons","Number of Protons",20,0,20);
TH1D* nucleons = new TH1D("nucleons","Number of Nucleons",20,0,20);
TH1D* primneutrons = new TH1D("primneutrons","Number of Primary Neutrons",20,0,20);
TH1D* primprotons = new TH1D("primprotons","Number of Primary Protons",20,0,20);
TH1D* primnucleons = new TH1D("primnucleons","Number of Primary Nucleons",20,0,20);

TH1D* nener = new TH1D("nener","nne",100,0,10);
Int_t ntankbranchval=0;
TBranch *nTankBranch=0;
Int_t* pdgbranchval=0;
Double_t* nenerval=0;
TBranch* pdgBranch=0;
Int_t* nuprimarybranchval=0;
TBranch* nuprimaryBranch=0;
TBranch* primene=0;
double totalplots=0;
double pots;
c->LoadTree(0);
cm->LoadTree(0);
TBranch* potsbranch;
Int_t treeNumber = -1;
Int_t numprimaryevents=0;
for(Int_t inputEntry=0; inputEntry<numents; inputEntry++){
    Long64_t localEntry = c->LoadTree(inputEntry);
    Int_t nextTreeNumber = c->GetTreeNumber();
    if(treeNumber!=nextTreeNumber){
        c->SetBranchAddress("ntank",&ntankbranchval,&nTankBranch);
        pdgBranch=c->GetBranch("pdgtank");
        nuprimaryBranch=c->GetBranch("primary");
        primene=c->GetBranch("kE");
        cm->SetBranchAddress("inputTotalPOTs", &pots, &potsbranch);
        potsbranch->GetEntry(0);
        totalplots+=pots;
        treeNumber=nextTreeNumber;
    }
    nTankBranch->GetEntry(localEntry);
    if(pdgbranchval){delete[] pdgbranchval;}
    if(nuprimarybranchval){delete[] nuprimarybranchval;}
    pdgbranchval = new Int_t[ntankbranchval];
    nuprimarybranchval = new Int_t[ntankbranchval];
    if(nenerval) delete[] nenerval;
    nenerval = new Double_t[ntankbranchval];
    primene->SetAddress(nenerval);
    pdgBranch->SetAddress(pdgbranchval);
    nuprimaryBranch->SetAddress(nuprimarybranchval);
    pdgBranch->GetEntry(localEntry);
    nuprimaryBranch->GetEntry(localEntry);
    primene->GetEntry(localEntry);
    Int_t protonsthisevent=0;
    Int_t neutronsthisevent=0;
    Int_t primprotonsinthisevent=0;
    Int_t primneutronsinthisevent=0;
    bool primaryevent=false;
    for(int i=0;i<ntankbranchval;i++){
        if(pdgbranchval[i]==2212){
            protonsthisevent++;
            if(nuprimarybranchval[i]==1){ primprotonsinthisevent++; primaryevent=true; }
        }
        else if(pdgbranchval[i]==2112){
            neutronsthisevent++;
            if(nuprimarybranchval[i]==1){ primneutronsinthisevent++; primaryevent=true; nener->Fill(nenerval[i]);}
        }
    }
    protons->Fill(protonsthisevent);
    neutrons->Fill(neutronsthisevent);
    nucleons->Fill(protonsthisevent+neutronsthisevent);
    primprotons->Fill(primprotonsinthisevent);
    primneutrons->Fill(primneutronsinthisevent);
    primnucleons->Fill(primprotonsinthisevent+primneutronsinthisevent);
    if(primaryevent) numprimaryevents++;
}
double normalisation=1.;
TCanvas c1;
nucleons->SetTitle("Number of Protons and Neutrons from Genie v2.8");
nucleons->GetXaxis()->SetTitle("multiplicity");
nucleons->GetYaxis()->SetTitle("num events");
nucleons->SetLineColor(kViolet);
nucleons->Scale(normalisation/nucleons->Integral());
nucleons->Draw();
protons->SetLineColor(kBlue);
protons->Scale(normalisation/protons->Integral());
protons->Draw("same");
neutrons->SetLineColor(kRed);
neutrons->Scale(normalisation/neutrons->Integral());
neutrons->Draw("same");
TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(neutrons,"neutrons","l");
leg->AddEntry(protons,"protons","l");
leg->AddEntry(nucleons,"nucleons","l");
leg->Draw();
TCanvas c2;
primnucleons->SetTitle("Number of Primary Protons and Neutrons from Genie v2.8");
primnucleons->GetXaxis()->SetTitle("multiplicity");
primnucleons->GetYaxis()->SetTitle("num events");
primnucleons->SetLineColor(kViolet);
primnucleons->Scale(normalisation/primnucleons->Integral());
primnucleons->Draw();
primprotons->SetLineColor(kBlue);
primprotons->Scale(normalisation/primprotons->Integral());
primprotons->Draw("same");
primneutrons->SetLineColor(kRed);
primneutrons->Scale(normalisation/primneutrons->Integral());
primneutrons->Draw("same");
leg->Draw();
nener->Draw();

cout<<"plots are for "<<numprimaryevents<<" primary events"<<endl;


}
