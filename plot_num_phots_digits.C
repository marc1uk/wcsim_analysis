{
//TFile* f=TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/ANNIEtest_10MeV_Uni_Iso_0.root");
//TFile* f=TFile::Open("/home/marc/LinuxSystemFiles/Bonsai/createlike_nickwp/analysis_wcsim/ANNIEtest_10MeV_Uni_Iso_0.root");
//TFile* f=TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in/wcsim_tankonly_13-04-17.root");
TFile* f=TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in/wcsim_10MeV_iso_e_wDN_13-04-17.root");
TTree* t = (TTree*)f->Get("wcsimT");
WCSimRootEvent* ev;
TBranch* evb =0;
t->SetBranchAddress("wcsimrootevent",&ev,&evb);
WCSimRootTrigger* tr;
TH1D phhist=TH1D("phhist","Number Of Photon Hits per Event",100,0,6000);
TH1D dhhist=TH1D("dhhist","Number Of Digits per Event",100,0,200);
for(int i=0; i<t->GetEntries(); i++){evb->GetEntry(i); tr=ev->GetTrigger(0); Int_t numphots=tr->GetNcherenkovhittimes(); Int_t numdigs = tr->GetNcherenkovdigihits(); cout<<"event "<<i<<" had "<<numphots<<" photons and "<<numdigs<<" digits"<<endl; phhist.Fill(numphots); dhhist.Fill(numdigs);}
TCanvas c1("c1","canv1",600,800);
phhist.Draw();
TCanvas c2("c2","canv2",600,800);
dhhist.Draw();
}
