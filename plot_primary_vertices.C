{
//TFile* f=TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/ANNIEtest_10MeV_Uni_Iso_0.root");
//TFile* f=TFile::Open("/home/marc/LinuxSystemFiles/Bonsai/createlike_nickwp/analysis_wcsim/ANNIEtest_10MeV_Uni_Iso_0.root");
TFile* f=TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in/wcsim_10MeV_iso_e_wDN_13-04-17.root");
TTree* t = (TTree*)f->Get("wcsimT");
WCSimRootEvent* ev;
TBranch* evb =0;
t->SetBranchAddress("wcsimrootevent",&ev,&evb);
WCSimRootTrigger* tr;
std::vector<TVector3> allvs;
for(int i=0; i<t->GetEntries(); i++){
    evb->GetEntry(i);
    tr=ev->GetTrigger(0);
    TVector3 avec(tr->GetVtx(0),tr->GetVtx(1),tr->GetVtx(2));
    allvs.push_back(avec);
}
TH3D ahist = TH3D("ahist2","My Hist",100,-300,300,100,-400,400,100,-100,500);
TH1D xplot = TH1D("xplot","My XPlot",100,-400,400);
TH1D yplot = TH1D("yplot","My YPlot",100,-300,300);
TH1D zplot = TH1D("zplot","My ZPlot",100,-100,500);
for(int i=0;i<allvs.size();i++){
  TVector3 avec=allvs.at(i); ahist.Fill(avec.X(), avec.Y(), avec.Z());
  xplot.Fill(avec.X()); yplot.Fill(avec.Y()); zplot.Fill(avec.Z());
}
TCanvas c1("c1","Ac1",600,400);
ahist.Draw();
TCanvas c2("c2","Ac2",600,400);
xplot.Draw();
TCanvas c3("c3","Ac3",600,400);
yplot.Draw();
TCanvas c4("c4","Ac4",600,400);
zplot.Draw();
}
