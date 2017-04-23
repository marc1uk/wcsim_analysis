//TFile* f = TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_0.root");
TFile* f = TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in/wcsim_10MeV_iso_e_wDN_13-04-17.root");
TTree* t= (TTree*)f->Get("wcsimGeoT")
WCSimRootGeom* geo =0
t->SetBranchAddress("wcsimrootgeom",&geo)
t->Draw("wcsimrootgeom.fPMTArray[].fPosition[0]:wcsimrootgeom.fPMTArray[].fPosition[1]:wcsimrootgeom.fPMTArray[].fPosition[2]")
t->Draw("wcsimrootgeom.fLAPPDArray[].fPosition[0]:wcsimrootgeom.fLAPPDArray[].fPosition[1]:wcsimrootgeom.fLAPPDArray[].fPosition[2]","","same")
TCanvas c1
c1.cd()
t->Draw("wcsimrootgeom.fLAPPDArray[].fPosition[0]")

// plot the positions as a different Draw("same") with markercolor set by number of photons hit?
// or superk event display: plot markers in a "box" option - i.e. with the size of the box set by number of PMTs hit, and with the colouring set by the time of the digit.

