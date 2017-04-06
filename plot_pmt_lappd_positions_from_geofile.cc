TFile* f = TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_0.root");
TTree* t= (TTree*)f->Get("wcsimGeoT")
WCSimRootGeom* geo =0
t->SetBranchAddress("wcsimrootgeom",&geo)
t->Draw("wcsimrootgeom.fPMTArray[].fPosition[0]:wcsimrootgeom.fPMTArray[].fPosition[1]:wcsimrootgeom.fPMTArray[].fPosition[2]")
t->Draw("wcsimrootgeom.fLAPPDArray[].fPosition[0]:wcsimrootgeom.fLAPPDArray[].fPosition[1]:wcsimrootgeom.fLAPPDArray[].fPosition[2]","","same")


// plot the positions as a different Draw("same") with markercolor set by number of photons hit?
// or superk event display: plot markers in a "box" option - i.e. with the size of the box set by number of PMTs hit, and with the colouring set by the time of the digit.
