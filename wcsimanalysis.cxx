//#include <TROOT.h>
//#include <TSystem.h>
//#include <TApplication.h>
//#include "TFile.h"
//#include "TH1.h"
//#include "TH2.h"
//#include "TNtuple.h"
//#include "TChain.h"
//#include "TProfile.h"
//#include "TProfile2D.h"
//#include "TCanvas.h"
//#include "TString.h"
//#include <exception>	// for stdexcept
//#include <vector>
//#include <map>
//#include <string>
//#include <algorithm> // remove and remove_if

//void analyseout(){
{

Int_t numtankpmts=200;
Int_t nummrdpmts=336;
Int_t numvetopmts=26;
Int_t caparraysize=12;   // pmts on the cap form an nxn grid where caparraysize=n
Int_t pmtsperring=18;   // pmts around each ring of the main walls
Int_t numpmtrings=7;    // num rings around the main walls - 5x16=200, and 16%8=0 so 2 per octagonal face
//gROOT->Reset();
gROOT->ProcessLine("#include <regex>");
gSystem->Load("/home/marc/LinuxSystemFiles/WCSim/gitver/build/libWCSimRoot.so");	//also include TSystem.h
//R__LOAD_LIBRARY(libWCSimRoot.so);
cout<<"Loading files..."<<endl;
//TFile* f = TFile::Open("wcsim6.root");
//TTree* t = (TTree*)f->Get("wcsimT");

TChain* t = new TChain("wcsimT");
// loop over subdirectories and add all the files in each
const char* dirname = "/home/marc/anniegpvm/stats10k/";
const char* ext = ".root";
TSystemDirectory dir(dirname, dirname);
TList *subfolders = dir.GetListOfFiles(); // despite name returns files and folders
if(subfolders) {
  TSystemDirectory *subfolder;
  TIter nextsf(subfolders);
  while (subfolder=(TSystemDirectory*)nextsf()) {
    if (subfolder->IsDirectory()){
      TString sfname = subfolder->GetName();
      //cout<<"found subdirectory "<<sfname<<endl;
      if(sfname=="."||sfname==".."){ continue; }
/*    // don't actually need this information
      TList *files = subfolder->GetListOfFiles();
      TSystemFile *file;
      TIter nextf(files);
      while (file=(TSystemFile*)nextf()) {
        TString fname;
        fname = file->GetName();
        if(fname.EndsWith(ext)) { int i=0; }
      }
*/    
      TString nextfilepattern = dirname + sfname + "/wcsim_*";
      cout<<"adding "<<nextfilepattern<<endl;
      t->Add(nextfilepattern.Data());
      // cf. t->Add("/home/marc/anniegpvm/stats10k/stats_1_79283/wcsim_*");
    } // end if directory
  } // end loop over subdirectories
} // end if subfolders exist.
cout<<"Loaded "<<t->GetEntriesFast()<<" entries"<<endl;

// first let's just get the geotree now, since we only need to retrieve it once
// can be used to get PMTs by ID: note GetPMT(i) returns a PMT object not a pointer!
cout<<"Loading geometry...";
TFile* fgeo = TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_geom.root");
// currently using a separate file because the geometry tree was broken in the full stats version
TTree *geotree = (TTree*)fgeo->Get("wcsimGeoT");
WCSimRootGeom *geo = 0; 
geotree->SetBranchAddress("wcsimrootgeom", &geo);
if (geotree->GetEntries() == 0) { cout<<"geotree has no entries!"<<endl; exit(9); }
geotree->GetEntry(0);
cout<<"Geometry loaded"<<endl;
// get the extent of z,x for the caps and x,y,z for the walls. Use the array of x,y,z values to produce a mapping of tubeID vs x-y for an unrolled view of the tank.
std::vector<double> topcapxvalsall, topcapzvalsall, bottomcapxvalsall, bottomcapzvalsall, wallyvalsall, wallthetavalsall;
std::vector<double> topcapxvals, topcapzvals, bottomcapxvals, bottomcapzvals, wallyvals, wallthetavals;
std::vector<int> topcaptubeids, bottomcaptubeids, walltubeids;

cout<<"Producing map of tank PMT positions"<<endl;
for(int i=0; i<geo->GetWCNumPMT();i++){
  WCSimRootPMT p = geo->GetPMT(i);
  double thexval = p.GetPosition(0);
  double theyval = p.GetPosition(1);
  double thezval = p.GetPosition(2);
  // we have some spurious entries? their positions are effectively 0,0,0.. 
  if(TMath::Abs(thexval)<1.||TMath::Abs(theyval)<1.||TMath::Abs(thezval)<1.){ continue; }
  switch(p.GetCylLoc()){
    case 0: { // top cap pmt
      topcaptubeids.push_back(p.GetTubeNo());
      // make a vector of all the pmt positions
      topcapxvalsall.push_back(thexval);
      topcapzvalsall.push_back(thezval);
      // add the position to the vectors if it's a new one
      std::vector<double>::iterator it = std::find(topcapxvals.begin(), topcapxvals.end(), thexval);
      if(it==topcapxvals.end()){ topcapxvals.push_back(thexval); }
      it = std::find(topcapzvals.begin(), topcapzvals.end(), thezval);
      if(it==topcapzvals.end()){ topcapzvals.push_back(thezval); }
      break;
    }
    case 1: { // wall pmt
      walltubeids.push_back(p.GetTubeNo());
      // make a vector of all the pmt positions
      thezval -= geo->GetWCOffset(2); // z offset of the tank origin
      double thethetaval = TMath::ATan(thexval/thezval);
      if(thezval<0.){ (thexval<0.) ? thethetaval+=(2*TMath::Pi()) : thethetaval-=(2*TMath::Pi()); }
      wallthetavalsall.push_back(thethetaval);
      wallyvalsall.push_back(theyval);
      // add the position to the vectors if it's a new one
      std::vector<double>::iterator it = std::find(wallyvals.begin(), wallyvals.end(), theyval);
      if(it==wallyvals.end()){ wallyvals.push_back(theyval); }
      it = std::find(wallthetavals.begin(), wallthetavals.end(), thethetaval);
      if(it==wallthetavals.end()){ wallthetavals.push_back(thethetaval); break; }
    }
    case 2: { // bottom cap pmt
      bottomcaptubeids.push_back(p.GetTubeNo());
      // make a vector of all the pmt positions
      bottomcapxvalsall.push_back(thexval);
      bottomcapzvalsall.push_back(thezval);
      // add the position to the vectors if it's a new one
      std::vector<double>::iterator it = std::find(bottomcapxvals.begin(), bottomcapxvals.end(), thexval);
      if(it==bottomcapxvals.end()){ bottomcapxvals.push_back(thexval); }
      it = std::find(bottomcapzvals.begin(), bottomcapzvals.end(), thezval);
      if(it==bottomcapzvals.end()){ bottomcapzvals.push_back(thezval); }
      break;
    }
  }
}
cout<<"    wall PMTs have "<<wallyvals.size()<<" unique y values and "<<wallthetavals.size()<<" unique theta vals"<<endl
    <<"    top cap PMTs have "<<topcapxvals.size()<<" unique x values and "<<topcapzvals.size()<<" unique z values"<<endl
    <<"    bottom cap PMTs have "<<bottomcapxvals.size()<<" unique x values and "<<bottomcapzvals.size()<<" unique z values"<<endl;
// sort the unique positions into ascending order
std::sort(topcapxvals.begin(),topcapxvals.end());
std::sort(topcapzvals.begin(),topcapzvals.end());
std::sort(bottomcapxvals.begin(),bottomcapxvals.end());
std::sort(bottomcapzvals.begin(),bottomcapzvals.end());
std::sort(wallyvals.begin(),wallyvals.end());
std::sort(wallthetavals.begin(),wallthetavals.end());

std::ofstream mapfile;
mapfile.open("mapfile.txt", std::ios::out);
//mapfile<<"wall map ordering:\n"<<endl;
//for(int i=0; i<wallthetavals.size(); i++){
//	mapfile<<wallthetavals.at(i)<<endl;
//}
//mapfile<<endl;
// create a map (for each pmt, via it's tubeid) of pairs (for the 2 positions), but with positions now ints representing bin
std::map<int, std::pair<int,int> > topcappositionmap;
for(int i=0; i<topcaptubeids.size();i++){
  double thexval=topcapxvalsall.at(i);
  double thezval=topcapzvalsall.at(i);
  std::vector<double>::iterator it = std::find(topcapxvals.begin(), topcapxvals.end(), thexval);
  int xbin = std::distance(topcapxvals.begin(),it);
  it = std::find(topcapzvals.begin(), topcapzvals.end(), thezval);
  int zbin = std::distance(topcapzvals.begin(),it);
  topcappositionmap.emplace(topcaptubeids.at(i), std::pair<int,int>(xbin,zbin));
  mapfile << "topcapmap entry "<<i<<" tubeID: "<<topcaptubeids.at(i)<<", x: "<<xbin<<", z: "<<zbin<<endl;
}
cout<<"    top cap map has "<<topcappositionmap.size()<<" entries"<<endl;
// repeat for bottom cap
std::map<int, std::pair<int,int> > bottomcappositionmap;
for(int i=0; i<bottomcaptubeids.size();i++){
  double thexval=bottomcapxvalsall.at(i);
  double thezval=bottomcapzvalsall.at(i);
  std::vector<double>::iterator it = std::find(bottomcapxvals.begin(), bottomcapxvals.end(), thexval);
  int xbin = std::distance(bottomcapxvals.begin(),it);
  it = std::find(bottomcapzvals.begin(), bottomcapzvals.end(), thezval);
  int zbin = std::distance(bottomcapzvals.begin(),it);
  bottomcappositionmap.emplace(bottomcaptubeids.at(i), std::pair<int,int>(xbin,zbin));
  mapfile << "bottomcapmap entry "<<i<<" tubeID: "<<bottomcaptubeids.at(i)<<", x: "<<xbin<<", z: "<<zbin<<endl;
}
cout<<"    bottom cap map has "<<bottomcappositionmap.size()<<" entries"<<endl;
// once more for the walls
std::map<int, std::pair<int,int> > wallpositionmap;
for(int i=0; i<walltubeids.size();i++){
  double theyval=wallyvalsall.at(i);
  double thethetaval=wallthetavalsall.at(i);
  std::vector<double>::iterator it = std::find(wallthetavals.begin(), wallthetavals.end(), thethetaval);
  int thetabin = std::distance(wallthetavals.begin(),it);
  it = std::find(wallyvals.begin(), wallyvals.end(), theyval);
  int ybin = std::distance(wallyvals.begin(),it);
  wallpositionmap.emplace(walltubeids.at(i), std::pair<int,int>(thetabin,ybin));
  mapfile << "wallmap entry "<<i<<" tubeID: "<<walltubeids.at(i)<<", theta: "<<thetabin<<", y: "<<ybin<<endl;
}
cout<<"    wall map has "<<wallpositionmap.size()<<" entries"<<endl;
mapfile.close();
// incidentally, topcapxvalsall.size() should be the same as caparraysize
//f->cd();

cout<<"Retrieving event branches"<<endl;
WCSimRootEvent* b = new WCSimRootEvent();
WCSimRootEvent* m = new WCSimRootEvent();
WCSimRootEvent* v = new WCSimRootEvent();
TBranch* bp=0, *mp=0, *vp=0;
t->LoadTree(0);
t->SetBranchAddress("wcsimrootevent",&b, &bp);
t->SetBranchAddress("wcsimrootevent_mrd",&m, &mp);
t->SetBranchAddress("wcsimrootevent_facc",&v, &vp);
if(bp==0||mp==0||vp==0){ cout<<"first set of branches are zombies!"<<endl; 
	cout<<"bp = "<<bp<<", mp = "<<mp<<", vp = "<<vp<<endl; 
	return; 
}
t->GetEntry(0);
b->GetNumberOfEvents(); 	// number of sub events (delayed decays) within event 0.
WCSimRootTrigger* atrigt = b->GetTrigger(0);
WCSimRootTrigger* atrigm = m->GetTrigger(0);
WCSimRootTrigger* atrigv = v->GetTrigger(0);
Int_t eventnum=0;

Int_t eventswithtankhit=0, eventswithmrdhit=0, eventswithvetohit=0;
//=======
cout<<"Making histograms"<<endl;
//tank
TH1F* tubeshithisttank = new TH1F("tubeshithist","Number of Tank PMTs Hit;Number of tubes;Frequency", 100, 1., numtankpmts);
TH1F* truehitcounthisttank = new TH1F("truehitcounthisttank","Number of True Tank PMT Hits per Event;Number of true hits;Frequency", 100, 1., 200.);
TH1F* digihitcountshisttank = new TH1F("digihitcountshisttank","Number of Digitized Tank PMT Hits per Event;Number of digitized hits;Frequency", 100, 1., 200.);
digihitcountshisttank->SetLineColor(kRed);
TH1D* PEdisttank = new TH1D("PEdisttank","Tank Num PEs Distribution;Num PEs;Frequency", 16,-0.5,30.5);
TH1D* PMThitfrequencytank = new TH1D("PMThitfrequencytank","Num Hits vs Tank PMT;PMT ID;Num Hits", numtankpmts,0.5,numtankpmts+0.5);
TH2D* QvsTtank = new TH2D("QvsTtank","Tank charge vs. time;time;charge", 40, 900, 1400, 40, -0.5, 30.5);
//mrd
//these plots should be some kind of 2D histogram? makes more sense given pmt layout. layers?
TH1F* tubeshithistmrd = new TH1F("tubeshithistmrd","Number of MRD PMTs Hit;Number of tubes;Frequency", 100, 1., nummrdpmts);
TH1F* truehitcounthistmrd = new TH1F("truehitcounthistmrd","Number of True MRD PMT Hits per Event;Number of true hits;Frequency", 100, 1., 250.);
TH1F* digihitcountshistmrd = new TH1F("digihitcountshistmrd","Number of Digitized MRD PMT Hits per Event;Number of digitized hits;Frequency", 100, 1., 250.);
digihitcountshistmrd->SetLineColor(kRed);
TH1D* PEdistmrd = new TH1D("PEdistmrd","MRD Num PEs Distribution;Num PEs;Frequency", 16,-0.5,30.5);
TH1D* PMThitfrequencymrd = new TH1D("PMThitfrequencymrd","Num Hits vs MRD PMT;PMT ID;Num Hits", nummrdpmts,0.5,nummrdpmts+0.5);
TH2D* QvsTmrd = new TH2D("QvsTmrd","MRD charge vs. time;time;charge", 40, 900, 1400, 40, -0.5, 30.5);
TH1D* hitTimeDistmrd = new TH1D("hitTimeDistmrd","MRD Hit Time Distribution;Hit Time (ns);Frequency", 100,900,3400);
TH1D* digitTimeDistmrd = new TH1D("digitTimeDistmrd","MRD Digitized Hit Time Distribution;Digit Time (ns);Frequency", 100,900,3400);
//veto
TH1F* tubeshithistveto = new TH1F("tubeshithistveto","Number of Veto PMTs Hit;Number of tubes;Frequency", numvetopmts, 1., numvetopmts);
TH1F* truehitcounthistveto = new TH1F("truehitcounthistveto","Number of True Veto PMT Hits per Event;Number of true hits;Frequency", numvetopmts, 1., numvetopmts);
TH1F* digihitcountshistveto = new TH1F("digihitcountshistveto","Number of Digitized Veto PMT Hits per Event;Number of digitized hits;Frequency", numvetopmts, 1., numvetopmts);
digihitcountshistveto->SetLineColor(kRed);
TH1D* PEdistveto = new TH1D("PEdistveto","Veto Num PEs Distribution;Num PEs;Frequency", 16,-0.5,30.5);
TH1D* PMThitfrequencyveto = new TH1D("PMThitfrequencyveto","Num Hits vs Veto PMT;PMT ID;Num Hits", numvetopmts,0.5,numvetopmts+0.5);
TH2D* QvsTveto = new TH2D("QvsTveto","Veto charge vs. time;time;charge", 40, 900, 1400, 40, -0.5, 30.5);
TH1D* hitTimeDistveto = new TH1D("hitTimeDistveto","Veto Hit Time Distribution;Hit Time (ns);Frequency", 100,900,3400);
TH1D* digitTimeDistveto = new TH1D("digitTimeDistveto","Veto Digitized Hit Time Distribution;Digit Time (ns);Frequency", 100,900,3400);
TH2D* PMTsvDigitTimeveto = new TH2D("PMTsvDigitTimeveto","Num Veto PMTs hit vs Last Digit Time;Digit Time (ns);Num PMTs Hit", 100,900,3400, numvetopmts, 0, numvetopmts);

TH2D* topcaphist = new TH2D("topcaphist","Map of Hit Frequency for PMTs on the Top Cap",caparraysize,-1,caparraysize-1,caparraysize,-1,caparraysize-1);
TH2D* bottomcaphist = new TH2D("bottomcaphist","Map of Hit Frequency for PMTs on the Bottom Cap",caparraysize,-1,caparraysize-1,caparraysize,-1,caparraysize-1);
TH2D* wallhist = new TH2D("wallhist","Map of Hit Frequency for PMTs on the Tank Wall",pmtsperring,-1,pmtsperring-1,numpmtrings,-1,numpmtrings-1);
TH2D* mrdhist = new TH2D("mrdhist","Map of MRD Hits", 15, 0, 14, 13, 0, 12);
TH2D* facchist = new TH2D("facchist","Map of FACC Hits", 30, 0, 29, 2, 0, 1);
//=======

cout<<"Looping over entries"<<endl;
//for(eventnum=0; eventnum<bp->GetEntries();eventnum++){
int treeNumber=0;
do {
	Long64_t localEntry = t->LoadTree(eventnum);
	if (localEntry<0){ cout<<"breaking at entry "<<eventnum<<" as localEntry<0"<<endl; break; }
	Int_t nextTreeNumber = t->GetTreeNumber();
	if(treeNumber!=nextTreeNumber){
		cout<< "Reached end of Tree. Last entries' tree number was "
					<< treeNumber <<", this entries' tree number is "<< nextTreeNumber <<endl;
		cout<<"entries in this tree: "<<bp->GetEntries()<<endl; 
		// do need to re-set branch addresses
		t->SetBranchAddress("wcsimrootevent",&b, &bp);
		t->SetBranchAddress("wcsimrootevent_mrd",&m, &mp);
		t->SetBranchAddress("wcsimrootevent_facc",&v, &vp);
		if(bp==0||mp==0||vp==0){ cout<<"branches are zombies!"<<endl; break; }
		treeNumber=nextTreeNumber;
	}
//	G4cout<<"Loading data from entry "<<inputEntry<<", localentry "<<localEntry<<"/"<<entriesInThisTree<<G4endl;
	t->GetEntry(eventnum);
	atrigt = b->GetTrigger(0);
	atrigm = m->GetTrigger(0);
	atrigv = v->GetTrigger(0);
	
	//tank
	int numtankpmtshit = atrigt->GetNumTubesHit();
	tubeshithisttank->Fill(numtankpmtshit);
	int totnumtankhits = atrigt->GetCherenkovHits()->GetEntries();
	truehitcounthisttank->Fill(totnumtankhits);
	int totnumtankdigits = atrigt->GetCherenkovDigiHits()->GetEntries();
	digihitcountshisttank->Fill(totnumtankdigits);
	// photon hit information
	for (int i = 0; i<totnumtankhits; i++){
		WCSimRootCherenkovHit* hit = (WCSimRootCherenkovHit*)atrigt->GetCherenkovHits()->At(i);
		PMThitfrequencytank->Fill(hit->GetTubeID());
		//let's do this better with a 2D map with positions
		int tubeID = hit->GetTubeID();
		WCSimRootPMT pmt = geo->GetPMT(tubeID - 1);	// TUBE ID NEEDS -1 IN GEO FILE
		// WCSimRootPMT has members GetTubeNo(), GetCylLoc(), GetPosition(j), GetOrientation(j)
		// GetCylLoc(): 0=top cap, 2=bottom cap, 1=wall, 4=mrd, 5=veto, 3=obselete outer veto (shouldnt come up)
		// GetPosition(j), j=0..2: Returns x,y,z coordinates of the center of the sphere that forms the PMT.
		// GetOrientation(j), j=0..2: Returns the x,y,z components of a vector of the direction the PMT faces.
		// GetPMT(j) Returns a pmt object - NOT a pointer to a PMT object.
		//cout<<"Filling histogram for cylloc "<<pmt.GetCylLoc()<<" for tubeID "<<tubeID<<endl;
		switch(pmt.GetCylLoc()){
			case 0: {
				if(topcappositionmap.count(tubeID)){
					std::pair<int,int> thebins = topcappositionmap.at(tubeID);
					topcaphist->Fill(thebins.first, thebins.second);
				} else {cout<<"bad pmt: ID "<<tubeID<<" in CylLoc "<<pmt.GetCylLoc()<<endl;}
				break;
			}
			case 1: {
				if(wallpositionmap.count(tubeID)){
					std::pair<int,int> thebins = wallpositionmap.at(tubeID);
					wallhist->Fill(thebins.first, thebins.second);
				} else {cout<<"bad pmt: ID "<<tubeID<<" in CylLoc "<<pmt.GetCylLoc()<<endl;}
				break;
			}
			case 2: {
				if(bottomcappositionmap.count(tubeID)){
					std::pair<int,int> thebins = bottomcappositionmap.at(tubeID);
					bottomcaphist->Fill(thebins.first, thebins.second);
				} else {cout<<"bad pmt: ID "<<tubeID<<" in CylLoc "<<pmt.GetCylLoc()<<endl;}
				break;
			}
			case 4: {
//				std::pair<int,int> thebins = mrdpositionmap.at(tubeID);
//				mrdhist->Fill(thebins.first, thebins.second);
				break;
			}
			case 5: {
//				std::pair<int,int> thebins = faccpositionmap.at(tubeID);
//				facchist->Fill(thebins.first, thebins.second);
				break;
			}
			default: {
				//cout<<"PMT "<<tubeID<<" has unknown location "<<pmt.GetCylLoc()<<"!"<<endl; 
				break;
			}
		}
		//WCSimRootCherenkovHit has methods GetTubeId(), GetTotalPe(int)
		PEdisttank->Fill(hit->GetTotalPe(1));	//(num 'pe's = number of true photon hits after QE)
	}
	// digit information
	for (int i = 0; i<totnumtankdigits; i++){
		WCSimRootCherenkovDigiHit* DigiHit = (WCSimRootCherenkovDigiHit*)atrigt->GetCherenkovDigiHits()->At(i);
		//WCSimRootChernkovDigiHit has methods GetTubeId(), GetT(), GetQ()
		QvsTtank->Fill(DigiHit->GetT(), DigiHit->GetQ());
	}
	
	//mrd
	int nummrdpmtshit = atrigm->GetNumTubesHit();
	tubeshithistmrd->Fill(nummrdpmtshit);
	int totnummrdhits = atrigm->GetCherenkovHits()->GetEntries();
	truehitcounthistmrd->Fill(totnummrdhits);
	int totnummrddigits = atrigm->GetCherenkovDigiHits()->GetEntries();
	digihitcountshistmrd->Fill(totnummrddigits);
	// photon hit information
	for (int i = 0; i<totnummrdhits; i++){
		WCSimRootCherenkovHit* hit = (WCSimRootCherenkovHit*)atrigm->GetCherenkovHits()->At(i);
		PMThitfrequencymrd->Fill(hit->GetTubeID());
		//WCSimRootCherenkovHit has methods GetTubeId(), GetTotalPe(int)
		PEdistmrd->Fill(hit->GetTotalPe(1));	//(num 'pe's = number of true photon hits after QE)
		
		WCSimRootCherenkovHitTime *cHitTime = (WCSimRootCherenkovHitTime*)atrigm->GetCherenkovHitTimes()->At(i);
		//WCSimRootCherenkovHitTime has methods GetTubeId(), GetTruetime()
		hitTimeDistmrd->Fill(((cHitTime->GetTruetime())/500)+950);
	}
	// digit information
	for (int i = 0; i<totnummrddigits; i++){
		WCSimRootCherenkovDigiHit* DigiHit = (WCSimRootCherenkovDigiHit*)atrigm->GetCherenkovDigiHits()->At(i);
		//WCSimRootChernkovDigiHit has methods GetTubeId(), GetT(), GetQ()
		QvsTmrd->Fill(DigiHit->GetT(), DigiHit->GetQ());
		
		digitTimeDistmrd->Fill(DigiHit->GetT());
	}
	
	//veto 
	int numvetopmtshit = atrigv->GetNumTubesHit();
	tubeshithistveto->Fill(numvetopmtshit);
	int totnumvetohits = atrigv->GetCherenkovHits()->GetEntries();
	truehitcounthistveto->Fill(totnumvetohits);
	int totnumvetodigits = atrigv->GetCherenkovDigiHits()->GetEntries();
	digihitcountshistveto->Fill(totnumvetodigits);
	// photon hit information
	for (int i = 0; i<totnumvetohits; i++){
		WCSimRootCherenkovHit* hit = (WCSimRootCherenkovHit*)atrigv->GetCherenkovHits()->At(i);
		PMThitfrequencyveto->Fill(hit->GetTubeID());
		//WCSimRootCherenkovHit has methods GetTubeId(), GetTotalPe(int)
		PEdistveto->Fill(hit->GetTotalPe(1));	//(num 'pe's = number of true photon hits after QE)
		
		WCSimRootCherenkovHitTime *cHitTime = (WCSimRootCherenkovHitTime*)atrigv->GetCherenkovHitTimes()->At(i);
		//WCSimRootCherenkovHitTime has methods GetTubeId(), GetTruetime()
		hitTimeDistveto->Fill(((cHitTime->GetTruetime())/500)+950);
	}
	// digit information
	double lastdigittime=0;
	for (int i = 0; i<totnumvetodigits; i++){
		WCSimRootCherenkovDigiHit* DigiHit = (WCSimRootCherenkovDigiHit*)atrigv->GetCherenkovDigiHits()->At(i);
		//WCSimRootChernkovDigiHit has methods GetTubeId(), GetT(), GetQ()
		QvsTveto->Fill(DigiHit->GetT(), DigiHit->GetQ());
		digitTimeDistveto->Fill(DigiHit->GetT());
		lastdigittime=DigiHit->GetT();
	}
	PMTsvDigitTimeveto->Fill(lastdigittime,totnumvetodigits);
	
	eventnum++;
} while (1); //(!(t->LoadTree(eventnum)<0) );	// at start of loop. Need to remove loadtree at front if so. 

float win_scale=0.75;
int n_wide=3;
int n_high=3;

TCanvas* wallmapcanv = new TCanvas("wallmapcanv","Map of Wall Hits",700*n_wide*win_scale,500*n_high*win_scale);
wallmapcanv->cd();
wallhist->Draw("colz");
TCanvas* topcapmapcanv = new TCanvas("topcapmapcanv","Map of Top Cap Hits",700*n_wide*win_scale,500*n_high*win_scale);
topcapmapcanv->cd();
topcaphist->Draw("colz");
TCanvas* bottomcapmapcanv = new TCanvas("bottomcapmapcanv","Map of Wall Hits",700*n_wide*win_scale,500*n_high*win_scale);
bottomcapmapcanv->cd();
bottomcaphist->Draw("colz");

TCanvas* hitCountCanv = new TCanvas("vectCanv","Title",700*n_wide*win_scale,500*n_high*win_scale);
hitCountCanv->Divide(3,3);
//tank plots
hitCountCanv->cd(1);
tubeshithisttank->Draw();
hitCountCanv->cd(2);
truehitcounthisttank->Draw();
//hitCountCanv->cd(3);
digihitcountshisttank->Draw("same");
hitCountCanv->cd(3);
PMThitfrequencytank->Draw();
//mrd plots
hitCountCanv->cd(4);
tubeshithistmrd->Draw();
hitCountCanv->cd(5);
truehitcounthistmrd->Draw();
//hitCountCanv->cd(6);
digihitcountshistmrd->Draw("same");
hitCountCanv->cd(6);
PMThitfrequencymrd->Draw();
//veto plots
hitCountCanv->cd(7);
tubeshithistveto->Draw();
hitCountCanv->cd(8);
truehitcounthistveto->Draw();
//hitCountCanv->cd(9);
digihitcountshistveto->Draw("same");
hitCountCanv->cd(9);
PMThitfrequencyveto->Draw();

// charge distributions for tank =========================
TH1 *temp;
n_wide=2;
n_high=3;
TCanvas *QvsTtankCanv = new TCanvas("QvsTtankCanv","QvsTtankCanv",700*n_wide*win_scale,500*n_high*win_scale);
QvsTtankCanv->Divide(n_wide,n_high);
QvsTtankCanv->cd(1);
QvsTtank->Draw("colz");

QvsTtankCanv->cd(2);
temp=QvsTtank->ProjectionY();
temp->SetTitle("charge");
temp->Draw();
QvsTtankCanv->GetPad(2)->SetLogy();

QvsTtankCanv->cd(3);
temp=QvsTtank->ProjectionX();
temp->SetTitle("hits vs time");
temp->Draw();
QvsTtankCanv->GetPad(3)->SetLogy();

QvsTtankCanv->cd(4);
temp=QvsTtank->ProfileX();
TH1F* temp2 = new TH1F("temp2","average charge vs time",temp->GetNbinsX(),temp->GetXaxis()->GetBinLowEdge(0),temp->GetXaxis()->GetBinUpEdge(temp->GetNbinsX()));
for(int i=-1; i<(temp->GetNbinsX()+1); i++){temp2->SetBinContent(i,temp->GetBinContent(i)); }
temp->SetTitle("average charge vs time");
temp->Draw();

//temp=QvsTtank->ProfileX();
//temp->SetTitle("average charge vs time");
//temp->SetLineColor(kRed);
//temp->SetMarkerColor(kRed);
//temp->Draw("same");

QvsTtankCanv->cd(5);
temp=PEdisttank;
temp->Draw();
QvsTtankCanv->GetPad(5)->SetLogy();

//QvsTtankCanv->cd(6);
//temp=PMThitfrequencytank;
//temp->SetTitle("hit frequency distribution;pmt number;hit frequency");
//temp->Draw();

// charge distributions for mrd =========================
TCanvas *QvsTmrdCanv = new TCanvas("QvsTmrdCanv","QvsTmrdCanv",700*n_wide*win_scale,500*n_high*win_scale);
QvsTmrdCanv->Divide(n_wide,n_high);
QvsTmrdCanv->cd(1);
QvsTmrd->Draw("colz");

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

// charge distributions for veto =========================
TCanvas *QvsTvetoCanv = new TCanvas("QvsTvetoCanv","QvsTvetoCanv",700*n_wide*win_scale,500*n_high*win_scale);
QvsTvetoCanv->Divide(n_wide,n_high);
QvsTvetoCanv->cd(1);
QvsTveto->Draw("colz");

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

n_wide=1;
n_high=1;
TCanvas *vetodistcheckcanv = new TCanvas("vetodistcheckcanv","vetodistcheckcanv",700*n_wide*win_scale,500*n_high*win_scale);
vetodistcheckcanv->cd();
temp=PMTsvDigitTimeveto->ProjectionX();
temp->Draw();

//=======================================
// TODO: read genie info
/*
// do we need to do makeclass before hand so the compiler knows about the class? or not compile...
TFile* geniefile = TFile::Open(genieprimsfile);	//check if open, open if new
TTree* genietree = (TTree*)geniefile->Get("gst");	//check if we already have a genietree? or new tchain...whatv
// makes the files "gst.C" and "gst.h" which describe an appropriate class
genietree->MakeClass();	// check if we alread have gst.C, gst.h, or load them from somewhere
.L gst.C	//do we need enclosing "'s? not in interactive session, but maybe in built...
gst geniedat;	// maybe pointer better? who owns this? does it matter?
// now the loop
geniedat.GetEntry(entry);
std::string inttype = intxnumtotype(geniedat);
*/

// ok, let's look at the primary neutrino interaction and try'n get something useful:
/* we want to know:
GENERALLY
to be able to put numbers in context, read the number of POTs from all the files, scale by 
(number of events used / number of events in file), then / number of POTs per bnb dump (??)
-> this should give number of bnb dumps that correspond to the used events - /bnb dumps per day
to get number of days the analysed events correspond to.

NEUTRONS
========
total number of events with one or more final state neutrons. 
make a tree subset of events with these events in. 
number of final state neutrons vs various neutrino properties - i.e. MC neutron multiplicity prediction
e.g. against neutrino E, neutrino Q^2
for each event with neutrons, how many captured on Gd = true efficiency of Gd capture. 
(perhaps events with at least 1 Gd capture, or events where all FS neutrons Gd capture)
average amount of light (pe's, hits, digits, hitpmts) from a neutron capture (on Gd).
as a function of? position? anything else relevant?
time spread of hits from ncapture = characteristic time signal of Gd ncapture 
 - use 2d intensity distribution charge vs time! intensity is num events!
average spatial distribution?? ...dunno.

based on num hit pmts / pe's, need to make a 'selection' criteria for detecting neutron capture
efficiency of detection = num Gd capture events that pass this criteria vs total. 
efficiency vs various numbers of pmts? i.e. how calculate the appropriate selection criteria with hits/digits only on a particular subset of pmts?

MUONS
=====
make tree of subset of events with muons.
number of events with muons vs number of total events, just out of interest - should match the ratio of number of CCQE events vs non-CCQE. 

retrieve starting and stopping point. 
If stopping point is IN THE MRD - not just z, but x,y, use these events as 'contained'.
If a straight line from starting point to stopping point projected to the mrd front and back passes within the mrd area, this is a completely penetrating track. use these events as 'penetrated'. 

number of events contained, penetrated, vs number of events (efficiency)
efficiency (contained) vs lepton angle (lab frame? neutrino frame?)
efficiency (penetrated+contained) vs lepton angle
efficienc vs *neutrino* energy

for contained events:
(fractional) number of events that penetrated into layer 1,2,3...11 of the MRD -> do an integral, so the
fraction of events that penetrated ito layer 11 is the stopping power of the MRD (efficiency)
plot penetrating tracks in the overflow bin.

plot penetration depth (cm) vs energy of muon

//  
*/


// ================ End and Close
t->ResetBranchAddresses();
//f->Close();
// ================ Cleanup
// events
delete b;
delete m;
delete v;
//// canvases
//delete hitCountCanv;
//delete QvsTtankCanv;
//delete // tank histograms
//delete tubeshithisttank;
//delete truehitcounthisttank;
//delete digihitcountshisttank;
//delete PEdisttank;
//delete PMThitfrequencytank;
//delete QvsTtank;
//// mrd histograms
//delete tubeshithistmrd;
//delete truehitcounthistmrd;
//delete digihitcountshistmrd;
//delete PEdistmrd;
//delete PMThitfrequencymrd;
//delete QvsTmrd;
//// veto
//delete tubeshithistveto;
//delete truehitcounthistveto;
//delete digihitcountshistveto;
//delete PEdistveto;
//delete PMThitfrequencyveto;
//delete QvsTveto;

}

//TODO
/*
for(int i=0; i<100;i++){t->GetEntry(i);atrigm=m->GetTrigger(0); int numdigis = atrigm->GetNcherenkovdigihits(); cout<<numdigis<<", ";if(i%10==0){cout<<endl;}if(numdigis!=0){break;}}
*/
//TH1F* histpoint = new TH1F("histpoint",TString::Format("Number of %s PMTs Hit;Number of tubes;Frequency",detectorElement), 100, 0., numpmts);

//std::string intxnumtotype(gst genieeventasclass){ //TODO

/* we have the following branch members:
  int    brIev         = 0;      // Event number 
  int    brNeutrino    = 0;      // Neutrino pdg code
  int    brFSPrimLept  = 0;      // Final state primary lepton pdg code
  int    brTarget      = 0;      // Nuclear target pdg code (10LZZZAAAI)
  int    brTargetZ     = 0;      // Nuclear target Z (extracted from pdg code above)
  int    brTargetA     = 0;      // Nuclear target A (extracted from pdg code above)
  int    brHitNuc      = 0;      // Hit nucleon pdg code      (not set for COH,IMD and NuEL events)
  int    brHitQrk      = 0;      // Hit quark pdg code        (set for DIS events only)
  bool   brFromSea     = false;  // Hit quark is from sea     (set for DIS events only)
  int    brResId       = 0;      // Produced baryon resonance (set for resonance events only)
  bool   brIsQel       = false;  // Is QEL?
  bool   brIsRes       = false;  // Is RES?
  bool   brIsDis       = false;  // Is DIS?
  bool   brIsCoh       = false;  // Is Coherent?
  bool   brIsMec       = false;  // Is MEC?
  bool   brIsDfr       = false;  // Is Diffractive?
  bool   brIsImd       = false;  // Is IMD?
  bool   brIsSingleK   = false;  // Is single kaon?  
  bool   brIsImdAnh    = false;  // Is IMD annihilation?
  bool   brIsNuEL      = false;  // Is ve elastic?
  bool   brIsEM        = false;  // Is EM process?
  bool   brIsCC        = false;  // Is Weak CC process?
  bool   brIsNC        = false;  // Is Weak NC process?
  bool   brIsCharmPro  = false;  // Produces charm?
  int    brCodeNeut    = 0;      // The equivalent NEUT reaction code (if any)
  int    brCodeNuance  = 0;      // The equivalent NUANCE reaction code (if any)
  double brWeight      = 0;      // Event weight
  double brKineXs      = 0;      // Bjorken x as was generated during kinematical selection; takes fermi momentum / off-shellness into account
  double brKineYs      = 0;      // Inelasticity y as was generated during kinematical selection; takes fermi momentum / off-shellness into account
  double brKineTs      = 0;      // Energy transfer to nucleus at COH events as was generated during kinematical selection
  double brKineQ2s     = 0;      // Momentum transfer Q^2 as was generated during kinematical selection; takes fermi momentum / off-shellness into account
  double brKineWs      = 0;      // Hadronic invariant mass W as was generated during kinematical selection; takes fermi momentum / off-shellness into account
  double brKineX       = 0;      // Experimental-like Bjorken x; neglects fermi momentum / off-shellness 
  double brKineY       = 0;      // Experimental-like inelasticity y; neglects fermi momentum / off-shellness 
  double brKineT       = 0;      // Experimental-like energy transfer to nucleus at COH events 
  double brKineQ2      = 0;      // Experimental-like momentum transfer Q^2; neglects fermi momentum / off-shellness
  double brKineW       = 0;      // Experimental-like hadronic invariant mass W; neglects fermi momentum / off-shellness 
  double brEvRF        = 0;      // Neutrino energy @ the rest-frame of the hit-object (eg nucleon for CCQE, e- for ve- elastic,...)
  double brEv          = 0;      // Neutrino energy @ LAB
  double brPxv         = 0;      // Neutrino px @ LAB
  double brPyv         = 0;      // Neutrino py @ LAB
  double brPzv         = 0;      // Neutrino pz @ LAB
  double brEn          = 0;      // Initial state hit nucleon energy @ LAB
  double brPxn         = 0;      // Initial state hit nucleon px @ LAB
  double brPyn         = 0;      // Initial state hit nucleon py @ LAB
  double brPzn         = 0;      // Initial state hit nucleon pz @ LAB
  double brEl          = 0;      // Final state primary lepton energy @ LAB
  double brPxl         = 0;      // Final state primary lepton px @ LAB
  double brPyl         = 0;      // Final state primary lepton py @ LAB
  double brPzl         = 0;      // Final state primary lepton pz @ LAB
  double brPl          = 0;      // Final state primary lepton p  @ LAB
  double brCosthl      = 0;      // Final state primary lepton cos(theta) wrt to neutrino direction
  int    brNfP         = 0;      // Nu. of final state p's + \bar{p}'s (after intranuclear rescattering)
  int    brNfN         = 0;      // Nu. of final state n's + \bar{n}'s
  int    brNfPip       = 0;      // Nu. of final state pi+'s
  int    brNfPim       = 0;      // Nu. of final state pi-'s
  int    brNfPi0       = 0;      // Nu. of final state pi0's (
  int    brNfKp        = 0;      // Nu. of final state K+'s
  int    brNfKm        = 0;      // Nu. of final state K-'s
  int    brNfK0        = 0;      // Nu. of final state K0's + \bar{K0}'s
  int    brNfEM        = 0;      // Nu. of final state gammas and e-/e+ 
  int    brNfOther     = 0;      // Nu. of heavier final state hadrons (D+/-,D0,Ds+/-,Lamda,Sigma,Lamda_c,Sigma_c,...)
  int    brNiP         = 0;      // Nu. of `primary' (: before intranuclear rescattering) p's + \bar{p}'s  
  int    brNiN         = 0;      // Nu. of `primary' n's + \bar{n}'s  
  int    brNiPip       = 0;      // Nu. of `primary' pi+'s 
  int    brNiPim       = 0;      // Nu. of `primary' pi-'s 
  int    brNiPi0       = 0;      // Nu. of `primary' pi0's 
  int    brNiKp        = 0;      // Nu. of `primary' K+'s  
  int    brNiKm        = 0;      // Nu. of `primary' K-'s  
  int    brNiK0        = 0;      // Nu. of `primary' K0's + \bar{K0}'s 
  int    brNiEM        = 0;      // Nu. of `primary' gammas and e-/e+ 
  int    brNiOther     = 0;      // Nu. of other `primary' hadron shower particles
  int    brNf          = 0;      // Nu. of final state particles in hadronic system
  int    brPdgf  [kNPmax];       // Pdg code of k^th final state particle in hadronic system
  double brEf    [kNPmax];       // Energy     of k^th final state particle in hadronic system @ LAB
  double brPxf   [kNPmax];       // Px         of k^th final state particle in hadronic system @ LAB
  double brPyf   [kNPmax];       // Py         of k^th final state particle in hadronic system @ LAB
  double brPzf   [kNPmax];       // Pz         of k^th final state particle in hadronic system @ LAB
  double brPf    [kNPmax];       // P          of k^th final state particle in hadronic system @ LAB
  double brCosthf[kNPmax];       // cos(theta) of k^th final state particle in hadronic system @ LAB wrt to neutrino direction
  int    brNi          = 0;      // Nu. of particles in 'primary' hadronic system (before intranuclear rescattering)
  int    brPdgi[kNPmax];         // Pdg code of k^th particle in 'primary' hadronic system 
  int    brResc[kNPmax];         // FSI code of k^th particle in 'primary' hadronic system 
  double brEi  [kNPmax];         // Energy   of k^th particle in 'primary' hadronic system @ LAB
  double brPxi [kNPmax];         // Px       of k^th particle in 'primary' hadronic system @ LAB
  double brPyi [kNPmax];         // Py       of k^th particle in 'primary' hadronic system @ LAB
  double brPzi [kNPmax];         // Pz       of k^th particle in 'primary' hadronic system @ LAB
  double brVtxX;                 // Vertex x in detector coord system (SI)
  double brVtxY;                 // Vertex y in detector coord system (SI)
  double brVtxZ;                 // Vertex z in detector coord system (SI)
  double brVtxT;                 // Vertex t in detector coord system (SI)
  double brSumKEf;               // Sum of kinetic energies of all final state particles
  double brCalResp0;             // Approximate calorimetric response to the hadronic system computed as sum of
				 //  - (kinetic energy) for pi+, pi-, p, n 
                                 //  - (energy + 2*mass) for antiproton, antineutron
                                 //  - ((e/h) * energy)   for pi0, gamma, e-, e+, where e/h is set to 1.3
                                 //  - (kinetic energy) for other particles
*/


//TODO
/*
scattering type
interaction type
neutrino energy (RfLab) 
neutrino type
hit nucleon type
hit nucleus type (Z,A, name (O16) or PDG)
hit nucleus energy
primary fs lepton energy
primary fs lepton type
Q^2
angle or Cos(θ)
# additional FS protons
# additional FS neutrons
# additonal FS pi+/0/-

  int    brNeutrino    = 0;      // Neutrino pdg code
  int    brFSPrimLept  = 0;      // Final state primary lepton pdg code
  int    brTarget      = 0;      // Nuclear target pdg code (10LZZZAAAI)
  int    brTargetZ     = 0;      // Nuclear target Z (extracted from pdg code above)
  int    brTargetA     = 0;      // Nuclear target A (extracted from pdg code above)
  int    brHitNuc      = 0;      // Hit nucleon pdg code      (not set for COH,IMD and NuEL events)
  int    brHitQrk      = 0;      // Hit quark pdg code        (set for DIS events only)
  bool   brFromSea     = false;  // Hit quark is from sea     (set for DIS events only)
  int    brResId       = 0;      // Produced baryon resonance (set for resonance events only)
  bool   brIsQel       = false;  // Is QEL?
  bool   brIsRes       = false;  // Is RES?
  bool   brIsDis       = false;  // Is DIS?
  bool   brIsCoh       = false;  // Is Coherent?
  bool   brIsMec       = false;  // Is MEC?
  bool   brIsDfr       = false;  // Is Diffractive?
  bool   brIsImd       = false;  // Is IMD?
  bool   brIsSingleK   = false;  // Is single kaon?  
  bool   brIsImdAnh    = false;  // Is IMD annihilation?
  bool   brIsNuEL      = false;  // Is ve elastic?
  bool   brIsEM        = false;  // Is EM process?
  bool   brIsCC        = false;  // Is Weak CC process?
  bool   brIsNC        = false;  // Is Weak NC process?
  bool   brIsCharmPro  = false;  // Produces charm?
  int    brCodeNeut    = 0;      // The equivalent NEUT reaction code (if any)
  int    brCodeNuance  = 0;      // The equivalent NUANCE reaction code (if any)
  double brWeight      = 0;      // Event weight
  double brKineXs      = 0;      // Bjorken x as was generated during kinematical selection; takes fermi momentum / off-shellness into account
  double brKineYs      = 0;      // Inelasticity y as was generated during kinematical selection; takes fermi momentum / off-shellness into account
  double brKineTs      = 0;      // Energy transfer to nucleus at COH events as was generated during kinematical selection
  double brKineQ2s     = 0;      // Momentum transfer Q^2 as was generated during kinematical selection; takes fermi momentum / off-shellness into account
  double brKineWs      = 0;      // Hadronic invariant mass W as was generated during kinematical selection; takes fermi momentum / off-shellness into account
  double brKineX       = 0;      // Experimental-like Bjorken x; neglects fermi momentum / off-shellness 
  double brKineY       = 0;      // Experimental-like inelasticity y; neglects fermi momentum / off-shellness 
  double brKineT       = 0;      // Experimental-like energy transfer to nucleus at COH events 
  double brKineQ2      = 0;      // Experimental-like momentum transfer Q^2; neglects fermi momentum / off-shellness
  double brKineW       = 0;      // Experimental-like hadronic invariant mass W; neglects fermi momentum / off-shellness 
  double brEvRF        = 0;      // Neutrino energy @ the rest-frame of the hit-object (eg nucleon for CCQE, e- for ve- elastic,...)
  double brEv          = 0;      // Neutrino energy @ LAB
  double brPxv         = 0;      // Neutrino px @ LAB
  double brPyv         = 0;      // Neutrino py @ LAB
  double brPzv         = 0;      // Neutrino pz @ LAB
  double brEn          = 0;      // Initial state hit nucleon energy @ LAB
  double brPxn         = 0;      // Initial state hit nucleon px @ LAB
  double brPyn         = 0;      // Initial state hit nucleon py @ LAB
  double brPzn         = 0;      // Initial state hit nucleon pz @ LAB
  double brEl          = 0;      // Final state primary lepton energy @ LAB
  double brPxl         = 0;      // Final state primary lepton px @ LAB
  double brPyl         = 0;      // Final state primary lepton py @ LAB
  double brPzl         = 0;      // Final state primary lepton pz @ LAB
  double brPl          = 0;      // Final state primary lepton p  @ LAB
  double brCosthl      = 0;      // Final state primary lepton cos(theta) wrt to neutrino direction
  int    brNfP         = 0;      // Nu. of final state p's + \bar{p}'s (after intranuclear rescattering)
  int    brNfN         = 0;      // Nu. of final state n's + \bar{n}'s
  int    brNfPip       = 0;      // Nu. of final state pi+'s
  int    brNfPim       = 0;      // Nu. of final state pi-'s
  int    brNfPi0       = 0;      // Nu. of final state pi0's (
  int    brNfKp        = 0;      // Nu. of final state K+'s
  int    brNfKm        = 0;      // Nu. of final state K-'s
  int    brNfK0        = 0;      // Nu. of final state K0's + \bar{K0}'s
  int    brNfEM        = 0;      // Nu. of final state gammas and e-/e+ 
  int    brNfOther     = 0;      // Nu. of heavier final state hadrons (D+/-,D0,Ds+/-,Lamda,Sigma,Lamda_c,Sigma_c,...)
  int    brNiP         = 0;      // Nu. of `primary' (: before intranuclear rescattering) p's + \bar{p}'s  
  int    brNiN         = 0;      // Nu. of `primary' n's + \bar{n}'s  
  int    brNiPip       = 0;      // Nu. of `primary' pi+'s 
  int    brNiPim       = 0;      // Nu. of `primary' pi-'s 
  int    brNiPi0       = 0;      // Nu. of `primary' pi0's 
  int    brNiKp        = 0;      // Nu. of `primary' K+'s  
  int    brNiKm        = 0;      // Nu. of `primary' K-'s  
  int    brNiK0        = 0;      // Nu. of `primary' K0's + \bar{K0}'s 
  int    brNiEM        = 0;      // Nu. of `primary' gammas and e-/e+ 
  int    brNiOther     = 0;      // Nu. of other `primary' hadron shower particles
  int    brNf          = 0;      // Nu. of final state particles in hadronic system
  int    brPdgf  [kNPmax];       // Pdg code of k^th final state particle in hadronic system
  double brEf    [kNPmax];       // Energy     of k^th final state particle in hadronic system @ LAB
  double brPxf   [kNPmax];       // Px         of k^th final state particle in hadronic system @ LAB
  double brPyf   [kNPmax];       // Py         of k^th final state particle in hadronic system @ LAB
  double brPzf   [kNPmax];       // Pz         of k^th final state particle in hadronic system @ LAB
  double brPf    [kNPmax];       // P          of k^th final state particle in hadronic system @ LAB
  double brCosthf[kNPmax];       // cos(theta) of k^th final state particle in hadronic system @ LAB wrt to neutrino direction
  int    brNi          = 0;      // Nu. of particles in 'primary' hadronic system (before intranuclear rescattering)
  int    brPdgi[kNPmax];         // Pdg code of k^th particle in 'primary' hadronic system 
  int    brResc[kNPmax];         // FSI code of k^th particle in 'primary' hadronic system 
  double brEi  [kNPmax];         // Energy   of k^th particle in 'primary' hadronic system @ LAB
  double brPxi [kNPmax];         // Px       of k^th particle in 'primary' hadronic system @ LAB
  double brPyi [kNPmax];         // Py       of k^th particle in 'primary' hadronic system @ LAB
  double brPzi [kNPmax];         // Pz       of k^th particle in 'primary' hadronic system @ LAB
  double brVtxX;                 // Vertex x in detector coord system (SI)
  double brVtxY;                 // Vertex y in detector coord system (SI)
  double brVtxZ;                 // Vertex z in detector coord system (SI)
  double brVtxT;                 // Vertex t in detector coord system (SI)
  double brSumKEf;               // Sum of kinetic energies of all final state particles
  double brCalResp0;             // Approximate calorimetric response to the hadronic system computed as sum of
				 //  - (kinetic energy) for pi+, pi-, p, n 
                                 //  - (energy + 2*mass) for antiproton, antineutron
                                 //  - ((e/h) * energy)   for pi0, gamma, e-, e+, where e/h is set to 1.3
                                 //  - (kinetic energy) for other particles
*/
