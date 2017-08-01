/* vim:set noexpandtab tabstop=4 wrap */
//############################################################################################

// INITIALISE ENVIRONMENT
// ======================
void WCSimAnalysis::InitEnvironment(){
	//gROOT->Reset();
	gROOT->ProcessLine("#include <regex>");
	TString pwd = gSystem->Getenv("PWD");
	TString librarypath = pwd + "libWCSimRoot.so";
	//"/home/marc/LinuxSystemFiles/WCSim/gitver/wcsim/libWCSimRoot.so"
	//gSystem->Load(librarypath);	//also include TSystem.h
	//R__LOAD_LIBRARY(libWCSimRoot.so);
	//gStyle->SetOptStat(0); 				// disable the stats box, or use hist->SetBit(TH1::kNoStats);
	//gStyle->SetOptStat(111111);			// show overflow and underflow contents in the stats box
	ColourPlotStyle();
}

//############################################################################################

// LOAD INPUT FILES
// ================
void WCSimAnalysis::LoadInputFiles(){
	cout<<"Setting output directory to "<<outputdir<<endl;
	cout<<"Loading files from"<<inputdir<<"..."<<endl;
	//TFile* f = TFile::Open("wcsim6.root");
	//TTree* t = (TTree*)f->Get("wcsimT");

	/*TChain* */ t = new TChain("wcsimT");
	// loop over subdirectories and add all the files in each
	//const char* inputdir = "/home/marc/anniegpvm/stats10k/";	// inputdir is class member 
	const char* ext = ".root";
	TString nextfilepattern;
	TSystemDirectory dir(inputdir, inputdir);
	TList *subfolders = dir.GetListOfFiles(); // despite name returns files and folders
	if(subfolders&&(subfolders->GetEntries()>2)) { // always at least 2 subfolders: '.' and '..'
		cout<<"looping over subfolders"<<endl;
		TSystemDirectory *subfolder;
		TIter nextsf(subfolders);
		while ((subfolder=(TSystemDirectory*)nextsf())) {
			if (subfolder->IsDirectory()){
				TString sfname = subfolder->GetName();
				if(sfname=="."||sfname==".."){ continue; }
				cout<<"found subfolder '"<<sfname<<"'"<<endl;
				/*
				// don't actually need this information
				TList *files = subfolder->GetListOfFiles();
				TSystemFile *file;
				TIter nextf(files);
				while (file=(TSystemFile*)nextf()) {
					TString fname;
					fname = file->GetName();
					if(fname.EndsWith(ext)) { int i=0; }
				}
				*/
				nextfilepattern = inputdir + sfname + "/wcsim_*"; // XXX
				cout<<"adding "<<nextfilepattern<<endl;
				t->Add(nextfilepattern.Data());
				// cf. t->Add("/home/marc/anniegpvm/stats10k/stats_1_79283/wcsim_*");
			} // end if directory
		} // end loop over subdirectories
	} // end if subfolders exist.
	if(1 /*t->GetEntriesFast()==0*/) {	// add files in parent directory
		// can't '+' const chars, convert to string
		nextfilepattern = std::string(inputdir) + std::string("/wcsim_0.*.root"); // XXX /wcsim_0.*
		cout<<"adding "<<nextfilepattern<<endl;
		t->Add(nextfilepattern.Data());
	}
	cout<<"Loaded "<<t->GetEntriesFast()<<" entries"<<endl;	// fast doesn't work with a tchain

	// first let's just get the geotree now, since we only need to retrieve it once
	// can be used to get PMTs by ID: note GetPMT(i) returns a PMT object not a pointer!
	Long64_t localEntry = t->LoadTree(0);
	TTree* currenttree = t->GetTree();
	//TFile* firstfile = currenttree->GetCurrentFile();
	//TFile* firstfile = TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_0.root");
	TFile* firstfile = TFile::Open("/pnfs/annie/persistent/users/moflaher/wcsim_wdirt_17-06-17/wcsim_0.1000.root");
	cout<<"Loading geometry from file "<<firstfile->GetName()<<endl;
	if(firstfile==0){ cout<<"NO GEOMETRY FILE?"<<endl; assert(false);}
	TTree* geotree = (TTree*)firstfile->Get("wcsimGeoT");
	if(geotree==0){ cout<<"NO GEOMETRY IN FIRST FILE?"<<endl; assert(false); }
	//WCSimRootGeom* geo = 0; 
	geotree->SetBranchAddress("wcsimrootgeom", &geo);
	if (geotree->GetEntries() == 0) { cout<<"geotree has no entries!"<<endl; exit(9); }
	geotree->GetEntry(0);
	cout<<"Geometry loaded"<<endl;
}

//############################################################################################

// GET TREE BRANCHES
// =================
void WCSimAnalysis::GetTreeData(){
	b = new WCSimRootEvent();
	m = new WCSimRootEvent();
	v = new WCSimRootEvent();
	bp=0, mp=0, vp=0;
	Long64_t localEntry = t->LoadTree(0);
	if(localEntry<0){cout<<"no trees in the chain!"<<endl; assert(false);}
	t->SetBranchAddress("wcsimrootevent",&b, &bp);
	t->SetBranchAddress("wcsimrootevent_mrd",&m, &mp);
	t->SetBranchAddress("wcsimrootevent_facc",&v, &vp);
	if(bp==0||mp==0||vp==0){ cout<<"first set of branches are zombies!"<<endl; 
		cout<<"bp = "<<bp<<", mp = "<<mp<<", vp = "<<vp<<endl; 
		return; 
	}
	t->GetEntry(0);
	Int_t firstentries = bp->GetEntries();
	if(firstentries==0){cout<<"no entries in the first tree!"<<endl; assert(false);}
}

//############################################################################################

// LOAD NEXT TCHAIN ENTRY
// ======================
int WCSimAnalysis::LoadTchainEntry(Int_t &eventnum){
	while(1){  // keep loading entries until the file number is > the start file offset
		Long64_t localEntry = t->LoadTree(eventnum);
		if (localEntry<0){ cout<<"breaking at entry "<<eventnum<<" as localEntry<0"<<endl; return 0; }
		Int_t nextTreeNumber = t->GetTreeNumber();
		if(treeNumber!=nextTreeNumber){
			cout<< "Reached end of Tree. Last entries' tree number was "
				<< treeNumber <<", this entries' tree number is "<< nextTreeNumber <<endl;
			cout<<"entries in this tree: "<<bp->GetEntries()<<endl; 
			currenttree = t->GetTree();
			currentfile = currenttree->GetCurrentFile();
			currentfilestring = std::string(currentfile->GetName());
		
			// filename is of the form "wcsim_0.####.root"  //annie_tank_flux.####.root
			// #### is input file num. Need this to match against genie/wcsim file names
			std::match_results<string::const_iterator> submatches;
			std::regex theexpression (".*/?[^\\.]+\\.([0-9]+)\\.root");
			cout<<"matching regex for filename "<<currentfilestring<<endl;
			std::regex_match (currentfilestring, submatches, theexpression);
			std::string submatch = (std::string)submatches[0];	// match 0 is 'whole match' or smthg
			if(submatch==""){ cout<<"unrecognised input file pattern: "<<currentfilestring<<endl; return -1; }
			submatch = (std::string)submatches[1];
			cout<<"extracted submatch is "<<submatch<<endl;
			wcsimfilenum = atoi(submatch.c_str());
		
			if(wcsimfilenum<firstfile){
				cout<<"skipping file "<<wcsimfilenum<<endl;
				eventnum+=bp->GetEntries(); 
				continue; 
			}
		
			// do need to re-set branch addresses
			t->SetBranchAddress("wcsimrootevent",&b, &bp);
			t->SetBranchAddress("wcsimrootevent_mrd",&m, &mp);
			t->SetBranchAddress("wcsimrootevent_facc",&v, &vp);
			if(bp==0||mp==0||vp==0){ cout<<"branches are zombies!"<<endl; }
	
			// open new MRD Track and Veto Event files
			OpenMRDtrackOutfile(wcsimfilenum);
			OpenFACCtrackOutfile(wcsimfilenum);
			treeNumber=nextTreeNumber;
			break;
		}
		break;
	} // end of do while - break after skipping selected entries.
//	G4cout<<"Loading data from entry "<<inputEntry<<", localentry "<<localEntry<<"/"<<entriesInThisTree<<G4endl;
	t->GetEntry(eventnum);
	return 1;
}

void WCSimAnalysis::ColourPlotStyle(){
	const Int_t NRGBs = 5;
	const Int_t NCont = 255;

	Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
	Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
	Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
	Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
	gStyle->SetNumberContours(NCont);
}
