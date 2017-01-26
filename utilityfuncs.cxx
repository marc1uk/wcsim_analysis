/* vim:set noexpandtab tabstop=4 wrap */
//############################################################################################

// INITIALISE ENVIRONMENT
// ======================
void WCSimAnalysis::InitEnvironment(){
	//gROOT->Reset();
	gROOT->ProcessLine("#include <regex>");
	gSystem->Load("/home/marc/LinuxSystemFiles/WCSim/gitver/wcsim/libWCSimRoot.so");	//also include TSystem.h
	//R__LOAD_LIBRARY(libWCSimRoot.so);
}

//############################################################################################

// LOAD INPUT FILES
// ================
void WCSimAnalysis::LoadInputFiles(){
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
		while (subfolder=(TSystemDirectory*)nextsf()) {
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
				nextfilepattern = inputdir + sfname + "/wcsim_*";
				cout<<"adding "<<nextfilepattern<<endl;
				t->Add(nextfilepattern.Data());
				// cf. t->Add("/home/marc/anniegpvm/stats10k/stats_1_79283/wcsim_*");
			} // end if directory
		} // end loop over subdirectories
	} // end if subfolders exist.
	if(1 /*t->GetEntriesFast()==0*/) {	// add files in parent directory
		// can't '+' const chars, convert to string
		nextfilepattern = std::string(inputdir) + std::string("/wcsim_*");
		cout<<"adding "<<nextfilepattern<<endl;
		t->Add(nextfilepattern.Data());
	}
	cout<<"Loaded "<<t->GetEntriesFast()<<" entries"<<endl;

	// first let's just get the geotree now, since we only need to retrieve it once
	// can be used to get PMTs by ID: note GetPMT(i) returns a PMT object not a pointer!
	Long64_t localEntry = t->LoadTree(0);
	TTree* currenttree = t->GetTree();
	//TFile* firstfile = currenttree->GetCurrentFile();
	TFile* firstfile = TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_0.root");
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
	cout<<"Retrieving event branches"<<endl;
	b = new WCSimRootEvent();
	m = new WCSimRootEvent();
	v = new WCSimRootEvent();
	bp=0, mp=0, vp=0;
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
	atrigt = b->GetTrigger(0);
	atrigm = m->GetTrigger(0);
	atrigv = v->GetTrigger(0);
}

//############################################################################################

// LOAD NEXT TCHAIN ENTRY
// ======================
int WCSimAnalysis::LoadTchainEntry(Int_t eventnum){
	Long64_t localEntry = t->LoadTree(eventnum);
	if (localEntry<0){ cout<<"breaking at entry "<<eventnum<<" as localEntry<0"<<endl; return 0; }
	Int_t nextTreeNumber = t->GetTreeNumber();
	if(treeNumber!=nextTreeNumber){
		cout<< "Reached end of Tree. Last entries' tree number was "
		    << treeNumber <<", this entries' tree number is "<< nextTreeNumber <<endl;
		cout<<"entries in this tree: "<<bp->GetEntries()<<endl; 
		// do need to re-set branch addresses
		t->SetBranchAddress("wcsimrootevent",&b, &bp);
		t->SetBranchAddress("wcsimrootevent_mrd",&m, &mp);
		t->SetBranchAddress("wcsimrootevent_facc",&v, &vp);
		if(bp==0||mp==0||vp==0){ cout<<"branches are zombies!"<<endl; }
		treeNumber=nextTreeNumber;
	}
//	G4cout<<"Loading data from entry "<<inputEntry<<", localentry "<<localEntry<<"/"<<entriesInThisTree<<G4endl;
	t->GetEntry(eventnum);
	return 1;
}
