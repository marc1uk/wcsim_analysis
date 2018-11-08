/* vim:set noexpandtab tabstop=4 wrap */
//############################################################################################
#include "type_name_as_string.hh" // this is some code from StackExchange that defines a function that can print a variables type
#include <sys/types.h> // for stat() test to see if file or folder
#include <sys/stat.h>
#include <unistd.h>
#include <fstream>

// INITIALISE ENVIRONMENT
// ======================
void WCSimAnalysis::InitEnvironment(){
	//gROOT->Reset();
	gROOT->ProcessLine("#include <regex>");
	TString pwd = gSystem->Getenv("PWD");
	//TString librarypath = pwd + "libWCSimRoot.so";
//	TString libraryloc = pwd + "../wcsim/";
//	TString librarypath = libraryloc + "libWCSimRoot.so";
//	gSystem->Load(librarypath);
//	gInterpreter->AddIncludePath(libraryloc);
	//R__LOAD_LIBRARY(../wcsim/libWCSimRoot.so);
	//gStyle->SetOptStat(0); 				// disable the stats box, or use hist->SetBit(TH1::kNoStats);
	//gStyle->SetOptStat(111111);			// show overflow and underflow contents in the stats box
	ColourPlotStyle();
}

//############################################################################################

// LOAD INPUT FILES
// ================
void WCSimAnalysis::LoadInputFiles(){
	cout<<"Setting output directory to "<<outputdir<<endl;
	cout<<"Loading files from "<<inputdir<<"..."<<endl;
	//TFile* f = TFile::Open("wcsim6.root");
	//TTree* t = (TTree*)f->Get("wcsimT");
	
	/*TChain* */ t = new TChain("wcsimT");
	

	bool isdir, addsubfolders=true;
	struct stat s;
	if(stat(inputdir,&s)==0){
		if(s.st_mode & S_IFDIR){        // mask to extract if it's a directory?? how does this work?
			isdir=true;  //it's a directory
		} else if(s.st_mode & S_IFREG){ // mask to check if it's a file??
			isdir=false; //it's a file
		} else {
			assert(false&&"Check input path: stat says it's neither file nor directory..?");
		}
	} else {
		//assert(false&&"stat failed on input path! Is it valid?"); // error
		// errors could be because this is a file pattern: e.g. wcsim_0.4*.root Treat as a file.
		isdir=false;
	}
	
	if(isdir){
		// loop over subdirectories and add all the files in each
		//const char* inputdir = "/home/marc/anniegpvm/stats10k/";	// inputdir is class member 
		const char* ext = ".root";
		TString nextfilepattern;
		TSystemDirectory dir(inputdir, inputdir);
		TList *subfolders = dir.GetListOfFiles(); // despite name returns files and folders
		if(addsubfolders&&subfolders&&(subfolders->GetEntries()>2)) { // always 2 subfolders: '.' and '..'
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
					nextfilepattern = TString::Format("%s/%s/%s",inputdir,sfname.Data(),"wcsim_*"); // XXX
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
	} else {
		// path passed is just a file
		cout<<"loading single file "<<inputdir<<endl;
		t->Add(inputdir);
	}
	
	// first let's just get the geotree now, since we only need to retrieve it once
	// can be used to get PMTs by ID: note GetPMT(i) returns a PMT object not a pointer!
	Long64_t localEntry = t->LoadTree(0);
	TTree* currenttree = t->GetTree();
	TFile* firstfile = currenttree->GetCurrentFile();
	//TFile* firstfile = TFile::Open("/home/marc/LinuxSystemFiles/WCSim/gitver/build/wcsim_0.root");
	//TFile* firstfile = TFile::Open("/pnfs/annie/persistent/users/moflaher/wcsim_wdirt_17-06-17/wcsim_0.1000.root");
	cout<<"Loading geometry from file "<<firstfile->GetName()<<endl;
	if(firstfile==0){ cout<<"NO GEOMETRY FILE?"<<endl; assert(false);}
	TTree* geotree = (TTree*)firstfile->Get("wcsimGeoT");
	if(geotree==0){ cout<<"NO GEOMETRY IN FIRST FILE?"<<endl; assert(false); }
	//WCSimRootGeom* geo = 0; 
	geotree->SetBranchAddress("wcsimrootgeom", &geo);
	if (geotree->GetEntries() == 0) { cout<<"geotree has no entries!"<<endl; exit(9); }
	geotree->GetEntry(0);
	cout<<"Geometry loaded"<<endl;
	
	cout<<"Loading WCSim options from file "<<firstfile->GetName()<<endl;
	if(firstfile==0){ cout<<"NO ROOT OPTIONS FILE?"<<endl; assert(false);}
	TTree *opttree = (TTree*)firstfile->Get("wcsimRootOptionsT");
	if(opttree==0){ cout<<"NO ROOT OPTIONS TREE IN FIRST FILE?"<<endl; assert(false); }
	//WCSimRootOptions *opt = 0; 
	opttree->SetBranchAddress("wcsimrootoptions", &opt);
	if (opttree->GetEntries() == 0) { cout<<"opttree has no entries!"<<endl; exit(9); }
	opttree->GetEntry(0);
	//opt->Print();
	pre_trigger_window_ns = -1.*opt->GetNDigitsPreTriggerWindow(); // NOTE: OPTIONS STORE AS NEGATIVE. 
	post_trigger_window_ns = opt->GetNDigitsPostTriggerWindow();
	
	// Load in configuration options for WCSimAnalysis
	cout<<"Loading WCSimAnalysis options from "<<optionsfile<<endl;
	std::ifstream file(optionsfile.c_str());
	std::string line;
	if(file.is_open()){
		while (getline(file,line)){
			if (line.size()>0){
				if (line.at(0)=='#')continue;
				std::string key;
				std::string value="";
				std::stringstream stream(line);
				if(stream>>key){
					std::string valuetemp; // combine everything after the key - allow spaces in values
					while(stream>>valuetemp){ value+=valuetemp; value+=" "; }
					//m_variables[key]=value;  // put into a map? must all be same type.
					if(key=="startDate") startDate = value;
					if(key=="triggerOffset") triggeroffset = stoi(value);
				}
			}
		}
	}
	file.close();
	
	// ***********************************
	// TIME TO SET THE CONSTANTS
	// ***********************************	
	
	// geometry constants
	numtankpmts=geo->GetWCNumPMT();
	nummrdpmts=geo->GetWCNumMRDPMT();
	numvetopmts=geo->GetWCNumFACCPMT();
	// FIXME DO NOT HARD CODE THESE FIXME
	caparraysize=8;             // pmts on the cap form an nxn grid where caparraysize=n FIXME wrong! Not square!
	pmtsperring=16;             // pmts around each ring of the main walls
	numpmtrings=6;              // num rings around the main walls
	
	// MRD track finding constants
	MAXTRACKSPEREVENT=50;       //
	maxsubeventduration=30.;    // in ns?
	
	
	// RAW file DAQ constants
	ADC_NS_PER_SAMPLE=2;
	MRD_NS_PER_SAMPLE=4;
	MRD_TIMEOUT_NS=4200;
	MRD_TIMESTAMP_DELAY = static_cast<unsigned long long>(MRD_TIMEOUT_NS);
	ADC_INPUT_RESISTANCE = 50.;  // Ohm
	ADC_TO_VOLT = 2.415 / std::pow(2., 12);// * by this constant converts ADC counts to Volts
	PULSE_HEIGHT_FUDGE_FACTOR = (1./300.); // WHAT UNITS ARE DIGIT Q's IN?!?
	MAXEVENTSIZE=10;
	MAXTRIGGERSIZE=10;
	channels_per_tdc_card = 32;
	channels_per_adc_card = 4;
	num_adc_cards = (numtankpmts-1)/channels_per_adc_card + 1; // round up, requiring numtankpmts!=0
	// n.b. variant that doesn't require !=0, but could overflow:
	// num_adc_cards = (numtankpmts + channels_per_adc_card - 1) / channels_per_adc_card;
	
	// WCSim Note: EventSize = (datapoints per chan per minibuf / 4) => event duration must be a multiple of 8
	minibuffers_per_fullbuffer = 40;
	minibuffer_datapoints_per_channel = (pre_trigger_window_ns+post_trigger_window_ns) / ADC_NS_PER_SAMPLE;
	buffer_size = minibuffer_datapoints_per_channel * minibuffers_per_fullbuffer;
	full_buffer_size = buffer_size * channels_per_adc_card;
	emulated_event_size = minibuffer_datapoints_per_channel / 4.; // the 4 comes from ADC firmware stuff.
}

//############################################################################################

// GET TREE BRANCHES
// =================
void WCSimAnalysis::GetTreeData(){
	b = new WCSimRootEvent();
	b->Initialize();
	m = new WCSimRootEvent();
	m->Initialize();
	v = new WCSimRootEvent();
	v->Initialize();
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
			//std::regex theexpression (".*/?[^\\.]+\\.([0-9]+)\\.?([0-9]+)?\\.root");
			cout<<"matching regex for filename "<<currentfilestring<<endl;
			std::regex_match (currentfilestring, submatches, theexpression);
			// match 0 is 'whole match' or smthg
			std::string submatch = (submatches.size()) ? (std::string)submatches[0] : "";
			if(submatch==""){
				cout<<"unrecognised input file pattern: "<<currentfilestring
					<<", will set wcsimfilenum=0"<<endl;
					wcsimfilenum=0;
				//return -1;
			} else {
				submatch = (std::string)submatches[1];
				cout<<"extracted submatch is "<<submatch<<endl;
				wcsimfilenum = atoi(submatch.c_str());
				
				// check if we're to skip this file
				if(wcsimfilenum<firstfilenum){
					cout<<"skipping file "<<wcsimfilenum<<endl;
					eventnum+=bp->GetEntries(); 
					continue; 
				}
			}
			
			// do need to re-set branch addresses
			cout<<"setting branch addresses"<<endl;
			t->SetBranchAddress("wcsimrootevent",&b, &bp);
			t->SetBranchAddress("wcsimrootevent_mrd",&m, &mp);
			t->SetBranchAddress("wcsimrootevent_facc",&v, &vp);
			if(bp==0||mp==0||vp==0){ cout<<"branches are zombies!"<<endl; }
			
			// open new MRD Track and Veto Event files
			cout<<"opening MRD and FACC track output files"<<endl;
			OpenMRDtrackOutfile(wcsimfilenum);
			OpenFACCtrackOutfile(wcsimfilenum);
			treeNumber=nextTreeNumber;
			break;
		}
		break;
	} // end of do while - break after skipping selected entries.
	//cout<<"Loading entry "<<eventnum<<endl;
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
