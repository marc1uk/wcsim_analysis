{
const char* dirtpath="/home/marc/anniegpvm/pnfs/annie/persistent/users/moflaher/g4dirt";
const char* geniepath="/home/marc/anniegpvm/pnfs/annie/persistent/users/rhatcher/genie";
const char* wcsimpath="/home/marc/anniegpvm/pnfs/annie/persistent/users/moflaher/wcsim";
const char* wcsimlibrarypath="/home/marc/anniegpvm/annie/app/users/moflaher/wcsim/wcsim/libWCSimRoot.so";
const char* outpath="/home/marc/anniegpvm/annie/app/users/moflaher/wcsim/root_work";

const bool VERBOSE=false;

// TChain for dirt files - this will be the main driver of the loop - all it's events will be processed.
TChain* inputdata =  new TChain("tankflux");
TString chainpattern = TString::Format("%s/annie_tank_flux.100*.root",dirtpath);   // replaced *->2000
cout<<"loading TChain entries from "<<chainpattern<<endl;
inputdata->Add(chainpattern);
Int_t numents = inputdata->GetEntries();
cout<<"loaded "<<numents<<" entries in the chain"<<endl;
if(numents<1){ exit(-1); }

Long64_t localEntry;
Int_t inputEntry;
Int_t entriesInThisTree;
TBranch* runBranch=0, *vtxxBranch=0, *vtxyBranch=0, *vtxzBranch=0, *vtxtBranch=0, *pxBranch=0, *pyBranch=0, *pzBranch=0, *EBranch=0, *KEBranch=0, *pdgBranch=0, *nTankBranch=0, *nupdgBranch=0, *nuvtxxBranch=0, *nuvtxyBranch=0, *nuvtxzBranch=0, *nuvtxtBranch=0, *nuPVBranch=0, *nuvtxmatBranch=0, *nuprimaryBranch=0, *nufluxfilenameBranch=0, *genieentryBranch=0, *genierecordBranch=0;
Int_t runbranchval, entrybranchval, ntankbranchval, nupdgval, genieentrybranchval, pdgval, nuprimaryval;
Double_t vtxxval, vtxyval, vtxzval, vtxtval, pxval, pyval, pzval, eval, keval, nuvtxxval, nuvtxyval, nuvtxzval, nuvtxtval;
Int_t* pdgbranchval=0, *nuprimarybranchval=0;
Double_t* vtxxbranchval=0, *vtxybranchval=0, *vtxzbranchval=0, *vtxtbranchval=0, *pxbranchval=0, *pybranchval=0, *pzbranchval=0, *ebranchval=0, *kebranchval=0;
Char_t nupvval[100];
Char_t numatval[100];
Char_t nufluxfilenameval[100];

TH1F* vertexxvals = new TH1F("vertexxvals","X Vertices of Primaries",100,-160, 160);
TH1F* vertexyvals = new TH1F("vertexyvals","Y Vertices of Primaries",100,-300, 300);
TH1F* vertexzvals = new TH1F("vertexzvals","Z Vertices of Primaries",100,0, 350);

inputdata->LoadTree(0);
Int_t treeNumber = -1;
TTree* tankflux = inputdata->GetTree();
Int_t thistreesentries = tankflux->GetEntries();
cout<<thistreesentries<<" entries in the first tree"<<endl;

int tankeventnum=0;

for(Int_t inputEntry=0; inputEntry<numents; inputEntry++){
	/* 1. Load next g4dirt entry */ 
if(VERBOSE){
cout<<"loading entry "<<inputEntry<<endl;
}
	Long64_t localEntry = inputdata->LoadTree(inputEntry);
	if( localEntry<0){ cout<<"end of tchain"<<endl; break; }
	Int_t nextTreeNumber = inputdata->GetTreeNumber();
	if(treeNumber!=nextTreeNumber){
		if(VERBOSE){
		cout<< "Reached end of Tree. Last entries' tree number was "
		      << treeNumber <<", this entries' tree number is "<< nextTreeNumber <<endl;
		cout<<"Getting new tree branches"<<endl;
		}
		// this means we've switched file - need to load the new meta tree and genie tree.
		tankflux = inputdata->GetTree();
		thistreesentries = tankflux->GetEntries();
		if(VERBOSE){
		cout<<"tankflux has "<<thistreesentries<<" entries in this file"<<endl;
		}
		/* Set the branch addresses for the new trees */
		inputdata->SetBranchAddress("run",&runbranchval,&runBranch);
		inputdata->SetBranchAddress("ntank",&ntankbranchval,&nTankBranch);
		inputdata->SetBranchAddress("nupdg",&nupdgval,&nupdgBranch);
		inputdata->SetBranchAddress("nuvtxx",&nuvtxxval,&nuvtxxBranch);
		inputdata->SetBranchAddress("nuvtxy",&nuvtxyval,&nuvtxyBranch);
		inputdata->SetBranchAddress("nuvtxz",&nuvtxzval,&nuvtxzBranch);
		inputdata->SetBranchAddress("nuvtxt",&nuvtxtval,&nuvtxtBranch);
		inputdata->SetBranchAddress("vtxvol",&nupvval,&nuPVBranch);
		inputdata->SetBranchAddress("vtxmat",&numatval,&nuvtxmatBranch);
		inputdata->SetBranchAddress("entry",&genieentrybranchval,&genieentryBranch);
		
		vtxxBranch=inputdata->GetBranch("vx");
		vtxyBranch=inputdata->GetBranch("vy");
		vtxzBranch=inputdata->GetBranch("vz");
		vtxtBranch=inputdata->GetBranch("vt");
		pxBranch=inputdata->GetBranch("px");
		pyBranch=inputdata->GetBranch("py");
		pzBranch=inputdata->GetBranch("pz");
		EBranch=inputdata->GetBranch("E");
		KEBranch=inputdata->GetBranch("kE");
		pdgBranch=inputdata->GetBranch("pdgtank");
		nuprimaryBranch=inputdata->GetBranch("primary");
		entriesInThisTree = runBranch->GetEntries();
		
		
		if(runBranch==0||nTankBranch==0||vtxxBranch==0||vtxyBranch==0||vtxzBranch==0||vtxtBranch==0||pxBranch==0||pyBranch==0||pzBranch==0||EBranch==0||KEBranch==0||pdgBranch==0||nupdgBranch==0||nuvtxxBranch==0||nuvtxyBranch==0||nuvtxzBranch==0||nuvtxtBranch==0||nuPVBranch==0||nuvtxmatBranch==0||nuprimaryBranch==0){
			cout<<"BRANCHES ARE ZOMBIES ARGH!"<<endl;
		} else if (VERBOSE) { cout<<"entries in this tree: "<<vtxxBranch->GetEntries()<<endl; }
		treeNumber=nextTreeNumber;
	}
	if(VERBOSE){
	cout<<"Loading primaries from entry "<<inputEntry<<", localentry "
	      <<localEntry<<"/"<<entriesInThisTree<<endl;
	}
	
	nTankBranch->GetEntry(localEntry);
	nupdgBranch->GetEntry(localEntry);
	nuvtxxBranch->GetEntry(localEntry);
	nuvtxyBranch->GetEntry(localEntry);
	nuvtxzBranch->GetEntry(localEntry);
	nuvtxtBranch->GetEntry(localEntry);
	nuPVBranch->GetEntry(localEntry);
	nuvtxmatBranch->GetEntry(localEntry);
	genieentryBranch->GetEntry(localEntry);
	
	if(vtxxbranchval){delete[] vtxxbranchval;}
	if(vtxybranchval){delete[] vtxybranchval;}
	if(vtxzbranchval){delete[] vtxzbranchval;}
	if(vtxtbranchval){delete[] vtxtbranchval;}
	if(pxbranchval){delete[] pxbranchval;}
	if(pybranchval){delete[] pybranchval;}
	if(pzbranchval){delete[] pzbranchval;}
	if(ebranchval){delete[] ebranchval;}
	if(kebranchval){delete[] kebranchval;}
	if(pdgbranchval){delete[] pdgbranchval;}
	if(nuprimarybranchval){delete[] nuprimarybranchval;}
	
	vtxxbranchval = new Double_t[ntankbranchval];
	vtxybranchval = new Double_t[ntankbranchval];
	vtxzbranchval = new Double_t[ntankbranchval];
	vtxtbranchval = new Double_t[ntankbranchval];
	pxbranchval = new Double_t[ntankbranchval];
	pybranchval = new Double_t[ntankbranchval];
	pzbranchval = new Double_t[ntankbranchval];
	ebranchval = new Double_t[ntankbranchval];
	kebranchval = new Double_t[ntankbranchval];
	pdgbranchval = new Int_t[ntankbranchval];
	nuprimarybranchval = new Int_t[ntankbranchval];
	
	if(vtxxbranchval==0||vtxybranchval==0||vtxzbranchval==0||vtxtbranchval==0||pxbranchval==0||pybranchval==0||pzbranchval==0||ebranchval==0||kebranchval==0||pdgbranchval==0||nuprimarybranchval==0){
		cout<<"Arrays are zombies!"<<endl;
	}
	
	//cout<<"Setting branch addresses"<<endl;
	vtxxBranch->SetAddress(vtxxbranchval);
	vtxyBranch->SetAddress(vtxybranchval);
	vtxzBranch->SetAddress(vtxzbranchval);
	vtxtBranch->SetAddress(vtxtbranchval);
	pxBranch->SetAddress(pxbranchval);
	pyBranch->SetAddress(pybranchval);
	pzBranch->SetAddress(pzbranchval);
	EBranch->SetAddress(ebranchval);
	KEBranch->SetAddress(kebranchval);
	pdgBranch->SetAddress(pdgbranchval);
	nuprimaryBranch->SetAddress(nuprimarybranchval);
	
	//cout<<"Getting primary arrays"<<endl;
	vtxxBranch->GetEntry(localEntry);
	vtxyBranch->GetEntry(localEntry);
	vtxzBranch->GetEntry(localEntry);
	vtxtBranch->GetEntry(localEntry);
	pxBranch->GetEntry(localEntry);
	pyBranch->GetEntry(localEntry);
	pzBranch->GetEntry(localEntry);
	EBranch->GetEntry(localEntry);
	KEBranch->GetEntry(localEntry);
	pdgBranch->GetEntry(localEntry);
	nuprimaryBranch->GetEntry(localEntry);
	
	if(VERBOSE){
	cout<<"processing inputEntry "<<inputEntry<<", localEntry "<<localEntry
	    <<"/"<<thistreesentries<<" in tree "<<treeNumber<<endl;
	}
	
	/* 2. Check if genie primary, and volume is in tank - if not, continue */
	if(strcmp(numatval,"TankWater")!=0){ continue; }	// neutrino intx wasn't in tank water
	
	tankeventnum++;
	for(int i=0;i<ntankbranchval;i++){
		//cout<<"Loading details of primary "<<i<<endl;
		vtxxval=vtxxbranchval[i];	// vtx vals are in cm
		vtxyval=vtxybranchval[i];	// note: nuvtx vals are in m!
		vtxzval=vtxzbranchval[i];
		vtxtval=vtxtbranchval[i];	// times are in nanoseconds
		pxval=pxbranchval[i];		// energies and momenta are in GeV
		pyval=pybranchval[i];
		pzval=pzbranchval[i];
		eval=ebranchval[i];
		keval=kebranchval[i];
		pdgval=pdgbranchval[i];
		nuprimaryval=nuprimarybranchval[i];
		TVector3 thevtx = TVector3(vtxxval, vtxyval, vtxzval);
		TVector3 thepdir = TVector3(pxval, pyval, pzval);
		thepdir.Unit();		// normalise to unit vector
		
		if(keval==0.){keval+=(pow(10.,-6.)); eval+=(pow(10.,-6.)); thepdir=TVector3(0.,0.,1.);}
		if(VERBOSE){
		if(nuprimaryval==1){cout<<"genie primary ";} else {cout<<"genie secondary ";}
		cout<<eval<<" GeV ";
		cout<<"PDG: "<<pdgval;
		cout<<" At: ("<<vtxxval<<", "<<vtxyval<<", "<<vtxzval<<")"<<endl;
		}
		
		vertexxvals->Fill(vtxxval);
		vertexyvals->Fill(vtxyval);
		vertexzvals->Fill(vtxzval);
		if(true){//std::find(goodevents.begin(), goodevents.end(), tankeventnum)!=goodevents.end()){ 
			if(i==0){ cout<<"TANK EVENT "<<tankeventnum<<", NU ("<<nuvtxxval<<", "<<nuvtxyval<<", "<<nuvtxzval<<"), PRIMARIES @"<<endl; }
			cout<<"("<<vtxxval<<", "<<vtxyval<<", "<<vtxzval<<")"<<endl;
		}
	}
	
} // end of loop over TChain entries

TCanvas c1;
c1.cd();
vertexxvals->Draw();
TCanvas c2;
c2.cd();
vertexyvals->Draw();
TCanvas c3;
c3.cd();
vertexzvals->Draw();
}

