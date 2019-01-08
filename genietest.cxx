/* vim:set noexpandtab tabstop=2 wrap */
//C++
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
//ROOT
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TString.h>
//GENIE
#include <Ntuple/NtpMCEventRecord.h>
#include <FluxDrivers/GSimpleNtpFlux.h>
#include <FluxDrivers/GNuMIFlux.h>
#include <GHEP/GHepUtils.h>               // neut reaction codes
#include <PDG/PDGLibrary.h>
#include <Conventions/Constants.h>

void genietest(){
	TChain* oldflux = new TChain("gtree");
	oldflux->Add("/pnfs/annie/persistent/users/vfischer/genie/BNB_Water_10k_22-05-17/gntp.10??.ghep.root");
	long long oldfluxentries = oldflux->GetEntries();
	cout<<"old flux has "<<oldfluxentries<<" entries"<<endl;
	
	cout<<"setting oldflux branch addresses"<<endl;
	
	genie::flux::GNuMIFluxPassThroughInfo* gnumipassthruentry  = nullptr;
	oldflux->SetBranchAddress("flux",&gnumipassthruentry);
	oldflux->GetBranch("flux")->SetAutoDelete(kTRUE);
	
	genie::NtpMCEventRecord* genieintx = nullptr; //new genie::NtpMCEventRecord;
	oldflux->SetBranchAddress("gmcrec",&genieintx);
	oldflux->GetBranch("gmcrec")->SetAutoDelete(kTRUE);
	
	TFile* curf=nullptr;
	TFile* curflast=nullptr;
	long long maxentriestoanalyse = 2000000000;
	for(int i=0; i<std::min(oldfluxentries,maxentriestoanalyse); i++){
		long long currentevtnum = oldflux->LoadTree(i);
		curf = oldflux->GetCurrentFile();
		if(curf!=curflast || curflast==nullptr){
			TString curftstring = curf->GetName();
			curflast=curf;
			cout<<"reading file "<<curftstring.Data()<<endl;
		}
		if((i%1000)==0) cout<<"i="<<i<<" = "<<(((double)i/(double)oldfluxentries)*100.)<<"%"<<endl;
		oldflux->GetEntry(i);
		
		// neutrino interaction info
		genie::EventRecord* gevtRec = genieintx->event;
		genie::Interaction* genieint = gevtRec->Summary();
		
		genieintx->Clear(); // REQUIRED TO PREVENT MEMORY LEAK
	}
}
