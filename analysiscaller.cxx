/* vim:set noexpandtab tabstop=4 wrap */
/*
TO RUN: 
cd to this folder
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:/annie/app/users/moflaher/wcsim/wcsim/include
root analysiscaller.cxx+g
WCSimAnalysis* theana = analysiscaller("/annie/app/users/moflaher/wcsim/root_work/in/MRD_muon_sample")
// This performs the analaysis as `DoAnalysis()` is called in the constructor.
delete theana;
UPDATE:
or just run `root dotheanalysis.C`
*/
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
#include <time.h>         // clock_t, clock, CLOCKS_PER_SEC
#include "TString.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "wcsimanalysis.hh"

WCSimAnalysis* analysiscaller(TString directoryin, TString configfilein, TString directoryout){

	//TString theincludepath = gSystem->Getenv("ROOT_INCLUDE_PATH");
	//gInterpreter->AddIncludePath(theincludepath);
	WCSimAnalysis* theanalysis = new WCSimAnalysis(directoryin.Data(), configfilein.Data(), directoryout.Data());
	//WCSimAnalysis* theanalysis = new WCSimAnalysis("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in");
	// no trailing slash in directory
	cout<<"calling DoAnalysis()"<<endl;
	theanalysis->DoAnalysis();
	cout<<"end of DoAnalysis() call"<<endl;
	//std::this_thread::sleep_for (std::chrono::seconds(5));	// a little wait so we can look at histos
	//delete theanalysis;
	return theanalysis;											// keeps the histograms open

}
