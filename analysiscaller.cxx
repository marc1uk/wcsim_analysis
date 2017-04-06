/* vim:set noexpandtab tabstop=4 wrap */
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
#include <time.h>         // clock_t, clock, CLOCKS_PER_SEC
#include "TString.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "wcsimanalysis.hh"

WCSimAnalysis* analysiscaller(TString directoryin){

	//TString theincludepath = gSystem->Getenv("ROOT_INCLUDE_PATH");
	//gInterpreter->AddIncludePath(theincludepath);
	WCSimAnalysis* theanalysis = new WCSimAnalysis(directoryin.Data());
	//WCSimAnalysis* theanalysis = new WCSimAnalysis("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in");
	// no trailing slash in directory
	theanalysis->DoAnalysis();
	cout<<"end of DoAnalysis() call"<<endl;
	//std::this_thread::sleep_for (std::chrono::seconds(5));	// a little wait so we can look at histos
	//delete theanalysis;
	return theanalysis;											// keeps the histograms open

}
