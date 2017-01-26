#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
#include <time.h>         // clock_t, clock, CLOCKS_PER_SEC
#include "TSystem.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "wcsimanalysis.hh"

WCSimAnalysis* analysiscaller(){

	//TString theincludepath = gSystem->Getenv("ROOT_INCLUDE_PATH");
	//gInterpreter->AddIncludePath(theincludepath);
	WCSimAnalysis* theanalysis = new WCSimAnalysis("/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in");	//no trailing slash
	theanalysis->DoAnalysis();
	cout<<"end of DoAnalysis() call"<<endl;
	//std::this_thread::sleep_for (std::chrono::seconds(5));	// a little wait so we can look at histos
	//delete theanalysis;
	return theanalysis;	// keeps the histograms open

}
